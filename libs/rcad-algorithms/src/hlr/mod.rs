//! Hidden-Line Removal (HLR).
//!
//! Projects a rcad_kernel::BRep's edges onto a view plane and classifies each edge segment
//! as **visible** or **hidden** by testing against the silhouette of all faces.
//!
//! Analytic silhouette curves are generated for curved surfaces (cylinder,
//! sphere, cone, torus) and processed through the same visibility pipeline as wire edges.
//! For general surfaces (BSpline, Bezier, etc.), numerical silhouette extraction
//! is performed using adaptive sampling with curvature-based refinement.
//!
//! Analogous to OCCT `HLRrcad_kernel::BRep_Algo` / `HLRrcad_kernel::BRep_HLRToShape`.
//!
//! # Algorithm
//!
//! For each edge (and silhouette curve):
//! 1. Project both endpoints onto the screen plane.
//! 2. Sample `N` points along the edge in 3D (adaptively if `curvature_adaptive` is true).
//! 3. For each sample, cast a ray from that point toward the camera.
//! 4. If any face triangle blocks the ray **closer** to the camera than the
//!    edge sample, that sample is hidden.
//! 5. Classify runs of consecutive samples → visible/hidden segments.
//!
//! The result is a set of `HlrSegment`s — 2D projected line segments labeled
//! visible or hidden.
//!
//! # Silhouette Classification
//!
//! Segments are classified as:
//! - **Visible**: Not occluded by any face
//! - **Hidden**: Occluded by at least one face
//! - **Contour**: Edge of a face's silhouette (marked via `is_contour`)
//!
//! # Curved Surface Enhancements
//!
//! This implementation includes several enhancements for better handling of curved geometry:
//! - **Curvature-adaptive sampling**: Uses surface curvature to concentrate samples in high-curvature regions
//! - **Marching silhouette extraction**: Robust detection of silhouette curves on general parametric surfaces
//! - **B-spline fitting**: Converts dense silhouette points to smooth B-spline curves
//! - **BVH acceleration**: Spatial acceleration structure for efficient ray casting
//! - **Grazing angle handling**: Special treatment for near-silhouette regions
//! - **Thread edge classification**: Support for helical edges on cylinders and cones
//! - **Seam edge detection**: Proper handling of seam edges on closed surfaces
//! - **Parallel processing**: Multi-threaded processing for large models

use crate::tolerance::*;
use glam::{DAffine3, DMat4, DVec2, DVec3, DVec4};
use rayon::prelude::*;
use rcad_kernel::geom::{Circle3, CurveEval, Surface3, any_perpendicular};
use rcad_kernel::{SurfaceEval, topods};
use std::collections::{HashMap, HashSet};

// ── Public types ──────────────────────────────────────────────────────────────

/// Configuration options for HLR computation.
#[derive(Debug, Clone)]
pub struct HlrOptions {
    /// Number of samples per edge for occlusion testing.
    /// Higher values give more accurate results but slower computation.
    /// Default: 8.
    pub edge_samples: usize,
    /// Base number of samples for silhouette curve generation.
    /// Default: 32.
    pub silhouette_samples: usize,
    /// Enable curvature-adaptive sampling for silhouette curves.
    /// When true, high-curvature regions receive more samples.
    /// Default: true.
    pub curvature_adaptive: bool,
    /// Tolerance for tangent alignment when computing silhouettes.
    /// Points where |normal · view_dir| < tangent_tolerance are considered silhouette candidates.
    /// Default: TOLERANCE_MESH_LEGACY.
    pub tangent_tolerance: f64,
    /// Maximum angle deviation (in radians) for adaptive sampling subdivision.
    /// Smaller values produce smoother curves at higher cost.
    /// Default: 0.05 (about 3 degrees).
    pub angular_tolerance: f64,
    /// Minimum number of subdivision iterations for adaptive sampling.
    /// Default: 2.
    pub min_subdivisions: usize,
    /// Maximum number of subdivision iterations for adaptive sampling.
    /// Default: 8.
    pub max_subdivisions: usize,
    /// Enable BVH acceleration for ray casting.
    /// Default: true.
    pub use_bvh: bool,
    /// Maximum curvature for adaptive sampling (higher = more samples in curved regions).
    /// Default: 100.0.
    pub max_curvature: f64,
    /// Minimum curvature for adaptive sampling (lower = fewer samples in flat regions).
    /// Default: 0.001.
    pub min_curvature: f64,
    /// Enable B-spline fitting for silhouette curves.
    /// Default: true.
    pub fit_bspline: bool,
    /// Tolerance for B-spline fitting (maximum deviation from original points).
    /// Default: 0.001.
    pub bspline_tolerance: f64,
    /// Grazing angle threshold (in radians). Points closer to silhouette receive special handling.
    /// Default: 0.1 (about 6 degrees).
    pub grazing_angle_threshold: f64,
    /// Enable smooth silhouette approximation.
    /// Default: true.
    pub smooth_silhouettes: bool,
    /// Enable parallel processing for multi-face models.
    /// Default: true.
    pub parallel: bool,
    /// Minimum number of faces to trigger parallel processing.
    /// Default: 4.
    pub parallel_threshold: usize,
    /// Enable surface property caching for improved performance.
    /// Default: true.
    pub cache_surface_properties: bool,
    /// Silhouette proximity factor for increased sampling density.
    /// Samples within this factor of the silhouette receive more refinement.
    /// Default: 0.1 (10% of local feature size).
    pub silhouette_proximity_factor: f64,
    /// Enable thread edge detection for helical geometry.
    /// Default: true.
    pub detect_thread_edges: bool,
    /// Enable seam edge detection for closed surfaces.
    /// Default: true.
    pub detect_seam_edges: bool,
    /// Maximum depth complexity for curve-surface intersection.
    /// Default: 16.
    pub max_depth_complexity: usize,
}

impl Default for HlrOptions {
    fn default() -> Self {
        Self {
            edge_samples: 8,
            silhouette_samples: 32,
            curvature_adaptive: true,
            tangent_tolerance: TOLERANCE_MESH_LEGACY,
            angular_tolerance: 0.05,
            min_subdivisions: 2,
            max_subdivisions: 8,
            use_bvh: true,
            max_curvature: 100.0,
            min_curvature: 0.001,
            fit_bspline: true,
            bspline_tolerance: 0.001,
            grazing_angle_threshold: 0.1,
            smooth_silhouettes: true,
            parallel: true,
            parallel_threshold: 4,
            cache_surface_properties: true,
            silhouette_proximity_factor: 0.1,
            detect_thread_edges: true,
            detect_seam_edges: true,
            max_depth_complexity: 16,
        }
    }
}

impl HlrOptions {
    /// Create options with a specific edge sample count.
    pub fn with_edge_samples(mut self, n: usize) -> Self {
        self.edge_samples = n.max(2);
        self
    }

    /// Create options with a specific silhouette sample count.
    pub fn with_silhouette_samples(mut self, n: usize) -> Self {
        self.silhouette_samples = n.max(8);
        self
    }

    /// Enable or disable curvature-adaptive sampling.
    pub fn with_curvature_adaptive(mut self, adaptive: bool) -> Self {
        self.curvature_adaptive = adaptive;
        self
    }

    /// Set the tangent tolerance for silhouette detection.
    pub fn with_tangent_tolerance(mut self, tol: f64) -> Self {
        self.tangent_tolerance = tol.abs().max(TOLERANCE_LEN_MIN);
        self
    }

    /// Enable or disable BVH acceleration.
    pub fn with_bvh(mut self, use_bvh: bool) -> Self {
        self.use_bvh = use_bvh;
        self
    }

    /// Set the maximum curvature for adaptive sampling.
    pub fn with_max_curvature(mut self, curv: f64) -> Self {
        self.max_curvature = curv.abs().max(0.1);
        self
    }

    /// Enable or disable B-spline fitting for silhouettes.
    pub fn with_bspline_fitting(mut self, fit: bool) -> Self {
        self.fit_bspline = fit;
        self
    }

    /// Set the grazing angle threshold.
    pub fn with_grazing_angle(mut self, angle: f64) -> Self {
        self.grazing_angle_threshold = angle.abs().min(std::f64::consts::FRAC_PI_2);
        self
    }

    /// Enable or disable smooth silhouette approximation.
    pub fn with_smooth_silhouettes(mut self, smooth: bool) -> Self {
        self.smooth_silhouettes = smooth;
        self
    }

    /// Enable or disable parallel processing.
    pub fn with_parallel(mut self, parallel: bool) -> Self {
        self.parallel = parallel;
        self
    }

    /// Set the parallel processing threshold (minimum faces to trigger parallelism).
    pub fn with_parallel_threshold(mut self, threshold: usize) -> Self {
        self.parallel_threshold = threshold.max(1);
        self
    }

    /// Enable or disable surface property caching.
    pub fn with_surface_caching(mut self, cache: bool) -> Self {
        self.cache_surface_properties = cache;
        self
    }

    /// Set the silhouette proximity factor.
    pub fn with_silhouette_proximity(mut self, factor: f64) -> Self {
        self.silhouette_proximity_factor = factor.abs().max(0.01).min(1.0);
        self
    }

    /// Enable or disable thread edge detection.
    pub fn with_thread_edge_detection(mut self, detect: bool) -> Self {
        self.detect_thread_edges = detect;
        self
    }

    /// Enable or disable seam edge detection.
    pub fn with_seam_edge_detection(mut self, detect: bool) -> Self {
        self.detect_seam_edges = detect;
        self
    }
}

/// Hint about the geometric type of the original 3D edge curve.
/// Used by consumers (e.g. SVG exporter) to emit arcs instead of polylines.
#[derive(Debug, Clone, PartialEq)]
pub enum CurveHint {
    /// Edge is a full or partial circle in 3D.
    Circle {
        /// Projected 2D center of the circle.
        center: DVec2,
        /// Projected radius (approximate — perspective not applied).
        radius: f64,
    },
    /// Any other non-straight curve (ellipse, spline, …).
    Other,
}

/// Classification of an HLR segment type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentType {
    /// Regular edge (part of the rcad_kernel::BRep wire).
    Edge,
    /// Silhouette curve (contour of a curved face).
    Silhouette,
    /// Thread edge (helical edge on cylinders/cones).
    Thread,
    /// Seam edge (closed surface seam).
    Seam,
}

/// Classification of edge visibility type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeClassification {
    /// Edge is fully visible.
    Visible,
    /// Edge is fully hidden.
    Hidden,
    /// Edge is a contour/silhouette edge.
    Contour,
    /// Edge is partially visible (some segments visible, some hidden).
    Partial,
    /// Edge is a thread edge (helical).
    Thread,
    /// Edge is a seam edge on a closed surface.
    Seam,
}

/// Information about an edge's classification for HLR.
#[derive(Debug, Clone)]
pub struct EdgeClassInfo {
    /// Edge index in the rcad_kernel::BRep.
    pub edge_idx: usize,
    /// Classification type.
    pub classification: EdgeClassification,
    /// Number of visible segments.
    pub visible_segments: usize,
    /// Number of hidden segments.
    pub hidden_segments: usize,
    /// Whether this edge is on a curved surface.
    pub on_curved_surface: bool,
    /// Surface index if on a curved surface.
    pub surface_idx: Option<usize>,
}

/// A projected edge segment labeled as visible or hidden.
#[derive(Debug, Clone, PartialEq)]
pub struct HlrSegment {
    /// Start point in 2D screen space.
    pub start: DVec2,
    /// End point in 2D screen space.
    pub end: DVec2,
    /// Whether this segment is visible from the camera.
    pub visible: bool,
    /// Optional hint about the underlying curve type (None for straight lines).
    pub curve_hint: Option<CurveHint>,
    /// Segment type (edge or silhouette).
    pub segment_type: SegmentType,
}

impl HlrSegment {
    /// Returns true if this segment is a silhouette/contour curve.
    pub fn is_contour(&self) -> bool {
        self.segment_type == SegmentType::Silhouette
    }

    /// Returns true if this segment is a thread edge.
    pub fn is_thread(&self) -> bool {
        self.segment_type == SegmentType::Thread
    }

    /// Returns true if this segment is a seam edge.
    pub fn is_seam(&self) -> bool {
        self.segment_type == SegmentType::Seam
    }
}

// ── Surface Normal Analysis ─────────────────────────────────────────────────────

/// Cached surface properties for efficient silhouette computation.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceProperties {
    /// Surface point at (u, v).
    pub point: DVec3,
    /// Unit normal at (u, v).
    pub normal: DVec3,
    /// Principal curvatures (k1, k2).
    pub curvatures: (f64, f64),
    /// Gaussian curvature.
    pub gaussian: f64,
    /// Mean curvature.
    pub mean: f64,
}

impl SurfaceProperties {
    /// Compute the dot product of the normal with a view direction.
    #[inline]
    pub fn normal_dot_view(&self, view_dir: DVec3) -> f64 {
        self.normal.dot(view_dir)
    }

    /// Check if this point is near a silhouette (normal nearly perpendicular to view).
    #[inline]
    pub fn is_near_silhouette(&self, view_dir: DVec3, threshold: f64) -> bool {
        self.normal_dot_view(view_dir).abs() < threshold
    }

    /// Get the maximum principal curvature magnitude.
    #[inline]
    pub fn max_curvature(&self) -> f64 {
        self.curvatures.0.abs().max(self.curvatures.1.abs())
    }

    /// Check if the surface is locally flat (low curvature).
    #[inline]
    pub fn is_flat(&self, tolerance: f64) -> bool {
        self.max_curvature() < tolerance
    }
}

/// Cache for surface property evaluations.
#[derive(Debug, Clone)]
pub struct SurfacePropertyCache {
    /// Cached properties keyed by (u, v) discretized to grid cells.
    cache: HashMap<(usize, usize), SurfaceProperties>,
    /// Grid resolution for cache.
    resolution: usize,
    /// UV domain of the surface.
    domain: [f64; 4],
}

impl SurfacePropertyCache {
    /// Create a new cache with given resolution.
    pub fn new(resolution: usize, domain: [f64; 4]) -> Self {
        Self {
            cache: HashMap::new(),
            resolution,
            domain,
        }
    }

    /// Get or compute surface properties at (u, v).
    pub fn get_or_compute(&mut self, surface: &Surface3, u: f64, v: f64) -> SurfaceProperties {
        let [u0, u1, v0, v1] = self.domain;
        let i = ((u - u0) / (u1 - u0) * self.resolution as f64).min(self.resolution as f64 - 1.0)
            as usize;
        let j = ((v - v0) / (v1 - v0) * self.resolution as f64).min(self.resolution as f64 - 1.0)
            as usize;

        if let Some(&props) = self.cache.get(&(i, j)) {
            return props;
        }

        let props = compute_surface_properties(surface, u, v);
        self.cache.insert((i, j), props);
        props
    }

    /// Get cached properties if available.
    pub fn get(&self, u: f64, v: f64) -> Option<&SurfaceProperties> {
        let [u0, u1, v0, v1] = self.domain;
        let i = ((u - u0) / (u1 - u0) * self.resolution as f64).min(self.resolution as f64 - 1.0)
            as usize;
        let j = ((v - v0) / (v1 - v0) * self.resolution as f64).min(self.resolution as f64 - 1.0)
            as usize;
        self.cache.get(&(i, j))
    }

    /// Clear the cache.
    pub fn clear(&mut self) {
        self.cache.clear();
    }

    /// Get the number of cached entries.
    pub fn len(&self) -> usize {
        self.cache.len()
    }

    /// Check if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }
}

/// Compute surface properties at a given parameter location.
pub fn compute_surface_properties(surface: &Surface3, u: f64, v: f64) -> SurfaceProperties {
    let point = surface.point_at(u, v);
    let normal = surface.normal_at(u, v);
    let curvatures = rcad_kernel::curvature::principal_curvatures(surface, u, v);
    let gaussian = rcad_kernel::curvature::gaussian_curvature(surface, u, v);
    let mean = rcad_kernel::curvature::mean_curvature(surface, u, v);

    SurfaceProperties {
        point,
        normal,
        curvatures,
        gaussian,
        mean,
    }
}

// ── Adaptive Silhouette Sampling ────────────────────────────────────────────────

/// Sample point along a silhouette curve with additional metadata.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveSample {
    /// Parameter space location (u, v).
    pub uv: (f64, f64),
    /// World space position.
    pub point: DVec3,
    /// Surface normal at this point.
    pub normal: DVec3,
    /// Maximum curvature at this point.
    pub curvature: f64,
    /// Distance to the exact silhouette curve (0 = on silhouette).
    pub silhouette_distance: f64,
    /// Sampling weight (higher = more samples nearby).
    pub weight: f64,
}

/// Adaptive sampling configuration for silhouette curves.
#[derive(Debug, Clone)]
pub struct AdaptiveSamplingConfig {
    /// Base number of samples.
    pub base_samples: usize,
    /// Maximum number of samples after adaptive refinement.
    pub max_samples: usize,
    /// Curvature threshold for refinement (higher curvature = more samples).
    pub curvature_threshold: f64,
    /// Proximity threshold for refinement (closer to silhouette = more samples).
    pub proximity_threshold: f64,
    /// Minimum chord length between samples.
    pub min_chord_length: f64,
    /// Maximum angle deviation between consecutive samples (radians).
    pub max_angle_deviation: f64,
}

impl Default for AdaptiveSamplingConfig {
    fn default() -> Self {
        Self {
            base_samples: 32,
            max_samples: 256,
            curvature_threshold: 10.0,
            proximity_threshold: 0.05,
            min_chord_length: TOLERANCE_RETRY_LADDER_COARSE,
            max_angle_deviation: 0.1,
        }
    }
}

/// Compute adaptive samples along a silhouette curve.
pub fn compute_adaptive_samples(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    config: &AdaptiveSamplingConfig,
    opts: &HlrOptions,
) -> Vec<AdaptiveSample> {
    let [u0, u1, v0, v1] = domain;

    // Find silhouette seed points
    let seeds = find_silhouette_seeds(
        surface,
        view_dir,
        domain,
        config.base_samples,
        opts.tangent_tolerance,
    );

    if seeds.is_empty() {
        return Vec::new();
    }

    // Trace silhouette curves from seeds
    let mut all_samples: Vec<AdaptiveSample> = Vec::new();
    let mut visited: HashSet<(usize, usize)> = HashSet::new();

    for (_, _, u, v) in seeds {
        // Check if this cell was already visited
        let cell_i = ((u - u0) / (u1 - u0) * config.base_samples as f64) as usize;
        let cell_j = ((v - v0) / (v1 - v0) * config.base_samples as f64) as usize;

        if visited.contains(&(cell_i, cell_j)) {
            continue;
        }

        // Trace the silhouette curve
        let curve_samples =
            trace_adaptive_silhouette(surface, view_dir, domain, u, v, config, opts);

        // Mark visited cells
        for sample in &curve_samples {
            let ci = ((sample.uv.0 - u0) / (u1 - u0) * config.base_samples as f64) as usize;
            let cj = ((sample.uv.1 - v0) / (v1 - v0) * config.base_samples as f64) as usize;
            visited.insert((
                ci.min(config.base_samples - 1),
                cj.min(config.base_samples - 1),
            ));
        }

        all_samples.extend(curve_samples);
    }

    // Refine samples in high-curvature and near-silhouette regions
    if opts.curvature_adaptive {
        refine_adaptive_samples(surface, view_dir, &mut all_samples, config, opts);
    }

    all_samples
}

/// Trace a silhouette curve with adaptive sampling.
fn trace_adaptive_silhouette(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    u_start: f64,
    v_start: f64,
    config: &AdaptiveSamplingConfig,
    opts: &HlrOptions,
) -> Vec<AdaptiveSample> {
    let mut samples: Vec<AdaptiveSample> = Vec::new();
    let [u0, u1, v0, v1] = domain;

    // Add the starting point
    if let Some(sample) = create_adaptive_sample(surface, view_dir, u_start, v_start) {
        samples.push(sample);
    }

    // March in both directions from the seed
    for direction in &[-1.0_f64, 1.0] {
        let mut u = u_start;
        let mut v = v_start;
        let mut prev_sample = samples.first().copied();

        for _ in 0..opts.max_subdivisions * 100 {
            // Compute the tangent direction to the silhouette curve
            let tangent = compute_silhouette_tangent(surface, view_dir, u, v);

            if tangent.length_squared() < 1e-16 {
                break;
            }

            // Choose direction along the tangent
            let step_dir = *direction * tangent.normalize_or_zero();

            // Compute adaptive step size based on curvature
            let props = compute_surface_properties(surface, u, v);
            let max_k = props.max_curvature().max(opts.min_curvature);
            let curvature_factor = (opts.max_curvature / max_k).min(4.0).max(0.25);
            let step_size = opts.angular_tolerance * curvature_factor;

            // Take a step
            let u_new = u + step_dir.x * step_size;
            let v_new = v + step_dir.y * step_size;

            // Check bounds
            if u_new < u0 || u_new > u1 || v_new < v0 || v_new > v1 {
                break;
            }

            // Project back onto the silhouette curve
            if let Some((u_proj, v_proj)) =
                project_to_silhouette(surface, view_dir, u_new, v_new, opts.tangent_tolerance)
            {
                u = u_proj;
                v = v_proj;

                // Create a new sample
                if let Some(sample) = create_adaptive_sample(surface, view_dir, u, v) {
                    // Check if we've moved enough to add a new sample
                    if let Some(prev) = prev_sample {
                        let dist = (sample.point - prev.point).length();
                        if dist < config.min_chord_length {
                            continue;
                        }

                        // Check angular deviation
                        let dir_new = (sample.point - prev.point).normalize_or_zero();
                        let dir_prev = if samples.len() >= 2 {
                            (samples[samples.len() - 1].point - samples[samples.len() - 2].point)
                                .normalize_or_zero()
                        } else {
                            dir_new
                        };
                        let angle = dir_new.dot(dir_prev).acos().abs();
                        if angle > config.max_angle_deviation {
                            // Add intermediate samples for high angular deviation
                            add_intermediate_samples(
                                surface,
                                view_dir,
                                prev,
                                sample,
                                &mut samples,
                                config,
                            );
                        }
                    }

                    samples.push(sample);
                    prev_sample = Some(sample);
                }
            } else {
                break;
            }

            // Check for closed loop
            if samples.len() > 10 {
                let first = samples[0];
                let dist = ((first.uv.0 - u).powi(2) + (first.uv.1 - v).powi(2)).sqrt();
                if dist < step_size * 2.0 {
                    break;
                }
            }
        }
    }

    samples
}

/// Create an adaptive sample at a parameter location.
fn create_adaptive_sample(
    surface: &Surface3,
    view_dir: DVec3,
    u: f64,
    v: f64,
) -> Option<AdaptiveSample> {
    let point = surface.point_at(u, v);
    let normal = surface.normal_at(u, v);
    let curvatures = rcad_kernel::curvature::principal_curvatures(surface, u, v);
    let curvature = curvatures.0.abs().max(curvatures.1.abs());

    // Compute distance to exact silhouette (absolute value of normal dot view)
    let silhouette_distance = normal.dot(view_dir).abs();

    // Compute sampling weight based on curvature and proximity
    let weight = (curvature + 1.0) * (silhouette_distance + 0.1).recip();

    Some(AdaptiveSample {
        uv: (u, v),
        point,
        normal,
        curvature,
        silhouette_distance,
        weight,
    })
}

/// Add intermediate samples between two samples for smooth curves.
fn add_intermediate_samples(
    surface: &Surface3,
    view_dir: DVec3,
    start: AdaptiveSample,
    end: AdaptiveSample,
    samples: &mut Vec<AdaptiveSample>,
    config: &AdaptiveSamplingConfig,
) {
    let num_intermediate =
        ((end.point - start.point).length() / config.min_chord_length).ceil() as usize;
    let num_intermediate = num_intermediate.min(4);

    for i in 1..num_intermediate {
        let t = i as f64 / num_intermediate as f64;
        let u = start.uv.0 + t * (end.uv.0 - start.uv.0);
        let v = start.uv.1 + t * (end.uv.1 - start.uv.1);

        if let Some(sample) = create_adaptive_sample(surface, view_dir, u, v) {
            samples.push(sample);
        }
    }
}

/// Refine adaptive samples based on curvature and silhouette proximity.
fn refine_adaptive_samples(
    surface: &Surface3,
    view_dir: DVec3,
    samples: &mut Vec<AdaptiveSample>,
    config: &AdaptiveSamplingConfig,
    _opts: &HlrOptions,
) {
    if samples.len() < 2 {
        return;
    }

    let mut refined: Vec<AdaptiveSample> = Vec::with_capacity(samples.len() * 2);
    refined.push(samples[0]);

    for i in 1..samples.len() {
        let prev = &samples[i - 1];
        let curr = &samples[i];

        // Determine if refinement is needed based on curvature and proximity
        let chord_len = (curr.point - prev.point).length();
        let avg_curvature = (prev.curvature + curr.curvature) * 0.5;
        let avg_proximity = (prev.silhouette_distance + curr.silhouette_distance) * 0.5;

        let needs_refinement = avg_curvature > config.curvature_threshold
            || avg_proximity < config.proximity_threshold
            || chord_len > config.min_chord_length * 4.0;

        if needs_refinement {
            let num_subdivisions =
                (avg_curvature * chord_len / config.curvature_threshold).ceil() as usize;
            let num_subdivisions = num_subdivisions.min(4).max(1);

            for j in 1..num_subdivisions {
                let t = j as f64 / num_subdivisions as f64;
                let u = prev.uv.0 + t * (curr.uv.0 - prev.uv.0);
                let v = prev.uv.1 + t * (curr.uv.1 - prev.uv.1);

                if let Some(sample) = create_adaptive_sample(surface, view_dir, u, v) {
                    refined.push(sample);
                }
            }
        }

        refined.push(*curr);
    }

    *samples = refined;
}

/// Output of an HLR computation.
#[derive(Debug, Clone, Default)]
pub struct HlrResult {
    pub segments: Vec<HlrSegment>,
}

impl HlrResult {
    /// Return only visible segments.
    pub fn visible(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible)
    }

    /// Return only hidden segments.
    pub fn hidden(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| !s.visible)
    }

    /// Return only silhouette/contour segments.
    pub fn silhouettes(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.is_contour())
    }

    /// Return only visible silhouette segments.
    pub fn visible_silhouettes(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible && s.is_contour())
    }
}

/// Camera / view specification for HLR.
#[derive(Debug, Clone)]
pub struct HlrCamera {
    /// Camera position in world space.
    pub eye: DVec3,
    /// Target point (look-at).
    pub target: DVec3,
    /// Up direction.
    pub up: DVec3,
}

impl HlrCamera {
    pub fn new(eye: DVec3, target: DVec3) -> Self {
        Self {
            eye,
            target,
            up: DVec3::Y,
        }
    }

    pub fn with_up(mut self, up: DVec3) -> Self {
        self.up = up;
        self
    }

    /// Isometric-style view from the +X+Y+Z octant.
    pub fn isometric(distance: f64) -> Self {
        let d = distance / 3.0_f64.sqrt();
        Self::new(DVec3::splat(d), DVec3::ZERO)
    }

    /// Front view (looking along +Y, up = +Z).
    pub fn front(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, -distance, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }

    /// Top view (looking down -Z).
    pub fn top(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, 0.0, distance), DVec3::ZERO).with_up(DVec3::Y)
    }

    /// Right-side view (looking along -X, up = +Z).
    pub fn right(distance: f64) -> Self {
        Self::new(DVec3::new(distance, 0.0, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Build a right-handed view matrix (world → camera space).
fn look_at(eye: DVec3, target: DVec3, up: DVec3) -> DMat4 {
    let forward = (target - eye).normalize_or_zero();
    let right = forward.cross(up).normalize_or_zero();
    let up = right.cross(forward);

    DMat4::from_cols(
        DVec4::new(right.x, right.y, right.z, -right.dot(eye)),
        DVec4::new(up.x, up.y, up.z, -up.dot(eye)),
        DVec4::new(-forward.x, -forward.y, -forward.z, forward.dot(eye)),
        DVec4::new(0.0, 0.0, 0.0, 1.0),
    )
    .transpose()
}

/// Project a world-space point to 2D screen space using the view matrix.
/// Returns (x, y) in camera space (z is depth; ignored for 2D output).
fn project(p: DVec3, view: &DMat4) -> (DVec2, f64) {
    let hp = view.mul_vec4(DVec4::new(p.x, p.y, p.z, 1.0));
    (DVec2::new(hp.x, hp.y), hp.z)
}

/// Collect all triangles from a rcad_kernel::BRep (fan-triangulate faces without pre-triangulated data).
fn collect_triangles(brep: &rcad_kernel::BRep) -> Vec<[DVec3; 3]> {
    let mut tris = Vec::new();
    for solid in &brep.solids() {
        for shell in &solid.shells {
            for face in &shell.faces {
                if !face.triangles.is_empty() {
                    for &[i, j, k] in &face.triangles {
                        if let (Some(a), Some(b), Some(c)) = (
                            brep.vertices().get(i),
                            brep.vertices().get(j),
                            brep.vertices().get(k),
                        ) {
                            tris.push([a.point, b.point, c.point]);
                        }
                    }
                } else {
                    // Fan-triangulate from wire
                    let edges = brep.edges();
                    let vertices = brep.vertices();
                    let pts: Vec<DVec3> = face
                        .outer_wire
                        .edges
                        .iter()
                        .filter_map(|we| {
                            let edge = edges.get(we.idx)?;
                            let vi = if we.forward { edge.start } else { edge.end };
                            vertices.get(vi).map(|v| v.point)
                        })
                        .collect();
                    if pts.len() >= 3 {
                        let origin = pts[0];
                        for i in 1..pts.len() - 1 {
                            tris.push([origin, pts[i], pts[i + 1]]);
                        }
                    }
                }
            }
        }
    }
    tris
}

/// Ray-triangle intersection (Möller–Trumbore). Returns `Some(t)` if the ray
/// `origin + t*dir` hits the triangle (t > epsilon, front-face only).
fn ray_triangle_intersect(origin: DVec3, dir: DVec3, tri: &[DVec3; 3]) -> Option<f64> {
    const EPS: f64 = TOLERANCE_LINEAR_RELAX_8;
    let edge1 = tri[1] - tri[0];
    let edge2 = tri[2] - tri[0];
    let h = dir.cross(edge2);
    let a = edge1.dot(h);
    if a.abs() < EPS {
        return None;
    }
    let f = 1.0 / a;
    let s = origin - tri[0];
    let u = f * s.dot(h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(edge1);
    let v = f * dir.dot(q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * edge2.dot(q);
    if t > EPS { Some(t) } else { None }
}

/// Test if a world-space point is occluded by any triangle when viewed from `eye`.
fn is_occluded(point: DVec3, eye: DVec3, triangles: &[[DVec3; 3]], dist_to_eye: f64) -> bool {
    let dir = (eye - point).normalize_or_zero();
    let origin = point + dir * TOLERANCE_RETRY_LADDER_MID; // push off surface
    for tri in triangles {
        if let Some(t) = ray_triangle_intersect(origin, dir, tri)
            && t < dist_to_eye - TOLERANCE_RETRY_LADDER_COARSE
        {
            return true;
        }
    }
    false
}

// ── Triangle BVH for accelerated ray casting ───────────────────────────────────

/// Axis-aligned bounding box for triangle BVH.
#[derive(Debug, Clone, Copy)]
struct TriAabb {
    min: DVec3,
    max: DVec3,
}

impl TriAabb {
    fn empty() -> Self {
        Self {
            min: DVec3::splat(f64::INFINITY),
            max: DVec3::splat(f64::NEG_INFINITY),
        }
    }

    fn from_triangle(tri: &[DVec3; 3]) -> Self {
        let mut aabb = Self::empty();
        for &p in tri {
            aabb.expand_point(p);
        }
        aabb
    }

    fn expand_point(&mut self, p: DVec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    fn expand_aabb(&mut self, other: &TriAabb) {
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
    }

    fn center(&self) -> DVec3 {
        (self.min + self.max) * 0.5
    }

    fn surface_area(&self) -> f64 {
        let d = self.max - self.min;
        if d.x < 0.0 || d.y < 0.0 || d.z < 0.0 {
            return 0.0;
        }
        2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
    }

    fn ray_intersect(&self, origin: DVec3, inv_dir: DVec3) -> Option<f64> {
        let t1 = (self.min - origin) * inv_dir;
        let t2 = (self.max - origin) * inv_dir;

        let t_min = t1.min(t2);
        let t_max = t1.max(t2);

        let t_enter = t_min.x.max(t_min.y).max(t_min.z);
        let t_exit = t_max.x.min(t_max.y).min(t_max.z);

        if t_exit >= t_enter.max(0.0) {
            Some(t_enter.max(0.0))
        } else {
            None
        }
    }
}

/// BVH node for triangle-level acceleration.
#[derive(Debug, Clone)]
enum TriBvhNode {
    Leaf {
        aabb: TriAabb,
        /// Triangle indices (index into the original triangle array).
        tris: Vec<usize>,
    },
    Internal {
        aabb: TriAabb,
        left: usize,
        right: usize,
    },
}

/// Triangle-level BVH for efficient ray casting.
#[derive(Debug, Clone)]
pub struct TriBvh {
    nodes: Vec<TriBvhNode>,
    triangle_aabbs: Vec<TriAabb>,
    triangle_centers: Vec<DVec3>,
}

const MAX_TRIS_PER_LEAF: usize = 8;
const SAH_BUCKETS: usize = 8;

impl TriBvh {
    /// Build a triangle BVH from a list of triangles.
    pub fn build(triangles: &[[DVec3; 3]]) -> Self {
        if triangles.is_empty() {
            return Self {
                nodes: Vec::new(),
                triangle_aabbs: Vec::new(),
                triangle_centers: Vec::new(),
            };
        }

        let triangle_aabbs: Vec<TriAabb> = triangles.iter().map(TriAabb::from_triangle).collect();
        let triangle_centers: Vec<DVec3> = triangle_aabbs.iter().map(|a| a.center()).collect();

        let tri_indices: Vec<usize> = (0..triangles.len()).collect();

        let mut bvh = Self {
            nodes: Vec::new(),
            triangle_aabbs,
            triangle_centers,
        };

        bvh.build_recursive(&tri_indices);
        bvh
    }

    fn build_recursive(&mut self, tri_indices: &[usize]) -> usize {
        let count = tri_indices.len();
        if count == 0 {
            return usize::MAX;
        }

        // Compute AABB for this node
        let mut aabb = TriAabb::empty();
        for &ti in tri_indices {
            aabb.expand_aabb(&self.triangle_aabbs[ti]);
        }

        // Leaf condition
        if count <= MAX_TRIS_PER_LEAF {
            let node_idx = self.nodes.len();
            self.nodes.push(TriBvhNode::Leaf {
                aabb,
                tris: tri_indices.to_vec(),
            });
            return node_idx;
        }

        // SAH split
        let (split_axis, split_pos) = self.sah_split(tri_indices, &aabb);

        // Partition triangles
        let (left_tris, right_tris): (Vec<usize>, Vec<usize>) =
            tri_indices.iter().copied().partition(|&ti| {
                let center = match split_axis {
                    0 => self.triangle_centers[ti].x,
                    1 => self.triangle_centers[ti].y,
                    _ => self.triangle_centers[ti].z,
                };
                center < split_pos
            });

        // Handle degenerate split
        let (left_tris, right_tris) = if left_tris.is_empty() || right_tris.is_empty() {
            let mid = count / 2;
            let mut sorted = tri_indices.to_vec();
            sorted.sort_by(|&a, &b| {
                let ca = match split_axis {
                    0 => self.triangle_centers[a].x,
                    1 => self.triangle_centers[a].y,
                    _ => self.triangle_centers[a].z,
                };
                let cb = match split_axis {
                    0 => self.triangle_centers[b].x,
                    1 => self.triangle_centers[b].y,
                    _ => self.triangle_centers[b].z,
                };
                ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
            });
            (sorted[..mid].to_vec(), sorted[mid..].to_vec())
        } else {
            (left_tris, right_tris)
        };

        let node_idx = self.nodes.len();
        self.nodes.push(TriBvhNode::Internal {
            aabb: TriAabb::empty(),
            left: 0,
            right: 0,
        });

        let left = self.build_recursive(&left_tris);
        let right = self.build_recursive(&right_tris);

        self.nodes[node_idx] = TriBvhNode::Internal { aabb, left, right };
        node_idx
    }

    fn sah_split(&self, tri_indices: &[usize], parent_aabb: &TriAabb) -> (usize, f64) {
        let parent_sa = parent_aabb.surface_area().max(TOLERANCE_LEN_SQ_DIV_SAFE);
        let mut best_cost = f64::INFINITY;
        let mut best_axis = 0usize;
        let mut best_pos = 0.0f64;

        for axis in 0..3usize {
            let axis_min = match axis {
                0 => parent_aabb.min.x,
                1 => parent_aabb.min.y,
                _ => parent_aabb.min.z,
            };
            let axis_max = match axis {
                0 => parent_aabb.max.x,
                1 => parent_aabb.max.y,
                _ => parent_aabb.max.z,
            };
            let span = axis_max - axis_min;
            if span < TOLERANCE_FLOAT_LOOSE {
                continue;
            }

            for b in 1..SAH_BUCKETS {
                let split = axis_min + span * b as f64 / SAH_BUCKETS as f64;

                let mut left_aabb = TriAabb::empty();
                let mut right_aabb = TriAabb::empty();
                let mut left_count = 0usize;
                let mut right_count = 0usize;

                for &ti in tri_indices {
                    let center_val = match axis {
                        0 => self.triangle_centers[ti].x,
                        1 => self.triangle_centers[ti].y,
                        _ => self.triangle_centers[ti].z,
                    };
                    if center_val < split {
                        left_aabb.expand_aabb(&self.triangle_aabbs[ti]);
                        left_count += 1;
                    } else {
                        right_aabb.expand_aabb(&self.triangle_aabbs[ti]);
                        right_count += 1;
                    }
                }

                if left_count == 0 || right_count == 0 {
                    continue;
                }

                let cost = (left_count as f64 * left_aabb.surface_area()
                    + right_count as f64 * right_aabb.surface_area())
                    / parent_sa;

                if cost < best_cost {
                    best_cost = cost;
                    best_axis = axis;
                    best_pos = split;
                }
            }
        }

        if best_cost.is_infinite() {
            let d = parent_aabb.max - parent_aabb.min;
            best_axis = if d.x >= d.y && d.x >= d.z {
                0
            } else if d.y >= d.z {
                1
            } else {
                2
            };
            best_pos = parent_aabb.center()[best_axis];
        }

        (best_axis, best_pos)
    }

    /// Test if a point is occluded by any triangle in the BVH.
    pub fn is_occluded(
        &self,
        point: DVec3,
        eye: DVec3,
        triangles: &[[DVec3; 3]],
        dist_to_eye: f64,
    ) -> bool {
        if self.nodes.is_empty() {
            return false;
        }

        let dir = (eye - point).normalize_or_zero();
        let origin = point + dir * TOLERANCE_RETRY_LADDER_MID;
        let inv_dir = DVec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);

        self.is_occluded_node(0, origin, dir, inv_dir, triangles, dist_to_eye)
    }

    fn is_occluded_node(
        &self,
        node_idx: usize,
        origin: DVec3,
        dir: DVec3,
        inv_dir: DVec3,
        triangles: &[[DVec3; 3]],
        dist_to_eye: f64,
    ) -> bool {
        let node = &self.nodes[node_idx];

        // Check AABB intersection
        let t_aabb = match node.aabb().ray_intersect(origin, inv_dir) {
            Some(t) => t,
            None => return false,
        };

        // Early exit if AABB is beyond the eye
        if t_aabb > dist_to_eye {
            return false;
        }

        match node {
            TriBvhNode::Leaf { tris, .. } => {
                for &ti in tris {
                    if let Some(t) = ray_triangle_intersect(origin, dir, &triangles[ti])
                        && t < dist_to_eye - TOLERANCE_RETRY_LADDER_COARSE
                    {
                        return true;
                    }
                }
                false
            }
            TriBvhNode::Internal { left, right, .. } => {
                self.is_occluded_node(*left, origin, dir, inv_dir, triangles, dist_to_eye)
                    || self.is_occluded_node(*right, origin, dir, inv_dir, triangles, dist_to_eye)
            }
        }
    }
}

impl TriBvhNode {
    fn aabb(&self) -> &TriAabb {
        match self {
            TriBvhNode::Leaf { aabb, .. } => aabb,
            TriBvhNode::Internal { aabb, .. } => aabb,
        }
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

// ── Silhouette generation ─────────────────────────────────────────────────────

/// Internal: one silhouette curve to process through the HLR pipeline.
struct SilhouetteCurve {
    /// World-space sample points (at least 2).
    world_pts: Vec<DVec3>,
    /// Optional curve hint for SVG output.
    curve_hint: Option<CurveHint>,
    /// If true, emit one segment per consecutive point pair instead of merging
    /// runs.  Used for dense polyline approximations (e.g. sphere silhouette).
    dense: bool,
}

/// A 3D silhouette curve extracted from a curved surface.
#[derive(Debug, Clone)]
pub struct SilhouetteCurve3 {
    /// World-space sample points along the silhouette curve.
    pub points: Vec<DVec3>,
    /// The surface index from which this silhouette was extracted.
    pub surface_index: usize,
}

/// Extract silhouette curves from a rcad_kernel::BRep for a given view direction.
///
/// This function computes the visible contour lines (silhouettes) of curved surfaces
/// as seen from a specific viewing direction. For analytic surfaces (cylinder, sphere,
/// cone, torus), exact silhouette curves are computed. For general surfaces (BSpline,
/// Bezier, etc.), numerical methods with adaptive sampling are used.
///
/// # Arguments
/// * `brep` - The rcad_kernel::BRep model to extract silhouettes from.
/// * `view_dir` - The normalized view direction (from target to eye).
/// * `opts` - Configuration options for sampling and tolerance.
///
/// # Returns
/// A vector of 3D silhouette curves, each represented as a series of world-space points.
pub fn extract_silhouette_curves(
    brep: &rcad_kernel::BRep,
    view_dir: DVec3,
    opts: &HlrOptions,
) -> Vec<SilhouetteCurve3> {
    let mut curves: Vec<SilhouetteCurve3> = Vec::new();

    if brep.solids().is_empty() {
        return curves;
    }

    let line_samples = opts.silhouette_samples.max(16);
    let dense_curve_samples = (opts.silhouette_samples * 4).max(64);

    use rcad_kernel::topods::TShape;
    let mut face_idx = 0usize;
    for ts in &brep.tshapes {
        if let TShape::Face(fd) = ts.as_ref() {
            let (surface, domain) = match (&fd.surface, fd.uv_domain) {
                (Some(surf), Some(dom)) => (surf.clone(), dom),
                (Some(surf), None) => (surf.clone(), surf.default_domain()),
                (None, _) => {
                    face_idx += 1;
                    continue;
                }
            };
            let [_u0, _u1, _v0, _v1] = domain;

            // Extract silhouettes based on surface type
            let face_curves = extract_surface_silhouettes(
                &surface,
                view_dir,
                domain,
                brep,
                opts,
                line_samples,
                dense_curve_samples,
            );

            for pts in face_curves {
                if pts.len() >= 2 {
                    curves.push(SilhouetteCurve3 {
                        points: pts,
                        surface_index: face_idx,
                    });
                }
            }

            face_idx += 1;
        }
    }

    curves
}

/// Extract silhouette curves from a single surface.
fn extract_surface_silhouettes(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    brep: &rcad_kernel::BRep,
    opts: &HlrOptions,
    line_samples: usize,
    dense_curve_samples: usize,
) -> Vec<Vec<DVec3>> {
    let [_u0, _u1, v0, v1] = domain;
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    match surface {
        Surface3::Cylinder(cyl) => {
            curves.extend(extract_cylinder_silhouettes(
                cyl,
                view_dir,
                brep,
                line_samples,
                v0,
                v1,
            ));
        }

        Surface3::Sphere(sph) => {
            curves.push(extract_sphere_silhouette(
                sph,
                view_dir,
                dense_curve_samples,
            ));
        }

        Surface3::Cone(con) => {
            curves.extend(extract_cone_silhouettes(
                con,
                view_dir,
                brep,
                line_samples,
                v0,
                v1,
            ));
        }

        Surface3::Torus(tor) => {
            curves.extend(extract_torus_silhouettes(
                tor,
                view_dir,
                dense_curve_samples,
            ));
        }

        Surface3::Ellipsoid(ell) => {
            curves.extend(extract_ellipsoid_silhouettes(
                ell,
                view_dir,
                opts,
                dense_curve_samples,
            ));
        }

        // For general surfaces, use numerical silhouette extraction
        Surface3::BSpline(_)
        | Surface3::Bezier(_)
        | Surface3::TriBezier(_)
        | Surface3::Offset(_)
        | Surface3::LinearExtrusion(_)
        | Surface3::Revolution(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::Pipe(_) => {
            curves.extend(extract_numerical_silhouettes(
                surface, view_dir, domain, opts, brep,
            ));
        }

        // Planes have no silhouette curves
        Surface3::Plane(_) | Surface3::Trimmed(_) | Surface3::Helicoid(_) => {}
    }

    curves
}

/// Extract silhouette lines from a cylinder.
fn extract_cylinder_silhouettes(
    cyl: &rcad_kernel::geom::CylindricalSurface,
    view_dir: DVec3,
    brep: &rcad_kernel::BRep,
    line_samples: usize,
    v0: f64,
    v1: f64,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    // Project view direction onto the plane perpendicular to the axis.
    let d_perp = view_dir - view_dir.dot(cyl.axis) * cyl.axis;
    if d_perp.length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT {
        // Viewing along the axis — no silhouette lines.
        return curves;
    }

    // Direction from axis to silhouette (perpendicular to both axis and d_perp).
    let sil_dir = cyl.axis.cross(d_perp).normalize_or_zero();

    // Resolve v range (height along axis).
    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
        (v0, v1)
    } else {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for vert in &brep.vertices() {
            let proj = (vert.point - cyl.origin).dot(cyl.axis);
            lo = lo.min(proj);
            hi = hi.max(proj);
        }
        if lo.is_finite() && hi.is_finite() {
            (lo, hi)
        } else {
            return curves;
        }
    };

    for &sign in &[1.0_f64, -1.0] {
        let offset = sil_dir * sign * cyl.radius;
        let world_pts: Vec<DVec3> = (0..line_samples)
            .map(|i| {
                let t = i as f64 / (line_samples - 1) as f64;
                let v = v0_eff + (v1_eff - v0_eff) * t;
                cyl.origin + v * cyl.axis + offset
            })
            .collect();
        curves.push(world_pts);
    }

    curves
}

/// Extract silhouette curve from a sphere (great circle perpendicular to view direction).
fn extract_sphere_silhouette(
    sph: &rcad_kernel::geom::SphericalSurface,
    view_dir: DVec3,
    samples: usize,
) -> Vec<DVec3> {
    let x_ax = any_perpendicular(view_dir);
    let y_ax = view_dir.cross(x_ax).normalize_or_zero();

    (0..samples)
        .map(|i| {
            let t = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
            sph.center + sph.radius * (t.cos() * x_ax + t.sin() * y_ax)
        })
        .collect()
}

/// Extract silhouette lines from a cone (two generators from apex).
fn extract_cone_silhouettes(
    con: &rcad_kernel::geom::ConicalSurface,
    view_dir: DVec3,
    brep: &rcad_kernel::BRep,
    line_samples: usize,
    v0: f64,
    v1: f64,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    let d_perp = view_dir - view_dir.dot(con.axis) * con.axis;
    if d_perp.length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT {
        return curves;
    }

    let sil_dir = con.axis.cross(d_perp).normalize_or_zero();

    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
        (v0, v1)
    } else {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for vert in &brep.vertices() {
            let proj = (vert.point - con.apex).dot(con.axis);
            lo = lo.min(proj);
            hi = hi.max(proj);
        }
        if lo.is_finite() && hi.is_finite() {
            (lo.max(0.0), hi.max(0.0))
        } else {
            return curves;
        }
    };

    let tan_a = con.half_angle_rad.tan();
    for &sign in &[1.0_f64, -1.0] {
        let world_pts: Vec<DVec3> = (0..line_samples)
            .map(|i| {
                let t = i as f64 / (line_samples - 1) as f64;
                let v = v0_eff + (v1_eff - v0_eff) * t;
                con.apex + v * con.axis + v * tan_a * sil_dir * sign
            })
            .collect();

        if world_pts
            .first()
            .zip(world_pts.last())
            .map(|(a, b)| (*b - *a).length_squared() > TOLERANCE_LEN_MIN)
            .unwrap_or(false)
        {
            curves.push(world_pts);
        }
    }

    curves
}

/// Extract silhouette curves from a torus.
fn extract_torus_silhouettes(
    tor: &rcad_kernel::geom::ToroidalSurface,
    view_dir: DVec3,
    samples: usize,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    let x_ax = any_perpendicular(tor.axis);
    let y_ax = tor.axis.cross(x_ax).normalize_or_zero();
    let axis_dot = tor.axis.dot(view_dir);

    for &offset in &[0.0_f64, std::f64::consts::PI] {
        let pts: Vec<DVec3> = (0..samples)
            .map(|i| {
                let u = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
                let radial = u.cos() * x_ax + u.sin() * y_ax;
                let radial_dot = radial.dot(view_dir);
                let v = (-radial_dot).atan2(axis_dot) + offset;
                let tube_center = tor.center + tor.major_radius * radial;
                tube_center + tor.minor_radius * (v.cos() * radial + v.sin() * tor.axis)
            })
            .collect();
        curves.push(pts);
    }

    curves
}

/// Extract silhouette curves from an ellipsoid using analytic methods.
///
/// The silhouette of an ellipsoid is the intersection of the ellipsoid surface
/// with a plane passing through the center. The plane normal is proportional to
/// (vx/a², vy/b², vz/c²) where v is the view direction and a, b, c are the radii.
///
/// This intersection is an ellipse, which we parameterize by:
/// 1. Finding two orthogonal directions in the silhouette plane
/// 2. For each angle, computing where a ray in that direction intersects the ellipsoid
fn extract_ellipsoid_silhouettes(
    ell: &rcad_kernel::geom::EllipsoidalSurface,
    view_dir: DVec3,
    opts: &HlrOptions,
    samples: usize,
) -> Vec<Vec<DVec3>> {
    use std::f64::consts::PI;

    // Build the orthonormal frame of the ellipsoid
    let (axis, x_axis, y_axis) = orthonormal_frame(ell.axis, ell.ref_dir);

    // Transform view direction into the ellipsoid's local coordinate frame
    // Local coordinates: x along x_axis, y along y_axis, z along axis
    let vx = view_dir.dot(x_axis);
    let vy = view_dir.dot(y_axis);
    let vz = view_dir.dot(axis);
    let view_local = DVec3::new(vx, vy, vz);

    // Handle degenerate case: view direction is zero
    if view_local.length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
        return Vec::new();
    }

    // The silhouette plane normal in local coordinates is proportional to
    // (vx/a², vy/b², vz/c²)
    let a = ell.radius_x;
    let b = ell.radius_y;
    let c = ell.radius_z;

    let plane_normal_local = DVec3::new(vx / (a * a), vy / (b * b), vz / (c * c));

    // Handle the degenerate case where the plane normal is zero
    // This happens when all components are zero, which shouldn't occur for valid view direction
    let plane_normal_len = plane_normal_local.length();
    if plane_normal_len < TOLERANCE_METRIC_SQ_NEAR_ZERO {
        // View direction is exactly perpendicular to all scaled axes
        // This is a degenerate case - return empty
        return Vec::new();
    }
    let plane_normal_local = plane_normal_local.normalize();

    // Check if view is along a principal axis (plane normal is near a coordinate axis)
    // In this case, the silhouette is an ellipse in the perpendicular plane
    let is_view_along_axis = plane_normal_local.z.abs() > 1.0 - TOLERANCE_MESH_LEGACY;
    let is_view_along_x = plane_normal_local.x.abs() > 1.0 - TOLERANCE_MESH_LEGACY;
    let is_view_along_y = plane_normal_local.y.abs() > 1.0 - TOLERANCE_MESH_LEGACY;

    // Find two orthogonal directions in the silhouette plane
    // These will be used to parameterize the ellipse
    let (u_dir, v_dir) = if is_view_along_axis {
        // View is along Z axis (ellipsoid's axis)
        // Silhouette plane is XY plane, silhouette is ellipse x²/a² + y²/b² = 1
        (DVec3::X, DVec3::Y)
    } else if is_view_along_x {
        // View is along X axis
        // Silhouette plane is YZ plane
        (DVec3::Y, DVec3::Z)
    } else if is_view_along_y {
        // View is along Y axis
        // Silhouette plane is XZ plane
        (DVec3::X, DVec3::Z)
    } else {
        // General case: find two orthogonal vectors in the silhouette plane
        // Use any_perpendicular to get a vector perpendicular to the plane normal
        let u = any_perpendicular(plane_normal_local);
        let v = plane_normal_local.cross(u).normalize_or_zero();
        (u, v)
    };

    // Parameterize the ellipse by sampling angles
    // For each angle θ, the ray direction is u*cos(θ) + v*sin(θ)
    // The intersection parameter t is: t = 1 / sqrt((dx/a)² + (dy/b)² + (dz/c)²)
    let actual_samples = samples.max(opts.silhouette_samples).max(32);
    let points: Vec<DVec3> = (0..actual_samples)
        .map(|i| {
            let theta = 2.0 * PI * i as f64 / actual_samples as f64;
            let cos_t = theta.cos();
            let sin_t = theta.sin();

            // Ray direction in local coordinates
            let dir_local = u_dir * cos_t + v_dir * sin_t;

            // Compute intersection parameter t
            // The ray is: p = t * dir_local
            // On ellipsoid: (tx/a)² + (ty/b)² + (tz/c)² = 1
            // t² * ((dx/a)² + (dy/b)² + (dz/c)²) = 1
            let dx = dir_local.x;
            let dy = dir_local.y;
            let dz = dir_local.z;
            let t_squared_recip = (dx * dx) / (a * a) + (dy * dy) / (b * b) + (dz * dz) / (c * c);

            let t = 1.0 / t_squared_recip.sqrt();

            // Local point on ellipsoid
            let local_point = dir_local * t;

            // Transform back to world coordinates
            ell.center + local_point.x * x_axis + local_point.y * y_axis + local_point.z * axis
        })
        .collect();

    // Verify that the silhouette points satisfy the silhouette condition
    // (normal · view_dir ≈ 0)
    let mut valid = true;
    for pt in &points {
        // Compute the point in local coordinates
        let p_local = *pt - ell.center;
        let x = p_local.dot(x_axis);
        let y = p_local.dot(y_axis);
        let z = p_local.dot(axis);

        // Normal direction (gradient of implicit equation)
        let grad_local = DVec3::new(x / (a * a), y / (b * b), z / (c * c));
        let normal = (grad_local.x * x_axis + grad_local.y * y_axis + grad_local.z * axis)
            .normalize_or_zero();
        let dot = normal.dot(view_dir);

        // Check if silhouette condition is approximately satisfied
        if dot.abs() > opts.tangent_tolerance.max(0.1) {
            valid = false;
            break;
        }
    }

    if valid && !points.is_empty() {
        vec![points]
    } else {
        // Fallback to numerical method if analytic result seems invalid
        extract_ellipsoid_silhouettes_numerical(ell, view_dir, opts, samples)
    }
}

/// Fallback numerical silhouette extraction for ellipsoids.
///
/// Used when the analytic method produces questionable results.
fn extract_ellipsoid_silhouettes_numerical(
    ell: &rcad_kernel::geom::EllipsoidalSurface,
    view_dir: DVec3,
    opts: &HlrOptions,
    samples: usize,
) -> Vec<Vec<DVec3>> {
    use rcad_kernel::geom::SurfaceEval;

    let domain = ell.default_domain(); // [0, 2π, 0, π]
    let grid_size = (samples / 4).max(16);

    // Find silhouette seed points on a grid
    let mut silhouette_points: Vec<DVec3> = Vec::new();

    let [u0, u1, v0, v1] = domain;
    for i in 0..grid_size {
        for j in 0..grid_size {
            let u = u0 + (u1 - u0) * i as f64 / (grid_size - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (grid_size - 1) as f64;

            let normal = ell.normal_at(u, v);
            let dot = normal.dot(view_dir);

            // Check if this is a silhouette point
            if dot.abs() < opts.tangent_tolerance.max(0.01) {
                let point = ell.point_at(u, v);
                silhouette_points.push(point);
            }
        }
    }

    // Sort points by angle around the silhouette center (approximation)
    if silhouette_points.len() >= 3 {
        // Compute centroid
        let _centroid = silhouette_points.iter().sum::<DVec3>() / silhouette_points.len() as f64;

        // Build the orthonormal frame
        let (_axis, x_axis, y_axis) = orthonormal_frame(ell.axis, ell.ref_dir);

        // Sort by angle in the local XY plane (projected)
        silhouette_points.sort_by(|a, b| {
            let a_local = *a - ell.center;
            let b_local = *b - ell.center;
            let ax = a_local.dot(x_axis);
            let ay = a_local.dot(y_axis);
            let bx = b_local.dot(x_axis);
            let by = b_local.dot(y_axis);
            let angle_a = ay.atan2(ax);
            let angle_b = by.atan2(bx);
            angle_a
                .partial_cmp(&angle_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        vec![silhouette_points]
    } else {
        Vec::new()
    }
}
include!("e1.rs");
