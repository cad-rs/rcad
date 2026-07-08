//! Point cloud analysis tools, analogous to OCCT 8.0 PointSetLib.
//!
//! Provides:
//! - Principal Component Analysis (PCA)
//! - Inertia tensor computation
//! - Dimensionality estimation
//! - Outlier detection and point cloud simplification
//! - Normal estimation
//! - Shape fitting (plane, sphere, cylinder)
//! - ICP registration (point-to-point, point-to-plane)
//! - Segmentation (region growing, Euclidean clustering, shape segmentation)
//! - Surface reconstruction (Poisson, Ball pivoting, Delaunay)
//! - Advanced sampling (curvature-aware, Poisson disk)
//! - BRep integration

use crate::tolerance::*;
use glam::DVec3;
use std::cmp::Ordering;

/// A collection of 3D points.
#[derive(Debug, Clone, Default)]
pub struct PointCloud {
    pub points: Vec<DVec3>,
}

impl PointCloud {
    /// Creates an empty point cloud.
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }

    /// Creates a point cloud from a slice of points.
    pub fn from_points(points: &[DVec3]) -> Self {
        Self { points: points.to_vec() }
    }

    /// Creates a point cloud from a vector of points.
    pub fn from_vec(points: Vec<DVec3>) -> Self {
        Self { points }
    }

    /// Returns the number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Returns true if the point cloud is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Computes the axis-aligned bounding box.
    pub fn bounding_box(&self) -> Option<(DVec3, DVec3)> {
        if self.points.is_empty() {
            return None;
        }
        let mut min = DVec3::splat(f64::INFINITY);
        let mut max = DVec3::splat(f64::NEG_INFINITY);
        for &p in &self.points {
            min = min.min(p);
            max = max.max(p);
        }
        Some((min, max))
    }

    /// Computes the centroid (mean) of all points.
    pub fn centroid(&self) -> Option<DVec3> {
        if self.points.is_empty() {
            return None;
        }
        let sum: DVec3 = self.points.iter().sum();
        Some(sum / self.points.len() as f64)
    }
}

/// Classification of point cloud dimensionality based on PCA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dimensionality {
    /// All points are at or very near a single location.
    Point,
    /// Points lie approximately along a line.
    Linear,
    /// Points lie approximately on a plane.
    Planar,
    /// Points have significant extent in all three dimensions.
    Volumetric,
}

/// Result of point cloud analysis.
#[derive(Debug, Clone)]
pub struct PointCloudAnalysis {
    /// Centroid (mean) of all points.
    pub centroid: DVec3,
    /// Principal axes sorted by eigenvalue (largest to smallest).
    /// - `principal_axes[0]`: direction of maximum variance
    /// - `principal_axes[2]`: direction of minimum variance (normal for planar data)
    pub principal_axes: [DVec3; 3],
    /// Principal values (eigenvalues) corresponding to each axis.
    pub principal_values: [f64; 3],
    /// Axis-aligned bounding box as (min, max).
    pub bounding_box: (DVec3, DVec3),
    /// Inertia tensor about the centroid.
    pub inertia_tensor: [[f64; 3]; 3],
    /// Estimated dimensionality.
    pub dimensionality: Dimensionality,
}

/// Performs comprehensive analysis on a point cloud.
///
/// Computes centroid, PCA, bounding box, inertia tensor, and dimensionality.
pub fn analyze_point_cloud(points: &[DVec3]) -> Option<PointCloudAnalysis> {
    if points.is_empty() {
        return None;
    }

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;

    let (principal_axes, principal_values) = compute_pca(points);

    let bounding_box = {
        let mut min = DVec3::splat(f64::INFINITY);
        let mut max = DVec3::splat(f64::NEG_INFINITY);
        for &p in points {
            min = min.min(p);
            max = max.max(p);
        }
        (min, max)
    };

    let inertia_tensor = compute_inertia_centroid(points, centroid);

    let dimensionality = estimate_dimensionality(principal_values, 0.01);

    Some(PointCloudAnalysis {
        centroid,
        principal_axes,
        principal_values,
        bounding_box,
        inertia_tensor,
        dimensionality,
    })
}

/// Computes Principal Component Analysis (PCA) on a point set.
///
/// Returns:
/// - Principal axes (eigenvectors) sorted by eigenvalue (largest first)
/// - Principal values (eigenvalues) sorted largest first
///
/// The principal axes form an orthonormal basis. For planar data,
/// `principal_axes[2]` is the normal of the best-fit plane.
pub fn compute_pca(points: &[DVec3]) -> ([DVec3; 3], [f64; 3]) {
    if points.is_empty() {
        return ([DVec3::X, DVec3::Y, DVec3::Z], [0.0; 3]);
    }

    let n = points.len() as f64;
    let centroid = points.iter().sum::<DVec3>() / n;

    // Compute covariance matrix
    let mut cov = [[0.0; 3]; 3];
    for &p in points {
        let d = p - centroid;
        cov[0][0] += d.x * d.x;
        cov[0][1] += d.x * d.y;
        cov[0][2] += d.x * d.z;
        cov[1][1] += d.y * d.y;
        cov[1][2] += d.y * d.z;
        cov[2][2] += d.z * d.z;
    }
    cov[0][0] /= n;
    cov[0][1] /= n;
    cov[0][2] /= n;
    cov[1][1] /= n;
    cov[1][2] /= n;
    cov[2][2] /= n;
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    // Compute eigenvalues and eigenvectors using power iteration
    let (eigenvalues, eigenvectors) = compute_eigendecomposition_3x3(&cov);

    // Sort by eigenvalue descending
    let mut indexed: [(usize, f64, DVec3); 3] = [
        (0, eigenvalues[0], eigenvectors[0]),
        (1, eigenvalues[1], eigenvectors[1]),
        (2, eigenvalues[2], eigenvectors[2]),
    ];
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));

    let mut axes = [DVec3::ZERO; 3];
    let mut values = [0.0; 3];
    for i in 0..3 {
        axes[i] = indexed[i].2;
        values[i] = indexed[i].1.max(0.0);
    }

    // Ensure orthonormal right-handed basis
    axes[2] = axes[0].cross(axes[1]).normalize_or(DVec3::Z);
    axes[1] = axes[2].cross(axes[0]).normalize_or(DVec3::Y);

    (axes, values)
}

/// Compute eigenvalues and eigenvectors of a symmetric 3x3 matrix.
fn compute_eigendecomposition_3x3(m: &[[f64; 3]; 3]) -> ([f64; 3], [DVec3; 3]) {
    // Use Jacobi eigenvalue algorithm for symmetric matrices
    let mut a = *m;
    let mut v = [
        DVec3::X,
        DVec3::Y,
        DVec3::Z,
    ];

    const MAX_ITERATIONS: usize = 100;
    const TOLERANCE: f64 = TOLERANCE_LEN_MIN;

    for _ in 0..MAX_ITERATIONS {
        // Find the largest off-diagonal element
        let mut max_val = 0.0;
        let (mut p, mut q) = (0, 1);

        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < TOLERANCE {
            break;
        }

        // Compute rotation angle
        let theta = if (a[p][p] - a[q][q]).abs() < TOLERANCE {
            std::f64::consts::FRAC_PI_4 * a[p][q].signum()
        } else {
            0.5 * (2.0 * a[p][q] / (a[p][p] - a[q][q])).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Update matrix A
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];

        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for i in 0..3 {
            if i != p && i != q {
                let aip = a[i][p];
                let aiq = a[i][q];
                a[i][p] = c * aip - s * aiq;
                a[p][i] = a[i][p];
                a[i][q] = s * aip + c * aiq;
                a[q][i] = a[i][q];
            }
        }

        // Update eigenvector matrix
        for i in 0..3 {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
        }
    }

    // Normalize eigenvectors
    for i in 0..3 {
        v[i] = v[i].normalize_or(match i {
            0 => DVec3::X,
            1 => DVec3::Y,
            _ => DVec3::Z,
        });
    }

    ([a[0][0], a[1][1], a[2][2]], v)
}

/// Computes the inertia tensor of a point set about the origin.
///
/// The inertia tensor is a symmetric 3x3 matrix defined as:
/// ```text
/// Ixx = Σ(y²+z²),  Iyy = Σ(x²+z²),  Izz = Σ(x²+y²)
/// Ixy = -Σxy,       Ixz = -Σxz,       Iyz = -Σyz
/// ```
pub fn compute_inertia(points: &[DVec3]) -> [[f64; 3]; 3] {
    if points.is_empty() {
        return [[0.0; 3]; 3];
    }

    compute_inertia_centroid(points, DVec3::ZERO)
}

/// Computes the inertia tensor of a point set about a given centroid.
fn compute_inertia_centroid(points: &[DVec3], centroid: DVec3) -> [[f64; 3]; 3] {
    let mut ixx = 0.0_f64;
    let mut iyy = 0.0_f64;
    let mut izz = 0.0_f64;
    let mut ixy = 0.0_f64;
    let mut ixz = 0.0_f64;
    let mut iyz = 0.0_f64;

    for &p in points {
        let d = p - centroid;
        let x = d.x;
        let y = d.y;
        let z = d.z;

        ixx += y * y + z * z;
        iyy += x * x + z * z;
        izz += x * x + y * y;
        ixy -= x * y;
        ixz -= x * z;
        iyz -= y * z;
    }

    [
        [ixx, ixy, ixz],
        [ixy, iyy, iyz],
        [ixz, iyz, izz],
    ]
}

/// Estimates the dimensionality of a point cloud from PCA eigenvalues.
///
/// The threshold is the relative tolerance for considering an eigenvalue
/// as "negligible" compared to the largest eigenvalue.
///
/// Classification:
/// - Point: all eigenvalues are negligible (total variance near zero)
/// - Linear: only one significant eigenvalue
/// - Planar: two significant eigenvalues, one negligible
/// - Volumetric: all three eigenvalues are significant
pub fn estimate_dimensionality(pca_values: [f64; 3], threshold: f64) -> Dimensionality {
    // Sort eigenvalues descending
    let mut sorted = pca_values;
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));

    let total: f64 = sorted.iter().sum();
    // If total variance is negligible, it's a point
    if total < TOLERANCE_LINEAR_ULTRA_STRICT {
        return Dimensionality::Point;
    }

    // Normalize by largest eigenvalue
    let max_val = sorted[0].max(TOLERANCE_METRIC_SQ_NEAR_ZERO);
    let rel1 = sorted[1] / max_val;
    let rel2 = sorted[2] / max_val;

    // Count how many eigenvalues are significant (relative to max)
    // First eigenvalue is always significant if total > 0
    let sig1 = true;
    let sig2 = rel1 > threshold;
    let sig3 = rel2 > threshold;

    let count = [sig1, sig2, sig3].iter().filter(|&&x| x).count();

    match count {
        0 => Dimensionality::Point,
        1 => Dimensionality::Linear,
        2 => Dimensionality::Planar,
        3 => Dimensionality::Volumetric,
        _ => Dimensionality::Volumetric,
    }
}

// ============================================================================
// Point Cloud Processing
// ============================================================================

/// Detected outlier point with its outlier score.
#[derive(Debug, Clone)]
pub struct OutlierPoint {
    /// Index of the outlier point in the original point cloud.
    pub index: usize,
    /// Outlier score (higher = more likely an outlier).
    pub score: f64,
}

/// Detects outlier points using the Local Outlier Factor (LOF) algorithm.
///
/// Parameters:
/// - `points`: the point cloud
/// - `k`: number of nearest neighbors to consider (default: 20)
/// - `threshold`: LOF score threshold for outliers (default: 2.0)
///
/// Returns a list of outlier points sorted by score (highest first).
pub fn detect_outliers(points: &[DVec3], k: usize, threshold: f64) -> Vec<OutlierPoint> {
    if points.len() <= k + 1 {
        return Vec::new();
    }

    let k = k.min(points.len() - 1).max(1);
    let n = points.len();

    // Compute k-distances and reachability distances
    let mut lof_scores = vec![0.0; n];

    for i in 0..n {
        // Find k nearest neighbors
        let mut distances: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, (points[j] - points[i]).length_squared()))
            .collect();
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

        let k_dist = distances[k - 1].1.sqrt();
        let neighbors: Vec<usize> = distances.iter().take(k).map(|&(j, _)| j).collect();

        // Compute local reachability density
        let mut lrd_sum = 0.0;
        for &j in &neighbors {
            let dist_ij = (points[j] - points[i]).length();
            let reach_dist = dist_ij.max(k_dist);
            lrd_sum += reach_dist;
        }
        let lrd_i = k as f64 / lrd_sum;

        // Compute LOF
        let mut lof_sum = 0.0;
        for &j in &neighbors {
            // Compute lrd of neighbor j (simplified - use same k_dist approximation)
            let mut j_dists: Vec<f64> = (0..n)
                .filter(|&l| l != j)
                .map(|l| (points[l] - points[j]).length())
                .collect();
            j_dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
            let j_k_dist = j_dists[k - 1];
            let mut j_lrd_sum = 0.0;
            for &l in j_dists.iter().take(k) {
                j_lrd_sum += l.max(j_k_dist);
            }
            let lrd_j = k as f64 / j_lrd_sum;
            lof_sum += lrd_j / lrd_i;
        }
        lof_scores[i] = lof_sum / k as f64;
    }

    // Collect outliers
    let mut outliers: Vec<OutlierPoint> = lof_scores
        .iter()
        .enumerate()
        .filter(|&(_, &score)| score > threshold)
        .map(|(i, &score)| OutlierPoint { index: i, score })
        .collect();

    outliers.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal));
    outliers
}

/// Removes outliers from a point cloud based on LOF detection.
///
/// Returns a new point cloud with outliers removed.
pub fn remove_outliers(points: &[DVec3], k: usize, threshold: f64) -> Vec<DVec3> {
    let outliers = detect_outliers(points, k, threshold);
    let outlier_set: std::collections::HashSet<usize> = outliers.iter().map(|o| o.index).collect();

    points
        .iter()
        .enumerate()
        .filter(|(i, _)| !outlier_set.contains(i))
        .map(|(_, &p)| p)
        .collect()
}

/// Sampling strategy for point cloud simplification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SamplingStrategy {
    /// Uniform random sampling.
    Random,
    /// Grid-based voxel sampling (one point per voxel).
    Voxel,
    /// Farthest point sampling (maximizes coverage).
    FarthestPoint,
}

/// Simplifies a point cloud by reducing the number of points.
///
/// Parameters:
/// - `points`: the input point cloud
/// - `target_count`: target number of output points
/// - `strategy`: sampling strategy to use
///
/// Returns the simplified point cloud.
pub fn simplify_point_cloud(
    points: &[DVec3],
    target_count: usize,
    strategy: SamplingStrategy,
) -> Vec<DVec3> {
    if points.len() <= target_count {
        return points.to_vec();
    }

    match strategy {
        SamplingStrategy::Random => random_sample(points, target_count),
        SamplingStrategy::Voxel => voxel_sample(points, target_count),
        SamplingStrategy::FarthestPoint => farthest_point_sample(points, target_count),
    }
}

fn random_sample(points: &[DVec3], target_count: usize) -> Vec<DVec3> {
    use std::collections::HashSet;
    let n = points.len();
    let mut indices: HashSet<usize> = HashSet::new();
    let mut rng = SimpleRng::new(12345);

    while indices.len() < target_count {
        indices.insert((rng.next() as usize) % n);
    }

    indices.iter().map(|&i| points[i]).collect()
}

fn voxel_sample(points: &[DVec3], target_count: usize) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    // Compute bounding box
    let (min, max) = points.iter().fold(
        (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
        |(min, max), &p| (min.min(p), max.max(p)),
    );

    // Estimate voxel size
    let volume = (max.x - min.x) * (max.y - min.y) * (max.z - min.z);
    let voxel_size = (volume / target_count as f64).cbrt().max(TOLERANCE_LINEAR_ULTRA_STRICT);

    // Group points by voxel
    let mut voxels: std::collections::HashMap<[i64; 3], Vec<DVec3>> = std::collections::HashMap::new();

    for &p in points {
        let key = [
            ((p.x - min.x) / voxel_size).floor() as i64,
            ((p.y - min.y) / voxel_size).floor() as i64,
            ((p.z - min.z) / voxel_size).floor() as i64,
        ];
        voxels.entry(key).or_default().push(p);
    }

    // Take centroid of each voxel
    voxels
        .values()
        .map(|pts| {
            let sum: DVec3 = pts.iter().sum();
            sum / pts.len() as f64
        })
        .collect()
}

fn farthest_point_sample(points: &[DVec3], target_count: usize) -> Vec<DVec3> {
    if points.len() <= target_count {
        return points.to_vec();
    }

    let n = points.len();
    let mut selected = Vec::with_capacity(target_count);
    let mut distances = vec![f64::INFINITY; n];

    // Start with centroid or first point
    let centroid = points.iter().sum::<DVec3>() / n as f64;
    let first_idx = points
        .iter()
        .enumerate()
        .map(|(i, &p)| (i, (p - centroid).length_squared()))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    selected.push(points[first_idx]);

    // Greedy farthest point selection
    while selected.len() < target_count {
        let last = selected.last().unwrap();
        let mut farthest_idx = 0;
        let mut farthest_dist = 0.0;

        for i in 0..n {
            let d = (points[i] - *last).length_squared();
            distances[i] = distances[i].min(d);
            if distances[i] > farthest_dist {
                farthest_dist = distances[i];
                farthest_idx = i;
            }
        }

        selected.push(points[farthest_idx]);
    }

    selected
}

/// Simple deterministic RNG for reproducible sampling.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        // xorshift64
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }
}

/// Estimates normals for each point using local PCA.
///
/// For each point, computes PCA on its k nearest neighbors.
/// The normal is the eigenvector corresponding to the smallest eigenvalue.
///
/// Returns a vector of unit normals (one per point).
pub fn estimate_normals(points: &[DVec3], k: usize) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    let k = k.min(points.len() - 1).max(2);
    let n = points.len();
    let mut normals = Vec::with_capacity(n);

    for i in 0..n {
        // Find k nearest neighbors
        let mut distances: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, (points[j] - points[i]).length_squared()))
            .collect();
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

        let neighbor_pts: Vec<DVec3> = distances
            .iter()
            .take(k)
            .map(|&(j, _)| points[j])
            .collect();

        // PCA on neighbors
        let (axes, values) = compute_pca(&neighbor_pts);

        // Normal is direction of minimum variance
        // Check if the smallest eigenvalue is small enough for a planar fit
        let max_val = values[0].max(TOLERANCE_METRIC_SQ_NEAR_ZERO);
        if values[2] / max_val < 0.1 {
            normals.push(axes[2]);
        } else {
            // Not planar enough, still use smallest variance direction
            normals.push(axes[2]);
        }
    }

    normals
}

// ============================================================================
// Point Cloud Fitting
// ============================================================================

/// Result of fitting a plane to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedPlane {
    /// A point on the plane.
    pub point: DVec3,
    /// Unit normal of the plane.
    pub normal: DVec3,
    /// RMS distance of points to the fitted plane.
    pub rms_error: f64,
}

/// Fits a plane to a point cloud using least squares.
///
/// Uses PCA: the normal is the eigenvector corresponding to the smallest eigenvalue.
pub fn fit_plane(points: &[DVec3]) -> Option<FittedPlane> {
    if points.len() < 3 {
        return None;
    }

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;
    let (axes, _values) = compute_pca(points);

    // Compute RMS error
    let normal = axes[2];
    let mut sum_sq = 0.0;
    for &p in points {
        let d = (p - centroid).dot(normal);
        sum_sq += d * d;
    }
    let rms_error = (sum_sq / points.len() as f64).sqrt();

    Some(FittedPlane {
        point: centroid,
        normal,
        rms_error,
    })
}

/// Result of fitting a sphere to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedSphere {
    /// Center of the sphere.
    pub center: DVec3,
    /// Radius of the sphere.
    pub radius: f64,
    /// RMS distance of points to the fitted sphere surface.
    pub rms_error: f64,
}

/// Fits a sphere to a point cloud using least squares.
///
/// Uses an algebraic fit followed by geometric refinement.
pub fn fit_sphere(points: &[DVec3]) -> Option<FittedSphere> {
    if points.len() < 4 {
        return None;
    }

    // Algebraic fit using linear least squares
    // Fit: (x - cx)^2 + (y - cy)^2 + (z - cz)^2 = r^2
    // Rewrite: x^2 + y^2 + z^2 - 2*cx*x - 2*cy*y - 2*cz*z + cx^2 + cy^2 + cz^2 - r^2 = 0
    // Let: a = -2*cx, b = -2*cy, c = -2*cz, d = cx^2 + cy^2 + cz^2 - r^2
    // Then: x^2 + y^2 + z^2 + a*x + b*y + c*z + d = 0
    // Solve for a, b, c, d using linear least squares

    let n = points.len();
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_y2 = 0.0;
    let mut sum_z2 = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_xz = 0.0;
    let mut sum_yz = 0.0;
    let mut sum_xyz = 0.0;  // x*(x^2+y^2+z^2)
    let mut sum_x2yz = 0.0; // y*(x^2+y^2+z^2)
    let mut sum_xz2 = 0.0;  // z*(x^2+y^2+z^2)
    let mut sum_r2 = 0.0;   // x^2 + y^2 + z^2

    for &p in points {
        let x = p.x;
        let y = p.y;
        let z = p.z;
        let r2 = x * x + y * y + z * z;

        sum_x += x;
        sum_y += y;
        sum_z += z;
        sum_x2 += x * x;
        sum_y2 += y * y;
        sum_z2 += z * z;
        sum_xy += x * y;
        sum_xz += x * z;
        sum_yz += y * z;
        sum_xyz += x * r2;
        sum_x2yz += y * r2;
        sum_xz2 += z * r2;
        sum_r2 += r2;
    }

    // Solve 4x4 linear system: A * [a, b, c, d]^T = B
    // Where A is the matrix of sums, B is the RHS
    let a = [
        [sum_x2, sum_xy, sum_xz, sum_x],
        [sum_xy, sum_y2, sum_yz, sum_y],
        [sum_xz, sum_yz, sum_z2, sum_z],
        [sum_x, sum_y, sum_z, n as f64],
    ];
    let b = [-sum_xyz, -sum_x2yz, -sum_xz2, -sum_r2];

    // Solve using Gaussian elimination
    let coeffs = solve_linear_4x4(&a, &b)?;

    let cx = -coeffs[0] / 2.0;
    let cy = -coeffs[1] / 2.0;
    let cz = -coeffs[2] / 2.0;
    let center = DVec3::new(cx, cy, cz);
    let radius = (cx * cx + cy * cy + cz * cz - coeffs[3]).sqrt().max(0.0);

    // Compute RMS error
    let mut sum_sq = 0.0;
    for &p in points {
        let d = (p - center).length() - radius;
        sum_sq += d * d;
    }
    let rms_error = (sum_sq / n as f64).sqrt();

    Some(FittedSphere {
        center,
        radius,
        rms_error,
    })
}

/// Solve a 4x4 linear system using Gaussian elimination with partial pivoting.
fn solve_linear_4x4(a: &[[f64; 4]; 4], b: &[f64; 4]) -> Option<[f64; 4]> {
    const N: usize = 4;
    let mut m = *a;
    let mut v = *b;

    // Forward elimination
    for col in 0..N {
        // Find pivot
        let mut max_row = col;
        let mut max_val = m[col][col].abs();
        for row in (col + 1)..N {
            if m[row][col].abs() > max_val {
                max_val = m[row][col].abs();
                max_row = row;
            }
        }

        if max_val < TOLERANCE_FLOAT_LOOSE {
            return None; // Singular matrix
        }

        // Swap rows
        m.swap(col, max_row);
        v.swap(col, max_row);

        // Eliminate
        for row in (col + 1)..N {
            let factor = m[row][col] / m[col][col];
            for j in col..N {
                m[row][j] -= factor * m[col][j];
            }
            v[row] -= factor * v[col];
        }
    }

    // Back substitution
    let mut x = [0.0; N];
    for i in (0..N).rev() {
        let mut sum = v[i];
        for j in (i + 1)..N {
            sum -= m[i][j] * x[j];
        }
        x[i] = sum / m[i][i];
    }

    Some(x)
}

/// Result of fitting a cylinder to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedCylinder {
    /// A point on the cylinder axis.
    pub axis_point: DVec3,
    /// Unit direction of the cylinder axis.
    pub axis_direction: DVec3,
    /// Radius of the cylinder.
    pub radius: f64,
    /// RMS distance of points to the fitted cylinder surface.
    pub rms_error: f64,
}

/// Fits a cylinder to a point cloud using iterative optimization.
///
/// The algorithm:
/// 1. Estimate the axis direction using PCA on differences from centroid
/// 2. Project points onto the plane perpendicular to axis
/// 3. Fit a circle to the projected points
pub fn fit_cylinder(points: &[DVec3]) -> Option<FittedCylinder> {
    if points.len() < 5 {
        return None;
    }

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;
    let (axes, values) = compute_pca(points);

    // For a cylinder, we expect two large eigenvalues and one small
    // The cylinder axis is the direction of minimum variance
    let max_val = values[0].max(TOLERANCE_METRIC_SQ_NEAR_ZERO);
    if values[2] / max_val > 0.3 {
        // Not cylindrical enough
        // Try alternative: axis might be along maximum variance direction
        // This happens for short cylinders
    }

    // Try both possible axis directions and pick the better fit
    let axis_candidates = [axes[2], axes[0]];

    let mut best_fit: Option<FittedCylinder> = None;
    let mut best_error = f64::INFINITY;

    for axis in axis_candidates {
        if let Some(cyl) = fit_cylinder_with_axis(points, centroid, axis)
            && cyl.rms_error < best_error {
                best_error = cyl.rms_error;
                best_fit = Some(cyl);
            }
    }

    best_fit
}

fn fit_cylinder_with_axis(points: &[DVec3], centroid: DVec3, axis: DVec3) -> Option<FittedCylinder> {
    let axis = axis.normalize_or(DVec3::Z);

    // Build orthonormal basis with axis as Z
    let u = if axis.x.abs() < 0.9 {
        axis.cross(DVec3::X).normalize()
    } else {
        axis.cross(DVec3::Y).normalize()
    };
    let v = axis.cross(u);

    // Project points onto the plane perpendicular to axis
    let projected: Vec<DVec2> = points
        .iter()
        .map(|&p| {
            let d = p - centroid;
            DVec2::new(d.dot(u), d.dot(v))
        })
        .collect();

    // Fit a circle to projected points
    let circle = fit_circle_2d(&projected)?;

    // Transform back to 3D
    let center_2d = circle.center;
    let center_3d = centroid + center_2d.x * u + center_2d.y * v;

    // Compute RMS error in 3D
    let mut sum_sq = 0.0;
    for &p in points {
        let to_center = p - center_3d;
        let axial_dist = to_center.dot(axis);
        let radial = to_center - axial_dist * axis;
        let d = radial.length() - circle.radius;
        sum_sq += d * d;
    }
    let rms_error = (sum_sq / points.len() as f64).sqrt();

    Some(FittedCylinder {
        axis_point: center_3d,
        axis_direction: axis,
        radius: circle.radius,
        rms_error,
    })
}

/// 2D vector for circle fitting.
#[derive(Debug, Clone, Copy)]
struct DVec2 {
    x: f64,
    y: f64,
}

impl DVec2 {
    fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    fn length(self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}

/// Fitted circle in 2D.
struct FittedCircle {
    center: DVec2,
    radius: f64,
}

/// Fit a circle to 2D points using least squares.
fn fit_circle_2d(points: &[DVec2]) -> Option<FittedCircle> {
    if points.len() < 3 {
        return None;
    }

    // Use algebraic fit similar to sphere fitting
    let n = points.len();
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_y2 = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_xr = 0.0;
    let mut sum_yr = 0.0;
    let mut sum_r = 0.0;

    for &p in points {
        let x = p.x;
        let y = p.y;
        let r = x * x + y * y;

        sum_x += x;
        sum_y += y;
        sum_x2 += x * x;
        sum_y2 += y * y;
        sum_xy += x * y;
        sum_xr += x * r;
        sum_yr += y * r;
        sum_r += r;
    }

    // Solve 3x3: [sum_x2, sum_xy, sum_x] [a]   [-sum_xr]
    //            [sum_xy, sum_y2, sum_y] [b] = [-sum_yr]
    //            [sum_x,  sum_y,  n   ] [d]   [-sum_r ]
    let a = [
        [sum_x2, sum_xy, sum_x],
        [sum_xy, sum_y2, sum_y],
        [sum_x, sum_y, n as f64],
    ];
    let b = [-sum_xr, -sum_yr, -sum_r];

    let coeffs = solve_linear_3x3(&a, &b)?;

    let cx = -coeffs[0] / 2.0;
    let cy = -coeffs[1] / 2.0;
    let radius = (cx * cx + cy * cy - coeffs[2]).sqrt().max(0.0);

    Some(FittedCircle {
        center: DVec2::new(cx, cy),
        radius,
    })
}

/// Solve a 3x3 linear system using Gaussian elimination.
fn solve_linear_3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    const N: usize = 3;
    let mut m = *a;
    let mut v = *b;

    // Forward elimination
    for col in 0..N {
        let mut max_row = col;
        let mut max_val = m[col][col].abs();
        for row in (col + 1)..N {
            if m[row][col].abs() > max_val {
                max_val = m[row][col].abs();
                max_row = row;
            }
        }

        if max_val < TOLERANCE_FLOAT_LOOSE {
            return None;
        }

        m.swap(col, max_row);
        v.swap(col, max_row);

        for row in (col + 1)..N {
            let factor = m[row][col] / m[col][col];
            for j in col..N {
                m[row][j] -= factor * m[col][j];
            }
            v[row] -= factor * v[col];
        }
    }

    // Back substitution
    let mut x = [0.0; N];
    for i in (0..N).rev() {
        let mut sum = v[i];
        for j in (i + 1)..N {
            sum -= m[i][j] * x[j];
        }
        x[i] = sum / m[i][i];
    }

    Some(x)
}

/// Result of fitting a convex polygon to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedPolygon {
    /// Vertices of the fitted polygon.
    pub vertices: Vec<DVec3>,
    /// Plane of the polygon.
    pub plane_point: DVec3,
    pub plane_normal: DVec3,
    /// Area of the polygon.
    pub area: f64,
}

/// Fits a convex polygon to a planar point cloud.
///
/// Projects points to the best-fit plane and computes the 2D convex hull.
pub fn fit_polygon(points: &[DVec3]) -> Option<FittedPolygon> {
    if points.len() < 3 {
        return None;
    }

    let plane = fit_plane(points)?;

    // Build orthonormal basis on the plane
    let normal = plane.normal;
    let u = if normal.x.abs() < 0.9 {
        normal.cross(DVec3::X).normalize()
    } else {
        normal.cross(DVec3::Y).normalize()
    };
    let v = normal.cross(u);

    // Project to 2D
    let projected_2d: Vec<DVec2> = points
        .iter()
        .map(|&p| {
            let d = p - plane.point;
            DVec2::new(d.dot(u), d.dot(v))
        })
        .collect();

    // Compute 2D convex hull
    let hull_2d = convex_hull_2d(&projected_2d);

    if hull_2d.len() < 3 {
        return None;
    }

    // Transform back to 3D
    let nodes: Vec<DVec3> = hull_2d
        .iter()
        .map(|&p2d| plane.point + p2d.x * u + p2d.y * v)
        .collect();

    // Compute area using shoelace formula in 3D
    let mut area = 0.0;
    let n = nodes.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let cross = nodes[i].cross(nodes[j]);
        area += cross.dot(normal);
    }
    area = area.abs() / 2.0;

    Some(FittedPolygon {
        vertices: nodes,
        plane_point: plane.point,
        plane_normal: plane.normal,
        area,
    })
}

/// Compute the convex hull of 2D points using Andrew's monotone chain algorithm.
/// This is more robust than Graham scan for numerical stability.
fn convex_hull_2d(points: &[DVec2]) -> Vec<DVec2> {
    if points.len() < 3 {
        return points.to_vec();
    }

    // Sort points by x, then by y
    let mut sorted: Vec<usize> = (0..points.len()).collect();
    sorted.sort_by(|&a, &b| {
        let pa = points[a];
        let pb = points[b];
        if (pa.x - pb.x).abs() > TOLERANCE_FLOAT_LOOSE {
            pa.x.partial_cmp(&pb.x).unwrap_or(Ordering::Equal)
        } else {
            pa.y.partial_cmp(&pb.y).unwrap_or(Ordering::Equal)
        }
    });

    // Cross product of OA and OB vectors
    let cross = |o: DVec2, a: DVec2, b: DVec2| -> f64 {
        (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
    };

    // Build lower hull
    let mut lower: Vec<DVec2> = Vec::new();
    for &i in &sorted {
        while lower.len() >= 2 {
            let n = lower.len();
            if cross(lower[n - 2], lower[n - 1], points[i]) <= 0.0 {
                lower.pop();
            } else {
                break;
            }
        }
        lower.push(points[i]);
    }

    // Build upper hull
    let mut upper: Vec<DVec2> = Vec::new();
    for &i in sorted.iter().rev() {
        while upper.len() >= 2 {
            let n = upper.len();
            if cross(upper[n - 2], upper[n - 1], points[i]) <= 0.0 {
                upper.pop();
            } else {
                break;
            }
        }
        upper.push(points[i]);
    }

    // Remove last point of each half because it's repeated at the beginning of the other half
    lower.pop();
    upper.pop();

    // Concatenate lower and upper hulls
    lower.extend(upper);
    lower
}

impl std::ops::Sub for DVec2 {
    type Output = DVec2;

    fn sub(self, other: DVec2) -> DVec2 {
        DVec2::new(self.x - other.x, self.y - other.y)
    }
}

// ============================================================================
// ICP Registration
// ============================================================================

/// Result of ICP registration.
#[derive(Debug, Clone)]
pub struct IcpResult {
    /// Rotation matrix (3x3) to transform source to target.
    pub rotation: [[f64; 3]; 3],
    /// Translation vector to transform source to target.
    pub translation: DVec3,
    /// Final RMS error after convergence.
    pub rms_error: f64,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the algorithm converged within tolerance.
    pub converged: bool,
}

impl IcpResult {
    /// Applies the transformation to a point.
    pub fn transform_point(&self, point: DVec3) -> DVec3 {
        let r = &self.rotation;
        DVec3::new(
            r[0][0] * point.x + r[0][1] * point.y + r[0][2] * point.z + self.translation.x,
            r[1][0] * point.x + r[1][1] * point.y + r[1][2] * point.z + self.translation.y,
            r[2][0] * point.x + r[2][1] * point.y + r[2][2] * point.z + self.translation.z,
        )
    }

    /// Returns the transformation as a 4x4 homogeneous matrix.
    pub fn to_matrix(&self) -> [[f64; 4]; 4] {
        let r = &self.rotation;
        let t = &self.translation;
        [
            [r[0][0], r[0][1], r[0][2], t.x],
            [r[1][0], r[1][1], r[1][2], t.y],
            [r[2][0], r[2][1], r[2][2], t.z],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
}

/// ICP algorithm variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IcpVariant {
    /// Standard point-to-point ICP.
    PointToPoint,
    /// Point-to-plane ICP (requires normals on target).
    PointToPlane,
}

/// ICP configuration parameters.
#[derive(Debug, Clone)]
pub struct IcpConfig {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance for RMS error change.
    pub tolerance: f64,
    /// Maximum correspondence distance (points beyond this are ignored).
    pub max_correspondence_distance: f64,
    /// Whether to use reciprocal correspondence (both directions).
    pub use_reciprocal: bool,
}

impl Default for IcpConfig {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: TOLERANCE_MESH_LEGACY,
            max_correspondence_distance: f64::INFINITY,
            use_reciprocal: false,
        }
    }
}

/// Performs Iterative Closest Point (ICP) registration.
///
/// Aligns the source point cloud to the target point cloud.
///
/// # Arguments
/// * `source` - Source point cloud to be transformed
/// * `target` - Target point cloud (reference)
/// * `variant` - ICP variant to use
/// * `config` - Configuration parameters
///
/// # Returns
/// * `IcpResult` containing the transformation and convergence info
pub fn icp_registration(
    source: &[DVec3],
    target: &[DVec3],
    variant: IcpVariant,
    config: &IcpConfig,
) -> Option<IcpResult> {
    if source.is_empty() || target.is_empty() {
        return None;
    }

    match variant {
        IcpVariant::PointToPoint => icp_point_to_point(source, target, config),
        IcpVariant::PointToPlane => {
            // Estimate normals for target if not provided
            let target_normals = estimate_normals(target, 10);
            icp_point_to_plane(source, target, &target_normals, config)
        }
    }
}

/// Performs ICP registration with pre-computed normals.
pub fn icp_registration_with_normals(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    config: &IcpConfig,
) -> Option<IcpResult> {
    if source.is_empty() || target.is_empty() || target_normals.len() != target.len() {
        return None;
    }
    icp_point_to_plane(source, target, target_normals, config)
}

fn icp_point_to_point(
    source: &[DVec3],
    target: &[DVec3],
    config: &IcpConfig,
) -> Option<IcpResult> {
    let mut transformed: Vec<DVec3> = source.to_vec();
    let mut cumulative_rotation = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut cumulative_translation = DVec3::ZERO;

    // Build KD-tree for target (simple brute-force for now)
    let mut prev_error = f64::INFINITY;
    let mut converged = false;

    for _iteration in 0..config.max_iterations {
        // Find correspondences
        let correspondences = find_correspondences(
            &transformed,
            target,
            config.max_correspondence_distance,
        );

        if correspondences.is_empty() {
            break;
        }

        // Compute transformation using SVD
        let (rotation, translation) = compute_transformation_svd(&transformed, target, &correspondences)?;

        // Apply transformation
        for p in &mut transformed {
            *p = apply_transform(*p, &rotation, &translation);
        }

        // Update cumulative transformation
        cumulative_rotation = multiply_matrices(&rotation, &cumulative_rotation);
        cumulative_translation = apply_transform_to_vector(&cumulative_translation, &rotation, &translation);

        // Compute error
        let rms_error = compute_rms_error(&transformed, target, &correspondences);

        // Check convergence
        if (prev_error - rms_error).abs() < config.tolerance {
            converged = true;
            break;
        }
        prev_error = rms_error;
    }

    Some(IcpResult {
        rotation: cumulative_rotation,
        translation: cumulative_translation,
        rms_error: prev_error,
        iterations: config.max_iterations,
        converged,
    })
}

fn icp_point_to_plane(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    config: &IcpConfig,
) -> Option<IcpResult> {
    let mut transformed: Vec<DVec3> = source.to_vec();
    let mut cumulative_rotation = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut cumulative_translation = DVec3::ZERO;

    let mut prev_error = f64::INFINITY;
    let mut converged = false;

    for _iteration in 0..config.max_iterations {
        let correspondences = find_correspondences(
            &transformed,
            target,
            config.max_correspondence_distance,
        );

        if correspondences.is_empty() {
            break;
        }

        // Compute point-to-plane transformation using linear least squares
        let (rotation, translation) = compute_point_to_plane_transformation(
            &transformed,
            target,
            target_normals,
            &correspondences,
        )?;

        // Apply transformation
        for p in &mut transformed {
            *p = apply_transform(*p, &rotation, &translation);
        }

        // Update cumulative transformation
        cumulative_rotation = multiply_matrices(&rotation, &cumulative_rotation);
        cumulative_translation = apply_transform_to_vector(&cumulative_translation, &rotation, &translation);

        // Compute point-to-plane error
        let rms_error = compute_point_to_plane_error(&transformed, target, target_normals, &correspondences);

        if (prev_error - rms_error).abs() < config.tolerance {
            converged = true;
            break;
        }
        prev_error = rms_error;
    }

    Some(IcpResult {
        rotation: cumulative_rotation,
        translation: cumulative_translation,
        rms_error: prev_error,
        iterations: config.max_iterations,
        converged,
    })
}

/// Correspondence between source and target points.
struct Correspondence {
    source_idx: usize,
    target_idx: usize,
    distance: f64,
}

fn find_correspondences(
    source: &[DVec3],
    target: &[DVec3],
    max_distance: f64,
) -> Vec<Correspondence> {
    let mut correspondences = Vec::with_capacity(source.len());

    for (i, &s) in source.iter().enumerate() {
        let mut best_dist = f64::INFINITY;
        let mut best_j = 0;

        for (j, &t) in target.iter().enumerate() {
            let d = (s - t).length_squared();
            if d < best_dist {
                best_dist = d;
                best_j = j;
            }
        }

        let dist = best_dist.sqrt();
        if dist <= max_distance {
            correspondences.push(Correspondence {
                source_idx: i,
                target_idx: best_j,
                distance: dist,
            });
        }
    }

    correspondences
}

fn compute_transformation_svd(
    source: &[DVec3],
    target: &[DVec3],
    correspondences: &[Correspondence],
) -> Option<([[f64; 3]; 3], DVec3)> {
    if correspondences.is_empty() {
        return None;
    }

    let n = correspondences.len();

    // Compute centroids
    let mut source_centroid = DVec3::ZERO;
    let mut target_centroid = DVec3::ZERO;

    for c in correspondences {
        source_centroid += source[c.source_idx];
        target_centroid += target[c.target_idx];
    }

    source_centroid /= n as f64;
    target_centroid /= n as f64;

    // Compute cross-covariance matrix
    let mut h = [[0.0_f64; 3]; 3];
    for c in correspondences {
        let s = source[c.source_idx] - source_centroid;
        let t = target[c.target_idx] - target_centroid;
        h[0][0] += s.x * t.x;
        h[0][1] += s.x * t.y;
        h[0][2] += s.x * t.z;
        h[1][0] += s.y * t.x;
        h[1][1] += s.y * t.y;
        h[1][2] += s.y * t.z;
        h[2][0] += s.z * t.x;
        h[2][1] += s.z * t.y;
        h[2][2] += s.z * t.z;
    }

    // SVD of H
    let (u, _, v) = svd_3x3(&h)?;

    // R = V * U^T
    let mut rotation = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            rotation[i][j] = v[i][0] * u[j][0] + v[i][1] * u[j][1] + v[i][2] * u[j][2];
        }
    }

    // Handle reflection case
    let det = rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
            - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
            + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0]);

    if det < 0.0 {
        // Flip sign of last column of V
        let mut v_corrected = v;
        for i in 0..3 {
            v_corrected[i][2] = -v_corrected[i][2];
        }
        for i in 0..3 {
            for j in 0..3 {
                rotation[i][j] = v_corrected[i][0] * u[j][0] + v_corrected[i][1] * u[j][1] + v_corrected[i][2] * u[j][2];
            }
        }
    }

    // Translation = target_centroid - R * source_centroid
    let translation = DVec3::new(
        target_centroid.x - (rotation[0][0] * source_centroid.x + rotation[0][1] * source_centroid.y + rotation[0][2] * source_centroid.z),
        target_centroid.y - (rotation[1][0] * source_centroid.x + rotation[1][1] * source_centroid.y + rotation[1][2] * source_centroid.z),
        target_centroid.z - (rotation[2][0] * source_centroid.x + rotation[2][1] * source_centroid.y + rotation[2][2] * source_centroid.z),
    );

    Some((rotation, translation))
}

/// Compute point-to-plane transformation using linear least squares.
/// Uses the linearized rotation approximation for small angles.
fn compute_point_to_plane_transformation(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    correspondences: &[Correspondence],
) -> Option<([[f64; 3]; 3], DVec3)> {
    if correspondences.len() < 6 {
        return None;
    }

    // Build linear system: A * x = b
    // Where x = [alpha, beta, gamma, tx, ty, tz]^T (rotation angles and translation)
    // For each correspondence: n_i^T * (R(p_i) + t - q_i) = 0
    // Linearized: n_i^T * (p_i + r x p_i + t - q_i) = 0
    // n_i^T * (p_i - q_i) + n_i^T * (r x p_i) + n_i^T * t = 0
    // n_i^T * (p_i - q_i) + (p_i x n_i)^T * r + n_i^T * t = 0

    let _n = correspondences.len();
    let mut ata = [[0.0_f64; 6]; 6];
    let mut atb = [0.0_f64; 6];

    for c in correspondences {
        let p = source[c.source_idx];
        let q = target[c.target_idx];
        let n = target_normals[c.target_idx];

        let cross = DVec3::new(
            p.y * n.z - p.z * n.y,
            p.z * n.x - p.x * n.z,
            p.x * n.y - p.y * n.x,
        );

        let diff = p - q;
        let rhs = -(n.x * diff.x + n.y * diff.y + n.z * diff.z);

        // Row: [cross.x, cross.y, cross.z, n.x, n.y, n.z]
        let row = [cross.x, cross.y, cross.z, n.x, n.y, n.z];

        for i in 0..6 {
            for j in 0..6 {
                ata[i][j] += row[i] * row[j];
            }
            atb[i] += row[i] * rhs;
        }
    }

    // Solve 6x6 system using Gaussian elimination
    let solution = solve_linear_6x6(&ata, &atb)?;

    // Convert angles to rotation matrix
    let (alpha, beta, gamma) = (solution[0], solution[1], solution[2]);
    let rotation = angles_to_rotation_matrix(alpha, beta, gamma);

    let translation = DVec3::new(solution[3], solution[4], solution[5]);

    Some((rotation, translation))
}

fn angles_to_rotation_matrix(alpha: f64, beta: f64, gamma: f64) -> [[f64; 3]; 3] {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let cg = gamma.cos();
    let sg = gamma.sin();

    // R = Rz(gamma) * Ry(beta) * Rx(alpha)
    [
        [cg * cb, cg * sb * sa - sg * ca, cg * sb * ca + sg * sa],
        [sg * cb, sg * sb * sa + cg * ca, sg * sb * ca - cg * sa],
        [-sb, cb * sa, cb * ca],
    ]
}

fn compute_rms_error(
    source: &[DVec3],
    target: &[DVec3],
    correspondences: &[Correspondence],
) -> f64 {
    let mut sum_sq = 0.0;
    for c in correspondences {
        let d = source[c.source_idx] - target[c.target_idx];
        sum_sq += d.length_squared();
    }
    (sum_sq / correspondences.len() as f64).sqrt()
}

fn compute_point_to_plane_error(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    correspondences: &[Correspondence],
) -> f64 {
    let mut sum_sq = 0.0;
    for c in correspondences {
        let diff = source[c.source_idx] - target[c.target_idx];
        let n = target_normals[c.target_idx];
        let dist = diff.dot(n);
        sum_sq += dist * dist;
    }
    (sum_sq / correspondences.len() as f64).sqrt()
}

fn apply_transform(point: DVec3, rotation: &[[f64; 3]; 3], translation: &DVec3) -> DVec3 {
    DVec3::new(
        rotation[0][0] * point.x + rotation[0][1] * point.y + rotation[0][2] * point.z + translation.x,
        rotation[1][0] * point.x + rotation[1][1] * point.y + rotation[1][2] * point.z + translation.y,
        rotation[2][0] * point.x + rotation[2][1] * point.y + rotation[2][2] * point.z + translation.z,
    )
}

fn apply_transform_to_vector(vec: &DVec3, rotation: &[[f64; 3]; 3], translation: &DVec3) -> DVec3 {
    DVec3::new(
        rotation[0][0] * vec.x + rotation[0][1] * vec.y + rotation[0][2] * vec.z + translation.x,
        rotation[1][0] * vec.x + rotation[1][1] * vec.y + rotation[1][2] * vec.z + translation.y,
        rotation[2][0] * vec.x + rotation[2][1] * vec.y + rotation[2][2] * vec.z + translation.z,
    )
}

fn multiply_matrices(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    result
}

/// Simplified SVD for 3x3 matrices using Jacobi iteration.
fn svd_3x3(a: &[[f64; 3]; 3]) -> Option<([[f64; 3]; 3], [f64; 3], [[f64; 3]; 3])> {
    // Compute A^T * A for eigenvalue decomposition
    let mut ata = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                ata[i][j] += a[k][i] * a[k][j];
            }
        }
    }

    // Compute eigenvalues and eigenvectors of A^T * A
    let (eigenvalues, v) = jacobi_eigen(&ata);

    // Compute singular values
    let mut sigma = [0.0; 3];
    for i in 0..3 {
        sigma[i] = eigenvalues[i].max(0.0).sqrt();
    }

    // Compute U = A * V * Sigma^(-1)
    let mut u = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            if sigma[j] > TOLERANCE_LINEAR_ULTRA_STRICT {
                for k in 0..3 {
                    u[i][j] += a[i][k] * v[k][j];
                }
                u[i][j] /= sigma[j];
            }
        }
    }

    // Orthonormalize U
    for j in 0..3 {
        let mut norm = 0.0;
        for i in 0..3 {
            norm += u[i][j] * u[i][j];
        }
        norm = norm.sqrt();
        if norm > TOLERANCE_LINEAR_ULTRA_STRICT {
            for i in 0..3 {
                u[i][j] /= norm;
            }
        }
    }

    Some((u, sigma, v))
}

fn jacobi_eigen(a: &[[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut m = *a;
    let mut v = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    const MAX_ITER: usize = 100;
    const TOL: f64 = TOLERANCE_LEN_MIN;

    for _ in 0..MAX_ITER {
        // Find largest off-diagonal element
        let mut max_val = 0.0;
        let (mut p, mut q) = (0, 1);

        for i in 0..3 {
            for j in (i + 1)..3 {
                if m[i][j].abs() > max_val {
                    max_val = m[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < TOL {
            break;
        }

        // Compute rotation angle
        let theta = if (m[p][p] - m[q][q]).abs() < TOL {
            std::f64::consts::FRAC_PI_4 * m[p][q].signum()
        } else {
            0.5 * (2.0 * m[p][q] / (m[p][p] - m[q][q])).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Apply rotation
        let app = m[p][p];
        let aqq = m[q][q];
        let apq = m[p][q];

        m[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        m[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        m[p][q] = 0.0;
        m[q][p] = 0.0;

        for i in 0..3 {
            if i != p && i != q {
                let aip = m[i][p];
                let aiq = m[i][q];
                m[i][p] = c * aip - s * aiq;
                m[p][i] = m[i][p];
                m[i][q] = s * aip + c * aiq;
                m[q][i] = m[i][q];
            }
        }

        // Update eigenvectors
        for i in 0..3 {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
        }
    }

    // Sort eigenvalues descending
    let mut indexed = [(m[0][0], 0), (m[1][1], 1), (m[2][2], 2)];
    indexed.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));

    let mut eigenvalues = [0.0; 3];
    let mut v_sorted = [[0.0; 3]; 3];
    for (i, &(val, idx)) in indexed.iter().enumerate() {
        eigenvalues[i] = val;
        for j in 0..3 {
            v_sorted[j][i] = v[j][idx];
        }
    }

    (eigenvalues, v_sorted)
}
include!("e1.rs");

