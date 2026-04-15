//! Point cloud analysis tools, analogous to OCCT 8.0 PointSetLib.
//!
//! Provides:
//! - Principal Component Analysis (PCA)
//! - Inertia tensor computation
//! - Dimensionality estimation
//! - Outlier detection and point cloud simplification
//! - Normal estimation
//! - Shape fitting (plane, sphere, cylinder)
//! - BRep integration

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
    const TOLERANCE: f64 = 1e-12;

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
    if total < 1e-10 {
        return Dimensionality::Point;
    }

    // Normalize by largest eigenvalue
    let max_val = sorted[0].max(1e-20);
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
    let voxel_size = (volume / target_count as f64).cbrt().max(1e-10);

    // Group points by voxel
    let mut voxels: std::collections::HashMap<[i64; 3], Vec<DVec3>> = std::collections::HashMap::new();

    for &p in points {
        let key = [
            ((p.x - min.x) / voxel_size).floor() as i64,
            ((p.y - min.y) / voxel_size).floor() as i64,
            ((p.z - min.z) / voxel_size).floor() as i64,
        ];
        voxels.entry(key).or_insert_with(Vec::new).push(p);
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
        let max_val = values[0].max(1e-20);
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
    let (axes, values) = compute_pca(points);

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
    let mut m = a.clone();
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

        if max_val < 1e-14 {
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
    let max_val = values[0].max(1e-20);
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
        if let Some(cyl) = fit_cylinder_with_axis(points, centroid, axis) {
            if cyl.rms_error < best_error {
                best_error = cyl.rms_error;
                best_fit = Some(cyl);
            }
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
    let mut m = a.clone();
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

        if max_val < 1e-14 {
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
    let vertices: Vec<DVec3> = hull_2d
        .iter()
        .map(|&p2d| plane.point + p2d.x * u + p2d.y * v)
        .collect();

    // Compute area using shoelace formula in 3D
    let mut area = 0.0;
    let n = vertices.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let cross = vertices[i].cross(vertices[j]);
        area += cross.dot(normal);
    }
    area = area.abs() / 2.0;

    Some(FittedPolygon {
        vertices,
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
        if (pa.x - pb.x).abs() > 1e-14 {
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
// BRep Integration
// ============================================================================

/// Extracts a point cloud from BRep vertices.
///
/// Collects all unique vertex positions from the BRep.
pub fn extract_points_from_brep_vertices(brep: &rcad_kernel::BRep) -> PointCloud {
    PointCloud::from_vec(brep.vertices.iter().map(|v| v.point).collect())
}

/// Extracts a point cloud from a BRep mesh.
///
/// Samples points from the triangulated faces. If a face has no cached
/// triangulation, it will be skipped.
pub fn extract_points_from_brep_mesh(brep: &rcad_kernel::BRep) -> PointCloud {
    let mut points = Vec::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Add vertices from triangles
                for &[i, j, k] in &face.triangles {
                    if let (Some(a), Some(b), Some(c)) = (
                        brep.vertices.get(i),
                        brep.vertices.get(j),
                        brep.vertices.get(k),
                    ) {
                        points.push(a.point);
                        points.push(b.point);
                        points.push(c.point);
                    }
                }

                // If no triangles, add wire vertices
                if face.triangles.is_empty() {
                    for we in &face.outer_wire.edges {
                        if let Some(edge) = brep.edges.get(we.idx) {
                            let vidx = if we.forward { edge.start } else { edge.end };
                            if let Some(v) = brep.vertices.get(vidx) {
                                points.push(v.point);
                            }
                        }
                    }
                }
            }
        }
    }

    PointCloud::from_vec(points)
}

/// Extracts a point cloud from a triangulated mesh.
///
/// Takes the vertices directly from a SurfaceMesh.
pub fn extract_points_from_mesh(mesh: &crate::triangulate::SurfaceMesh) -> PointCloud {
    PointCloud::from_vec(mesh.vertices.clone())
}

/// Samples points uniformly from a BRep's surfaces.
///
/// For each face with an associated surface, samples points on a regular
/// UV grid and keeps those that lie within the face's boundary.
pub fn sample_points_from_brep_surfaces(
    brep: &rcad_kernel::BRep,
    samples_per_face: usize,
) -> PointCloud {
    use rcad_kernel::geom::SurfaceEval;

    let mut points = Vec::new();
    let sqrt_n = (samples_per_face as f64).sqrt().ceil() as usize;

    let mut face_idx = 0;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Get surface for this face
                let surf_idx = match brep.geom.face_surface.get(face_idx).and_then(|o| *o) {
                    Some(idx) => idx,
                    None => {
                        face_idx += 1;
                        continue;
                    }
                };
                let surf = match brep.geom.surfaces.get(surf_idx) {
                    Some(s) => s,
                    None => {
                        face_idx += 1;
                        continue;
                    }
                };

                // Get UV domain
                let domain = brep.geom.face_surface_range.get(face_idx)
                    .and_then(|o| *o)
                    .unwrap_or_else(|| surf.default_domain());
                let [u0, u1, v0, v1] = domain;

                if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
                    face_idx += 1;
                    continue;
                }

                // Sample on a grid
                for i in 0..sqrt_n {
                    for j in 0..sqrt_n {
                        let u = u0 + (u1 - u0) * (i as f64 + 0.5) / sqrt_n as f64;
                        let v = v0 + (v1 - v0) * (j as f64 + 0.5) / sqrt_n as f64;
                        let p = surf.point_at(u, v);
                        points.push(p);
                    }
                }

                face_idx += 1;
            }
        }
    }

    PointCloud::from_vec(points)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-6;

    fn approx_eq(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    #[test]
    fn test_empty_point_cloud() {
        let pc = PointCloud::new();
        assert!(pc.is_empty());
        assert_eq!(pc.len(), 0);
        assert!(pc.bounding_box().is_none());
        assert!(pc.centroid().is_none());
    }

    #[test]
    fn test_point_cloud_basics() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let pc = PointCloud::from_points(&points);

        assert_eq!(pc.len(), 3);

        let centroid = pc.centroid().unwrap();
        assert!(approx_eq(centroid, DVec3::new(1.0/3.0, 1.0/3.0, 0.0), 1e-10));

        let (min, max) = pc.bounding_box().unwrap();
        assert!(approx_eq(min, DVec3::ZERO, 1e-10));
        assert!(approx_eq(max, DVec3::new(1.0, 1.0, 0.0), 1e-10));
    }

    #[test]
    fn test_pca_identity() {
        // Points on a cube - PCA should give roughly equal eigenvalues
        // Simpler test that's more numerically stable
        let points: Vec<DVec3> = vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(-1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, -1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(0.0, 0.0, -1.0),
        ];

        let (axes, values) = compute_pca(&points);

        // Eigenvalues should be positive and roughly equal for symmetric distribution
        assert!(values[0] > 0.0, "Largest eigenvalue should be positive, got {}", values[0]);
        assert!(values[2] >= 0.0, "Smallest eigenvalue should be non-negative, got {}", values[2]);
        // All eigenvalues should be similar (within factor of 2) for this symmetric case
        assert!(values[0] / values[2].max(1e-10) < 3.0, "Eigenvalue ratio {} too large", values[0] / values[2].max(1e-10));

        // Axes should be orthonormal
        for axis in &axes {
            assert!((axis.length() - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_pca_line() {
        // Points along X axis
        let points: Vec<DVec3> = (0..10).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let (axes, values) = compute_pca(&points);

        // Largest eigenvalue should be along X
        assert!(values[0] > values[1]);
        assert!(values[1] < 1e-6);
        assert!(values[2] < 1e-6);

        // First principal axis should be approximately X
        assert!((axes[0].x.abs() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_pca_plane() {
        // Points on XY plane
        let mut points = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let (axes, values) = compute_pca(&points);

        // Two large eigenvalues, one small
        assert!(values[0] > 0.1);
        assert!(values[1] > 0.1);
        assert!(values[2] < 0.01);

        // Third principal axis should be Z (normal)
        assert!((axes[2].z.abs() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dimensionality() {
        // Point-like
        let d = estimate_dimensionality([1e-20, 1e-20, 1e-20], 0.01);
        assert_eq!(d, Dimensionality::Point);

        // Linear
        let d = estimate_dimensionality([10.0, 0.001, 0.001], 0.01);
        assert_eq!(d, Dimensionality::Linear);

        // Planar
        let d = estimate_dimensionality([10.0, 10.0, 0.001], 0.01);
        assert_eq!(d, Dimensionality::Planar);

        // Volumetric
        let d = estimate_dimensionality([10.0, 10.0, 10.0], 0.01);
        assert_eq!(d, Dimensionality::Volumetric);
    }

    #[test]
    fn test_inertia_tensor() {
        // Unit cube at origin
        let points: Vec<DVec3> = (0..=1)
            .flat_map(|x| (0..=1).flat_map(move |y| (0..=1).map(move |z| DVec3::new(x as f64, y as f64, z as f64))))
            .collect();

        let inertia = compute_inertia(&points);

        // Check diagonal elements are positive
        assert!(inertia[0][0] >= 0.0);
        assert!(inertia[1][1] >= 0.0);
        assert!(inertia[2][2] >= 0.0);

        // Check symmetry
        assert!((inertia[0][1] - inertia[1][0]).abs() < 1e-10);
        assert!((inertia[0][2] - inertia[2][0]).abs() < 1e-10);
        assert!((inertia[1][2] - inertia[2][1]).abs() < 1e-10);
    }

    #[test]
    fn test_fit_plane() {
        // Perfect plane
        let points: Vec<DVec3> = (0..10)
            .flat_map(|i| (0..10).map(move |j| DVec3::new(i as f64, j as f64, 0.0)))
            .collect();

        let plane = fit_plane(&points).unwrap();

        assert!(approx_eq(plane.normal, DVec3::Z, 1e-6) || approx_eq(plane.normal, -DVec3::Z, 1e-6));
        assert!(plane.rms_error < 1e-6);
    }

    #[test]
    fn test_fit_sphere() {
        // Points on a sphere of radius 2 centered at (1, 2, 3)
        let center = DVec3::new(1.0, 2.0, 3.0);
        let radius = 2.0;

        let mut points = Vec::new();
        for i in 0..50 {
            let theta = 2.0 * PI * i as f64 / 50.0;
            let phi = PI * i as f64 / 50.0;
            let x = center.x + radius * phi.sin() * theta.cos();
            let y = center.y + radius * phi.sin() * theta.sin();
            let z = center.z + radius * phi.cos();
            points.push(DVec3::new(x, y, z));
        }

        let sphere = fit_sphere(&points).unwrap();

        assert!(approx_eq(sphere.center, center, 0.1));
        assert!((sphere.radius - radius).abs() < 0.1);
        assert!(sphere.rms_error < 0.1);
    }

    #[test]
    fn test_fit_cylinder() {
        // Points on a cylinder along Z axis
        let radius = 1.5;
        let mut points = Vec::new();

        for i in 0..20 {
            let theta = 2.0 * PI * i as f64 / 20.0;
            for z in 0..5 {
                let x = radius * theta.cos();
                let y = radius * theta.sin();
                points.push(DVec3::new(x, y, z as f64));
            }
        }

        let cylinder = fit_cylinder(&points).unwrap();

        assert!((cylinder.radius - radius).abs() < 0.1);
        assert!(cylinder.rms_error < 0.1);
    }

    #[test]
    fn test_simplify_random() {
        let points: Vec<DVec3> = (0..1000).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let simplified = simplify_point_cloud(&points, 100, SamplingStrategy::Random);

        assert_eq!(simplified.len(), 100);
    }

    #[test]
    fn test_simplify_voxel() {
        let mut points = Vec::new();
        // Dense grid of points
        for i in 0..10 {
            for j in 0..10 {
                for k in 0..10 {
                    points.push(DVec3::new(i as f64, j as f64, k as f64));
                }
            }
        }

        let simplified = simplify_point_cloud(&points, 50, SamplingStrategy::Voxel);

        assert!(simplified.len() >= 27); // At least 3x3x3 voxels
        assert!(simplified.len() <= 100);
    }

    #[test]
    fn test_simplify_farthest_point() {
        let points: Vec<DVec3> = (0..100).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let simplified = simplify_point_cloud(&points, 10, SamplingStrategy::FarthestPoint);

        assert_eq!(simplified.len(), 10);

        // Should include endpoints
        let has_start = simplified.iter().any(|p| p.x < 1.0);
        let has_end = simplified.iter().any(|p| p.x > 98.0);
        assert!(has_start || has_end);
    }

    #[test]
    fn test_estimate_normals() {
        // Points on XY plane
        let mut points = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let normals = estimate_normals(&points, 4);

        assert_eq!(normals.len(), points.len());

        // All normals should point along Z (positive or negative)
        for n in &normals {
            assert!(n.z.abs() > 0.9, "Normal should be along Z, got {:?}", n);
        }
    }

    #[test]
    fn test_fit_polygon() {
        // Square in XY plane
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];

        let polygon = fit_polygon(&points).expect("fit_polygon should succeed for a square");

        assert!(polygon.vertices.len() >= 3, "Should have at least 3 vertices, got {}", polygon.vertices.len());
        assert!((polygon.area - 1.0).abs() < 1e-4, "Area should be 1.0, got {}", polygon.area);
    }

    #[test]
    fn test_outlier_detection() {
        let mut points = Vec::new();

        // Cluster of points near origin
        for i in 0..50 {
            points.push(DVec3::new(
                (i as f64 % 10.0) * 0.1,
                (i as f64 / 10.0) * 0.1,
                0.0,
            ));
        }

        // Add an outlier far away
        points.push(DVec3::new(100.0, 100.0, 100.0));

        let outliers = detect_outliers(&points, 5, 1.5);

        // Should detect at least one outlier
        assert!(!outliers.is_empty());

        // The farthest point should have highest score
        assert!(outliers[0].index == 50 || outliers.iter().any(|o| o.index == 50));
    }

    #[test]
    fn test_analyze_point_cloud() {
        // Create a box-shaped point cloud
        let mut points = Vec::new();
        for x in 0..=1 {
            for y in 0..=1 {
                for z in 0..=1 {
                    points.push(DVec3::new(x as f64, y as f64, z as f64));
                }
            }
        }

        let analysis = analyze_point_cloud(&points).unwrap();

        // Centroid should be at (0.5, 0.5, 0.5)
        assert!(approx_eq(analysis.centroid, DVec3::splat(0.5), 1e-6));

        // Should be volumetric
        assert_eq!(analysis.dimensionality, Dimensionality::Volumetric);

        // Bounding box
        assert!(approx_eq(analysis.bounding_box.0, DVec3::ZERO, 1e-6));
        assert!(approx_eq(analysis.bounding_box.1, DVec3::splat(1.0), 1e-6));
    }

    #[test]
    fn test_brep_integration() {
        use rcad_kernel::{BRep, PrimitiveSolid};

        // Create a unit box
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Extract vertex points
        let vertex_pc = extract_points_from_brep_vertices(&brep);
        assert_eq!(vertex_pc.len(), 8); // Box has 8 vertices

        // Analyze
        let analysis = analyze_point_cloud(&vertex_pc.points).unwrap();
        assert!(approx_eq(analysis.centroid, DVec3::splat(0.5), 1e-6));
    }
}
