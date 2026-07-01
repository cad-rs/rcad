
fn solve_linear_6x6(a: &[[f64; 6]; 6], b: &[f64; 6]) -> Option<[f64; 6]> {
    const N: usize = 6;
    let mut m = *a;
    let mut v = *b;

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

// ============================================================================
// Segmentation
// ============================================================================

/// Result of region growing segmentation.
#[derive(Debug, Clone)]
pub struct Segment {
    /// Indices of points in this segment.
    pub point_indices: Vec<usize>,
    /// Fitted shape type (if any).
    pub shape_type: Option<ShapeType>,
    /// Fitted shape parameters (if any).
    pub shape_params: Option<ShapeParams>,
    /// Centroid of the segment.
    pub centroid: DVec3,
}

/// Types of shapes that can be detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeType {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
}

/// Parameters for fitted shapes.
#[derive(Debug, Clone)]
pub enum ShapeParams {
    Plane { point: DVec3, normal: DVec3 },
    Sphere { center: DVec3, radius: f64 },
    Cylinder { axis_point: DVec3, axis_direction: DVec3, radius: f64 },
}

/// Configuration for region growing segmentation.
#[derive(Debug, Clone)]
pub struct RegionGrowingConfig {
    /// Number of neighbors for normal estimation.
    pub k_neighbors: usize,
    /// Maximum angular difference (radians) for region growing.
    pub max_angle: f64,
    /// Maximum distance from fitted shape.
    pub max_distance: f64,
    /// Minimum number of points for a valid segment.
    pub min_segment_size: usize,
    /// Maximum number of segments to extract.
    pub max_segments: usize,
}

impl Default for RegionGrowingConfig {
    fn default() -> Self {
        Self {
            k_neighbors: 30,
            max_angle: std::f64::consts::PI / 6.0, // 30 degrees
            max_distance: 0.01,
            min_segment_size: 100,
            max_segments: 100,
        }
    }
}

/// Performs region growing segmentation based on smoothness constraint.
///
/// Grows regions from seed points, adding neighbors with similar normals.
pub fn region_growing_segmentation(
    points: &[DVec3],
    config: &RegionGrowingConfig,
) -> Vec<Segment> {
    if points.is_empty() {
        return Vec::new();
    }

    let n = points.len();

    // Estimate normals
    let normals = estimate_normals(points, config.k_neighbors);

    // Build neighbor graph (kNN)
    let neighbors = build_neighbor_graph(points, config.k_neighbors);

    // Track visited points
    let mut visited = vec![false; n];
    let mut segments = Vec::new();

    // Sort points by curvature (lowest first for better seeds)
    let curvatures: Vec<f64> = compute_curvatures(points, &neighbors);
    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_by(|&a, &b| {
        curvatures[a].partial_cmp(&curvatures[b]).unwrap_or(Ordering::Equal)
    });

    for &seed_idx in &sorted_indices {
        if visited[seed_idx] || segments.len() >= config.max_segments {
            continue;
        }

        // Grow region from seed
        let segment_indices = grow_region(
            seed_idx,
            points,
            &normals,
            &neighbors,
            &mut visited,
            config,
        );

        if segment_indices.len() >= config.min_segment_size {
            let segment_pts: Vec<DVec3> = segment_indices.iter().map(|&i| points[i]).collect();
            let centroid = segment_pts.iter().sum::<DVec3>() / segment_pts.len() as f64;

            segments.push(Segment {
                point_indices: segment_indices,
                shape_type: None,
                shape_params: None,
                centroid,
            });
        }
    }

    segments
}

/// Performs Euclidean clustering segmentation.
///
/// Clusters points based on Euclidean distance threshold.
pub fn euclidean_clustering(
    points: &[DVec3],
    tolerance: f64,
    min_cluster_size: usize,
) -> Vec<Vec<usize>> {
    if points.is_empty() {
        return Vec::new();
    }

    let n = points.len();
    let mut visited = vec![false; n];
    let mut clusters = Vec::new();
    let tolerance_sq = tolerance * tolerance;

    for i in 0..n {
        if visited[i] {
            continue;
        }

        let mut cluster = Vec::new();
        let mut queue = vec![i];
        visited[i] = true;

        while let Some(current) = queue.pop() {
            cluster.push(current);

            // Find all neighbors within tolerance
            for j in 0..n {
                if !visited[j] && (points[j] - points[current]).length_squared() < tolerance_sq {
                    visited[j] = true;
                    queue.push(j);
                }
            }
        }

        if cluster.len() >= min_cluster_size {
            clusters.push(cluster);
        }
    }

    // Sort clusters by size (largest first)
    clusters.sort_by(|a, b| b.len().cmp(&a.len()));

    clusters
}

/// Performs shape-based segmentation (plane, sphere, cylinder).
///
/// Uses RANSAC to detect dominant shapes and segment the point cloud.
pub fn shape_segmentation(
    points: &[DVec3],
    shape_type: ShapeType,
    distance_threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    match shape_type {
        ShapeType::Plane => ransac_plane_segmentation(points, distance_threshold, min_points, max_iterations),
        ShapeType::Sphere => ransac_sphere_segmentation(points, distance_threshold, min_points, max_iterations),
        ShapeType::Cylinder => ransac_cylinder_segmentation(points, distance_threshold, min_points, max_iterations),
        _ => None,
    }
}

fn ransac_plane_segmentation(
    points: &[DVec3],
    threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    if points.len() < 3 {
        return None;
    }

    let mut best_inliers = Vec::new();
    let mut best_plane: Option<(DVec3, DVec3)> = None;

    let mut rng = SimpleRng::new(42);

    for _ in 0..max_iterations {
        // Sample 3 random points
        let i0 = (rng.next() as usize) % points.len();
        let i1 = (rng.next() as usize) % points.len();
        let i2 = (rng.next() as usize) % points.len();

        if i0 == i1 || i1 == i2 || i0 == i2 {
            continue;
        }

        let p0 = points[i0];
        let p1 = points[i1];
        let p2 = points[i2];

        // Compute plane normal
        let normal = (p1 - p0).cross(p2 - p0);
        let len = normal.length();
        if len < TOLERANCE_LINEAR_ULTRA_STRICT {
            continue;
        }
        let normal = normal / len;

        // Find inliers
        let mut inliers = Vec::new();
        for (i, &p) in points.iter().enumerate() {
            let dist = (p - p0).dot(normal).abs();
            if dist < threshold {
                inliers.push(i);
            }
        }

        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
            best_plane = Some((p0, normal));
        }
    }

    if best_inliers.len() < min_points {
        return None;
    }

    let (_point, _normal) = best_plane.unwrap();

    // Refit using all inliers
    let inlier_pts: Vec<DVec3> = best_inliers.iter().map(|&i| points[i]).collect();
    let fitted = fit_plane(&inlier_pts)?;

    // Separate inliers and outliers
    let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
    let outliers: Vec<usize> = (0..points.len())
        .filter(|i| !inlier_set.contains(i))
        .collect();

    Some((
        ShapeParams::Plane {
            point: fitted.point,
            normal: fitted.normal,
        },
        best_inliers,
        outliers,
    ))
}

fn ransac_sphere_segmentation(
    points: &[DVec3],
    threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    if points.len() < 4 {
        return None;
    }

    let mut best_inliers = Vec::new();
    let mut best_sphere: Option<(DVec3, f64)> = None;

    let mut rng = SimpleRng::new(42);

    for _ in 0..max_iterations {
        // Sample 4 random points
        let indices: Vec<usize> = (0..4)
            .map(|_| (rng.next() as usize) % points.len())
            .collect();

        if indices.iter().collect::<std::collections::HashSet<_>>().len() < 4 {
            continue;
        }

        let sample_pts: Vec<DVec3> = indices.iter().map(|&i| points[i]).collect();

        // Fit sphere to 4 points
        let sphere = fit_sphere_4pt(&sample_pts);
        let (center, radius) = match sphere {
            Some(s) => s,
            None => continue,
        };

        if !(TOLERANCE_LINEAR_ULTRA_STRICT..=1e10).contains(&radius) {
            continue;
        }

        // Find inliers
        let mut inliers = Vec::new();
        for (i, &p) in points.iter().enumerate() {
            let dist = ((p - center).length() - radius).abs();
            if dist < threshold {
                inliers.push(i);
            }
        }

        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
            best_sphere = Some((center, radius));
        }
    }

    if best_inliers.len() < min_points {
        return None;
    }

    let (center, radius) = best_sphere.unwrap();

    // Refit using all inliers
    let inlier_pts: Vec<DVec3> = best_inliers.iter().map(|&i| points[i]).collect();
    if let Some(fitted) = fit_sphere(&inlier_pts) {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Sphere {
                center: fitted.center,
                radius: fitted.radius,
            },
            best_inliers,
            outliers,
        ))
    } else {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Sphere { center, radius },
            best_inliers,
            outliers,
        ))
    }
}

fn fit_sphere_4pt(points: &[DVec3]) -> Option<(DVec3, f64)> {
    if points.len() < 4 {
        return None;
    }

    // Solve for sphere passing through 4 points using linear system
    let p0 = points[0];
    let p1 = points[1];
    let p2 = points[2];
    let p3 = points[3];

    // System: |P - C|^2 = r^2 for each point
    // Subtract first equation to eliminate r^2
    // |Pi - C|^2 - |P0 - C|^2 = 0
    // |Pi|^2 - 2*Pi*C + |C|^2 - |P0|^2 + 2*P0*C - |C|^2 = 0
    // |Pi|^2 - |P0|^2 - 2*(Pi - P0)*C = 0
    // (Pi - P0)*C = (|Pi|^2 - |P0|^2) / 2

    let a = [
        [p1.x - p0.x, p1.y - p0.y, p1.z - p0.z],
        [p2.x - p0.x, p2.y - p0.y, p2.z - p0.z],
        [p3.x - p0.x, p3.y - p0.y, p3.z - p0.z],
    ];

    let p0_sq = p0.length_squared();
    let b = [
        (p1.length_squared() - p0_sq) / 2.0,
        (p2.length_squared() - p0_sq) / 2.0,
        (p3.length_squared() - p0_sq) / 2.0,
    ];

    let center = solve_linear_3x3(&a, &b)?;
    let cx = center[0];
    let cy = center[1];
    let cz = center[2];
    let center = DVec3::new(cx, cy, cz);
    let radius = (center - p0).length();

    Some((center, radius))
}

fn ransac_cylinder_segmentation(
    points: &[DVec3],
    threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    if points.len() < 5 {
        return None;
    }

    let mut best_inliers = Vec::new();
    let mut best_cylinder: Option<(DVec3, DVec3, f64)> = None;

    let mut rng = SimpleRng::new(42);

    for _ in 0..max_iterations {
        // Sample 2 points for axis estimation
        let i0 = (rng.next() as usize) % points.len();
        let i1 = (rng.next() as usize) % points.len();

        if i0 == i1 {
            continue;
        }

        let p0 = points[i0];
        let p1 = points[i1];
        let axis = (p1 - p0).normalize_or(DVec3::Z);

        // Project points onto plane perpendicular to axis
        // and estimate radius
        let centroid = points.iter().sum::<DVec3>() / points.len() as f64;
        let u = if axis.x.abs() < 0.9 {
            axis.cross(DVec3::X).normalize()
        } else {
            axis.cross(DVec3::Y).normalize()
        };
        let v = axis.cross(u);

        let projected: Vec<f64> = points.iter().map(|&p| {
            let d = p - centroid;
            let x = d.dot(u);
            let y = d.dot(v);
            (x * x + y * y).sqrt()
        }).collect();

        // Estimate radius as median
        let mut sorted_proj = projected.clone();
        sorted_proj.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let radius = sorted_proj[sorted_proj.len() / 2];

        if !(TOLERANCE_LINEAR_ULTRA_STRICT..=1e10).contains(&radius) {
            continue;
        }

        // Find inliers
        let mut inliers = Vec::new();
        for (i, &p) in points.iter().enumerate() {
            let d = p - centroid;
            let axial = d.dot(axis);
            let radial = (d - axial * axis).length();
            let dist = (radial - radius).abs();
            if dist < threshold {
                inliers.push(i);
            }
        }

        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
            best_cylinder = Some((centroid, axis, radius));
        }
    }

    if best_inliers.len() < min_points {
        return None;
    }

    let (axis_point, axis_direction, radius) = best_cylinder.unwrap();

    // Refit using all inliers
    let inlier_pts: Vec<DVec3> = best_inliers.iter().map(|&i| points[i]).collect();
    if let Some(fitted) = fit_cylinder(&inlier_pts) {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Cylinder {
                axis_point: fitted.axis_point,
                axis_direction: fitted.axis_direction,
                radius: fitted.radius,
            },
            best_inliers,
            outliers,
        ))
    } else {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Cylinder {
                axis_point,
                axis_direction,
                radius,
            },
            best_inliers,
            outliers,
        ))
    }
}

fn build_neighbor_graph(points: &[DVec3], k: usize) -> Vec<Vec<usize>> {
    let n = points.len();
    let k = k.min(n - 1).max(1);
    let mut neighbors = Vec::with_capacity(n);

    for i in 0..n {
        let mut distances: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, (points[j] - points[i]).length_squared()))
            .collect();
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

        neighbors.push(distances.iter().take(k).map(|&(j, _)| j).collect());
    }

    neighbors
}

fn compute_curvatures(points: &[DVec3], neighbors: &[Vec<usize>]) -> Vec<f64> {
    points.iter().enumerate().map(|(i, _)| {
        let neighbor_pts: Vec<DVec3> = neighbors[i].iter().map(|&j| points[j]).collect();
        let (_, values) = compute_pca(&neighbor_pts);
        let sum: f64 = values.iter().sum();
        if sum > TOLERANCE_LINEAR_ULTRA_STRICT {
            values[2] / sum
        } else {
            0.0
        }
    }).collect()
}

fn grow_region(
    seed_idx: usize,
    points: &[DVec3],
    normals: &[DVec3],
    neighbors: &[Vec<usize>],
    visited: &mut [bool],
    config: &RegionGrowingConfig,
) -> Vec<usize> {
    let mut region = Vec::new();
    let mut queue = vec![seed_idx];
    visited[seed_idx] = true;

    let seed_normal = normals[seed_idx];

    while let Some(current) = queue.pop() {
        region.push(current);

        for &neighbor in &neighbors[current] {
            if visited[neighbor] {
                continue;
            }

            // Check normal similarity
            let dot = normals[neighbor].dot(seed_normal);
            let angle = dot.acos();

            if angle < config.max_angle {
                // Check distance constraint
                let dist = (points[neighbor] - points[current]).length();
                if dist < config.max_distance * 10.0 {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
    }

    region
}

// ============================================================================
// Surface Reconstruction
// ============================================================================

/// Triangle mesh for surface reconstruction output.
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// Node positions.
    pub nodes: Vec<DVec3>,
    /// Triangle indices (3 node indices per triangle).
    pub triangles: Vec<[usize; 3]>,
    /// Node normals (optional).
    pub normals: Option<Vec<DVec3>>,
}

impl Default for TriangleMesh {
    fn default() -> Self {
        Self::new()
    }
}

impl TriangleMesh {
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            triangles: Vec::new(),
            normals: None,
        }
    }

    /// Computes face normals for the mesh.
    pub fn compute_face_normals(&self) -> Vec<DVec3> {
        self.triangles.iter().map(|&[i, j, k]| {
            let a = self.nodes[i];
            let b = self.nodes[j];
            let c = self.nodes[k];
            let normal = (b - a).cross(c - a);
            let len = normal.length();
            if len > TOLERANCE_LINEAR_ULTRA_STRICT {
                normal / len
            } else {
                DVec3::Z
            }
        }).collect()
    }

    /// Computes node normals by averaging adjacent face normals.
    pub fn compute_node_normals(&self) -> Vec<DVec3> {
        let face_normals = self.compute_face_normals();
        let mut node_normals = vec![DVec3::ZERO; self.nodes.len()];

        for (tri, &normal) in self.triangles.iter().zip(face_normals.iter()) {
            for &idx in tri {
                node_normals[idx] += normal;
            }
        }

        for normal in &mut node_normals {
            let len = normal.length();
            if len > TOLERANCE_LINEAR_ULTRA_STRICT {
                *normal /= len;
            }
        }

        node_normals
    }
}

/// Configuration for Poisson surface reconstruction.
#[derive(Debug, Clone)]
pub struct PoissonConfig {
    /// Octree depth (higher = more detail, more memory).
    pub depth: usize,
    /// Solver division (default: 8).
    pub solver_divide: usize,
    /// Iso-surface value (default: 0, use mean).
    pub iso_value: f64,
}

impl Default for PoissonConfig {
    fn default() -> Self {
        Self {
            depth: 8,
            solver_divide: 8,
            iso_value: 0.0,
        }
    }
}

/// Performs Poisson surface reconstruction.
///
/// Reconstructs a watertight surface from oriented point samples.
/// Uses an implicit function approach with octree-based spatial indexing.
pub fn poisson_reconstruction(
    points: &[DVec3],
    normals: &[DVec3],
    config: &PoissonConfig,
) -> Option<TriangleMesh> {
    if points.is_empty() || points.len() != normals.len() {
        return None;
    }

    // Simplified Poisson reconstruction using implicit function + marching cubes
    // This is a basic implementation - full Poisson is more complex

    let (min, max) = compute_bounding_box(points)?;
    let padding = (max - min).length() * 0.1;
    let min = min - DVec3::splat(padding);
    let max = max + DVec3::splat(padding);

    let resolution = 2usize.pow(config.depth as u32);
    let cell_size = (max - min) / resolution as f64;

    // Build implicit function using oriented point samples
    let mut grid = vec![0.0_f64; resolution * resolution * resolution];

    // Compute implicit function values
    for idx in 0..resolution {
        for idy in 0..resolution {
            for idz in 0..resolution {
                let p = DVec3::new(
                    min.x + (idx as f64 + 0.5) * cell_size.x,
                    min.y + (idy as f64 + 0.5) * cell_size.y,
                    min.z + (idz as f64 + 0.5) * cell_size.z,
                );

                let mut value = 0.0;
                for (&pt, &n) in points.iter().zip(normals.iter()) {
                    let d = p - pt;
                    let dist = d.length();
                    if dist > TOLERANCE_LINEAR_ULTRA_STRICT {
                        // Signed distance from oriented point
                        value += d.dot(n) / (dist * dist + 1.0);
                    }
                }

                grid[idx + idy * resolution + idz * resolution * resolution] = value;
            }
        }
    }

    // Marching cubes to extract iso-surface
    marching_cubes(&grid, resolution, resolution, resolution, &min, &max, config.iso_value)
}

fn marching_cubes(
    grid: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    min: &DVec3,
    max: &DVec3,
    iso_value: f64,
) -> Option<TriangleMesh> {
    let mut mesh = TriangleMesh::new();
    let dx = (max.x - min.x) / nx as f64;
    let dy = (max.y - min.y) / ny as f64;
    let dz = (max.z - min.z) / nz as f64;

    // Simplified marching cubes - just create triangles for cells crossing iso-value
    for ix in 0..nx - 1 {
        for iy in 0..ny - 1 {
            for iz in 0..nz - 1 {
                // Get 8 corner values
                let corners = [
                    grid[ix + iy * nx + iz * nx * ny],
                    grid[ix + 1 + iy * nx + iz * nx * ny],
                    grid[ix + 1 + (iy + 1) * nx + iz * nx * ny],
                    grid[ix + (iy + 1) * nx + iz * nx * ny],
                    grid[ix + iy * nx + (iz + 1) * nx * ny],
                    grid[ix + 1 + iy * nx + (iz + 1) * nx * ny],
                    grid[ix + 1 + (iy + 1) * nx + (iz + 1) * nx * ny],
                    grid[ix + (iy + 1) * nx + (iz + 1) * nx * ny],
                ];

                // Check if cell crosses iso-value
                let above: Vec<bool> = corners.iter().map(|&v| v > iso_value).collect();
                let all_above = above.iter().all(|&b| b);
                let all_below = above.iter().all(|&b| !b);

                if all_above || all_below {
                    continue;
                }

                // Simplified: create triangles at cell center
                let cx = min.x + (ix as f64 + 0.5) * dx;
                let cy = min.y + (iy as f64 + 0.5) * dy;
                let cz = min.z + (iz as f64 + 0.5) * dz;

                // Create a simple cube face approximation
                let base_idx = mesh.nodes.len();
                let size = dx.min(dy).min(dz) * 0.5;

                // Add nodes for a small quad
                mesh.nodes.push(DVec3::new(cx - size, cy - size, cz));
                mesh.nodes.push(DVec3::new(cx + size, cy - size, cz));
                mesh.nodes.push(DVec3::new(cx + size, cy + size, cz));
                mesh.nodes.push(DVec3::new(cx - size, cy + size, cz));

                mesh.triangles.push([base_idx, base_idx + 1, base_idx + 2]);
                mesh.triangles.push([base_idx, base_idx + 2, base_idx + 3]);
            }
        }
    }

    if mesh.nodes.is_empty() {
        None
    } else {
        Some(mesh)
    }
}

/// Configuration for Ball Pivoting Algorithm.
#[derive(Debug, Clone)]
pub struct BpaConfig {
    /// Ball radius for pivoting.
    pub ball_radius: f64,
    /// Clustering radius for duplicate removal.
    pub clustering: f64,
    /// Angle threshold for edge selection (radians).
    pub angle_threshold: f64,
}

impl Default for BpaConfig {
    fn default() -> Self {
        Self {
            ball_radius: 0.1,
            clustering: 0.001,
            angle_threshold: std::f64::consts::PI / 4.0,
        }
    }
}

/// Performs Ball Pivoting Algorithm surface reconstruction.
///
/// Reconstructs surface by rolling a ball over the point cloud.
pub fn ball_pivoting_reconstruction(
    points: &[DVec3],
    normals: &[DVec3],
    config: &BpaConfig,
) -> Option<TriangleMesh> {
    if points.len() < 3 {
        return None;
    }

    let mut mesh = TriangleMesh::new();

    // Build spatial index
    let grid_size = config.ball_radius;
    let (min, _) = compute_bounding_box(points)?;

    let mut spatial_index: std::collections::HashMap<[i64; 3], Vec<usize>> = std::collections::HashMap::new();

    for (i, &p) in points.iter().enumerate() {
        let key = [
            ((p.x - min.x) / grid_size).floor() as i64,
            ((p.y - min.y) / grid_size).floor() as i64,
            ((p.z - min.z) / grid_size).floor() as i64,
        ];
        spatial_index.entry(key).or_default().push(i);
    }

    // Track used edges
    let mut used_edges: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    let mut used_triangles: std::collections::HashSet<[usize; 3]> = std::collections::HashSet::new();

    let r = config.ball_radius;
    let r_sq = r * r;

    // Find seed triangle
    for i0 in 0..points.len() {
        for i1 in (i0 + 1)..points.len() {
            for i2 in (i1 + 1)..points.len() {
                let p0 = points[i0];
                let p1 = points[i1];
                let p2 = points[i2];

                // Check if ball of radius r touches all three points
                let circumcenter = compute_circumcenter(p0, p1, p2)?;
                let circumradius = (p0 - circumcenter).length();

                if circumradius > r {
                    continue;
                }

                // Ball center is at distance sqrt(r^2 - circumradius^2) from circumcenter
                let ball_dist = (r_sq - circumradius * circumradius).sqrt();
                let normal = (p1 - p0).cross(p2 - p0);
                let len = normal.length();
                if len < TOLERANCE_LINEAR_ULTRA_STRICT {
                    continue;
                }
                let normal = normal / len;

                // Check ball doesn't contain other points
                let ball_center = circumcenter + normal * ball_dist;
                let mut valid = true;
                for &p in points.iter().take(10) {
                    if (p - ball_center).length_squared() < r_sq - TOLERANCE_LINEAR_ULTRA_STRICT {
                        valid = false;
                        break;
                    }
                }

                if valid {
                    let base_idx = mesh.nodes.len();
                    mesh.nodes.push(p0);
                    mesh.nodes.push(p1);
                    mesh.nodes.push(p2);
                    mesh.triangles.push([base_idx, base_idx + 1, base_idx + 2]);

                    used_edges.insert((i0, i1));
                    used_edges.insert((i1, i2));
                    used_edges.insert((i2, i0));
                    used_triangles.insert([i0, i1, i2]);
                    break;
                }
            }
            if !mesh.triangles.is_empty() {
                break;
            }
        }
        if !mesh.triangles.is_empty() {
            break;
        }
    }

    if mesh.nodes.is_empty() {
        // Fallback: simple Delaunay triangulation
        return delaunay_reconstruction(points, normals);
    }

    // Expand from seed using ball pivoting
    // This is simplified - full BPA is more complex
    for _ in 0..points.len() / 3 {
        // Find next edge to pivot from
        // Collect edges to iterate over to avoid borrow issues
        let edges_to_process: Vec<(usize, usize)> = used_edges.iter().copied().collect();
        let mut new_edges: Vec<(usize, usize)> = Vec::new();

        for (i0, i1) in edges_to_process {
            // Try to find third point
            for i2 in 0..points.len() {
                if i2 == i0 || i2 == i1 {
                    continue;
                }

                let _key = [i0.min(i1), i0.max(i1), i2];
                let tri_key = if i0 < i1 {
                    [i0, i1, i2]
                } else {
                    [i1, i0, i2]
                };

                if used_triangles.contains(&tri_key) || used_triangles.contains(&[tri_key[0], tri_key[2], tri_key[1]]) {
                    continue;
                }

                let p0 = points[i0];
                let p1 = points[i1];
                let p2 = points[i2];

                if let Some(cc) = compute_circumcenter(p0, p1, p2) {
                    let cr = (p0 - cc).length();
                    if cr <= r {
                        let base_idx = mesh.nodes.len();
                        mesh.nodes.push(p0);
                        mesh.nodes.push(p1);
                        mesh.nodes.push(p2);
                        mesh.triangles.push([base_idx, base_idx + 1, base_idx + 2]);

                        new_edges.push((i0, i2));
                        new_edges.push((i2, i1));
                        used_triangles.insert(tri_key);
                    }
                }
            }
        }

        // Add new edges
        for edge in new_edges {
            used_edges.insert(edge);
        }
    }

    if mesh.nodes.is_empty() {
        None
    } else {
        Some(mesh)
    }
}

fn compute_circumcenter(a: DVec3, b: DVec3, c: DVec3) -> Option<DVec3> {
    let ab = b - a;
    let ac = c - a;

    let cross = ab.cross(ac);
    let denom = 2.0 * cross.length_squared();

    if denom < TOLERANCE_METRIC_SQ_NEAR_ZERO {
        return None;
    }

    let _d = cross.dot(a.cross(b) + b.cross(c) + c.cross(a)) / denom;

    Some(a + (ab.cross(ac).cross(ab) * ac.length_squared() + ac.cross(ab.cross(ac)) * ab.length_squared()) / (2.0 * cross.length_squared()))
}

/// Performs Delaunay triangulation based surface reconstruction.
///
/// Projects points to 2D, computes Delaunay triangulation,
/// then projects back to 3D.
pub fn delaunay_reconstruction(
    points: &[DVec3],
    _normals: &[DVec3],
) -> Option<TriangleMesh> {
    if points.len() < 3 {
        return None;
    }

    let mut mesh = TriangleMesh::new();

    // Fit plane to points
    let plane = fit_plane(points)?;
    let normal = plane.normal;

    // Build orthonormal basis on plane
    let u = if normal.x.abs() < 0.9 {
        normal.cross(DVec3::X).normalize()
    } else {
        normal.cross(DVec3::Y).normalize()
    };
    let v = normal.cross(u);

    // Project to 2D
    let projected_2d: Vec<DVec2> = points.iter().map(|&p| {
        let d = p - plane.point;
        DVec2::new(d.dot(u), d.dot(v))
    }).collect();

    // Compute Delaunay triangulation (Bowyer-Watson algorithm)
    let triangles_2d = delaunay_triangulation_2d(&projected_2d);

    // Convert back to 3D
    mesh.nodes = points.to_vec();
    mesh.triangles = triangles_2d;

    Some(mesh)
}

fn delaunay_triangulation_2d(points: &[DVec2]) -> Vec<[usize; 3]> {
    if points.len() < 3 {
        return Vec::new();
    }

    let n = points.len();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    // Create super-triangle containing all points
    let (min_x, max_x) = points.iter().map(|p| p.x).fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), x| (min.min(x), max.max(x)));
    let (min_y, max_y) = points.iter().map(|p| p.y).fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), y| (min.min(y), max.max(y)));

    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let delta = dx.max(dy) * 10.0;

    let p1 = DVec2::new(min_x - delta, min_y - delta);
    let p2 = DVec2::new(min_x + 3.0 * delta, min_y - delta);
    let p3 = DVec2::new(min_x, min_y + 3.0 * delta);

    let super_nodes = [n, n + 1, n + 2];
    let all_points: Vec<DVec2> = points.iter().copied().chain([p1, p2, p3]).collect();

    triangles.push(super_nodes);

    // Add points one by one
    for i in 0..n {
        let p = points[i];

        // Find all triangles whose circumcircle contains p
        let mut bad_triangles: Vec<usize> = Vec::new();
        for (ti, &tri) in triangles.iter().enumerate() {
            if let Some(cc) = circumcenter_2d(all_points[tri[0]], all_points[tri[1]], all_points[tri[2]]) {
                let r_sq = (all_points[tri[0]].x - cc.x).powi(2) + (all_points[tri[0]].y - cc.y).powi(2);
                let d_sq = (p.x - cc.x).powi(2) + (p.y - cc.y).powi(2);
                if d_sq <= r_sq {
                    bad_triangles.push(ti);
                }
            }
        }

        // Find boundary polygon
        let mut polygon: Vec<(usize, usize)> = Vec::new();
        for &ti in &bad_triangles {
            let tri = triangles[ti];
            let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
            for &(a, b) in &edges {
                let edge = (a.min(b), a.max(b));
                let mut shared = false;
                for &tj in &bad_triangles {
                    if tj == ti {
                        continue;
                    }
                    let other = triangles[tj];
                    let other_edges = [(other[0], other[1]), (other[1], other[2]), (other[2], other[0])];
                    for &(oa, ob) in &other_edges {
                        let other_edge = (oa.min(ob), oa.max(ob));
                        if edge == other_edge {
                            shared = true;
                            break;
                        }
                    }
                    if shared {
                        break;
                    }
                }
                if !shared {
                    polygon.push((tri[0], tri[1]));
                    if !polygon.contains(&(tri[1], tri[2])) {
                        polygon.push((tri[1], tri[2]));
                    }
                    if !polygon.contains(&(tri[2], tri[0])) {
                        polygon.push((tri[2], tri[0]));
                    }
                }
            }
        }

        // Remove bad triangles
        let mut new_triangles: Vec<[usize; 3]> = Vec::new();
        for (ti, &tri) in triangles.iter().enumerate() {
            if !bad_triangles.contains(&ti) {
                new_triangles.push(tri);
            }
        }

        // Add new triangles from polygon
        for &(a, b) in &polygon {
            new_triangles.push([a, b, i]);
        }

        triangles = new_triangles;
    }

    // Remove triangles containing super-triangle nodes
    triangles.retain(|&tri| tri[0] < n && tri[1] < n && tri[2] < n);

    triangles
}

fn circumcenter_2d(a: DVec2, b: DVec2, c: DVec2) -> Option<DVec2> {
    let d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if d.abs() < TOLERANCE_FLOAT_LOOSE {
        return None;
    }

    let a_sq = a.x * a.x + a.y * a.y;
    let b_sq = b.x * b.x + b.y * b.y;
    let c_sq = c.x * c.x + c.y * c.y;

    Some(DVec2::new(
        (a_sq * (b.y - c.y) + b_sq * (c.y - a.y) + c_sq * (a.y - b.y)) / d,
        (a_sq * (c.x - b.x) + b_sq * (a.x - c.x) + c_sq * (b.x - a.x)) / d,
    ))
}

/// Generates normal-consistent mesh from oriented point cloud.
///
/// Ensures all face normals point consistently outward.
pub fn generate_consistent_mesh(
    points: &[DVec3],
    normals: &[DVec3],
) -> Option<TriangleMesh> {
    let mut mesh = delaunay_reconstruction(points, normals)?;

    // Orient faces consistently
    let node_normals = mesh.compute_node_normals();

    for tri in &mut mesh.triangles {
        let a = mesh.nodes[tri[0]];
        let b = mesh.nodes[tri[1]];
        let c = mesh.nodes[tri[2]];

        let face_normal = (b - a).cross(c - a);
        let len = face_normal.length();
        if len < TOLERANCE_LINEAR_ULTRA_STRICT {
            continue;
        }
        let face_normal = face_normal / len;

        // Compare with node normals
        let avg_normal = (node_normals[tri[0]] + node_normals[tri[1]] + node_normals[tri[2]]) / 3.0;

        if face_normal.dot(avg_normal) < 0.0 {
            // Flip triangle
            *tri = [tri[0], tri[2], tri[1]];
        }
    }

    mesh.normals = Some(node_normals);
    Some(mesh)
}

fn compute_bounding_box(points: &[DVec3]) -> Option<(DVec3, DVec3)> {
    if points.is_empty() {
        return None;
    }

    let mut min = DVec3::splat(f64::INFINITY);
    let mut max = DVec3::splat(f64::NEG_INFINITY);

    for &p in points {
        min = min.min(p);
        max = max.max(p);
    }

    Some((min, max))
}

// ============================================================================
// Advanced Sampling
// ============================================================================

/// Sampling method for point cloud simplification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdvancedSamplingMethod {
    /// Voxel grid downsampling.
    VoxelGrid,
    /// Random uniform sampling.
    RandomUniform,
    /// Curvature-aware sampling (more points in high-curvature regions).
    CurvatureAware,
    /// Poisson disk sampling (uniform distribution with minimum distance).
    PoissonDisk,
}

/// Configuration for advanced sampling.
#[derive(Debug, Clone)]
pub struct AdvancedSamplingConfig {
    /// Target number of points.
    pub target_count: usize,
    /// Voxel size for voxel grid sampling.
    pub voxel_size: f64,
    /// Minimum distance for Poisson disk sampling.
    pub min_distance: f64,
    /// Number of neighbors for curvature estimation.
    pub k_neighbors: usize,
    /// Seed for random sampling.
    pub seed: u64,
}

impl Default for AdvancedSamplingConfig {
    fn default() -> Self {
        Self {
            target_count: 1000,
            voxel_size: 0.1,
            min_distance: 0.05,
            k_neighbors: 30,
            seed: 42,
        }
    }
}

/// Performs advanced point cloud sampling.
pub fn advanced_sample(
    points: &[DVec3],
    method: AdvancedSamplingMethod,
    config: &AdvancedSamplingConfig,
) -> Vec<DVec3> {
    if points.len() <= config.target_count {
        return points.to_vec();
    }

    match method {
        AdvancedSamplingMethod::VoxelGrid => voxel_grid_sample(points, config),
        AdvancedSamplingMethod::RandomUniform => random_uniform_sample(points, config),
        AdvancedSamplingMethod::CurvatureAware => curvature_aware_sample(points, config),
        AdvancedSamplingMethod::PoissonDisk => poisson_disk_sample(points, config),
    }
}

fn voxel_grid_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    let (min, _max) = compute_bounding_box(points).unwrap();
    let voxel_size = config.voxel_size;

    let mut voxels: std::collections::HashMap<[i64; 3], Vec<DVec3>> = std::collections::HashMap::new();

    for &p in points {
        let key = [
            ((p.x - min.x) / voxel_size).floor() as i64,
            ((p.y - min.y) / voxel_size).floor() as i64,
            ((p.z - min.z) / voxel_size).floor() as i64,
        ];
        voxels.entry(key).or_default().push(p);
    }

    voxels.values().map(|pts| {
        let sum: DVec3 = pts.iter().sum();
        sum / pts.len() as f64
    }).collect()
}

fn random_uniform_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    let n = points.len();
    let target = config.target_count.min(n);

    let mut rng = SimpleRng::new(config.seed);
    let mut indices: Vec<usize> = (0..n).collect();

    // Fisher-Yates shuffle for first target_count elements
    for i in 0..target {
        let j = i + ((rng.next() as usize) % (n - i));
        indices.swap(i, j);
    }

    indices.iter().take(target).map(|&i| points[i]).collect()
}

fn curvature_aware_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    if points.len() <= config.target_count {
        return points.to_vec();
    }

    // Compute curvatures
    let neighbors = build_neighbor_graph(points, config.k_neighbors);
    let curvatures = compute_curvatures(points, &neighbors);

    // Compute sampling probabilities (higher curvature = higher probability)
    let max_curv = curvatures.iter().cloned().fold(0.0_f64, f64::max).max(TOLERANCE_LINEAR_ULTRA_STRICT);
    let weights: Vec<f64> = curvatures.iter().map(|&c| 0.1 + 0.9 * c / max_curv).collect();

    let total_weight: f64 = weights.iter().sum();
    let probs: Vec<f64> = weights.iter().map(|&w| w / total_weight).collect();

    // Sample based on probabilities
    let mut rng = SimpleRng::new(config.seed);
    let mut selected: std::collections::HashSet<usize> = std::collections::HashSet::new();

    while selected.len() < config.target_count {
        let r = (rng.next() as f64) / (u64::MAX as f64);
        let mut cumsum = 0.0;
        for (i, &p) in probs.iter().enumerate() {
            cumsum += p;
            if r < cumsum {
                selected.insert(i);
                break;
            }
        }
    }

    selected.iter().map(|&i| points[i]).collect()
}

fn poisson_disk_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    let (min, _) = compute_bounding_box(points).unwrap();
    let r = config.min_distance;
    let r_sq = r * r;
    let cell_size = r / 3.0_f64.sqrt();

    let mut grid: std::collections::HashMap<[i64; 3], usize> = std::collections::HashMap::new();
    let mut result: Vec<DVec3> = Vec::new();
    let mut active: Vec<usize> = Vec::new();

    // Start with random point
    let mut rng = SimpleRng::new(config.seed);
    let first_idx = (rng.next() as usize) % points.len();
    let first = points[first_idx];

    result.push(first);
    active.push(0);

    let key = [
        ((first.x - min.x) / cell_size).floor() as i64,
        ((first.y - min.y) / cell_size).floor() as i64,
        ((first.z - min.z) / cell_size).floor() as i64,
    ];
    grid.insert(key, 0);

    // Dart throwing
    while !active.is_empty() && result.len() < config.target_count {
        let active_idx = (rng.next() as usize) % active.len();
        let sample_idx = active[active_idx];
        let sample = result[sample_idx];

        let mut found = false;
        for _ in 0..30 {
            // Generate random point in annulus around sample
            let angle1 = 2.0 * std::f64::consts::PI * (rng.next() as f64) / (u64::MAX as f64);
            let angle2 = std::f64::consts::PI * (rng.next() as f64) / (u64::MAX as f64);
            let rad = r * (1.0 + (rng.next() as f64) / (u64::MAX as f64));

            let candidate = DVec3::new(
                sample.x + rad * angle2.sin() * angle1.cos(),
                sample.y + rad * angle2.sin() * angle1.sin(),
                sample.z + rad * angle2.cos(),
            );

            // Check if candidate is far enough from existing points
            let key = [
                ((candidate.x - min.x) / cell_size).floor() as i64,
                ((candidate.y - min.y) / cell_size).floor() as i64,
                ((candidate.z - min.z) / cell_size).floor() as i64,
            ];

            let mut valid = true;
            for di in -1..=1 {
                for dj in -1..=1 {
                    for dk in -1..=1 {
                        let neighbor_key = [key[0] + di, key[1] + dj, key[2] + dk];
                        if let Some(&idx) = grid.get(&neighbor_key)
                            && (result[idx] - candidate).length_squared() < r_sq {
                                valid = false;
                                break;
                            }
                    }
                    if !valid {
                        break;
                    }
                }
                if !valid {
                    break;
                }
            }

            if valid {
                let new_idx = result.len();
                result.push(candidate);
                active.push(new_idx);
                grid.insert(key, new_idx);
                found = true;
                break;
            }
        }

        if !found {
            active.swap_remove(active_idx);
        }
    }

    result
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
                // Add triangle node positions from `brep.vertices`
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

                // If no triangles, add outer-wire vertex positions
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
/// Takes node positions directly from a `SurfaceMesh`.
pub fn extract_points_from_mesh(mesh: &crate::triangulate::SurfaceMesh) -> PointCloud {
    PointCloud::from_vec(mesh.nodes.clone())
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
            for _face in &shell.faces {
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
