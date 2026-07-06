//! Medial axis extraction and wall thickness analysis.
//!
//! This module provides algorithms for computing the medial axis (also known as
//! the skeleton) of 2D profiles and 3D solids. The medial axis represents the
//! set of points that have at least two closest points on the boundary.
//!
//! Applications include:
//! - Wall thickness analysis for injection molding and casting
//! - Detection of thin regions that may cause manufacturing defects
//! - Generation of rib/stiffener paths for structural reinforcement
//! - Shape simplification and feature recognition
//! - Mid-surface extraction for FEA shell meshing
//!
//! # OCCT Equivalents
//!
//! This module provides functionality similar to:
//! - `GeomAPI_PointsToBSpline` for medial curve approximation
//! - `BRepExtrema_DistShapeShape` for distance computations
//! - `BRepAdaptor_Surface` for surface analysis
//!
//! # Examples
//!
//! ```
//! use rcad_algorithms::medial_axis::{compute_medial_axis_2d, MedialAxisOptions};
//! use glam::dvec3;
//!
//! let polygon = vec![
//!     dvec3(0.0, 0.0, 0.0),
//!     dvec3(2.0, 0.0, 0.0),
//!     dvec3(2.0, 1.0, 0.0),
//!     dvec3(1.0, 1.0, 0.0),
//!     dvec3(1.0, 2.0, 0.0),
//!     dvec3(0.0, 2.0, 0.0),
//! ];
//! let opts = MedialAxisOptions::default();
//! let axis = compute_medial_axis_2d(&polygon, &opts);
//! println!("Found {} medial points", axis.all_points.len());
//! ```

use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::{topods, BRep, Curve3, Surface3, Face, Shell, Solid, SurfaceEval, CurveEval, Wire, WireEdge};
use rcad_kernel::geom::{Line3, Plane};
use std::collections::{HashMap, HashSet};

/// Options for medial axis computation.
#[derive(Debug, Clone)]
pub struct MedialAxisOptions {
    /// Geometric tolerance for numerical operations.
    pub tolerance: f64,
    /// Minimum thickness threshold for filtering medial points.
    pub min_thickness: f64,
    /// Whether to simplify the result by removing close points.
    pub simplify: bool,
    /// Number of samples per direction for surface sampling.
    pub sample_density: usize,
    /// Maximum recursion depth for Voronoi subdivision.
    pub voronoi_depth: usize,
    /// Angle tolerance for detecting sharp corners (radians).
    pub corner_angle_tol: f64,
    /// Maximum distance for clustering medial points (3D).
    pub cluster_distance: f64,
    /// Number of refinement iterations for medial surface extraction.
    pub refinement_iterations: usize,
    /// Enable chordal axis transform for better thin feature detection.
    pub use_chordal_axis: bool,
    /// Minimum feature size to detect (for thin region analysis).
    pub min_feature_size: f64,
    /// Angular resolution for ray casting (3D distance field).
    pub angular_resolution: f64,
}

impl Default for MedialAxisOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_MESH_LEGACY,
            min_thickness: 0.001,
            simplify: true,
            sample_density: 100,
            voronoi_depth: 10,
            corner_angle_tol: 0.1,
            cluster_distance: 0.01,
            refinement_iterations: 3,
            use_chordal_axis: true,
            min_feature_size: 0.01,
            angular_resolution: std::f64::consts::PI / 36.0, // 5 degrees
        }
    }
}

/// Options for mid-surface extraction.
#[derive(Debug, Clone)]
pub struct MidSurfaceOptions {
    /// Base computation options.
    pub base: MedialAxisOptions,
    /// Maximum thickness ratio for treating as thin-walled.
    pub max_thickness_ratio: f64,
    /// Minimum aspect ratio for thin wall detection.
    pub min_aspect_ratio: f64,
    /// Target surface continuity.
    pub continuity: ContinuityLevel,
    /// Whether to preserve sharp features.
    pub preserve_features: bool,
    /// Feature angle threshold (radians).
    pub feature_angle: f64,
}

impl Default for MidSurfaceOptions {
    fn default() -> Self {
        Self {
            base: MedialAxisOptions::default(),
            max_thickness_ratio: 0.1,
            min_aspect_ratio: 10.0,
            continuity: ContinuityLevel::C0,
            preserve_features: true,
            feature_angle: std::f64::consts::PI / 6.0, // 30 degrees
        }
    }
}

/// Surface continuity levels for mid-surface extraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuityLevel {
    /// Position continuity only.
    C0,
    /// Tangent continuity.
    C1,
    /// Curvature continuity.
    C2,
}

/// Options for rib/stiffener generation.
#[derive(Debug, Clone)]
pub struct RibGenerationOptions {
    /// Base computation options.
    pub base: MedialAxisOptions,
    /// Minimum rib height.
    pub min_height: f64,
    /// Maximum rib height.
    pub max_height: f64,
    /// Rib draft angle (radians).
    pub draft_angle: f64,
    /// Minimum rib length.
    pub min_length: f64,
    /// Spacing between parallel ribs.
    pub spacing: f64,
    /// Whether to optimize for structural stiffness.
    pub optimize_stiffness: bool,
    /// Weight for thickness uniformity in optimization.
    pub thickness_weight: f64,
}

impl Default for RibGenerationOptions {
    fn default() -> Self {
        Self {
            base: MedialAxisOptions::default(),
            min_height: 2.0,
            max_height: 20.0,
            draft_angle: std::f64::consts::PI / 36.0, // 5 degrees
            min_length: 10.0,
            spacing: 20.0,
            optimize_stiffness: true,
            thickness_weight: 0.5,
        }
    }
}

// ============================================================================
// 2D Data Structures
// ============================================================================

/// A point on a 2D medial axis with associated radius.
#[derive(Debug, Clone, Copy)]
pub struct MedialPoint2d {
    /// Position in 2D space.
    pub point: DVec2,
    /// Radius of the maximal inscribed circle at this point.
    pub radius: f64,
    /// Whether this is a branch point (3+ touching boundary points).
    pub is_branch: bool,
    /// Whether this is an end point (touches a convex vertex).
    pub is_end: bool,
}

/// A branch of the 2D medial axis.
///
/// Each branch traces the locus of centers of inscribed circles
/// that touch the boundary at exactly two points.
#[derive(Debug, Clone)]
pub struct MedialBranch2d {
    /// Points along the branch in order.
    pub points: Vec<MedialPoint2d>,
    /// Index of the parent branch (-1 for root branches).
    pub parent: Option<usize>,
    /// Indices of child branches.
    pub children: Vec<usize>,
    /// Source edge indices on the original polygon.
    pub source_edges: (usize, usize),
}

/// The complete 2D medial axis (skeleton) of a polygon.
#[derive(Debug, Clone, Default)]
pub struct MedialAxis2d {
    /// All branches of the medial axis.
    pub branches: Vec<MedialBranch2d>,
    /// All unique points across all branches.
    pub all_points: Vec<MedialPoint2d>,
    /// Branch points where 3+ branches meet.
    pub branch_points: Vec<usize>,
    /// End points (leaves) of the medial axis.
    pub end_points: Vec<usize>,
    /// Maximum inscribed circle info.
    pub max_inscribed_circle: Option<(DVec2, f64)>,
}

/// A Voronoi vertex in 2D.
#[derive(Debug, Clone)]
pub struct VoronoiVertex2d {
    /// Position of the vertex.
    pub point: DVec2,
    /// Index of the input site this vertex is equidistant to.
    pub sites: Vec<usize>,
}

/// A Voronoi edge in 2D.
#[derive(Debug, Clone)]
pub struct VoronoiEdge2d {
    /// Start vertex index (or None for unbounded).
    pub start: Option<usize>,
    /// End vertex index (or None for unbounded).
    pub end: Option<usize>,
    /// The two sites this edge bisects.
    pub sites: (usize, usize),
    /// Whether this is a finite edge.
    pub is_finite: bool,
}

/// A Voronoi diagram in 2D.
#[derive(Debug, Clone, Default)]
pub struct VoronoiDiagram2d {
    /// Input sites (points).
    pub sites: Vec<DVec2>,
    /// Voronoi vertices.
    pub vertices: Vec<VoronoiVertex2d>,
    /// Voronoi edges.
    pub edges: Vec<VoronoiEdge2d>,
    /// Cells: for each site, the indices of its cell edges.
    pub cells: Vec<Vec<usize>>,
}

// ============================================================================
// 3D Data Structures
// ============================================================================

/// A vertex on the medial axis/surface.
///
/// Each vertex represents a point where the inscribed sphere
/// touches the boundary at two or more points.
#[derive(Debug, Clone)]
pub struct MedialVertex {
    /// Position of the medial vertex.
    pub point: DVec3,
    /// Radius of the inscribed sphere at this point.
    pub radius: f64,
    /// Indices of the boundary elements this point is closest to.
    pub boundary_elements: Vec<usize>,
}

/// An edge on the medial axis.
///
/// Edges connect vertices and trace the locus of centers of
/// inscribed spheres that touch the boundary at two points.
#[derive(Debug, Clone)]
pub struct MedialEdge {
    /// Index of the start vertex.
    pub start_vertex: usize,
    /// Index of the end vertex.
    pub end_vertex: usize,
    /// The curve geometry (if representable).
    pub curve: Option<Curve3>,
    /// Radius at the start of the edge.
    pub start_radius: f64,
    /// Radius at the end of the edge.
    pub end_radius: f64,
}

/// A face on the medial axis (3D case).
///
/// Faces represent regions where the inscribed sphere touches
/// the boundary at three or more points.
#[derive(Debug, Clone)]
pub struct MedialFace {
    /// Indices of the vertices forming the face boundary.
    pub vertices: Vec<usize>,
    /// The surface geometry (if representable).
    pub surface: Option<Surface3>,
    /// Minimum inscribed radius within this face.
    pub min_radius: f64,
    /// Maximum inscribed radius within this face.
    pub max_radius: f64,
}

/// The computed medial axis/surface for 3D solids.
#[derive(Debug, Clone, Default)]
pub struct MedialSurface {
    /// Medial vertices.
    pub vertices: Vec<MedialVertex>,
    /// Medial edges (centerlines).
    pub edges: Vec<MedialEdge>,
    /// Medial faces (surface patches).
    pub faces: Vec<MedialFace>,
    /// Overall thickness statistics.
    pub thickness_stats: ThicknessStats,
}

/// Statistics about thickness distribution.
#[derive(Debug, Clone, Copy, Default)]
pub struct ThicknessStats {
    /// Minimum thickness found.
    pub min: f64,
    /// Maximum thickness found.
    pub max: f64,
    /// Mean thickness.
    pub mean: f64,
    /// Standard deviation.
    pub std_dev: f64,
}

/// A detected thin-walled region.
#[derive(Debug, Clone)]
pub struct ThinRegion {
    /// Center point of the thin region.
    pub center: DVec3,
    /// Thickness at this region.
    pub thickness: f64,
    /// Approximate area affected.
    pub area: f64,
    /// Indices of associated faces in the original model.
    pub face_indices: Vec<usize>,
    /// Severity level (0-1, where 1 is critically thin).
    pub severity: f64,
}

/// Result of wall thickness analysis.
#[derive(Debug, Clone)]
pub struct ThicknessMap {
    /// Thickness values at sample points.
    pub samples: Vec<ThicknessSample>,
    /// Overall statistics.
    pub stats: ThicknessStats,
    /// Detected thin regions.
    pub thin_regions: Vec<ThinRegion>,
}

/// A sample point in the thickness map.
#[derive(Debug, Clone, Copy)]
pub struct ThicknessSample {
    /// Sample point position.
    pub point: DVec3,
    /// Thickness at this point (2x the medial radius).
    pub thickness: f64,
    /// Direction to nearest boundary point.
    pub normal: DVec3,
    /// Index of the nearest face.
    pub nearest_face: usize,
}

/// Result of wall thickness analysis.
#[derive(Debug, Clone)]
pub struct WallThicknessResult {
    /// Minimum wall thickness found.
    pub min_thickness: f64,
    /// Maximum wall thickness found.
    pub max_thickness: f64,
    /// Average wall thickness.
    pub avg_thickness: f64,
    /// Detected thin regions below threshold.
    pub thin_regions: Vec<ThinRegion>,
}

/// Result of mid-surface extraction for FEA shell meshing.
#[derive(Debug, Clone)]
pub struct MidSurfaceResult {
    /// The extracted mid-surface as a B-Rep.
    pub brep: BRep,
    /// Thickness at each face.
    pub face_thickness: Vec<f64>,
    /// Mapping from mid-surface face to original solid faces.
    pub face_mapping: Vec<(usize, usize)>,
}

// ============================================================================
// Enhanced 3D Data Structures
// ============================================================================

/// A chordal axis vertex for thin feature detection.
///
/// The chordal axis is a simplified version of the medial axis that
/// focuses on the centerlines of thin-walled regions.
#[derive(Debug, Clone)]
pub struct ChordalVertex {
    /// Position of the vertex.
    pub point: DVec3,
    /// Local thickness at this point.
    pub thickness: f64,
    /// Principal direction of the thin feature.
    pub direction: DVec3,
    /// Normal to the mid-surface.
    pub normal: DVec3,
    /// Associated boundary face pairs.
    pub face_pairs: Vec<(usize, usize)>,
}

/// A chordal axis edge connecting two vertices.
#[derive(Debug, Clone)]
pub struct ChordalEdge {
    /// Start vertex index.
    pub start: usize,
    /// End vertex index.
    pub end: usize,
    /// Approximate curve geometry.
    pub curve: Option<Curve3>,
    /// Average thickness along this edge.
    pub avg_thickness: f64,
    /// Length of the edge.
    pub length: f64,
}

/// The chordal axis of a thin-walled solid.
#[derive(Debug, Clone, Default)]
pub struct ChordalAxis {
    /// Vertices of the chordal axis.
    pub vertices: Vec<ChordalVertex>,
    /// Edges connecting vertices.
    pub edges: Vec<ChordalEdge>,
    /// Identified thin sheets.
    pub sheets: Vec<ThinSheet>,
}

/// A thin sheet region in the solid.
#[derive(Debug, Clone)]
pub struct ThinSheet {
    /// Index of the chordal edge forming the sheet spine.
    pub spine_edge: usize,
    /// Face indices on one side of the sheet.
    pub side_a_faces: Vec<usize>,
    /// Face indices on the other side.
    pub side_b_faces: Vec<usize>,
    /// Average thickness of the sheet.
    pub avg_thickness: f64,
    /// Area of the sheet region.
    pub area: f64,
    /// Quality of the thin sheet (0-1).
    pub quality: f64,
}

/// Classification of wall thickness regions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThicknessClass {
    /// Very thin region (< 25% of target).
    VeryThin,
    /// Thin region (25-50% of target).
    Thin,
    /// Normal thickness (50-150% of target).
    Normal,
    /// Thick region (150-200% of target).
    Thick,
    /// Very thick region (> 200% of target).
    VeryThick,
}

/// Detailed thin region analysis result.
#[derive(Debug, Clone)]
pub struct ThinRegionAnalysis {
    /// All detected thin regions.
    pub regions: Vec<ThinRegion>,
    /// Overall classification of wall thickness.
    pub classification: ThicknessClass,
    /// Recommended minimum wall thickness.
    pub recommended_min: f64,
    /// Regions grouped by severity.
    pub severity_groups: HashMap<ThinRegionSeverity, Vec<usize>>,
    /// Histogram of thickness values.
    pub thickness_histogram: Vec<ThicknessHistogramBin>,
}

/// Severity level for thin regions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ThinRegionSeverity {
    /// Critical: immediate manufacturing risk.
    Critical,
    /// Warning: may cause issues.
    Warning,
    /// Acceptable: within tolerance but notable.
    Acceptable,
}

/// A bin in the thickness histogram.
#[derive(Debug, Clone, Copy)]
pub struct ThicknessHistogramBin {
    /// Lower bound of thickness values in this bin.
    pub lower: f64,
    /// Upper bound of thickness values in this bin.
    pub upper: f64,
    /// Number of samples in this bin.
    pub count: usize,
}

/// A rib/stiffener placement recommendation.
#[derive(Debug, Clone)]
pub struct RibPlacement {
    /// Centerline curve of the rib.
    pub centerline: Curve3,
    /// Start point of the rib.
    pub start: DVec3,
    /// End point of the rib.
    pub end: DVec3,
    /// Recommended height.
    pub height: f64,
    /// Recommended width (at base).
    pub width: f64,
    /// Draft angle.
    pub draft_angle: f64,
    /// Structural efficiency score (0-1).
    pub efficiency: f64,
    /// Associated medial axis edge.
    pub medial_edge: Option<usize>,
    /// Index of the face the rib attaches to.
    pub attached_face: usize,
}

/// Result of rib/stiffener generation.
#[derive(Debug, Clone)]
pub struct RibGenerationResult {
    /// Generated rib placements.
    pub ribs: Vec<RibPlacement>,
    /// Total rib volume added.
    pub total_volume: f64,
    /// Estimated stiffness improvement.
    pub stiffness_improvement: f64,
    /// Weight increase percentage.
    pub weight_increase: f64,
    /// Quality of the rib layout.
    pub quality_score: f64,
}

/// An octree node for distance field computation.
#[derive(Debug, Clone)]
struct OctreeNode {
    /// Bounding box minimum.
    min: DVec3,
    /// Bounding box maximum.
    max: DVec3,
    /// Distance value at the center.
    distance: f64,
    /// Children (8 for internal nodes, 0 for leaves).
    children: Vec<OctreeNode>,
    /// Whether this is a medial point (local maximum).
    is_medial: bool,
    /// Depth in the octree.
    depth: usize,
}

/// A voxel grid for distance field representation.
#[derive(Debug, Clone)]
pub struct VoxelGrid {
    /// Origin of the grid.
    pub origin: DVec3,
    /// Size of each voxel.
    pub voxel_size: f64,
    /// Number of voxels in each dimension.
    pub dimensions: [usize; 3],
    /// Distance values at each voxel.
    pub distances: Vec<f64>,
    /// Gradient vectors at each voxel.
    pub gradients: Vec<DVec3>,
    /// Whether each voxel is inside the solid.
    pub inside: Vec<bool>,
}

impl VoxelGrid {
    /// Create a new voxel grid.
    pub fn new(origin: DVec3, voxel_size: f64, dimensions: [usize; 3]) -> Self {
        let total = dimensions[0] * dimensions[1] * dimensions[2];
        Self {
            origin,
            voxel_size,
            dimensions,
            distances: vec![0.0; total],
            gradients: vec![DVec3::ZERO; total],
            inside: vec![false; total],
        }
    }

    /// Get the index for a voxel position.
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + j * self.dimensions[0] + k * self.dimensions[0] * self.dimensions[1]
    }

    /// Get the world position of a voxel center.
    pub fn voxel_center(&self, i: usize, j: usize, k: usize) -> DVec3 {
        self.origin + DVec3::new(
            (i as f64 + 0.5) * self.voxel_size,
            (j as f64 + 0.5) * self.voxel_size,
            (k as f64 + 0.5) * self.voxel_size,
        )
    }

    /// Get the distance at a voxel.
    pub fn get_distance(&self, i: usize, j: usize, k: usize) -> f64 {
        self.distances[self.index(i, j, k)]
    }

    /// Set the distance at a voxel.
    pub fn set_distance(&mut self, i: usize, j: usize, k: usize, d: f64) {
        let idx = self.index(i, j, k);
        self.distances[idx] = d;
    }

    /// Check if a voxel is inside the solid.
    pub fn is_inside(&self, i: usize, j: usize, k: usize) -> bool {
        self.inside[self.index(i, j, k)]
    }

    /// Find local maxima in the distance field (medial axis candidates).
    pub fn find_local_maxima(&self, threshold: f64) -> Vec<(usize, usize, usize, f64)> {
        let mut maxima = Vec::new();

        for k in 1..self.dimensions[2] - 1 {
            for j in 1..self.dimensions[1] - 1 {
                for i in 1..self.dimensions[0] - 1 {
                    if !self.is_inside(i, j, k) {
                        continue;
                    }

                    let d = self.get_distance(i, j, k);
                    if d < threshold {
                        continue;
                    }

                    // Check if this is a local maximum
                    let mut is_max = true;
                    for di in -1..=1 {
                        for dj in -1..=1 {
                            for dk in -1..=1 {
                                if di == 0 && dj == 0 && dk == 0 {
                                    continue;
                                }
                                let ni = (i as isize + di) as usize;
                                let nj = (j as isize + dj) as usize;
                                let nk = (k as isize + dk) as usize;
                                if self.get_distance(ni, nj, nk) > d {
                                    is_max = false;
                                    break;
                                }
                            }
                            if !is_max {
                                break;
                            }
                        }
                        if !is_max {
                            break;
                        }
                    }

                    if is_max {
                        maxima.push((i, j, k, d));
                    }
                }
            }
        }

        maxima
    }
}

/// Mid-surface extraction with enhanced geometry.
#[derive(Debug, Clone)]
pub struct EnhancedMidSurfaceResult {
    /// The extracted mid-surface as a B-Rep.
    pub brep: BRep,
    /// Thickness at each face.
    pub face_thickness: Vec<f64>,
    /// Mapping from mid-surface face to original solid faces.
    pub face_mapping: Vec<(usize, usize)>,
    /// Chordal axis of the thin-walled solid.
    pub chordal_axis: ChordalAxis,
    /// Quality metrics for the extraction.
    pub quality: MidSurfaceQuality,
}

/// Quality metrics for mid-surface extraction.
#[derive(Debug, Clone, Copy, Default)]
pub struct MidSurfaceQuality {
    /// Percentage of the solid successfully represented.
    pub coverage: f64,
    /// Average deviation from true mid-surface.
    pub avg_deviation: f64,
    /// Maximum deviation from true mid-surface.
    pub max_deviation: f64,
    /// Thickness accuracy (correlation coefficient).
    pub thickness_accuracy: f64,
    /// Number of discontinuities in the mid-surface.
    pub discontinuities: usize,
    /// Overall quality score (0-1).
    pub overall_score: f64,
}

// ============================================================================
// 2D Medial Axis Computation
// ============================================================================

/// Compute the medial axis of a 2D polygon.
///
/// Uses a Voronoi-based approach:
/// 1. Sample points on the polygon boundary
/// 2. Compute constrained Voronoi diagram
/// 3. Extract internal Voronoi edges as the medial axis
///
/// # Arguments
/// * `polygon` - Ordered boundary points (closed polygon, Z-coordinate ignored)
/// * `opts` - Computation options
///
/// # Returns
/// The computed 2D medial axis.
pub fn compute_medial_axis_2d(polygon: &[DVec3], opts: &MedialAxisOptions) -> MedialAxis2d {
    let n = polygon.len();
    if n < 3 {
        return MedialAxis2d::default();
    }

    // Convert to 2D points
    let points2d: Vec<DVec2> = polygon
        .iter()
        .map(|p| DVec2::new(p.x, p.y))
        .collect();

    compute_medial_axis_2d_from_points(&points2d, opts)
}

/// Compute the medial axis from 2D points.
pub fn compute_medial_axis_2d_from_points(points: &[DVec2], opts: &MedialAxisOptions) -> MedialAxis2d {
    let n = points.len();
    if n < 3 {
        return MedialAxis2d::default();
    }

    // Step 1: Sample points on the polygon edges
    let sampled_points = sample_polygon_boundary(points, opts);

    // Step 2: Compute Voronoi diagram
    let voronoi = compute_voronoi_2d(&sampled_points, opts);

    // Step 3: Extract medial axis as internal Voronoi edges
    extract_medial_axis_from_voronoi(&voronoi, points, opts)
}

/// Sample points densely on the polygon boundary.
fn sample_polygon_boundary(polygon: &[DVec2], opts: &MedialAxisOptions) -> Vec<DVec2> {
    let n = polygon.len();
    if n == 0 {
        return vec![];
    }

    let mut samples = Vec::new();

    for i in 0..n {
        let p0 = polygon[i];
        let p1 = polygon[(i + 1) % n];
        let edge_len = (p1 - p0).length();

        // Sample based on edge length and tolerance
        let num_samples = (edge_len / opts.tolerance).max(2.0).ceil() as usize;

        for j in 0..num_samples {
            let t = j as f64 / num_samples as f64;
            samples.push(p0 + t * (p1 - p0));
        }
    }

    samples
}

/// Compute a 2D Voronoi diagram using a simple incremental approach.
///
/// For robustness, this uses the Bowyer-Watson algorithm for Delaunay
/// triangulation, then extracts the dual Voronoi graph.
pub fn compute_voronoi_2d(sites: &[DVec2], opts: &MedialAxisOptions) -> VoronoiDiagram2d {
    let n = sites.len();
    if n < 2 {
        return VoronoiDiagram2d {
            sites: sites.to_vec(),
            vertices: vec![],
            edges: vec![],
            cells: vec![],
        };
    }

    // Compute Delaunay triangulation
    let triangles = compute_delaunay_2d(sites, opts);

    // Extract Voronoi nodes and edges from Delaunay triangles
    let mut voronoi_nodes: Vec<VoronoiVertex2d> = Vec::new();
    let mut edges: Vec<VoronoiEdge2d> = Vec::new();
    let mut cells: Vec<Vec<usize>> = vec![vec![]; n];

    // Map from edge (sorted pair of sites) to Voronoi edge
    let mut edge_map: HashMap<(usize, usize), (usize, Option<usize>)> = HashMap::new();

    // Each Delaunay triangle gives one Voronoi node (circumcenter)
    for tri in &triangles {
        let p0 = sites[tri[0]];
        let p1 = sites[tri[1]];
        let p2 = sites[tri[2]];

        // Compute circumcenter
        if let Some((center, _radius)) = circumcenter(p0, p1, p2) {
            let v_idx = voronoi_nodes.len();
            voronoi_nodes.push(VoronoiVertex2d {
                point: center,
                sites: tri.to_vec(),
            });

            // Create Voronoi edges for each triangle edge
            for k in 0..3 {
                let i1 = tri[k];
                let i2 = tri[(k + 1) % 3];
                let key = if i1 < i2 { (i1, i2) } else { (i2, i1) };

                if let Some((prev_v, _prev_t)) = edge_map.get(&key).copied() {
                    // Connect the two Voronoi nodes
                    edges.push(VoronoiEdge2d {
                        start: Some(prev_v),
                        end: Some(v_idx),
                        sites: key,
                        is_finite: true,
                    });
                    let e_idx = edges.len() - 1;
                    cells[i1].push(e_idx);
                    cells[i2].push(e_idx);
                } else {
                    edge_map.insert(key, (v_idx, None));
                }
            }
        }
    }

    VoronoiDiagram2d {
        sites: sites.to_vec(),
        vertices: voronoi_nodes,
        edges,
        cells,
    }
}

/// Compute Delaunay triangulation using Bowyer-Watson algorithm.
fn compute_delaunay_2d(points: &[DVec2], opts: &MedialAxisOptions) -> Vec<[usize; 3]> {
    let n = points.len();
    if n < 3 {
        return vec![];
    }

    // Find bounding box
    let mut min_pt = points[0];
    let mut max_pt = points[0];
    for &p in points {
        min_pt = min_pt.min(p);
        max_pt = max_pt.max(p);
    }

    // Create super-triangle that contains all points
    let margin = (max_pt - min_pt).length() * 10.0;
    let super_p0 = DVec2::new(min_pt.x - margin, min_pt.y - margin);
    let super_p1 = DVec2::new(max_pt.x + margin, min_pt.y - margin);
    let super_p2 = DVec2::new((min_pt.x + max_pt.x) / 2.0, max_pt.y + margin);

    // Extended points array with super-triangle nodes
    let mut all_points = points.to_vec();
    all_points.push(super_p0);
    all_points.push(super_p1);
    all_points.push(super_p2);
    let super_idx = n;

    // Initial triangulation: just the super-triangle
    let mut triangles: Vec<[usize; 3]> = vec![[super_idx, super_idx + 1, super_idx + 2]];

    // Insert each point
    for i in 0..n {
        let p = points[i];
        let mut bad_triangles: Vec<usize> = Vec::new();

        // Find all triangles whose circumcircle contains this point
        for (t_idx, tri) in triangles.iter().enumerate() {
            let c0 = all_points[tri[0]];
            let c1 = all_points[tri[1]];
            let c2 = all_points[tri[2]];

            if let Some((center, radius)) = circumcenter(c0, c1, c2)
                && (p - center).length() < radius + opts.tolerance {
                    bad_triangles.push(t_idx);
                }
        }

        // Find the boundary polygon of the cavity
        let mut polygon: Vec<(usize, usize)> = Vec::new();
        for &t_idx in &bad_triangles {
            let tri = triangles[t_idx];
            for k in 0..3 {
                let e1 = tri[k];
                let e2 = tri[(k + 1) % 3];

                // Check if this edge is shared by another bad triangle
                let mut is_shared = false;
                for &other_idx in &bad_triangles {
                    if other_idx != t_idx {
                        let other = triangles[other_idx];
                        for j in 0..3 {
                            if (other[j] == e1 && other[(j + 1) % 3] == e2)
                                || (other[j] == e2 && other[(j + 1) % 3] == e1)
                            {
                                is_shared = true;
                                break;
                            }
                        }
                    }
                    if is_shared {
                        break;
                    }
                }

                if !is_shared {
                    polygon.push((e1, e2));
                }
            }
        }

        // Remove bad triangles
        let mut new_triangles: Vec<[usize; 3]> = Vec::new();
        for (t_idx, tri) in triangles.iter().enumerate() {
            if !bad_triangles.contains(&t_idx) {
                new_triangles.push(*tri);
            }
        }
        triangles = new_triangles;

        // Create new triangles from the polygon boundary
        for (e1, e2) in polygon {
            triangles.push([e1, e2, i]);
        }
    }

    // Remove triangles that contain super-triangle nodes
    triangles.retain(|tri| tri[0] < n && tri[1] < n && tri[2] < n);

    triangles
}

/// Compute the circumcenter and circumradius of a triangle.
fn circumcenter(p0: DVec2, p1: DVec2, p2: DVec2) -> Option<(DVec2, f64)> {
    let d0 = p1 - p0;
    let d1 = p2 - p0;

    let cross = d0.x * d1.y - d0.y * d1.x;
    if cross.abs() < TOLERANCE_FLOAT_DEDUP {
        return None; // Degenerate triangle
    }

    let len0_sq = d0.length_squared();
    let len1_sq = d1.length_squared();

    let s = (len0_sq * d1.y - len1_sq * d0.y) / (2.0 * cross);
    let t = (len0_sq * d1.x - len1_sq * d0.x) / (2.0 * cross);

    let center = p0 + DVec2::new(s, -t);
    let radius = (center - p0).length();

    Some((center, radius))
}

/// Extract the medial axis from a Voronoi diagram by filtering internal edges.
fn extract_medial_axis_from_voronoi(
    voronoi: &VoronoiDiagram2d,
    polygon: &[DVec2],
    _opts: &MedialAxisOptions,
) -> MedialAxis2d {
    let mut result = MedialAxis2d::default();

    // Find Voronoi vertices that are inside the polygon
    let mut inside_vertices: HashSet<usize> = HashSet::new();
    for (i, v) in voronoi.vertices.iter().enumerate() {
        if point_in_polygon_2d(v.point, polygon) {
            inside_vertices.insert(i);
        }
    }

    // Collect internal Voronoi edges as medial axis edges
    let mut medial_points: Vec<MedialPoint2d> = Vec::new();
    let mut medial_edges: Vec<(usize, usize)> = Vec::new();
    let mut point_index_map: HashMap<usize, usize> = HashMap::new();

    for edge in &voronoi.edges {
        if !edge.is_finite {
            continue;
        }

        if let (Some(start_idx), Some(end_idx)) = (edge.start, edge.end) {
            // Both vertices must be inside the polygon
            if inside_vertices.contains(&start_idx) && inside_vertices.contains(&end_idx) {
                // Add start vertex
                let s_idx = if let Some(&idx) = point_index_map.get(&start_idx) {
                    idx
                } else {
                    let idx = medial_points.len();
                    let v = &voronoi.vertices[start_idx];
                    let radius = compute_distance_to_boundary(v.point, polygon);
                    medial_points.push(MedialPoint2d {
                        point: v.point,
                        radius,
                        is_branch: false,
                        is_end: false,
                    });
                    point_index_map.insert(start_idx, idx);
                    idx
                };

                // Add end vertex
                let e_idx = if let Some(&idx) = point_index_map.get(&end_idx) {
                    idx
                } else {
                    let idx = medial_points.len();
                    let v = &voronoi.vertices[end_idx];
                    let radius = compute_distance_to_boundary(v.point, polygon);
                    medial_points.push(MedialPoint2d {
                        point: v.point,
                        radius,
                        is_branch: false,
                        is_end: false,
                    });
                    point_index_map.insert(end_idx, idx);
                    idx
                };

                medial_edges.push((s_idx, e_idx));
            }
        }
    }

    // Identify branch points and end points
    let mut degree = vec![0usize; medial_points.len()];
    for (s, e) in &medial_edges {
        degree[*s] += 1;
        degree[*e] += 1;
    }

    for (i, &deg) in degree.iter().enumerate() {
        if deg > 2 {
            medial_points[i].is_branch = true;
            result.branch_points.push(i);
        } else if deg == 1 {
            medial_points[i].is_end = true;
            result.end_points.push(i);
        }
    }

    // Find maximum inscribed circle
    if !medial_points.is_empty() {
        let max_pt = medial_points
            .iter()
            .max_by(|a, b| a.radius.partial_cmp(&b.radius).unwrap_or(std::cmp::Ordering::Equal));
        if let Some(pt) = max_pt {
            result.max_inscribed_circle = Some((pt.point, pt.radius));
        }
    }

    // Build branches
    result.all_points = medial_points;
    result.branches = build_medial_branches(&result.all_points, &medial_edges, &result.branch_points);

    result
}

/// Build branch structures from the medial axis graph.
fn build_medial_branches(
    points: &[MedialPoint2d],
    edges: &[(usize, usize)],
    branch_points: &[usize],
) -> Vec<MedialBranch2d> {
    if points.is_empty() {
        return vec![];
    }

    let branch_set: HashSet<usize> = branch_points.iter().cloned().collect();

    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![vec![]; points.len()];
    for &(s, e) in edges {
        adj[s].push(e);
        adj[e].push(s);
    }

    let mut visited: HashSet<usize> = HashSet::new();
    let mut branches: Vec<MedialBranch2d> = Vec::new();

    // Start from end points or branch points
    for &start in branch_points {
        if visited.contains(&start) {
            continue;
        }

        // Trace each branch from this branch point
        for &next in &adj[start] {
            if visited.contains(&next) {
                continue;
            }

            let mut branch_pts = vec![start, next];
            visited.insert(start);

            let mut current = next;
            loop {
                visited.insert(current);

                // Find next unvisited neighbor
                let mut found_next = false;
                for &neighbor in &adj[current] {
                    if !visited.contains(&neighbor) && !branch_pts.contains(&neighbor) {
                        branch_pts.push(neighbor);
                        current = neighbor;
                        found_next = true;
                        break;
                    }
                }

                if !found_next {
                    break;
                }

                // Stop at branch points
                if branch_set.contains(&current) {
                    break;
                }
            }

            let branch_points_data: Vec<MedialPoint2d> = branch_pts
                .iter()
                .map(|&i| points[i])
                .collect();

            branches.push(MedialBranch2d {
                points: branch_points_data,
                parent: None,
                children: vec![],
                source_edges: (0, 0), // Would need more info to determine
            });
        }
    }

    // Also trace branches starting from end points
    for (i, pt) in points.iter().enumerate() {
        if pt.is_end && !visited.contains(&i) {
            let mut branch_pts = vec![i];
            visited.insert(i);

            let mut current = i;
            loop {
                let mut found_next = false;
                for &neighbor in &adj[current] {
                    if !visited.contains(&neighbor) {
                        branch_pts.push(neighbor);
                        current = neighbor;
                        visited.insert(current);
                        found_next = true;

                        // Stop at branch points
                        if branch_set.contains(&current) {
                            break;
                        }
                    }
                }

                if !found_next || branch_set.contains(&current) {
                    break;
                }
            }

            let branch_points_data: Vec<MedialPoint2d> = branch_pts
                .iter()
                .map(|&idx| points[idx])
                .collect();

            branches.push(MedialBranch2d {
                points: branch_points_data,
                parent: None,
                children: vec![],
                source_edges: (0, 0),
            });
        }
    }

    branches
}

/// Compute distance from a point to the polygon boundary.
fn compute_distance_to_boundary(point: DVec2, polygon: &[DVec2]) -> f64 {
    let n = polygon.len();
    let mut min_dist = f64::MAX;

    for i in 0..n {
        let p0 = polygon[i];
        let p1 = polygon[(i + 1) % n];
        let d = distance_point_to_segment_2d(point, p0, p1);
        min_dist = min_dist.min(d);
    }

    min_dist
}

/// Distance from a point to a line segment in 2D.
fn distance_point_to_segment_2d(p: DVec2, a: DVec2, b: DVec2) -> f64 {
    let ab = b - a;
    let len_sq = ab.length_squared();
    if len_sq < TOLERANCE_FLOAT_DEDUP {
        return (p - a).length();
    }

    let t = ((p - a).dot(ab) / len_sq).clamp(0.0, 1.0);
    let closest = a + t * ab;
    (p - closest).length()
}

/// Check if a point is inside a 2D polygon using ray casting.
pub fn point_in_polygon_2d(point: DVec2, polygon: &[DVec2]) -> bool {
    if polygon.len() < 3 {
        return false;
    }

    let mut inside = false;
    let n = polygon.len();

    for i in 0..n {
        let j = (i + 1) % n;
        let pi = polygon[i];
        let pj = polygon[j];

        if ((pi.y > point.y) != (pj.y > point.y))
            && (point.x < (pj.x - pi.x) * (point.y - pi.y) / (pj.y - pi.y) + pi.x)
        {
            inside = !inside;
        }
    }

    inside
}

// ============================================================================
// 3D Medial Surface Computation
// ============================================================================

/// Compute the medial axis of a 3D solid (approximate).
///
/// Uses distance field sampling:
/// - Sample points within each face
/// - Compute distance to nearest boundary
/// - Points with local maxima in distance are medial axis candidates
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `opts` - Computation options
///
/// # Returns
/// The computed medial surface.
pub fn compute_medial_surface(brep: &BRep, opts: &MedialAxisOptions) -> MedialSurface {
    let mut result = MedialSurface::default();

    for solid in &brep.solids {
        for shell in &solid.shells {
            compute_shell_medial_surface(shell, brep, opts, &mut result);
        }
    }

    if opts.simplify {
        simplify_medial_surface(&mut result, opts.tolerance);
    }

    // Compute thickness statistics
    compute_thickness_stats(&mut result);

    result
}

fn compute_shell_medial_surface(
    shell: &Shell,
    brep: &BRep,
    opts: &MedialAxisOptions,
    result: &mut MedialSurface,
) {
    // Collect all face surfaces
    let mut face_idx = 0;
    for face in &shell.faces {
        if let Some(&Some(surf_idx)) = brep.geom.face_surface.get(face_idx)
            && let Some(surf) = brep.geom.surfaces.get(surf_idx) {
                sample_surface_medial_points(surf, face, face_idx, brep, opts, result);
            }
        face_idx += 1;
    }
}

fn sample_surface_medial_points(
    surf: &Surface3,
    face: &Face,
    face_idx: usize,
    brep: &BRep,
    opts: &MedialAxisOptions,
    result: &mut MedialSurface,
) {
    let [u_min, u_max, v_min, v_max] = surf.default_domain();

    // Skip unbounded surfaces
    if !u_min.is_finite() || !u_max.is_finite() || !v_min.is_finite() || !v_max.is_finite() {
        return;
    }

    let du = (u_max - u_min) / opts.sample_density as f64;
    let dv = (v_max - v_min) / opts.sample_density as f64;

    let mut samples: Vec<(DVec3, f64)> = Vec::new();

    for i in 0..opts.sample_density {
        for j in 0..opts.sample_density {
            let u = u_min + (i as f64 + 0.5) * du;
            let v = v_min + (j as f64 + 0.5) * dv;

            let point = surf.point_at(u, v);
            let dist = distance_to_boundary_3d(&point, face, brep);

            if dist > opts.min_thickness {
                samples.push((point, dist));
            }
        }
    }

    // Find local maxima in distance field
    let local_maxima = find_local_maxima(&samples, opts.tolerance * 10.0);

    for &idx in &local_maxima {
        let (point, radius) = samples[idx];
        result.vertices.push(MedialVertex {
            point,
            radius,
            boundary_elements: vec![face_idx],
        });
    }
}

/// Find local maxima in a set of distance samples.
fn find_local_maxima(samples: &[(DVec3, f64)], radius: f64) -> Vec<usize> {
    let n = samples.len();
    if n == 0 {
        return vec![];
    }

    let mut maxima = Vec::new();

    for i in 0..n {
        let (p_i, d_i) = samples[i];
        let mut is_max = true;

        for j in 0..n {
            if i == j {
                continue;
            }
            let (p_j, d_j) = samples[j];

            if (p_i - p_j).length() < radius && d_j > d_i {
                is_max = false;
                break;
            }
        }

        if is_max {
            maxima.push(i);
        }
    }

    maxima
}

fn distance_to_boundary_3d(point: &DVec3, face: &Face, brep: &BRep) -> f64 {
    let mut min_dist = f64::MAX;

    // Check distance to outer wire edges
    for we in &face.outer_wire.edges {
        if let Some(&Some(curve_idx)) = brep.geom.edge_curve.get(we.idx)
            && let Some(curve) = brep.geom.curves.get(curve_idx) {
                let [t0, t1] = curve.default_domain();
                if t0.is_finite() && t1.is_finite() {
                    // Sample curve points
                    for k in 0..20 {
                        let t = t0 + (k as f64 / 19.0) * (t1 - t0);
                        let cp = curve.point_at(t);
                        let d = (*point - cp).length();
                        min_dist = min_dist.min(d);
                    }
                }
            }
    }

    // Check distance to inner wire edges
    for wire in &face.inner_wires {
        for we in &wire.edges {
            if let Some(&Some(curve_idx)) = brep.geom.edge_curve.get(we.idx)
                && let Some(curve) = brep.geom.curves.get(curve_idx) {
                    let [t0, t1] = curve.default_domain();
                    if t0.is_finite() && t1.is_finite() {
                        for k in 0..20 {
                            let t = t0 + (k as f64 / 19.0) * (t1 - t0);
                            let cp = curve.point_at(t);
                            let d = (*point - cp).length();
                            min_dist = min_dist.min(d);
                        }
                    }
                }
        }
    }

    min_dist
}

fn simplify_medial_surface(surface: &mut MedialSurface, tolerance: f64) {
    let n = surface.vertices.len();
    if n == 0 {
        return;
    }

    let mut keep = vec![true; n];

    for i in 0..n {
        for j in (i + 1)..n {
            if keep[i] && keep[j] {
                let d = (surface.vertices[i].point - surface.vertices[j].point).length();
                if d < tolerance {
                    // Keep the one with larger radius
                    if surface.vertices[i].radius >= surface.vertices[j].radius {
                        keep[j] = false;
                    } else {
                        keep[i] = false;
                    }
                }
            }
        }
    }

    // Build vertex index mapping
    let mut old_to_new: HashMap<usize, usize> = HashMap::new();
    let mut new_vertices = Vec::new();
    let mut new_idx = 0;

    for (i, v) in surface.vertices.drain(..).enumerate() {
        if keep[i] {
            old_to_new.insert(i, new_idx);
            new_vertices.push(v);
            new_idx += 1;
        }
    }

    surface.vertices = new_vertices;

    // Update edge vertex indices
    for edge in &mut surface.edges {
        if let Some(&new_start) = old_to_new.get(&edge.start_vertex) {
            edge.start_vertex = new_start;
        }
        if let Some(&new_end) = old_to_new.get(&edge.end_vertex) {
            edge.end_vertex = new_end;
        }
    }

    // Remove edges with invalid vertices
    surface.edges.retain(|e| {
        e.start_vertex < surface.vertices.len() && e.end_vertex < surface.vertices.len()
    });

    // Update face vertex indices
    for face in &mut surface.faces {
        face.vertices = face
            .vertices
            .iter()
            .filter_map(|&v| old_to_new.get(&v).copied())
            .collect();
    }
}

fn compute_thickness_stats(surface: &mut MedialSurface) {
    if surface.vertices.is_empty() {
        return;
    }

    let radii: Vec<f64> = surface.vertices.iter().map(|v| v.radius * 2.0).collect();
    let n = radii.len();

    let min = radii.iter().cloned().fold(f64::MAX, f64::min);
    let max = radii.iter().cloned().fold(0.0, f64::max);
    let mean = radii.iter().sum::<f64>() / n as f64;

    let variance = radii.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / n as f64;
    let std_dev = variance.sqrt();

    surface.thickness_stats = ThicknessStats { min, max, mean, std_dev };
}

// ============================================================================
// Enhanced 3D Medial Axis Computation
// ============================================================================

/// Compute the medial axis of a 3D solid using voxel-based distance field.
///
/// This method provides more accurate medial axis extraction for complex 3D geometries
/// by using a voxelized distance field and detecting local maxima.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `opts` - Computation options
///
/// # Returns
/// The computed medial surface with vertices, edges, and faces.
pub fn compute_medial_surface_voxel(brep: &BRep, opts: &MedialAxisOptions) -> MedialSurface {
    let mut result = MedialSurface::default();

    // Compute bounding box of the solid
    let bbox = compute_brep_bbox(brep);
    if !bbox.is_valid() {
        return result;
    }

    // Determine voxel size based on tolerance and min thickness
    let voxel_size = (opts.tolerance * 10.0).max(opts.min_feature_size / 2.0);

    // Create voxel grid
    let dimensions = [
        ((bbox.max.x - bbox.min.x) / voxel_size).ceil() as usize + 2,
        ((bbox.max.y - bbox.min.y) / voxel_size).ceil() as usize + 2,
        ((bbox.max.z - bbox.min.z) / voxel_size).ceil() as usize + 2,
    ];

    let mut grid = VoxelGrid::new(
        bbox.min - DVec3::splat(voxel_size),
        voxel_size,
        dimensions,
    );

    // Compute signed distance field
    compute_signed_distance_field(brep, &mut grid, opts);

    // Find local maxima (medial axis candidates)
    let maxima = grid.find_local_maxima(opts.min_thickness);

    // Convert maxima to medial vertices
    for (i, j, k, distance) in maxima {
        let point = grid.voxel_center(i, j, k);
        result.vertices.push(MedialVertex {
            point,
            radius: distance,
            boundary_elements: vec![],
        });
    }

    // Connect nearby vertices with edges
    connect_medial_vertices(&mut result, opts.cluster_distance);

    // Build medial faces from edge loops
    build_medial_faces(&mut result);

    if opts.simplify {
        simplify_medial_surface(&mut result, opts.tolerance);
    }

    compute_thickness_stats(&mut result);
    result
}

/// Compute the chordal axis of a thin-walled solid.
///
/// The chordal axis is a simplified representation of the medial axis
/// specifically designed for thin-walled parts, capturing the centerlines
/// of sheet-like regions.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `opts` - Computation options
///
/// # Returns
/// The chordal axis with vertices, edges, and thin sheet information.
pub fn compute_chordal_axis(brep: &BRep, opts: &MedialAxisOptions) -> ChordalAxis {
    let mut result = ChordalAxis::default();

    // First compute the medial surface
    let medial = if opts.use_chordal_axis {
        compute_medial_surface_voxel(brep, opts)
    } else {
        compute_medial_surface(brep, opts)
    };

    // Extract face pairs for chordal axis computation
    let face_pairs = compute_opposing_face_pairs(brep, opts);

    // Convert medial vertices to chordal vertices
    for vertex in &medial.vertices {
        // Find associated face pairs for this vertex
        let associated_pairs = find_associated_face_pairs(&vertex.point, &face_pairs, opts.cluster_distance);

        if !associated_pairs.is_empty() {
            // Compute the direction along the thin feature
            let direction = compute_chordal_direction(&vertex.point, &associated_pairs, brep);

            // Compute the normal to the mid-surface
            let normal = compute_mid_surface_normal(&vertex.point, &associated_pairs, brep);

            result.vertices.push(ChordalVertex {
                point: vertex.point,
                thickness: vertex.radius * 2.0,
                direction,
                normal,
                face_pairs: associated_pairs,
            });
        }
    }

    // Connect chordal vertices with edges
    connect_chordal_vertices(&mut result, opts.cluster_distance);

    // Identify thin sheets
    result.sheets = identify_thin_sheets(&result, brep, opts);

    result
}

/// Compute enhanced mid-surface extraction for FEA shell meshing.
///
/// This function extracts the mid-surface from thin-walled solids with
/// improved accuracy and quality metrics suitable for FEA analysis.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `opts` - Mid-surface extraction options
///
/// # Returns
/// Enhanced mid-surface result with quality metrics.
pub fn compute_enhanced_mid_surface(brep: &BRep, opts: &MidSurfaceOptions) -> EnhancedMidSurfaceResult {
    // Compute chordal axis for better thin feature detection
    let chordal_axis = compute_chordal_axis(brep, &opts.base);

    // Create mid-surface B-Rep
    let mut mid_brep = BRep::default();
    let mut face_thickness: Vec<f64> = Vec::new();
    let mut face_mapping: Vec<(usize, usize)> = Vec::new();

    // Create mid-surface faces from chordal sheets
    for sheet in &chordal_axis.sheets {
        let edge_idx = sheet.spine_edge;
        if let Some(edge) = chordal_axis.edges.get(edge_idx) {
            let start_idx = edge.start;
            let end_idx = edge.end;
            if let (Some(start_v), Some(end_v)) = (
                chordal_axis.vertices.get(start_idx),
                chordal_axis.vertices.get(end_idx),
            ) {
                // Create a surface patch between the two vertices
                create_mid_surface_patch(
                    start_v,
                    end_v,
                    sheet,
                    &mut mid_brep,
                    &mut face_thickness,
                    &mut face_mapping,
                    opts,
                );
            }
        }
    }

    // Also create faces for isolated chordal vertices
    for vertex in &chordal_axis.vertices {
        create_mid_surface_point(vertex, &mut mid_brep, &mut face_thickness, &mut face_mapping, opts);
    }

    // Compute quality metrics
    let quality = compute_mid_surface_quality(&mid_brep, brep, &chordal_axis, opts);

    EnhancedMidSurfaceResult {
        brep: mid_brep,
        face_thickness,
        face_mapping,
        chordal_axis,
        quality,
    }
}

/// Detect thin-walled regions with detailed analysis.
///
/// Performs comprehensive thin region detection including clustering,
/// severity classification, and histogram analysis.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `target_thickness` - Target wall thickness for comparison
/// * `opts` - Computation options
///
/// # Returns
/// Detailed thin region analysis with classifications.
pub fn analyze_thin_regions(brep: &BRep, target_thickness: f64, opts: &MedialAxisOptions) -> ThinRegionAnalysis {
    let medial = compute_medial_surface_voxel(brep, opts);

    // Compute basic thin regions
    let mut regions: Vec<ThinRegion> = medial
        .vertices
        .iter()
        .filter(|v| v.radius * 2.0 < target_thickness)
        .map(|v| {
            let thickness = v.radius * 2.0;
            let severity = 1.0 - (thickness / target_thickness).min(1.0);
            ThinRegion {
                center: v.point,
                thickness,
                area: 0.0,
                face_indices: v.boundary_elements.clone(),
                severity,
            }
        })
        .collect();

    // Cluster nearby thin regions
    cluster_thin_regions(&mut regions, opts.cluster_distance);

    // Compute areas for each region
    for region in &mut regions {
        region.area = estimate_region_area(&region.center, region.thickness, &medial);
    }

    // Classify overall thickness
    let classification = classify_thickness(&medial.thickness_stats, target_thickness);

    // Group by severity
    let mut severity_groups: HashMap<ThinRegionSeverity, Vec<usize>> = HashMap::new();
    for (i, region) in regions.iter().enumerate() {
        let severity = if region.severity > 0.75 {
            ThinRegionSeverity::Critical
        } else if region.severity > 0.5 {
            ThinRegionSeverity::Warning
        } else {
            ThinRegionSeverity::Acceptable
        };
        severity_groups.entry(severity).or_default().push(i);
    }

    // Build thickness histogram
    let thickness_histogram = build_thickness_histogram(&medial, 20);

    // Compute recommended minimum thickness
    let recommended_min = compute_recommended_min_thickness(&medial, target_thickness);

    ThinRegionAnalysis {
        regions,
        classification,
        recommended_min,
        severity_groups,
        thickness_histogram,
    }
}

/// Generate optimal rib/stiffener placements along the medial axis.
///
/// Analyzes the medial axis to determine optimal rib placement for
/// structural reinforcement, considering load paths and thickness distribution.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `opts` - Rib generation options
///
/// # Returns
/// Rib generation result with placement recommendations.
pub fn generate_ribs(brep: &BRep, opts: &RibGenerationOptions) -> RibGenerationResult {
    // Compute medial surface
    let medial = compute_medial_surface_voxel(brep, &opts.base);

    // Find candidate rib paths (medial edges with low thickness)
    let candidates = find_rib_candidates(&medial, opts);

    // Generate rib placements
    let mut ribs: Vec<RibPlacement> = Vec::new();

    for candidate in candidates {
        if let Some(placement) = create_rib_placement(&candidate, &medial, brep, opts) {
            ribs.push(placement);
        }
    }

    // Optimize rib layout if requested
    if opts.optimize_stiffness {
        optimize_rib_layout(&mut ribs, brep, opts);
    }

    // Compute statistics
    let total_volume: f64 = ribs.iter().map(|r| {
        // Approximate volume as trapezoidal cross-section
        let length = (r.end - r.start).length();
        let avg_width = r.width;
        let avg_height = r.height;
        length * avg_width * avg_height * 0.5 // Triangular-ish cross-section
    }).sum();

    let stiffness_improvement = estimate_stiffness_improvement(&ribs, &medial);
    let weight_increase = compute_weight_increase(&ribs, brep);
    let quality_score = compute_rib_quality_score(&ribs, &medial, opts);

    RibGenerationResult {
        ribs,
        total_volume,
        stiffness_improvement,
        weight_increase,
        quality_score,
    }
}
include!("e1.rs");
include!("tests_inc.rs");
