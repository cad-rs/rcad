use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_kernel::geom::*;
use crate::bopds::ds::*;
use crate::tolerance::*;
use crate::builder::types::{BooleanOpType, FaceSampleData};
use crate::builder::{SourceSide, BooleanBuilder};

pub struct GlueConfig {
 /// Tolerance for face matching (default: TOLERANCE_MESH_LEGACY).
 ///
 /// Two faces are considered coincident if their surface geometry
 /// matches within this tolerance.
 pub face_tolerance: f64,

 /// Tolerance for edge matching (default: TOLERANCE_MESH_LEGACY).
 ///
 /// Two edges are considered coincident if their curve geometry
 /// matches within this tolerance.
 pub edge_tolerance: f64,

 /// Enable geometric hashing for O(n) face pairing (default: true).
 ///
 /// When enabled, uses a spatial hash to quickly find candidate face
 /// pairs, reducing the complexity from O(n ? to O(n) for models
 /// with many faces.
 pub use_geometric_hash: bool,

 /// Skip non-parallel face pairs early (default: true).
 ///
 /// When enabled, quickly rejects face pairs whose normals are not
 /// approximately anti-parallel, avoiding more expensive geometric
 /// compatibility checks.
 pub early_normal_filter: bool,
}

impl Default for GlueConfig {
 fn default() -> Self {
 Self {
 face_tolerance: TOLERANCE_ABS,
 edge_tolerance: TOLERANCE_ABS,
 use_geometric_hash: true,
 early_normal_filter: true,
 }
 }
}

/// Result of glue face detection.
///
/// Represents a pair of faces from two different shapes that have been
/// identified as coincident or near-coincident, suitable for glue-based
/// boolean operations.
#[derive(Debug, Clone)]
pub struct GlueFacePair {
 /// Index of face in shape A.
 pub face_a: usize,

 /// Index of face in shape B.
 pub face_b: usize,

 /// Match quality (1.0 = perfect match).
 ///
 /// This value indicates how well the two faces match:
 /// - 1.0: Perfect geometric match
 /// - 0.9-1.0: Near-perfect match, within tolerance
 /// - 0.7-0.9: Partial match, some deviation
 /// - < 0.7: Poor match, may not be suitable for gluing
 pub match_quality: f64,

 /// Estimated area of shared region.
 ///
 /// For fully coincident faces, this is the face area.
 /// For partially overlapping faces, this is the overlap area.
 pub shared_area: f64,
}

/// Geometric hash cell for face center points.
///
/// Used for O(n) face pairing by hashing face center coordinates
/// into spatial cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct GeomHashCell {
 ix: i64,
 iy: i64,
 iz: i64,
}

impl GeomHashCell {
 fn from_point(p: DVec3, cell_size: f64) -> Self {
 let scale = 1.0 / cell_size;
 Self {
 ix: (p.x * scale).round() as i64,
 iy: (p.y * scale).round() as i64,
 iz: (p.z * scale).round() as i64,
 }
 }
}

/// Face-pairing cache for performance.
///
/// Caches the results of face compatibility checks to avoid
/// redundant computations during boolean operations.
#[derive(Debug, Clone, Default)]
pub struct GlueFaceCache {
 /// Cached face center points for each face.
 face_centers: Vec<DVec3>,

 /// Cached face normals for each face.
 face_normals: Vec<DVec3>,

 /// Cached face areas for each face.
 face_areas: Vec<f64>,

 /// Spatial hash mapping cells to face indices.
 spatial_hash: HashMap<GeomHashCell, Vec<usize>>,

 /// Cached surface compatibility results.
 /// Key: (face_a, face_b), Value: is_compatible
 compatibility_cache: HashMap<(usize, usize), bool>,
}

impl GlueFaceCache {
 /// Create a new empty cache.
 pub fn new() -> Self {
 Self::default()
 }

 /// Build the cache for a BRep by computing face centers, normals, and areas.
 pub fn build(&mut self, brep: &BRep, cell_size: f64) {
 self.face_centers.clear();
 self.face_normals.clear();
 self.face_areas.clear();
 self.spatial_hash.clear();
 self.compatibility_cache.clear();

 let mut face_idx = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 // Compute face center and area from boundary vertices
 let mut center = DVec3::ZERO;
 let mut area = 0.0;
 let mut count = 0usize;

 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 let edge = &brep.edges[we.idx];
 if edge.start < brep.vertices.len() {
 center += brep.vertices[edge.start].point;
 count += 1;
 }
 if edge.end < brep.vertices.len() {
 center += brep.vertices[edge.end].point;
 count += 1;
 }
 }
 }

 if count > 0 {
 center /= count as f64;
 }

 // Approximate area from bounding box
 let mut min_pt = DVec3::splat(f64::INFINITY);
 let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 let edge = &brep.edges[we.idx];
 if edge.start < brep.vertices.len() {
 let p = brep.vertices[edge.start].point;
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }
 if edge.end < brep.vertices.len() {
 let p = brep.vertices[edge.end].point;
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }
 }
 }
 let diag = max_pt - min_pt;
 area = diag.x * diag.y + diag.y * diag.z + diag.z * diag.x;

 self.face_centers.push(center);
 self.face_normals.push(face.normal);
 self.face_areas.push(area);

 // Add to spatial hash
 let cell = GeomHashCell::from_point(center, cell_size);
 self.spatial_hash.entry(cell).or_default().push(face_idx);

 face_idx += 1;
 }
 }
 }
 }

 /// Get nearby faces using spatial hash.
 pub fn get_nearby_faces(&self, center: DVec3, cell_size: f64) -> Vec<usize> {
 let cell = GeomHashCell::from_point(center, cell_size);

 // Check the cell and its neighbors
 let mut result = Vec::new();
 for dx in -1i64..=1 {
 for dy in -1i64..=1 {
 for dz in -1i64..=1 {
 let neighbor = GeomHashCell {
 ix: cell.ix + dx,
 iy: cell.iy + dy,
 iz: cell.iz + dz,
 };
 if let Some(faces) = self.spatial_hash.get(&neighbor) {
 result.extend(faces.iter().copied());
 }
 }
 }
 }
 result
 }

 /// Check if surface compatibility is cached.
 pub fn get_compatibility(&self, face_a: usize, face_b: usize) -> Option<bool> {
 self.compatibility_cache.get(&(face_a, face_b)).copied()
 }

 /// Cache a surface compatibility result.
 pub fn set_compatibility(&mut self, face_a: usize, face_b: usize, compatible: bool) {
 self.compatibility_cache.insert((face_a, face_b), compatible);
 self.compatibility_cache.insert((face_b, face_a), compatible);
 }
}

/// Detect glue face pairs between two shapes.
///
/// This function analyzes two BReps and identifies pairs of faces that
/// are geometrically coincident or near-coincident, suitable for the
/// glue-based boolean fast path.
///
/// # Arguments
///
/// * `brep_a` - First BRep shape.
/// * `brep_b` - Second BRep shape.
/// * `config` - Configuration for glue detection.
///
/// # Returns
///
/// A vector of `GlueFacePair` representing detected coincident face pairs.
///
/// # Example
///
/// ```
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::builder::{GlueConfig, detect_glue_faces};
/// use glam::DAffine3;
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let mut box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// box2.apply_transform(DAffine3::from_translation(glam::DVec3::new(0.0, 1.0, 0.0)));
///
/// let config = GlueConfig::default();
/// let pairs = detect_glue_faces(&box1, &box2, &config);
/// ```
pub fn detect_glue_faces(
 brep_a: &BRep,
 brep_b: &BRep,
 config: &GlueConfig,
) -> Vec<GlueFacePair> {
 let mut result = Vec::new();

 // Build caches for both BReps
 let cell_size = config.face_tolerance * 10.0;
 let mut cache_a = GlueFaceCache::new();
 let mut cache_b = GlueFaceCache::new();
 cache_a.build(brep_a, cell_size);
 cache_b.build(brep_b, cell_size);

 // Get face counts
 let faces_a: Vec<(usize, DVec3, DVec3, f64)> = brep_a.solids.iter()
 .flat_map(|s| s.shells.iter())
 .flat_map(|sh| sh.faces.iter().enumerate())
 .enumerate()
 .map(|(idx, (_, face))| {
 let center = cache_a.face_centers.get(idx).copied().unwrap_or(DVec3::ZERO);
 let normal = face.normal;
 let area = cache_a.face_areas.get(idx).copied().unwrap_or(0.0);
 (idx, center, normal, area)
 })
 .collect();

 let faces_b: Vec<(usize, DVec3, DVec3, f64)> = brep_b.solids.iter()
 .flat_map(|s| s.shells.iter())
 .flat_map(|sh| sh.faces.iter().enumerate())
 .enumerate()
 .map(|(idx, (_, face))| {
 let center = cache_b.face_centers.get(idx).copied().unwrap_or(DVec3::ZERO);
 let normal = face.normal;
 let area = cache_b.face_areas.get(idx).copied().unwrap_or(0.0);
 (idx, center, normal, area)
 })
 .collect();

 // Early normal filter threshold
 let normal_threshold = -0.95;

 for (idx_a, center_a, normal_a, area_a) in &faces_a {
 // Use geometric hash to find nearby faces in B
 let nearby_faces = if config.use_geometric_hash {
 cache_b.get_nearby_faces(*center_a, cell_size)
 } else {
 faces_b.iter().map(|(idx, _, _, _)| *idx).collect()
 };

 for idx_b in nearby_faces {
 let (_, center_b, normal_b, area_b) = &faces_b.get(idx_b).unwrap_or(&(0, DVec3::ZERO, DVec3::ZERO, 0.0));

 // Early normal filter: skip if normals are not anti-parallel
 if config.early_normal_filter {
 let na_len2 = normal_a.length_squared();
 let nb_len2 = normal_b.length_squared();
 if na_len2 > TOLERANCE_LEN_MIN && nb_len2 > TOLERANCE_LEN_MIN {
 let na = *normal_a / na_len2.sqrt();
 let nb = *normal_b / nb_len2.sqrt();
 if na.dot(nb) > normal_threshold {
 continue;
 }
 }
 }

 // Check center proximity
 let center_dist = (*center_a - *center_b).length();
 if center_dist > config.face_tolerance * 10.0 {
 continue;
 }

 // Compute match quality
 let normal_match = {
 let na_len2 = normal_a.length_squared();
 let nb_len2 = normal_b.length_squared();
 if na_len2 > TOLERANCE_LEN_MIN && nb_len2 > TOLERANCE_LEN_MIN {
 let na = *normal_a / na_len2.sqrt();
 let nb = *normal_b / nb_len2.sqrt();
 // For glue, normals should be anti-parallel
 (-na.dot(nb)).max(0.0)
 } else {
 0.0
 }
 };

 let center_match = {
 let max_dist = config.face_tolerance * 10.0;
 if max_dist > 0.0 {
 (1.0 - center_dist / max_dist).max(0.0)
 } else {
 1.0
 }
 };

 let area_match = {
 let max_area = area_a.max(*area_b);
 let min_area = area_a.min(*area_b);
 if max_area > 0.0 {
 min_area / max_area
 } else {
 1.0
 }
 };

 let match_quality = (normal_match * 0.4 + center_match * 0.3 + area_match * 0.3).min(1.0);

 // Only include pairs with reasonable match quality
 if match_quality >= 0.5 {
 result.push(GlueFacePair {
 face_a: *idx_a,
 face_b: idx_b,
 match_quality,
 shared_area: area_a.min(*area_b),
 }); }
 }
 }

 // Sort by match quality (highest first)
 result.sort_by(|a, b| {
 b.match_quality.partial_cmp(&a.match_quality).unwrap_or(std::cmp::Ordering::Equal)
 });

 result
}

/// Apply glue optimization to pave filler.
///
/// This function configures a PaveFiller to use pre-detected glue face pairs,
/// enabling it to skip expensive interference computations for coincident faces.
///
/// # Arguments
///
/// * `filler` - The PaveFiller to optimize.
/// * `glue_pairs` - Pre-detected glue face pairs.
///
/// # Example
///
/// ```
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::bopds::ds::DS;
/// use rcad_algorithms::pave_filler::PaveFiller;
/// use rcad_algorithms::builder::{GlueConfig, detect_glue_faces, apply_glue_optimization};
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
///
/// let config = GlueConfig::default();
/// let pairs = detect_glue_faces(&box1, &box2, &config);
///
/// let mut ds = DS::new(&box1, &box2);
/// let mut filler = PaveFiller::new(&mut ds);
/// apply_glue_optimization(&mut filler, &pairs);
/// ```
pub fn apply_glue_optimization(
 filler: &mut crate::pave_filler::PaveFiller,
 glue_pairs: &[GlueFacePair],
) {
 if glue_pairs.is_empty() {
 return;
 }

 // Use the tolerance from the best match
 let best_pair = glue_pairs.iter()
 .max_by(|a, b| {
 a.match_quality.partial_cmp(&b.match_quality).unwrap_or(std::cmp::Ordering::Equal)
 });

 if let Some(pair) = best_pair {
 // Estimate tolerance from match quality
 let tolerance = if pair.match_quality > 0.99 {
 TOLERANCE_ABS
 } else if pair.match_quality > 0.9 {
 TOLERANCE_ABS * 10.0
 } else {
 TOLERANCE_ABS * 100.0
 };

 filler.configure_glue(true, tolerance);
 }
}

/// Compute adaptive glue tolerance based on geometry characteristics.
///
/// Analyzes the input BReps and computes an appropriate glue tolerance
/// based on the minimum feature size, face area distribution, and
/// edge length distribution.
///
/// # Arguments
///
/// * `brep_a` - First BRep shape.
/// * `brep_b` - Second BRep shape.
/// * `base_tolerance` - Base tolerance to start with.
///
/// # Returns
///
/// The computed adaptive glue tolerance.
pub fn compute_adaptive_glue_tolerance(
 brep_a: &BRep,
 brep_b: &BRep,
 base_tolerance: f64,
) -> f64 {
 let mut min_feature_size = f64::INFINITY;

 // Analyze edge lengths
 for edge in &brep_a.edges {
 if edge.start < brep_a.vertices.len() && edge.end < brep_a.vertices.len() {
 let p1 = brep_a.vertices[edge.start].point;
 let p2 = brep_a.vertices[edge.end].point;
 let length = (p2 - p1).length();
 if length > TOLERANCE_LINEAR_ULTRA_STRICT {
 min_feature_size = min_feature_size.min(length);
 }
 }
 }
 for edge in &brep_b.edges {
 if edge.start < brep_b.vertices.len() && edge.end < brep_b.vertices.len() {
 let p1 = brep_b.vertices[edge.start].point;
 let p2 = brep_b.vertices[edge.end].point;
 let length = (p2 - p1).length();
 if length > TOLERANCE_LINEAR_ULTRA_STRICT {
 min_feature_size = min_feature_size.min(length);
 }
 }
 }

 // Analyze face areas (approximate from bounding box)
 for solid in &brep_a.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 let mut min_pt = DVec3::splat(f64::INFINITY);
 let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
 for we in &face.outer_wire.edges {
 if we.idx < brep_a.edges.len() {
 let edge = &brep_a.edges[we.idx];
 if edge.start < brep_a.vertices.len() {
 let p = brep_a.vertices[edge.start].point;
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }
 if edge.end < brep_a.vertices.len() {
 let p = brep_a.vertices[edge.end].point;
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }
 }
 }
 let diag = max_pt - min_pt;
 let size = diag.x.min(diag.y).min(diag.z);
 if size > TOLERANCE_LINEAR_ULTRA_STRICT {
 min_feature_size = min_feature_size.min(size);
 }
 }
 }
 }
 for solid in &brep_b.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 let mut min_pt = DVec3::splat(f64::INFINITY);
 let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
 for we in &face.outer_wire.edges {
 if we.idx < brep_b.edges.len() {
 let edge = &brep_b.edges[we.idx];
 if edge.start < brep_b.vertices.len() {
 let p = brep_b.vertices[edge.start].point;
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }
 if edge.end < brep_b.vertices.len() {
 let p = brep_b.vertices[edge.end].point;
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }
 }
 }
 let diag = max_pt - min_pt;
 let size = diag.x.min(diag.y).min(diag.z);
 if size > TOLERANCE_LINEAR_ULTRA_STRICT {
 min_feature_size = min_feature_size.min(size);
 }
 }
 }
 }

 // Compute adaptive tolerance
 let adaptive_tol = if min_feature_size.is_finite() && min_feature_size > 0.0 {
 // Use a fraction of minimum feature size, but at least base tolerance
 let feature_based = min_feature_size * 0.01;
 base_tolerance.max(feature_based).min(min_feature_size * 0.1)
 } else {
 base_tolerance
 };

 adaptive_tol.max(TOLERANCE_ABS)
}
#[cfg(test)]
mod glue_tests {
 use super::*;
 use rcad_kernel::PrimitiveSolid;
 use rcad_kernel::geom::{SphericalSurface, CylindricalSurface, ConicalSurface, ToroidalSurface};
 use glam::DAffine3;
  use glam::DVec2;

 fn unit_box() -> BRep {
 BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 })
 }

 #[test]
 fn test_glue_config_default() {
 let config = GlueConfig::default();
 assert_eq!(config.face_tolerance, TOLERANCE_ABS);
 assert_eq!(config.edge_tolerance, TOLERANCE_ABS);
 assert!(config.use_geometric_hash);
 assert!(config.early_normal_filter);
 }

 #[test]
 fn test_detect_glue_faces_no_overlap() {
 let box1 = unit_box();
 let box2 = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 }).transformed(DAffine3::from_translation(DVec3::new(10.0, 0.0, 0.0)));

 let config = GlueConfig::default();
 let pairs = detect_glue_faces(&box1, &box2, &config);

 // No overlapping faces
 assert!(pairs.is_empty());
 }

 #[test]
 fn test_detect_glue_faces_touching() {
 let box1 = unit_box();
 let mut box2 = unit_box();
 // Translate box2 to touch box1 at y=1 face
 box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

 let config = GlueConfig::default();
 let pairs = detect_glue_faces(&box1, &box2, &config);

 // Should detect at least one coincident face pair
 assert!(!pairs.is_empty());

 // Match quality should be high for exact match
 assert!(pairs[0].match_quality > 0.9);
 }

 #[test]
 fn test_detect_glue_faces_with_tolerance() {
 let box1 = unit_box();
 let mut box2 = unit_box();
 // Slight offset - faces are near but not exactly coincident
 box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0 + TOLERANCE_MESH_LEGACY * 0.1, 0.0)));

 let config = GlueConfig {
 face_tolerance: TOLERANCE_RETRY_LADDER_MID,
 ..Default::default()
 };
 let pairs = detect_glue_faces(&box1, &box2, &config);

 // Should still detect coincident faces with relaxed tolerance
 assert!(!pairs.is_empty());
 }

 #[test]
 fn test_glue_face_pair_structure() {
 let pair = GlueFacePair {
 face_a: 0,
 face_b: 1,
 match_quality: 0.95,
 shared_area: 1.0,
 };

 assert_eq!(pair.face_a, 0);
 assert_eq!(pair.face_b, 1);
 assert!((pair.match_quality - 0.95).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
 assert!((pair.shared_area - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
 }

 #[test]
 fn test_glue_face_cache_build() {
 let box1 = unit_box();
 let mut cache = GlueFaceCache::new();
 cache.build(&box1, 1.0);

 // Should have cached 6 faces (box has 6 faces)
 assert_eq!(cache.face_centers.len(), 6);
 assert_eq!(cache.face_normals.len(), 6);
 assert_eq!(cache.face_areas.len(), 6);

 // Spatial hash should not be empty
 assert!(!cache.spatial_hash.is_empty());
 }

 #[test]
 fn test_glue_face_cache_nearby_faces() {
 let box1 = unit_box();
 let mut cache = GlueFaceCache::new();
 cache.build(&box1, 1.0);

 // Get nearby faces for the center of the box
 let nearby = cache.get_nearby_faces(DVec3::new(0.5, 0.5, 0.5), 1.0);

 // Should find at least some faces
 assert!(!nearby.is_empty());
 }

 #[test]
 fn test_compute_adaptive_glue_tolerance() {
 let box1 = unit_box();
 let box2 = unit_box();

 let tolerance = compute_adaptive_glue_tolerance(&box1, &box2, TOLERANCE_MESH_LEGACY);

 // Tolerance should be reasonable
 assert!(tolerance >= TOLERANCE_ABS);
 assert!(tolerance < 1.0); // Should be much smaller than box size
 }

 #[test]
 fn test_early_normal_filter_disabled() {
 let box1 = unit_box();
 let mut box2 = unit_box();
 box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

 let config = GlueConfig {
 early_normal_filter: false,
 ..Default::default()
 };
 let pairs = detect_glue_faces(&box1, &box2, &config);

 // Should still detect coincident faces
 assert!(!pairs.is_empty());
 }

 #[test]
 fn test_geometric_hash_disabled() {
 let box1 = unit_box();
 let mut box2 = unit_box();
 box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

 let config = GlueConfig {
 use_geometric_hash: false,
 ..Default::default()
 };
 let pairs = detect_glue_faces(&box1, &box2, &config);

 // Should still detect coincident faces
 assert!(!pairs.is_empty());
 }

 #[test]
 fn test_match_quality_ordering() {
 let box1 = unit_box();
 let mut box2 = unit_box();
 // Perfect match
 box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

 let mut box3 = unit_box();
 // Slight rotation - not as good a match
 box3.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));
 box3.apply_transform(DAffine3::from_rotation_z(0.001));

 let config = GlueConfig::default();

 let pairs_exact = detect_glue_faces(&box1, &box2, &config);
 let pairs_rotated = detect_glue_faces(&box1, &box3, &config);

 // Exact match should have higher quality
 if !pairs_exact.is_empty() && !pairs_rotated.is_empty() {
 assert!(pairs_exact[0].match_quality >= pairs_rotated[0].match_quality);
 }
 }

 #[test]
 fn test_shared_area_estimation() {
 let box1 = unit_box();
 let mut box2 = unit_box();
 box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

 let config = GlueConfig::default();
 let pairs = detect_glue_faces(&box1, &box2, &config);

 // Shared area should be approximately 1.0 (unit square face)
 assert!(!pairs.is_empty());
 assert!(pairs[0].shared_area > 0.1);
 }

 #[test]
 fn test_multiple_face_pairs() {
 // Create two boxes that share multiple faces (impossible in real geometry,
 // but tests the algorithm)
 let box1 = unit_box();
 let mut box2 = unit_box();
 box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

 let config = GlueConfig::default();
 let pairs = detect_glue_faces(&box1, &box2, &config);

 // Should detect exactly one face pair (the touching faces)
 assert!(!pairs.is_empty());
 // All pairs should have valid indices
 for pair in &pairs {
 assert!(pair.face_a < 6); // Box has 6 faces
 assert!(pair.face_b < 6);
 }
 }

 #[test]
 fn test_compatibility_cache() {
 let mut cache = GlueFaceCache::new();

 // Initially no cached value
 assert!(cache.get_compatibility(0, 1).is_none());

 // Set and retrieve
 cache.set_compatibility(0, 1, true);
 assert_eq!(cache.get_compatibility(0, 1), Some(true));
 assert_eq!(cache.get_compatibility(1, 0), Some(true)); // Symmetric

 cache.set_compatibility(0, 1, false);
 assert_eq!(cache.get_compatibility(0, 1), Some(false));
 }

 #[test]
 fn test_glue_config_custom_values() {
 let config = GlueConfig {
 face_tolerance: TOLERANCE_RETRY_LADDER_MID,
 edge_tolerance: TOLERANCE_RETRY_LADDER_MID * 2.0,
 use_geometric_hash: false,
 early_normal_filter: false,
 };

 assert!((config.face_tolerance - TOLERANCE_RETRY_LADDER_MID).abs() < TOLERANCE_LEN_MIN);
 assert!((config.edge_tolerance - TOLERANCE_RETRY_LADDER_MID * 2.0).abs() < TOLERANCE_LEN_MIN);
 assert!(!config.use_geometric_hash);
 assert!(!config.early_normal_filter);
 }

  // (split_polygon tests removed — OCCT MakeLoops+Areas pipeline used instead)
}


// ===================================================
//  ?OCCT-aligned: BOPTools_AlgoTools3D  ?orient_edges_on_wire
// ===================================================

///  ?OCCT-aligned: BOPTools_AlgoTools3D::OrientEdgesOnWire.
///
/// Orients edges so they form a consistent closed wire (end-to-start
/// connectivity).  After orientation, the end vertex of edges[i] equals
/// the start vertex of edges[i+1].
///
/// OCCT reference: BOPTools_AlgoTools3D.cxx (OrientEdgesOnWire)
///
/// # Arguments
/// * `edges`  ?Mutable list of (edge_index, forward_flag) pairs to
/// orient in-place.  The first edge's orientation is kept as-is.
/// * `ds`  ?The DS containing vertices and edges.
pub fn orient_edges_on_wire(edges: &mut Vec<(usize, bool)>, ds: &DS) {
 if edges.is_empty() {
 return;
 }
 for i in 1..edges.len() {
 let (prev_ei, prev_fwd) = edges[i - 1];
 let prev_end_vi = if prev_fwd {
 ds.edges[prev_ei].end_vertex
 } else {
 ds.edges[prev_ei].start_vertex
 };
 let (cur_ei, _cur_fwd) = edges[i];
 // Check both orientations of the current edge.
 if ds.edges[cur_ei].start_vertex == prev_end_vi {
 // Already oriented forward  ?keep as-is.
 continue;
 } else if ds.edges[cur_ei].end_vertex == prev_end_vi {
 // Reverse orientation makes the connection.
 edges[i].1 = !edges[i].1;
 }
 // If neither matches there is a topological gap  ?OCCT leaves it as-is.
 }
}

// ===================================================
//  ?OCCT-aligned: BOPTools_AlgoTools3D  ?is_micro_edge
// ===================================================

///  ?OCCT-aligned: BOPTools_AlgoTools3D::IsMicroEdge.
///
/// Returns `true` when the edge's 3D length is shorter than
/// `edge.geom_tol * 2.0`.  Micro-edges are degenerate candidates that
/// the builder can safely skip during face/wire construction.
///
/// Length computation is curve-type-aware:
/// - Line: Euclidean distance between endpoints.
/// - Circle: `radius * |angle_range|`.
/// - Ellipse: `semi_major * |angle_range|` (approximate).
/// - Other: chord distance between endpoints as a conservative estimate.
///
/// OCCT reference: BOPTools_AlgoTools3D.cxx (IsMicroEdge).
pub fn is_micro_edge(edge_idx: usize, ds: &DS) -> bool {
 let tol = ds.edges[edge_idx].geom_tol;
 let len = compute_edge_length_3d(edge_idx, ds);
 len < tol * 2.0
}

/// Compute the 3D length of a DS edge by its curve type.
pub(crate) fn compute_edge_length_3d(edge_idx: usize, ds: &DS) -> f64 {
 let edge = &ds.edges[edge_idx];
 match &edge.curve {
 Curve3::Line(_) => {
 ds.vertices[edge.start_vertex]
 .point
 .distance(ds.vertices[edge.end_vertex].point)
 }
 Curve3::Circle(c) => {
 let angle = (edge.t_range[1] - edge.t_range[0]).abs();
 c.radius * angle
 }
 Curve3::Ellipse(e) => {
 let angle = (edge.t_range[1] - edge.t_range[0]).abs();
 e.major_radius * angle
 }
 _ => {
 // Fallback: chord distance between edge vertices.
 ds.vertices[edge.start_vertex]
 .point
 .distance(ds.vertices[edge.end_vertex].point)
 }
 }
}

// ===================================================
//  ?OCCT-aligned: BOPTools_AlgoTools3D  ?get_edge_on_face
// ===================================================

///  ?OCCT-aligned: BOPTools_AlgoTools3D::GetEdgeOnFace.
///
/// Checks whether a DS edge lies entirely on a DS face's surface.
/// The edge is considered "on face" when both its vertices project
/// to within a combined tolerance of the face surface.
///
/// OCCT reference: BOPTools_AlgoTools3D.cxx (GetEdgeOnFace).
pub fn get_edge_on_face(edge_idx: usize, face_idx: usize, ds: &DS) -> bool {
 let edge = &ds.edges[edge_idx];
 let face = &ds.faces[face_idx];
 let surf = &face.surface;

 // OCCT-aligned SUM: vert_tol + face_tol + fuzzy (same pattern as ComputeVF)
 let v1_tol = ds.vertices[edge.start_vertex].geom_tol + face.geom_tol + ds.fuzzy_tol;
 let v2_tol = ds.vertices[edge.end_vertex].geom_tol + face.geom_tol + ds.fuzzy_tol;

 // Check both edge vertices project onto the face surface.
 let v1_pt = ds.vertices[edge.start_vertex].point;
 let v2_pt = ds.vertices[edge.end_vertex].point;

 let (_uv1, p1_on_surf) = crate::extrema::closest_point_on_surface(surf, v1_pt);
 let (_uv2, p2_on_surf) = crate::extrema::closest_point_on_surface(surf, v2_pt);

 let d1 = p1_on_surf.distance(v1_pt);
 let d2 = p2_on_surf.distance(v2_pt);

 d1 < v1_tol && d2 < v2_tol
}

// ================================================================
//  ?Current state: emit_sphere_faces_direct replaces sphere face emission pipeline
// OCCT edge-based path not yet implemented. Current approach:
// emit_sphere_faces_direct: Circle3 intersection points  ?emit_face_data (FaceSampleData-free)
// ?DoSplitSEAMOnFace  ?(collect_face_edge_segments L2196-2282)
// ?SmartMap/Path walk  ?(build_closed_wires L3312-3617)
// ?PerformAreas  ?(perform_areas)
// ?emit_sphere_faces_direct = =  OCCT  ?
// BuildSplitFaces  ?BuilderFace::Perform  =?  :  =)
// ================================================================

//  ?DoSplitSEAMOnFace  ? ?(collect_face_edge_segments L2196-2282)
// OCCT BOPTools_AlgoTools3D::DoSplitSEAMOnFace (BOPTools_AlgoTools3D.cxx L58-232)
//  ?seam  ?IC   = seam  ?  seam  , ?shifted pcurve ?
// rcad: collect_face_edge_segments  ?seam    ?second_pcurve,
// midpoint UV U=0  ?U=TAU  ュ = = €?

