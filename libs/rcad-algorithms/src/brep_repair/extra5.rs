/// Repair a gap at a periodic surface boundary.
fn repair_periodic_seam_gap(
 result: &mut BRep,
 gap: &PeriodicGap,
 surface_idx: usize,
 surface: &rcad_kernel::geom::Surface3,
 domain: &[f64; 4],
 config: &UvGapRepairConfig,
) -> Result<bool, GapRepairFailureReason> {
 let _ = (result, gap, surface_idx, surface, domain, config);
 // Periodic seam handling is complex and may require:
 // 1. Adjusting the PCurve to wrap correctly
 // 2. Creating a seam edge representation
 // 3. Ensuring continuity across the seam

 // For now, return success without modification
 // A full implementation would adjust PCurve parameters
 Ok(false)
}

/// Repair all UV bounds violations in a BRep.
///
/// This function analyzes all faces in the BRep and attempts to repair
/// any UV bounds violations detected.
///
/// # Arguments
///
/// * `brep` - The BRep to repair.
/// * `config` - Configuration for the repair operations.
///
/// # Returns
///
/// A tuple of (repaired BRep, repair report).
pub fn fix_all_uv_gaps(brep: &BRep, config: &UvGapRepairConfig) -> (BRep, UvGapRepairReport) {
 let mut result = brep.clone();
 let mut total_report = UvGapRepairReport::default();

 // Iterate through all faces
 for (si, solid) in brep.solids.iter().enumerate() {
 for (shi, shell) in solid.shells.iter().enumerate() {
 for (fi, _) in shell.faces.iter().enumerate() {
 let (new_brep, face_report) = fix_uv_gaps(si, shi, fi, &result, config);
 result = new_brep;

 total_report.faces_processed += face_report.faces_processed;
 total_report.gaps_repaired += face_report.gaps_repaired;
 total_report.pcurves_extended += face_report.pcurves_extended;
 total_report.pcurves_trimmed += face_report.pcurves_trimmed;
 total_report.seam_edges_adjusted += face_report.seam_edges_adjusted;
 total_report.unrepaired_gaps.extend(face_report.unrepaired_gaps);
 }
 }
 }

 (result, total_report)
}

/// Repair UV bounds for a specific edge's PCurve.
///
/// This is a more targeted repair function that fixes the PCurve
/// for a specific edge on a specific surface.
///
/// # Arguments
///
/// * `edge_idx` - Index of the edge to repair.
/// * `surface_idx` - Index of the surface for the PCurve.
/// * `brep` - The BRep structure.
/// * `config` - Configuration for the repair operation.
///
/// # Returns
///
/// A tuple of (repaired BRep, whether repair was performed).
pub fn fix_edge_pcurve_uv_bounds(
 edge_idx: usize,
 surface_idx: usize,
 brep: &BRep,
 config: &UvGapRepairConfig,
) -> (BRep, bool) {
 let mut result = brep.clone();
 let mut repaired = false;

 let Some(surface) = brep.geom.surfaces.get(surface_idx) else {
 return (result, repaired);
 };

 let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
 return (result, repaired);
 };

 let domain = surface.default_domain();

 for (pc_idx, pc) in pcurves.iter().enumerate() {
 if pc.surface_idx != surface_idx {
 continue;
 }

 let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else {
 continue;
 };

 let range = brep.geom.curve2d_range.get(pc.curve2d_idx)
 .and_then(|r| *r)
 .unwrap_or([0.0, 1.0]);

 // Sample the PCurve to find bounds
 let mut u_min = f64::INFINITY;
 let mut u_max = f64::NEG_INFINITY;
 let mut v_min = f64::INFINITY;
 let mut v_max = f64::NEG_INFINITY;

 for i in 0..=32 {
 let t = range[0] + (range[1] - range[0]) * i as f64 / 32.0;
 let uv = curve2d.point_at(t);
 u_min = u_min.min(uv.x);
 u_max = u_max.max(uv.x);
 v_min = v_min.min(uv.y);
 v_max = v_max.max(uv.y);
 }

 // Check for violations
 let u_violation_low = domain[0] - u_min;
 let u_violation_high = u_max - domain[1];
 let v_violation_low = domain[2] - v_min;
 let v_violation_high = v_max - domain[3];

 if u_violation_low > config.closure_tolerance ||
 u_violation_high > config.closure_tolerance ||
 v_violation_low > config.closure_tolerance ||
 v_violation_high > config.closure_tolerance {
 // Attempt to wrap or adjust the PCurve
 if let Some(wrapped) = wrap_pcurve_to_domain(curve2d, &range, &domain, config) {
 let new_idx = result.geom.curve2ds.len();
 result.geom.curve2ds.push(wrapped);

 if let Some(pcs) = result.geom.edge_pcurves.get_mut(edge_idx) {
 pcs[pc_idx].curve2d_idx = new_idx;
 }

 repaired = true;
 }
 }
 }

 (result, repaired)
}

/// Wrap a PCurve to fit within the surface domain.
fn wrap_pcurve_to_domain(
 curve2d: &rcad_kernel::Curve2d,
 range: &[f64; 2],
 domain: &[f64; 4],
 config: &UvGapRepairConfig,
) -> Option<rcad_kernel::Curve2d> {
 use rcad_kernel::Curve2d;

 match curve2d {
 Curve2d::Line(line) => {
 let mut new_line = *line;

 // Wrap origin to be within domain
 let u_period = domain[1] - domain[0];
 let v_period = domain[3] - domain[2];

 // Wrap U coordinate
 if new_line.origin.x < domain[0] - config.closure_tolerance {
 new_line.origin.x += u_period;
 } else if new_line.origin.x > domain[1] + config.closure_tolerance {
 new_line.origin.x -= u_period;
 }

 // Wrap V coordinate
 if new_line.origin.y < domain[2] - config.closure_tolerance {
 new_line.origin.y += v_period;
 } else if new_line.origin.y > domain[3] + config.closure_tolerance {
 new_line.origin.y -= v_period;
 }

 Some(Curve2d::Line(new_line))
 }
 Curve2d::BSpline(_) | Curve2d::Circle(_) | Curve2d::Ellipse(_) |
 Curve2d::CircleInvolute(_) | Curve2d::ArchimedeanSpiral(_) |
 Curve2d::LogarithmicSpiral(_) | Curve2d::SineWave(_) | Curve2d::Bezier(_) |
 Curve2d::Trimmed(_) => {
 let _ = range;
 None
 }
 Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => None,
 }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Internal Face Detection and Removal (Post-Boolean Cleanup)
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Classification of duplicate face types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DuplicateFaceKind {
 /// Faces are geometrically identical (same surface, same bounds).
 GeometricallyIdentical,
 /// Faces share topology (same edges, opposite orientation).
 TopologicallyShared,
 /// Faces are coincident but have different geometry representations.
 CoincidentDifferentGeometry,
 /// Faces share the same surface but have different parameter bounds.
 SameSurfaceDifferentBounds,
}

/// Information about a pair of duplicate faces.
#[derive(Debug, Clone)]
pub struct DuplicateFacePair {
 /// Flattened index of the first face.
 pub face_a: usize,
 /// Flattened index of the second face.
 pub face_b: usize,
 /// Classification of the duplication.
 pub kind: DuplicateFaceKind,
 /// Whether the faces have opposite normals.
 pub opposite_orientation: bool,
 /// Maximum geometric deviation between the faces.
 pub max_deviation: f64,
 /// Indices of shared edges (if any).
 pub shared_edges: Vec<usize>,
 /// Whether one face is internal (should be removed).
 pub is_internal: bool,
}

/// Report from duplicate face detection.
#[derive(Debug, Clone, Default)]
pub struct DuplicateFaceReport {
 /// All detected duplicate face pairs.
 pub duplicate_pairs: Vec<DuplicateFacePair>,
 /// Number of faces that are internal candidates for removal.
 pub internal_face_count: usize,
 /// Indices of faces identified as internal.
 pub internal_face_indices: Vec<usize>,
 /// Summary string for debugging.
 pub summary: String,
}

/// Detect duplicate faces in a BRep using geometric and topological comparison.
///
/// This function identifies faces that are geometrically or topologically
/// duplicated, which commonly occurs after boolean operations.
///
/// # Arguments
/// * `brep` - The BRep to analyze.
/// * `tolerance` - Maximum distance for considering geometry coincident.
///
/// # Returns
/// A `DuplicateFaceReport` containing all detected duplicate pairs.
///
/// # Example
/// ```
/// use rcad_algorithms::brep_repair::detect_duplicate_faces;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
/// width: 1.0,
/// height: 1.0,
/// depth: 1.0,
/// });
///
/// let report = detect_duplicate_faces(&brep, TOLERANCE_MESH_LEGACY);
/// println!("Found {} duplicate pairs", report.duplicate_pairs.len());
/// ```
pub fn detect_duplicate_faces(brep: &BRep, tolerance: f64) -> DuplicateFaceReport {
 let tol = tolerance.max(TOLERANCE_ABS);
 let mut report = DuplicateFaceReport::default();

 // Collect all faces with their flattened indices
 let faces: Vec<(usize, usize, usize, &Face)> = brep
 .solids
 .iter()
 .enumerate()
 .flat_map(|(si, solid)| {
 solid.shells.iter().enumerate().flat_map(move |(shi, shell)| {
 shell.faces.iter().enumerate().map(move |(fi, face)| (si, shi, fi, face))
 })
 })
 .collect();

 let n_faces = faces.len();
 if n_faces < 2 {
 report.summary = "No faces to compare".to_string();
 return report;
 }

 // Build surface compatibility map
 let surface_map = build_surface_compatibility_map(brep, &faces, tol);

 // Compare each pair of faces
 let mut processed: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();

 for i in 0..n_faces {
 for j in (i + 1)..n_faces {
 if processed.contains(&(i, j)) {
 continue;
 }

 let (si1, shi1, fi1, face1) = faces[i];
 let (si2, shi2, fi2, face2) = faces[j];

 // Skip faces in the same shell at the same position
 if si1 == si2 && shi1 == shi2 && fi1 == fi2 {
 continue;
 }

 if let Some(pair) = analyze_face_duplication(
 brep,
 face1,
 face2,
 i,
 j,
 &surface_map,
 tol,
 ) {
 processed.insert((i, j));

 // Check if this face is internal
 let is_internal = check_if_internal(brep, &faces, i, j, &pair, tol);
 let mut pair = pair;
 pair.is_internal = is_internal;

 if is_internal {
 report.internal_face_indices.push(j); // Remove the second face
 }

 report.duplicate_pairs.push(pair);
 }
 }
 }

 report.internal_face_count = report.internal_face_indices.len();
 report.summary = format!(
 "DuplicateFaceReport: {} pairs found, {} internal faces",
 report.duplicate_pairs.len(),
 report.internal_face_count
 );

 report
}

/// Build a map of surface compatibility between faces.
fn build_surface_compatibility_map(
 brep: &BRep,
 faces: &[(usize, usize, usize, &Face)],
 tolerance: f64,
) -> std::collections::HashMap<(usize, usize), bool> {
 let mut map: std::collections::HashMap<(usize, usize), bool> = std::collections::HashMap::new();

 for (i, (_, _, _, face1)) in faces.iter().enumerate() {
 for (j, (_, _, _, face2)) in faces.iter().enumerate() {
 if i >= j {
 continue;
 }

 // Check if faces have compatible surfaces
 let compatible = check_surface_compatibility(brep, face1, face2, tolerance);
 map.insert((i, j), compatible);
 }
 }

 map
}

/// Check if two faces have compatible surfaces.
fn check_surface_compatibility(
 brep: &BRep,
 face1: &Face,
 face2: &Face,
 tolerance: f64,
) -> bool {
 // First check normal compatibility - duplicate faces should have parallel or anti-parallel normals
 let normal_dot = face1.normal.dot(face2.normal);
 if normal_dot.abs() < 0.99 {
 return false;
 }

 // Check geometric bounds compatibility
 let pts1: Vec<DVec3> = face1
 .outer_wire
 .edges
 .iter()
 .filter_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 let vidx = if we.forward { edge.start } else { edge.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 .collect();

 let pts2: Vec<DVec3> = face2
 .outer_wire
 .edges
 .iter()
 .filter_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 let vidx = if we.forward { edge.start } else { edge.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 .collect();

 if pts1.is_empty() || pts2.is_empty() {
 return false;
 }

 // Check bounding box overlap
 let (min1, max1) = compute_bounding_box(&pts1);
 let (min2, max2) = compute_bounding_box(&pts2);

 // Allow some tolerance for bounding box comparison
 

 (min1.x - tolerance <= max2.x && max1.x + tolerance >= min2.x) &&
 (min1.y - tolerance <= max2.y && max1.y + tolerance >= min2.y) &&
 (min1.z - tolerance <= max2.z && max1.z + tolerance >= min2.z)
}

/// Compute bounding box of a set of points.
fn compute_bounding_box(points: &[DVec3]) -> (DVec3, DVec3) {
 if points.is_empty() {
 return (DVec3::ZERO, DVec3::ZERO);
 }

 let mut min_pt = points[0];
 let mut max_pt = points[0];

 for &p in points.iter().skip(1) {
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }

 (min_pt, max_pt)
}

/// Analyze two faces for duplication.
fn analyze_face_duplication(
 brep: &BRep,
 face1: &Face,
 face2: &Face,
 flat_idx1: usize,
 flat_idx2: usize,
 surface_map: &std::collections::HashMap<(usize, usize), bool>,
 tolerance: f64,
) -> Option<DuplicateFacePair> {
 // Check surface compatibility
 let surface_compatible = surface_map
 .get(&(flat_idx1.min(flat_idx2), flat_idx1.max(flat_idx2)))
 .copied()
 .unwrap_or(false);

 if !surface_compatible {
 return None;
 }

 // Collect boundary vertices for both faces
 let pts1: Vec<DVec3> = face1
 .outer_wire
 .edges
 .iter()
 .filter_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 let vidx = if we.forward { edge.start } else { edge.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 .collect();

 let pts2: Vec<DVec3> = face2
 .outer_wire
 .edges
 .iter()
 .filter_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 let vidx = if we.forward { edge.start } else { edge.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 .collect();

 // Compare vertex positions
 let _tol_sq = tolerance * tolerance;
 let mut matched_vertices = 0;
 let mut max_deviation = 0.0f64;

 for &p1 in &pts1 {
 let mut best_dist = f64::INFINITY;
 for &p2 in &pts2 {
 let dist_sq = (p1 - p2).length_squared();
 if dist_sq < best_dist {
 best_dist = dist_sq;
 }
 }
 let dist = best_dist.sqrt();
 max_deviation = max_deviation.max(dist);
 if dist <= tolerance {
 matched_vertices += 1;
 }
 }

 // Require most vertices to match
 let match_ratio = matched_vertices as f64 / pts1.len().max(1) as f64;
 if match_ratio < 0.8 {
 return None;
 }

 // Check for shared edges
 let edges1: std::collections::HashSet<usize> =
 face1.outer_wire.edges.iter().map(|we| we.idx).collect();
 let edges2: std::collections::HashSet<usize> =
 face2.outer_wire.edges.iter().map(|we| we.idx).collect();

 let shared_edges: Vec<usize> = edges1.intersection(&edges2).copied().collect();

 // Determine duplication kind
 let kind = if shared_edges.len() == edges1.len() && shared_edges.len() == edges2.len() {
 // All edges are shared - topologically identical
 if max_deviation < tolerance * 0.1 {
 DuplicateFaceKind::GeometricallyIdentical
 } else {
 DuplicateFaceKind::CoincidentDifferentGeometry
 }
 } else if !shared_edges.is_empty() {
 // Some edges shared
 DuplicateFaceKind::TopologicallyShared
 } else {
 // No shared edges but geometrically close
 DuplicateFaceKind::SameSurfaceDifferentBounds
 };

 // Check orientation
 let normal_dot = face1.normal.dot(face2.normal);
 let opposite_orientation = normal_dot < -0.99;

 Some(DuplicateFacePair {
 face_a: flat_idx1,
 face_b: flat_idx2,
 kind,
 opposite_orientation,
 max_deviation,
 shared_edges,
 is_internal: false, // Will be set later
 })
}

/// Check if a face pair indicates one face is internal.
fn check_if_internal(
 brep: &BRep,
 faces: &[(usize, usize, usize, &Face)],
 flat_idx1: usize,
 flat_idx2: usize,
 pair: &DuplicateFacePair,
 _tolerance: f64,
) -> bool {
 // A face is considered internal if:
 // 1. It's a duplicate with opposite orientation
 // 2. It's inside another solid
 // 3. It belongs to a void shell (internal shell in a solid)

 let (si1, shi1, _, _) = faces[flat_idx1];
 let (si2, shi2, _, _) = faces[flat_idx2];

 // If faces are in different solids, check for containment
 if si1 != si2 {
 // For now, consider the second face as potentially internal
 // A more sophisticated check would do ray casting
 return pair.opposite_orientation;
 }

 // If in the same solid but different shells
 if shi1 != shi2 {
 // Check if one shell is internal (void)
 // Shell index > 0 in a solid typically indicates a void
 let solid = &brep.solids[si1];
 if shi2 > 0 && shi2 < solid.shells.len() {
 // Second shell is likely a void shell
 return true;
 }
 }

 // If faces have opposite orientation and are geometrically identical
 pair.opposite_orientation && matches!(
 pair.kind,
 DuplicateFaceKind::GeometricallyIdentical | DuplicateFaceKind::CoincidentDifferentGeometry
 )
}

/// Identify internal faces in a BRep using geometric analysis.
///
/// Internal faces are faces that are completely contained within the solid
/// and do not contribute to the outer boundary. These typically arise from
/// boolean operations where internal separator faces are not removed.
///
/// # Arguments
/// * `brep` - The BRep to analyze.
///
/// # Returns
/// A vector of flattened face indices that are identified as internal.
///
/// # Detection Methods
/// 1. Faces with zero outward normal contribution (sandwiched between other faces)
/// 2. Faces in void shells (shell index > 0 in a solid)
/// 3. Duplicate faces with opposite orientation
/// 4. Faces completely inside other solids (via ray casting)
pub fn identify_internal_faces(brep: &BRep) -> Vec<usize> {
 let mut internal_faces = Vec::new();

 // Method 1: Check for void shells (internal cavities)
 for (si, solid) in brep.solids.iter().enumerate() {
 if solid.shells.len() > 1 {
 // First shell is typically the outer shell
 // Subsequent shells are voids (internal cavities)
 // Faces in void shells with inverted normals are internal separators
 for shi in 1..solid.shells.len() {
 let mut flat_idx = 0usize;
 for (prev_si, prev_solid) in brep.solids.iter().enumerate() {
 for (prev_shi, prev_shell) in prev_solid.shells.iter().enumerate() {
 if prev_si == si && prev_shi == shi {
 // This is a void shell - add all its faces
 for fi in 0..prev_shell.faces.len() {
 internal_faces.push(flat_idx + fi);
 }
 }
 flat_idx += prev_shell.faces.len();
 }
 }
 }
 }
 }

 // Method 2: Check for duplicate faces with opposite orientation
 let duplicate_report = detect_duplicate_faces(brep, TOLERANCE_MESH_LEGACY);
 for pair in &duplicate_report.duplicate_pairs {
 if pair.opposite_orientation && pair.is_internal {
 // Add the second face (the one that should be removed)
 if !internal_faces.contains(&pair.face_b) {
 internal_faces.push(pair.face_b);
 }
 }
 }

 // Method 3: Check for faces with no volume contribution using ray casting
 let ray_internal = identify_internal_faces_by_raycast(brep);
 for idx in ray_internal {
 if !internal_faces.contains(&idx) {
 internal_faces.push(idx);
 }
 }

 // Sort and deduplicate
 internal_faces.sort();
 internal_faces.dedup();

 internal_faces
}

/// Identify internal faces using ray casting.
fn identify_internal_faces_by_raycast(brep: &BRep) -> Vec<usize> {
 let mut internal_faces = Vec::new();

 // Collect all faces with their flattened indices and centroids
 let faces: Vec<(usize, &Face)> = brep
 .solids
 .iter()
 .flat_map(|solid| solid.shells.iter())
 .flat_map(|shell| shell.faces.iter())
 .enumerate()
 .collect();

 if faces.is_empty() {
 return internal_faces;
 }

 // For each face, cast a ray along its normal and check if it's inside the solid
 for (flat_idx, face) in &faces {
 // Compute face centroid
 let centroid = compute_face_centroid_from_wire(brep, face);
 if centroid.is_nan() {
 continue;
 }

 // Cast ray along the face normal
 let ray_origin = centroid + face.normal * TOLERANCE_RETRY_LADDER_COARSE; // Offset slightly
 let ray_dir = face.normal;

 // Count intersections with other faces
 let mut intersection_count = 0;
 for (other_idx, other_face) in &faces {
 if *other_idx == *flat_idx {
 continue;
 }

 if ray_intersects_face(brep, other_face, ray_origin, ray_dir) {
 intersection_count += 1;
 }
 }

 // If odd number of intersections in the direction of the normal,
 // the face is likely internal
 if intersection_count > 0 && intersection_count % 2 == 1 {
 internal_faces.push(*flat_idx);
 }
 }

 internal_faces
}

/// Compute the centroid of a face from its wire vertices.
fn compute_face_centroid_from_wire(brep: &BRep, face: &Face) -> DVec3 {
 let pts: Vec<DVec3> = face
 .outer_wire
 .edges
 .iter()
 .filter_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 let vidx = if we.forward { edge.start } else { edge.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 .collect();

 if pts.is_empty() {
 return DVec3::NAN;
 }

 pts.iter().sum::<DVec3>() / pts.len() as f64
}

/// Check if a ray intersects a face.
fn ray_intersects_face(
 brep: &BRep,
 face: &Face,
 ray_origin: DVec3,
 ray_dir: DVec3,
) -> bool {
 // Get face vertices
 let pts: Vec<DVec3> = face
 .outer_wire
 .edges
 .iter()
 .filter_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 let vidx = if we.forward { edge.start } else { edge.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 .collect();

 if pts.len() < 3 {
 return false;
 }

 // Use M ler= rumbore algorithm for ray-triangle intersection
 // Triangulate the face using fan triangulation
 for i in 1..pts.len() - 1 {
 let v0 = pts[0];
 let v1 = pts[i];
 let v2 = pts[i + 1];

 if ray_triangle_intersection(ray_origin, ray_dir, v0, v1, v2) {
 return true;
 }
 }

 false
}

/// M ler= rumbore ray-triangle intersection.
fn ray_triangle_intersection(
 origin: DVec3,
 dir: DVec3,
 v0: DVec3,
 v1: DVec3,
 v2: DVec3,
) -> bool {
 const EPSILON: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

 let edge1 = v1 - v0;
 let edge2 = v2 - v0;

 let h = dir.cross(edge2);
 let a = edge1.dot(h);

 if a.abs() < EPSILON {
 return false;
 }

 let f = 1.0 / a;
 let s = origin - v0;
 let u = f * s.dot(h);

 if !(0.0..=1.0).contains(&u) {
 return false;
 }

 let q = s.cross(edge1);
 let v = f * dir.dot(q);

 if v < 0.0 || u + v > 1.0 {
 return false;
 }

 let t = f * edge2.dot(q);

 t > EPSILON
}

/// Report from internal face removal.
#[derive(Debug, Clone, Default)]
pub struct InternalFaceRemovalReport {
 /// Number of faces removed.
 pub faces_removed: usize,
 /// Indices of faces that were removed.
 pub removed_indices: Vec<usize>,
 /// Number of edges that became orphaned and were removed.
 pub edges_removed: usize,
 /// Number of vertices that became orphaned and were removed.
 pub vertices_removed: usize,
 /// Whether the result is valid.
 pub is_valid: bool,
}

/// Remove internal faces from a BRep while maintaining topology consistency.
///
/// This function safely removes specified internal faces, updating shell
/// references and handling edge sharing correctly.
///
/// # Arguments
/// * `brep` - The BRep to modify.
/// * `face_indices` - Flattened indices of faces to remove.
///
/// # Returns
/// A new BRep with the internal faces removed and a report of changes.
///
/// # Topology Handling
/// - Removes faces from shells
/// - Removes orphaned edges (edges no longer referenced by any face)
/// - Removes orphaned vertices (vertices no longer referenced by any edge)
/// - Updates geometric data arrays to match new topology
pub fn remove_internal_faces(brep: &BRep, face_indices: &[usize]) -> (BRep, InternalFaceRemovalReport) {
 let mut report = InternalFaceRemovalReport::default();
 let remove_set: std::collections::HashSet<usize> = face_indices.iter().copied().collect();

 if remove_set.is_empty() {
 report.is_valid = true;
 return (brep.clone(), report);
 }

 // Build a map from flat face index to (solid_idx, shell_idx, face_idx)
 let mut flat_to_local: std::collections::HashMap<usize, (usize, usize, usize)> =
 std::collections::HashMap::new();
 let mut flat_idx = 0usize;

 for (si, solid) in brep.solids.iter().enumerate() {
 for (shi, shell) in solid.shells.iter().enumerate() {
 for fi in 0..shell.faces.len() {
 flat_to_local.insert(flat_idx, (si, shi, fi));
 flat_idx += 1;
 }
 }
 }

 // Identify edges to keep (edges referenced by faces NOT being removed)
 let mut edges_to_keep: std::collections::HashSet<usize> = std::collections::HashSet::new();

 for (flat_idx, (_, _, _face)) in flat_to_local.iter().flat_map(|(idx, &(si, shi, fi))| {
 let face = &brep.solids[si].shells[shi].faces[fi];
 Some((idx, (si, shi, face)))
 }) {
 if !remove_set.contains(flat_idx) {
 // Collect all edges from this face's wires
 let face = &brep.solids[flat_to_local[flat_idx].0]
 .shells[flat_to_local[flat_idx].1]
 .faces[flat_to_local[flat_idx].2];

 for we in &face.outer_wire.edges {
 edges_to_keep.insert(we.idx);
 }
 for inner in &face.inner_wires {
 for we in &inner.edges {
 edges_to_keep.insert(we.idx);
 }
 }
 }
 }

 // Also collect edges from faces being kept
 flat_idx = 0;
 for solid in brep.solids.iter() {
 for shell in solid.shells.iter() {
 for face in shell.faces.iter() {
 if !remove_set.contains(&flat_idx) {
 for we in &face.outer_wire.edges {
 edges_to_keep.insert(we.idx);
 }
 for inner in &face.inner_wires {
 for we in &inner.edges {
 edges_to_keep.insert(we.idx);
 }
 }
 }
 flat_idx += 1;
 }
 }
 }

 // Build new solids with faces removed
 let mut new_solids: Vec<Solid> = Vec::new();
 flat_idx = 0;

 for solid in &brep.solids {
 let mut new_shells: Vec<Shell> = Vec::new();

 for shell in &solid.shells {
 let mut new_faces: Vec<Face> = Vec::new();

 for face in &shell.faces {
 if remove_set.contains(&flat_idx) {
 report.faces_removed += 1;
 report.removed_indices.push(flat_idx);
 } else {
 new_faces.push(face.clone());
 }
 flat_idx += 1;
 }

 // Only add shell if it has faces
 if !new_faces.is_empty() {
 new_shells.push(Shell { faces: new_faces });
 }
 }

 // Only add solid if it has shells
 if !new_shells.is_empty() {
 new_solids.push(Solid { shells: new_shells });
 }
 }

 // Create result BRep
 let mut result = BRep::new();
 result.vertices = brep.vertices.clone();
 result.edges = brep.edges.clone();
 result.solids = new_solids;
 result.geom = brep.geom.clone();

 // Remove orphaned edges
 let old_edge_count = result.edges.len();
 let (cleaned_brep, edge_remap) = remove_orphaned_edges(&result, &edges_to_keep);
 result = cleaned_brep;
 report.edges_removed = old_edge_count - result.edges.len();

 // Remove orphaned vertices
 let old_vertex_count = result.vertices.len();
 let cleaned_brep = remove_orphaned_vertices(&result);
 result = cleaned_brep;
 report.vertices_removed = old_vertex_count - result.vertices.len();

 // Update geometric data arrays
 result = update_geom_after_removal(&result, &edge_remap);

 report.is_valid = true;
 (result, report)
}

/// Remove edges that are no longer referenced by any face.
fn remove_orphaned_edges(
 brep: &BRep,
 edges_to_keep: &std::collections::HashSet<usize>,
) -> (BRep, std::collections::HashMap<usize, usize>) {
 let _n_edges = brep.edges.len();

 // Build remap: old_idx -> new_idx
 let mut remap: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 let mut new_edges: Vec<Edge> = Vec::new();

 for (old_idx, edge) in brep.edges.iter().enumerate() {
 if edges_to_keep.contains(&old_idx) {
 let new_idx = new_edges.len();
 new_edges.push(*edge);
 remap.insert(old_idx, new_idx);
 }
 }

 // Update wires to use new edge indices
 let mut result = brep.clone();
 result.edges = new_edges;

 // Update face wires with remapped edge indices
 for solid in &mut result.solids {
 for shell in &mut solid.shells {
 for face in &mut shell.faces {
 // Update outer wire
 for we in &mut face.outer_wire.edges {
 if let Some(&new_idx) = remap.get(&we.idx) {
 we.idx = new_idx;
 }
 }
 // Update inner wires
 for inner in &mut face.inner_wires {
 for we in &mut inner.edges {
 if let Some(&new_idx) = remap.get(&we.idx) {
 we.idx = new_idx;
 }
 }
 }
 }
 }
 }

 // Update geom store edge-related arrays
 // We keep the shared pools (curves, surfaces, curve2ds) intact
 // and only remap edge-level indices

 // Rebuild edge_curve with remapped indices (these index into the shared curves pool)
 let mut new_edge_curve: Vec<Option<usize>> = vec![None; result.edges.len()];
 for (&old_idx, &new_idx) in remap.iter() {
 if new_idx < new_edge_curve.len() {
 new_edge_curve[new_idx] = brep.geom.edge_curve.get(old_idx).copied().flatten();
 }
 }

 // Rebuild edge_curve_range
 let mut new_edge_curve_range: Vec<Option<[f64; 2]>> = vec![None; result.edges.len()];
 for (&old_idx, &new_idx) in remap.iter() {
 if new_idx < new_edge_curve_range.len() {
 new_edge_curve_range[new_idx] = brep.geom.edge_curve_range.get(old_idx).copied().flatten();
 }
 }

 // Rebuild edge_tolerance
 let mut new_edge_tolerance: Vec<f64> = vec![0.0; result.edges.len()];
 for (&old_idx, &new_idx) in remap.iter() {
 if new_idx < new_edge_tolerance.len() {
 new_edge_tolerance[new_idx] = brep.geom.edge_tolerance.get(old_idx).copied().unwrap_or(0.0);
 }
 }

 // Rebuild edge_pcurves
 let mut new_edge_pcurves: Vec<Vec<rcad_kernel::PCurve>> = vec![Vec::new(); result.edges.len()];
 for (&old_idx, &new_idx) in remap.iter() {
 if new_idx < new_edge_pcurves.len()
 && let Some(pcurves) = brep.geom.edge_pcurves.get(old_idx) {
 new_edge_pcurves[new_idx] = pcurves.clone();
 }
 }

 // Rebuild edge_degenerated
 let mut new_edge_degenerated: Vec<bool> = vec![false; result.edges.len()];
 for (&old_idx, &new_idx) in remap.iter() {
 if new_idx < new_edge_degenerated.len() {
 new_edge_degenerated[new_idx] = brep.geom.edge_degenerated.get(old_idx).copied().unwrap_or(false);
 }
 }

 // Rebuild edge_same_parameter
 let mut new_edge_same_parameter: Vec<bool> = vec![true; result.edges.len()];
 for (&old_idx, &new_idx) in remap.iter() {
 if new_idx < new_edge_same_parameter.len() {
 new_edge_same_parameter[new_idx] = brep.geom.edge_same_parameter.get(old_idx).copied().unwrap_or(true);
 }
 }

 // Rebuild edge_same_range
 let mut new_edge_same_range: Vec<bool> = vec![true; result.edges.len()];
 for (&old_idx, &new_idx) in remap.iter() {
 if new_idx < new_edge_same_range.len() {
 new_edge_same_range[new_idx] = brep.geom.edge_same_range.get(old_idx).copied().unwrap_or(true);
 }
 }

 result.geom.edge_curve = new_edge_curve;
 result.geom.edge_curve_range = new_edge_curve_range;
 result.geom.edge_tolerance = new_edge_tolerance;
 result.geom.edge_pcurves = new_edge_pcurves;
 result.geom.edge_degenerated = new_edge_degenerated;
 result.geom.edge_same_parameter = new_edge_same_parameter;
 result.geom.edge_same_range = new_edge_same_range;

 (result, remap)
}

/// Remove vertices that are no longer referenced by any edge.
fn remove_orphaned_vertices(brep: &BRep) -> BRep {
 // Find all vertices that are referenced by edges
 let mut vertices_used: std::collections::HashSet<usize> = std::collections::HashSet::new();

 for edge in &brep.edges {
 vertices_used.insert(edge.start);
 vertices_used.insert(edge.end);
 }

 // Build remap
 let n_verts = brep.vertices.len();
 let mut remap: Vec<usize> = vec![0; n_verts];
 let mut new_vertices: Vec<Vertex> = Vec::new();

 for (old_idx, vertex) in brep.vertices.iter().enumerate() {
 if vertices_used.contains(&old_idx) {
 let new_idx = new_vertices.len();
 new_vertices.push(*vertex);
 remap[old_idx] = new_idx;
 }
 }

 // Update edges with new vertex indices
 let mut result = brep.clone();
 result.vertices = new_vertices;

 for edge in &mut result.edges {
 edge.start = remap[edge.start];
 edge.end = remap[edge.end];
 }

 // Update vertex tolerance array
 let mut new_vertex_tolerance: Vec<f64> = vec![0.0; result.vertices.len()];
 for (old_idx, &new_idx) in remap.iter().enumerate() {
 if let Some(&tol) = brep.geom.vertex_tolerance.get(old_idx)
 && new_idx < new_vertex_tolerance.len() {
 new_vertex_tolerance[new_idx] = tol;
 }
 }
 result.geom.vertex_tolerance = new_vertex_tolerance;

 result
}

/// Update geometric data arrays after edge removal.
fn update_geom_after_removal(
 brep: &BRep,
 edge_remap: &std::collections::HashMap<usize, usize>,
) -> BRep {
 let mut result = brep.clone();

 // Update pcurve references to use new edge indices
 for (old_idx, &new_idx) in edge_remap {
 if let Some(pcurves) = brep.geom.edge_pcurves.get(*old_idx).cloned()
 && new_idx < result.geom.edge_pcurves.len() {
 result.geom.edge_pcurves[new_idx] = pcurves;
 }
 }

 result
}

/// Report from boolean cleanup.
#[derive(Debug, Clone, Default)]
pub struct BooleanCleanupReport {
 /// Number of internal faces removed.
 pub internal_faces_removed: usize,
 /// Number of duplicate faces merged.
 pub duplicate_faces_merged: usize,
 /// Number of vertices merged.
 pub vertices_merged: usize,
 /// Number of degenerate faces removed.
 pub degenerate_faces_removed: usize,
 /// Number of edges sewn.
 pub edges_sewn: usize,
 /// Whether the result is valid.
 pub is_valid: bool,
 /// Summary string.
 pub summary: String,
}

/// Clean up a BRep after boolean operations.
///
/// This function applies a comprehensive cleanup pipeline designed to
/// remove artifacts commonly produced by boolean operations:
///
/// 1. Remove internal faces (separator faces between merged volumes)
/// 2. Merge duplicate faces
/// 3. Remove degenerate faces
/// 4. Merge close vertices
/// 5. Sew close edges
/// 6. Fix tolerances
///
/// # Arguments
/// * `brep` - The BRep to clean up.
/// * `tolerance` - Tolerance for geometric comparisons.
///
/// # Returns
/// A cleaned BRep and a report of all changes made.
///
/// # Example
/// ```
/// use rcad_algorithms::brep_repair::cleanup_boolean_result;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
///
/// // After a boolean operation, clean up the result
/// fn process_boolean_result(result: &BRep) -> BRep {
/// let (cleaned, report) = cleanup_boolean_result(result, TOLERANCE_MESH_LEGACY);
/// println!("Cleaned: {} internal faces removed", report.internal_faces_removed);
/// cleaned
/// }
/// ```
pub fn cleanup_boolean_result(brep: &BRep, tolerance: f64) -> (BRep, BooleanCleanupReport) {
 let mut report = BooleanCleanupReport::default();
 let tol = tolerance.max(TOLERANCE_ABS);

 // Step 1: Detect and remove internal faces
 let internal_faces = identify_internal_faces(brep);
 let (brep, removal_report) = remove_internal_faces(brep, &internal_faces);
 report.internal_faces_removed = removal_report.faces_removed;

 // Step 2: Merge duplicate faces
 let duplicate_report = detect_duplicate_faces(&brep, tol);
 let mut faces_to_merge: Vec<usize> = Vec::new();
 for pair in &duplicate_report.duplicate_pairs {
 if pair.opposite_orientation {
 faces_to_merge.push(pair.face_b);
 }
 }
 let (brep, merge_report) = remove_internal_faces(&brep, &faces_to_merge);
 report.duplicate_faces_merged = merge_report.faces_removed;

 // Step 3: Remove degenerate faces
 let (brep, degenerate_removed) = remove_degenerate_faces(&brep);
 report.degenerate_faces_removed = degenerate_removed;

 // Step 4: Merge close vertices
 let (brep, vertices_merged) = merge_close_vertices(&brep, tol);
 report.vertices_merged = vertices_merged;

 // Step 5: Sew close edges
 let (brep, sew_report) = sew_close_edges(&brep, tol);
 report.edges_sewn = sew_report.edges_sewn;

 // Step 6: Fix tolerances
 let brep = propagate_tolerances(&brep, tol, ToleranceFlowDirection::BottomUp);

 // Validate result
 report.is_valid = !brep.solids.is_empty();
 report.summary = format!(
 "BooleanCleanup: {} internal faces, {} duplicates merged, {} degenerate removed, {} vertices merged, {} edges sewn",
 report.internal_faces_removed,
 report.duplicate_faces_merged,
 report.degenerate_faces_removed,
 report.vertices_merged,
 report.edges_sewn
 );

 (brep, report)
}

// = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =
// Boolean Operation Type for Tolerance Propagation
// = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =

/// Type of boolean operation that was performed.
///
/// Used by tolerance propagation to apply operation-specific rules.
/// This is distinct from `builder::BooleanOpTypeForTolerance` to avoid naming conflicts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BooleanOpTypeForTolerance {
 /// Union (fuse) operation.
 #[default]
 Union,
 /// Intersection operation.
 Intersection,
 /// Difference (cut) operation.
 Difference,
 /// General boolean operation (unknown type).
 General,
}

/// Configuration for post-boolean tolerance propagation.
#[derive(Debug, Clone)]
pub struct PostBooleanToleranceConfig {
 /// Base tolerance floor for entities without explicit tolerance.
 pub tolerance_floor: f64,
 /// Multiplier applied to intersection edge tolerances.
 pub intersection_edge_factor: f64,
 /// Maximum allowed edge tolerance after propagation.
 pub max_edge_tolerance: f64,
 /// Maximum allowed face tolerance after propagation.
 pub max_face_tolerance: f64,
 /// Whether to propagate from intersection vertices to edges.
 pub propagate_vertex_to_edge: bool,
 /// Whether to propagate from edges to faces.
 pub propagate_edge_to_face: bool,
 /// Whether to detect and handle tolerance conflicts.
 pub handle_conflicts: bool,
}

impl Default for PostBooleanToleranceConfig {
 fn default() -> Self {
 Self {
 tolerance_floor: TOLERANCE_ABS,
 intersection_edge_factor: 1.0,
 max_edge_tolerance: 1.0,
 max_face_tolerance: 1.0,
 propagate_vertex_to_edge: true,
 propagate_edge_to_face: true,
 handle_conflicts: true,
 }
 }
}

impl PostBooleanToleranceConfig {
 /// Create a config for high-precision boolean operations.
 pub fn high_precision() -> Self {
 Self {
 tolerance_floor: TOLERANCE_COORD_SUB,
 intersection_edge_factor: 1.0,
 max_edge_tolerance: 0.01,
 max_face_tolerance: 0.01,
 ..Default::default()
 }
 }

 /// Create a config for standard CAD operations.
 pub fn standard() -> Self {
 Self::default()
 }

 /// Create a config for relaxed tolerance (e.g., visualization, 3D printing).
 pub fn relaxed() -> Self {
 Self {
 tolerance_floor: TOLERANCE_RETRY_LADDER_MID,
 intersection_edge_factor: 2.0,
 max_edge_tolerance: 1.0,
 max_face_tolerance: 1.0,
 ..Default::default()
 }
 }
}

/// Report from post-boolean tolerance propagation.
#[derive(Debug, Clone, Default)]
pub struct PostBooleanToleranceReport {
 /// Number of vertices whose tolerance was increased.
 pub vertices_updated: usize,
 /// Number of edges whose tolerance was increased.
 pub edges_updated: usize,
 /// Number of faces whose tolerance was increased.
 pub faces_updated: usize,
 /// Number of tolerance conflicts detected.
 pub conflicts_detected: usize,
 /// Number of tolerance conflicts resolved.
 pub conflicts_resolved: usize,
 /// Maximum vertex tolerance after propagation.
 pub max_vertex_tolerance: f64,
 /// Maximum edge tolerance after propagation.
 pub max_edge_tolerance: f64,
 /// Maximum face tolerance after propagation.
 pub max_face_tolerance: f64,
}

/// Propagate tolerances after a boolean operation.
///
/// This function applies OCCT-style tolerance propagation rules tailored to
/// the type of boolean operation performed. It handles:
///
/// 1. Intersection vertices: New vertices created at curve/surface intersections
/// receive tolerances based on the geometric precision of the intersection.
/// 2. Edge propagation: Edge tolerance >= max(vertex tolerances at endpoints).
/// 3. Face propagation: Face tolerance >= max(edge tolerances on boundary).
/// 4. Conflict resolution: Detects and resolves cases where vertex tolerance
/// exceeds edge tolerance, etc.
///
/// # Arguments
///
/// * `brep` - The BRep after boolean operation.
/// * `operation_type` - The type of boolean operation performed.
/// * `intersection_edge_indices` - Indices of edges created during intersection.
/// * `intersection_vertex_indices` - Indices of vertices created during intersection.
///
/// # Returns
///
/// A tuple of (updated BRep, propagation report).
pub fn propagate_tolerances_post_boolean_op(
 brep: &BRep,
 operation_type: BooleanOpTypeForTolerance,
 intersection_edge_indices: &[usize],
 intersection_vertex_indices: &[usize],
) -> (BRep, PostBooleanToleranceReport) {
 propagate_tolerances_post_boolean_op_with_config(
 brep,
 operation_type,
 intersection_edge_indices,
 intersection_vertex_indices,
 &PostBooleanToleranceConfig::default(),
 )
}

/// Propagate tolerances after a boolean operation with custom configuration.
pub fn propagate_tolerances_post_boolean_op_with_config(
 brep: &BRep,
 operation_type: BooleanOpTypeForTolerance,
 intersection_edge_indices: &[usize],
 intersection_vertex_indices: &[usize],
 config: &PostBooleanToleranceConfig,
) -> (BRep, PostBooleanToleranceReport) {
 let floor = config.tolerance_floor.max(TOLERANCE_ABS);
 let mut result = brep.clone();
 let mut report = PostBooleanToleranceReport::default();

 let n_verts = result.vertices.len();
 let n_edges = result.edges.len();
 let n_faces: usize = result.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();

 // Ensure tolerance arrays are sized
 if result.geom.vertex_tolerance.len() < n_verts {
 result.geom.vertex_tolerance.resize(n_verts, floor);
 }
 if result.geom.edge_tolerance.len() < n_edges {
 result.geom.edge_tolerance.resize(n_edges, floor);
 }
 if result.geom.face_tolerance.len() < n_faces {
 result.geom.face_tolerance.resize(n_faces, floor);
 }

 // Step 1: Set initial tolerances for intersection entities
 // OCCT-style: intersection edges get a tolerance based on operation type
 let base_intersection_tol = match operation_type {
 BooleanOpTypeForTolerance::Intersection => floor * 10.0,
 BooleanOpTypeForTolerance::Union => floor * 5.0,
 BooleanOpTypeForTolerance::Difference => floor * 8.0,
 BooleanOpTypeForTolerance::General => floor * 10.0,
 };

 // Apply intersection edge tolerances
 for &ei in intersection_edge_indices {
 if ei < result.geom.edge_tolerance.len() {
 let new_tol = base_intersection_tol * config.intersection_edge_factor;
 let old_tol = result.geom.edge_tolerance[ei];
 if new_tol > old_tol {
 result.geom.edge_tolerance[ei] = new_tol.min(config.max_edge_tolerance);
 report.edges_updated += 1;
 }
 }
 }

 // Apply intersection vertex tolerances
 for &vi in intersection_vertex_indices {
 if vi < result.geom.vertex_tolerance.len() {
 let new_tol = base_intersection_tol;
 let old_tol = result.geom.vertex_tolerance[vi];
 if new_tol > old_tol {
 result.geom.vertex_tolerance[vi] = new_tol;
 report.vertices_updated += 1;
 }
 }
 }

 // Step 2: Propagate vertex -> edge (OCCT BRepLib::UpdateEdgeTol rule)
 if config.propagate_vertex_to_edge {
 for ei in 0..n_edges {
 let edge = &result.edges[ei];
 let vtol_start = result.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
 let vtol_end = result.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
 let max_vtol = vtol_start.max(vtol_end);

 let cur_etol = result.geom.edge_tolerance[ei];
 let new_etol = cur_etol.max(max_vtol).min(config.max_edge_tolerance);

 if new_etol > cur_etol {
 result.geom.edge_tolerance[ei] = new_etol;
 report.edges_updated += 1;
 }
 }
 }

 // Step 3: Propagate edge -> face
 if config.propagate_edge_to_face {
 let mut flat_fi = 0usize;
 for solid in &result.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 let mut max_etol = floor;
 for we in &face.outer_wire.edges {
 if we.idx < result.geom.edge_tolerance.len() {
 max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
 }
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 if we.idx < result.geom.edge_tolerance.len() {
 max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
 }
 }
 }

 let cur_ftol = result.geom.face_tolerance.get(flat_fi).copied().unwrap_or(floor);
 let new_ftol = cur_ftol.max(max_etol).min(config.max_face_tolerance);

 if new_ftol > cur_ftol
 && flat_fi < result.geom.face_tolerance.len() {
 result.geom.face_tolerance[flat_fi] = new_ftol;
 report.faces_updated += 1;
 }
 flat_fi += 1;
 }
 }
 }
 }

 // Step 4: Detect and handle tolerance conflicts
 if config.handle_conflicts {
 let (conflicts, resolved) = detect_and_resolve_tolerance_conflicts(&mut result, floor);
 report.conflicts_detected = conflicts;
 report.conflicts_resolved = resolved;
 }

 // Compute max tolerances for report
 if !result.geom.vertex_tolerance.is_empty() {
 report.max_vertex_tolerance = result.geom.vertex_tolerance.iter()
 .cloned()
 .fold(0.0_f64, f64::max);
 }
 if !result.geom.edge_tolerance.is_empty() {
 report.max_edge_tolerance = result.geom.edge_tolerance.iter()
 .cloned()
 .fold(0.0_f64, f64::max);
 }
 if !result.geom.face_tolerance.is_empty() {
 report.max_face_tolerance = result.geom.face_tolerance.iter()
 .cloned()
 .fold(0.0_f64, f64::max);
 }

 (result, report)
}

/// Detect and resolve tolerance conflicts in a BRep.
///
/// A conflict occurs when:
/// - A vertex tolerance exceeds the tolerance of an edge it belongs to
/// - An edge tolerance exceeds the tolerance of a face it bounds
///
/// Returns (conflicts_detected, conflicts_resolved).
fn detect_and_resolve_tolerance_conflicts(brep: &mut BRep, floor: f64) -> (usize, usize) {
 let mut conflicts = 0usize;
 let mut resolved = 0usize;

 // Check vertex > edge conflicts
 for ei in 0..brep.edges.len() {
 let edge = &brep.edges[ei];
 let vtol_start = brep.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
 let vtol_end = brep.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
 let etol = brep.geom.edge_tolerance.get(ei).copied().unwrap_or(floor);

 if vtol_start > etol + TOLERANCE_FLOAT_DEDUP || vtol_end > etol + TOLERANCE_FLOAT_DEDUP {
 conflicts += 1;
 // Resolve: increase edge tolerance
 if ei < brep.geom.edge_tolerance.len() {
 let new_etol = etol.max(vtol_start).max(vtol_end);
 brep.geom.edge_tolerance[ei] = new_etol;
 resolved += 1;
 }
 }
 }

 // Check edge > face conflicts
 let mut flat_fi = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 let ftol = brep.geom.face_tolerance.get(flat_fi).copied().unwrap_or(floor);

 let mut max_etol = floor;
 let mut has_conflict = false;
 for we in &face.outer_wire.edges {
 if we.idx < brep.geom.edge_tolerance.len() {
 let etol = brep.geom.edge_tolerance[we.idx];
 max_etol = max_etol.max(etol);
 if etol > ftol + TOLERANCE_FLOAT_DEDUP {
 has_conflict = true;
 }
 }
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 if we.idx < brep.geom.edge_tolerance.len() {
 let etol = brep.geom.edge_tolerance[we.idx];
 max_etol = max_etol.max(etol);
 if etol > ftol + TOLERANCE_FLOAT_DEDUP {
 has_conflict = true;
 }
 }
 }
 }

 if has_conflict {
 conflicts += 1;
 // Resolve: increase face tolerance
 if flat_fi < brep.geom.face_tolerance.len() {
 brep.geom.face_tolerance[flat_fi] = max_etol;
 resolved += 1;
 }
 }
 flat_fi += 1;
 }
 }
 }

 (conflicts, resolved)
}

// = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =
// Post-Sew Tolerance Propagation
// = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =

/// Configuration for post-sew tolerance propagation.
#[derive(Debug, Clone)]
pub struct PostSewToleranceConfig {
 /// Base tolerance floor for entities without explicit tolerance.
 pub tolerance_floor: f64,
 /// Factor to multiply sewing tolerance by for seam edges.
 pub seam_tolerance_factor: f64,
 /// Whether to ensure consistency across sewn edges.
 pub ensure_seam_consistency: bool,
 /// Maximum allowed tolerance growth ratio.
 pub max_growth_ratio: f64,
}

impl Default for PostSewToleranceConfig {
 fn default() -> Self {
 Self {
 tolerance_floor: TOLERANCE_ABS,
 seam_tolerance_factor: 1.5,
 ensure_seam_consistency: true,
 max_growth_ratio: 100.0,
 }
 }
}

/// Report from post-sew tolerance propagation.
#[derive(Debug, Clone, Default)]
pub struct PostSewToleranceReport {
 /// Number of seam edges whose tolerance was updated.
 pub seam_edges_updated: usize,
 /// Number of faces whose tolerance was updated for seam consistency.
 pub faces_updated: usize,
 /// Maximum tolerance among seam edges.
 pub max_seam_tolerance: f64,
 /// Number of edges that required tolerance harmonization.
 pub edges_harmonized: usize,
}

/// Propagate tolerances after a sewing operation.
///
/// After sewing, edges that were joined together (seam edges) need their
/// tolerances updated to ensure geometric consistency. This function:
///
/// 1. Updates seam edge tolerances to be at least the sewing tolerance
/// 2. Ensures consistency across both sides of a seam
/// 3. Propagates tolerance updates to adjacent faces
///
/// # Arguments
///
/// * `brep` - The BRep after sewing.
/// * `sewing_tolerance` - The tolerance used during sewing.
/// * `seam_edge_pairs` - Pairs of edge indices that were sewn together.
///
/// # Returns
///
/// A tuple of (updated BRep, propagation report).
pub fn propagate_tolerances_post_sew(
 brep: &BRep,
 sewing_tolerance: f64,
 seam_edge_pairs: &[(usize, usize)],
) -> (BRep, PostSewToleranceReport) {
 propagate_tolerances_post_sew_with_config(
 brep,
 sewing_tolerance,
 seam_edge_pairs,
 &PostSewToleranceConfig::default(),
 )
}

/// Propagate tolerances after a sewing operation with custom configuration.
pub fn propagate_tolerances_post_sew_with_config(
 brep: &BRep,
 sewing_tolerance: f64,
 seam_edge_pairs: &[(usize, usize)],
 config: &PostSewToleranceConfig,
) -> (BRep, PostSewToleranceReport) {
 let floor = config.tolerance_floor.max(TOLERANCE_ABS);
 let seam_tol = sewing_tolerance.max(floor) * config.seam_tolerance_factor;

 let mut result = brep.clone();
 let mut report = PostSewToleranceReport::default();

 let n_verts = result.vertices.len();
 let n_edges = result.edges.len();
 let n_faces: usize = result.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();

 // Ensure tolerance arrays are sized
 if result.geom.vertex_tolerance.len() < n_verts {
 result.geom.vertex_tolerance.resize(n_verts, floor);
 }
 if result.geom.edge_tolerance.len() < n_edges {
 result.geom.edge_tolerance.resize(n_edges, floor);
 }
 if result.geom.face_tolerance.len() < n_faces {
 result.geom.face_tolerance.resize(n_faces, floor);
 }

 // Step 1: Harmonize seam edge tolerances
 let mut edge_tol_updates: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();

 for &(e1, e2) in seam_edge_pairs {
 let tol1 = result.geom.edge_tolerance.get(e1).copied().unwrap_or(floor);
 let tol2 = result.geom.edge_tolerance.get(e2).copied().unwrap_or(floor);
 let harmonized_tol = tol1.max(tol2).max(seam_tol);

 // Check growth ratio
 let growth = harmonized_tol / floor;
 let final_tol = if growth > config.max_growth_ratio {
 floor * config.max_growth_ratio
 } else {
 harmonized_tol
 };

 edge_tol_updates.insert(e1, edge_tol_updates.get(&e1).copied().unwrap_or(floor).max(final_tol));
 edge_tol_updates.insert(e2, edge_tol_updates.get(&e2).copied().unwrap_or(floor).max(final_tol));
 report.edges_harmonized += 1;
 }

 // Apply edge tolerance updates
 for (&ei, &new_tol) in &edge_tol_updates {
 if ei < result.geom.edge_tolerance.len() {
 let old_tol = result.geom.edge_tolerance[ei];
 if new_tol > old_tol {
 result.geom.edge_tolerance[ei] = new_tol;
 report.seam_edges_updated += 1;
 }
 }
 }

 // Step 2: Update vertex tolerances at seam endpoints
 for &(e1, e2) in seam_edge_pairs {
 if e1 < result.edges.len() && e2 < result.edges.len() {
 let edge1 = &result.edges[e1];
 let edge2 = &result.edges[e2];
 let seam_etol = edge_tol_updates.get(&e1).copied().unwrap_or(seam_tol);

 // Update vertices at seam edge endpoints
 for &vi in &[edge1.start, edge1.end, edge2.start, edge2.end] {
 if vi < result.geom.vertex_tolerance.len() {
 let old_vtol = result.geom.vertex_tolerance[vi];
 if seam_etol > old_vtol {
 result.geom.vertex_tolerance[vi] = seam_etol;
 }
 }
 }
 }
 }

 // Step 3: Ensure face tolerance consistency
 if config.ensure_seam_consistency {
 let mut flat_fi = 0usize;
 for solid in &result.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 let mut max_etol = floor;
 let mut has_seam_edge = false;

 for we in &face.outer_wire.edges {
 if we.idx < result.geom.edge_tolerance.len() {
 let etol = result.geom.edge_tolerance[we.idx];
 max_etol = max_etol.max(etol);
 if edge_tol_updates.contains_key(&we.idx) {
 has_seam_edge = true;
 }
 }
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 if we.idx < result.geom.edge_tolerance.len() {
 let etol = result.geom.edge_tolerance[we.idx];
 max_etol = max_etol.max(etol);
 if edge_tol_updates.contains_key(&we.idx) {
 has_seam_edge = true;
 }
 }
 }
 }

 if has_seam_edge {
 let old_ftol = result.geom.face_tolerance.get(flat_fi).copied().unwrap_or(floor);
 if max_etol > old_ftol
 && flat_fi < result.geom.face_tolerance.len() {
 result.geom.face_tolerance[flat_fi] = max_etol;
 report.faces_updated += 1;
 }
 }
 flat_fi += 1;
 }
 }
 }
 }

 // Compute max seam tolerance
 report.max_seam_tolerance = edge_tol_updates.values()
 .cloned()
 .fold(0.0_f64, f64::max);

 (result, report)
}

// = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =
// Tolerance Rules Engine
// = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =ㄦ = =

/// Rules for tolerance propagation.
///
/// These rules determine how tolerances propagate through the BRep topology
/// and how conflicts are resolved.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ToleranceRule {
 /// OCCT standard: vertex =edge =face propagation.
 /// Edge tolerance >= max(vertex tolerances at endpoints).
 /// Face tolerance >= max(edge tolerances on boundary).
 #[default]
 OcctStandard,

 /// Conservative: only propagate when absolutely necessary.
 /// Maintains minimum tolerances required for geometric validity.
 Conservative,

 /// Aggressive: propagate all tolerances upward.
 /// Useful for ensuring geometric operations succeed.
 Aggressive,

 /// Harmonized: ensure all connected entities have consistent tolerances.
 /// Propagates the maximum tolerance through connected topology.
 Harmonized,

 /// Bounded: propagate but cap at a maximum value.
 /// Prevents tolerances from growing unboundedly.
 Bounded,

 /// Model-scale: scale tolerances based on model bounding box.
 /// Useful for models at unusual scales.
 ModelScale,
}

/// Policy for handling tolerance conflicts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConflictResolutionPolicy {
 /// Do not modify tolerances when conflicts are detected.
 Ignore,
 /// Increase the lower-level tolerance to resolve conflicts.
 #[default]
 PropagateUp,
 /// Decrease the higher-level tolerance if safe to do so.
 ClampDown,
 /// Report conflicts but do not modify.
 ReportOnly,
}

/// Configuration for the tolerance propagation engine.
#[derive(Debug, Clone)]
pub struct TolerancePropagationConfig {
 /// Primary propagation rule to apply.
 pub rule: ToleranceRule,
 /// How to handle tolerance conflicts.
 pub conflict_policy: ConflictResolutionPolicy,
 /// Base tolerance floor.
 pub tolerance_floor: f64,
 /// Maximum allowed tolerance.
 pub max_tolerance: f64,
 /// For Bounded rule: the cap value.
 pub bound_value: f64,
 /// For ModelScale rule: the model scale factor.
 pub model_scale: f64,
 /// Number of propagation passes to run.
 pub propagation_passes: usize,
 /// Whether to validate after propagation.
 pub validate_result: bool,
}

impl Default for TolerancePropagationConfig {
 fn default() -> Self {
 Self {
 rule: ToleranceRule::OcctStandard,
 conflict_policy: ConflictResolutionPolicy::PropagateUp,
 tolerance_floor: TOLERANCE_ABS,
 max_tolerance: 1.0,
 bound_value: 0.01,
 model_scale: 1.0,
 propagation_passes: 3,
 validate_result: true,
 }
 }
}

impl TolerancePropagationConfig {
 /// Create config for OCCT-standard propagation.
 pub fn occt_standard() -> Self {
 Self::default()
 }

 /// Create config for conservative propagation.
 pub fn conservative() -> Self {
 Self {
 rule: ToleranceRule::Conservative,
 propagation_passes: 1,
 ..Default::default()
 }
 }

 /// Create config for aggressive propagation.
 pub fn aggressive() -> Self {
 Self {
 rule: ToleranceRule::Aggressive,
 propagation_passes: 5,
 ..Default::default()
 }
 }

 /// Create config for harmonized propagation.
 pub fn harmonized() -> Self {
 Self {
 rule: ToleranceRule::Harmonized,
 propagation_passes: 3,
 ..Default::default()
 }
 }

 /// Create config for bounded propagation.
 pub fn bounded(max_tol: f64) -> Self {
 Self {
 rule: ToleranceRule::Bounded,
 bound_value: max_tol,
 max_tolerance: max_tol,
 ..Default::default()
 }
 }

 /// Create config for model-scale propagation.
 pub fn model_scale(scale: f64) -> Self {
 Self {
 rule: ToleranceRule::ModelScale,
 model_scale: scale,
 tolerance_floor: TOLERANCE_ABS * scale,
 max_tolerance: 1.0 * scale,
 ..Default::default()
 }
 }
}

/// Engine for applying tolerance propagation rules.
///
/// This engine provides configurable tolerance propagation following
/// OCCT-style rules with additional customization options.
#[derive(Debug, Clone)]
pub struct TolerancePropagationEngine {
 /// Configuration for the engine.
 pub config: TolerancePropagationConfig,
}

impl Default for TolerancePropagationEngine {
 fn default() -> Self {
 Self::new()
 }
}
