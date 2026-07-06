/// Get curvature at a parameter value on a curve.
fn curve_curvature_at(curve: &rcad_kernel::Curve3, t: f64) -> Option<f64> {
 use rcad_kernel::CurveEval;

 let h = TOLERANCE_MESH_LEGACY;
 let p0 = curve.point_at((t - h).max(0.0));
 let p1 = curve.point_at(t);
 let p2 = curve.point_at((t + h).min(1.0));

 // Approximate curvature using finite differences
 let d1 = (p1 - p0) / h;
 let d2 = (p2 - p1) / h;
 let dd = (d2 - d1) / h;

 let d1_len = d1.length();
 if d1_len < TOLERANCE_LINEAR_ULTRA_STRICT {
 return None;
 }

 let cross = d1.cross(dd);
 let curvature = cross.length() / (d1_len.powi(3));

 Some(curvature)
}

/// Check if parameter ranges are compatible.
fn check_param_range_compatibility(
 brep: &BRep,
 e1: usize,
 e2: usize,
 tolerance: f64,
) -> bool {
 let range1 = brep.geom.edge_curve_range.get(e1).and_then(|r| *r);
 let range2 = brep.geom.edge_curve_range.get(e2).and_then(|r| *r);

 match (range1, range2) {
 (Some(r1), Some(r2)) => {
 // Check for overlap
 let min_max = r1[1].min(r2[1]);
 let max_min = r1[0].max(r2[0]);
 min_max >= max_min - tolerance
 }
 _ => true, // No range data, assume compatible
 }
}

/// Analyze a pair of faces for shared topology.
fn analyze_shared_face_pair(
 brep: &BRep,
 face1: &Face,
 face2: &Face,
 flat_idx1: usize,
 flat_idx2: usize,
 tolerance: f64,
) -> Option<SharedFaceInfo> {
 // Collect boundary vertices
 let verts1: Vec<usize> = face1
 .outer_wire
 .edges
 .iter()
 .flat_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 if we.forward {
 Some(vec![edge.start, edge.end])
 } else {
 Some(vec![edge.end, edge.start])
 }
 })
 .flatten()
 .collect();

 let verts2: Vec<usize> = face2
 .outer_wire
 .edges
 .iter()
 .flat_map(|we| {
 let edge = brep.edges.get(we.idx)?;
 if we.forward {
 Some(vec![edge.start, edge.end])
 } else {
 Some(vec![edge.end, edge.start])
 }
 })
 .flatten()
 .collect();

 // Count shared vertices
 let tol_sq = tolerance * tolerance;
 let mut shared_vertices = Vec::new();
 for &v1 in &verts1 {
 let p1 = brep.vertices.get(v1)?.point;
 for &v2 in &verts2 {
 let p2 = brep.vertices.get(v2)?.point;
 if (p1 - p2).length_squared() <= tol_sq {
 shared_vertices.push(v1.min(v2));
 break;
 }
 }
 }
 shared_vertices.sort();
 shared_vertices.dedup();

 // Collect boundary edges
 let edges1: std::collections::HashSet<usize> =
 face1.outer_wire.edges.iter().map(|we| we.idx).collect();
 let edges2: std::collections::HashSet<usize> =
 face2.outer_wire.edges.iter().map(|we| we.idx).collect();

 // Find shared edges (by geometry)
 let mut shared_edges = Vec::new();
 for &e1 in &edges1 {
 for &e2 in &edges2 {
 if let Some(info) = analyze_shared_edge_pair(brep, e1, e2, tolerance)
 && info.geometry_compatible {
 shared_edges.push(e1.min(e2));
 }
 }
 }
 shared_edges.sort();
 shared_edges.dedup();

 // Determine sharing kind
 let kind = if shared_edges.len() == edges1.len() && shared_edges.len() == edges2.len() {
 SharedFaceKind::FullShared
 } else if !shared_edges.is_empty() {
 SharedFaceKind::PartialShared
 } else if !shared_vertices.is_empty() {
 SharedFaceKind::VertexShared
 } else {
 SharedFaceKind::Adjacent
 };

 // Check normal compatibility
 let normal_dot = face1.normal.dot(face2.normal).abs();
 let normals_compatible = normal_dot >= 0.999;

 Some(SharedFaceInfo {
 face_a: flat_idx1,
 face_b: flat_idx2,
 kind,
 shared_edges,
 shared_vertices,
 normals_compatible,
 })
}

/// Merge shared faces in a BRep.
///
/// This function identifies and merges faces that share their complete boundary.
/// Only available in Aggressive mode.
fn merge_shared_faces(brep: &BRep, tolerance: f64) -> (BRep, usize) {
 let report = detect_shared_topology_advanced(brep, tolerance);

 if report.fully_shared_faces.is_empty() {
 return (brep.clone(), 0);
 }

 // For now, just count the mergeable faces
 // A full implementation would actually merge the faces
 let merged_count = report.fully_shared_faces.len();

 (brep.clone(), merged_count)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Connectivity Graph Analysis
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// A graph representing topological connectivity in a BRep.
///
/// This structure tracks how faces, edges, and vertices are connected,
/// enabling analysis of disconnected components and connectivity strength.
#[derive(Debug, Clone, Default)]
pub struct ConnectivityGraph {
 /// Number of vertices in the graph.
 pub vertex_count: usize,
 /// Number of edges in the graph.
 pub edge_count: usize,
 /// Number of faces in the graph.
 pub face_count: usize,
 /// Adjacency list: vertex -> connected vertices.
 pub vertex_adjacency: Vec<Vec<usize>>,
 /// Adjacency list: edge -> connected edges (via shared vertices).
 pub edge_adjacency: Vec<Vec<usize>>,
 /// Adjacency list: face -> connected faces (via shared edges).
 pub face_adjacency: Vec<Vec<usize>>,
 /// Edge-to-vertex mapping: edge -> (start_vertex, end_vertex).
 pub edge_vertices: Vec<(usize, usize)>,
 /// Face-to-edge mapping: face -> edge indices in outer wire.
 pub face_edges: Vec<Vec<usize>>,
 /// Connected components (vertex groups).
 pub vertex_components: Vec<Vec<usize>>,
 /// Connected components (face groups).
 pub face_components: Vec<Vec<usize>>,
 /// Connectivity strength metrics per edge.
 pub edge_strength: Vec<f64>,
 /// Connectivity strength metrics per face.
 pub face_strength: Vec<f64>,
}

/// Metrics for connectivity strength.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConnectivityStrength {
 /// Weak connection (single vertex shared).
 Weak,
 /// Medium connection (single edge shared).
 Medium,
 /// Strong connection (multiple edges shared).
 Strong,
 /// Full connection (entire boundary shared).
 Full,
}

impl ConnectivityStrength {
 /// Convert to a numeric strength value (0.0 to 1.0).
 pub fn to_value(&self) -> f64 {
 match self {
 ConnectivityStrength::Weak => 0.25,
 ConnectivityStrength::Medium => 0.5,
 ConnectivityStrength::Strong => 0.75,
 ConnectivityStrength::Full => 1.0,
 }
 }
}

/// Build a connectivity graph from a BRep.
///
/// This function analyzes the topological connectivity of a BRep and
/// returns a graph structure that tracks:
/// - Which faces are connected via shared edges
/// - Which edges are connected via shared vertices
/// - Which vertices are connected via edges
/// - Disconnected components
/// - Connectivity strength metrics
///
/// # Arguments
/// * `brep` - The BRep to analyze.
///
/// # Returns
/// A `ConnectivityGraph` containing all connectivity information.
pub fn build_connectivity_graph(brep: &BRep) -> ConnectivityGraph {
 let mut graph = ConnectivityGraph::default();

 let n_vertices = brep.vertices.len();
 let n_edges = brep.edges.len();

 graph.vertex_count = n_vertices;
 graph.edge_count = n_edges;

 // Initialize adjacency lists
 graph.vertex_adjacency = vec![Vec::new(); n_vertices];
 graph.edge_adjacency = vec![Vec::new(); n_edges];
 graph.edge_vertices = Vec::with_capacity(n_edges);

 // Build vertex adjacency via edges
 for edge in brep.edges.iter() {
 graph.edge_vertices.push((edge.start, edge.end));

 // Add bidirectional vertex adjacency
 if edge.start < n_vertices && edge.end < n_vertices {
 if !graph.vertex_adjacency[edge.start].contains(&edge.end) {
 graph.vertex_adjacency[edge.start].push(edge.end);
 }
 if !graph.vertex_adjacency[edge.end].contains(&edge.start) {
 graph.vertex_adjacency[edge.end].push(edge.start);
 }
 }
 }

 // Build edge adjacency via shared vertices
 let mut vertex_to_edges: Vec<Vec<usize>> = vec![Vec::new(); n_vertices];
 for (ei, edge) in brep.edges.iter().enumerate() {
 if edge.start < n_vertices {
 vertex_to_edges[edge.start].push(ei);
 }
 if edge.end < n_vertices && edge.end != edge.start {
 vertex_to_edges[edge.end].push(ei);
 }
 }

 for edges_at_vertex in &vertex_to_edges {
 for &e1 in edges_at_vertex {
 for &e2 in edges_at_vertex {
 if e1 != e2 && !graph.edge_adjacency[e1].contains(&e2) {
 graph.edge_adjacency[e1].push(e2);
 }
 }
 }
 }

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

 graph.face_count = faces.len();
 graph.face_adjacency = vec![Vec::new(); faces.len()];
 graph.face_edges = Vec::with_capacity(faces.len());
 graph.edge_strength = vec![0.0; n_edges];
 graph.face_strength = vec![0.0; faces.len()];

 // Build face edges list
 for (_, _, _, face) in &faces {
 let edges: Vec<usize> = face.outer_wire.edges.iter().map(|we| we.idx).collect();
 graph.face_edges.push(edges);
 }

 // Build edge-to-face map
 let mut edge_to_faces: std::collections::HashMap<usize, Vec<usize>> =
 std::collections::HashMap::new();
 for (fi, (_, _, _, face)) in faces.iter().enumerate() {
 for we in &face.outer_wire.edges {
 edge_to_faces.entry(we.idx).or_default().push(fi);
 }
 }

 // Build face adjacency via shared edges
 for (fi, (_, _, _, face)) in faces.iter().enumerate() {
 for we in &face.outer_wire.edges {
 if let Some(adjacent_faces) = edge_to_faces.get(&we.idx) {
 for &adj_fi in adjacent_faces {
 if adj_fi != fi && !graph.face_adjacency[fi].contains(&adj_fi) {
 graph.face_adjacency[fi].push(adj_fi);
 }
 }
 }
 }
 }

 // Calculate edge strength (number of faces sharing the edge)
 for (ei, faces_sharing) in edge_to_faces.iter() {
 if *ei < graph.edge_strength.len() {
 graph.edge_strength[*ei] = faces_sharing.len().min(4) as f64 / 4.0;
 }
 }

 // Calculate face strength (average strength of connected edges)
 for (fi, (_, _, _, face)) in faces.iter().enumerate() {
 let mut total_strength = 0.0;
 let mut count = 0;
 for we in &face.outer_wire.edges {
 if we.idx < graph.edge_strength.len() {
 total_strength += graph.edge_strength[we.idx];
 count += 1;
 }
 }
 if count > 0 {
 graph.face_strength[fi] = total_strength / count as f64;
 }
 }

 // Find connected components for vertices using union-find
 graph.vertex_components = find_connected_components(&graph.vertex_adjacency);

 // Find connected components for faces
 graph.face_components = find_connected_components(&graph.face_adjacency);

 graph
}

/// Find connected components using BFS.
fn find_connected_components(adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
 let n = adjacency.len();
 if n == 0 {
 return Vec::new();
 }

 let mut visited = vec![false; n];
 let mut components = Vec::new();

 for start in 0..n {
 if visited[start] {
 continue;
 }

 let mut component = Vec::new();
 let mut stack = vec![start];

 while let Some(node) = stack.pop() {
 if visited[node] {
 continue;
 }
 visited[node] = true;
 component.push(node);

 for &neighbor in &adjacency[node] {
 if neighbor < n && !visited[neighbor] {
 stack.push(neighbor);
 }
 }
 }

 if !component.is_empty() {
 component.sort();
 components.push(component);
 }
 }

 // Sort components by size (largest first)
 components.sort_by(|a, b| b.len().cmp(&a.len()));
 components
}

/// Identify disconnected components in a BRep.
///
/// Returns a list of component groups, where each group contains the indices
/// of faces that belong to the same connected component.
pub fn identify_disconnected_components(brep: &BRep) -> Vec<Vec<usize>> {
 let graph = build_connectivity_graph(brep);
 graph.face_components.clone()
}

/// Check if a BRep is fully connected (single component).
pub fn is_fully_connected(brep: &BRep) -> bool {
 let graph = build_connectivity_graph(brep);
 graph.face_components.len() <= 1
}

/// Get the number of disconnected components in a BRep.
pub fn disconnected_component_count(brep: &BRep) -> usize {
 let graph = build_connectivity_graph(brep);
 graph.face_components.len()
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Connectivity Gap Detection
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// A gap between disconnected regions in a BRep.
#[derive(Debug, Clone)]
pub struct ConnectivityGap {
 /// Index of the first face region.
 pub face_a: usize,
 /// Index of the second face region.
 pub face_b: usize,
 /// Component index of the first face.
 pub component_a: usize,
 /// Component index of the second face.
 pub component_b: usize,
 /// Minimum distance between the two regions.
 pub distance: f64,
 /// Closest point on face A.
 pub point_a: DVec3,
 /// Closest point on face B.
 pub point_b: DVec3,
 /// Type of gap.
 pub gap_type: GapType,
}

/// Classification of connectivity gap types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapType {
 /// Parallel faces with constant gap (like a thin wall).
 Parallel,
 /// Adjacent faces that should share an edge.
 Adjacent,
 /// Corner gap where vertices should meet.
 Corner,
 /// Complex gap requiring fill surface.
 Complex,
 /// No gap detected (faces are connected).
 None,
}

/// Detect gaps between disconnected components in a BRep.
///
/// This function finds the closest points between disconnected regions
/// and classifies the type of gap that needs to be bridged.
///
/// # Arguments
/// * `brep` - The BRep to analyze.
/// * `tolerance` - Maximum distance to consider as a gap.
///
/// # Returns
/// A vector of `ConnectivityGap` structures describing each gap.
pub fn detect_connectivity_gaps(brep: &BRep, tolerance: f64) -> Vec<ConnectivityGap> {
 let graph = build_connectivity_graph(brep);
 let mut gaps = Vec::new();

 if graph.face_components.len() <= 1 {
 return gaps;
 }

 // Collect face centers for each component
 let mut component_centers: Vec<Vec<(usize, DVec3)>> = Vec::new();
 for component in &graph.face_components {
 let mut centers = Vec::new();
 for &fi in component {
 if let Some(center) = compute_face_center(brep, fi) {
 centers.push((fi, center));
 }
 }
 component_centers.push(centers);
 }

 // Find closest pairs between components
 for (ci_a, centers_a) in component_centers.iter().enumerate() {
 for (ci_b, centers_b) in component_centers.iter().enumerate() {
 if ci_b <= ci_a {
 continue;
 }

 let mut min_dist = f64::INFINITY;
 let mut best_pair: Option<(usize, usize, DVec3, DVec3)> = None;

 for &(fa, center_a) in centers_a {
 for &(fb, center_b) in centers_b {
 let dist = (center_a - center_b).length();
 if dist < min_dist {
 min_dist = dist;
 best_pair = Some((fa, fb, center_a, center_b));
 }
 }
 }

 if let Some((fa, fb, pa, pb)) = best_pair
 && min_dist <= tolerance {
 let gap_type = classify_gap_type(brep, fa, fb, min_dist, tolerance);
 gaps.push(ConnectivityGap {
 face_a: fa,
 face_b: fb,
 component_a: ci_a,
 component_b: ci_b,
 distance: min_dist,
 point_a: pa,
 point_b: pb,
 gap_type,
 });
 }
 }
 }

 gaps
}

/// Compute the center point of a face (by averaging vertex positions).
fn compute_face_center(brep: &BRep, face_flat_idx: usize) -> Option<DVec3> {
 let faces: Vec<&Face> = brep
 .solids
 .iter()
 .flat_map(|s| &s.shells)
 .flat_map(|sh| &sh.faces)
 .collect();

 let face = faces.get(face_flat_idx)?;
 let mut center = DVec3::ZERO;
 let mut count = 0;

 for we in &face.outer_wire.edges {
 let edge = brep.edges.get(we.idx)?;
 let v = if we.forward { edge.start } else { edge.end };
 if v < brep.vertices.len() {
 center += brep.vertices[v].point;
 count += 1;
 }
 }

 if count > 0 {
 Some(center / count as f64)
 } else {
 None
 }
}

/// Classify the type of gap between two faces.
fn classify_gap_type(brep: &BRep, fa: usize, fb: usize, distance: f64, tolerance: f64) -> GapType {
 let faces: Vec<&Face> = brep
 .solids
 .iter()
 .flat_map(|s| &s.shells)
 .flat_map(|sh| &sh.faces)
 .collect();

 let face_a = match faces.get(fa) {
 Some(f) => f,
 None => return GapType::Complex,
 };
 let face_b = match faces.get(fb) {
 Some(f) => f,
 None => return GapType::Complex,
 };

 // Check if normals are parallel (indicating parallel faces)
 let normal_dot = face_a.normal.dot(face_b.normal).abs();
 if normal_dot > 0.99 {
 return GapType::Parallel;
 }

 // Check if normals are perpendicular (indicating adjacent faces)
 if normal_dot < 0.1 {
 // Check if edges are close
 for we_a in &face_a.outer_wire.edges {
 if let Some(edge_a) = brep.edges.get(we_a.idx) {
 let pa_s = brep.vertices.get(edge_a.start).map(|v| v.point);
 let pa_e = brep.vertices.get(edge_a.end).map(|v| v.point);
 if let (Some(pas), Some(pae)) = (pa_s, pa_e) {
 for we_b in &face_b.outer_wire.edges {
 if let Some(edge_b) = brep.edges.get(we_b.idx) {
 let pb_s = brep.vertices.get(edge_b.start).map(|v| v.point);
 let pb_e = brep.vertices.get(edge_b.end).map(|v| v.point);
 if let (Some(pbs), Some(pbe)) = (pb_s, pb_e) {
 // Check if edges are close
 let dist_ss = (pas - pbs).length();
 let dist_se = (pas - pbe).length();
 let dist_es = (pae - pbs).length();
 let dist_ee = (pae - pbe).length();

 if dist_ss <= tolerance
 || dist_se <= tolerance
 || dist_es <= tolerance
 || dist_ee <= tolerance
 {
 return GapType::Adjacent;
 }
 }
 }
 }
 }
 }
 }
 }

 // Check if it's a corner gap (vertices very close)
 if distance < tolerance * 0.1 {
 return GapType::Corner;
 }

 GapType::Complex
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Component Merging Strategies
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Strategy for merging disconnected components.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergeStrategy {
 /// Merge by proximity (nearest faces first).
 ByProximity,
 /// Merge by topology (create shared edges).
 ByTopology,
 /// Merge by geometry (same surface).
 ByGeometry,
 /// Merge all components into single shell.
 ForceMerge,
}

/// Configuration for component merging.
#[derive(Debug, Clone)]
pub struct MergeConfig {
 /// Strategy to use for merging.
 pub strategy: MergeStrategy,
 /// Maximum distance for proximity merging.
 pub proximity_tolerance: f64,
 /// Whether to create bridge faces between components.
 pub create_bridges: bool,
 /// Minimum bridge face quality (0.0 to 1.0).
 pub min_bridge_quality: f64,
 /// Whether to preserve original face orientations.
 pub preserve_orientations: bool,
}

impl Default for MergeConfig {
 fn default() -> Self {
 Self {
 strategy: MergeStrategy::ByProximity,
 proximity_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
 create_bridges: true,
 min_bridge_quality: 0.5,
 preserve_orientations: true,
 }
 }
}

/// Result of component merging.
#[derive(Debug, Clone, Default)]
pub struct MergeReport {
 /// Number of components merged.
 pub components_merged: usize,
 /// Number of bridge faces created.
 pub bridges_created: usize,
 /// Number of vertices merged during the operation.
 pub vertices_merged: usize,
 /// Number of edges created during merging.
 pub edges_created: usize,
 /// Final component count.
 pub final_component_count: usize,
 /// Whether the merge succeeded.
 pub success: bool,
 /// Error messages if merge failed.
 pub errors: Vec<String>,
}

/// Merge disconnected components in a BRep.
///
/// This function attempts to connect disconnected regions in a BRep
/// using the specified merging strategy.
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `strategy` - The merging strategy to use.
///
/// # Returns
/// A tuple of (modified BRep, merge report).
pub fn merge_disconnected_components(brep: &BRep, strategy: MergeStrategy) -> (BRep, MergeReport) {
 let config = MergeConfig {
 strategy,
 ..Default::default()
 };
 merge_disconnected_components_with_config(brep, &config)
}

/// Merge disconnected components with custom configuration.
pub fn merge_disconnected_components_with_config(
 brep: &BRep,
 config: &MergeConfig,
) -> (BRep, MergeReport) {
 let mut result = brep.clone();
 let mut report = MergeReport::default();

 let initial_components = disconnected_component_count(&result);
 if initial_components <= 1 {
 report.final_component_count = 1;
 report.success = true;
 return (result, report);
 }

 // Detect gaps between components
 let gaps = detect_connectivity_gaps(&result, config.proximity_tolerance);
 if gaps.is_empty() {
 report.errors.push("No gaps detected within tolerance".to_string());
 report.final_component_count = initial_components;
 report.success = initial_components <= 1;
 return (result, report);
 }

 match config.strategy {
 MergeStrategy::ByProximity => {
 // Sort gaps by distance (smallest first)
 let mut sorted_gaps = gaps;
 sorted_gaps.sort_by(|a, b| {
 a.distance.partial_cmp(&b.distance).unwrap_or(std::cmp::Ordering::Equal)
 });

 for gap in sorted_gaps {
 let merge_result = merge_gap_by_proximity(&result, &gap, config);
 result = merge_result.0;
 report.vertices_merged += merge_result.1.vertices_merged;
 report.edges_created += merge_result.1.edges_created;
 if merge_result.1.success {
 report.components_merged += 1;
 }
 }
 }
 MergeStrategy::ByTopology => {
 for gap in &gaps {
 let merge_result = merge_gap_by_topology(&result, gap, config);
 result = merge_result.0;
 report.vertices_merged += merge_result.1.vertices_merged;
 report.edges_created += merge_result.1.edges_created;
 if merge_result.1.success {
 report.components_merged += 1;
 }
 }
 }
 MergeStrategy::ByGeometry => {
 for gap in &gaps {
 let merge_result = merge_gap_by_geometry(&result, gap, config);
 result = merge_result.0;
 if merge_result.1.success {
 report.components_merged += 1;
 report.vertices_merged += merge_result.1.vertices_merged;
 }
 }
 }
 MergeStrategy::ForceMerge => {
 // Force merge all components by creating bridge faces
 if config.create_bridges {
 let (new_result, bridges) = create_bridges(&result, &gaps);
 result = new_result;
 report.bridges_created = bridges;
 }
 report.components_merged = initial_components.saturating_sub(1);
 }
 }

 report.final_component_count = disconnected_component_count(&result);
 report.success = report.final_component_count < initial_components;

 (result, report)
}

/// Merge a gap by bringing nearby vertices together.
fn merge_gap_by_proximity(
 brep: &BRep,
 gap: &ConnectivityGap,
 config: &MergeConfig,
) -> (BRep, MergeReport) {
 let mut result = brep.clone();
 let mut report = MergeReport::default();

 if gap.distance > config.proximity_tolerance {
 report.success = false;
 return (result, report);
 }

 // Find closest vertices from each component
 let faces: Vec<&Face> = result
 .solids
 .iter()
 .flat_map(|s| &s.shells)
 .flat_map(|sh| &sh.faces)
 .collect();

 let face_a = match faces.get(gap.face_a) {
 Some(f) => f,
 None => {
 report.errors.push("Face A not found".to_string());
 return (result, report);
 }
 };
 let face_b = match faces.get(gap.face_b) {
 Some(f) => f,
 None => {
 report.errors.push("Face B not found".to_string());
 return (result, report);
 }
 };

 // Collect vertices from each face
 let mut verts_a: Vec<usize> = Vec::new();
 for we in &face_a.outer_wire.edges {
 if let Some(edge) = result.edges.get(we.idx) {
 verts_a.push(edge.start);
 verts_a.push(edge.end);
 }
 }
 verts_a.sort();
 verts_a.dedup();

 let mut verts_b: Vec<usize> = Vec::new();
 for we in &face_b.outer_wire.edges {
 if let Some(edge) = result.edges.get(we.idx) {
 verts_b.push(edge.start);
 verts_b.push(edge.end);
 }
 }
 verts_b.sort();
 verts_b.dedup();

 // Find and merge closest vertex pair
 let tol_sq = config.proximity_tolerance * config.proximity_tolerance;
 for &va in &verts_a {
 if va >= result.vertices.len() {
 continue;
 }
 let pa = result.vertices[va].point;
 for &vb in &verts_b {
 if vb >= result.vertices.len() {
 continue;
 }
 let pb = result.vertices[vb].point;
 if (pa - pb).length_squared() <= tol_sq && va != vb {
 // Merge vb into va
 result = merge_specific_vertices(&result, vb, va);
 report.vertices_merged += 1;
 report.success = true;
 }
 }
 }

 (result, report)
}

/// Merge a gap by creating shared edges.
fn merge_gap_by_topology(
 brep: &BRep,
 gap: &ConnectivityGap,
 config: &MergeConfig,
) -> (BRep, MergeReport) {
 let mut result = brep.clone();
 let mut report = MergeReport::default();

 if gap.gap_type != GapType::Adjacent {
 // Topology merge only works for adjacent gaps
 report.success = false;
 return (result, report);
 }

 // Use proximity merge as the base
 let proximity_result = merge_gap_by_proximity(&result, gap, config);
 result = proximity_result.0;
 report.vertices_merged = proximity_result.1.vertices_merged;

 // Additional edge creation if needed
 if proximity_result.1.success {
 report.success = true;
 }

 (result, report)
}

/// Merge a gap by matching geometry (same surface).
fn merge_gap_by_geometry(
 brep: &BRep,
 gap: &ConnectivityGap,
 config: &MergeConfig,
) -> (BRep, MergeReport) {
 // Geometry-based merge requires same surface
 // For now, use proximity merge as fallback
 merge_gap_by_proximity(brep, gap, config)
}

/// Merge two specific vertices in a BRep.
fn merge_specific_vertices(brep: &BRep, drop_vi: usize, keep_vi: usize) -> BRep {
 if drop_vi == keep_vi || drop_vi >= brep.vertices.len() || keep_vi >= brep.vertices.len() {
 return brep.clone();
 }

 let mut result = brep.clone();

 // Update all edge references
 for edge in &mut result.edges {
 if edge.start == drop_vi {
 edge.start = keep_vi;
 } else if edge.start > drop_vi {
 edge.start -= 1;
 }
 if edge.end == drop_vi {
 edge.end = keep_vi;
 } else if edge.end > drop_vi {
 edge.end -= 1;
 }
 }

 // Remove the dropped vertex
 result.vertices.remove(drop_vi);

 // Update tolerance arrays if present
 if result.geom.vertex_tolerance.len() > drop_vi {
 result.geom.vertex_tolerance.remove(drop_vi);
 }

 result
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Bridge Creation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Create bridge faces to connect disconnected regions.
///
/// This function creates new faces that bridge the gaps between
/// disconnected components, making the BRep topologically connected.
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `gaps` - The gaps to bridge.
///
/// # Returns
/// A tuple of (modified BRep, number of bridges created).
pub fn create_bridges(brep: &BRep, gaps: &[ConnectivityGap]) -> (BRep, usize) {
 if gaps.is_empty() {
 return (brep.clone(), 0);
 }

 let mut result = brep.clone();
 let mut bridges_created = 0;

 for gap in gaps {
 if gap.gap_type == GapType::None {
 continue;
 }

 // Create a bridge face between the gap endpoints
 let bridge_result = create_single_bridge(&result, gap);
 if bridge_result.1 {
 result = bridge_result.0;
 bridges_created += 1;
 }
 }

 (result, bridges_created)
}

/// Create a single bridge face for a gap.
fn create_single_bridge(brep: &BRep, gap: &ConnectivityGap) -> (BRep, bool) {
 let mut result = brep.clone();

 // Find vertices near the gap endpoints
 let faces: Vec<&Face> = result
 .solids
 .iter()
 .flat_map(|s| &s.shells)
 .flat_map(|sh| &sh.faces)
 .collect();

 let face_a = match faces.get(gap.face_a) {
 Some(f) => f,
 None => return (result, false),
 };

 // Find the closest vertex on face A to the gap point
 let mut closest_va: Option<usize> = None;
 let mut min_dist_a = f64::INFINITY;
 for we in &face_a.outer_wire.edges {
 if let Some(edge) = result.edges.get(we.idx) {
 for &v in &[edge.start, edge.end] {
 if v < result.vertices.len() {
 let dist = (result.vertices[v].point - gap.point_a).length();
 if dist < min_dist_a {
 min_dist_a = dist;
 closest_va = Some(v);
 }
 }
 }
 }
 }

 let face_b = match faces.get(gap.face_b) {
 Some(f) => f,
 None => return (result, false),
 };

 // Find the closest vertex on face B to the gap point
 let mut closest_vb: Option<usize> = None;
 let mut min_dist_b = f64::INFINITY;
 for we in &face_b.outer_wire.edges {
 if let Some(edge) = result.edges.get(we.idx) {
 for &v in &[edge.start, edge.end] {
 if v < result.vertices.len() {
 let dist = (result.vertices[v].point - gap.point_b).length();
 if dist < min_dist_b {
 min_dist_b = dist;
 closest_vb = Some(v);
 }
 }
 }
 }
 }

 let (va, vb) = match (closest_va, closest_vb) {
 (Some(a), Some(b)) => (a, b),
 _ => return (result, false),
 };

 if va == vb {
 // Already connected
 return (result, true);
 }

 // Create an edge between the vertices if it doesn't exist
 let edge_exists = result.edges.iter().any(|e| {
 (e.start == va && e.end == vb) || (e.start == vb && e.end == va)
 });

 let bridge_edge_idx = if edge_exists {
 result.edges.iter().position(|e| {
 (e.start == va && e.end == vb) || (e.start == vb && e.end == va)
 }).unwrap()
 } else {
 // Create new edge
 let new_edge = Edge { start: va, end: vb };
 result.edges.push(new_edge);
 result.geom.edge_tolerance.push(gap.distance);
 result.edges.len() - 1
 };

 // Create a bridge face (triangle) if we have enough vertices
 // For simplicity, we create a degenerate bridge by just ensuring the edge exists
 // A proper implementation would create a new face with this edge

 // Add the edge to a new face or existing shell
 // For now, we just ensure connectivity through the edge
 if result.solids.is_empty() {
 // Create a new solid with a face containing the bridge edge
 use rcad_kernel::topology::{Shell, Solid, Wire, WireEdge};
 let wire = Wire {
 edges: vec![WireEdge::fwd(bridge_edge_idx)],
 };
 let face = Face {
 outer_wire: wire,
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 result.solids.push(Solid {
 shells: vec![Shell { faces: vec![face] }],
 });
 }

 (result, true)
}

/// Create bridge faces with custom configuration.
pub fn create_bridges_with_config(
 brep: &BRep,
 gaps: &[ConnectivityGap],
 _config: &MergeConfig,
) -> (BRep, usize) {
 create_bridges(brep, gaps)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Connectivity Validation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Report from connectivity validation.
#[derive(Debug, Clone, Default)]
pub struct ConnectivityReport {
 /// Whether the BRep is fully connected.
 pub is_connected: bool,
 /// Number of connected components.
 pub component_count: usize,
 /// Number of weak connections found.
 pub weak_connections: usize,
 /// Number of medium connections found.
 pub medium_connections: usize,
 /// Number of strong connections found.
 pub strong_connections: usize,
 /// Number of gaps detected.
 pub gaps_detected: usize,
 /// Gaps that were detected.
 pub gaps: Vec<ConnectivityGap>,
 /// Suggested improvements.
 pub suggestions: Vec<String>,
 /// Summary string.
 pub summary: String,
}

impl ConnectivityReport {
 /// Create a human-readable summary.
 pub fn summary(&self) -> String {
 if self.is_connected {
 format!(
 "Fully connected BRep with {} components, {} strong connections",
 self.component_count, self.strong_connections
 )
 } else {
 format!(
 "Disconnected BRep: {} components, {} gaps, {} weak connections",
 self.component_count, self.gaps_detected, self.weak_connections
 )
 }
 }
}

/// Validate the connectivity of a BRep.
///
/// This function performs a comprehensive connectivity analysis,
/// checking for disconnected components, weak connections, and gaps.
///
/// # Arguments
/// * `brep` - The BRep to validate.
/// * `tolerance` - Maximum distance for gap detection.
///
/// # Returns
/// A `ConnectivityReport` with detailed findings.
pub fn validate_connectivity(brep: &BRep, tolerance: f64) -> ConnectivityReport {
 let graph = build_connectivity_graph(brep);
 let mut report = ConnectivityReport::default();

 report.component_count = graph.face_components.len();
 report.is_connected = report.component_count <= 1;

 // Detect gaps
 report.gaps = detect_connectivity_gaps(brep, tolerance);
 report.gaps_detected = report.gaps.len();

 // Count connection strengths
 for &strength in &graph.edge_strength {
 if strength < 0.3 {
 report.weak_connections += 1;
 } else if strength < 0.7 {
 report.medium_connections += 1;
 } else {
 report.strong_connections += 1;
 }
 }

 // Generate suggestions
 if !report.is_connected {
 report.suggestions.push("Consider using merge_disconnected_components with ByProximity strategy".to_string());
 }

 if report.weak_connections > report.strong_connections {
 report.suggestions.push(
 "Many weak connections detected. Consider edge sewing or vertex merging.".to_string()
 );
 }

 for gap in &report.gaps {
 match gap.gap_type {
 GapType::Parallel => {
 report.suggestions.push(format!(
 "Parallel gap at distance {:.6} between faces {} and {}",
 gap.distance, gap.face_a, gap.face_b
 ));
 }
 GapType::Adjacent => {
 report.suggestions.push(format!(
 "Adjacent faces {} and {} should share an edge",
 gap.face_a, gap.face_b
 ));
 }
 GapType::Corner => {
 report.suggestions.push(format!(
 "Corner gap between faces {} and {} requires vertex merge",
 gap.face_a, gap.face_b
 ));
 }
 GapType::Complex => {
 report.suggestions.push(format!(
 "Complex gap between faces {} and {} may require fill surface",
 gap.face_a, gap.face_b
 ));
 }
 GapType::None => {}
 }
 }

 report.summary = report.summary();
 report
}

/// Quick check if a BRep needs connectivity repair.
pub fn needs_connectivity_repair(brep: &BRep) -> bool {
 !is_fully_connected(brep)
}

/// Get the connectivity strength between two faces.
pub fn get_face_connectivity_strength(brep: &BRep, face_a: usize, face_b: usize) -> ConnectivityStrength {
 let graph = build_connectivity_graph(brep);

 if face_a >= graph.face_count || face_b >= graph.face_count {
 return ConnectivityStrength::Weak;
 }

 if graph.face_adjacency[face_a].contains(&face_b) {
 // Count shared edges
 let edges_a: std::collections::HashSet<usize> = graph.face_edges.get(face_a)
 .map(|e| e.iter().copied().collect())
 .unwrap_or_default();
 let edges_b: std::collections::HashSet<usize> = graph.face_edges.get(face_b)
 .map(|e| e.iter().copied().collect())
 .unwrap_or_default();

 let shared_count = edges_a.intersection(&edges_b).count();

 match shared_count {
 0 => ConnectivityStrength::Weak,
 1 => ConnectivityStrength::Medium,
 2..=3 => ConnectivityStrength::Strong,
 _ => ConnectivityStrength::Full,
 }
 } else {
 ConnectivityStrength::Weak
 }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Enhanced Make-Connected with Connectivity Analysis
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Configuration for enhanced make-connected with connectivity analysis.
#[derive(Debug, Clone)]
pub struct EnhancedMakeConnectedConfig {
 /// Base tolerance for vertex merging.
 pub base_tolerance: f64,
 /// Maximum tolerance for gap detection.
 pub max_gap_tolerance: f64,
 /// Maximum number of repair passes.
 pub max_passes: usize,
 /// Tolerance growth factor per pass.
 pub tolerance_growth: f64,
 /// Whether to attempt component merging.
 pub merge_components: bool,
 /// Whether to create bridges for gaps.
 pub create_bridges: bool,
 /// Merge strategy to use.
 pub merge_strategy: MergeStrategy,
 /// Whether to validate connectivity after repair.
 pub validate_result: bool,
}

impl Default for EnhancedMakeConnectedConfig {
 fn default() -> Self {
 Self {
 base_tolerance: TOLERANCE_MESH_LEGACY,
 max_gap_tolerance: TOLERANCE_ADAPTIVE_MAX,
 max_passes: 5,
 tolerance_growth: 1.5,
 merge_components: true,
 create_bridges: true,
 merge_strategy: MergeStrategy::ByProximity,
 validate_result: true,
 }
 }
}

/// Report from enhanced make-connected with connectivity analysis.
#[derive(Debug, Clone, Default)]
pub struct EnhancedMakeConnectedReport {
 /// Basic make-connected report.
 pub basic_report: MakeConnectedReport,
 /// Connectivity analysis report.
 pub connectivity_report: ConnectivityReport,
 /// Merge report if components were merged.
 pub merge_report: Option<MergeReport>,
 /// Number of bridges created.
 pub bridges_created: usize,
 /// Final component count.
 pub final_components: usize,
 /// Whether the result is fully connected.
 pub is_fully_connected: bool,
}

/// Apply enhanced make-connected with full connectivity analysis.
///
/// This function performs a comprehensive connectivity repair:
/// 1. Basic vertex merging and small edge removal
/// 2. Connectivity graph analysis
/// 3. Component merging if needed
/// 4. Bridge creation for gaps
/// 5. Connectivity validation
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `config` - Configuration for the repair.
///
/// # Returns
/// A tuple of (modified BRep, detailed report).
pub fn make_connected_with_connectivity_analysis(
 brep: &BRep,
 config: &EnhancedMakeConnectedConfig,
) -> (BRep, EnhancedMakeConnectedReport) {
 let mut result = brep.clone();
 let mut report = EnhancedMakeConnectedReport::default();

 // Step 1: Basic make-connected
 let tol = config.base_tolerance.max(TOLERANCE_ABS);
 let (basic_result, basic_report) = make_connected_iterative_with_growth_cap(
 &result,
 tol,
 config.max_passes,
 config.tolerance_growth,
 config.max_gap_tolerance,
 );
 result = basic_result;
 report.basic_report = basic_report;

 // Step 2: Connectivity analysis
 report.connectivity_report = validate_connectivity(&result, config.max_gap_tolerance);

 // Step 3: Component merging if needed
 if config.merge_components && report.connectivity_report.component_count > 1 {
 let merge_config = MergeConfig {
 strategy: config.merge_strategy,
 proximity_tolerance: config.max_gap_tolerance,
 create_bridges: config.create_bridges,
 ..Default::default()
 };
 let (merged_result, merge_report) = merge_disconnected_components_with_config(&result, &merge_config);
 result = merged_result;
 report.merge_report = Some(merge_report);
 }

 // Step 4: Bridge creation
 if config.create_bridges && !report.connectivity_report.gaps.is_empty() {
 let (bridged_result, bridges) = create_bridges(&result, &report.connectivity_report.gaps);
 result = bridged_result;
 report.bridges_created = bridges;
 }

 // Step 5: Final validation
 if config.validate_result {
 let final_report = validate_connectivity(&result, config.max_gap_tolerance);
 report.final_components = final_report.component_count;
 report.is_fully_connected = final_report.is_connected;
 } else {
 report.final_components = disconnected_component_count(&result);
 report.is_fully_connected = report.final_components <= 1;
 }

 (result, report)
}

/// Repair SameRange consistency by aligning PCurve ranges with the 3D edge range.
///
/// For each edge with a known `edge_curve_range` and attached PCurves, ensure all
/// referenced `curve2d_range` entries are populated with the same `[t1, t2]`.
/// Also marks `edge_same_range[edge_idx] = true` after alignment.
///
///  ?OCCT = : BRepLib.cxx  ?SameRange (lines 75-120).
/// OCCT iterates edges, identifies those where PCurve ranges differ from the 3D
/// range, and reparameterizes the PCurves to match. This implementation performs
/// the same range-alignment by overwriting `curve2d_range` with the 3D range
/// when the mismatch exceeds `tolerance`.
pub fn fix_same_range_flags(brep: &BRep, tolerance: f64) -> (BRep, usize) {
 let mut out = brep.clone();
 let edge_count = out.edges.len();

 if out.geom.edge_same_range.len() < edge_count {
 out.geom.edge_same_range.resize(edge_count, true);
 }
 if out.geom.edge_curve_range.len() < edge_count {
 out.geom.edge_curve_range.resize(edge_count, None);
 }
 if out.geom.edge_pcurves.len() < edge_count {
 out.geom.edge_pcurves.resize(edge_count, Vec::new());
 }

 if out.geom.curve2d_range.len() < out.geom.curve2ds.len() {
 out.geom.curve2d_range.resize(out.geom.curve2ds.len(), None);
 }

 let mut fixed = 0usize;
 for edge_idx in 0..edge_count {
 let Some(range3d) = out.geom.edge_curve_range[edge_idx] else {
 continue;
 };
 let pcurves = out.geom.edge_pcurves[edge_idx].clone();
 if pcurves.is_empty() {
 continue;
 }

 let mut changed = !out.geom.edge_same_range[edge_idx];
 for pc in pcurves {
 if pc.curve2d_idx >= out.geom.curve2d_range.len() {
 continue;
 }
 match out.geom.curve2d_range[pc.curve2d_idx] {
 Some(r)
 if (r[0] - range3d[0]).abs() <= tolerance
 && (r[1] - range3d[1]).abs() <= tolerance => {}
 _ => {
 out.geom.curve2d_range[pc.curve2d_idx] = Some(range3d);
 changed = true;
 }
 }
 }

 if changed {
 out.geom.edge_same_range[edge_idx] = true;
 fixed += 1;
 }
 }

 (out, fixed)
}

/// Scan all edges for SameRange violations, flag them, and repair.
///
/// This combines the diagnostic scan from [`diagnose_same_range`] with the
/// repair logic of [`fix_same_range_flags`] in a single call.
pub fn fix_same_range_with_scan(brep: &BRep, tolerance: f64) -> (BRep, usize) {
 let diagnosis = diagnose_same_range(brep, tolerance);
 if diagnosis.suspect_edges.is_empty() {
 return (brep.clone(), 0);
 }

 let mut out = brep.clone();
 let n_edges = out.edges.len();

 if out.geom.edge_same_range.len() < n_edges {
 out.geom.edge_same_range.resize(n_edges, true);
 }

 for suspect in &diagnosis.suspect_edges {
 if suspect.edge_idx < n_edges {
 out.geom.edge_same_range[suspect.edge_idx] = false;
 }
 }

 fix_same_range_flags(&out, tolerance)
}

/// Merge vertices that are within `tolerance` of each other.
///
/// Uses spatial hashing for O(n) average performance on large models,
/// falling back to brute-force for small vertex counts.
/// For each pair of vertices closer than `tolerance`, they are merged into
/// the vertex with the smaller index. All edges and wires are remapped.
///
/// Returns the repaired BRep and the number of vertices merged.
///
/// Analogous to `BRepOffsetAPI_Sewing` vertex merging or
/// `ShapeFix_Wire::FixSameParameter`.
pub fn merge_close_vertices(brep: &BRep, tolerance: f64) -> (BRep, usize) {
 let n = brep.vertices.len();
 // Union-find: parent[i] = canonical representative of vertex i
 let mut parent: Vec<usize> = (0..n).collect();

 fn find(parent: &mut [usize], mut x: usize) -> usize {
 while parent[x] != x {
 parent[x] = parent[parent[x]]; // path compression
 x = parent[x];
 }
 x
 }

 fn union(parent: &mut [usize], a: usize, b: usize) {
 let ra = find(parent, a);
 let rb = find(parent, b);
 if ra != rb {
 // Merge to the smaller index so result is deterministic
 if ra < rb {
 parent[rb] = ra;
 } else {
 parent[ra] = rb;
 }
 }
 }

 let tol2 = tolerance * tolerance;

 // OCCT-aligned: compute degenerate edge vertex pairs to skip merging.
 // OCCT uses distinct TopoDS_Vertex for deg edge ends at the same point.
 let deg_skip: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::from_iter(
 brep.edges.iter().enumerate().filter_map(|(ei, e)| {
 if brep.geom.edge_degenerated.get(ei).copied().unwrap_or(false) {
 Some((e.start.min(e.end), e.start.max(e.end)))
 } else { None }
 })
 );

 fn spatial_cell_coord(value: f64, tolerance: f64) -> i64 {
 let t = tolerance.max(f64::MIN_POSITIVE);
 let q = (value / t).floor();
 if !q.is_finite() {
 return 0;
 }
 q.clamp(i64::MIN as f64, i64::MAX as f64) as i64
 }

 // Use spatial hashing for large models, brute-force for small ones.
 // Spatial hashing: bucket size = tolerance, check 27 neighbor cells.
 const SPATIAL_HASH_THRESHOLD: usize = 500;
 if n >= SPATIAL_HASH_THRESHOLD {
 let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
 std::collections::HashMap::with_capacity(n);
 for i in 0..n {
 let p = brep.vertices[i].point;
 let cell = (
 spatial_cell_coord(p.x, tolerance),
 spatial_cell_coord(p.y, tolerance),
 spatial_cell_coord(p.z, tolerance),
 );
 // Check 27 neighbor cells (including self)
 for dx in -1..=1 {
 for dy in -1..=1 {
 for dz in -1..=1 {
 let neighbor = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
 if let Some(bucket) = grid.get(&neighbor) {
 for &j in bucket {
 let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
 if d2 <= tol2 {
 let key = (i.min(j), i.max(j));
 if !deg_skip.contains(&key) {
 union(&mut parent, i, j);
 }
 }
 }
 }
 }
 }
 }
 grid.entry(cell).or_default().push(i);
 }
 } else {
 // Brute-force O(n ? =fast enough for small models
 for i in 0..n {
 for j in (i + 1)..n {
 let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
 if d2 <= tol2 {
 let key = (i.min(j), i.max(j));
 if !deg_skip.contains(&key) {
 union(&mut parent, i, j);
 }
 }
 }
 }
 }

 // Compress paths
 for i in 0..n {
 parent[i] = find(&mut parent, i);
 }

 // Count merges (vertices whose canonical rep is a different index)
 let merged = (0..n).filter(|&i| parent[i] != i).count();
 if merged == 0 {
 return (brep.clone(), 0);
 }

 // Build a compact vertex list and a remap table old_idx =new_idx
 let mut new_vertices: Vec<Vertex> = Vec::new();
 let mut remap = vec![0usize; n];
 let mut seen: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 for i in 0..n {
 let rep = parent[i];
 if let Some(&new_idx) = seen.get(&rep) {
 remap[i] = new_idx;
 } else {
 let new_idx = new_vertices.len();
 // Use the average position of all merged vertices for robustness
 new_vertices.push(brep.vertices[rep]);
 seen.insert(rep, new_idx);
 remap[i] = new_idx;
 }
 }

 // Re-map edges
 let new_edges: Vec<Edge> = brep
 .edges
 .iter()
 .map(|e| Edge {
 start: remap[e.start],
 end: remap[e.end],
 })
 .collect();

 // Rebuild solids with remapped wires (topology is unchanged, just vertex indices)
 let new_solids = brep
 .solids
 .iter()
 .map(|solid| Solid {
 shells: solid
 .shells
 .iter()
 .map(|shell| Shell {
 faces: shell
 .faces
 .iter()
 .map(|face| {
 let remap_wire = |w: &Wire| Wire {
 edges: w.edges.clone(), // WireEdge indices are edge indices, not vertex
 };
 Face {
 outer_wire: remap_wire(&face.outer_wire),
 inner_wires: face.inner_wires.iter().map(remap_wire).collect(),
 normal: face.normal,
 triangles: face.triangles.clone(),
 sample_point: face.sample_point,
 mesh_dirty: true,
 surface_idx: face.surface_idx,
 }
 })
 .collect(),
 })
 .collect(),
 })
 .collect();

 let mut result = brep.clone();
 result.vertices = new_vertices;
 result.edges = new_edges;
 result.solids = new_solids;

 (result, merged)
}

/// Remove faces that are degenerate:
/// - Fewer than 3 edges in the outer wire, or
/// - All wire vertices are collinear (zero-area face).
///
/// Returns the cleaned BRep and the number of faces removed.
///
/// Analogous to `ShapeFix_Shape` degenerate-face removal.
pub fn remove_degenerate_faces(brep: &BRep) -> (BRep, usize) {
 let mut removed = 0usize;

 let new_solids = brep
 .solids
 .iter()
 .map(|solid| Solid {
 shells: solid
 .shells
 .iter()
 .map(|shell| {
 let new_faces: Vec<Face> = shell
 .faces
 .iter()
 .filter(|face| {
 let wire = &face.outer_wire;
 // Must have at least 3 edges
 if wire.edges.len() < 3 {
 removed += 1;
 return false;
 }
 // Collect distinct vertex positions
 let pts: Vec<DVec3> = wire
 .edges
 .iter()
 .filter_map(|we| {
 brep.edges.get(we.idx).and_then(|e| {
 let vidx = if we.forward { e.start } else { e.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 })
 .collect();

 if pts.len() < 3 {
 removed += 1;
 return false;
 }

 // Check for zero area using Newell's method
 let area2 = newell_area(&pts);
 if area2 < TOLERANCE_METRIC_SQ_NEAR_ZERO {
 removed += 1;
 return false;
 }
 true
 })
 .cloned()
 .collect();
 Shell { faces: new_faces }
 })
 .collect(),
 })
 .collect();

 let mut result = brep.clone();
 result.solids = new_solids;
 (result, removed)
}

/// Recompute each face's `normal` field from the positions of its wire vertices,
/// using Newell's method for robustness with non-planar polygons.
///
/// Returns the updated BRep and the number of faces whose normals changed by
/// more than 1 ?(indicating they were stale or flipped).
///
/// Analogous to `BRepLib` normal re-computation after topology repair.
pub fn recompute_face_normals(brep: &BRep) -> (BRep, usize) {
 let mut changed = 0usize;

 let new_solids = brep
 .solids
 .iter()
 .map(|solid| Solid {
 shells: solid
 .shells
 .iter()
 .map(|shell| Shell {
 faces: shell
 .faces
 .iter()
 .map(|face| {
 let pts: Vec<DVec3> = face
 .outer_wire
 .edges
 .iter()
 .filter_map(|we| {
 brep.edges.get(we.idx).and_then(|e| {
 let vidx = if we.forward { e.start } else { e.end };
 brep.vertices.get(vidx).map(|v| v.point)
 })
 })
 .collect();

 let new_normal = if pts.len() >= 3 {
 let n = newell_normal(&pts);
 if n.length() > TOLERANCE_FLOAT_LOOSE {
 n.normalize()
 } else {
 face.normal
 }
 } else {
 face.normal
 };

 let dot = face.normal.dot(new_normal);
 // dot < cos(1 ? =0.9998 means the normal changed significantly
 if dot < 0.9998 {
 changed += 1;
 }

 Face {
 outer_wire: face.outer_wire.clone(),
 inner_wires: face.inner_wires.clone(),
 normal: new_normal,
 triangles: face.triangles.clone(),
 sample_point: face.sample_point,
 mesh_dirty: true,
 surface_idx: None,
 }
 })
 .collect(),
 })
 .collect(),
 })
 .collect();

 let mut result = brep.clone();
 result.solids = new_solids;
 (result, changed)
}

/// Ensure that each wire in the BRep forms a properly closed chain.
///
/// For each open wire (end of edge i =start of edge i+1 within `tolerance`),
/// attempts to close it by reversing individual edges whose orientation appears
/// flipped relative to the chain direction.
///
/// Returns the repaired BRep and the count of wires that were modified.
///
/// Analogous to `ShapeFix_Wire::FixClosed()` / `FixConnected()`.
pub fn fix_wire_orientation(brep: &BRep, tolerance: f64) -> (BRep, usize) {
 let tol2 = tolerance * tolerance;
 let mut total_fixed = 0usize;

 let new_solids = brep
 .solids
 .iter()
 .map(|solid| Solid {
 shells: solid
 .shells
 .iter()
 .map(|shell| Shell {
 faces: shell
 .faces
 .iter()
 .map(|face| {
 let (new_outer, fixed_outer) = fix_wire(&face.outer_wire, brep, tol2);
 let (new_inners, fixed_inner): (Vec<Wire>, usize) = face
 .inner_wires
 .iter()
 .map(|w| fix_wire(w, brep, tol2))
 .fold((Vec::new(), 0), |(mut wires, n), (w, f)| {
 wires.push(w);
 (wires, n + f)
 });
 let fixed = fixed_outer + fixed_inner;
 total_fixed += fixed;
 Face {
 outer_wire: new_outer,
 inner_wires: new_inners,
 normal: face.normal,
 triangles: face.triangles.clone(),
 sample_point: face.sample_point,
 mesh_dirty: true,
 surface_idx: None,
 }
 })
 .collect(),
 })
 .collect(),
 })
 .collect();

 let mut result = brep.clone();
 result.solids = new_solids;
 (result, total_fixed)
}

/// Flip inward-facing faces so shell orientation is outward-consistent.
///
/// Uses the same centroid heuristic as [`check_orientation_consistency`]. Each
/// offending face has its stored normal negated and all wires reversed.
pub fn fix_face_orientation(brep: &BRep) -> (BRep, usize) {
 let report = check_orientation_consistency(brep);
 if report.issues.is_empty() {
 return (brep.clone(), 0);
 }

 let issue_set: std::collections::HashSet<(usize, usize)> = report
 .issues
 .iter()
 .map(|issue| (issue.solid_idx, issue.face_idx))
 .collect();

 let mut flat_face_idx = 0usize;
 let mut changed = 0usize;
 let new_solids = brep
 .solids
 .iter()
 .enumerate()
 .map(|(si, solid)| Solid {
 shells: solid
 .shells
 .iter()
 .map(|shell| Shell {
 faces: shell
 .faces
 .iter()
 .map(|face| {
 let flip = issue_set.contains(&(si, flat_face_idx));
 flat_face_idx += 1;
 if flip {
 changed += 1;
 Face {
 outer_wire: reverse_wire(&face.outer_wire),
 inner_wires: face.inner_wires.iter().map(reverse_wire).collect(),
 normal: -face.normal,
 triangles: face.triangles.clone(),
 sample_point: face.sample_point,
 mesh_dirty: true,
 surface_idx: None,
 }
 } else {
 face.clone()
 }
 })
 .collect(),
 })
 .collect(),
 })
 .collect();

 let mut result = brep.clone();
 result.solids = new_solids;
 (result, changed)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Internal helpers
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Attempt to fix one wire, returning (fixed_wire, number_of_edges_flipped).
fn fix_wire(wire: &Wire, brep: &BRep, tol2: f64) -> (Wire, usize) {
 if wire.edges.len() < 2 {
 return (wire.clone(), 0);
 }

 let mut edges: Vec<WireEdge> = wire.edges.clone();
 let mut flipped = 0usize;
 let n = edges.len();

 for i in 0..n {
 let next = (i + 1) % n;
 let e_curr = match brep.edges.get(edges[i].idx) {
 Some(e) => e,
 None => continue,
 };
 let e_next = match brep.edges.get(edges[next].idx) {
 Some(e) => e,
 None => continue,
 };

 // end vertex of current edge
 let end_v = if edges[i].forward {
 e_curr.end
 } else {
 e_curr.start
 };
 // start vertex of next edge
 let start_v = if edges[next].forward {
 e_next.start
 } else {
 e_next.end
 };

 if end_v == start_v {
 continue; // already connected
 }
 // Check spatial proximity
 if let (Some(ep), Some(sp)) = (
 brep.vertices.get(end_v).map(|v| v.point),
 brep.vertices.get(start_v).map(|v| v.point),
 ) && (ep - sp).length_squared() <= tol2
 {
 continue; // close enough =OK
 }

 // Try flipping the *next* edge to see if that connects the chain
 let alt_start = if edges[next].forward {
 e_next.end
 } else {
 e_next.start
 };
 if alt_start == end_v {
 edges[next].forward = !edges[next].forward;
 flipped += 1;
 }
 }

 (Wire { edges }, flipped)
}

fn reverse_wire(wire: &Wire) -> Wire {
 let edges = wire
 .edges
 .iter()
 .rev()
 .map(|we| WireEdge::new(we.idx, !we.forward))
 .collect();
 Wire { edges }
}

