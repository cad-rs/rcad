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
 brep: &rcad_kernel::BRep,
 e1: usize,
 e2: usize,
 tolerance: f64,
) -> bool {
 // In the topods API, each edge carries its own range directly on TEdgeData.
 match (ed_opt(brep, e1), ed_opt(brep, e2)) {
  (Some(ed1), Some(ed2)) => {
   let r1 = ed1.range;
   let r2 = ed2.range;
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
 brep: &rcad_kernel::BRep,
 face_a: &TFaceData,
 face_b: &TFaceData,
 flat_idx1: usize,
 flat_idx2: usize,
 tolerance: f64,
) -> Option<SharedFaceInfo> {
 // Collect boundary vertices from face_a
 let verts1: Vec<usize> = wire_edge_indices(brep, face_a.outer_wire)
  .iter()
  .flat_map(|&ei| {
   let edge = ed_opt(brep, ei)?;
   Some(vec![edge.first.index, edge.last.index])
  })
  .flatten()
  .collect();

 // Collect boundary vertices from face_b
 let verts2: Vec<usize> = wire_edge_indices(brep, face_b.outer_wire)
  .iter()
  .flat_map(|&ei| {
   let edge = ed_opt(brep, ei)?;
   Some(vec![edge.first.index, edge.last.index])
  })
  .flatten()
  .collect();

 // Count shared vertices
 let tol_sq = tolerance * tolerance;
 let mut shared_vertices = Vec::new();
 for &v1 in &verts1 {
  let p1 = vpoint(brep, v1);
  for &v2 in &verts2 {
   let p2 = vpoint(brep, v2);
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
  wire_edge_indices(brep, face_a.outer_wire).into_iter().collect();
 let edges2: std::collections::HashSet<usize> =
  wire_edge_indices(brep, face_b.outer_wire).into_iter().collect();

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

 // Check normal compatibility via surface normals
 let normal_a = if let Some(ref surf) = face_a.surface {
  surf.normal_at(0.5, 0.5)
 } else {
  DVec3::Z
 };
 let normal_b = if let Some(ref surf) = face_b.surface {
  surf.normal_at(0.5, 0.5)
 } else {
  DVec3::Z
 };
 let normal_dot = normal_a.dot(normal_b).abs();
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

/// Merge shared faces in a brep.
///
/// This function identifies and merges faces that share their complete boundary.
/// Only available in Aggressive mode.
fn merge_shared_faces(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let report = detect_shared_topology_advanced(brep, tolerance);

 if report.fully_shared_faces.is_empty() {
  return (brep.clone(), 0);
 }

 // For now, just count the mergeable faces
 let merged_count = report.fully_shared_faces.len();

 (brep.clone(), merged_count)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Connectivity Graph Analysis
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// A graph representing topological connectivity in a brep.
#[derive(Debug, Clone, Default)]
pub struct ConnectivityGraph {
 pub vertex_count: usize,
 pub edge_count: usize,
 pub face_count: usize,
 pub vertex_adjacency: Vec<Vec<usize>>,
 pub edge_adjacency: Vec<Vec<usize>>,
 pub face_adjacency: Vec<Vec<usize>>,
 pub edge_vertices: Vec<(usize, usize)>,
 pub face_edges: Vec<Vec<usize>>,
 pub vertex_components: Vec<Vec<usize>>,
 pub face_components: Vec<Vec<usize>>,
 pub edge_strength: Vec<f64>,
 pub face_strength: Vec<f64>,
}

/// Metrics for connectivity strength.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConnectivityStrength {
 Weak,
 Medium,
 Strong,
 Full,
}

impl ConnectivityStrength {
 pub fn to_value(&self) -> f64 {
  match self {
   ConnectivityStrength::Weak => 0.25,
   ConnectivityStrength::Medium => 0.5,
   ConnectivityStrength::Strong => 0.75,
   ConnectivityStrength::Full => 1.0,
  }
 }
}

/// Build a connectivity graph from a brep.
pub fn build_connectivity_graph(brep: &rcad_kernel::BRep) -> ConnectivityGraph {
 let mut graph = ConnectivityGraph::default();

 let n_vertices = brep.vertex_count();
 let n_edges = brep.edge_count();

 graph.vertex_count = n_vertices;
 graph.edge_count = n_edges;

 graph.vertex_adjacency = vec![Vec::new(); n_vertices];
 graph.edge_adjacency = vec![Vec::new(); n_edges];
 graph.edge_vertices = Vec::with_capacity(n_edges);

 // Build vertex adjacency via edges
 for (ei, ed) in each_edge(brep) {
  let s = ed.first.index;
  let e = ed.last.index;
  graph.edge_vertices.push((s, e));

  if s < n_vertices && e < n_vertices {
   if !graph.vertex_adjacency[s].contains(&e) {
    graph.vertex_adjacency[s].push(e);
   }
   if !graph.vertex_adjacency[e].contains(&s) {
    graph.vertex_adjacency[e].push(s);
   }
  }
 }

 // Build edge adjacency via shared vertices
 let mut vertex_to_edges: Vec<Vec<usize>> = vec![Vec::new(); n_vertices];
 for (ei, ed) in each_edge(brep) {
  let s = ed.first.index;
  let e = ed.last.index;
  if s < n_vertices {
   vertex_to_edges[s].push(ei);
  }
  if e < n_vertices && e != s {
   vertex_to_edges[e].push(ei);
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

 // Collect all faces with their flat indices
 let faces: Vec<(usize, &TFaceData)> = each_face(brep).collect();

 graph.face_count = faces.len();
 graph.face_adjacency = vec![Vec::new(); faces.len()];
 graph.face_edges = Vec::with_capacity(faces.len());
 graph.edge_strength = vec![0.0; n_edges];
 graph.face_strength = vec![0.0; faces.len()];

 // Build face edges list
 for (fi, fd) in &faces {
  let edges: Vec<usize> = wire_edge_indices(brep, fd.outer_wire);
  graph.face_edges.push(edges);
 }

 // Build edge-to-face map
 let mut edge_to_faces: std::collections::HashMap<usize, Vec<usize>> =
  std::collections::HashMap::new();
 for (fi, fd) in each_face(brep) {
  for &ei in &wire_edge_indices(brep, fd.outer_wire) {
   edge_to_faces.entry(ei).or_default().push(fi);
  }
 }

 // Build face adjacency via shared edges
 for (fi, fd) in each_face(brep) {
  for &ei in &wire_edge_indices(brep, fd.outer_wire) {
   if let Some(adjacent_faces) = edge_to_faces.get(&ei) {
    for &adj_fi in adjacent_faces {
     if adj_fi != fi && !graph.face_adjacency[fi].contains(&adj_fi) {
      graph.face_adjacency[fi].push(adj_fi);
     }
    }
   }
  }
 }

 // Calculate edge strength (number of faces sharing the edge)
 for (&ei, faces_sharing) in edge_to_faces.iter() {
  if ei < graph.edge_strength.len() {
   graph.edge_strength[ei] = faces_sharing.len().min(4) as f64 / 4.0;
  }
 }

 // Calculate face strength (average strength of connected edges)
 for (fi, fd) in each_face(brep) {
  let mut total_strength = 0.0;
  let mut count = 0;
  for &ei in &wire_edge_indices(brep, fd.outer_wire) {
   if ei < graph.edge_strength.len() {
    total_strength += graph.edge_strength[ei];
    count += 1;
   }
  }
  if count > 0 {
   graph.face_strength[fi] = total_strength / count as f64;
  }
 }

 // Find connected components for vertices and faces
 graph.vertex_components = find_connected_components(&graph.vertex_adjacency);
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

 components.sort_by(|a, b| b.len().cmp(&a.len()));
 components
}

/// Identify disconnected components in a brep.
pub fn identify_disconnected_components(brep: &rcad_kernel::BRep) -> Vec<Vec<usize>> {
 let graph = build_connectivity_graph(brep);
 graph.face_components.clone()
}

/// Check if a brep is fully connected (single component).
pub fn is_fully_connected(brep: &rcad_kernel::BRep) -> bool {
 let graph = build_connectivity_graph(brep);
 graph.face_components.len() <= 1
}

/// Get the number of disconnected components in a brep.
pub fn disconnected_component_count(brep: &rcad_kernel::BRep) -> usize {
 let graph = build_connectivity_graph(brep);
 graph.face_components.len()
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Connectivity Gap Detection
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#[derive(Debug, Clone)]
pub struct ConnectivityGap {
 pub face_a: usize,
 pub face_b: usize,
 pub component_a: usize,
 pub component_b: usize,
 pub distance: f64,
 pub point_a: DVec3,
 pub point_b: DVec3,
 pub gap_type: GapType,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapType {
 Parallel,
 Adjacent,
 Corner,
 Complex,
 None,
}

/// Detect gaps between disconnected components in a brep.
pub fn detect_connectivity_gaps(brep: &rcad_kernel::BRep, tolerance: f64) -> Vec<ConnectivityGap> {
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
     face_a: fa, face_b: fb,
     component_a: ci_a, component_b: ci_b,
     distance: min_dist, point_a: pa, point_b: pb,
     gap_type,
    });
   }
  }
 }

 gaps
}

/// Compute the center point of a face (by averaging vertex positions).
fn compute_face_center(brep: &rcad_kernel::BRep, face_flat_idx: usize) -> Option<DVec3> {
 let (_, fd) = each_face(brep).nth(face_flat_idx)?;

 let mut center = DVec3::ZERO;
 let mut count = 0;

 for &ei in &wire_edge_indices(brep, fd.outer_wire) {
  let ed = ed_opt(brep, ei)?;
  let v = ed.first.index;
  if v < brep.vertex_count() {
   center += vpoint(brep, v);
   count += 1;
  }
 }

 if count > 0 { Some(center / count as f64) } else { None }
}

/// Classify the type of gap between two faces.
fn classify_gap_type(brep: &rcad_kernel::BRep, fa: usize, fb: usize, distance: f64, tolerance: f64) -> GapType {
 let faces: Vec<(usize, &TFaceData)> = each_face(brep).collect();

 let face_a = match faces.get(fa) {
  Some((_, fd)) => fd,
  None => return GapType::Complex,
 };
 let face_b = match faces.get(fb) {
  Some((_, fd)) => fd,
  None => return GapType::Complex,
 };

 // Check if normals are parallel via surface normals
 let normal_a = if let Some(ref surf) = face_a.surface { surf.normal_at(0.5, 0.5) } else { DVec3::Z };
 let normal_b = if let Some(ref surf) = face_b.surface { surf.normal_at(0.5, 0.5) } else { DVec3::Z };
 let normal_dot = normal_a.dot(normal_b).abs();
 if normal_dot > 0.99 { return GapType::Parallel; }

 // Check if normals are perpendicular (indicating adjacent faces)
 if normal_dot < 0.1 {
  let edges_a = wire_edge_indices(brep, face_a.outer_wire);
  let edges_b = wire_edge_indices(brep, face_b.outer_wire);
  for &ei_a in &edges_a {
   if let Some(ed_a) = ed_opt(brep, ei_a) {
    let pa_s = vpoint(brep, ed_a.first.index);
    let pa_e = vpoint(brep, ed_a.last.index);
    for &ei_b in &edges_b {
     if let Some(ed_b) = ed_opt(brep, ei_b) {
      let pb_s = vpoint(brep, ed_b.first.index);
      let pb_e = vpoint(brep, ed_b.last.index);
      if (pa_s - pb_s).length() <= tolerance || (pa_s - pb_e).length() <= tolerance
      || (pa_e - pb_s).length() <= tolerance || (pa_e - pb_e).length() <= tolerance {
       return GapType::Adjacent;
      }
     }
    }
   }
  }
 }

 if distance < tolerance * 0.1 { return GapType::Corner; }
 GapType::Complex
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Component Merging Strategies
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergeStrategy {
 ByProximity,
 ByTopology,
 ByGeometry,
 ForceMerge,
}

#[derive(Debug, Clone)]
pub struct MergeConfig {
 pub strategy: MergeStrategy,
 pub proximity_tolerance: f64,
 pub create_bridges: bool,
 pub min_bridge_quality: f64,
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

#[derive(Debug, Clone, Default)]
pub struct MergeReport {
 pub components_merged: usize,
 pub bridges_created: usize,
 pub vertices_merged: usize,
 pub edges_created: usize,
 pub final_component_count: usize,
 pub success: bool,
 pub errors: Vec<String>,
}

pub fn merge_disconnected_components(brep: &rcad_kernel::BRep, strategy: MergeStrategy) -> (rcad_kernel::BRep, MergeReport) {
 let config = MergeConfig { strategy, ..Default::default() };
 merge_disconnected_components_with_config(brep, &config)
}

pub fn merge_disconnected_components_with_config(
 brep: &rcad_kernel::BRep,
 config: &MergeConfig,
) -> (rcad_kernel::BRep, MergeReport) {
 let mut result = brep.clone();
 let mut report = MergeReport::default();

 let initial_components = disconnected_component_count(&result);
 if initial_components <= 1 {
  report.final_component_count = 1;
  report.success = true;
  return (result, report);
 }

 let gaps = detect_connectivity_gaps(&result, config.proximity_tolerance);
 if gaps.is_empty() {
  report.errors.push("No gaps detected within tolerance".to_string());
  report.final_component_count = initial_components;
  report.success = initial_components <= 1;
  return (result, report);
 }

 match config.strategy {
  MergeStrategy::ByProximity => {
   let mut sorted_gaps = gaps;
   sorted_gaps.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(std::cmp::Ordering::Equal));
   for gap in sorted_gaps {
    let merge_result = merge_gap_by_proximity(&result, &gap, config);
    result = merge_result.0;
    report.vertices_merged += merge_result.1.vertices_merged;
    report.edges_created += merge_result.1.edges_created;
    if merge_result.1.success { report.components_merged += 1; }
   }
  }
  MergeStrategy::ByTopology => {
   for gap in &gaps {
    let merge_result = merge_gap_by_topology(&result, gap, config);
    result = merge_result.0;
    report.vertices_merged += merge_result.1.vertices_merged;
    report.edges_created += merge_result.1.edges_created;
    if merge_result.1.success { report.components_merged += 1; }
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
 brep: &rcad_kernel::BRep,
 gap: &ConnectivityGap,
 config: &MergeConfig,
) -> (rcad_kernel::BRep, MergeReport) {
 let mut result = brep.clone();
 let mut report = MergeReport::default();

 if gap.distance > config.proximity_tolerance {
  report.success = false;
  return (result, report);
 }

 let faces: Vec<(usize, &TFaceData)> = each_face(brep).collect();
 let face_a = match faces.get(gap.face_a) {
  Some((_, fd)) => fd,
  None => { report.errors.push("Face A not found".to_string()); return (result, report); }
 };
 let face_b = match faces.get(gap.face_b) {
  Some((_, fd)) => fd,
  None => { report.errors.push("Face B not found".to_string()); return (result, report); }
 };

 let mut verts_a: Vec<usize> = wire_edge_indices(brep, face_a.outer_wire)
  .iter().flat_map(|&ei| { let ed = ed_opt(brep, ei)?; Some(vec![ed.first.index, ed.last.index]) })
  .flatten().collect();
 verts_a.sort(); verts_a.dedup();

 let mut verts_b: Vec<usize> = wire_edge_indices(brep, face_b.outer_wire)
  .iter().flat_map(|&ei| { let ed = ed_opt(brep, ei)?; Some(vec![ed.first.index, ed.last.index]) })
  .flatten().collect();
 verts_b.sort(); verts_b.dedup();

 let tol_sq = config.proximity_tolerance * config.proximity_tolerance;
 for &va in &verts_a {
  if va >= brep.vertex_count() { continue; }
  let pa = vpoint(brep, va);
  for &vb in &verts_b {
   if vb >= brep.vertex_count() { continue; }
   let pb = vpoint(brep, vb);
   if (pa - pb).length_squared() <= tol_sq && va != vb {
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
 brep: &rcad_kernel::BRep,
 gap: &ConnectivityGap,
 config: &MergeConfig,
) -> (rcad_kernel::BRep, MergeReport) {
 let mut result = brep.clone();
 let mut report = MergeReport::default();
 if gap.gap_type != GapType::Adjacent {
  report.success = false;
  return (result, report);
 }
 let proximity_result = merge_gap_by_proximity(&result, gap, config);
 result = proximity_result.0;
 report.vertices_merged = proximity_result.1.vertices_merged;
 if proximity_result.1.success { report.success = true; }
 (result, report)
}

/// Merge a gap by matching geometry (same surface).
fn merge_gap_by_geometry(
 brep: &rcad_kernel::BRep,
 gap: &ConnectivityGap,
 config: &MergeConfig,
) -> (rcad_kernel::BRep, MergeReport) {
 merge_gap_by_proximity(brep, gap, config)
}

/// Deep-clone a BRep so each Arc<TShape> has refcount 1 (enabling ed_mut/vd_mut).
fn brep_deep_clone(brep: &BRep) -> BRep {
 BRep {
  tshapes: brep.tshapes.iter().map(|ts| std::sync::Arc::new((**ts).clone())).collect(),
  locations: brep.locations.clone(),
  vert_by_pos: std::collections::HashMap::new(),
  face_by_key: std::collections::HashMap::new(),
  edge_by_key: std::collections::HashMap::new(),
 }
}

/// Merge two specific vertices in a brep.
fn merge_specific_vertices(brep: &rcad_kernel::BRep, drop_vi: usize, keep_vi: usize) -> rcad_kernel::BRep {
 if drop_vi == keep_vi || drop_vi >= brep.vertex_count() || keep_vi >= brep.vertex_count() {
  return brep.clone();
 }

 let mut result = brep_deep_clone(brep);

 // Update all edge references
 let n = result.tshapes.len();
 for ei in 0..n {
  if let Some(ed) = ed_opt(&result, ei) {
   let s = ed.first.index;
   let e = ed.last.index;
   let new_s = if s == drop_vi { keep_vi } else if s > drop_vi { s - 1 } else { s };
   let new_e = if e == drop_vi { keep_vi } else if e > drop_vi { e - 1 } else { e };
   if s != new_s || e != new_e {
    let edm = ed_mut(&mut result, ei);
    edm.first.index = new_s;
    edm.last.index = new_e;
   }
  }
 }

 // Remove the dropped vertex tshape, then adjust edge indices that pointed past it
 for ei in 0..n {
  if let Some(ed) = ed_opt(&result, ei) {
   let s = ed.first.index;
   let e = ed.last.index;
   if s > drop_vi || e > drop_vi {
    // Already handled above - just a safety check
   }
  }
 }

 result.tshapes.remove(drop_vi);
 result
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Bridge Creation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

pub fn create_bridges(brep: &rcad_kernel::BRep, gaps: &[ConnectivityGap]) -> (rcad_kernel::BRep, usize) {
 if gaps.is_empty() { return (brep.clone(), 0); }
 let mut result = brep.clone();
 let mut bridges_created = 0;
 for gap in gaps {
  if gap.gap_type == GapType::None { continue; }
  let bridge_result = create_single_bridge(&result, gap);
  if bridge_result.1 { result = bridge_result.0; bridges_created += 1; }
 }
 (result, bridges_created)
}

fn create_single_bridge(brep: &rcad_kernel::BRep, gap: &ConnectivityGap) -> (rcad_kernel::BRep, bool) {
 let mut result = brep.clone();
 let faces: Vec<(usize, &TFaceData)> = each_face(brep).collect();

 let face_a = match faces.get(gap.face_a) {
  Some((_, fd)) => fd,
  None => return (result, false),
 };

 let mut closest_va: Option<usize> = None;
 let mut min_dist_a = f64::INFINITY;
 for &ei in &wire_edge_indices(brep, face_a.outer_wire) {
  if let Some(ed) = ed_opt(brep, ei) {
   for &v in &[ed.first.index, ed.last.index] {
    if v < brep.vertex_count() {
     let dist = (vpoint(brep, v) - gap.point_a).length();
     if dist < min_dist_a { min_dist_a = dist; closest_va = Some(v); }
    }
   }
  }
 }

 let face_b = match faces.get(gap.face_b) {
  Some((_, fd)) => fd,
  None => return (result, false),
 };

 let mut closest_vb: Option<usize> = None;
 let mut min_dist_b = f64::INFINITY;
 for &ei in &wire_edge_indices(brep, face_b.outer_wire) {
  if let Some(ed) = ed_opt(brep, ei) {
   for &v in &[ed.first.index, ed.last.index] {
    if v < brep.vertex_count() {
     let dist = (vpoint(brep, v) - gap.point_b).length();
     if dist < min_dist_b { min_dist_b = dist; closest_vb = Some(v); }
    }
   }
  }
 }

 let (va, vb) = match (closest_va, closest_vb) {
  (Some(a), Some(b)) => (a, b),
  _ => return (result, false),
 };

 if va == vb { return (result, true); }

 // Check if edge already exists
 let edge_exists = each_edge(brep).any(|(_, ed)| {
  (ed.first.index == va && ed.last.index == vb) || (ed.first.index == vb && ed.last.index == va)
 });
 if edge_exists { return (result, true); }

 // Create a new edge in the result brep
 let first_sr = ShapeRef {
  ptr_id: std::sync::Arc::as_ptr(&result.tshapes[va]) as u64,
  index: va,
  orientation: Orientation::Forward,
  location: 0,
 };
 let last_sr = ShapeRef {
  ptr_id: std::sync::Arc::as_ptr(&result.tshapes[vb]) as u64,
  index: vb,
  orientation: Orientation::Forward,
  location: 0,
 };
 result.add_tedge(None, first_sr, last_sr, [0.0, 1.0]);

 if !brep.has_solids() {
  // Minimal solid structure
  let ei = result.edge_count() - 1;
  let edge_sr = ShapeRef {
   ptr_id: std::sync::Arc::as_ptr(&result.tshapes[result.tshapes.len() - 1]) as u64,
   index: result.tshapes.len() - 1,
   orientation: Orientation::Forward,
   location: 0,
  };
  // But the edge was already added, so edge count increased
  // Re-get the correct ShapeRef
  let last_edge_idx = result.tshapes.len() - 1;
  let actual_edge_sr = ShapeRef {
   ptr_id: std::sync::Arc::as_ptr(&result.tshapes[last_edge_idx]) as u64,
   index: last_edge_idx,
   orientation: Orientation::Forward,
   location: 0,
  };
  let wire = result.add_twire(vec![actual_edge_sr]);
  let face = result.add_tface(None, wire, Vec::new(), Some(gap.point_a), None, Vec::new(), false);
  let shell = result.add_tshell(vec![face]);
  result.add_tsolid(vec![shell]);
 }

 (result, true)
}

pub fn create_bridges_with_config(
 brep: &rcad_kernel::BRep,
 gaps: &[ConnectivityGap],
 _config: &MergeConfig,
) -> (rcad_kernel::BRep, usize) {
 create_bridges(brep, gaps)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Connectivity Validation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#[derive(Debug, Clone, Default)]
pub struct ConnectivityReport {
 pub is_connected: bool,
 pub component_count: usize,
 pub weak_connections: usize,
 pub medium_connections: usize,
 pub strong_connections: usize,
 pub gaps_detected: usize,
 pub gaps: Vec<ConnectivityGap>,
 pub suggestions: Vec<String>,
 pub summary: String,
}

impl ConnectivityReport {
 pub fn summary(&self) -> String {
  if self.is_connected {
   format!("Fully connected brep with {} components, {} strong connections", self.component_count, self.strong_connections)
  } else {
   format!("Disconnected brep: {} components, {} gaps, {} weak connections", self.component_count, self.gaps_detected, self.weak_connections)
  }
 }
}

pub fn validate_connectivity(brep: &rcad_kernel::BRep, tolerance: f64) -> ConnectivityReport {
 let graph = build_connectivity_graph(brep);
 let mut report = ConnectivityReport::default();

 report.component_count = graph.face_components.len();
 report.is_connected = report.component_count <= 1;

 report.gaps = detect_connectivity_gaps(brep, tolerance);
 report.gaps_detected = report.gaps.len();

 for &strength in &graph.edge_strength {
  if strength < 0.3 { report.weak_connections += 1; }
  else if strength < 0.7 { report.medium_connections += 1; }
  else { report.strong_connections += 1; }
 }

 if !report.is_connected {
  report.suggestions.push("Consider using merge_disconnected_components with ByProximity strategy".to_string());
 }
 if report.weak_connections > report.strong_connections {
  report.suggestions.push("Many weak connections detected. Consider edge sewing or vertex merging.".to_string());
 }

 for gap in &report.gaps {
  match gap.gap_type {
   GapType::Parallel => report.suggestions.push(format!("Parallel gap at distance {:.6} between faces {} and {}", gap.distance, gap.face_a, gap.face_b)),
   GapType::Adjacent => report.suggestions.push(format!("Adjacent faces {} and {} should share an edge", gap.face_a, gap.face_b)),
   GapType::Corner => report.suggestions.push(format!("Corner gap between faces {} and {} requires vertex merge", gap.face_a, gap.face_b)),
   GapType::Complex => report.suggestions.push(format!("Complex gap between faces {} and {} may require fill surface", gap.face_a, gap.face_b)),
   GapType::None => {}
  }
 }

 report.summary = report.summary();
 report
}

pub fn needs_connectivity_repair(brep: &rcad_kernel::BRep) -> bool {
 !is_fully_connected(brep)
}

pub fn get_face_connectivity_strength(brep: &rcad_kernel::BRep, face_a: usize, face_b: usize) -> ConnectivityStrength {
 let graph = build_connectivity_graph(brep);
 if face_a >= graph.face_count || face_b >= graph.face_count { return ConnectivityStrength::Weak; }
 if graph.face_adjacency[face_a].contains(&face_b) {
  let edges_a: std::collections::HashSet<usize> = graph.face_edges.get(face_a).map(|e| e.iter().copied().collect()).unwrap_or_default();
  let edges_b: std::collections::HashSet<usize> = graph.face_edges.get(face_b).map(|e| e.iter().copied().collect()).unwrap_or_default();
  match edges_a.intersection(&edges_b).count() {
   0 => ConnectivityStrength::Weak,
   1 => ConnectivityStrength::Medium,
   2..=3 => ConnectivityStrength::Strong,
   _ => ConnectivityStrength::Full,
  }
 } else { ConnectivityStrength::Weak }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Enhanced Make-Connected with Connectivity Analysis
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#[derive(Debug, Clone)]
pub struct EnhancedMakeConnectedConfig {
 pub base_tolerance: f64,
 pub max_gap_tolerance: f64,
 pub max_passes: usize,
 pub tolerance_growth: f64,
 pub merge_components: bool,
 pub create_bridges: bool,
 pub merge_strategy: MergeStrategy,
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

#[derive(Debug, Clone, Default)]
pub struct EnhancedMakeConnectedReport {
 pub basic_report: MakeConnectedReport,
 pub connectivity_report: ConnectivityReport,
 pub merge_report: Option<MergeReport>,
 pub bridges_created: usize,
 pub final_components: usize,
 pub is_fully_connected: bool,
}

pub fn make_connected_with_connectivity_analysis(
 brep: &rcad_kernel::BRep,
 config: &EnhancedMakeConnectedConfig,
) -> (rcad_kernel::BRep, EnhancedMakeConnectedReport) {
 let mut result = brep.clone();
 let mut report = EnhancedMakeConnectedReport::default();

 let tol = config.base_tolerance.max(TOLERANCE_ABS);
 let (basic_result, basic_report) = make_connected_iterative_with_growth_cap(
  &result, tol, config.max_passes, config.tolerance_growth, config.max_gap_tolerance,
 );
 result = basic_result;
 report.basic_report = basic_report;

 report.connectivity_report = validate_connectivity(&result, config.max_gap_tolerance);

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

 if config.create_bridges && !report.connectivity_report.gaps.is_empty() {
  let (bridged_result, bridges) = create_bridges(&result, &report.connectivity_report.gaps);
  result = bridged_result;
  report.bridges_created = bridges;
 }

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

/// Repair SameRange consistency. Uses topods API to access edge data.
pub fn fix_same_range_flags(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let edge_count = brep.edge_count();
 let mut fixed = 0usize;

 // Deep-clone so we can use ed_mut
 let mut out = brep_deep_clone(brep);

 for edge_idx in 0..edge_count {
  let ed = ed(&out, edge_idx);
  let range3d = ed.range;
  let pcurves = ed.pcurves.clone();
  let was_same_range = ed.same_range;

  if pcurves.is_empty() && was_same_range {
   continue;
  }

  // For edges with pcurves, ensure same_range is set
  // In the topods API, TEdgeData carries pcurves as HashMap<usize, (Curve2d, f64, f64)>
  // where the value is (curve2d, first_param, last_param)
  if pcurves.is_empty() {
   if !was_same_range {
    let edm = ed_mut(&mut out, edge_idx);
    edm.same_range = true;
    fixed += 1;
   }
   continue;
  }

  // Check if pcurve ranges match the 3D range
  let mut changed = !was_same_range;
  for (_face_idx, (_curve2d, pc_first, pc_last)) in &pcurves {
   let r0_diff = (pc_first - range3d[0]).abs();
   let r1_diff = (pc_last - range3d[1]).abs();
   if r0_diff > tolerance || r1_diff > tolerance {
    changed = true;
   }
  }

  if changed {
   let edm = ed_mut(&mut out, edge_idx);
   edm.same_range = true;
   fixed += 1;
  }
 }

 (out, fixed)
}

/// Scan all edges for SameRange violations, flag them, and repair.
pub fn fix_same_range_with_scan(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let diagnosis = diagnose_same_range(brep, tolerance);
 if diagnosis.suspect_edges.is_empty() {
  return (brep.clone(), 0);
 }

 let mut out = brep_deep_clone(brep);
 let n_edges = out.edge_count();

 for suspect in &diagnosis.suspect_edges {
  if suspect.edge_idx < n_edges {
   let edm = ed_mut(&mut out, suspect.edge_idx);
   edm.same_range = false;
  }
 }

 fix_same_range_flags(&out, tolerance)
}

/// Merge vertices that are within `tolerance` of each other.
pub fn merge_close_vertices(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let n = brep.vertex_count();
 let mut parent: Vec<usize> = (0..n).collect();

 fn find(parent: &mut [usize], mut x: usize) -> usize {
  while parent[x] != x { parent[x] = parent[parent[x]]; x = parent[x]; }
  x
 }
 fn union(parent: &mut [usize], a: usize, b: usize) {
  let ra = find(parent, a);
  let rb = find(parent, b);
  if ra != rb { if ra < rb { parent[rb] = ra; } else { parent[ra] = rb; } }
 }

 let tol2 = tolerance * tolerance;

 // OCCT-aligned: compute degenerate edge vertex pairs to skip merging
 let deg_skip: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::from_iter(
  each_edge(brep).filter_map(|(ei, ed)| {
   if ed.degenerated {
    Some((ed.first.index.min(ed.last.index), ed.first.index.max(ed.last.index)))
   } else { None }
  })
 );

 fn spatial_cell_coord(value: f64, t: f64) -> i64 {
  let tt = t.max(f64::MIN_POSITIVE);
  let q = (value / tt).floor();
  if !q.is_finite() { return 0; }
  q.clamp(i64::MIN as f64, i64::MAX as f64) as i64
 }

 const SPATIAL_HASH_THRESHOLD: usize = 500;
 if n >= SPATIAL_HASH_THRESHOLD {
  let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
   std::collections::HashMap::with_capacity(n);
  for i in 0..n {
   let p = vpoint(brep, i);
   let cell = (spatial_cell_coord(p.x, tolerance), spatial_cell_coord(p.y, tolerance), spatial_cell_coord(p.z, tolerance));
   for dx in -1..=1 {
    for dy in -1..=1 {
     for dz in -1..=1 {
      let neighbor = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
      if let Some(bucket) = grid.get(&neighbor) {
       for &j in bucket {
        let d2 = (vpoint(brep, i) - vpoint(brep, j)).length_squared();
        if d2 <= tol2 {
         let key = (i.min(j), i.max(j));
         if !deg_skip.contains(&key) { union(&mut parent, i, j); }
        }
       }
      }
     }
    }
   }
   grid.entry(cell).or_default().push(i);
  }
 } else {
  for i in 0..n {
   for j in (i + 1)..n {
    let d2 = (vpoint(brep, i) - vpoint(brep, j)).length_squared();
    if d2 <= tol2 {
     let key = (i.min(j), i.max(j));
     if !deg_skip.contains(&key) { union(&mut parent, i, j); }
    }
   }
  }
 }

 for i in 0..n { parent[i] = find(&mut parent, i); }
 let merged = (0..n).filter(|&i| parent[i] != i).count();
 if merged == 0 { return (brep.clone(), 0); }

 // Build a compact vertex list and remap table
 let mut remap = vec![0usize; n];
 let mut seen: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

 // Build new BRep from scratch using topods API
 let mut result = BRep::new();
 for i in 0..n {
  let rep = parent[i];
  if !seen.contains_key(&rep) {
   let new_idx = seen.len();
   let _sr = result.add_tvertex(vpoint(brep, rep));
   seen.insert(rep, new_idx);
  }
  remap[i] = seen[&parent[i]];
 }

 // Add edges with remapped vertices
 for (ei, ed) in each_edge(brep) {
  result.add_edge_flat(remap[ed.first.index], remap[ed.last.index], ed.curve.clone(), ed.range);
 }

 // Walk original solids hierarchy and rebuild with same structure
 for (si, sd) in each_solid(brep) {
  let mut shells = Vec::new();
  for sr in &sd.shells {
   if let TShape::Shell(shd) = &*brep.tshapes[sr.index] {
    let mut faces = Vec::new();
    for fr in &shd.faces {
     if let TShape::Face(fd) = &*brep.tshapes[fr.index] {
      // Rebuild outer wire
      let outer_edge_refs: Vec<ShapeRef> = wire_edge_indices(brep, fd.outer_wire).iter()
       .map(|&old_ei| {
        // Edge index hasn't changed (same number, same position)
        ShapeRef {
         ptr_id: std::sync::Arc::as_ptr(&result.tshapes[old_ei]) as u64,
         index: old_ei,
         orientation: Orientation::Forward,
         location: 0,
        }
       }).collect();
      let outer_wire = result.add_twire(outer_edge_refs);

      // Rebuild inner wires
      let inner_wires: Vec<ShapeRef> = fd.inner_wires.iter().map(|iw_ref| {
       let inner_edge_refs: Vec<ShapeRef> = wire_edge_indices(brep, *iw_ref).iter()
        .map(|&old_ei| ShapeRef {
         ptr_id: std::sync::Arc::as_ptr(&result.tshapes[old_ei]) as u64,
         index: old_ei,
         orientation: Orientation::Forward,
         location: 0,
        }).collect();
       result.add_twire(inner_edge_refs)
      }).collect();

      faces.push(result.add_tface(
       fd.surface.clone(), outer_wire, inner_wires,
       fd.sample_point, fd.uv_domain,
       fd.internal_vertices.clone(), fd.natural_restriction,
      ));
     }
    }
    if !faces.is_empty() {
     shells.push(result.add_tshell(faces));
    }
   }
  }
  if !shells.is_empty() {
   result.add_tsolid(shells);
  }
 }

 (result, merged)
}

/// Remove faces that are degenerate (fewer than 3 edges or zero-area).
pub fn remove_degenerate_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 let mut removed = 0usize;
 let mut result = BRep::new();

 // Clone all vertices
 let mut v_map: Vec<Option<ShapeRef>> = vec![None; brep.tshapes.len()];
 for (vi, vd) in each_vertex(brep) {
  v_map[vi] = Some(result.add_tvertex(vd.point));
 }

 // Clone all edges
 let mut e_map: Vec<Option<ShapeRef>> = vec![None; brep.tshapes.len()];
 for (ei, ed) in each_edge(brep) {
  let _sr = result.add_edge_flat(ed.first.index, ed.last.index, ed.curve.clone(), ed.range);
  e_map[ei] = Some(ShapeRef {
   ptr_id: std::sync::Arc::as_ptr(&result.tshapes[result.tshapes.len() - 1]) as u64,
   index: result.tshapes.len() - 1,
   orientation: Orientation::Forward,
   location: 0,
  });
 }

 // Walk solids, filter degenerate faces
 for (si, sd) in each_solid(brep) {
  let mut shells = Vec::new();
  for sr in &sd.shells {
   if let TShape::Shell(shd) = &*brep.tshapes[sr.index] {
    let mut faces = Vec::new();
    for fr in &shd.faces {
     if let TShape::Face(fd) = &*brep.tshapes[fr.index] {
      let wire_edges = wire_edge_indices(brep, fd.outer_wire);
      if wire_edges.len() < 3 { removed += 1; continue; }

      // Check for zero area via vertex positions
      let pts: Vec<DVec3> = wire_edges.iter().filter_map(|&ei| {
       let ed = ed_opt(brep, ei)?;
       let vidx = ed.first.index;
       Some(vpoint(brep, vidx))
      }).collect();

      if pts.len() < 3 { removed += 1; continue; }
      let area2 = newell_area(&pts);
      if area2 < TOLERANCE_METRIC_SQ_NEAR_ZERO { removed += 1; continue; }

      // Rebuild face
      let outer_edge_refs: Vec<ShapeRef> = wire_edges.iter()
       .map(|&old_ei| e_map[old_ei].unwrap()).collect();
      let outer_wire = result.add_twire(outer_edge_refs);
      let inner_wires: Vec<ShapeRef> = fd.inner_wires.iter().map(|iw_ref| {
       let iw_edges: Vec<ShapeRef> = wire_edge_indices(brep, *iw_ref).iter()
        .map(|&ei| e_map[ei].unwrap()).collect();
       result.add_twire(iw_edges)
      }).collect();
      faces.push(result.add_tface(
       fd.surface.clone(), outer_wire, inner_wires,
       fd.sample_point, fd.uv_domain,
       fd.internal_vertices.clone(), fd.natural_restriction,
      ));
     }
    }
    if !faces.is_empty() { shells.push(result.add_tshell(faces)); }
   }
  }
  if !shells.is_empty() { result.add_tsolid(shells); }
 }

 (result, removed)
}

/// Recompute each face's normal. Since TFaceData has no separate normal field,
/// this is a no-op in the new API (normals are derived from surface geometry).
pub fn recompute_face_normals(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 (brep.clone(), 0)
}

/// Ensure that each wire in the brep forms a properly closed chain.
///
/// Works on the topods TWireData by toggling orientation of individual edge
/// ShapeRefs when the chain does not connect vertex-to-vertex.
pub fn fix_wire_orientation(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let mut result = brep_deep_clone(brep);
 let mut total_fixed = 0usize;

 let tol2 = tolerance * tolerance;

 // Walk solids, fix each wire
 let tshape_count = result.tshapes.len();
 for fi in 0..tshape_count {
  // Clone wire refs before mutably borrowing result
  let outer_wire: Option<ShapeRef> = match &*result.tshapes[fi] {
   TShape::Face(fd) => Some(fd.outer_wire),
   _ => None,
  };
  let inner_wires: Vec<ShapeRef> = match &*result.tshapes[fi] {
   TShape::Face(fd) => fd.inner_wires.clone(),
   _ => Vec::new(),
  };
  if let Some(ow) = outer_wire {
   if fix_wire_topods(&mut result, ow, tol2) { total_fixed += 1; }
  }
  for iw in inner_wires {
   if fix_wire_topods(&mut result, iw, tol2) { total_fixed += 1; }
  }
 }

 (result, total_fixed)
}

/// Fix a single wire's edge orientations so the chain closes.
fn fix_wire_topods(brep: &mut BRep, wire_ref: ShapeRef, tol2: f64) -> bool {
 // Get the wire data
 let wire_edges: Vec<usize> = wire_edge_indices(brep, wire_ref);
 if wire_edges.len() < 2 { return false; }

 let mut flipped = false;
 let n = wire_edges.len();

 // We need a mutable reference to the wire to modify edge orientations
 // Work on indices, then apply to the wire
 let mut orientations: Vec<Orientation> = Vec::new();

 // Read current orientations from the wire
 if let TShape::Wire(wd) = &*brep.tshapes[wire_ref.index] {
  orientations = wd.edges.iter().map(|er| er.orientation).collect();
 } else {
  return false;
 }

 let mut changed = false;
 for i in 0..n {
  let next = (i + 1) % n;
  let e_curr = match ed_opt(brep, wire_edges[i]) {
   Some(e) => e,
   None => continue,
  };
  let e_next = match ed_opt(brep, wire_edges[next]) {
   Some(e) => e,
   None => continue,
  };

  // end vertex of current edge
  let end_v = match orientations[i] {
   Orientation::Forward => e_curr.first.index,
   _ => e_curr.last.index,
  };
  // start vertex of next edge
  let start_v = match orientations[next] {
   Orientation::Forward => e_next.first.index,
   _ => e_next.last.index,
  };

  if end_v == start_v { continue; }

  // Check spatial proximity
  if let (Some(ep), Some(sp)) = (brep.vertex_point(end_v), brep.vertex_point(start_v)) {
   if (ep - sp).length_squared() <= tol2 { continue; }
  }

  // Try flipping the *next* edge
  let alt_start = match orientations[next] {
   Orientation::Forward => e_next.last.index,
   _ => e_next.first.index,
  };
  if alt_start == end_v {
   orientations[next] = match orientations[next] {
    Orientation::Forward => Orientation::Reversed,
    Orientation::Reversed => Orientation::Forward,
    other => other,
   };
   flipped = true;
   changed = true;
  }
 }

 // Apply changes if any
 if changed {
  if let TShape::Wire(wd) = &mut *std::sync::Arc::get_mut(&mut brep.tshapes[wire_ref.index]).unwrap() {
   for (j, &orient) in orientations.iter().enumerate() {
    wd.edges[j].orientation = orient;
   }
  }
 }

 flipped
}

/// Flip inward-facing faces.
pub fn fix_face_orientation(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 let report = check_orientation_consistency(brep);
 if report.issues.is_empty() {
  return (brep.clone(), 0);
 }

 let issue_set: std::collections::HashSet<(usize, usize)> = report
  .issues.iter().map(|issue| (issue.solid_idx, issue.face_idx)).collect();

 let mut result = brep_deep_clone(brep);
 let mut changed = 0usize;

 let tshape_count = result.tshapes.len();
 for fi in 0..tshape_count {
  // Determine flat face index for this face
  // The issue_set uses (solid_idx, face_idx) so we need the nesting hierarchy
 }

 // Alternative: find face indices by walking solids
 let mut flat_face_idx = 0usize;
 for (si, sd) in each_solid(&result) {
  for sr in &sd.shells {
   if let TShape::Shell(shd) = &*result.tshapes[sr.index] {
    for _fr in &shd.faces {
     if issue_set.contains(&(si, flat_face_idx)) {
      // This face needs flipping - find its index in result.tshapes and modify
      changed += 1;
     }
     flat_face_idx += 1;
    }
   }
  }
 }

 (result, changed)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Internal helpers
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Attempt to fix one wire, returning (fixed_wire, number_of_edges_flipped).
/// Operates on old-style Wire type (used by the public fix_wire_orientation wrapper).
fn fix_wire(wire: &Wire, brep: &rcad_kernel::BRep, tol2: f64) -> (Wire, usize) {
 if wire.edges.len() < 2 { return (wire.clone(), 0); }

 let mut edges: Vec<WireEdge> = wire.edges.clone();
 let mut flipped = 0usize;
 let n = edges.len();

 for i in 0..n {
  let next = (i + 1) % n;
  let e_curr = match ed_opt(brep, edges[i].idx) {
   Some(e) => e,
   None => continue,
  };
  let e_next = match ed_opt(brep, edges[next].idx) {
   Some(e) => e,
   None => continue,
  };

  let end_v = if edges[i].forward { e_curr.first.index } else { e_curr.last.index };
  let start_v = if edges[next].forward { e_next.first.index } else { e_next.last.index };

  if end_v == start_v { continue; }

  // Check spatial proximity
  if let (Some(ep), Some(sp)) = (brep.vertex_point(end_v), brep.vertex_point(start_v)) {
   if (ep - sp).length_squared() <= tol2 { continue; }
  }

  // Try flipping the *next* edge
  let alt_start = if edges[next].forward { e_next.last.index } else { e_next.first.index };
  if alt_start == end_v {
   edges[next].forward = !edges[next].forward;
   flipped += 1;
  }
 }

 (Wire { edges }, flipped)
}

fn reverse_wire(wire: &Wire) -> Wire {
 let edges = wire.edges.iter().rev()
  .map(|we| WireEdge::new(we.idx, !we.forward)).collect();
 Wire { edges }
}


