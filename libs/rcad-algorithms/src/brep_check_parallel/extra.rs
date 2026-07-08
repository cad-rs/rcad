pub fn validate_solids_parallel(brep: &rcad_kernel::BRep) -> Vec<SolidValidationResult> {
 brep.solids.par_iter()
 .enumerate()
 .map(|(si, _)| validate_single_solid(brep, si))
 .collect()
}

/// Validate a single solid.
fn validate_single_solid(brep: &rcad_kernel::BRep, si: usize) -> SolidValidationResult {
 use std::collections::HashSet;

 let solid = &brep.solids[si];
 let n_edges = brep.edges.len();

 let mut errors = Vec::new();
 let mut warnings = Vec::new();

 // Get shell results
 let shell_results: Vec<ShellValidationResult> = solid.shells
 .iter()
 .enumerate()
 .map(|(shi, _)| validate_single_shell(brep, si, shi))
 .collect();

 // Aggregate counts
 let face_count: usize = solid.shells.iter().map(|s| s.faces.len()).sum();
 let edge_count: usize;
 let vertex_count: usize;

 {
 let mut edges: HashSet<usize> = HashSet::new();
 let mut verts: HashSet<usize> = HashSet::new();

 for shell in &solid.shells {
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if we.idx < n_edges {
 edges.insert(we.idx);
 let edge = &brep.edges[we.idx];
 verts.insert(edge.start);
 verts.insert(edge.end);
 }
 }
 }
 }

 edge_count = edges.len();
 vertex_count = verts.len();
 }

 // Compute Euler characteristic
 let euler_characteristic = vertex_count as i64 - edge_count as i64 + face_count as i64;

 // Check if all shells are closed and manifold
 let is_closed = shell_results.iter().all(|s| s.is_closed);
 let is_manifold = shell_results.iter().all(|s| s.is_manifold);
 let orientation_valid = shell_results.iter().all(|s| s.orientation_consistent);

 // Compute volume (approximate using shell volumes)
 let volume: f64 = solid.shells.iter()
 .map(|shell| compute_shell_volume(shell, brep))
 .sum();

 let has_positive_volume = volume > 0.0;

 // Compute genus
 let genus = if is_closed && is_manifold {
 let g = (2 - euler_characteristic) / 2;
 if (2 - euler_characteristic) % 2 == 0 && g >= 0 { Some(g) } else { None }
 } else {
 None
 };

 // Generate errors
 if !is_closed {
 errors.push("Solid has unclosed shells".to_string());
 }
 if !is_manifold {
 errors.push("Solid has non-manifold topology".to_string());
 }
 if !has_positive_volume {
 warnings.push("Solid has zero or negative volume".to_string());
 }

 let is_valid = errors.is_empty() && shell_results.iter().all(|s| s.is_valid);

 SolidValidationResult {
 solid_idx: si,
 is_valid,
 shell_count: solid.shells.len(),
 face_count,
 edge_count,
 vertex_count,
 euler_characteristic,
 is_closed,
 is_manifold,
 orientation_valid,
 has_positive_volume,
 volume,
 genus,
 shell_results,
 errors,
 warnings,
 }
}

/// Compute the volume of a shell using signed volume method.
fn compute_shell_volume(shell: &rcad_kernel::topology::Shell, brep: &rcad_kernel::BRep) -> f64 {
 

 let n_edges = brep.edges.len();
 let mut volume = 0.0_f64;

 for face in &shell.faces {
 // Get vertices of the outer wire
 let mut verts: Vec<DVec3> = Vec::new();
 for we in &face.outer_wire.edges {
 if we.idx < n_edges {
 let edge = &brep.edges[we.idx];
 let vi = if we.forward { edge.start } else { edge.end };
 if vi < brep.vertices.len() {
 verts.push(brep.vertices[vi].point);
 }
 }
 }

 // Compute signed volume contribution using triangulation
 if verts.len() >= 3 {
 let origin = verts[0];
 for i in 1..verts.len() - 1 {
 let v1 = verts[i] - origin;
 let v2 = verts[i + 1] - origin;
 let signed_vol = v1.cross(v2).dot(face.normal) / 6.0;
 volume += signed_vol;
 }
 }
 }

 volume.abs()
}

// ===========================================================?
// Comprehensive Parallel Check
// ===========================================================?

/// Perform a comprehensive parallel check of a brep.
///
/// This function runs all configured checks in parallel and returns a detailed
/// report including timing information for each phase.
///
/// # Arguments
///
/// * `brep` - The brep to check.
/// * `config` - Configuration for the check.
///
/// # Returns
///
/// A `ParallelCheckReport` containing all results and timing information.
pub fn check_brep_parallel(brep: &rcad_kernel::BRep, config: &ParallelCheckConfig) -> ParallelCheckReport {
 let start_time = Instant::now();
 let mut phase_timings: Vec<CheckPhaseTiming> = Vec::new();
 let mut structural_issues = Vec::new();
 let mut parallel_issues = Vec::new();

 // Configure thread pool
 let threads_used = if config.num_threads > 0 {
 config.num_threads
 } else {
 rayon::current_num_threads()
 };

 // Count totals
 let total_solids = brep.solids.len();
 let total_shells: usize = brep.solids.iter().map(|s| s.shells.len()).sum();
 let total_faces: usize = brep.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();
 let total_edges = brep.edges.len();
 let total_vertices = brep.vertices.len();

 let use_parallel = total_faces >= config.parallel_threshold
 || total_edges >= config.parallel_threshold
 || total_vertices >= config.parallel_threshold;

 // Face checking
 let mut face_results = Vec::new();
 if config.check_faces {
 let phase_start = Instant::now();
 face_results = if use_parallel {
 check_faces_parallel(brep, threads_used)
 } else {
 check_faces_sequential(brep)
 };
 phase_timings.push(CheckPhaseTiming {
 phase: "faces".to_string(),
 duration_ms: phase_start.elapsed().as_millis() as u64,
 items_processed: total_faces,
 });

 // Collect issues from face results
 for fr in &face_results {
 if !fr.is_valid {
 structural_issues.push(CheckIssue::DegenerateFace {
 solid: fr.solid_idx,
 shell: fr.shell_idx,
 face: fr.face_idx,
 });
 }
 }
 }

 // Edge checking
 let mut edge_results = Vec::new();
 if config.check_edges {
 let phase_start = Instant::now();
 edge_results = if use_parallel {
 check_edges_parallel(brep, threads_used)
 } else {
 check_edges_sequential(brep)
 };
 phase_timings.push(CheckPhaseTiming {
 phase: "edges".to_string(),
 duration_ms: phase_start.elapsed().as_millis() as u64,
 items_processed: total_edges,
 });

 // Collect issues from edge results
 for er in &edge_results {
 for issue in &er.issues {
 match issue {
 EdgeCheckIssue::InvalidVertexIndex { vertex_idx } => {
 structural_issues.push(CheckIssue::InvalidVertexIndex {
 edge: er.edge_idx,
 vertex_idx: *vertex_idx,
 });
 }
 EdgeCheckIssue::NonManifold { face_count } => {
 structural_issues.push(CheckIssue::NonManifoldEdge {
 edge_idx: er.edge_idx,
 face_count: *face_count,
 });
 }
 _ => {}
 }
 }
 }
 }

 // Vertex checking
 if config.check_vertices {
 let phase_start = Instant::now();

 // Check for non-finite vertices
 if config.check_finite_vertices {
 for (vidx, v) in brep.vertices.iter().enumerate() {
 if !v.point.is_finite() {
 parallel_issues.push(ParallelCheckIssue::NonFiniteVertex { vertex_idx: vidx });
 }
 }
 }

 // Check for isolated vertices
 if config.check_isolated_vertices {
 let mut referenced = vec![false; brep.vertices.len()];
 for edge in &brep.edges {
 if edge.start < brep.vertices.len() {
 referenced[edge.start] = true;
 }
 if edge.end < brep.vertices.len() {
 referenced[edge.end] = true;
 }
 }
 for (vidx, &is_ref) in referenced.iter().enumerate() {
 if !is_ref {
 parallel_issues.push(ParallelCheckIssue::IsolatedVertex { vertex_idx: vidx });
 }
 }
 }

 // Check for duplicate vertices
 if config.check_duplicate_vertices {
 let duplicates = find_duplicate_vertices_parallel(&brep.vertices, config.tolerance);
 parallel_issues.extend(duplicates);
 }

 phase_timings.push(CheckPhaseTiming {
 phase: "vertices".to_string(),
 duration_ms: phase_start.elapsed().as_millis() as u64,
 items_processed: total_vertices,
 });
 }

 // Shell validation
 let mut shell_results = Vec::new();
 if config.check_shells {
 let phase_start = Instant::now();
 shell_results = validate_shells_parallel(brep);
 phase_timings.push(CheckPhaseTiming {
 phase: "shells".to_string(),
 duration_ms: phase_start.elapsed().as_millis() as u64,
 items_processed: total_shells,
 });
 }

 // Solid validation
 let mut solid_results = Vec::new();
 if config.check_solids {
 let phase_start = Instant::now();
 solid_results = validate_solids_parallel(brep);
 phase_timings.push(CheckPhaseTiming {
 phase: "solids".to_string(),
 duration_ms: phase_start.elapsed().as_millis() as u64,
 items_processed: total_solids,
 });
 }

 let total_duration_ms = start_time.elapsed().as_millis() as u64;

 // Determine overall validity
 let is_valid = structural_issues.is_empty()
 && parallel_issues.is_empty()
 && shell_results.iter().all(|s| s.is_valid)
 && solid_results.iter().all(|s| s.is_valid);

 // Build stats
 let stats = ParallelCheckStats {
 face_count: total_faces,
 edge_count: total_edges,
 vertex_count: total_vertices,
 issue_count: structural_issues.len() + parallel_issues.len(),
 is_valid,
 was_parallel: use_parallel,
 thread_count: threads_used,
 };

 ParallelCheckReport {
 is_valid,
 total_faces,
 total_edges,
 total_vertices,
 total_solids,
 total_shells,
 threads_used,
 was_parallel: use_parallel,
 total_duration_ms,
 phase_timings,
 face_results,
 edge_results,
 shell_results,
 solid_results,
 structural_issues,
 parallel_issues,
 stats,
 }
}

/// Sequential face checking fallback.
fn check_faces_sequential(brep: &rcad_kernel::BRep) -> Vec<FaceCheckResult> {
 let n_edges = brep.edges.len();
 let tolerance = TOLERANCE_MESH_LEGACY;

 let mut results = Vec::new();
 for (si, solid) in brep.solids.iter().enumerate() {
 for (shi, shell) in solid.shells.iter().enumerate() {
 for fi in 0..shell.faces.len() {
 results.push(check_single_face_detailed(brep, si, shi, fi, n_edges, tolerance));
 }
 }
 }
 results
}

/// Sequential edge checking fallback.
fn check_edges_sequential(brep: &rcad_kernel::BRep) -> Vec<EdgeCheckResult> {
 let n_verts = brep.vertices.len();
 let tolerance = TOLERANCE_MESH_LEGACY;

 // Compute edge face counts
 let mut edge_face_counts = vec![0usize; brep.edges.len()];
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 edge_face_counts[we.idx] += 1;
 }
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if we.idx < brep.edges.len() {
 edge_face_counts[we.idx] += 1;
 }
 }
 }
 }
 }
 }

 brep.edges.iter()
 .enumerate()
 .map(|(eidx, edge)| check_single_edge(brep, eidx, edge, n_verts, edge_face_counts[eidx], tolerance))
 .collect()
}

/// Perform parallel check and return detailed statistics.
pub fn check_parallel_with_stats(brep: &rcad_kernel::BRep) -> (CheckResult, ParallelCheckStats) {
 let face_count: usize = brep.solids.iter()
 .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum::<usize>())
 .sum();
 let edge_count = brep.edges.len();
 let vertex_count = brep.vertices.len();

 let options = ParallelCheckOptions::default();
 let result = check_parallel_with_options(brep, &options);

 let stats = ParallelCheckStats {
 face_count,
 edge_count,
 vertex_count,
 issue_count: result.issues.len() + result.parallel_issues.len(),
 is_valid: result.is_valid(),
 was_parallel: result.was_parallel,
 thread_count: result.thread_count,
 };

 (result.to_check_result(), stats)
}

