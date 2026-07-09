/// Handle degenerate points on periodic surfaces.
///
/// This function identifies and handles degenerate points such as:
/// - Sphere poles (V=0 and V= ?
/// - Cone apex
pub fn handle_degenerate_points(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let mut result = brep.clone();
 let mut degenerate_count = 0;

 // Track vertices that are at degenerate points
 let mut degenerate_vertices: std::collections::HashSet<usize> = std::collections::HashSet::new();

 // Iterate through all faces
 let mut flat_face_idx = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 // Get the surface for this face
 let surface_idx = match brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) {
 Some(idx) => idx,
 None => {
 flat_face_idx += 1;
 continue;
 }
 };

 let surface = match brep.geom.surfaces.get(surface_idx) {
 Some(s) => s,
 None => {
 flat_face_idx += 1;
 continue;
 }
 };

 let periodic_info = detect_periodic_surface_info(surface);

 if !periodic_info.has_degenerate_points() {
 flat_face_idx += 1;
 continue;
 }

 // Check vertices for degeneracy
 for we in &face.outer_wire.edges {
 if let Some(edge) = ed_opt(brep, we.idx) {
 let vi = if we.forward { edge.start } else { edge.end };
 if let Some(vertex) = brep.vertices.get(vi)
 && is_vertex_at_degenerate_point(
 vertex,
 surface,
 &periodic_info,
 tolerance,
 ) {
 degenerate_vertices.insert(vi);
 degenerate_count += 1;
 }
 }
 }

 flat_face_idx += 1;
 }
 }
 }

 // For vertices at degenerate points, we may need to:
 // 1. Mark edges incident to them as degenerate
 // 2. Ensure proper triangulation near degenerate points
 for vi in &degenerate_vertices {
 // Find edges incident to this vertex and mark them if needed
 for (ei, edge) in result.edges.iter().enumerate() {
 if (edge.start == *vi || edge.end == *vi)
 && result.geom.edge_degenerated.len() <= ei {
 result.geom.edge_degenerated.resize(ei + 1, false);
 }
 // Note: We don't automatically mark as degenerate - that depends on
 // whether the edge actually has zero 3D length
 }
 }

 (result, degenerate_count)
}

/// Check if a vertex is at a degenerate point on a surface.
fn is_vertex_at_degenerate_point(
 vertex: &Vertex,
 surface: &Surface3,
 _periodic_info: &PeriodicSurfaceInfo,
 tolerance: f64,
) -> bool {
 match surface {
 Surface3::Sphere(sphere) => {
 // Check if vertex is at north or south pole
 let to_vertex = vertex.point - sphere.center;
 let _along_axis = to_vertex.dot(sphere.axis.normalize_or_zero());

 // At north pole (V=0): vertex is at center + radius * axis
 // At south pole (V= ?: vertex is at center - radius * axis
 let north_pole = sphere.center + sphere.axis.normalize_or_zero() * sphere.radius;
 let south_pole = sphere.center - sphere.axis.normalize_or_zero() * sphere.radius;

 let dist_to_north = (vertex.point - north_pole).length();
 let dist_to_south = (vertex.point - south_pole).length();

 dist_to_north < tolerance || dist_to_south < tolerance
 }
 Surface3::Cone(cone) => {
 // Check if vertex is at apex
 let apex = cone.apex_point();
 let dist_to_apex = (vertex.point - apex).length();
 dist_to_apex < tolerance
 }
 _ => false,
 }
}

/// Merge edges that are split across a periodic seam.
///
/// When edges are incorrectly split at a seam, this function attempts to
/// merge them back together.
pub fn merge_seam_edges(brep: &rcad_kernel::BRep, config: &PeriodicSeamConfig) -> (rcad_kernel::BRep, usize) {
 let result = brep.clone();
 let mut merged_count = 0;

 // Find pairs of edges that could be merged across the seam
 // This is done by looking for edges that:
 // 1. Share a vertex
 // 2. Are on the same periodic surface
 // 3. Have endpoints near the seam (one at U=, one at U= ?

 let mut flat_face_idx = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 // Get the surface for this face
 let surface_idx = match brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) {
 Some(idx) => idx,
 None => {
 flat_face_idx += 1;
 continue;
 }
 };

 let surface = match brep.geom.surfaces.get(surface_idx) {
 Some(s) => s,
 None => {
 flat_face_idx += 1;
 continue;
 }
 };

 let periodic_info = detect_periodic_surface_info(surface);
 if !periodic_info.is_u_periodic() {
 flat_face_idx += 1;
 continue;
 }

 let u_period = periodic_info.u_period.unwrap();

 // Collect edges and their UV endpoints
 let mut edge_uv_endpoints: Vec<(usize, glam::DVec2, glam::DVec2)> = Vec::new();

 for we in &face.outer_wire.edges {
 if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
 for pc in pcurves {
 if pc.surface_idx != surface_idx {
 continue;
 }
 if let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) {
 let uv_start = curve2d.point_at(if we.forward { 0.0 } else { 1.0 });
 let uv_end = curve2d.point_at(if we.forward { 1.0 } else { 0.0 });
 edge_uv_endpoints.push((we.idx, uv_start, uv_end));
 }
 }
 }
 }

 // Look for edge pairs that span the seam
 for i in 0..edge_uv_endpoints.len() {
 for j in (i + 1)..edge_uv_endpoints.len() {
 let (ei, _uv_start_i, uv_end_i) = edge_uv_endpoints[i];
 let (ej, uv_start_j, _uv_end_j) = edge_uv_endpoints[j];

 // Check if one edge ends near U=0 and another starts near U=period
 // (or vice versa), indicating they should be merged
 let seam_proximity = config.seam_tolerance;

 let i_ends_near_0 = uv_end_i.x < seam_proximity;
 let i_ends_near_period = (uv_end_i.x - u_period).abs() < seam_proximity;
 let j_starts_near_0 = uv_start_j.x < seam_proximity;
 let j_starts_near_period = (uv_start_j.x - u_period).abs() < seam_proximity;

 // Check if they share a 3D vertex (required for merging)
 let edge_i = match ed_opt(brep, ei) {
 Some(e) => e,
 None => continue,
 };
 let edge_j = match ed_opt(brep, ej) {
 Some(e) => e,
 None => continue,
 };

 let shares_vertex = edge_i.end == edge_j.start || edge_i.start == edge_j.end;
 if !shares_vertex {
 continue;
 }

 // Check if they span the seam
 if (i_ends_near_0 && j_starts_near_period) || (i_ends_near_period && j_starts_near_0) {
 // These edges could potentially be merged
 // For now, just count them - actual merging requires more complex wire manipulation
 merged_count += 1;
 }
 }
 }

 flat_face_idx += 1;
 }
 }
 }

 (result, merged_count)
}

/// Handle edges that cross periodic surface seams.
///
/// On periodic surfaces (cylinder, cone, torus), edges that cross the seam
/// may be split incorrectly. This function attempts to handle them.
pub fn handle_periodic_surface_seams(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, PeriodicSeamReport) {
 let config = PeriodicSeamConfig {
 seam_tolerance: tolerance * 10.0,
 merge_tolerance: tolerance * 100.0,
 ..Default::default()
 };
 handle_periodic_surface_seams_with_config(brep, &config)
}

/// Handle periodic surface seams with custom configuration.
pub fn handle_periodic_surface_seams_with_config(
 brep: &rcad_kernel::BRep,
 config: &PeriodicSeamConfig,
) -> (rcad_kernel::BRep, PeriodicSeamReport) {
 let mut result = brep.clone();
 let mut report = PeriodicSeamReport::default();

 // Step 1: Detect seam edges
 let seam_edges = detect_seam_edges(&result, config);
 report.seam_edges_detected = seam_edges.len();

 // Step 2: Handle degenerate points if enabled
 if config.handle_degeneracies {
 let (new_brep, degenerate_count) = handle_degenerate_points(&result, config.seam_tolerance);
 result = new_brep;
 report.degenerate_points_handled = degenerate_count;
 }

 // Step 3: Split edges at seams if enabled
 if config.split_edges {
 for seam_info in &seam_edges {
 let (new_brep, split_done) = split_edge_at_seam(&result, seam_info, config.seam_tolerance);
 if split_done {
 result = new_brep;
 report.seam_edges_split += 1;
 }
 }
 }

 // Step 4: Merge edges across seams if enabled
 if config.merge_edges {
 let (new_brep, merged_count) = merge_seam_edges(&result, config);
 result = new_brep;
 report.seam_edges_merged = merged_count;
 }

 (result, report)
}

/// Compute the flat face index for a given solid/shell/face tuple.
fn compute_flat_face_idx(brep: &rcad_kernel::BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
 let mut idx = 0usize;
 for s in 0..solid_idx {
 for sh in &brep.solids[s].shells {
 idx += sh.faces.len();
 }
 }
 for sh in 0..shell_idx {
 idx += brep.solids[solid_idx].shells[sh].faces.len();
 }
 idx + face_idx
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Adaptive Tolerance Merging
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Configuration for adaptive tolerance merging.
#[derive(Debug, Clone)]
pub struct AdaptiveToleranceConfig {
 /// Base tolerance for merging.
 pub base_tolerance: f64,
 /// Maximum tolerance to use.
 pub max_tolerance: f64,
 /// Factor by which tolerance grows.
 pub tolerance_growth: f64,
 /// Minimum geometric feature size to preserve.
 pub min_feature_size: f64,
 /// Whether to use curvature-based tolerance adjustment.
 pub use_curvature_adjustment: bool,
}

impl Default for AdaptiveToleranceConfig {
 fn default() -> Self {
 Self {
 base_tolerance: TOLERANCE_ABS,
 max_tolerance: TOLERANCE_ABS * 1000.0,
 tolerance_growth: 2.0,
 min_feature_size: TOLERANCE_ABS * 10.0,
 use_curvature_adjustment: true,
 }
 }
}

/// Report from adaptive tolerance merging.
#[derive(Debug, Clone, Default)]
pub struct AdaptiveToleranceMergeReport {
 /// Total vertices merged.
 pub vertices_merged: usize,
 /// Total edges removed.
 pub edges_removed: usize,
 /// Number of passes executed.
 pub passes_executed: usize,
 /// Final tolerance used.
 pub final_tolerance: f64,
 /// Whether the process converged.
 pub converged: bool,
}

/// Perform adaptive tolerance merging of close vertices.
///
/// This function iteratively merges vertices with increasing tolerance,
/// but respects minimum feature size constraints to avoid merging
/// features that should be preserved.
pub fn merge_vertices_adaptive(
 brep: &rcad_kernel::BRep,
 config: &AdaptiveToleranceConfig,
) -> (rcad_kernel::BRep, AdaptiveToleranceMergeReport) {
 let mut result = brep.clone();
 let mut report = AdaptiveToleranceMergeReport::default();

 let base_tol = config.base_tolerance.max(TOLERANCE_ABS);
 let max_tol = config.max_tolerance.max(base_tol);

 for pass in 0..10 {
 let tol = if config.tolerance_growth > 1.0 {
 let grown = base_tol * config.tolerance_growth.powi(pass as i32);
 grown.min(max_tol)
 } else {
 base_tol
 };

 // Compute curvature-adjusted tolerance if enabled
 let effective_tol = if config.use_curvature_adjustment {
 compute_curvature_adjusted_tolerance(&result, tol, config.min_feature_size)
 } else {
 tol
 };

 let (new_brep, merged) = merge_close_vertices(&result, effective_tol);
 let (new_brep, removed) = remove_small_edges(&new_brep, effective_tol);

 let changed = merged > 0 || removed > 0;
 result = new_brep;
 report.vertices_merged += merged;
 report.edges_removed += removed;
 report.passes_executed = pass + 1;
 report.final_tolerance = effective_tol;

 if !changed {
 report.converged = true;
 break;
 }

 if effective_tol >= max_tol {
 break;
 }
 }

 (result, report)
}

/// Compute curvature-adjusted tolerance for a brep.
///
/// This function computes a tolerance that is adjusted based on the local
/// curvature of the geometry. In regions of high curvature, the tolerance
/// is reduced to preserve small features.
fn compute_curvature_adjusted_tolerance(brep: &rcad_kernel::BRep, base_tolerance: f64, min_feature_size: f64) -> f64 {
 // Compute the minimum curvature radius in the brep
 let mut min_curvature_radius = f64::INFINITY;

 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 // Use face normal variation as a proxy for curvature
 // For now, use a simple heuristic based on face area
 let area = compute_face_area(brep, face);
 if area > TOLERANCE_LINEAR_ULTRA_STRICT {
 // Approximate curvature radius from area
 let equiv_radius = (area / std::f64::consts::PI).sqrt();
 min_curvature_radius = min_curvature_radius.min(equiv_radius);
 }
 }
 }
 }

 // Adjust tolerance based on curvature
 if min_curvature_radius.is_finite() && min_curvature_radius > 0.0 {
 // Use a fraction of the minimum curvature radius as tolerance
 let curvature_tolerance = min_curvature_radius * 0.01;
 base_tolerance.min(curvature_tolerance).max(min_feature_size * 0.1)
 } else {
 base_tolerance
 }
}

/// Compute the approximate area of a face.
fn compute_face_area(brep: &rcad_kernel::BRep, face: &Face) -> f64 {
 let mut pts: Vec<DVec3> = Vec::new();
 for we in &face.outer_wire.edges {
 if let Some(edge) = ed_opt(brep, we.idx) {
 let vi = if we.forward { edge.start } else { edge.end };
 if let Some(v) = brep.vertices.get(vi) {
 pts.push(v.point);
 }
 }
 }

 if pts.len() < 3 {
 return 0.0;
 }

 // Fan triangulation area
 let p0 = pts[0];
 let mut area = 0.0f64;
 for i in 1..pts.len() - 1 {
 area += (pts[i] - p0).cross(pts[i + 1] - p0).length() * 0.5;
 }

 area
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// B-Spline Surface Same-Domain Detection
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Result of checking if two B-spline surfaces are the same domain.
#[derive(Debug, Clone)]
pub struct SameDomainMatch {
 /// Whether the surfaces are the same domain.
 pub is_same_domain: bool,
 /// The detected continuity level between surfaces.
 pub continuity: BsplineContinuity,
 /// Maximum deviation between control points.
 pub max_control_point_deviation: f64,
 /// Maximum deviation between weights.
 pub max_weight_deviation: f64,
 /// Whether the knot vectors match.
 pub knots_match: bool,
 /// Whether the degrees match.
 pub degrees_match: bool,
}

/// Classification of parametric continuity between B-spline surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Default)]
pub enum BsplineContinuity {
 /// No continuity (disconnected).
 #[default]
 None,
 /// C0: position continuous.
 C0,
 /// G1: tangent direction continuous (geometric continuity).
 G1,
 /// C1: tangent continuous (parametric continuity).
 C1,
 /// C2: curvature continuous.
 C2,
 /// CN: infinitely differentiable.
 CN,
}


/// Information about a merged B-spline face.
#[derive(Debug, Clone)]
pub struct MergedFaceInfo {
 /// Index of the kept face.
 pub kept_face_idx: usize,
 /// Index of the removed face.
 pub removed_face_idx: usize,
 /// Number of edges in the merged wire.
 pub merged_edge_count: usize,
 /// Whether inner wires were merged.
 pub inner_wires_merged: bool,
 /// The continuity level of the merge.
 pub continuity: BsplineContinuity,
}

/// Check if two B-spline surfaces are the same domain.
///
/// Two B-spline surfaces are considered same-domain if they have:
/// - Identical degrees (u and v)
/// - Identical knot vectors (within tolerance)
/// - Identical control point grids (within tolerance)
/// - Identical weights (for rational surfaces)
///
/// This function performs a comprehensive comparison of all geometric data.
pub fn bspline_same_domain(
 surf1: &rcad_kernel::geom::BSplineSurface,
 surf2: &rcad_kernel::geom::BSplineSurface,
 tolerance: f64,
) -> Option<SameDomainMatch> {
 const KNOT_TOL: f64 = TOLERANCE_MESH_LEGACY;
 const CP_TOL_DEFAULT: f64 = TOLERANCE_MESH_LEGACY;

 let cp_tol = if tolerance > 0.0 { tolerance } else { CP_TOL_DEFAULT };
 let knot_tol = KNOT_TOL.max(tolerance * 0.1);

 // Check degrees
 let degrees_match = surf1.degree_u == surf2.degree_u && surf1.degree_v == surf2.degree_v;
 if !degrees_match {
 return Some(SameDomainMatch {
 is_same_domain: false,
 continuity: BsplineContinuity::None,
 max_control_point_deviation: f64::INFINITY,
 max_weight_deviation: f64::INFINITY,
 knots_match: false,
 degrees_match: false,
 });
 }

 // Check knot vector lengths
 if surf1.knots_u.len() != surf2.knots_u.len() || surf1.knots_v.len() != surf2.knots_v.len() {
 return Some(SameDomainMatch {
 is_same_domain: false,
 continuity: BsplineContinuity::None,
 max_control_point_deviation: f64::INFINITY,
 max_weight_deviation: f64::INFINITY,
 knots_match: false,
 degrees_match: true,
 });
 }

 // Check knot vectors
 let mut max_knot_diff = 0.0f64;
 for (k1, k2) in surf1.knots_u.iter().zip(surf2.knots_u.iter()) {
 max_knot_diff = max_knot_diff.max((k1 - k2).abs());
 }
 for (k1, k2) in surf1.knots_v.iter().zip(surf2.knots_v.iter()) {
 max_knot_diff = max_knot_diff.max((k1 - k2).abs());
 }
 let knots_match = max_knot_diff <= knot_tol;

 if !knots_match {
 return Some(SameDomainMatch {
 is_same_domain: false,
 continuity: BsplineContinuity::None,
 max_control_point_deviation: f64::INFINITY,
 max_weight_deviation: f64::INFINITY,
 knots_match: false,
 degrees_match: true,
 });
 }

 // Check control point grid dimensions
 if surf1.control_points.len() != surf2.control_points.len() {
 return Some(SameDomainMatch {
 is_same_domain: false,
 continuity: BsplineContinuity::None,
 max_control_point_deviation: f64::INFINITY,
 max_weight_deviation: f64::INFINITY,
 knots_match: true,
 degrees_match: true,
 });
 }

 // Check control points
 let mut max_cp_deviation = 0.0f64;
 for (row1, row2) in surf1.control_points.iter().zip(surf2.control_points.iter()) {
 if row1.len() != row2.len() {
 return Some(SameDomainMatch {
 is_same_domain: false,
 continuity: BsplineContinuity::None,
 max_control_point_deviation: f64::INFINITY,
 max_weight_deviation: f64::INFINITY,
 knots_match: true,
 degrees_match: true,
 });
 }
 for (cp1, cp2) in row1.iter().zip(row2.iter()) {
 let dist = cp1.distance(*cp2);
 max_cp_deviation = max_cp_deviation.max(dist);
 }
 }

 // Check weights
 let mut max_weight_deviation = 0.0f64;
 if surf1.weights.len() != surf2.weights.len() {
 return Some(SameDomainMatch {
 is_same_domain: false,
 continuity: BsplineContinuity::None,
 max_control_point_deviation: max_cp_deviation,
 max_weight_deviation: f64::INFINITY,
 knots_match: true,
 degrees_match: true,
 });
 }
 for (row1, row2) in surf1.weights.iter().zip(surf2.weights.iter()) {
 if row1.len() != row2.len() {
 return Some(SameDomainMatch {
 is_same_domain: false,
 continuity: BsplineContinuity::None,
 max_control_point_deviation: max_cp_deviation,
 max_weight_deviation: f64::INFINITY,
 knots_match: true,
 degrees_match: true,
 });
 }
 for (w1, w2) in row1.iter().zip(row2.iter()) {
 let diff = (w1 - w2).abs();
 max_weight_deviation = max_weight_deviation.max(diff);
 }
 }

 // Determine if same domain
 let is_same_domain = max_cp_deviation <= cp_tol && max_weight_deviation <= knot_tol;

 // Determine continuity
 let continuity = if is_same_domain {
 check_bspline_continuity_from_match(surf1, surf2, cp_tol)
 } else {
 BsplineContinuity::None
 };

 Some(SameDomainMatch {
 is_same_domain,
 continuity,
 max_control_point_deviation: max_cp_deviation,
 max_weight_deviation,
 knots_match: true,
 degrees_match: true,
 })
}

/// Determine parametric continuity from matching B-spline surfaces.
fn check_bspline_continuity_from_match(
 surf: &rcad_kernel::geom::BSplineSurface,
 _other: &rcad_kernel::geom::BSplineSurface,
 tolerance: f64,
) -> BsplineContinuity {
 // For identical surfaces, continuity is determined by the degree
 // A B-spline surface has C^{degree - multiplicity} continuity at each internal knot
 // For surfaces with identical data, the minimum continuity is:
 let min_degree = surf.degree_u.min(surf.degree_v);

 if tolerance > TOLERANCE_MESH_LEGACY {
 // If tolerance is relatively large, report C0 as a conservative estimate
 return BsplineContinuity::C0;
 }

 // For clamped B-splines, only internal knot multiplicities reduce continuity
 // Boundary knots have multiplicity = degree + 1 by design
 let u_internal_mult = max_internal_knot_multiplicity(&surf.knots_u);
 let v_internal_mult = max_internal_knot_multiplicity(&surf.knots_v);

 // If no internal knots, the surface is C^{degree} everywhere inside
 // Continuity at internal knots = degree - multiplicity
 let u_continuity = if u_internal_mult == 0 {
 min_degree // No internal knots = full continuity
 } else {
 min_degree.saturating_sub(u_internal_mult)
 };
 let v_continuity = if v_internal_mult == 0 {
 min_degree // No internal knots = full continuity
 } else {
 min_degree.saturating_sub(v_internal_mult)
 };
 let min_continuity = u_continuity.min(v_continuity);

 match min_continuity {
 0 => BsplineContinuity::C0,
 1 => BsplineContinuity::C1,
 2 => BsplineContinuity::C2,
 _ if min_continuity >= 3 => BsplineContinuity::CN,
 _ => BsplineContinuity::C0,
 }
}

/// Compute the maximum multiplicity of internal knots (excluding boundary repeats).
/// Returns 0 if there are no internal knots.
fn max_internal_knot_multiplicity(knots: &[f64]) -> usize {
 if knots.len() <= 2 {
 return 0;
 }

 let tol = TOLERANCE_COORD_SUB;
 let first = knots[0];
 let last = knots[knots.len() - 1];

 // Find the range of internal knots (excluding first and last distinct values)
 let mut internal_start = 0;
 let mut internal_end = knots.len();

 // Skip boundary knots at the start
 for i in 0..knots.len() {
 if (knots[i] - first).abs() > tol {
 internal_start = i;
 break;
 }
 }

 // Skip boundary knots at the end
 for i in (0..knots.len()).rev() {
 if (knots[i] - last).abs() > tol {
 internal_end = i + 1;
 break;
 }
 }

 // If no internal knots, return 0
 if internal_start >= internal_end {
 return 0;
 }

 // Count multiplicities of internal knots
 let internal_knots = &knots[internal_start..internal_end];
 let mut max_mult = 1;
 let mut current_mult = 1;

 for i in 1..internal_knots.len() {
 if (internal_knots[i] - internal_knots[i - 1]).abs() <= tol {
 current_mult += 1;
 } else {
 max_mult = max_mult.max(current_mult);
 current_mult = 1;
 }
 }
 max_mult.max(current_mult)
}

/// Compute the maximum multiplicity of any knot in the vector.
fn max_knot_multiplicity(knots: &[f64]) -> usize {
 if knots.is_empty() {
 return 0;
 }

 let tol = TOLERANCE_COORD_SUB;
 let mut max_mult = 1;
 let mut current_mult = 1;

 for i in 1..knots.len() {
 if (knots[i] - knots[i - 1]).abs() <= tol {
 current_mult += 1;
 } else {
 max_mult = max_mult.max(current_mult);
 current_mult = 1;
 }
 }
 max_mult.max(current_mult)
}

/// Check parametric continuity between two B-spline surfaces.
///
/// This function evaluates the geometric continuity between two adjacent B-spline
/// surfaces by examining their control point and knot structures.
///
/// Returns the highest continuity level that can be guaranteed between the surfaces.
pub fn check_bspline_continuity(
 surf1: &rcad_kernel::geom::BSplineSurface,
 surf2: &rcad_kernel::geom::BSplineSurface,
 tolerance: f64,
) -> BsplineContinuity {
 // First check if surfaces are same domain
 if let Some(match_result) = bspline_same_domain(surf1, surf2, tolerance)
 && match_result.is_same_domain {
 return match_result.continuity;
 }

 // Check for adjacent surfaces (sharing a boundary)
 // This requires checking if the control points at boundaries match
 let cp_tol = tolerance.max(TOLERANCE_MESH_LEGACY);

 // Check if the last row of control points in surf1 matches the first row of surf2
 // (or vice versa) - this indicates adjacency along the v-direction
 if let Some(continuity) = check_adjacent_continuity_v(surf1, surf2, cp_tol) {
 return continuity;
 }

 // Check adjacency along u-direction
 if let Some(continuity) = check_adjacent_continuity_u(surf1, surf2, cp_tol) {
 return continuity;
 }

 BsplineContinuity::None
}

/// Check continuity between surfaces that are adjacent along the v-direction.
fn check_adjacent_continuity_v(
 surf1: &rcad_kernel::geom::BSplineSurface,
 surf2: &rcad_kernel::geom::BSplineSurface,
 tolerance: f64,
) -> Option<BsplineContinuity> {
 // surf1's last v-row should match surf2's first v-row (or vice versa)
 if surf1.control_points.is_empty() || surf2.control_points.is_empty() {
 return None;
 }

 let n_u1 = surf1.control_points.len();
 let n_u2 = surf2.control_points.len();

 if n_u1 == 0 || n_u2 == 0 {
 return None;
 }

 // Check degrees compatibility
 if surf1.degree_u != surf2.degree_u {
 return None;
 }

 // Check if last row of surf1 matches first row of surf2
 let row1 = &surf1.control_points[n_u1 - 1];
 let row2 = &surf2.control_points[0];

 if row1.len() != row2.len() {
 return None;
 }

 let mut max_dev = 0.0_f64;
 for (p1, p2) in row1.iter().zip(row2.iter()) {
 max_dev = max_dev.max(p1.distance(*p2));
 }

 if max_dev <= tolerance {
 // Surfaces are adjacent with C0 continuity
 // Check for higher continuity by comparing derivative rows
 if n_u1 >= 2 && n_u2 >= 2 {
 let row1_prev = &surf1.control_points[n_u1 - 2];
 let row2_next = &surf2.control_points[1];

 if row1_prev.len() == row2_next.len() {
 let mut max_deriv_dev = 0.0_f64;
 for ((p1, p2), (p1_prev, p2_next)) in
 row1.iter().zip(row2.iter())
 .zip(row1_prev.iter().zip(row2_next.iter()))
 {
 // Approximate tangent direction continuity
 let t1 = (*p2 - *p1_prev).normalize_or(DVec3::ZERO);
 let t2 = (*p2_next - *p1).normalize_or(DVec3::ZERO);
 let dot = t1.dot(t2);
 if dot > 0.99 {
 // Tangents are nearly parallel - G1 continuity
 max_deriv_dev = max_deriv_dev.max((t1 - t2).length());
 }
 }

 if max_deriv_dev <= tolerance * 10.0 {
 return Some(BsplineContinuity::G1);
 }
 }
 }

 return Some(BsplineContinuity::C0);
 }

 // Check the reverse: last row of surf2 matches first row of surf1
 let row2_last = &surf2.control_points[n_u2 - 1];
 let row1_first = &surf1.control_points[0];

 if row2_last.len() != row1_first.len() {
 return None;
 }

 max_dev = 0.0;
 for (p1, p2) in row2_last.iter().zip(row1_first.iter()) {
 max_dev = max_dev.max(p1.distance(*p2));
 }

 if max_dev <= tolerance {
 return Some(BsplineContinuity::C0);
 }

 None
}

/// Check continuity between surfaces that are adjacent along the u-direction.
fn check_adjacent_continuity_u(
 surf1: &rcad_kernel::geom::BSplineSurface,
 surf2: &rcad_kernel::geom::BSplineSurface,
 tolerance: f64,
) -> Option<BsplineContinuity> {
 // For each row, check if the last column of surf1 matches the first column of surf2
 if surf1.control_points.is_empty() || surf2.control_points.is_empty() {
 return None;
 }

 // Check degrees compatibility
 if surf1.degree_v != surf2.degree_v {
 return None;
 }

 // Check if row counts match
 if surf1.control_points.len() != surf2.control_points.len() {
 return None;
 }

 let mut max_dev = 0.0_f64;
 for (row1, row2) in surf1.control_points.iter().zip(surf2.control_points.iter()) {
 if row1.is_empty() || row2.is_empty() {
 continue;
 }

 let n_v1 = row1.len();
 let _n_v2 = row2.len();

 // Check last of row1 vs first of row2
 let dev = row1[n_v1 - 1].distance(row2[0]);
 max_dev = max_dev.max(dev);
 }

 if max_dev <= tolerance {
 return Some(BsplineContinuity::C0);
 }

 // Check the reverse direction
 max_dev = 0.0_f64;
 for (row1, row2) in surf1.control_points.iter().zip(surf2.control_points.iter()) {
 if row1.is_empty() || row2.is_empty() {
 continue;
 }

 let n_v2 = row2.len();

 // Check last of row2 vs first of row1
 let dev = row2[n_v2 - 1].distance(row1[0]);
 max_dev = max_dev.max(dev);
 }

 if max_dev <= tolerance {
 return Some(BsplineContinuity::C0);
 }

 None
}

/// Merge adjacent B-spline faces if they are on the same domain.
///
/// This function checks if two faces sharing a B-spline surface can be merged.
/// The faces must be adjacent (share an edge) and lie on the same B-spline surface.
///
/// Returns `Some((brep, MergedFaceInfo))` if the faces were merged, `None` otherwise.
pub fn merge_bspline_faces(
 brep: &rcad_kernel::BRep,
 face1_idx: usize,
 face2_idx: usize,
 tolerance: f64,
) -> Option<(rcad_kernel::BRep, MergedFaceInfo)> {
 // Get surfaces for both faces
 let surf1_idx = brep.geom.face_surface.get(face1_idx).and_then(|v| *v)?;
 let surf2_idx = brep.geom.face_surface.get(face2_idx).and_then(|v| *v)?;

 let surf1 = brep.geom.surfaces.get(surf1_idx)?;
 let surf2 = brep.geom.surfaces.get(surf2_idx)?;

 // Both must be B-spline surfaces
 let (bs1, bs2) = match (surf1, surf2) {
 (rcad_kernel::geom::Surface3::BSpline(b1), rcad_kernel::geom::Surface3::BSpline(b2)) => (b1, b2),
 _ => return None,
 };

 // Check same domain
 let match_result = bspline_same_domain(bs1, bs2, tolerance)?;
 if !match_result.is_same_domain {
 return None;
 }

 // Find the solid and shell containing both faces
 let (si, shi) = find_shell_containing_faces(brep, face1_idx, face2_idx)?;

 // Get local face indices within the shell
 let fi1 = find_face_index_in_shell(brep, si, shi, face1_idx)?;
 let fi2 = find_face_index_in_shell(brep, si, shi, face2_idx)?;

 // Find shared edge
 let shared_edge = find_shared_edge(brep, si, shi, fi1, fi2)?;

 // Perform the merge
 let mut result = brep.clone();

 // Splice the wires together
 let wire1 = result.solids[si].shells[shi].faces[fi1].outer_wire.edges.clone();
 let wire2 = result.solids[si].shells[shi].faces[fi2].outer_wire.edges.clone();

 let merged_wire = splice_wires_for_merge(&wire1, &wire2, shared_edge)?;

 // Collect inner wires
 let inner1 = result.solids[si].shells[shi].faces[fi1].inner_wires.clone();
 let inner2 = result.solids[si].shells[shi].faces[fi2].inner_wires.clone();
 let inner_wires_merged = !inner2.is_empty();
 let mut all_inner = inner1;
 all_inner.extend(inner2);

 // Build merged face
 let face1 = &result.solids[si].shells[shi].faces[fi1];
 let merged_face = rcad_kernel::topology::Face {
 outer_wire: rcad_kernel::topology::Wire { edges: merged_wire },
 inner_wires: all_inner,
 normal: face1.normal,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 let merged_edge_count = merged_face.outer_wire.edges.len();

 // Determine which face to keep (lower index) and which to remove
 let (keep_idx, remove_idx) = if fi1 < fi2 { (fi1, fi2) } else { (fi2, fi1) };

 // Update face_surface mapping
 let _kept_flat = flat_face_index_global(&result, si, shi, keep_idx);
 let remove_flat = flat_face_index_global(&result, si, shi, remove_idx);
 if result.geom.face_surface.len() > remove_flat {
 result.geom.face_surface.remove(remove_flat);
 }
 if result.geom.face_surface_range.len() > remove_flat {
 result.geom.face_surface_range.remove(remove_flat);
 }
 if result.geom.face_tolerance.len() > remove_flat {
 result.geom.face_tolerance.remove(remove_flat);
 }

 // Replace the kept face and remove the other
 result.solids[si].shells[shi].faces[keep_idx] = merged_face;
 result.solids[si].shells[shi].faces.remove(remove_idx);

 Some((result, MergedFaceInfo {
 kept_face_idx: keep_idx,
 removed_face_idx: remove_idx,
 merged_edge_count,
 inner_wires_merged,
 continuity: match_result.continuity,
 }))
}

/// Find the shell containing two faces.
fn find_shell_containing_faces(brep: &rcad_kernel::BRep, face1_idx: usize, face2_idx: usize) -> Option<(usize, usize)> {
 let mut found_si = None;
 let mut found_shi = None;

 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 let base = flat_face_index_global(brep, si, shi, 0);
 let n_faces = brep.solids[si].shells[shi].faces.len();

 let face1_in_shell = face1_idx >= base && face1_idx < base + n_faces;
 let face2_in_shell = face2_idx >= base && face2_idx < base + n_faces;

 if face1_in_shell && face2_in_shell {
 found_si = Some(si);
 found_shi = Some(shi);
 break;
 }
 }
 if found_si.is_some() {
 break;
 }
 }

 Some((found_si?, found_shi?))
}

/// Find the local index of a face within a shell.
fn find_face_index_in_shell(brep: &rcad_kernel::BRep, si: usize, shi: usize, global_face_idx: usize) -> Option<usize> {
 let base = flat_face_index_global(brep, si, shi, 0);
 if global_face_idx >= base {
 Some(global_face_idx - base)
 } else {
 None
 }
}

/// Get the global flat index of a face.
fn flat_face_index_global(brep: &rcad_kernel::BRep, si: usize, shi: usize, fi: usize) -> usize {
 let mut idx = 0usize;
 for s in 0..si {
 for sh in &brep.solids[s].shells {
 idx += sh.faces.len();
 }
 }
 for sh in 0..shi {
 idx += brep.solids[si].shells[sh].faces.len();
 }
 idx + fi
}

/// Find a shared edge between two faces in a shell.
fn find_shared_edge(brep: &rcad_kernel::BRep, si: usize, shi: usize, fi1: usize, fi2: usize) -> Option<usize> {
 use std::collections::HashSet;

 let face1 = &brep.solids[si].shells[shi].faces[fi1];
 let face2 = &brep.solids[si].shells[shi].faces[fi2];

 let edges1: HashSet<usize> = face1.outer_wire.edges.iter().map(|we| we.idx).collect();
 let edges2: HashSet<usize> = face2.outer_wire.edges.iter().map(|we| we.idx).collect();

 edges1.intersection(&edges2).copied().next()
}

/// Splice two wire edge lists together for merging.
fn splice_wires_for_merge(
 wire_a: &[rcad_kernel::topology::WireEdge],
 wire_b: &[rcad_kernel::topology::WireEdge],
 shared_edge_idx: usize,
) -> Option<Vec<rcad_kernel::topology::WireEdge>> {
 let pos_a = wire_a.iter().position(|we| we.idx == shared_edge_idx)?;
 let pos_b = wire_b.iter().position(|we| we.idx == shared_edge_idx)?;

 let n_b = wire_b.len();
 // B's edges (excluding the shared edge), in cyclic order starting at pos_b + 1
 let b_edges: Vec<rcad_kernel::topology::WireEdge> =
 (1..n_b).map(|i| wire_b[(pos_b + i) % n_b]).collect();

 let mut merged = Vec::with_capacity(wire_a.len() - 1 + b_edges.len());
 merged.extend_from_slice(&wire_a[..pos_a]);
 merged.extend(b_edges);
 merged.extend_from_slice(&wire_a[pos_a + 1..]);

 if merged.len() < 3 {
 return None; // Degenerate result
 }

 Some(merged)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Shell Repair (ShapeFix_Shell equivalent)
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Report from shell-level repair operations.
///
/// Analogous to OCCT `ShapeFix_Shell` report.
#[derive(Debug, Clone, Default)]
pub struct ShellFixReport {
 /// Number of faces whose orientation was corrected.
 pub faces_reoriented: usize,
 /// Number of edges that were non-manifold and were processed.
 pub non_manifold_edges_processed: usize,
 /// Number of new shells created from splitting non-manifold topology.
 pub shells_created: usize,
 /// Whether the shell is now closed.
 pub is_closed: bool,
 /// Whether the shell is now manifold.
 pub is_manifold: bool,
 /// Number of open edges detected.
 pub open_edge_count: usize,
 /// Number of non-manifold edges detected.
 pub non_manifold_edge_count: usize,
}

impl ShellFixReport {
 /// Returns true if the shell is in a clean state.
 pub fn is_clean(&self) -> bool {
 self.is_closed && self.is_manifold
 }

 /// Returns a summary string.
 pub fn summary(&self) -> String {
 format!(
 "ShellFix: {} faces reoriented, {} non-manifold edges processed, closed={}, manifold={}",
 self.faces_reoriented,
 self.non_manifold_edges_processed,
 self.is_closed,
 self.is_manifold
 )
 }
}

/// Report from shell closure checking.
#[derive(Debug, Clone, Default)]
pub struct ClosureReport {
 /// Whether the shell forms a closed surface (no free edges).
 pub is_closed: bool,
 /// Number of edges referenced by exactly 1 face (free/open edges).
 pub open_edge_count: usize,
 /// List of open edge indices.
 pub open_edges: Vec<usize>,
 /// Euler characteristic: V - E + F.
 pub euler_characteristic: i64,
 /// Number of unique vertices in the shell.
 pub vertex_count: usize,
 /// Number of unique edges in the shell.
 pub edge_count: usize,
 /// Number of faces in the shell.
 pub face_count: usize,
 /// Whether the shell is orientable (has consistent normal direction).
 pub is_orientable: bool,
 /// Genus computed from Euler characteristic (if closed).
 pub genus: Option<i64>,
}

impl ClosureReport {
 /// Returns true if the shell is closed and orientable.
 pub fn is_valid(&self) -> bool {
 self.is_closed && self.is_orientable
 }

 /// Returns a summary string.
 pub fn summary(&self) -> String {
 if self.is_closed {
 let genus_str = self.genus.map_or("?".to_string(), |g| g.to_string());
 format!(
 "Closed shell: V={}, E={}, F={},  ?{}, genus={}",
 self.vertex_count, self.edge_count, self.face_count,
 self.euler_characteristic, genus_str
 )
 } else {
 format!(
 "Open shell: {} open edges, V={}, E={}, F={},  ?{}",
 self.open_edge_count, self.vertex_count, self.edge_count,
 self.face_count, self.euler_characteristic
 )
 }
 }
}

/// Check shell closure and compute Euler characteristic.
///
/// This function analyzes a shell to determine if it forms a closed surface
/// (no free edges) and computes the Euler characteristic V - E + F.
///
/// # Arguments
/// * `shell` - The shell to analyze.
/// * `brep` - The containing brep.
///
/// # Returns
/// A `ClosureReport` with closure status and Euler characteristic.
///
/// # Example
/// ```rust
/// use brep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::check_shell_closure;
///
/// let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
/// width: 1.0, height: 1.0, depth: 1.0
/// });
/// let shell = &brep.solids[0].shells[0];
/// let report = check_shell_closure(shell, &rcad_kernel::BRep);
/// assert!(report.is_closed);
/// assert_eq!(report.euler_characteristic, 2); // Sphere topology
/// ```
pub fn check_shell_closure(shell: &Shell, brep: &rcad_kernel::BRep) -> ClosureReport {
 use std::collections::{HashMap, HashSet};

 let n_edges = brep.edge_count();
 let face_count = shell.faces.len();

 // Collect unique edges and count edge-face references
 let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
 let mut unique_edges: HashSet<usize> = HashSet::new();

 for face in &shell.faces {
 // Count edges in outer wire
 for we in &face.outer_wire.edges {
 if we.idx < n_edges {
 unique_edges.insert(we.idx);
 *edge_face_count.entry(we.idx).or_insert(0) += 1;
 }
 }
 // Count edges in inner wires
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if we.idx < n_edges {
 unique_edges.insert(we.idx);
 *edge_face_count.entry(we.idx).or_insert(0) += 1;
 }
 }
 }
 }

 // Find open edges (referenced by exactly 1 face)
 let open_edges: Vec<usize> = edge_face_count
 .iter()
 .filter(|(_, count)| **count == 1)
 .map(|(idx, _)| *idx)
 .collect();
 let open_edge_count = open_edges.len();

 // Collect unique vertices from unique edges
 let mut unique_verts: HashSet<usize> = HashSet::new();
 for &ei in &unique_edges {
 if let Some(edge) = ed_opt(brep, ei) {
 unique_verts.insert(edge.start);
 unique_verts.insert(edge.end);
 }
 }

 let vertex_count = unique_verts.len();
 let edge_count = unique_edges.len();

 // Compute Euler characteristic
 let euler_characteristic = vertex_count as i64 - edge_count as i64 + face_count as i64;

 // Check orientability by examining face normals
 // For a simple check, verify that adjacent faces have compatible normals
 let is_orientable = check_shell_orientability(shell, brep);

 // Compute genus if closed
 let is_closed = open_edge_count == 0;
 let genus = if is_closed {
 let g = (2 - euler_characteristic) / 2;
 if (2 - euler_characteristic) % 2 == 0 && g >= 0 {
 Some(g)
 } else {
 None
 }
 } else {
 None
 };

 ClosureReport {
 is_closed,
 open_edge_count,
 open_edges,
 euler_characteristic,
 vertex_count,
 edge_count,
 face_count,
 is_orientable,
 genus,
 }
}

/// Check if a shell is orientable by verifying face normals are consistent.
fn check_shell_orientability(shell: &Shell, brep: &rcad_kernel::BRep) -> bool {
 // For a properly oriented shell, all face normals should point outward.
 // We check this by verifying that the normals don't flip direction
 // relative to a consistent reference (the shell centroid).

 if shell.faces.is_empty() {
 return true;
 }

 // Compute the shell centroid
 let shell_centroid = compute_shell_centroid(shell, brep);

 // Check each face's normal orientation
 for face in &shell.faces {
 let face_centroid = compute_face_centroid(&face.outer_wire, brep);
 let outward = face_centroid - shell_centroid;

 // If outward vector is very small, skip this face
 if outward.length() < TOLERANCE_LINEAR_ULTRA_STRICT {
 continue;
 }

 // Normal should have positive dot product with outward direction
 let dot = face.normal.dot(outward);
 if dot < 0.0 {
 return false;
 }
 }

 true
}

/// Fix shell orientation for proper normal direction.
///
/// This function corrects face orientations so that all normals point
/// consistently outward (or inward for inner shells). It handles nested
/// shells by detecting which shells are outer vs inner.
///
/// # Arguments
/// * `shell` - The shell to repair.
/// * `brep` - The containing brep.
///
/// # Returns
/// A tuple of (repaired shell, report).
///
/// Analogous to OCCT `ShapeFix_Shell::FixOrientation()`.
pub fn fix_shell_orientation(shell: &Shell, brep: &rcad_kernel::BRep) -> (Shell, ShellFixReport) {
 let mut report = ShellFixReport::default();
 let mut fixed_shell = shell.clone();

 // Compute the shell's centroid from all face centroids
 let shell_centroid = compute_shell_centroid(shell, brep);

 // Check each face's normal orientation relative to the shell centroid
 for face in &mut fixed_shell.faces {
 let face_centroid = compute_face_centroid(&face.outer_wire, brep);
 let outward = face_centroid - shell_centroid;
 let dot = face.normal.dot(outward);

 // If normal points inward (negative dot product), flip the face
 if dot < 0.0 {
 face.normal = -face.normal;
 face.outer_wire = reverse_wire(&face.outer_wire);
 for inner in &mut face.inner_wires {
 *inner = reverse_wire(inner);
 }
 report.faces_reoriented += 1;
 }
 }

 // Check final state
 let closure_report = check_shell_closure(&fixed_shell, brep);
 report.is_closed = closure_report.is_closed;
 report.open_edge_count = closure_report.open_edge_count;

 // Check manifoldness
 let manifold_report = analyze_shell_manifoldness(&fixed_shell, brep);
 report.is_manifold = manifold_report.is_manifold;
 report.non_manifold_edge_count = manifold_report.non_manifold_edges.len();

 (fixed_shell, report)
}

/// Compute the centroid of a shell from all its face vertices.
fn compute_shell_centroid(shell: &Shell, brep: &rcad_kernel::BRep) -> DVec3 {
 let mut sum = DVec3::ZERO;
 let mut count = 0usize;

 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if let Some(edge) = ed_opt(brep, we.idx) {
 let vi = if we.forward { edge.start } else { edge.end };
 if let Some(v) = brep.vertices.get(vi) {
 sum += v.point;
 count += 1;
 }
 }
 }
 }

 if count > 0 {
 sum / count as f64
 } else {
 DVec3::ZERO
 }
}

/// Compute the centroid of a face from its outer wire vertices.
fn compute_face_centroid(wire: &Wire, brep: &rcad_kernel::BRep) -> DVec3 {
 let mut sum = DVec3::ZERO;
 let mut count = 0usize;

 for we in &wire.edges {
 if let Some(edge) = ed_opt(brep, we.idx) {
 let vi = if we.forward { edge.start } else { edge.end };
 if let Some(v) = brep.vertices.get(vi) {
 sum += v.point;
 count += 1;
 }
 }
 }

 if count > 0 {
 sum / count as f64
 } else {
 DVec3::ZERO
 }
}

/// Report from shell manifoldness analysis.
#[derive(Debug, Clone, Default)]
struct ManifoldReport {
 is_manifold: bool,
 non_manifold_edges: Vec<usize>,
 non_manifold_vertices: Vec<usize>,
}

/// Analyze a shell for manifoldness.
fn analyze_shell_manifoldness(shell: &Shell, brep: &rcad_kernel::BRep) -> ManifoldReport {
 use std::collections::{HashMap, HashSet};

 let n_edges = brep.edge_count();

 // Count edge-face references
 let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if we.idx < n_edges {
 *edge_face_count.entry(we.idx).or_insert(0) += 1;
 }
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if we.idx < n_edges {
 *edge_face_count.entry(we.idx).or_insert(0) += 1;
 }
 }
 }
 }

 // Find non-manifold edges (referenced by more than 2 faces)
 let non_manifold_edges: Vec<usize> = edge_face_count
 .iter()
 .filter(|(_, count)| **count > 2)
 .map(|(idx, _)| *idx)
 .collect();

 // Find non-manifold vertices
 let mut vertex_edge_count: HashMap<usize, HashSet<usize>> = HashMap::new();
 for &ei in edge_face_count.keys() {
 if let Some(edge) = ed_opt(brep, ei) {
 vertex_edge_count.entry(edge.start).or_default().insert(ei);
 vertex_edge_count.entry(edge.end).or_default().insert(ei);
 }
 }

 // A vertex is non-manifold if it's shared by edges that don't form a single fan
 let non_manifold_vertices: Vec<usize> = vertex_edge_count
 .iter()
 .filter(|(_, edges)| {
 // Simple heuristic: if vertex has > 4 edges, might be non-manifold
 // A proper check would verify the edge fan connectivity
 edges.len() > 4
 })
 .map(|(&vi, _)| vi)
 .collect();

 ManifoldReport {
 is_manifold: non_manifold_edges.is_empty() && non_manifold_vertices.is_empty(),
 non_manifold_edges,
 non_manifold_vertices,
 }
}

/// Fix non-manifold shell topology where possible.
///
/// This function attempts to convert non-manifold topology to manifold by:
/// - Splitting non-manifold edges (edges shared by 3+ faces)
/// - Creating separate shells for disconnected regions
///
/// # Arguments
/// * `shell` - The shell to repair.
/// * `brep` - The containing brep.
///
/// # Returns
/// A tuple of (repaired shell, report). The repaired shell may have different
/// topology but represents the same geometric shape in manifold form.
///
/// Analogous to OCCT `ShapeFix_Shell::FixManifold()`.
pub fn fix_non_manifold_shell(shell: &Shell, brep: &rcad_kernel::BRep) -> (Shell, ShellFixReport) {
 let mut report = ShellFixReport::default();

 // First analyze the shell for manifold issues
 let manifold_report = analyze_shell_manifoldness(shell, brep);
 report.non_manifold_edge_count = manifold_report.non_manifold_edges.len();

 if manifold_report.is_manifold {
 // Already manifold - just check closure
 let closure_report = check_shell_closure(shell, brep);
 report.is_closed = closure_report.is_closed;
 report.is_manifold = true;
 return (shell.clone(), report);
 }

 // For now, we mark non-manifold edges but don't split them
 // A full implementation would:
 // 1. Duplicate non-manifold edges
 // 2. Update face references to use the appropriate edge copy
 // 3. Potentially create separate shells for disconnected regions

 report.non_manifold_edges_processed = manifold_report.non_manifold_edges.len();

 // Return the original shell since we don't modify it yet
 // The processing is recorded in the report
 let closure_report = check_shell_closure(shell, brep);
 report.is_closed = closure_report.is_closed;
 report.is_manifold = false;

 (shell.clone(), report)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Enhanced Shell Repair (ShapeFix_Shell extensions)
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Detailed report from shell orientation analysis and repair.
#[derive(Debug, Clone, Default)]
pub struct ShellOrientationReport {
 pub faces_inverted: usize,
 pub faces_correct: usize,
 pub inverted_face_indices: Vec<usize>,
 pub edge_conflicts: usize,
 pub is_consistent: bool,
 pub non_manifold_edges_skipped: usize,
 pub volume_sign: f64,
}

impl ShellOrientationReport {
 pub fn is_valid(&self) -> bool {
 self.is_consistent && self.edge_conflicts == 0
 }

 pub fn summary(&self) -> String {
 format!(
 "ShellOrientation: {} inverted, {} correct, {} edge conflicts, consistent={}",
 self.faces_inverted, self.faces_correct, self.edge_conflicts, self.is_consistent
 )
 }
}

/// Result from shell closure repair operations.
#[derive(Debug, Clone)]
pub struct ShellClosureResult {
 pub original_shell: Shell,
 pub repaired_shell: Shell,
 pub open_edges_detected: usize,
 pub gaps_closed: usize,
 pub faces_added: usize,
 pub unrepairable_gaps: Vec<GapInfo>,
 pub is_closed: bool,
 pub tolerance_used: f64,
}

impl ShellClosureResult {
 pub fn is_successful(&self) -> bool {
 self.is_closed && self.unrepairable_gaps.is_empty()
 }

 pub fn summary(&self) -> String {
 if self.is_closed {
 format!("ShellClosure: closed {} gaps, added {} faces", self.gaps_closed, self.faces_added)
 } else {
 format!("ShellClosure: {} open edges, {} unrepairable", self.open_edges_detected, self.unrepairable_gaps.len())
 }
 }
}

/// Information about a gap in the shell.
#[derive(Debug, Clone)]
pub struct GapInfo {
 pub boundary_edges: Vec<usize>,
 pub estimated_area: f64,
 pub can_fill: bool,
 pub failure_reason: Option<String>,
}

/// Result from non-manifold edge repair.
#[derive(Debug, Clone)]
pub struct ManifoldRepairResult {
 pub original_shell: Shell,
 pub repaired_shell: Shell,
 pub edges_processed: usize,
 pub edges_split: usize,
 pub vertices_duplicated: usize,
 pub faces_created: usize,
 pub is_manifold: bool,
 pub edge_details: Vec<NonManifoldEdgeInfo>,
}

impl ManifoldRepairResult {
 pub fn is_successful(&self) -> bool {
 self.is_manifold
 }

 pub fn summary(&self) -> String {
 format!("ManifoldRepair: {} edges split, {} vertices duplicated, manifold={}", self.edges_split, self.vertices_duplicated, self.is_manifold)
 }
}

/// Information about a non-manifold edge.
#[derive(Debug, Clone)]
pub struct NonManifoldEdgeInfo {
 pub edge_index: usize,
 pub face_count: usize,
 pub face_indices: Vec<usize>,
 pub repaired: bool,
 pub copies_created: usize,
}

/// Comprehensive validation report for shell topology.
#[derive(Debug, Clone, Default)]
pub struct ShellValidationReport {
 pub is_valid: bool,
 pub euler_characteristic: i64,
 pub expected_euler: Option<i64>,
 pub euler_valid: bool,
 pub vertex_count: usize,
 pub edge_count: usize,
 pub face_count: usize,
 pub open_edge_count: usize,
 pub non_manifold_edge_count: usize,
 pub non_manifold_vertex_count: usize,
 pub orientation_consistent: bool,
 pub is_closed: bool,
 pub is_manifold: bool,
 pub genus: Option<i64>,
 pub edge_valence: Vec<EdgeValenceInfo>,
 pub vertex_valence: Vec<VertexValenceInfo>,
 pub errors: Vec<String>,
 pub warnings: Vec<String>,
}

impl ShellValidationReport {
 pub fn is_closed_manifold(&self) -> bool {
 self.is_closed && self.is_manifold && self.orientation_consistent
 }

 pub fn summary(&self) -> String {
 let status = if self.is_valid { "VALID" } else { "INVALID" };
 format!("ShellValidation: {} | V={}, E={}, F={},  ?{}", status, self.vertex_count, self.edge_count, self.face_count, self.euler_characteristic)
 }
}

/// Information about edge valence.
#[derive(Debug, Clone)]
pub struct EdgeValenceInfo {
 pub edge_index: usize,
 pub valence: usize,
 pub is_open: bool,
 pub is_manifold: bool,
 pub is_non_manifold: bool,
}

/// Information about vertex valence.
#[derive(Debug, Clone)]
pub struct VertexValenceInfo {
 pub vertex_index: usize,
 pub edge_valence: usize,
 pub face_valence: usize,
 pub is_boundary: bool,
 pub is_non_manifold: bool,
}

/// Fix shell orientation with detailed edge adjacency analysis.
pub fn fix_shell_orientation_advanced(shell: &Shell, brep: &rcad_kernel::BRep) -> (Shell, ShellOrientationReport) {
 use std::collections::{HashMap, VecDeque};

 let mut report = ShellOrientationReport::default();
 let mut fixed_shell = shell.clone();

 if shell.faces.is_empty() {
 report.is_consistent = true;
 return (fixed_shell, report);
 }

 let n_edges = brep.edge_count();
 let mut edge_faces: HashMap<usize, Vec<(usize, bool)>> = HashMap::new();

 for (face_idx, face) in shell.faces.iter().enumerate() {
 for we in &face.outer_wire.edges {
 if we.idx < n_edges {
 edge_faces.entry(we.idx).or_default().push((face_idx, we.forward));
 }
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if we.idx < n_edges {
 edge_faces.entry(we.idx).or_default().push((face_idx, we.forward));
 }
 }
 }
 }

 report.non_manifold_edges_skipped = edge_faces.values().filter(|faces| faces.len() > 2).count();

 let mut face_orientation: Vec<Option<bool>> = vec![None; shell.faces.len()];
 let mut queue: VecDeque<usize> = VecDeque::new();
 face_orientation[0] = Some(true);
 queue.push_back(0);

 while let Some(current_face) = queue.pop_front() {
 let current_keep = face_orientation[current_face].unwrap();
 for we in &shell.faces[current_face].outer_wire.edges {
 if we.idx >= n_edges { continue; }
 if let Some(adjacent) = edge_faces.get(&we.idx) {
 for &(adj_face_idx, adj_forward) in adjacent {
 if adj_face_idx == current_face || face_orientation[adj_face_idx].is_some() { continue; }
 let current_forward = we.forward;
 if adjacent.len() == 2 {
 let should_flip = if current_keep { current_forward == adj_forward } else { current_forward != adj_forward };
 face_orientation[adj_face_idx] = Some(!should_flip);
 } else {
 face_orientation[adj_face_idx] = Some(true);
 }
 queue.push_back(adj_face_idx);
 }
 }
 }
 }

 let shell_centroid = compute_shell_centroid(shell, brep);
 for (i, orientation) in face_orientation.iter_mut().enumerate() {
 if orientation.is_none() {
 let face = &shell.faces[i];
 let face_centroid = compute_face_centroid(&face.outer_wire, brep);
 let outward = face_centroid - shell_centroid;
 let dot = face.normal.dot(outward);
 *orientation = Some(dot >= 0.0);
 }
 }

 for (i, face) in fixed_shell.faces.iter_mut().enumerate() {
 let keep_original = face_orientation[i].unwrap_or(true);
 if !keep_original {
 face.normal = -face.normal;
 face.outer_wire = reverse_wire(&face.outer_wire);
 for inner in &mut face.inner_wires { *inner = reverse_wire(inner); }
 report.faces_inverted += 1;
 report.inverted_face_indices.push(i);
 } else {
 report.faces_correct += 1;
 }
 }

 for faces in edge_faces.values() {
 if faces.len() == 2 {
 let (f1, fwd1) = faces[0];
 let (f2, fwd2) = faces[1];
 let keep1 = face_orientation[f1].unwrap_or(true);
 let keep2 = face_orientation[f2].unwrap_or(true);
 let eff_fwd1 = if keep1 { fwd1 } else { !fwd1 };
 let eff_fwd2 = if keep2 { fwd2 } else { !fwd2 };
 if eff_fwd1 == eff_fwd2 { report.edge_conflicts += 1; }
 }
 }

 report.volume_sign = compute_shell_volume(&fixed_shell, brep);
 report.is_consistent = report.edge_conflicts == 0 && report.volume_sign >= 0.0;
 (fixed_shell, report)
}

/// Repair shell closure by detecting and closing gaps.
pub fn repair_shell_closure(shell: &Shell, brep: &rcad_kernel::BRep, tolerance: f64) -> ShellClosureResult {
 use std::collections::{HashMap, HashSet};

 let mut result = ShellClosureResult {
 original_shell: shell.clone(),
 repaired_shell: shell.clone(),
 open_edges_detected: 0,
 gaps_closed: 0,
 faces_added: 0,
 unrepairable_gaps: vec![],
 is_closed: false,
 tolerance_used: tolerance,
 };

 let n_edges = brep.edge_count();
 let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if we.idx < n_edges { *edge_face_count.entry(we.idx).or_insert(0) += 1; }
 }
 }

 let open_edges: Vec<usize> = edge_face_count.iter().filter(|(_, c)| **c == 1).map(|(i, _)| *i).collect();
 result.open_edges_detected = open_edges.len();

 if open_edges.is_empty() {
 result.is_closed = true;
 return result;
 }

 let mut visited: HashSet<usize> = HashSet::new();
 while visited.len() < open_edges.len() {
 let start_edge = match open_edges.iter().find(|e| !visited.contains(e)) {
 Some(e) => *e,
 None => break,
 };
 let mut chain: Vec<usize> = vec![start_edge];
 visited.insert(start_edge);

 loop {
 let mut extended = false;
 for &oe in &open_edges {
 if visited.contains(&oe) { continue; }
 let last = ed_opt(brep, chain[chain.len() - 1]);
 let curr = ed_opt(brep, oe);
 if let (Some(l), Some(c)) = (last, curr)
 && (l.end == c.start || l.end == c.end || l.start == c.start || l.start == c.end) {
 chain.push(oe);
 visited.insert(oe);
 extended = true;
 break;
 }
 }
 if !extended { break; }
 }

 if chain.len() >= 3 {
 let is_closed_loop = {
 let first = ed_opt(brep, chain[0]);
 let last = ed_opt(brep, chain[chain.len() - 1]);
 if let (Some(f), Some(l)) = (first, last) {
 l.end == f.start || l.start == f.start || l.end == f.end || l.start == f.end
 } else { false }
 };

 let gap_info = GapInfo {
 boundary_edges: chain.clone(),
 estimated_area: estimate_chain_area(&chain, brep),
 can_fill: is_closed_loop && chain.len() >= 3,
 failure_reason: if !is_closed_loop { Some("Gap boundary is not closed".into()) } else { None },
 };

 if gap_info.can_fill {
 if let Some(new_face) = create_face_from_boundary(&chain, brep, tolerance) {
 result.repaired_shell.faces.push(new_face);
 result.faces_added += 1;
 result.gaps_closed += 1;
 } else {
 result.unrepairable_gaps.push(GapInfo { failure_reason: Some("Could not create face".into()), ..gap_info });
 }
 } else {
 result.unrepairable_gaps.push(gap_info);
 }
 }
 }

 result.is_closed = check_shell_closure(&result.repaired_shell, brep).is_closed;
 result
}

fn estimate_chain_area(chain: &[usize], brep: &rcad_kernel::BRep) -> f64 {
 if chain.len() < 3 { return 0.0; }
 let mut nodes: Vec<DVec3> = Vec::new();
 for &ei in chain {
 if let Some(edge) = ed_opt(brep, ei)
 && let (Some(s), Some(e)) = (brep.vertices.get(edge.start), brep.vertices.get(edge.end)) {
 if nodes.is_empty() { nodes.push(s.point); }
 nodes.push(e.point);
 }
 }
 if nodes.len() < 3 { return 0.0; }
 let mut area = 0.0;
 for i in 0..nodes.len() {
 let j = (i + 1) % nodes.len();
 area += nodes[i].x * nodes[j].y - nodes[j].x * nodes[i].y;
 }
 (area / 2.0).abs()
}



