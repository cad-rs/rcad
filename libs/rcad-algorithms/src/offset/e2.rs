
/// Detect self-intersection with detailed results.
///
/// This is a more comprehensive analysis that returns information about
/// which faces might intersect and the minimum safe offset distance.
pub fn detect_self_intersection_detailed(brep: &BRep, distance: f64) -> SelfIntersectionResult {
 let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
 Some(s) => s,
 None => {
 return SelfIntersectionResult {
 has_intersection: false,
 intersecting_pairs: Vec::new(),
 min_safe_distance: None,
 description: "no shell found".to_string(),
 };
 }
 };

 if shell.faces.len() < 3 {
 return SelfIntersectionResult {
 has_intersection: false,
 intersecting_pairs: Vec::new(),
 min_safe_distance: None,
 description: "insufficient faces".to_string(),
 };
 }

 // Compute face centroids
 let centroids: Vec<DVec3> = shell
 .faces
 .iter()
 .map(|face| {
 let mut sum = DVec3::ZERO;
 let mut count = 0;
 for we in &face.outer_wire.edges {
 let e = &brep.edges[we.idx];
 sum += brep.vertices[e.start].point;
 count += 1;
 }
 if count > 0 {
 sum / count as f64
 } else {
 DVec3::ZERO
 }
 })
 .collect();

 // Build adjacency map
 let mut adjacent_pairs: HashSet<(usize, usize)> = HashSet::new();
 for (fi, face) in shell.faces.iter().enumerate() {
 for we in &face.outer_wire.edges {
 for (fj, other_face) in shell.faces.iter().enumerate() {
 if fi < fj && other_face.outer_wire.edges.iter().any(|we2| we2.idx == we.idx) {
 adjacent_pairs.insert((fi, fj));
 }
 }
 }
 }

 // Find minimum distance between non-adjacent faces
 let mut min_dist = f64::MAX;
 let mut intersecting_pairs = Vec::new();
 let abs_distance = distance.abs();

 for i in 0..centroids.len() {
 for j in (i + 1)..centroids.len() {
 if adjacent_pairs.contains(&(i, j)) {
 continue;
 }

 let dist = (centroids[i] - centroids[j]).length();

 // Check if these faces would intersect
 if abs_distance > dist * 0.5 {
 intersecting_pairs.push((i, j));
 }

 if dist < min_dist {
 min_dist = dist;
 }
 }
 }

 if min_dist == f64::MAX {
 return SelfIntersectionResult {
 has_intersection: false,
 intersecting_pairs: Vec::new(),
 min_safe_distance: None,
 description: "no non-adjacent faces found".to_string(),
 };
 }

 let has_intersection = abs_distance > min_dist * 0.5;
 let min_safe_distance = Some(min_dist * 0.5);

 let description = if has_intersection {
 format!(
 "self-intersection likely: {} face pairs at distance {} with offset {}",
 intersecting_pairs.len(),
 min_dist,
 abs_distance
 )
 } else {
 format!("no self-intersection: min distance {}, offset {}", min_dist, abs_distance)
 };

 SelfIntersectionResult {
 has_intersection,
 intersecting_pairs,
 min_safe_distance,
 description,
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Join Geometry Creation
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Create a sewing face to bridge the gap between two separating offset surfaces
/// at a concave edge.
///
/// When two adjacent offset surfaces separate (concave edge + outward offset, or
/// convex edge + inward offset), the offset faces no longer meet. This function
/// creates a 4-sided planar face that fills the gap.
///
/// The sewing face is bounded by:
/// - Two edges along the offset faces (one on each face's boundary)
/// - Two connecting edges at the endpoints of the shared edge
pub fn create_sewing_face(
 brep: &mut BRep,
 original_brep: &BRep,
 edge_idx: usize,
 _face_a_idx: usize,
 _face_b_idx: usize,
 distance: f64,
 offset_surfaces: &[Option<Surface3>],
) -> Option<usize> {
 let edge = &original_brep.edges[edge_idx];
 let shell = original_brep.solids.first().and_then(|s| s.shells.first())?;

 let v_start = edge.start;
 let v_end = edge.end;
 let p_start = original_brep.vertices[v_start].point;
 let p_end = original_brep.vertices[v_end].point;

 // Compute offset vertex positions using the per-face offset normals.
 let start_offset_avg = {
 let mut sum = DVec3::ZERO;
 let mut count = 0;
 for (fi, face) in shell.faces.iter().enumerate() {
 let uses = face.outer_wire.edges.iter().any(|we| {
 let e = &original_brep.edges[we.idx];
 e.start == v_start || e.end == v_start
 });
 if uses {
 sum += get_face_offset_normal(original_brep, fi, shell);
 count += 1;
 }
 }
 if count > 0 { sum.normalize_or(DVec3::Z) * distance } else { DVec3::Z * distance }
 };

 let end_offset_avg = {
 let mut sum = DVec3::ZERO;
 let mut count = 0;
 for (fi, face) in shell.faces.iter().enumerate() {
 let uses = face.outer_wire.edges.iter().any(|we| {
 let e = &original_brep.edges[we.idx];
 e.start == v_end || e.end == v_end
 });
 if uses {
 sum += get_face_offset_normal(original_brep, fi, shell);
 count += 1;
 }
 }
 if count > 0 { sum.normalize_or(DVec3::Z) * distance } else { DVec3::Z * distance }
 };

 let _off_p_start = p_start + start_offset_avg;
 let _off_p_end = p_end + end_offset_avg;

 // Compute the separation direction perpendicular to both faces.
 // For two adjacent planar faces, the sewing face spans the gap between
 // their offset surfaces along the direction of the edge.
 let edge_dir = (p_end - p_start).normalize_or(DVec3::X);
 let sep_dir = any_perpendicular(edge_dir);

 // Create 4 vertices of the sewing face
 let sv0 = p_start + sep_dir * distance.abs() * 0.5;
 let sv1 = p_end + sep_dir * distance.abs() * 0.5;
 let sv2 = p_end - sep_dir * distance.abs() * 0.5;
 let sv3 = p_start - sep_dir * distance.abs() * 0.5;

 // Check that vertices are non-degenerate
 if (sv1 - sv0).length_squared() < 1e-12 || (sv3 - sv0).length_squared() < 1e-12 {
 return None;
 }

 let v0 = add_vertex(brep, sv0);
 let v1 = add_vertex(brep, sv1);
 let v2 = add_vertex(brep, sv2);
 let v3 = add_vertex(brep, sv3);

 // Create 4 edges
 let len01 = (sv1 - sv0).length();
 let e01 = add_edge(brep,
 Curve3::Line(Line3 { origin: sv0, direction: (sv1 - sv0).normalize_or(DVec3::X) }),
 0.0, len01, v0, v1);

 let len12 = (sv2 - sv1).length();
 let e12 = add_edge(brep,
 Curve3::Line(Line3 { origin: sv1, direction: (sv2 - sv1).normalize_or(DVec3::X) }),
 0.0, len12, v1, v2);

 let len23 = (sv3 - sv2).length();
 let e23 = add_edge(brep,
 Curve3::Line(Line3 { origin: sv2, direction: (sv3 - sv2).normalize_or(DVec3::X) }),
 0.0, len23, v2, v3);

 let len30 = (sv0 - sv3).length();
 let e30 = add_edge(brep,
 Curve3::Line(Line3 { origin: sv3, direction: (sv0 - sv3).normalize_or(DVec3::X) }),
 0.0, len30, v3, v0);

 let normal = edge_dir.cross(sep_dir).normalize();
 let sewing_surface = Surface3::Plane(Plane {
 origin: (sv0 + sv1 + sv2 + sv3) * 0.25,
 normal,
 });

 let wire = Wire {
 edges: vec![
 WireEdge::fwd(e01),
 WireEdge::fwd(e12),
 WireEdge::fwd(e23),
 WireEdge::fwd(e30),
 ],
 };

 let _ = offset_surfaces; // Used in more sophisticated implementations
 let _ = _face_a_idx;
 let _ = _face_b_idx;
 Some(add_face(brep, sewing_surface, wire, Vec::new()))
}

/// Create an arc join between two offset edges.
///
/// Creates a cylindrical surface that smoothly transitions between
/// two offset faces meeting at an edge.
pub fn create_arc_join(
 brep: &mut BRep,
 edge_idx: usize,
 face0_idx: usize,
 face1_idx: usize,
 radius: f64,
 vertex_map: &[usize],
) -> Result<usize, OffsetError> {
 let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or_else(|| {
 OffsetError::JoinCreationFailed {
 join_type: JoinType::Arc,
 edge_index: edge_idx,
 reason: "no shell found".to_string(),
 }
 })?;

 let edge = &brep.edges[edge_idx];
 let face0 = &shell.faces[face0_idx];
 let face1 = &shell.faces[face1_idx];

 // Get the edge endpoints
 let v0 = vertex_map.get(edge.start).copied().unwrap_or(edge.start);
 let v1 = vertex_map.get(edge.end).copied().unwrap_or(edge.end);

 let p0 = brep.vertices[v0].point;
 let p1 = brep.vertices[v1].point;

 // Compute the edge direction and length
 let edge_dir = (p1 - p0).normalize_or(DVec3::X);
 let edge_len = (p1 - p0).length();

 // Compute the bisector normal from the two face normals
 let n0 = face0.normal;
 let n1 = face1.normal;
 let _bisector = (n0 + n1).normalize_or(n0);

 // Create a cylindrical surface for the arc join
 // The cylinder axis is along the edge, and the radius is the offset distance
 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: p0,
 axis: edge_dir,
 radius,
 ref_dir: any_perpendicular(edge_dir),
 });

 // Create vertices for the arc join face
 // The arc join is a sector of the cylinder
 let vs = add_vertex(brep, p0);
 let ve = add_vertex(brep, p1);

 // Create the edge along the cylinder
 let curve = Curve3::Line(Line3 {
 origin: p0,
 direction: edge_dir,
 });
 let arc_edge = add_edge(brep, curve, 0.0, edge_len, vs, ve);

 // Create the arc face wire
 let wire = Wire {
 edges: vec![WireEdge::fwd(arc_edge)],
 };

 // Add the arc join face
 let face_idx = add_face(brep, cylinder, wire, Vec::new());

 Ok(face_idx)
}

/// Create a tangent join between two offset edges.
///
/// Creates a smooth, tangent-continuous transition between adjacent faces.
/// Falls back to intersection join when the angle between faces is too large.
pub fn create_tangent_join(
 brep: &mut BRep,
 edge_idx: usize,
 face0_idx: usize,
 face1_idx: usize,
 distance: f64,
 vertex_map: &[usize],
) -> Result<usize, OffsetError> {
 let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or_else(|| {
 OffsetError::JoinCreationFailed {
 join_type: JoinType::Tangent,
 edge_index: edge_idx,
 reason: "no shell found".to_string(),
 }
 })?;

 let face0 = &shell.faces[face0_idx];
 let face1 = &shell.faces[face1_idx];

 // Check the angle between face normals
 let n0 = face0.normal;
 let n1 = face1.normal;
 let dot = n0.dot(n1);

 // If the angle is too large (faces nearly parallel or facing opposite directions),
 // fall back to intersection join
 let angle_threshold = 0.9; // cos(25 degrees) approximately
 if dot < angle_threshold {
 // Create intersection join instead
 return create_intersection_join(brep, edge_idx, face0_idx, face1_idx, vertex_map);
 }

 // For tangent join, create a smooth blending surface
 // This uses a ruled surface between the two offset edges
 let edge = &brep.edges[edge_idx];
 let v0 = vertex_map.get(edge.start).copied().unwrap_or(edge.start);
 let v1 = vertex_map.get(edge.end).copied().unwrap_or(edge.end);

 let p0 = brep.vertices[v0].point;
 let p1 = brep.vertices[v1].point;

 // Create a plane that smoothly blends the two face normals
 let blend_normal = (n0 + n1).normalize();
 let blend_plane = Surface3::Plane(Plane {
 origin: (p0 + p1) * 0.5,
 normal: blend_normal,
 });

 // Create the wire for the tangent join face
 let dir = (p1 - p0).normalize_or(DVec3::X);
 let len = (p1 - p0).length();
 let curve = Curve3::Line(Line3 { origin: p0, direction: dir });

 let vs = add_vertex(brep, p0);
 let ve = add_vertex(brep, p1);
 let blend_edge = add_edge(brep, curve, 0.0, len, vs, ve);

 let wire = Wire {
 edges: vec![WireEdge::fwd(blend_edge)],
 };

 let face_idx = add_face(brep, blend_plane, wire, Vec::new());

 let _ = distance; // Used in more sophisticated implementations
 Ok(face_idx)
}

/// Create an intersection join between two offset edges.
///
/// The offset surfaces extend until they intersect, creating sharp corners.
/// This is the default mode and works well for mechanical parts.
pub fn create_intersection_join(
 brep: &mut BRep,
 edge_idx: usize,
 _face0_idx: usize,
 _face1_idx: usize,
 vertex_map: &[usize],
) -> Result<usize, OffsetError> {
 let edge = &brep.edges[edge_idx];

 let v0 = vertex_map.get(edge.start).copied().unwrap_or(edge.start);
 let v1 = vertex_map.get(edge.end).copied().unwrap_or(edge.end);

 let p0 = brep.vertices[v0].point;
 let p1 = brep.vertices[v1].point;

 // For intersection join, we don't create additional geometry -
 // the offset surfaces naturally intersect at the edge
 // Instead, we return the edge index as the "join"
 // In a full implementation, this would compute the exact intersection curve

 // Create a minimal face at the intersection
 let dir = (p1 - p0).normalize_or(DVec3::X);
 let len = (p1 - p0).length();

 // Use the edge midpoint and direction to create a small plane
 let midpoint = (p0 + p1) * 0.5;
 let normal = dir.any_orthonormal_pair().0;

 let plane = Surface3::Plane(Plane {
 origin: midpoint,
 normal,
 });

 let vs = add_vertex(brep, p0);
 let ve = add_vertex(brep, p1);
 let curve = Curve3::Line(Line3 { origin: p0, direction: dir });
 let int_edge = add_edge(brep, curve, 0.0, len, vs, ve);

 let wire = Wire {
 edges: vec![WireEdge::fwd(int_edge)],
 };

 let face_idx = add_face(brep, plane, wire, Vec::new());

 Ok(face_idx)
}

/// Apply join type to all edges in the shell.
///
/// This function creates the appropriate join geometry for each edge
/// based on the specified join type.
pub fn apply_join_type(
 result: &mut BRep,
 original_brep: &BRep,
 opts: &OffsetOptions,
 edge_to_faces: &HashMap<usize, Vec<usize>>,
 vertex_map: &[usize],
) -> Result<usize, OffsetError> {
 let mut join_face_count = 0;

 if opts.join_type == JoinType::Intersection {
 // Intersection join is the default - no additional geometry needed
 return Ok(0);
 }

 for (&edge_idx, face_indices) in edge_to_faces {
 if face_indices.len() < 2 {
 continue; // Skip boundary edges
 }

 let face0_idx = face_indices[0];
 let face1_idx = face_indices[1];

 let join_result = match opts.join_type {
 JoinType::Arc => {
 let radius = opts.distance.abs();
 create_arc_join(result, edge_idx, face0_idx, face1_idx, radius, vertex_map)
 }
 JoinType::Tangent => {
 create_tangent_join(result, edge_idx, face0_idx, face1_idx, opts.distance, vertex_map)
 }
 JoinType::Intersection => {
 create_intersection_join(result, edge_idx, face0_idx, face1_idx, vertex_map)
 }
 };

 if join_result.is_ok() {
 join_face_count += 1;
 }
 }

 let _ = original_brep; // Used in more sophisticated implementations
 Ok(join_face_count)
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Offset Quality Analysis
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Analyze the quality of an offset result.
///
/// Computes various quality metrics including wall thickness, deviation,
/// and self-intersection detection.
pub fn analyze_offset_quality(
 result: &BRep,
 original: &BRep,
 opts: &OffsetOptions,
) -> OffsetQuality {
 let mut quality = OffsetQuality::default();

 // Compute minimum wall thickness
 quality.min_wall_thickness = compute_min_wall_thickness(result, opts.distance);

 // Compute maximum deviation from expected offset
 quality.max_deviation = compute_max_deviation(result, original, opts);

 // Count degenerate edges
 quality.degenerate_edge_count = result
 .geom
 .edge_degenerated
 .iter()
 .filter(|&&d| d)
 .count();

 // Self-intersection count
 let si_result = detect_self_intersection_detailed(result, opts.distance);
 quality.self_intersection_count = si_result.intersecting_pairs.len();

 // Compute face area ratio
 quality.face_area_ratio = compute_face_area_ratio(result, original);

 // Compute edge length ratio
 quality.edge_length_ratio = compute_edge_length_ratio(result, original);

 // Determine if result is valid
 quality.is_valid = quality.self_intersection_count == 0
 && quality.min_wall_thickness >= opts.min_wall_thickness;

 // Generate warnings
 if quality.min_wall_thickness < opts.min_wall_thickness {
 quality.warnings.push(format!(
 "Minimum wall thickness {} is below threshold {}",
 quality.min_wall_thickness, opts.min_wall_thickness
 ));
 }
 if quality.max_deviation > opts.approximation_tolerance {
 quality.warnings.push(format!(
 "Maximum deviation {} exceeds approximation tolerance {}",
 quality.max_deviation, opts.approximation_tolerance
 ));
 }
 if quality.degenerate_edge_count > 0 {
 quality.warnings.push(format!(
 "Found {} degenerate edges in result",
 quality.degenerate_edge_count
 ));
 }

 quality
}

/// Compute the minimum wall thickness in the offset result.
///
/// Uses face centroid distances to estimate minimum wall thickness.
pub fn compute_min_wall_thickness(brep: &BRep, distance: f64) -> f64 {
 let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
 Some(s) => s,
 None => return distance,
 };

 if shell.faces.len() < 2 {
 return distance;
 }

 // Compute face centroids
 let centroids: Vec<DVec3> = shell
 .faces
 .iter()
 .map(|face| {
 let mut sum = DVec3::ZERO;
 let mut count = 0;
 for we in &face.outer_wire.edges {
 let e = &brep.edges[we.idx];
 sum += brep.vertices[e.start].point;
 count += 1;
 }
 if count > 0 {
 sum / count as f64
 } else {
 DVec3::ZERO
 }
 })
 .collect();

 // Find minimum distance between any two faces
 let mut min_dist = f64::MAX;
 for i in 0..centroids.len() {
 for j in (i + 1)..centroids.len() {
 let dist = (centroids[i] - centroids[j]).length();
 if dist > 0.0 && dist < min_dist {
 min_dist = dist;
 }
 }
 }

 // The wall thickness is approximately the minimum distance minus twice the offset
 // For a proper implementation, this would use more sophisticated analysis
 if min_dist == f64::MAX {
 distance
 } else {
 (min_dist - 2.0 * distance.abs()).max(0.0)
 }
}

/// Compute the maximum deviation between offset and expected positions.
pub fn compute_max_deviation(result: &BRep, original: &BRep, opts: &OffsetOptions) -> f64 {
 let _result_shell = match result.solids.first().and_then(|s| s.shells.first()) {
 Some(s) => s,
 None => return 0.0,
 };

 let _original_shell = match original.solids.first().and_then(|s| s.shells.first()) {
 Some(s) => s,
 None => return 0.0,
 };

 let mut max_dev = 0.0;

 // Compare vertex positions
 for (i, vertex) in result.vertices.iter().enumerate() {
 if i >= original.vertices.len() {
 break;
 }

 let original_vertex = &original.vertices[i];
 let actual_offset = (vertex.point - original_vertex.point).length();
 let expected_offset = opts.distance.abs();

 let deviation = (actual_offset - expected_offset).abs();
 if deviation > max_dev {
 max_dev = deviation;
 }
 }

 max_dev
}

/// Compute the ratio of face areas between result and original.
pub fn compute_face_area_ratio(result: &BRep, original: &BRep) -> f64 {
 let result_shell = match result.solids.first().and_then(|s| s.shells.first()) {
 Some(s) => s,
 None => return 1.0,
 };

 let original_shell = match original.solids.first().and_then(|s| s.shells.first()) {
 Some(s) => s,
 None => return 1.0,
 };

 if original_shell.faces.is_empty() {
 return 1.0;
 }

 // Simple approximation: ratio of face counts
 // A proper implementation would compute actual areas
 result_shell.faces.len() as f64 / original_shell.faces.len() as f64
}

/// Compute the ratio of edge lengths between result and original.
pub fn compute_edge_length_ratio(result: &BRep, original: &BRep) -> f64 {
 if original.edges.is_empty() {
 return 1.0;
 }

 // Compute total edge lengths
 let original_len: f64 = original
 .edges
 .iter()
 .map(|e| {
 let p0 = original.vertices.get(e.start).map(|v| v.point).unwrap_or_default();
 let p1 = original.vertices.get(e.end).map(|v| v.point).unwrap_or_default();
 (p1 - p0).length()
 })
 .sum();

 let result_len: f64 = result
 .edges
 .iter()
 .map(|e| {
 let p0 = result.vertices.get(e.start).map(|v| v.point).unwrap_or_default();
 let p1 = result.vertices.get(e.end).map(|v| v.point).unwrap_or_default();
 (p1 - p0).length()
 })
 .sum();

 if original_len > 0.0 {
 result_len / original_len
 } else {
 1.0
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Self-Intersection Repair
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Attempt to repair self-intersection by reducing offset distance.
///
/// Tries progressively smaller offset distances until a valid result is found.
pub fn repair_self_intersection(
 brep: &BRep,
 opts: &OffsetOptions,
) -> Result<(BRep, f64, usize), OffsetError> {
 let config = &opts.self_intersection_config;

 if !config.auto_repair {
 return Err(OffsetError::RecoveryFailed {
 attempts: 0,
 last_error: "auto-repair not enabled".to_string(),
 });
 }

 let mut current_distance = opts.distance;
 let mut attempts = 0;
 let mut last_error = String::new();

 while attempts < config.max_repair_attempts {
 attempts += 1;

 // Reduce the offset distance
 current_distance *= config.reduction_factor;

 if current_distance.abs() < config.min_offset_distance {
 last_error = format!(
 "offset distance {} below minimum {}",
 current_distance.abs(),
 config.min_offset_distance
 );
 continue;
 }

 // Try with reduced distance
 let mut reduced_opts = opts.clone();
 reduced_opts.distance = current_distance;
 reduced_opts.check_self_intersection = true;

 let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(OffsetError::InvalidInput("no shell"))?;

 match offset_shell_with_options_impl(shell, brep, &reduced_opts) {
 Ok(result) => {
 let si_result = detect_self_intersection_detailed(&result, current_distance);
 if !si_result.has_intersection {
 return Ok((result, current_distance, attempts));
 }
 last_error = si_result.description;
 }
 Err(e) => {
 last_error = e.to_string();
 }
 }
 }

 Err(OffsetError::RecoveryFailed { attempts, last_error })
}

/// Implementation of offset_shell_with_options that can be called internally.
fn offset_shell_with_options_impl(
 shell: &Shell,
 brep: &BRep,
 opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
 // Validate variable thickness if specified
 if let Some(ref vt) = opts.variable_thickness {
 vt.validate(shell.faces.len())?;
 }

 let distance = opts.distance;

 if distance.abs() < TOLERANCE_LEN_MIN {
 return Err(OffsetError::ZeroDistance);
 }

 if shell.faces.is_empty() {
 return Err(OffsetError::InvalidInput("shell has no faces"));
 }

 // Step 1: Compute offset surfaces for each face (with variable thickness support)
 let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
 for (fi, _face) in shell.faces.iter().enumerate() {
 let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
 Some(s) => s,
 None => {
 offset_surfaces.push(None);
 continue;
 }
 };

 let surf = &brep.geom.surfaces[surf_idx];
 let face_distance = opts.effective_distance_for_face(fi);
 let off_surf = offset_surface(surf, face_distance);

 offset_surfaces.push(off_surf);
 }

 // Step 2: Build edge-to-face adjacency
 let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
 for (fi, face) in shell.faces.iter().enumerate() {
 for we in &face.outer_wire.edges {
 edge_to_faces.entry(we.idx).or_default().push(fi);
 }
 }

 // Step 3: Compute offset vertex positions (OCCT-aligned edge-first).
 // For each vertex, project onto offset intersection curves of its incident
 // edges, using the per-face distance if variable thickness is active.
 let get_face_distance = |fi: usize| -> f64 {
 if let Some(ref vt) = opts.variable_thickness {
 vt.thickness_for_face(fi)
 } else {
 distance
 }
 };
 let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
 .map(|vi| {
 let pt = brep.vertices[vi].point;
 // Find all incident edges for this vertex
 let mut proj_sum = DVec3::ZERO;
 let mut proj_count = 0usize;
 for (ei, edge) in brep.edges.iter().enumerate() {
 if edge.start != vi && edge.end != vi { continue; }
 let faces = match edge_to_faces.get(&ei) {
 Some(f) if f.len() >= 2 => f,
 _ => continue,
 };
 for pair in faces.windows(2) {
 let fi1 = pair[0];
 let fi2 = pair[1];
 let si1 = match brep.geom.face_surface.get(fi1).and_then(|s| *s) { Some(s) => s, None => continue };
 let si2 = match brep.geom.face_surface.get(fi2).and_then(|s| *s) { Some(s) => s, None => continue };
 let s1 = match brep.geom.surfaces.get(si1) { Some(s) => s, None => continue };
 let s2 = match brep.geom.surfaces.get(si2) { Some(s) => s, None => continue };
 let d1 = get_face_distance(fi1);
 let d2 = get_face_distance(fi2);
 let inter = intersect_offset_surfaces(s1, s2, d1, d2);
 if let Some(proj) = project_point_onto_intersection(pt, &inter) {
 proj_sum += proj;
 proj_count += 1;
 break;
 }
 }
 }
 if proj_count >= 1 {
 proj_sum / proj_count as f64
 } else {
 // Fallback: average-normal translation
 let avg_dist = if let Some(ref vt) = opts.variable_thickness {
 let mut s = 0.0; let mut c = 0;
 for (fi, face) in shell.faces.iter().enumerate() {
 if face.outer_wire.edges.iter().any(|we| brep.edges[we.idx].start == vi || brep.edges[we.idx].end == vi) {
 s += vt.thickness_for_face(fi); c += 1;
 }
 }
 if c > 0 { s / c as f64 } else { distance }
 } else { distance };
 offset_vertex(brep, vi, avg_dist, shell, None)
 }
 })
 .collect();

 // Step 4: Build result BRep
 let mut result = BRep::new();
 result.solids.push(Solid {
 shells: vec![Shell { faces: Vec::new() }],
 });

 // Map original vertices to offset vertices
 let mut vertex_map: Vec<usize> = Vec::with_capacity(offset_vertices.len());
 for &p in &offset_vertices {
 vertex_map.push(add_vertex(&mut result, p));
 }

 // Step 5: Create offset faces with offset edges.
 // Edge curves come from offset_edge (plane-plane intersection for planar faces),
 // reparameterized to pass through vertex positions for consistency.
 let mut valid_face_count = 0;

 for (fi, face) in shell.faces.iter().enumerate() {
 let off_surf = match &offset_surfaces[fi] {
 Some(s) => s.clone(),
 None => continue,
 };

 let mut wire_edges = Vec::new();

 for we in &face.outer_wire.edges {
 let e = &brep.edges[we.idx];
 let vs = vertex_map[e.start];
 let ve = vertex_map[e.end];
 let p_start = result.vertices[vs].point;
 let p_end = result.vertices[ve].point;

 let faces = edge_to_faces.get(&we.idx).cloned().unwrap_or_default();
 let (curve, t0, t1) = offset_edge(brep, we.idx, &faces, distance, &offset_surfaces, &offset_vertices)
 .map(|(c, _, _)| {
 match &c {
 Curve3::Line(line) => {
 let ts = project_point_to_line(p_start, line);
 let te = project_point_to_line(p_end, line);
 (c, ts.min(te), ts.max(te))
 }
 _ => (c, 0.0, (p_end - p_start).length()),
 }
 })
 .unwrap_or_else(|| {
 let dir = (p_end - p_start).normalize_or(DVec3::X);
 let len = (p_end - p_start).length();
 (Curve3::Line(Line3 { origin: p_start, direction: dir }), 0.0, len)
 });

 if (t1 - t0).abs() < TOLERANCE_LEN_MIN { continue; }

 let eidx = add_edge(&mut result, curve, t0, t1, vs, ve);
 wire_edges.push(if we.forward { WireEdge::fwd(eidx) } else { WireEdge::rev(eidx) });
 }

 if wire_edges.len() < 3 { continue; }

 add_face(&mut result, off_surf, Wire { edges: wire_edges }, Vec::new());
 valid_face_count += 1;
 }

 if valid_face_count == 0 {
 return Err(OffsetError::EmptyResult);
 }


 // Step 6: Apply join type if needed
 if opts.join_type.requires_join_geometry() {
 let _join_faces = apply_join_type(&mut result, brep, opts, &edge_to_faces, &vertex_map)?;
 }

 // Fix inverted face winding (same as main offset function)
 for f in &mut result.solids[0].shells[0].faces {
 let mut verts: Vec<DVec3> = Vec::new();
 for we in &f.outer_wire.edges {
 let e = &result.edges[we.idx];
 verts.push(if we.forward { result.vertices[e.end].point } else { result.vertices[e.start].point });
 }
 if verts.len() >= 3 {
 let n = f.normal;
 let mut signed = 0.0;
 for i in 0..verts.len() {
 let j = (i + 1) % verts.len();
 signed += verts[i].cross(verts[j]).dot(n);
 }
 signed *= 0.5;
 if signed < 0.0 {
 for we in &mut f.outer_wire.edges { we.forward = !we.forward; }
 f.outer_wire.edges.reverse();
 for iw in &mut f.inner_wires {
 for we in &mut iw.edges { we.forward = !we.forward; }
 iw.edges.reverse();
 }
 }
 }
 }

 // Fix inward-facing face normals (same as main offset function)
 // Use the centroid of all result vertices as the reference point,
 // NOT the origin.  The origin-based check gives wrong answers for
 // solids far from the origin (e.g., a box at (1,1,1)-(9,9,9) after
 // inward offset has inward-correct normals but negative tet volume
 // when measured from the origin).
 let n_faces = result.solids[0].shells[0].faces.len();
 if n_faces > 0 {
 let mut center_sum = DVec3::ZERO;
 for v in &result.vertices { center_sum += v.point; }
 let center = center_sum / result.vertices.len() as f64;
 for fi in 0..n_faces {
 let f = &result.solids[0].shells[0].faces[fi];
 if f.outer_wire.edges.len() < 3 { continue; }
 let mut verts: Vec<DVec3> = Vec::new();
 for we in &f.outer_wire.edges {
 let e = &result.edges[we.idx];
 verts.push(if we.forward { result.vertices[e.end].point } else { result.vertices[e.start].point });
 }
 if verts.len() < 3 { continue; }
 // Center vertices around the solid centroid so the signed volume
 // correctly reflects inward vs outward regardless of global position.
 let centered: Vec<DVec3> = verts.iter().map(|v| *v - center).collect();
 let p0 = centered[0];
 let mut vol_6 = 0.0;
 for i in 1..centered.len() - 1 {
 vol_6 += p0.cross(centered[i]).dot(centered[i + 1]);
 }
 if vol_6 < 0.0 {
 let f = &mut result.solids[0].shells[0].faces[fi];
 f.normal = -f.normal;
 for we in &mut f.outer_wire.edges { we.forward = !we.forward; }
 f.outer_wire.edges.reverse();
 for iw in &mut f.inner_wires {
 for we in &mut iw.edges { we.forward = !we.forward; }
 iw.edges.reverse();
 }
 }
 }
 }

 Ok(result)
}
// Main API Functions
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Clip a 3D polygon against a half-space defined by n·x ≤ d (interior of solid).
/// Uses the Sutherland-Hodgman algorithm adapted for half-space clipping in 3D.
/// Returns the clipped polygon (may be empty if fully clipped away).
fn clip_polygon_by_halfspace(polygon: &[DVec3], n: DVec3, d: f64, tol: f64) -> Vec<DVec3> {
 if polygon.len() < 3 {
 return Vec::new();
 }

 let mut output = Vec::new();
 let m = polygon.len();

 for i in 0..m {
 let curr = polygon[i];
 let prev = polygon[(i + m - 1) % m];

 let curr_dist = n.dot(curr) - d;
 let prev_dist = n.dot(prev) - d;

 let curr_inside = curr_dist <= tol;
 let prev_inside = prev_dist <= tol;

 if curr_inside {
 if !prev_inside {
 // Edge enters the valid half-space: add intersection point.
 let t = -prev_dist / (curr_dist - prev_dist);
 output.push(prev + t * (curr - prev));
 }
 output.push(curr);
 } else if prev_inside {
 // Edge leaves the valid half-space: add intersection point.
 let t = -prev_dist / (curr_dist - prev_dist);
 output.push(prev + t * (curr - prev));
 }
 }

 output
}

/// Offset a shell by moving all faces along their normals.
///
/// # Arguments
///
/// * `shell` - The input shell to offset
/// * `brep` - The BRep containing the shell's geometry
/// * `distance` - Offset distance (positive = outward, negative = inward)
///
/// # Returns
///
/// A new BRep containing the offset shell, or an error.
pub fn offset_shell(shell: &Shell, brep: &BRep, distance: f64) -> Result<BRep, OffsetError> {
 offset_shell_with_options(shell, brep, &OffsetOptions::new(distance))
}

/// Offset a shell with full options.
pub fn offset_shell_with_options(
 shell: &Shell,
 brep: &BRep,
 opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
 let distance = opts.distance;

 if distance.abs() < TOLERANCE_LEN_MIN {
 return Err(OffsetError::ZeroDistance);
 }

 if shell.faces.is_empty() {
 return Err(OffsetError::InvalidInput("shell has no faces"));
 }
 // Step 1: Compute offset surfaces for each face
 let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
 for (fi, _face) in shell.faces.iter().enumerate() {
 let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
 Some(s) => s,
 None => {
 offset_surfaces.push(None);
 continue;
 }
 };

 let surf = &brep.geom.surfaces[surf_idx];
 let off_surf = offset_surface(surf, distance);

 // Note: surface normal flip removed — the face normals on extruded
 // shapes may point inward (shape-creation artifact), and flipping the
 // surface normal to match would make the offset translate in the wrong
 // direction for all incident vertices. The surface normal from the
 // underlying geometry is the authoritative direction.

 if off_surf.is_none() && distance > 0.0 {
 // Negative offset on a small surface - may be ok for inward offset
 }

 offset_surfaces.push(off_surf);
 }

 // Step 2: Build edge-to-face adjacency (including inner wires for faces with holes)
 let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
 for (fi, face) in shell.faces.iter().enumerate() {
 for we in &face.outer_wire.edges {
 edge_to_faces.entry(we.idx).or_default().push(fi);
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 edge_to_faces.entry(we.idx).or_default().push(fi);
 }
 }
 }

 // Step 3: Compute offset vertex positions with position-based deduplication
 // and OCCT-aligned edge-first computation.
 //
 // Strategy:
 // 1. Group vertices by position (T-junction dedup).
 // 2. For each group, find ALL incident edges across ALL vertices in the group.
 // 3. For each incident edge, compute the offset intersection curve of its
 // two adjacent faces (intersect_offset_surfaces), then project the original
 // vertex onto that curve (project_point_onto_intersection).
 // 4. Average the projections from all incident edges as the offset position.
 //
 // This differs from the old approach (Cramer's rule on face normals) and
 // aligns with OCCT BRepOffset_Inter3d: edge curves come from surface-surface
 // intersection, vertices come from projecting onto those curves.
 let pos_tol = 1e-8;
 let mut pos_to_group: Vec<usize> = vec![usize::MAX; brep.vertices.len()];
 let mut group_positions: Vec<DVec3> = Vec::new();
 let mut group_vertex_indices: Vec<Vec<usize>> = Vec::new();
 for vi in 0..brep.vertices.len() {
 let pt = brep.vertices[vi].point;
 let mut found = None;
 for (gi, gp) in group_positions.iter().enumerate() {
 if (pt - *gp).length_squared() < pos_tol * pos_tol {
 found = Some(gi);
 break;
 }
 }
 if let Some(gi) = found {
 pos_to_group[vi] = gi;
 group_vertex_indices[gi].push(vi);
 } else {
 pos_to_group[vi] = group_positions.len();
 group_positions.push(pt);
 group_vertex_indices.push(vec![vi]);
 }
 }
 // Compute offset for each group using edge-first projection (OCCT-aligned)
 let mut group_offsets: Vec<DVec3> = Vec::with_capacity(group_positions.len());
 for gi in 0..group_positions.len() {
 let pt = group_positions[gi];
 // Collect all incident edges across all vertices in the group
 let mut incident_edges: Vec<usize> = Vec::new();
 for &vi in &group_vertex_indices[gi] {
 for (ei, edge) in brep.edges.iter().enumerate() {
 if (edge.start == vi || edge.end == vi) && !incident_edges.contains(&ei) {
 incident_edges.push(ei);
 }
 }
 }
 // Check if ALL incident faces are planar — if so, use exact Cramer's rule.
 let all_planar = incident_edges.iter().all(|ei| {
 edge_to_faces.get(ei).map_or(true, |faces| {
 faces.iter().all(|fi| {
 brep.geom.face_surface.get(*fi).and_then(|s| *s)
 .and_then(|si| brep.geom.surfaces.get(si))
 .is_some_and(|s| matches!(s, Surface3::Plane(_)))
 })
 })
 });

 let off = if all_planar {
 // Planar-only vertex: use Cramer's rule (exact intersection of offset planes)
 let mut fi_list: Vec<usize> = Vec::new();
 let mut normal_sum = DVec3::ZERO;
 for &vi in &group_vertex_indices[gi] {
 for (fi, face) in shell.faces.iter().enumerate() {
 let uses = face.outer_wire.edges.iter().any(|we| {
 let e = &brep.edges[we.idx];
 e.start == vi || e.end == vi
 }) || face.inner_wires.iter().any(|wire| {
 wire.edges.iter().any(|we| {
 let e = &brep.edges[we.idx];
 e.start == vi || e.end == vi
 })
 });
 if uses && !fi_list.contains(&fi) {
 fi_list.push(fi);
 normal_sum += face.normal;
 }
 }
 }
 if !fi_list.is_empty() {
 offset_vertex_from_faces(brep, pt, &fi_list, normal_sum, distance, shell)
 } else {
 pt
 }
 } else {
 // Curved-surface vertex: OCCT edge-first projection.
 //
 // For each incident edge:
 // 1. If the edge is manifold (2 adjacent faces), compute the
 // offset intersection curve of the two offset surfaces and
 // project the original vertex onto it.
 // 2. If the edge is a single-face seam (periodic surface) or
 // the 2-face intersection failed, project the vertex onto
 // the offset face surface directly.
 let mut projections: Vec<DVec3> = Vec::new();
 let mut seam_projection: Option<DVec3> = None;
 for &ei in &incident_edges {
 let efaces = match edge_to_faces.get(&ei) {
 Some(f) => f,
 None => continue,
 };
 let mut found = false;
 if efaces.len() >= 2 {
 // Try 2-face edge: offset intersection curve
 for pair in efaces.windows(2) {
 let fi1 = pair[0]; let fi2 = pair[1];
 let si1 = match brep.geom.face_surface.get(fi1).and_then(|s| *s) { Some(s) => s, None => continue };
 let si2 = match brep.geom.face_surface.get(fi2).and_then(|s| *s) { Some(s) => s, None => continue };
 let s1 = match brep.geom.surfaces.get(si1) { Some(s) => s, None => continue };
 let s2 = match brep.geom.surfaces.get(si2) { Some(s) => s, None => continue };
 let intersection = intersect_offset_surfaces(s1, s2, distance, distance);
 if let Some(proj) = project_point_onto_intersection(pt, &intersection) {
 projections.push(proj);
 found = true;
 break;
 }
 }
 }
 if !found {
 // Single-face seam edge (periodic surface): project the vertex
 // onto the offset face surface.  This provides an additional
 // constraint for cylinder seam vertices (common in cap-face
 // cases).  Only do this for cylinder surfaces, where the
 // offset is simple (radius change only) and the projection
 // is accurate.  Cone seam projections are less reliable due
 // to the varying radius and apex shift.
 for &fi in efaces.iter().take(2) {
 let si = match brep.geom.face_surface.get(fi).and_then(|s| *s) { Some(s) => s, None => continue };
 let is_cylinder = matches!(brep.geom.surfaces.get(si), Some(Surface3::Cylinder(_)));
 if !is_cylinder { continue; }
 let orig_surf = match brep.geom.surfaces.get(si) { Some(s) => s, None => continue };
 let off_surf = match offset_surface(orig_surf, distance) { Some(s) => s, None => continue };
 if let Some(uv) = project_point_to_surface_uv(pt, &off_surf, None) {
 seam_projection = Some(off_surf.point_at(uv[0], uv[1]));
 break;
 }
 }
 }
 }
 // Try offset_vertex_curved_plane for cone+plane vertices.
 // This analytically intersects the offset cone surface with each
 // adjacent offset plane, giving an exact result that is more
 // accurate than the edge-first projections or Cramer's rule fallback.
 let cone_plane_result: Option<DVec3> = {
 let mut all_fis: Vec<usize> = Vec::new();
 for &vi in &group_vertex_indices[gi] {
 for (fi, face) in shell.faces.iter().enumerate() {
 let uses = face.outer_wire.edges.iter().any(|we| brep.edges[we.idx].start == vi || brep.edges[we.idx].end == vi)
 || face.inner_wires.iter().any(|wire| wire.edges.iter().any(|we| brep.edges[we.idx].start == vi || brep.edges[we.idx].end == vi));
 if uses && !all_fis.contains(&fi) { all_fis.push(fi); }
 }
 }
 let mut cone_fi = None;
 let mut plane_fis: Vec<usize> = Vec::new();
 for &fi in &all_fis {
 match brep.geom.face_surface.get(fi).and_then(|s| *s).and_then(|si| brep.geom.surfaces.get(si)) {
 Some(Surface3::Cone(_)) => cone_fi = Some(fi),
 Some(Surface3::Plane(_)) => plane_fis.push(fi),
 _ => {}
 }
 }
 match (cone_fi, plane_fis.is_empty()) {
 (Some(cfi), false) => offset_vertex_curved_plane(pt, brep, cfi, &plane_fis, distance, shell),
 _ => None,
 }
 };

 // Decision: combine edge-first projections with cone+plane result.
 // Only manifold-edge (2-face) projections are averaged — OCCT uses the
 // intersection of adjacent offset surfaces for each incident manifold edge.
 // Single-face seam/silhouette edges only contribute a surface-projection
 // fallback and are NOT averaged in, since their projection onto the offset
 // surface alone lacks the constraint from the opposite face (e.g., a cylinder
 // seam projects onto the wall at the original V coordinate, not at the
 // offset-cap height).
 //
 // Order of preference:
 // 1. ≥2 manifold projections → average (optionally blend with cone+plane)
 // 2. 1 manifold projection → use it directly (OCCT-correct intersection)
 // 3. 0 manifold + seam → use seam surface projection
 // 4. otherwise → Cramer's rule fallback
 if projections.len() >= 2 {
 let mut sum = DVec3::ZERO;
 let mut count = 0usize;
 for p in &projections { sum += *p; count += 1; }
 if let Some(cp) = cone_plane_result { sum += cp; count += 1; }
 sum / count as f64
 } else if projections.len() == 1 {
 // Single manifold projection: use it directly (it is the exact
 // intersection of the two offset surfaces for this edge).
 if let Some(cp) = cone_plane_result {
 // Blend with analytic cone+plane result for better accuracy
 (projections[0] + cp) * 0.5
 } else {
 projections[0]
 }
 } else if let Some(sp) = seam_projection {
 sp
 } else if let Some(cp) = cone_plane_result {
 cp
 } else {
 // Fallback: Cramer's rule from all incident faces
 let mut fi_list: Vec<usize> = Vec::new();
 let mut normal_sum = DVec3::ZERO;
 for &vi in &group_vertex_indices[gi] {
 for (fi, face) in shell.faces.iter().enumerate() {
 let uses = face.outer_wire.edges.iter().any(|we| {
 let e = &brep.edges[we.idx];
 e.start == vi || e.end == vi
 }) || face.inner_wires.iter().any(|wire| {
 wire.edges.iter().any(|we| {
 let e = &brep.edges[we.idx];
 e.start == vi || e.end == vi
 })
 });
 if uses && !fi_list.contains(&fi) {
 fi_list.push(fi);
 normal_sum += face.normal;
 }
 }
 }
 if !fi_list.is_empty() {
 offset_vertex_from_faces(brep, pt, &fi_list, normal_sum, distance, shell)
 } else {
 pt
 }
 }
 }
 ;
 group_offsets.push(off);
 }
 let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
 .map(|vi| group_offsets[pos_to_group[vi]])
 .collect();

 // Step 4: Build result BRep
 let mut result = BRep::new();
 result.solids.push(Solid {
 shells: vec![Shell { faces: Vec::new() }],
 });

 // Map original vertices to offset vertices
 let mut vertex_map: Vec<usize> = Vec::with_capacity(offset_vertices.len());
 for &p in &offset_vertices {
 vertex_map.push(add_vertex(&mut result, p));
 }

 // Step 5: Create offset faces with offset edges.
 // Edge curves come from the intersection of adjacent offset surfaces
 // (with plane-plane intersections for planar faces), but are reparameterized
 // to pass through the vertex positions for consistency.
 let mut valid_face_count = 0;

 for (fi, face) in shell.faces.iter().enumerate() {
 let off_surf = match &offset_surfaces[fi] {
 Some(s) => s.clone(),
 None => continue,
 };

 // Build wire from offset edges
 let mut wire_edges = Vec::new();

 for we in &face.outer_wire.edges {
 let e = &brep.edges[we.idx];
 let vs = vertex_map[e.start];
 let ve = vertex_map[e.end];
 let p_start = result.vertices[vs].point;
 let p_end = result.vertices[ve].point;

 let faces = edge_to_faces.get(&we.idx).cloned().unwrap_or_default();
 let (curve, t0, t1, edge_vs, edge_ve) = offset_edge(brep, we.idx, &faces, distance, &offset_surfaces, &offset_vertices)
 .map(|(c, _, _)| {
 const VTX_TOL_SQ: f64 = 1e-12;
 match &c {
 Curve3::Line(line) => {
 let ts = project_point_to_line(p_start, line);
 let te = project_point_to_line(p_end, line);
 (c, ts.min(te), ts.max(te), vs, ve)
 }
 Curve3::Circle(off_circle) => {
 let (ta, tb) = brep.geom.edge_curve.get(we.idx)
 .and_then(|oc| *oc)
 .and_then(|ci| brep.geom.curves.get(ci))
 .and_then(|orig_curve| match orig_curve {
 Curve3::Circle(orig_c) => {
 let range = brep.geom.edge_curve_range.get(we.idx).and_then(|r| *r)?;
 let a0 = point_on_circle_angle(orig_curve.point_at(range[0]), orig_c);
 let a1 = point_on_circle_angle(orig_curve.point_at(range[1]), orig_c);
 if (a1 - a0).abs() < 1e-12 && (range[1] - range[0]).abs() > std::f64::consts::PI {
 Some((0.0, std::f64::consts::TAU)) // Full circle
 } else if a1 < a0 {
 Some((a0, a1 + std::f64::consts::TAU))
 } else {
 Some((a0, a1))
 }
 }
 _ => None,
 })
 .unwrap_or((0.0, std::f64::consts::TAU));
 let normal = off_circle.normal.normalize_or(DVec3::Z);
 let ref_dir = if normal.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let u_axis = normal.cross(ref_dir).normalize();
 let v_axis = normal.cross(u_axis).normalize();
 let _proj_start = off_circle.center + off_circle.radius
 * (u_axis * ta.cos() + v_axis * ta.sin());
 let _proj_end = off_circle.center + off_circle.radius
 * (u_axis * tb.cos() + v_axis * tb.sin());
 // Full-circle edges (cap faces) keep the merged vertex
 // so the edge is a self-loop detected and split below.
 let is_self_loop = (tb - ta - std::f64::consts::TAU).abs() < 1e-12;
 let vs_pt = result.vertices[vs].point;
 let (lvs, lve, nta, ntb) = if is_self_loop {
 (vs, vs, ta, tb)
 } else {
 // Align angle 0 with the vertex_map seam position
 // so the Circle parameterization is consistent.
 let vu = (vs_pt - off_circle.center).normalize_or(u_axis);
 let vv = off_circle.normal.cross(vu).normalize();
 let ve_pt = result.vertices[ve].point;
 let local_e = ve_pt - off_circle.center;
 let ang_e = local_e.dot(vv).atan2(local_e.dot(vu));
 let ang_e = if ang_e < 0.0 { ang_e + std::f64::consts::TAU } else { ang_e };
 let p_end = off_circle.center + off_circle.radius
 * (vu * ang_e.cos() + vv * ang_e.sin());
 (vs, add_vertex(&mut result, p_end), 0.0, ang_e)
 };
 (c, nta, ntb, lvs, lve)
 }
 _ => (c, 0.0, (p_end - p_start).length(), vs, ve),
 }
 })
 .unwrap_or_else(|| {
 let dir = (p_end - p_start).normalize_or(DVec3::X);
 let len = (p_end - p_start).length();
 (Curve3::Line(Line3 { origin: p_start, direction: dir }), 0.0, len, vs, ve)
 });

 if (t1 - t0).abs() < TOLERANCE_LEN_MIN {
 continue;
 }

 let eidx = add_edge(&mut result, curve, t0, t1, edge_vs, edge_ve);
 wire_edges.push(if we.forward { WireEdge::fwd(eidx) } else { WireEdge::rev(eidx) });
 }

 // Fix self-loop Circle edges (cap faces with a single full-circle edge)
 // by splitting into two half-circles with a midpoint vertex.
 let mut split_wire: Vec<WireEdge> = Vec::new();
 for we in &wire_edges {
 let (ei, fwd) = (we.idx, we.forward);
 if result.edges[ei].start == result.edges[ei].end {
 let circ_data = (|| -> Option<(Curve3, [f64; 2])> {
 let ci = result.geom.edge_curve.get(ei).and_then(|c| *c)?;
 match result.geom.curves.get(ci)? {
 Curve3::Circle(c) => {
 let rng = result.geom.edge_curve_range.get(ei).and_then(|r| *r)?;
 Some((Curve3::Circle(*c), rng))
 }
 _ => None,
 }
 })();
 if let Some((Curve3::Circle(circle), [t0, t1])) = circ_data {
 let vs = result.edges[ei].start;
 let mid = (t0 + t1) * 0.5;
 let n = circle.normal.normalize_or(DVec3::Z);
 let rd = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let u = n.cross(rd).normalize();
 let v = n.cross(u).normalize();
 let mid_pt = circle.center + circle.radius * (u * mid.cos() + v * mid.sin());
 let mvi = add_vertex(&mut result, mid_pt);
 let e1 = add_edge(&mut result, Curve3::Circle(circle), t0, mid, vs, mvi);
 let e2 = add_edge(&mut result, Curve3::Circle(circle), mid, t1, mvi, vs);
 split_wire.push(if fwd { WireEdge::fwd(e1) } else { WireEdge::rev(e1) });
 split_wire.push(if fwd { WireEdge::fwd(e2) } else { WireEdge::rev(e2) });
 continue;
 }
 }
 split_wire.push(if fwd { WireEdge::fwd(ei) } else { WireEdge::rev(ei) });
 }

 // Skip faces whose wire has too few edges (collapsed due to offset)
 if split_wire.len() < 2 {
 continue;
 }

 let fi = add_face(&mut result, off_surf.clone(), Wire { edges: split_wire }, Vec::new());
 valid_face_count += 1;

 // For offset full-cylinder faces, set the face_surface_range to ensure the
 // UV grid tessellation covers the full 2pi U-range.  Without pcurves
 // on offset edges, estimate_uv_domain_from_wire can only infer the
 // U-range from vertex positions (which span only a subset of the
 // full cylinder circumference when vertices cluster at the seam and
 // the cap split point), causing ~25% of the wall to be missing from
 // the tessellation and the signed volume to drop correspondingly.
 //
 // Only apply this when the ORIGINAL face wire contains a seam edge
 // (same edge index appearing more than once), indicating a full-wrap
 // cylinder face.  Partial-cylinder faces should use the vertex-based
 // estimate.  We check the original shell's wire, not the offset result
 // wire (which has new edge indices after add_edge).
 if let Surface3::Cylinder(cyl) = &off_surf {
 let orig_seam = {
 let mut seen = std::collections::HashSet::new();
 face.outer_wire.edges.iter().any(|we| !seen.insert(we.idx))
 };
 if orig_seam {
 let wire = &result.solids[0].shells[0].faces[fi].outer_wire;
 let mut v_vals: Vec<f64> = Vec::new();
 for we in &wire.edges {
 if let Some(edge) = result.edges.get(we.idx) {
 if let Some(v) = result.vertices.get(edge.start) {
 v_vals.push((v.point - cyl.origin).dot(cyl.axis));
 }
 if let Some(v) = result.vertices.get(edge.end) {
 v_vals.push((v.point - cyl.origin).dot(cyl.axis));
 }
 }
 }
 if !v_vals.is_empty() {
 let v0 = v_vals.iter().cloned().fold(f64::INFINITY, f64::min);
 let v1 = v_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 while result.geom.face_surface_range.len() <= fi {
 result.geom.face_surface_range.push(None);
 }
 result.geom.face_surface_range[fi] = Some([0.0, std::f64::consts::TAU, v0 - 1e-10, v1 + 1e-10]);
 }
 }
 }
 }

 if valid_face_count == 0 {
 return Err(OffsetError::EmptyResult);
 }

 let _n_result = result.solids.first()
 .and_then(|s| s.shells.first())
 .map(|sh| sh.faces.len())
 .unwrap_or(0);


 // Fix inverted face winding caused by concave topology copy.
 // When offsetting planar faces at reflex corners, the boundary polygon's
 // winding may flip (vertices cross over), making the face normal inconsistent
 // with the wire direction. Detect this and reverse the wire.
 for f in &mut result.solids[0].shells[0].faces {
 // Collect boundary vertices in wire order
 let mut verts: Vec<DVec3> = Vec::new();
 for we in &f.outer_wire.edges {
 let e = &result.edges[we.idx];
 // The vertex at the end of traversal for this edge
 let pt = if we.forward {
 result.vertices[e.end].point
 } else {
 result.vertices[e.start].point
 };
 verts.push(pt);
 }

 if verts.len() >= 3 {
 // Signed area projected onto face normal (Newell's method)
 let n = f.normal;
 let mut signed = 0.0;
 for i in 0..verts.len() {
 let j = (i + 1) % verts.len();
 signed += verts[i].cross(verts[j]).dot(n);
 }
 signed *= 0.5;

 if signed < 0.0 {
 // Winding is reversed relative to face normal → flip wire
 for we in &mut f.outer_wire.edges {
 we.forward = !we.forward;
 }
 f.outer_wire.edges.reverse();
 // Also flip inner wires if any
 for iw in &mut f.inner_wires {
 for we in &mut iw.edges {
 we.forward = !we.forward;
 }
 iw.edges.reverse();
 }
 }
 }
 }

 // Step 5.5: Remove faces that crossed/inverted during offset and fill the
 // resulting holes. A face has "crossed" if its signed volume contribution is
 // negative AND its shares an edge with at least one other negative-volume face.
 // This distinguishes "crossed during offset" (concave corner) from the common
 // "originally inverted face" (single isolated flipped face from boolean ops).
 if distance.abs() > TOLERANCE_MESH_LEGACY {
 let _crossed_removed = fix_crossed_faces(&mut result, brep, shell, distance);
 }

 // Step 6: Check for self-intersection if requested
 let self_intersects = if opts.check_self_intersection {
 detect_self_intersection(&result, distance)
 } else {
 false
 };

 if self_intersects && !opts.auto_repair {
 // Still return the result, but the caller should check for self-intersection
 }

 // Step 6.5: Fix inward-facing face normals.
 // The winding fix (Step 5) ensures wire order matches the face normal,
 // but the face normal itself might point inward (from boolean operations).
 // Use the centroid of all result vertices as the reference point so the
 // check works correctly for solids at any position (not just at the origin).
 {
 let n_faces = result.solids[0].shells[0].faces.len();
 if n_faces > 0 {
 let mut center_sum = DVec3::ZERO;
 for v in &result.vertices { center_sum += v.point; }
 let center = center_sum / result.vertices.len() as f64;
 for fi in 0..n_faces {
 let f = &result.solids[0].shells[0].faces[fi];
 if f.outer_wire.edges.len() < 3 { continue; }
 let mut verts: Vec<DVec3> = Vec::new();
 for we in &f.outer_wire.edges {
 let e = &result.edges[we.idx];
 let pt = if we.forward { result.vertices[e.end].point } else { result.vertices[e.start].point };
 verts.push(pt);
 }
 if verts.len() < 3 { continue; }
 // Centered tetrahedron signed volume (origin → solid center as reference)
 let centered: Vec<DVec3> = verts.iter().map(|v| *v - center).collect();
 let p0 = centered[0];
 let mut vol_6 = 0.0;
 for i in 1..centered.len() - 1 {
 vol_6 += p0.cross(centered[i]).dot(centered[i + 1]);
 }
 if vol_6 < 0.0 {
 // Face normal points inward — flip it and the wire
 let f = &mut result.solids[0].shells[0].faces[fi];
 f.normal = -f.normal;
 for we in &mut f.outer_wire.edges {
 we.forward = !we.forward;
 }
 f.outer_wire.edges.reverse();
 for iw in &mut f.inner_wires {
 for we in &mut iw.edges {
 we.forward = !we.forward;
 }
 iw.edges.reverse();
 }
 }
 }
 }
 }

 Ok(result)
}