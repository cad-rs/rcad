
/// Offset a solid by moving all faces along their normals.
///
/// # Arguments
///
/// * `solid` - The input solid to offset
/// * `brep` - The BRep containing the solid's geometry
/// * `distance` - Offset distance
/// - Positive: outward expansion (thickening)
/// - Negative: inward contraction (shelling)
///
/// # Returns
///
/// A new BRep containing the offset solid, or an error.
pub fn offset_solid(solid: &Solid, brep: &rcad_kernel::BRep, distance: f64) -> Result<rcad_kernel::BRep, OffsetError> {
 offset_solid_with_options(solid, brep, &OffsetOptions::new(distance))
}

/// Offset a solid with full options.
pub fn offset_solid_with_options(
 solid: &Solid,
 brep: &rcad_kernel::BRep,
 opts: &OffsetOptions,
) -> Result<rcad_kernel::BRep, OffsetError> {
 let distance = opts.distance;

 if distance.abs() < TOLERANCE_LEN_MIN {
 return Err(OffsetError::ZeroDistance);
 }

 // For a solid, offset each shell
 let mut result = rcad_kernel::BRep::new();
 result.solids.push(Solid { shells: Vec::new() });

 for shell in &solid.shells {
 let offset_brep = offset_shell_with_options(shell, brep, opts)?;

 // Merge the offset shell into the result
 for offset_solid in offset_brep.solids {
 for offset_shell in offset_solid.shells {
 result.solids[0].shells.push(offset_shell);
 }
 }

 // Merge geometry
 let vertex_offset = result.vertices.len();
 result.vertices.extend(offset_brep.vertices);

 // Remap edge vertex indices
 for edge in offset_brep.edges {
 result.edges.push(Edge {
 start: edge.start + vertex_offset,
 end: edge.end + vertex_offset,
 });
 }

 // Merge geometry store
 let curve_offset = result.geom.curves.len();
 let surface_offset = result.geom.surfaces.len();

 result.geom.curves.extend(offset_brep.geom.curves);
 result.geom.surfaces.extend(offset_brep.geom.surfaces);

 for idx in offset_brep.geom.edge_curve {
 result.geom.edge_curve.push(idx.map(|i| i + curve_offset));
 }
 for range in offset_brep.geom.edge_curve_range {
 result.geom.edge_curve_range.push(range);
 }
 for deg in offset_brep.geom.edge_degenerated {
 result.geom.edge_degenerated.push(deg);
 }
 for idx in offset_brep.geom.face_surface {
 result.geom.face_surface.push(idx.map(|i| i + surface_offset));
 }
 }

 Ok(result)
}

/// Create a hollow solid by removing specified faces and offsetting remaining faces inward.
///
/// This is analogous to the "shell" or "hollow" operation in CAD systems.
///
/// # Arguments
///
/// * `solid` - The input solid
/// * `brep` - The BRep containing the solid's geometry
/// * `thickness` - Wall thickness (positive value)
/// * `open_faces` - Indices of faces to remove (creates openings)
///
/// # Returns
///
/// A new BRep containing the hollow solid with the specified faces removed,
/// or an error.
pub fn hollow_solid(
 solid: &Solid,
 brep: &rcad_kernel::BRep,
 thickness: f64,
 open_faces: &[usize],
) -> Result<rcad_kernel::BRep, OffsetError> {
 hollow_solid_with_options(solid, brep, thickness, open_faces, &OffsetOptions::new(-thickness))
}

/// Create a hollow solid with full options.
pub fn hollow_solid_with_options(
 solid: &Solid,
 brep: &rcad_kernel::BRep,
 thickness: f64,
 open_faces: &[usize],
 opts: &OffsetOptions,
) -> Result<rcad_kernel::BRep, OffsetError> {
 if thickness <= 0.0 {
 return Err(OffsetError::InvalidInput("thickness must be positive"));
 }

 let shell = match solid.shells.first() {
 Some(s) => s,
 None => return Err(OffsetError::InvalidInput("solid has no shells")),
 };

 if open_faces.len() >= shell.faces.len() {
 return Err(OffsetError::InvalidInput("cannot remove all faces"));
 }

 let open_set: HashSet<usize> = open_faces.iter().copied().collect();

 // Step 1: Find boundary edges of the open faces
 let mut edge_use: HashMap<usize, usize> = HashMap::new();
 for (fi, face) in shell.faces.iter().enumerate() {
 if open_set.contains(&fi) {
 continue;
 }
 for we in &face.outer_wire.edges {
 *edge_use.entry(we.idx).or_insert(0) += 1;
 }
 }

 // Boundary edges: edges that were used by removed faces but not by kept faces
 // These are edges where one adjacent face is removed and one is kept
 let mut boundary_edges: Vec<usize> = Vec::new();

 for (fi, face) in shell.faces.iter().enumerate() {
 if !open_set.contains(&fi) {
 continue;
 }
 for we in &face.outer_wire.edges {
 let _e = &brep.edges[we.idx];
 // Check if this edge is shared with a kept face
 let is_shared = shell.faces.iter().enumerate().any(|(fj, fj_face)| {
 !open_set.contains(&fj)
 && fj_face.outer_wire.edges.iter().any(|we2| we2.idx == we.idx)
 });
 if is_shared && !boundary_edges.contains(&we.idx) {
 boundary_edges.push(we.idx);
 }
 }
 }

 // Step 2: Create offset of kept faces (inward offset = negative distance)
 let inward_offset = -thickness;
 let mut offset_opts = opts.clone();
 offset_opts.distance = inward_offset;

 // Compute offset surfaces
 let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
 for (fi, _face) in shell.faces.iter().enumerate() {
 if open_set.contains(&fi) {
 offset_surfaces.push(None);
 continue;
 }

 let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
 Some(s) => s,
 None => {
 offset_surfaces.push(None);
 continue;
 }
 };

 let surf = &brep.geom.surfaces[surf_idx];
 offset_surfaces.push(offset_surface(surf, inward_offset));
 }

 // Step 3: Compute offset vertex positions
 // Exclude open-set (stopper) faces so their normals don't influence
 // boundary vertex positions (e.g., bottom-face vertex at z=0 stays at z=0).
 let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
 .map(|vi| offset_vertex(brep, vi, inward_offset, shell, Some(&open_set)))
 .collect();

 // Step 4: Build result BRep
 let mut result = rcad_kernel::BRep::new();
 result.solids.push(Solid {
 shells: vec![Shell { faces: Vec::new() }],
 });

 // Add original vertices
 let mut orig_vertex_map: Vec<usize> = Vec::new();
 for v in &brep.vertices {
 orig_vertex_map.push(add_vertex(&mut result, v.point));
 }

 // Add offset vertices
 let mut off_vertex_map: Vec<usize> = Vec::new();
 for &p in &offset_vertices {
 off_vertex_map.push(add_vertex(&mut result, p));
 }

 // Step 4.5: Create original kept faces (outer boundary of the thickened wall)
 for (fi, face) in shell.faces.iter().enumerate() {
 if open_set.contains(&fi) {
 continue;
 }
 let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
 Some(s) => s,
 None => continue,
 };
 let surf = &brep.geom.surfaces[surf_idx];

 let mut wire_edges = Vec::new();
 for we in &face.outer_wire.edges {
 let e = &brep.edges[we.idx];
 let vs = orig_vertex_map[e.start];
 let ve = orig_vertex_map[e.end];
 let p0 = result.vertices[vs].point;
 let p1 = result.vertices[ve].point;
 let dir = (p1 - p0).normalize_or(DVec3::X);
 let len = (p1 - p0).length();
 let curve = Curve3::Line(Line3 {
 origin: p0,
 direction: dir,
 });
 let eidx = add_edge(&mut result, curve, 0.0, len, vs, ve);
 wire_edges.push(if we.forward {
 WireEdge::fwd(eidx)
 } else {
 WireEdge::rev(eidx)
 });
 }
 add_face(&mut result, surf.clone(), Wire { edges: wire_edges }, Vec::new());
 }

 // Step 5: Create offset faces for kept faces
 let mut offset_face_count = 0;

 for (fi, face) in shell.faces.iter().enumerate() {
 if open_set.contains(&fi) {
 continue;
 }

 let off_surf = match &offset_surfaces[fi] {
 Some(s) => s.clone(),
 None => continue,
 };

 // Build wire from offset vertices
 let mut wire_edges = Vec::new();

 for we in &face.outer_wire.edges {
 let e = &brep.edges[we.idx];
 let vs = off_vertex_map[e.start];
 let ve = off_vertex_map[e.end];

 let p0 = result.vertices[vs].point;
 let p1 = result.vertices[ve].point;
 let dir = (p1 - p0).normalize_or(DVec3::X);
 let len = (p1 - p0).length();

 let curve = Curve3::Line(Line3 {
 origin: p0,
 direction: dir,
 });

 let eidx = add_edge(&mut result, curve, 0.0, len, vs, ve);
 wire_edges.push(if we.forward { WireEdge::fwd(eidx) } else { WireEdge::rev(eidx) });
 }

 let reversed_edges: Vec<WireEdge> = wire_edges.iter().rev()
 .map(|we| WireEdge { idx: we.idx, forward: !we.forward })
 .collect();
 let off_face_idx = add_face(&mut result, off_surf, Wire { edges: reversed_edges }, Vec::new());
 // Negate the stored normal so face_triangles' orient_tri produces
 // reversed (inward-pointing) triangle winding for the inner boundary.
 result.solids[0].shells[0].faces[off_face_idx].normal *= -1.0;
 offset_face_count += 1;
 }

 // Step 6: Create lateral faces along boundary edges
 let loops = chain_boundary_edges(&boundary_edges, &brep.edges);
 let mut lateral_count = 0;

 for loop_edges in &loops {
 for &eidx in loop_edges {
 let e = &brep.edges[eidx];
 let o_vs = orig_vertex_map[e.start];
 let o_ve = orig_vertex_map[e.end];
 let f_vs = off_vertex_map[e.start];
 let f_ve = off_vertex_map[e.end];

 let p0 = result.vertices[o_vs].point;
 let p1 = result.vertices[o_ve].point;
 let p3 = result.vertices[f_vs].point;

 let normal = (p1 - p0).cross(p3 - p0).normalize_or(DVec3::Z);
 if normal.length() < TOLERANCE_LINEAR_ULTRA_STRICT {
 continue;
 }

 let surf = Surface3::Plane(Plane {
 origin: p0,
 normal,
 });

 // Quad: orig_start -> orig_end -> off_end -> off_start
 let vseq = [o_vs, o_ve, f_ve, f_vs];
 let mut edges = Vec::new();

 for i in 0..4 {
 let s = vseq[i];
 let en = vseq[(i + 1) % 4];
 let dir = (result.vertices[en].point - result.vertices[s].point).normalize_or(DVec3::X);
 let len = (result.vertices[en].point - result.vertices[s].point).length();
 let curve = Curve3::Line(Line3 {
 origin: result.vertices[s].point,
 direction: dir,
 });
 edges.push(WireEdge::fwd(add_edge(&mut result, curve, 0.0, len, s, en)));
 }

 add_face(&mut result, surf, Wire { edges }, Vec::new());
 lateral_count += 1;
 }
 }

 if offset_face_count == 0 {
 return Err(OffsetError::EmptyResult);
 }

 // Triangulate the result
 crate::triangulate::mesh_brep(&mut result, &crate::triangulate::TessellationParams::default());

 Ok(result)
}

/// Offset any BRep shape (shell or solid).
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `opts` - Offset options
///
/// # Returns
///
/// A new BRep with offset geometry.
pub fn offset_shape(brep: &rcad_kernel::BRep, opts: OffsetOptions) -> Result<OffsetResult, OffsetError> {
 if opts.distance.abs() < TOLERANCE_LEN_MIN {
 return Err(OffsetError::ZeroDistance);
 }

 let solid = match brep.solids.first() {
 Some(s) => s,
 None => return Err(OffsetError::InvalidInput("BRep has no solids")),
 };

 let shell = match solid.shells.first() {
 Some(s) => s,
 None => return Err(OffsetError::InvalidInput("solid has no shells")),
 };

 let result_brep = offset_shell_with_options(shell, brep, &opts)?;

 let self_intersection = if opts.check_self_intersection {
 detect_self_intersection(&result_brep, opts.distance)
 } else {
 false
 };

 let face_count = result_brep
 .solids
 .first()
 .and_then(|s| s.shells.first())
 .map(|sh| sh.faces.len())
 .unwrap_or(0);

 Ok(OffsetResult {
 brep: result_brep,
 offset_faces: face_count,
 lateral_faces: 0,
 join_faces: 0,
 self_intersection,
 quality: OffsetQuality {
 min_wall_thickness: f64::INFINITY,
 max_deviation: 0.0,
 degenerate_edge_count: 0,
 self_intersection_count: if self_intersection { 1 } else { 0 },
 face_area_ratio: 1.0,
 edge_length_ratio: 1.0,
 is_valid: true,
 warnings: Vec::new(),
 },
 warnings: Vec::new(),
 effective_distance: opts.distance,
 repair_attempts: 0,
 })
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Tests
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
