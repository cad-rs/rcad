
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
///
/// For single-shell solids (the common case), delegates to offset_shell_with_options.
/// For multi-shell solids, offsets each shell and merges via tshape append.
pub fn offset_solid_with_options(
 solid: &Solid,
 brep: &rcad_kernel::BRep,
 opts: &OffsetOptions,
) -> Result<rcad_kernel::BRep, OffsetError> {
 let distance = opts.distance;

 if distance.abs() < TOLERANCE_LEN_MIN {
  return Err(OffsetError::ZeroDistance);
 }

 // For each shell, compute the offset, then merge tshapes into the result
 let mut result = rcad_kernel::BRep::new();
 for shell in &solid.shells {
  let offset_brep = offset_shell_with_options(shell, brep, opts)?;
  // Merge tshapes: extend the result with the offset tshapes.
  // ShapeRef indices within the offset BRep are self-consistent,
  // so we can just append all tshapes.
  result.tshapes.extend(offset_brep.tshapes);
 }

 if result.tshapes.is_empty() {
  return Err(OffsetError::EmptyResult);
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
 let mut boundary_edges: Vec<usize> = Vec::new();

 for (fi, face) in shell.faces.iter().enumerate() {
  if !open_set.contains(&fi) {
  continue;
  }
  for we in &face.outer_wire.edges {
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

  let surf = match &*brep.tshapes[fi] {
  rcad_kernel::topods::TShape::Face(fd) => fd.surface.as_ref(),
  _ => None,
  };
  let off_surf = surf.and_then(|s| offset_surface(s, inward_offset));

  offset_surfaces.push(off_surf);
 }

 // Step 3: Compute offset vertex positions
 // Exclude open-set (stopper) faces so their normals don't influence
 // boundary vertex positions (e.g., bottom-face vertex at z=0 stays at z=0).
 let offset_vertices: Vec<DVec3> = (0..brep.vertex_count())
 .map(|vi| offset_vertex(brep, vi, inward_offset, shell, Some(&open_set)))
 .collect();

 // Step 4: Build result BRep
 let mut result = rcad_kernel::BRep::new();

 // Add original vertices
 let mut orig_vertex_map: Vec<usize> = Vec::new();
 for vi in 0..brep.vertex_count() {
  let pt = brep.vertex_point(vi).unwrap_or(DVec3::ZERO);
  orig_vertex_map.push(add_vertex(&mut result, pt));
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
  let surf = match &*brep.tshapes[fi] {
  rcad_kernel::topods::TShape::Face(fd) => fd.surface.clone(),
  _ => continue,
  };
  let Some(surf) = surf else { continue };

  let mut wire_edges = Vec::new();
  for we in &face.outer_wire.edges {
  let ed = match &*brep.tshapes[we.idx] { rcad_kernel::topods::TShape::Edge(ed) => ed, _ => continue };
  let vs = orig_vertex_map[ed.first.index];
  let ve = orig_vertex_map[ed.last.index];
  let p0 = result.vertex_point(vs).unwrap_or(DVec3::ZERO);
  let p1 = result.vertex_point(ve).unwrap_or(DVec3::ZERO);
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
  if wire_edges.len() >= 3 {
  add_face(&mut result, surf, Wire { edges: wire_edges }, Vec::new());
  }
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
  let ed = match &*brep.tshapes[we.idx] { rcad_kernel::topods::TShape::Edge(ed) => ed, _ => continue };
  let vs = off_vertex_map[ed.first.index];
  let ve = off_vertex_map[ed.last.index];

  let p0 = result.vertex_point(vs).unwrap_or(DVec3::ZERO);
  let p1 = result.vertex_point(ve).unwrap_or(DVec3::ZERO);
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
  add_face(&mut result, off_surf, Wire { edges: reversed_edges }, Vec::new());
  offset_face_count += 1;
 }

 // Step 6: Create lateral faces along boundary edges
 // Build an edges vec for chain_boundary_edges
 let edges_vec: Vec<Edge> = (0..brep.edge_count()).map(|ei| {
  match &*brep.tshapes[ei] {
  rcad_kernel::topods::TShape::Edge(ed) => Edge { start: ed.first.index, end: ed.last.index },
  _ => Edge { start: 0, end: 0 },
  }
 }).collect();
 let loops = chain_boundary_edges(&boundary_edges, &edges_vec);
 let mut _lateral_count = 0;

 for loop_edges in &loops {
  for &eidx in loop_edges {
  let ed_ref = match &*brep.tshapes[eidx] { rcad_kernel::topods::TShape::Edge(ed) => ed, _ => continue };
  let o_vs = orig_vertex_map[ed_ref.first.index];
  let o_ve = orig_vertex_map[ed_ref.last.index];
  let f_vs = off_vertex_map[ed_ref.first.index];
  let f_ve = off_vertex_map[ed_ref.last.index];

  let p0 = result.vertex_point(o_vs).unwrap_or(DVec3::ZERO);
  let p1 = result.vertex_point(o_ve).unwrap_or(DVec3::ZERO);
  let p3 = result.vertex_point(f_vs).unwrap_or(DVec3::ZERO);

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
  let mut quad_edges = Vec::new();

  for i in 0..4 {
  let s = vseq[i];
  let en = vseq[(i + 1) % 4];
  let sp = result.vertex_point(s).unwrap_or(DVec3::ZERO);
  let ep = result.vertex_point(en).unwrap_or(DVec3::ZERO);
  let dir = (ep - sp).normalize_or(DVec3::X);
  let len = (ep - sp).length();
  let curve = Curve3::Line(Line3 {
  origin: sp,
  direction: dir,
  });
  quad_edges.push(WireEdge::fwd(add_edge(&mut result, curve, 0.0, len, s, en)));
  }

  add_face(&mut result, surf, Wire { edges: quad_edges }, Vec::new());
  _lateral_count += 1;
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

 // Find first face tshape in brep to extract a shell for offset_shell_with_options
 let face_indices: Vec<usize> = brep.tshapes.iter().enumerate()
  .filter(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Face(_)))
  .map(|(i, _)| i)
  .collect();

 if face_indices.is_empty() {
  return Err(OffsetError::InvalidInput("BRep has no faces"));
 }

 // Build an old-style Shell from the face tshape indices.
 // Each face's surface comes from the tshape data, and its wire edges
 // reference edge tshape indices directly.
 let mut shell_faces: Vec<Face> = Vec::new();
 for &fi in &face_indices {
  let fd = match &*brep.tshapes[fi] { rcad_kernel::topods::TShape::Face(fd) => fd, _ => continue };
  // Convert the outer wire from ShapeRef-based to WireEdge-based
  let outer_wire = {
  let wd = match &*brep.tshapes[fd.outer_wire.index] { rcad_kernel::topods::TShape::Wire(wd) => wd, _ => continue };
  let edges: Vec<WireEdge> = wd.edges.iter().map(|sr| {
   WireEdge { idx: sr.index, forward: sr.orientation == rcad_kernel::topods::Orientation::Forward }
  }).collect();
  Wire { edges }
  };
  let normal = fd.surface.as_ref().map(|s| s.normal_at(0.0, 0.0)).unwrap_or(DVec3::Z);
  shell_faces.push(Face {
  outer_wire,
  inner_wires: Vec::new(),
  normal,
  triangles: Vec::new(),
  sample_point: None,
  mesh_dirty: true,
  surface_idx: None,
  });
 }
 let shell = Shell { faces: shell_faces };

 let result_brep = offset_shell_with_options(&shell, brep, &opts)?;

 let self_intersection = if opts.check_self_intersection {
  detect_self_intersection(&result_brep, opts.distance)
 } else {
  false
 };

 let face_count = result_brep.face_count();

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
