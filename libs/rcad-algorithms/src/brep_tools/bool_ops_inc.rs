// =============================================================================
// Shell / Solid Extraction (explode equivalent)
// =============================================================================

/// Remove stale vertices/edges and rebuild with dense indexing.
///
/// After boolean operations, the result BRep may retain vertices from both
/// inputs that are not part of the result geometry.  These stale vertices
/// inflate [`bounding_box`] and can cause subsequent booleans to produce
/// wrong results (e.g. via [`try_containment`](crate::try_containment)).
///
/// This function rebuilds the BRep using only vertices and edges that are
/// referenced by at least one face wire, producing a minimal self-contained
/// copy with correct bounding box.
pub(crate) fn compact_brep(brep: &rcad_kernel::BRep) -> rcad_kernel::BRep {
 crate::brep_tools::compact_brep_topods(brep)
}

/// Create a new self-contained BRep containing only the specified flat face
/// indices from the source BRep.  Vertices, edges, and geometry referenced by
/// the selected faces are copied into the new BRep with dense index renumbering.
fn extract_brep_subset(source: &rcad_kernel::BRep, face_indices: &[usize]) -> rcad_kernel::BRep {
 compact_brep_face_subset(source, face_indices)
}

/// Extract each solid from a (possibly compound) BRep as a separate
/// self-contained BRep.  Equivalent to OCCT `explode ... so`.
///
/// Each returned BRep has only the vertices, edges, and geometry belonging
/// to that solid, with indices renumbered from 0.
pub fn extract_solids(brep: &rcad_kernel::BRep) -> Vec<rcad_kernel::BRep> {
 crate::brep_tools::extract_solids_topods(brep)
}

/// Extract each shell from a BRep as a separate self-contained BRep.
/// Equivalent to OCCT `explode ... Sh`.
///
/// Each returned BRep has only the vertices, edges, and geometry belonging
/// to that shell, with indices renumbered from 0.
pub fn extract_shells(brep: &rcad_kernel::BRep) -> Vec<rcad_kernel::BRep> {
 crate::brep_tools::extract_shells_topods(brep)
}

/// Partition objects by tools using boolean-subset decomposition.
///
/// For each object and each combination of tools (inside/outside per tool mask),
/// computes the corresponding cell using pairwise boolean operations.
/// Returns all non-empty cells as individual self-contained BReps (one solid each).
///
/// This is equivalent to OCCT's `BRepAlgoAPI_Splitter` / `BRepAlgoAPI_Partition`.
///
/// Face tools (planar faces acting as half-space dividers) are automatically
/// expanded into two half-space solids. Face objects (zero-volume faces) are
/// partitioned via `split_shape` + point classification.
///
/// The number of boolean operations per call is O(objects.len() 2^n_tools n_tools),
/// so this is suitable only for small numbers of tools ( ?10).
pub fn n_ary_partition(objects: &[rcad_kernel::BRep], tools: &[rcad_kernel::BRep]) -> Result<Vec<rcad_kernel::BRep>, crate::BooleanError> {
 let mut all_cells = Vec::new();

 for obj in objects {
 if is_face_like(obj) {
 all_cells.extend(partition_face_object(obj, tools)?);
 } else {
 all_cells.extend(partition_solid_object(obj, tools)?);
 }
 }

 all_cells.retain(|c| count_faces(c) > 0);
 Ok(all_cells)
}

/// Decompose the complement of `inner_bbox` within `outer_bbox` into up to 6
/// non-overlapping axis-aligned boxes.
///
/// Each returned tuple is `(origin, u_dir, v_dir, width, height, depth)` suitable
/// for passing to `make_box_brep`.  The boxes are disjoint and their union is
/// exactly `outer_bbox \ inner_bbox` (the region inside the outer box but outside
/// the inner box).
fn box_complement_of_bbox(
 inner: &[DVec3; 2],
 outer: &[DVec3; 2],
) -> Vec<(DVec3, DVec3, DVec3, f64, f64, f64)> {
 let (omin, omax) = (outer[0], outer[1]);
 let (imin, imax) = (inner[0], inner[1]);

 // Clamp inner bbox to outer bbox so the complement doesn't extend past the outer.
 let imin = imin.max(omin);
 let imax = imax.min(omax);

 let mut boxes = Vec::with_capacity(6);

 // Left: x < imin.x (full y,z range of outer)
 if imin.x > omin.x {
 boxes.push((
 DVec3::new(omin.x, omin.y, omin.z),
 DVec3::X,
 DVec3::Y,
 imin.x - omin.x,
 omax.y - omin.y,
 omax.z - omin.z,
 ));
 }

 // Right: x > imax.x (full y,z range of outer)
 if omax.x > imax.x {
 boxes.push((
 DVec3::new(imax.x, omin.y, omin.z),
 DVec3::X,
 DVec3::Y,
 omax.x - imax.x,
 omax.y - omin.y,
 omax.z - omin.z,
 ));
 }

 // Front: y < imin.y (within tool's x range, full z range)
 if imin.y > omin.y {
 boxes.push((
 DVec3::new(imin.x, omin.y, omin.z),
 DVec3::X,
 DVec3::Y,
 imax.x - imin.x,
 imin.y - omin.y,
 omax.z - omin.z,
 ));
 }

 // Back: y > imax.y (within tool's x range, full z range)
 if omax.y > imax.y {
 boxes.push((
 DVec3::new(imin.x, imax.y, omin.z),
 DVec3::X,
 DVec3::Y,
 imax.x - imin.x,
 omax.y - imax.y,
 omax.z - omin.z,
 ));
 }

 // Bottom: z < imin.z (within tool's x,y range)
 if imin.z > omin.z {
 boxes.push((
 DVec3::new(imin.x, imin.y, omin.z),
 DVec3::X,
 DVec3::Y,
 imax.x - imin.x,
 imax.y - imin.y,
 imin.z - omin.z,
 ));
 }

 // Top: z > imax.z (within tool's x,y range)
 if omax.z > imax.z {
 boxes.push((
 DVec3::new(imin.x, imin.y, imax.z),
 DVec3::X,
 DVec3::Y,
 imax.x - imin.x,
 imax.y - imin.y,
 omax.z - imax.z,
 ));
 }

 boxes
}

/// Check whether a BRep is a simple axis-aligned box (its volume matches its
/// bounding-box volume within tolerance).
fn is_box_like(brep: &rcad_kernel::BRep) -> bool {
 let vol = crate::total_volume(brep);
 if vol <= 0.0 {
 return false;
 }
 if let Some(bbox) = bounding_box(brep) {
 let diag = bbox[1] - bbox[0];
 let bbox_vol = diag.x * diag.y * diag.z;
 if bbox_vol <= 0.0 {
 return false;
 }
 (vol - bbox_vol).abs() < 1e-6
 } else {
 false
 }
}

/// Partition a solid object by tools, expanding face tools into half-space solids.
fn partition_solid_object(obj: &rcad_kernel::BRep, tools: &[rcad_kernel::BRep]) -> Result<Vec<rcad_kernel::BRep>, crate::BooleanError> {
 // Detect face-like tools (planar faces used as dividing surfaces).
 let face_tool_info: Vec<Option<rcad_kernel::geom::Plane>> = tools
 .iter()
 .map(|t| if is_face_like(t) { try_as_planar_face(t) } else { None })
 .collect();
 let has_face_tool = face_tool_info.iter().any(|p| p.is_some());

 // Expand tools + track complement indices for half-spaces.
 //
 // For a half-space pair (h_plus, h_minus), each half-space's "outside"
 // is simply Intersection with the OTHER half-space  ?no Diff needed.
 // expanded_complements[i] = Some(j) means bit-flip (outside) at index i
 // can be handled by Intersection with expanded_tools[j].
 //
 // For solid tools (non-face), expanded_complements[i] = None  ?those
 // may need the complement-box fallback or Diff.
 let mut expanded_tools: Vec<rcad_kernel::BRep> = Vec::new();
 #[allow(clippy::type_complexity)]
 let mut expanded_complements: Vec<Option<usize>> = Vec::new();

 if has_face_tool {
 let mut bbox = bounding_box(obj);
 for (ti, tool) in tools.iter().enumerate() {
 if face_tool_info[ti].is_some() {
 if let Some(tb) = bounding_box(tool) {
 bbox = match bbox {
 Some(b) => Some([b[0].min(tb[0]), b[1].max(tb[1])]),
 None => Some(tb),
 };
 }
 }
 }

 for (ti, tool) in tools.iter().enumerate() {
 if let Some(ref plane) = face_tool_info[ti] {
 if let Some(b) = bbox {
 let h_plus_idx = expanded_tools.len();
  let h_plus = make_face_half_space(plane, &b, true);
  let h_minus = make_face_half_space(plane, &b, false);
 expanded_tools.push(h_plus);
 expanded_tools.push(h_minus);
 // h_plus  ?h_minus: each is the "outside" of the other.
 expanded_complements.push(Some(h_plus_idx + 1));
 expanded_complements.push(Some(h_plus_idx));
 continue;
 }
 }
 expanded_complements.push(None);
 expanded_tools.push(tool.clone());
 }
 } else {
 for tool in tools {
 expanded_complements.push(None);
 expanded_tools.push(tool.clone());
 }
 }

 let n_tools = expanded_tools.len();
 let mut cells = Vec::new();
 let max_mask = if n_tools >= 32 { 1u32 << 31 } else { 1u32 << n_tools };

 // Detect complementary half-space pairs (h_plus/h_minus).
 // When both bits of a pair are equal (both 0 or both 1), the cell lies on
 // the plane between them and has zero volume  ?skip those masks.
 let mut comp_pairs: Vec<(usize, usize)> = Vec::new();
 for (i, comp_i) in expanded_complements.iter().enumerate() {
 if let Some(j) = comp_i {
 if *j > i && expanded_complements.get(*j) == Some(&Some(i)) {
 comp_pairs.push((i, *j));
 }
 }
 }

 for mask in 0..max_mask {
 // Skip masks where complementary half-spaces have the same bit
 // (both inside or both outside = intersection is just the plane).
 if comp_pairs.iter().any(|&(i, j)| {
 ((mask >> i) & 1) == ((mask >> j) & 1)
 }) {
 continue;
 }
 let mut cell = obj.clone();
 let mut failed = false;
 let mut first_tool = true;

 for i in 0..n_tools {
 let inside = (mask >> i) & 1 != 0;
 let tool = &expanded_tools[i];

 if inside {
 match crate::boolean_op_pave_fill_build(crate::BooleanOpType::Intersection, &cell, tool) {
 Ok(r) => { cell = compact_brep(&r); },
 Err(_) => { failed = true; break; }
 }
 } else if let Some(complement_idx) = expanded_complements[i] {
 // Half-space: "outside" = Intersection with complementary half-space.
 let complement = &expanded_tools[complement_idx];
 match crate::boolean_op_pave_fill_build(crate::BooleanOpType::Intersection, &cell, complement) {
 Ok(r) => { cell = compact_brep(&r); },
 Err(_) => { failed = true; break; }
 }
 } else {
 // Solid tool: for the first tool application, use complement box
 // decomposition to avoid coincident-face issues with Diff.
 // For subsequent tools, use Diff (works better when cell is already
 // a multi-solid compound from previous operations).
 if first_tool && is_box_like(tool) {
 // First tool: use complement box decomposition.
 if let (Some(tool_bbox), Some(cell_bbox)) =
 (bounding_box(tool), bounding_box(&cell))
 {
 let comp_boxes =
 box_complement_of_bbox(&tool_bbox, &cell_bbox);
 if comp_boxes.is_empty() {
 cell = rcad_kernel::BRep::new();
 break;
 }
 let cell_solids = extract_solids(&cell);
 let mut parts = Vec::new();
 for (origin, u_dir, v_dir, w, h, d) in comp_boxes.iter() {
  let Ok(comp_box) =
  rcad_modeling::make_box_brep(*origin, *u_dir, *v_dir, *w, *h, *d)
  else { continue; };
  let comp_box_old = comp_box;
  for cell_part in &cell_solids {
  if let Ok(part) =
  crate::boolean_op_pave_fill_build(crate::BooleanOpType::Intersection, cell_part, &comp_box_old)
 {
 let p = part;
 if count_faces(&p) > 0 {
 parts.push(p);
 }
 }
 }
 }
 cell = rcad_kernel::BRep::compound_from_shapes(&parts);
 first_tool = false;
 continue;
 }
 }
 // Subsequent tool or non-box tool: use Diff.
 match crate::boolean_op_pave_fill_build(crate::BooleanOpType::Difference, &cell, tool) {
 Ok(r) => { cell = compact_brep(&r); }
 Err(_) => { failed = true; break; }
 }
 }
 first_tool = false;
 }

 if failed {
 continue;
 }

 let face_indices = collect_flat_face_indices(&cell);
 if face_indices.is_empty() {
 continue;
 }

 let comps = connected_face_components(&cell, &face_indices);
 for component in comps {
 if !component.is_empty() {
 let subset = extract_brep_subset(&cell, &component);
 cells.push(subset);
 }
 }
 }

 // Filter out degenerate cells (zero volume from tangent-coincident masks).
 cells.retain(|c| crate::total_volume(c) > 1e-10);

 Ok(cells)
}

/// Partition a face-like object by tools using split_shape and centroid classification.
fn partition_face_object(obj: &rcad_kernel::BRep, tools: &[rcad_kernel::BRep]) -> Result<Vec<rcad_kernel::BRep>, crate::BooleanError> {
 // Collect solid (non-face) tools.
 let solid_tools: Vec<rcad_kernel::BRep> = tools.iter().filter(|t| !is_face_like(t)).cloned().collect();

 if solid_tools.is_empty() || count_faces(obj) == 0 {
 return Ok(vec![obj.clone()]);
 }

 let orig_surface_info = get_surface(obj, 0).ok().map(|s| match s {
 Surface3::Plane(p) => (p.origin, p.normal, "plane"),
 Surface3::Cylinder(_) => (DVec3::ZERO, DVec3::Z, "cylinder"),
 Surface3::Sphere(_) => (DVec3::ZERO, DVec3::Z, "sphere"),
 Surface3::Cone(_) => (DVec3::ZERO, DVec3::Z, "cone"),
 Surface3::Torus(_) => (DVec3::ZERO, DVec3::Z, "torus"),
 _ => (DVec3::ZERO, DVec3::Z, "other"),
 });

 let collect_on_plane = |brep: &rcad_kernel::BRep| -> Vec<usize> {
 let Some((plane_origin, plane_normal, _)) = orig_surface_info else { return vec![] };
 let mut out = Vec::new();
 let mut fi = 0usize;
 for ts in &brep.tshapes {
 if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
 for shell_sr in &sd.shells {
 if let rcad_kernel::topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
 for _ in &shd.faces {
 if let Ok(Surface3::Plane(p)) = get_surface(brep, fi) {
 let dist = (p.origin - plane_origin).dot(plane_normal).abs();
 if dist < 1e-6 && p.normal.dot(plane_normal) > 0.9999 {
 out.push(fi);
 }
 }
 fi += 1;
 }
 }
 }
 }
 }
 out
 };

 // Use boolean ops to carve the face into inside/outside per tool.
 let mut cells = Vec::new();
 let mut remaining = obj.clone();
 for tool in &solid_tools {
 let inside_t = crate::boolean_op(crate::BooleanOpType::Intersection, &remaining, tool)?;
 let inside = inside_t;
 let in_faces = collect_on_plane(&inside);
 if !in_faces.is_empty() {
 cells.push(extract_brep_subset(&inside, &in_faces));
 }
 let remaining_t = crate::boolean_op(crate::BooleanOpType::Difference, &remaining, tool)?;
 remaining = remaining_t;
 }
 let out_faces = collect_on_plane(&remaining);
 if !out_faces.is_empty() {
 cells.push(extract_brep_subset(&remaining, &out_faces));
 }

 Ok(cells)
}

/// Check if a BRep represents a single planar face and extract its plane.
fn try_as_planar_face(brep: &rcad_kernel::BRep) -> Option<rcad_kernel::geom::Plane> {
 if count_faces(brep) != 1 {
 return None;
 }
 match get_surface(brep, 0).ok()? {
 Surface3::Plane(plane) => Some(plane.clone()),
 _ => None,
 }
}

/// Check if a BRep is face-like (open surface, not a proper 3D solid).
///
/// A BRep is face-like if every shell contains exactly one planar face. Proper 3D
/// solids always have at least 4 faces per shell (minimum tetrahedron), except
/// analytic primitives like spheres/cones/cylinders which may have only 1-3 faces.
fn is_face_like(brep: &rcad_kernel::BRep) -> bool {
 if count_faces(brep) == 0 {
 return false;
 }
 let mut flat_idx = 0usize;
 for ts in &brep.tshapes {
 if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
 for shell_sr in &sd.shells {
 if let rcad_kernel::topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
 if shd.faces.len() != 1 {
 return false;
 }
 let surface = get_surface(brep, flat_idx).ok();
 if !matches!(surface, Some(Surface3::Plane(_))) {
 return false;
 }
 flat_idx += 1;
 }
 }
 }
 }
 true
}

/// Create a half-space solid extending from a plane in the normal (or opposite) direction.
///
/// The resulting box occupies a prism with one face exactly on the plane through
/// `plane.origin` (with the plane's normal pointing inward for `normal_side=true`
/// or outward for `normal_side=false`) and extending far enough along the normal
/// to fully contain the `bbox` extent.
pub fn make_face_half_space(plane: &rcad_kernel::geom::Plane, bbox: &[DVec3; 2], normal_side: bool) -> topods::BRep {
 let [bmin, bmax] = *bbox;
 let diag = bmax - bmin;
 let margin = diag.length().max(1.0) * 2.0;

 let n = if normal_side { plane.normal } else { -plane.normal };
 let n = n.normalize();

 // Build a tangent basis in the plane.
 let abs = n.abs();
 let candidate = if abs.x <= abs.y && abs.x <= abs.z {
 DVec3::X
 } else if abs.y <= abs.z {
 DVec3::Y
 } else {
 DVec3::Z
 };
 let u = n.cross(candidate).normalize();
 let v = n.cross(u);

 // Build the box positioned so it starts at the plane and extends `margin`
 // along +n (normal_side=true) or -n (normal_side=false).
 //
 // The box origin is the corner at (-u*margin/2, -v*margin/2),
 // extended to (+u*margin/2, +v*margin/2) in the plane, and from
 // the plane to  margin along n.
 let origin = if normal_side {
 plane.origin - u * (margin / 2.0) - v * (margin / 2.0)
 } else {
 plane.origin - u * (margin / 2.0) - v * (margin / 2.0) - n * margin
 };

 rcad_modeling::make_box_brep(origin, u, v, margin, margin, margin)
 .expect("make_face_half_space: box construction should not fail")
}

/// Compute the average vertex position of a face's boundary using topods access.
fn average_vertex_of_face(brep: &rcad_kernel::BRep, face_sr: &rcad_kernel::topods::ShapeRef) -> DVec3 {
 let mut sum = DVec3::ZERO;
 let mut count = 0usize;
 if let rcad_kernel::topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
 let walk_edges = |wire_sr: &rcad_kernel::topods::ShapeRef, sum: &mut DVec3, count: &mut usize| {
 if let rcad_kernel::topods::TShape::Wire(wd) = &*brep.tshapes[wire_sr.index] {
 for e_sr in &wd.edges {
 if let rcad_kernel::topods::TShape::Edge(ed) = &*brep.tshapes[e_sr.index] {
 if let rcad_kernel::topods::TShape::Vertex(vd) = &*brep.tshapes[ed.first.index] {
 *sum += vd.point;
 *count += 1;
 }
 if let rcad_kernel::topods::TShape::Vertex(vd) = &*brep.tshapes[ed.last.index] {
 *sum += vd.point;
 *count += 1;
 }
 }
 }
 }
 };
 walk_edges(&fd.outer_wire, &mut sum, &mut count);
 for iw_sr in &fd.inner_wires {
 walk_edges(iw_sr, &mut sum, &mut count);
 }
 }
 if count > 0 { sum / count as f64 } else { DVec3::ZERO }
}

/// Collect flat face indices for all faces in a BRep.
fn collect_flat_face_indices(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let mut indices = Vec::new();
 let mut flat_idx = 0;
 for ts in &brep.tshapes {
 if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
 for shell_sr in &sd.shells {
 if let rcad_kernel::topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
 for _ in &shd.faces {
 indices.push(flat_idx);
 flat_idx += 1;
 }
 }
 }
 }
 }
 indices
}

/// Find connected components of a set of flat face indices within a BRep.
/// Two faces are connected if they share at least one edge (same edge tshape index).
fn connected_face_components(brep: &rcad_kernel::BRep, face_indices: &[usize]) -> Vec<Vec<usize>> {
 use std::collections::{HashMap, HashSet};

 let face_set: HashSet<usize> = face_indices.iter().copied().collect();
 if face_set.is_empty() {
 return Vec::new();
 }

 // Build edge → face list for our faces of interest.
 let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
 let mut flat_idx: usize = 0;

 for ts in &brep.tshapes {
 if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
 for shell_sr in &sd.shells {
 if let rcad_kernel::topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
 for lfi in 0..shd.faces.len() {
 let global_fi = flat_idx + lfi;
 if face_set.contains(&global_fi) {
 let face_sr = &shd.faces[lfi];
 if let rcad_kernel::topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
 let collect_edges = |wire_sr: &rcad_kernel::topods::ShapeRef, e2f: &mut HashMap<usize, Vec<usize>>| {
 if let rcad_kernel::topods::TShape::Wire(wd) = &*brep.tshapes[wire_sr.index] {
 for e_sr in &wd.edges {
 e2f.entry(e_sr.index).or_default().push(global_fi);
 }
 }
 };
 collect_edges(&fd.outer_wire, &mut edge_to_faces);
 for iw_sr in &fd.inner_wires {
 collect_edges(iw_sr, &mut edge_to_faces);
 }
 }
 }
 }
 flat_idx += shd.faces.len();
 }
 }
 }
 }

 // Build adjacency: face A  ?[faces that share an edge with A].
 let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
 for faces in edge_to_faces.values() {
 if faces.len() >= 2 {
 for i in 0..faces.len() {
 for j in (i + 1)..faces.len() {
 adjacency.entry(faces[i]).or_default().push(faces[j]);
 adjacency.entry(faces[j]).or_default().push(faces[i]);
 }
 }
 }
 }

 // DFS over face indices to find connected components.
 let mut visited: HashSet<usize> = HashSet::new();
 let mut components: Vec<Vec<usize>> = Vec::new();

 for &fi in face_indices {
 if !visited.insert(fi) {
 continue;
 }

 let mut component: Vec<usize> = Vec::new();
 let mut stack: Vec<usize> = vec![fi];

 while let Some(current) = stack.pop() {
 component.push(current);
 if let Some(neighbors) = adjacency.get(&current) {
 for &neighbor in neighbors {
 if visited.insert(neighbor) {
 stack.push(neighbor);
 }
 }
 }
 }

 components.push(component);
 }

 components
}
