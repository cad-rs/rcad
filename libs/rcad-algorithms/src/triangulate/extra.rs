
fn weld_surface_mesh_nodes_with_exclusion(
 mesh: &SurfaceMesh,
 excluded_nodes: &std::collections::HashSet<usize>,
) -> SurfaceMesh {
 const WELD_TOLERANCE: f64 = TOLERANCE_COORD_SUB;

 let mut remap = vec![0usize; mesh.nodes.len()];
 let mut welded_nodes: Vec<DVec3> = Vec::new();
 let mut welded_normals: Vec<DVec3> = Vec::new();
 let mut normal_counts = Vec::new();
 let mut buckets: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
 let scale = 1.0 / WELD_TOLERANCE;

 for (index, point) in mesh.nodes.iter().enumerate() {
 // Pinned nodes bypass spatial hashing
 if excluded_nodes.contains(&index) {
 let new_index = welded_nodes.len();
 welded_nodes.push(*point);
 welded_normals.push(mesh.normals.get(index).copied().unwrap_or(DVec3::ZERO));
 normal_counts.push(1);
 remap[index] = new_index;
 continue;
 }

 let key = [
 (point.x * scale).round() as i64,
 (point.y * scale).round() as i64,
 (point.z * scale).round() as i64,
 ];

 let mut matched = None;
 if let Some(candidates) = buckets.get(&key) {
 for &candidate in candidates {
 if excluded_nodes.contains(&candidate) {
 continue;
 }
 if (welded_nodes[candidate] - *point).length_squared() <= WELD_TOLERANCE * WELD_TOLERANCE {
 matched = Some(candidate);
 break;
 }
 }
 }

 let target = if let Some(existing) = matched {
 existing
 } else {
 let new_index = welded_nodes.len();
 welded_nodes.push(*point);
 welded_normals.push(DVec3::ZERO);
 normal_counts.push(0);
 buckets.entry(key).or_default().push(new_index);
 new_index
 };

 remap[index] = target;
 if let Some(normal) = mesh.normals.get(index) {
 welded_normals[target] += *normal;
 normal_counts[target] += 1;
 }
 }

 let welded_triangles: Vec<[usize; 3]> = mesh
 .triangles
 .iter()
 .filter_map(|&[a, b, c]| {
 let ra = remap[a];
 let rb = remap[b];
 let rc = remap[c];
 if ra == rb || rb == rc || rc == ra {
 None
 } else {
 Some([ra, rb, rc])
 }
 })
 .collect();

 let welded_normals: Vec<DVec3> = welded_normals
 .into_iter()
 .zip(normal_counts)
 .map(|(normal, count)| {
 if count == 0 {
 DVec3::ZERO
 } else {
 normal.normalize_or_zero()
 }
 })
 .collect();

 SurfaceMesh {
 nodes: welded_nodes,
 triangles: welded_triangles,
 normals: welded_normals,
 dirty: mesh.dirty,
 }
}

// ============================================================================
// Incremental mesh updates
// ============================================================================

/// Lightweight change list for dirty tracking.
#[derive(Debug, Clone, Default)]
pub struct MeshDelta {
 /// Vertex indices that changed.
 pub modified_vertices: Vec<usize>,
 /// Edge indices that changed.
 pub modified_edges: Vec<usize>,
 /// Flattened global face indices that changed.
 pub modified_faces: Vec<usize>,
}

impl MeshDelta {
 /// Empty delta.
 pub fn new() -> Self {
 Self::default()
 }

 /// Delta containing only vertex edits.
 pub fn from_vertices(vertices: Vec<usize>) -> Self {
 Self {
 modified_vertices: vertices,
 ..Default::default()
 }
 }

 /// Delta containing only edge edits.
 pub fn from_edges(edges: Vec<usize>) -> Self {
 Self {
 modified_edges: edges,
 ..Default::default()
 }
 }

 /// Delta containing only face edits.
 pub fn from_faces(faces: Vec<usize>) -> Self {
 Self {
 modified_faces: faces,
 ..Default::default()
 }
 }

 /// `true` if no topology elements were touched.
 pub fn is_empty(&self) -> bool {
 self.modified_vertices.is_empty()
 && self.modified_edges.is_empty()
 && self.modified_faces.is_empty()
 }
}

/// Tracks which `BRep` faces need re-tessellation after local edits.
#[derive(Debug, Clone, Default)]
pub struct IncrementalMesher {
 /// Flattened face indices pending retessellation.
 pub dirty_faces: std::collections::HashSet<usize>,
 /// Edge indices whose incident faces should be refreshed.
 pub dirty_edges: std::collections::HashSet<usize>,
 /// Vertex indices whose incident faces should be refreshed.
 pub dirty_vertices: std::collections::HashSet<usize>,
}

impl IncrementalMesher {
 /// Empty dirty set.
 pub fn new() -> Self {
 Self::default()
 }

 /// Mark one flattened face index dirty.
 pub fn invalidate_face(&mut self, face_idx: usize) {
 self.dirty_faces.insert(face_idx);
 }

 /// Mark many face indices dirty.
 pub fn invalidate_faces(&mut self, face_indices: &[usize]) {
 for &idx in face_indices {
 self.dirty_faces.insert(idx);
 }
 }

 /// Mark an edge dirty (incident faces inferred separately).
 pub fn invalidate_edge(&mut self, edge_idx: usize) {
 self.dirty_edges.insert(edge_idx);
 }

 /// Mark a vertex dirty.
 pub fn invalidate_vertex(&mut self, vertex_idx: usize) {
 self.dirty_vertices.insert(vertex_idx);
 }

 /// Expand `dirty_faces` from explicit faces/edges/vertices in `delta`.
 pub fn infer_dirty_faces_from_delta(&mut self, brep: &BRep, delta: &MeshDelta) {
 // Faces named explicitly in the delta
 self.invalidate_faces(&delta.modified_faces);

 // Faces incident on modified edges
 for &edge_idx in &delta.modified_edges {
 if let Some(_edge) = brep.edges.get(edge_idx) {
 // Faces incident on this edge
 let mut face_idx = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if we.idx == edge_idx {
 self.dirty_faces.insert(face_idx);
 }
 }
 face_idx += 1;
 }
 }
 }
 }
 }

 // Faces touching modified vertices
 for &vertex_idx in &delta.modified_vertices {
 let mut face_idx = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if let Some(edge) = brep.edges.get(we.idx) {
 let start = if we.forward { edge.start } else { edge.end };
 let end = if we.forward { edge.end } else { edge.start };
 if start == vertex_idx || end == vertex_idx {
 self.dirty_faces.insert(face_idx);
 }
 }
 }
 face_idx += 1;
 }
 }
 }
 }
 }

 /// Flag dirty faces on `brep` then invoke `mesh_brep` to refresh tessellation.
 pub fn update_mesh_for_face_change(
 &self,
 brep: &mut BRep,
 params: &TessellationParams,
 ) {
 if self.dirty_faces.is_empty() {
 return;
 }

 let mut face_flat_idx = 0usize;

 for solid_idx in 0..brep.solids.len() {
 for shell_idx in 0..brep.solids[solid_idx].shells.len() {
 let n_faces = brep.solids[solid_idx].shells[shell_idx].faces.len();
 for face_idx in 0..n_faces {
 if self.dirty_faces.contains(&face_flat_idx) {
 brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
 .mesh_dirty = true;
 }
 face_flat_idx += 1;
 }
 }
 }

 mesh_brep(brep, params);
 }

 /// Clear all dirty flags.
 pub fn clear(&mut self) {
 self.dirty_faces.clear();
 self.dirty_edges.clear();
 self.dirty_vertices.clear();
 }

 /// `true` if any dirty set is non-empty.
 pub fn is_dirty(&self) -> bool {
 !self.dirty_faces.is_empty()
 || !self.dirty_edges.is_empty()
 || !self.dirty_vertices.is_empty()
 }
}

// ============================================================================
// Mesh simplification (edge collapse)
// ============================================================================

/// Candidate half-edge collapse with a cheap length-based metric.
#[derive(Debug, Clone)]
struct EdgeCollapseInfo {
 /// Canonical undirected edge `(min(v0,v1), max(v0,v1))`.
 edge: (usize, usize),
 /// Collapse priority (here: edge length).
 error: f64,
 /// Vertex position after collapsing `edge` onto its midpoint.
 new_position: DVec3,
}

/// Very small edge-collapse helper for decimating `SurfaceMesh` data.
#[derive(Debug, Clone)]
pub struct MeshSimplifier {
 /// Fraction of triangles to keep (`0.0`= 1.0`).
 pub target_ratio: f64,
 /// Skip collapses longer than this edge length.
 pub max_error: f64,
 /// When `true`, do not collapse edges on the mesh boundary.
 pub preserve_boundary: bool,
}

impl Default for MeshSimplifier {
 fn default() -> Self {
 Self {
 target_ratio: 0.5,
 max_error: 0.01,
 preserve_boundary: true,
 }
 }
}

impl MeshSimplifier {
 /// Default simplifier (retain ~50% of triangles, short-edge priority).
 pub fn new() -> Self {
 Self::default()
 }

 /// Builder: clamped [`Self::target_ratio`].
 pub fn with_target_ratio(mut self, ratio: f64) -> Self {
 self.target_ratio = ratio.clamp(0.0, 1.0);
 self
 }

 /// Builder: [`Self::max_error`] (max edge length eligible for collapse).
 pub fn with_max_error(mut self, error: f64) -> Self {
 self.max_error = error;
 self
 }

 /// Convenience wrapper that derives `target_ratio` from `target_count`.
 pub fn simplify_to_target_count(&self, mesh: &SurfaceMesh, target_count: usize) -> SurfaceMesh {
 if mesh.triangles.len() <= target_count {
 return mesh.clone();
 }

 let ratio = target_count as f64 / mesh.triangles.len() as f64;
 Self {
 target_ratio: ratio,
 ..self.clone()
 }
 .simplify_mesh(mesh)
 }

 /// Greedy edge collapses until the target triangle count is reached.
 pub fn simplify_mesh(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
 if mesh.triangles.is_empty() {
 return mesh.clone();
 }

 let target_triangle_count = (mesh.triangles.len() as f64 * self.target_ratio).max(4.0) as usize;

 let mut nodes = mesh.nodes.clone();
 let mut normals = mesh.normals.clone();
 let mut triangles = mesh.triangles.clone();

 // Boundary nodes (degree-1 edges) when preserving openings
 let boundary_nodes = if self.preserve_boundary {
 find_boundary_nodes(&triangles)
 } else {
 std::collections::HashSet::new()
 };

 // Greedy collapses
 while triangles.len() > target_triangle_count {
 let collapse = find_best_edge_collapse(
 &nodes,
 &triangles,
 &boundary_nodes,
 self.max_error,
 );

 let Some(collapse) = collapse else {
 break;
 };

 // Collapse shortest eligible edge
 apply_edge_collapse(
 &mut nodes,
 &mut normals,
 &mut triangles,
 collapse.edge,
 collapse.new_position,
 );

 if triangles.len() <= target_triangle_count {
 break;
 }
 }

 SurfaceMesh {
 nodes,
 triangles,
 normals,
 dirty: false,
 }
 }
}

fn find_boundary_nodes(triangles: &[[usize; 3]]) -> std::collections::HashSet<usize> {
 let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();

 for &tri in triangles {
 let edges = [
 (tri[0].min(tri[1]), tri[0].max(tri[1])),
 (tri[1].min(tri[2]), tri[1].max(tri[2])),
 (tri[2].min(tri[0]), tri[2].max(tri[0])),
 ];
 for edge in edges {
 *edge_count.entry(edge).or_insert(0) += 1;
 }
 }

 let mut boundary_nodes = std::collections::HashSet::new();
 for (edge, count) in edge_count {
 if count == 1 {
 boundary_nodes.insert(edge.0);
 boundary_nodes.insert(edge.1);
 }
 }

 boundary_nodes
}

fn find_best_edge_collapse(
 nodes: &[DVec3],
 triangles: &[[usize; 3]],
 boundary_nodes: &std::collections::HashSet<usize>,
 max_error: f64,
) -> Option<EdgeCollapseInfo> {
 let mut best: Option<EdgeCollapseInfo> = None;

 // Unique undirected edges
 let mut edges: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
 for &tri in triangles {
 edges.insert((tri[0].min(tri[1]), tri[0].max(tri[1])));
 edges.insert((tri[1].min(tri[2]), tri[1].max(tri[2])));
 edges.insert((tri[2].min(tri[0]), tri[2].max(tri[0])));
 }

 for edge in edges {
 // Skip collapsing true boundary edges when requested
 if boundary_nodes.contains(&edge.0) && boundary_nodes.contains(&edge.1) {
 continue;
 }

 let p0 = nodes.get(edge.0)?;
 let p1 = nodes.get(edge.1)?;

 // Length-based priority
 let error = (*p1 - *p0).length();

 if error > max_error {
 continue;
 }

 let new_position = (*p0 + *p1) * 0.5;

 match &best {
 None => best = Some(EdgeCollapseInfo { edge, error, new_position }),
 Some(current) if error < current.error => {
 best = Some(EdgeCollapseInfo { edge, error, new_position })
 }
 _ => {}
 }
 }

 best
}

fn apply_edge_collapse(
 nodes: &mut Vec<DVec3>,
 normals: &mut Vec<DVec3>,
 triangles: &mut Vec<[usize; 3]>,
 edge: (usize, usize),
 new_position: DVec3,
) {
 let (v0, v1) = edge;

 // Move `v0` to the collapsed position
 if v0 < nodes.len() {
 nodes[v0] = new_position;
 }

 // Average normals at the merged vertex
 if v0 < normals.len() && v1 < normals.len() {
 normals[v0] = (normals[v0] + normals[v1]).normalize_or_zero();
 }

 // Rewire triangles: `v1` -> `v0`
 for tri in triangles.iter_mut() {
 for i in 0..3 {
 if tri[i] == v1 {
 tri[i] = v0;
 }
 }
 }

 // Drop collapsed/degenerate faces
 triangles.retain(|&tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0]);
}

