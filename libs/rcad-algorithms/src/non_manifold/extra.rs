
/// Finds coincident faces between two shells based on geometric matching.
fn find_coincident_faces(
 shell1: &Shell,
 shell2: &Shell,
 brep1: &BRep,
 brep2: &BRep,
 tolerance: f64,
) -> Vec<(usize, usize)> {
 let mut pairs = Vec::new();

 for (f1_idx, face1) in shell1.faces.iter().enumerate() {
 let normal1 = face1.normal;
 let center1 = compute_face_center(face1, brep1);

 for (f2_idx, face2) in shell2.faces.iter().enumerate() {
 let normal2 = face2.normal;
 let center2 = compute_face_center(face2, brep2);

 // Check if normals are parallel (same or opposite direction)
 let dot = normal1.dot(normal2).abs();
 if dot < 1.0 - tolerance {
 continue; // Not parallel
 }

 // Check if centers are coincident
 if (center1 - center2).length() < tolerance {
 pairs.push((f2_idx, f1_idx));
 }
 }
 }

 pairs
}

/// Computes the geometric center of a face.
fn compute_face_center(face: &Face, brep: &BRep) -> DVec3 {
 let mut center = DVec3::ZERO;
 let mut count = 0;

 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 let edge = &brep.edges[we.idx];
 if edge.start < brep.vertices.len() {
 center += brep.vertices[edge.start].point;
 count += 1;
 }
 if edge.end < brep.vertices.len() {
 center += brep.vertices[edge.end].point;
 count += 1;
 }
 }
 }

 if count > 0 {
 center / count as f64
 } else {
 DVec3::ZERO
 }
}

/// Finds a corresponding edge in face1 that matches edge_idx from face2/brep2.
fn find_corresponding_edge(
 edge_idx: usize,
 face1: &Face,
 brep1: &BRep,
 brep2: &BRep,
 tolerance: f64,
) -> Option<usize> {
 let edge2 = &brep2.edges[edge_idx];
 let start2 = brep2.vertices.get(edge2.start)?.point;
 let end2 = brep2.vertices.get(edge2.end)?.point;

 for we in &face1.outer_wire.edges {
 let edge1 = &brep1.edges[we.idx];
 let start1 = brep1.vertices.get(edge1.start)?.point;
 let end1 = brep1.vertices.get(edge1.end)?.point;

 // Check if edges are coincident (same or reversed direction)
 let same_dir = (start1 - start2).length() < tolerance && (end1 - end2).length() < tolerance;
 let rev_dir = (start1 - end2).length() < tolerance && (end1 - start2).length() < tolerance;

 if same_dir || rev_dir {
 return Some(we.idx);
 }
 }

 None
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Non-manifold Traversal Extensions for BRepGraph
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Extension trait for non-manifold traversal on BRepGraph.
pub trait NonManifoldTraversal {
 /// Returns all faces that share a non-manifold edge with the given face.
 fn non_manifold_adjacent_faces(&self, face_idx: usize) -> Vec<usize>;

 /// Returns the number of manifold regions (connected components via manifold edges only).
 fn manifold_region_count(&self) -> usize;

 /// Returns faces grouped by manifold region.
 fn manifold_regions(&self) -> Vec<Vec<usize>>;

 /// Iterator over non-manifold edges with their adjacent faces.
 fn non_manifold_edge_info(&self) -> Vec<(usize, Vec<usize>)>;
}

impl NonManifoldTraversal for BRepGraph {
 /// Returns all faces that share a non-manifold edge with the given face.
 fn non_manifold_adjacent_faces(&self, face_idx: usize) -> Vec<usize> {
 let mut result = Vec::new();
 let edges = self.face_edges(face_idx).to_vec();

 for &ei in &edges {
 if self.edge_valence(ei) > 2 {
 // This is a non-manifold edge
 for &adj_face in self.edge_adjacent_faces(ei) {
 if adj_face != face_idx {
 result.push(adj_face);
 }
 }
 }
 }

 result.sort_unstable();
 result.dedup();
 result
 }

 /// Returns the number of manifold regions (connected components via manifold edges only).
 fn manifold_region_count(&self) -> usize {
 self.manifold_regions().len()
 }

 /// Returns faces grouped by manifold region.
 fn manifold_regions(&self) -> Vec<Vec<usize>> {
 let mut visited = vec![false; self.face_count];
 let mut regions = Vec::new();

 // Get non-manifold edges
 let nm_edges: HashSet<usize> = self.multi_face_edges().into_iter().collect();

 for start in 0..self.face_count {
 if visited[start] {
 continue;
 }

 let mut region = Vec::new();
 let mut stack = vec![start];

 while let Some(fi) = stack.pop() {
 if visited[fi] {
 continue;
 }
 visited[fi] = true;
 region.push(fi);

 // Traverse only through manifold edges
 for &ei in self.face_edges(fi) {
 if nm_edges.contains(&ei) {
 continue; // Skip non-manifold edges
 }

 for &adj in self.edge_adjacent_faces(ei) {
 if !visited[adj] {
 stack.push(adj);
 }
 }
 }
 }

 if !region.is_empty() {
 regions.push(region);
 }
 }

 regions
 }

 /// Returns information about each non-manifold edge and its adjacent faces.
 fn non_manifold_edge_info(&self) -> Vec<(usize, Vec<usize>)> {
 self.multi_face_edges()
 .into_iter()
 .map(|ei| (ei, self.edge_adjacent_faces(ei).to_vec()))
 .collect()
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Tests
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

