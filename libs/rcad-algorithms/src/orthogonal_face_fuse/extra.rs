fn clip_one_coplanar_pair(
 brep: &mut BRep,
 si: usize,
 shi: usize,
 fi: usize,
 fj: usize,
 tol: f64,
) -> Option<()> {
 // Skip curved (non-plane) faces =a cylinder wall classified with an
 // axis-aligned face normal (e.g. after boolean On= ace dedup) must not
 // be treated as a planar face.  Clipping its 2D axis projection corrupts
 // face_surface_range and destroys the curved= urface sub= ace.
 for &f in &[fi, fj] {
 let flat = flat_face_index(brep, si, shi, f);
 let sidx = brep.geom.face_surface.get(flat).copied().flatten()?;
 if !matches!(brep.geom.surfaces.get(sidx)?, rcad_kernel::geom::Surface3::Plane(_)) {
 return None;
 }
 }

 // = =  Read phase (shared borrows only) = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 let (face_i, face_j, n, axes, poly_i_uv, poly_j_uv, plane_origin, normal);

 {
 let shell = &brep.solids[si].shells[shi];
 if fi >= shell.faces.len() || fj >= shell.faces.len() {
 return None;
 }
 face_i = shell.faces[fi].clone();
 face_j = shell.faces[fj].clone();

 // No holes
 if !face_i.inner_wires.is_empty() || !face_j.inner_wires.is_empty() {
 return None;
 }

 // Both must snap to  xis
 let n_i = snap_almost_axis(face_i.normal.normalize_or_zero());
 let n_j = snap_almost_axis(face_j.normal.normalize_or_zero());
 let axes_i = axis_aligned_world_plane_uv_axes(n_i)?;
 let axes_j = axis_aligned_world_plane_uv_axes(n_j)?;
 // Must be on the same infinite plane
 if axes_i != axes_j {
 // Different axis families cannot be the same oriented plane.
 // (e.g. one Z-plane and one X-plane are different no matter what.)
 return None;
 }

 let p_i = face_first_point(brep, &face_i)?;
 let p_j = face_first_point(brep, &face_j)?;
 let d_i = n_i.dot(p_i);
 let d_j = n_j.dot(p_j);
 let (n_i_c, d_i_c) = canonicalize_plane_n_d(n_i, d_i);
 let (n_j_c, d_j_c) = canonicalize_plane_n_d(n_j, d_j);
 if plane_key(n_i_c, d_i_c, tol) != plane_key(n_j_c, d_j_c, tol) {
 return None;
 }

 // UV bbox overlap with positive area
 let bi = face_axis_world_bbox(brep, &face_i, n_i)?;
 let bj = face_axis_world_bbox(brep, &face_j, n_j)?;
 let gap = (tol * 1e2).max(TOLERANCE_MESH_LEGACY);
 if !rects_2d_bbox_positive_area_overlap(bi, bj, gap) {
 return None;
 }

 // Skip strict bbox subset =already handled by remove_axis_coplanar_redundant_child_faces.
 let scale = (bi.1 - bi.0)
 .abs()
 .max(bi.3 - bi.2)
 .abs()
 .max(bj.1 - bj.0)
 .abs()
 .max(bj.3 - bj.2)
 .abs()
 .max(1.0);
 let eps = (TOLERANCE_COORD_SUB * scale).max(tol * TOLERANCE_TOL_SCALE_MICRO);
 let subset = |a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)| -> bool {
 a.0 >= b.0 - eps && a.1 <= b.1 + eps && a.2 >= b.2 - eps && a.3 <= b.3 + eps
 };
 let s_ij = subset(bi, bj) && !subset(bj, bi);
 let s_ji = subset(bj, bi) && !subset(bi, bj);
 if s_ij || s_ji {
 return None;
 }
 if subset(bi, bj) && subset(bj, bi) {
 // Equal bboxes =also handled by the subset pass.
 return None;
 }

 // Project both face boundaries to world-axis UV
 let [i_axis, j_axis] = axes_i;
 poly_i_uv = face_outer_points(brep, &face_i)
 .iter()
 .map(|p| [p[i_axis], p[j_axis]])
 .collect::<Vec<_>>();
 poly_j_uv = face_outer_points(brep, &face_j)
 .iter()
 .map(|p| [p[i_axis], p[j_axis]])
 .collect::<Vec<_>>();

 if poly_i_uv.len() < 3 || poly_j_uv.len() < 3 {
 return None;
 }

 n = n_i;
 axes = axes_i;
 plane_origin = p_i;
 normal = face_i.normal;
 }

 // = =  Compute 2D polygon intersection = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // SH expects the clip polygon to be CCW (left-of-edge = inside).  Faces with negative
 // normals (e.g. = ) project to CW in world-axis UV, so we normalise the clip polygon.
 let poly_j_uv = ensure_ccw(&poly_j_uv);
 let overlap =
 crate::inttools::coplanar::sutherland_hodgman_clip(&poly_i_uv, &poly_j_uv);
 if overlap.len() < 3 {
 return None;
 }

 // Minimum area guard
 {
 let mut area2 = 0.0_f64;
 for k in 0..overlap.len() {
 let (x0, y0) = (overlap[k][0], overlap[k][1]);
 let (x1, y1) = (overlap[(k + 1) % overlap.len()][0], overlap[(k + 1) % overlap.len()][1]);
 area2 += x0 * y1 - x1 * y0;
 }
 let area = 0.5 * area2.abs();
 let min_area = (tol * tol).max(TOLERANCE_FLOAT_ULTRA);
 if area < min_area {
 return None;
 }
 }

 // = =  Write phase: create face from overlap polygon = = = = = = = = = = = = = = = = = = = = = 
 let [_i_ax, _j_ax] = axes;

 // Convert to (f64, f64) rings for add_vertices_for_rings_with_eval
 let overlap_uv: Vec<(f64, f64)> = overlap.iter().map(|&c| (c[0], c[1])).collect();
 let rings = vec![overlap_uv];

 let ring_vertices = add_vertices_for_rings_with_eval(brep, &rings, |u, v| {
 point_from_axis_plane_world_uv(n, plane_origin, u, v)
 }, tol);

 if ring_vertices.is_empty() || ring_vertices[0].len() < 3 {
 return None;
 }

 // Build edges
 let mut edge_pairs: Vec<(usize, usize)> = Vec::new();
 for rv in &ring_vertices {
 let nv = rv.len();
 for k in 0..nv {
 edge_pairs.push((rv[k], rv[(k + 1) % nv]));
 }
 }

 let base_ei = brep.edges.len();
 push_new_edges(brep, edge_pairs);

 // Outer wire
 let n0 = ring_vertices[0].len();
 let outer_wire = Wire {
 edges: (0..n0).map(|k| WireEdge::fwd(base_ei + k)).collect(),
 };

 let merged_face = Face {
 outer_wire,
 inner_wires: vec![],
 normal,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 // Plane surface
 let surf_idx = {
 let p0 = face_first_point(brep, &merged_face).unwrap_or(plane_origin);
 let idx = brep.geom.surfaces.len();
 brep.geom
 .surfaces
 .push(Surface3::Plane(Plane { origin: p0, normal }));
 idx
 };

 // Replace both faces with the overlap face
 let remove_fis = [fi, fj];
 let flat_indices: Vec<usize> = remove_fis
 .iter()
 .map(|&f| flat_face_index(brep, si, shi, f))
 .collect();
 replace_shell_faces_and_geom(brep, si, shi, &remove_fis, merged_face, surf_idx, &flat_indices);

 Some(())
}

