#[cfg(test)]
mod orth_union_tests {
 use crate::orthogonal_face_fuse::ring_corner_count_after_collinear_removal;
 use crate::orthogonal_face_fuse::rects_2d_bbox_positive_area_overlap;
 use crate::orthogonal_face_fuse::union_rects_to_rings_grid;
 use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_MESH_LEGACY, TOLERANCE_RETRY_LADDER_COARSE};

 #[test]
 fn union_rects_three_adjacent_strips_forms_outer_ring() {
 let rects = [
 (0.0, 5.0, 0.0, 10.0),
 (5.0, 10.0, 0.0, 10.0),
 (10.0, 15.0, 0.0, 10.0),
 ];
 let rings = union_rects_to_rings_grid(&rects, TOLERANCE_ABS);
 assert!(rings.is_some(), "expected grid union to succeed for three strips");
 let rings = rings.unwrap();
 assert!(!rings.is_empty(), "expected at least one ring");
 }

 /// Same bucket key as disjoint islands on one plane: no 2D area overlap in UV.
 #[test]
 fn bbox_positive_area_overlap_distinguishes_disjoint_corner_edge() {
 let gap = TOLERANCE_RETRY_LADDER_COARSE;
 let a = (0.0, 1.0, 0.0, 1.0);
 let b_corner = (2.0, 3.0, 2.0, 3.0);
 assert!(!rects_2d_bbox_positive_area_overlap(a, b_corner, gap));
 let b_edge = (1.0, 2.0, 0.0, 1.0);
 assert!(!rects_2d_bbox_positive_area_overlap(a, b_edge, gap));
 let c_overlap = (0.5, 1.5, 0.0, 1.0);
 assert!(rects_2d_bbox_positive_area_overlap(a, c_overlap, gap));
 }

 /// L-shaped outline keeps >4 corners; a 3 ? rectangle of samples collapses to 4 corners.
 #[test]
 fn ring_collinear_simplify_rect_vs_l() {
 let tol = TOLERANCE_MESH_LEGACY;
 let l_ring: Vec<(f64, f64)> = vec![
 (0.0, 0.0),
 (1.0, 0.0),
 (1.0, 1.0),
 (2.0, 1.0),
 (2.0, 2.0),
 (0.0, 2.0),
 ];
 assert!(ring_corner_count_after_collinear_removal(&l_ring, tol) >= 5);

 let rect_dense: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)];
 assert_eq!(ring_corner_count_after_collinear_removal(&rect_dense, tol), 4);
 }
}

/// OCCT `bcommon_simple/G1` intersection: document axis-UV bbox overlap vs strict containment between coplanar  xis faces.
///
/// Run with `cargo test -p rcad-algorithms g1_intersection_axis_bbox_relationship_probe -- --nocapture` to print counts.
/// Intended for diagnosing +1 `checkprops -s` gaps when duplicate caps do **not** satisfy strict bbox subset.
#[cfg(test)]
mod bcommon_g1_bbox_probe_tests {
 use glam::{DAffine3, DVec3};
 use rcad_kernel::BRep;
 use rcad_modeling::make_box_brep;

 use crate::boolean_op;
 use crate::tolerance::{
 TOLERANCE_ABS, TOLERANCE_COORD_SUB, TOLERANCE_MESH_LEGACY, TOLERANCE_TOL_SCALE_MICRO,
 };
 use crate::BooleanOpType;

 use crate::orthogonal_face_fuse::{
 axis_aligned_world_plane_uv_axes, canonicalize_plane_n_d, face_axis_world_bbox, face_first_point,
 plane_key, rects_2d_bbox_positive_area_overlap, snap_almost_axis,
 };

 fn g1_operands() -> (BRep, BRep) {
 let ba = make_box_brep(
 DVec3::new(0.0, 0.0, 0.0),
 DVec3::new(1.0, 0.0, 0.0).normalize(),
 DVec3::new(0.0, 1.0, 0.0).normalize(),
 1.0,
 1.0,
 1.0,
 )
 .expect("ba");
 let bb = make_box_brep(
 DVec3::new(0.0, 0.7071067811865476, 0.0),
 DVec3::new(0.0, 0.0, 1.0).normalize(),
 DVec3::new(0.0, -1.0, 0.0).normalize(),
 1.0,
 0.7071067811865476,
 1.4142135623730951,
 )
 .expect("bb");
 let bb = {
 let mut shape = bb;
 let pivot = DVec3::new(0.0, 0.0, 0.0);
 let axis = DVec3::new(0.0, 0.0, 1.0).normalize_or(DVec3::Z);
 let rot = DAffine3::from_axis_angle(axis, (45.0_f64).to_radians());
 let xf = DAffine3::from_translation(pivot) * rot * DAffine3::from_translation(-pivot);
 shape.apply_transform(xf);
 shape
 };
 (ba, bb)
 }

 #[test]
 fn g1_intersection_axis_bbox_relationship_probe() {
 let (ba, bb) = g1_operands();
 let brep = boolean_op(BooleanOpType::Intersection, &ba, &bb).expect("g1 intersection");
 let t = TOLERANCE_ABS;
 let mut strict_ij = 0usize;
 let mut strict_ji = 0usize;
 let mut overlap_only = 0usize;
 let gap = (t * 1e2).max(TOLERANCE_MESH_LEGACY);

 let subset_containment = |a: (f64, f64, f64, f64), b: (f64, f64, f64, f64), scale: f64| -> bool {
 let eps = (TOLERANCE_COORD_SUB * scale).max(t * TOLERANCE_TOL_SCALE_MICRO);
 let (au0, au1, av0, av1) = a;
 let (bu0, bu1, bv0, bv1) = b;
 au0 >= bu0 - eps && au1 <= bu1 + eps && av0 >= bv0 - eps && av1 <= bv1 + eps
 };

 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 let shell = &brep.solids[si].shells[shi];
 let n = shell.faces.len();
 for fi in 0..n {
 for fj in (fi + 1)..n {
 let fa = &shell.faces[fi];
 let fb = &shell.faces[fj];
 if !fa.inner_wires.is_empty() || !fb.inner_wires.is_empty() {
 continue;
 }
 let n_i = snap_almost_axis(fa.normal.normalize_or_zero());
 let n_j = snap_almost_axis(fb.normal.normalize_or_zero());
 if axis_aligned_world_plane_uv_axes(n_i).is_none()
 || axis_aligned_world_plane_uv_axes(n_j).is_none()
 {
 continue;
 }
 let Some(p_i) = face_first_point(&brep, fa) else {
 continue;
 };
 let Some(p_j) = face_first_point(&brep, fb) else {
 continue;
 };
 let d_i = n_i.dot(p_i);
 let d_j = n_j.dot(p_j);
 let (n_i_c, d_i_c) = canonicalize_plane_n_d(n_i, d_i);
 let (n_j_c, d_j_c) = canonicalize_plane_n_d(n_j, d_j);
 let key_i = plane_key(n_i_c, d_i_c, t);
 let key_j = plane_key(n_j_c, d_j_c, t);
 if key_i != key_j {
 continue;
 }
 let Some(bi) = face_axis_world_bbox(&brep, fa, n_i) else {
 continue;
 };
 let Some(bj) = face_axis_world_bbox(&brep, fb, n_j) else {
 continue;
 };
 let scale = (bi.1 - bi.0)
 .abs()
 .max((bi.3 - bi.2).abs())
 .max((bj.1 - bj.0).abs())
 .max((bj.3 - bj.2).abs())
 .max(1.0);

 let s_ij =
 subset_containment(bi, bj, scale) && !subset_containment(bj, bi, scale);
 let s_ji =
 subset_containment(bj, bi, scale) && !subset_containment(bi, bj, scale);
 if s_ij {
 strict_ij += 1;
 }
 if s_ji {
 strict_ji += 1;
 }
 if rects_2d_bbox_positive_area_overlap(bi, bj, gap)
 && !(subset_containment(bi, bj, scale) && subset_containment(bj, bi, scale))
 && !s_ij
 && !s_ji
 {
 overlap_only += 1;
 }
 }
 }
 }
 }

 eprintln!(
 "G1 intersection axis-plane face pairs (same shell): strict-subset one-way ij={strict_ij} ji={strict_ji}; overlap-not-mutual-subset={overlap_only}"
 );
 }
}

 
