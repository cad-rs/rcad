// NOTE: this module uses crate::* to access items re-exported at the crate root.
use crate::*;
use rcad_kernel::BRep;
use rcad_kernel::topods;
use super::booleans::{merge_pairwise_model_tol_into_boolean_options, run_make_connected_for_boolean_output, tune_boolean_options_for_retry_class, boolean_retry_followup_attempts};
use super::face_count_of;


/// Merge identical surface geometries (same plane, same cylinder etc.) into
/// a single surface entry.  The PaveFiller creates separate GeomStore entries
/// for each sub-face even when they share the same geometric surface.
fn deduplicate_surfaces(mut brep: BRep) -> BRep {
 use rcad_kernel::geom::Surface3;
 let n = brep.geom.surfaces.len();
 if n < 2 { return brep; }
 let ang_tol = 1e-6;  // TOLERANCE_ANG_HEURISTIC_RAD
 let lin_tol = crate::tolerance::TOLERANCE_PLANE_DIST_RELAX;  // 5e-6 — PaveFiller numerical noise exceeds 1e-6

 // Compute a canonical index for each surface.
 let mut canon: Vec<usize> = (0..n).collect();
 for i in 0..n {
 if canon[i] != i { continue; }  // already mapped
 for j in (i + 1)..n {
 let same = match (&brep.geom.surfaces[i], &brep.geom.surfaces[j]) {
 (Surface3::Plane(p1), Surface3::Plane(p2)) => {
 let cross = p1.normal.cross(p2.normal).length();
 cross <= ang_tol
 && (p2.origin - p1.origin).dot(p1.normal).abs() <= lin_tol
 }
 (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
 (c1.radius - c2.radius).abs() <= lin_tol
 && c1.axis.cross(c2.axis).length() <= ang_tol
 && (c2.origin - c1.origin).cross(c1.axis).length() <= lin_tol
 }
 (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
 (s1.radius - s2.radius).abs() <= lin_tol
 && (s1.center - s2.center).length() <= lin_tol
 }
 (Surface3::Cone(c1), Surface3::Cone(c2)) => {
 (c1.radius - c2.radius).abs() <= lin_tol
 && (c1.half_angle_rad - c2.half_angle_rad).abs() <= ang_tol
 && c1.axis.cross(c2.axis).length() <= ang_tol
 && (c1.apex - c2.apex).length() <= lin_tol
 }
 _ => false,  // different types or BSpline — keep separate
 };
 if same { canon[j] = i; }
 }
 }
 // Count unique surfaces.
 let mut unique: Vec<usize> = Vec::new();
 let mut old_to_new: Vec<usize> = vec![0; n];
 for i in 0..n {
 if canon[i] == i {
 old_to_new[i] = unique.len();
 unique.push(i);
 }
 }
 for i in 0..n {
 old_to_new[i] = old_to_new[canon[i]];
 }
 // Remap face_surface references.
 for s in &mut brep.solids {
 for sh in &mut s.shells {
 for _ in &sh.faces {
 // face_surface is indexed by flat face index, which is complex
 // to iterate here.  We use a different approach: rebuild surfaces.
 }
 }
 }
 // Actually, iterate via flat_face_index.
 // We need to access face_surface by flat index across all solids/shells.
 // Simpler approach: rebuild the surfaces array.
 let new_surfaces: Vec<Surface3> = unique.iter().map(|&i| brep.geom.surfaces[i].clone()).collect();
 // Now remap face_surface: for each face, find its surface in the new array.
 let mut fi = 0usize;
 for s in &mut brep.solids {
 for sh in &mut s.shells {
 for _ in &sh.faces {
 if let Some(Some(old_si)) = brep.geom.face_surface.get(fi).copied() {
 let new_si = old_to_new[old_si];
 brep.geom.face_surface[fi] = Some(new_si);
 }
 fi += 1;
 }
 }
 }
 brep.geom.surfaces = new_surfaces;
 brep
}

/// Post-process a boolean operation's BRep result: merge coplanar faces on
/// the same plane, share edges between adjacent faces, and detect holes
/// (inner wires) for faces with missing interior regions.
fn count_topo(brep: &BRep, label: &str) {
 use rcad_kernel::geom::Surface3;
 let nf = brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
 let ne = brep.edges.len();
 let sphere_fe = brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).filter(|f| {
 f.surface_idx.and_then(|si| brep.geom.surfaces.get(si)).is_some_and(|s| matches!(s, Surface3::Sphere(_)))
 }).map(|f| f.outer_wire.edges.len()).next().unwrap_or(0);
 eprintln!("[CNT] {}: {} faces, {} edges, sphere_fe={}", label, nf, ne, sphere_fe);
}

fn optimize_boolean_topology(mut brep: BRep) -> BRep {
 count_topo(&brep, "topo-enter");
 if brep.vertices.len() < 4 { return brep; }
 // Allow fast step-only mode: skip all topology passes.
 if std::env::var("RCAD_SKIP_TOPOLOGY").is_ok() { return brep; }
 use rcad_kernel::topology::{Face, Wire, WireEdge};
 use rcad_kernel::{Edge, Vertex};
 use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};

 let tol = tolerance::TOLERANCE_ABS.max(1e-8);
 // Pass 1: removed orthogonal_face_fuse (self-created, no OCCT equivalent).
 //          OCCT-aligned FillSameDomainFaces below handles same-domain merging.
 let m1 = deduplicate_edges(brep);
 // Surface deduplication: merge identical surface geometries (same plane,
 // same cylinder, etc.) into a single surface entry.  The PaveFiller creates
 // separate entries for each sub-face even when they share the same geometry.
 let m1 = deduplicate_surfaces(m1);
 // ✅ OCCT : FillSameDomainFaces — edge-set  (BOPAlgo_Builder_2.cxx L636-L796)
 let (m2, _) = crate::occt_fill_same_domain_faces(&m1);
 count_topo(&m2, "topo-pass2");
 let mut brep = m2;

 // Pass 3: detect remaining coplanar groups with hole patterns and merge
 // sub-faces into a single face with outer wire + inner wire, reusing
 // existing edges so the resulting shell has shared edge topology.
 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 let nf = brep.solids[si].shells[shi].faces.len();
 if nf < 2 { continue; }
 // Group faces by plane
 let mut pg: Vec<Vec<usize>> = Vec::new();
 let mut pk: Vec<(f64,f64,f64,f64)> = Vec::new();
 for fi in 0..nf {
 let face = &brep.solids[si].shells[shi].faces[fi];
 let n = face.normal;
 if !face.inner_wires.is_empty() { continue; }
 let pd = face.outer_wire.edges.first().and_then(|we|
 brep.edges.get(we.idx).and_then(|e|
 brep.vertices.get(e.start).map(|v| n.dot(v.point))
 )
 ).unwrap_or(0.0);
 let key = (n.x, n.y, n.z, pd);
 if let Some(pos) = pk.iter().position(|k| {
 (k.0-key.0).abs()<1e-8 && (k.1-key.1).abs()<1e-8
 && (k.2-key.2).abs()<1e-8 && (k.3-key.3).abs()<1e-8
 }) { pg[pos].push(fi); }
 else { pk.push(key); pg.push(vec![fi]); }
 }
 let mut to_remove: Vec<usize> = Vec::new();
 let mut new_faces: Vec<Face> = Vec::new();
 for group in &pg {
 if group.len() < 2 { continue; }
 let mut pts: Vec<glam::DVec3> = Vec::new();
 for &fi in group { for we in &brep.solids[si].shells[shi].faces[fi].outer_wire.edges {
 if let Some(e) = brep.edges.get(we.idx) {
 if let Some(v) = brep.vertices.get(e.start) { pts.push(v.point); }
 if let Some(v) = brep.vertices.get(e.end) { pts.push(v.point); }
 }
 }}
 let omin = pts.iter().copied().fold(glam::DVec3::splat(f64::MAX), glam::DVec3::min);
 let omax = pts.iter().copied().fold(glam::DVec3::splat(f64::NEG_INFINITY), glam::DVec3::max);
 let (u_idx, v_idx): (usize,usize) = if (omin.x-omax.x).abs()<1e-8 {(1,2)}
 else if (omin.y-omax.y).abs()<1e-8 {(0,2)} else {(0,1)};
 let w_idx = 3 - u_idx - v_idx;
 let u_min = pts.iter().map(|p|p[u_idx]).fold(f64::MAX,f64::min);
 let u_max = pts.iter().map(|p|p[u_idx]).fold(f64::NEG_INFINITY,f64::max);
 let v_min = pts.iter().map(|p|p[v_idx]).fold(f64::MAX,f64::min);
 let v_max = pts.iter().map(|p|p[v_idx]).fold(f64::NEG_INFINITY,f64::max);
 let mut h_umin = f64::MAX; let mut h_umax = f64::NEG_INFINITY;
 let mut h_vmin = f64::MAX; let mut h_vmax = f64::NEG_INFINITY;
 for &p in &pts {
 let on_outer = (p[u_idx]-u_min).abs()<1e-8||(p[u_idx]-u_max).abs()<1e-8
 ||(p[v_idx]-v_min).abs()<1e-8||(p[v_idx]-v_max).abs()<1e-8;
 if !on_outer { h_umin=h_umin.min(p[u_idx]); h_umax=h_umax.max(p[u_idx]);
 h_vmin=h_vmin.min(p[v_idx]); h_vmax=h_vmax.max(p[v_idx]); }
 }
 if h_umin.is_infinite() { continue; }
 let n = brep.solids[si].shells[shi].faces[group[0]].normal;
 let w_val = omin[w_idx];
 let mk_pt = |u: f64, v: f64| -> glam::DVec3 {
 let mut a=[0.0;3]; a[w_idx]=w_val; a[u_idx]=u; a[v_idx]=v; glam::DVec3::from_array(a)
 };
 // Find outer perimeter edges from non-removed faces (adjacent
 // faces that have clean corner-to-corner topology)
 let o_pts = [mk_pt(u_min,v_min), mk_pt(u_max,v_min), mk_pt(u_max,v_max), mk_pt(u_min,v_max)];
 let mut outer_we: Vec<WireEdge> = Vec::new();
 for k in 0..4 {
 let a = o_pts[k]; let b = o_pts[(k+1)%4];
 let mut found = false;
 for (fi, face) in brep.solids[si].shells[shi].faces.iter().enumerate() {
 if to_remove.contains(&fi) { continue; }
 for we in &face.outer_wire.edges {
 if let Some(e) = brep.edges.get(we.idx) {
 let sa = brep.vertices.get(e.start).map(|v|v.point);
 let sb = brep.vertices.get(e.end).map(|v|v.point);
 if let (Some(sa), Some(sb)) = (sa, sb) {
 if (sa-a).length()<1e-8 && (sb-b).length()<1e-8 {
 outer_we.push(WireEdge{idx:we.idx, forward:true});
 found = true; break;
 }
 if (sb-a).length()<1e-8 && (sa-b).length()<1e-8 {
 outer_we.push(WireEdge{idx:we.idx, forward:false});
 found = true; break;
 }
 }
 }
 }
 if found { break; }
 }
 if !found { break; }
 }
 if outer_we.len() != 4 { continue; }
 // Find inner perimeter edges from non-removed faces (channel walls)
 let i_pts = [mk_pt(h_umin,h_vmin), mk_pt(h_umax,h_vmin), mk_pt(h_umax,h_vmax), mk_pt(h_umin,h_vmax)];
 let mut inner_we: Vec<WireEdge> = Vec::new();
 for k in 0..4 {
 let a = i_pts[k]; let b = i_pts[(k+1)%4];
 let mut found = false;
 for (fi, face) in brep.solids[si].shells[shi].faces.iter().enumerate() {
 if to_remove.contains(&fi) { continue; }
 for we in &face.outer_wire.edges {
 if let Some(e) = brep.edges.get(we.idx) {
 let sa = brep.vertices.get(e.start).map(|v|v.point);
 let sb = brep.vertices.get(e.end).map(|v|v.point);
 if let (Some(sa), Some(sb)) = (sa, sb) {
 if (sa-a).length()<1e-8 && (sb-b).length()<1e-8 {
 inner_we.push(WireEdge{idx:we.idx, forward:false});
 found = true; break;
 }
 if (sb-a).length()<1e-8 && (sa-b).length()<1e-8 {
 inner_we.push(WireEdge{idx:we.idx, forward:true});
 found = true; break;
 }
 }
 }
 }
 if found { break; }
 }
 if !found { break; }
 }
 if inner_we.len() != 4 { continue; }
 // Create merged face with outer wire + inner wire
 for &fi in group { to_remove.push(fi); }
 let surf_idx = brep.geom.surfaces.len();
 brep.geom.surfaces.push(Surface3::Plane(Plane{origin: o_pts[0], normal: n}));
 new_faces.push(Face {
 outer_wire: Wire{edges: outer_we},
 inner_wires: vec![Wire{edges: inner_we}],
 normal: n, triangles: vec![], sample_point: None, mesh_dirty: true,
 surface_idx: None,
 });
 }
 if new_faces.is_empty() { continue; }
 let mut kept: Vec<Face> = Vec::new();
 for (fi, face) in brep.solids[si].shells[shi].faces.iter().enumerate() {
 if !to_remove.contains(&fi) { kept.push(face.clone()); }
 }
 kept.extend(new_faces);
 // Rebuild face_surface
 let mut nfs: Vec<Option<usize>> = Vec::with_capacity(kept.len());
 let mut nfsr: Vec<Option<[f64;4]>> = Vec::with_capacity(kept.len());
 for face in &kept {
 let origin = face.outer_wire.edges.first().and_then(|we|
 brep.edges.get(we.idx).and_then(|e|
 brep.vertices.get(e.start).map(|v| v.point)
 )
 ).unwrap_or(glam::DVec3::ZERO);
 let norm = face.normal;
 let si2 = brep.geom.surfaces.iter().position(|s| {
 if let Surface3::Plane(p) = s {
 (p.normal-norm).length()<1e-8 && (p.origin-origin).length()<1e-8
 } else { false }
 });
 match si2 {
 Some(idx) => { nfs.push(Some(idx)); nfsr.push(None); }
 None => {
 let idx = brep.geom.surfaces.len();
 brep.geom.surfaces.push(Surface3::Plane(Plane{origin, normal: norm}));
 nfs.push(Some(idx)); nfsr.push(None);
 }
 }
 }
 brep.geom.face_surface = nfs;
 brep.geom.face_surface_range = nfsr;
 brep.solids[si].shells[shi].faces = kept;
 }
 }

 // Pass 4: General cleanup for planar-heavy models (internal faces,
 // duplicates, degenerate faces, vertex merge, edge sewing).
 // Skip for curved-surface results (most faces are Cylinder/Cone/Sphere)
 // where the O(n^2) detection is expensive and rarely beneficial.
// let n_planar = brep.geom.surfaces.iter().filter(|s| matches!(s, Surface3::Plane(_))).count();
// let n_curved = brep.geom.surfaces.len().saturating_sub(n_planar);
// if n_planar > n_curved && brep.solids.iter().any(|s| s.shells.iter().any(|sh| sh.faces.len() > 4)) {
// // Pass 4 disabled — cleanup_boolean_result can incorrectly remove
// // faces from concave-extruded shapes (H1/H2), breaking the solid.
// // let (cleaned, _report) = crate::brep_repair::cleanup_boolean_result(&brep, tol);
// // brep = cleaned;
// }
// 
// // Pass 5: Advanced simplification for planar-heavy results.
// if n_planar > n_curved {
// let s_opts = crate::SimplifyOptions {
// remove_small_edges: true,
// ..Default::default()
// };
// let (simplified, _srep) = crate::simplify_brep_post_ops(&brep, s_opts);
// brep = simplified;
 if std::env::var("RCAD_DEBUG_BOX").is_ok() {
 let nv = brep.vertices.len();
 let ne = brep.edges.len();
 let nf = brep.solids.get(0).and_then(|s| s.shells.get(0)).map(|sh| sh.faces.len()).unwrap_or(0);
 eprintln!("[POST_OPT] V={} E={} F={}", nv, ne, nf);
 }
 brep
}

/// Like [`boolean_op`] but with conservative auto-retry for numerical-instability cases.
///
/// First tries the standard [`boolean_op_pave_fill_build`] path (identical to
/// [`boolean_op`]'s first attempt).  On failure, delegates to
/// [`boolean::boolean_op_with_retry_policy`] with [`RetryPolicy::conservative`]
/// and default [`BooleanOptions`] to escalate fuzzy tolerance, glue mode, and
/// make-connected passes.
pub fn boolean_op_with_retry(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
) -> Result<BRep, BooleanError> {
 // First attempt: standard path including fast-paths (containment, box-box, etc.).
 let brep = if let Ok(mut t) = boolean_op(op, &a.to_topods(), &b.to_topods()) {
 // ✅ OCCT-aligned: correct_tolerances on topods::BRep directly.
 rcad_kernel::tolerance::correct_tolerances_topods(&mut t, 23, 0.05);
 rcad_kernel::tolerance::correct_shape_tolerances_topods(&mut t);
 rcad_kernel::BRep::from_topods(&t)
 } else {
 // Fallback: retry with escalating tolerance.
 let mut brep = boolean::boolean_op_with_retry_policy(
 op, a, b, &RetryPolicy::conservative(), BooleanOptions::default(),
 )
 .map(|(brep, _report)| brep)?;
 // ✅ OCCT-aligned: correct_tolerances on old BRep fallback.
 rcad_kernel::tolerance::correct_tolerances(&mut brep, 23, 0.05);
 rcad_kernel::tolerance::correct_shape_tolerances(&mut brep);
 brep
 };
 Ok(brep)
}

/// Perform a boolean operation with advanced execution options and report.
pub fn boolean_op_with_options(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
 mut options: BooleanOptions,
) -> Result<(BRep, BooleanExecutionReport), BooleanError> {
 merge_pairwise_model_tol_into_boolean_options(&mut options, a, b);

 let input_faces_a = face_count_of(a);
 let input_faces_b = face_count_of(b);
 let used_bvh = options.use_bvh && has_faces(a) && has_faces(b);

 let (mut out, mut report, history_opt) = if options.include_history {
 let (result, history) = if options.use_bvh {
 if options.fuzzy_tol <= 0.0 && !options.use_glue {
 let (r, h) = boolean_op_with_history(op, a, b)?; (r.to_topods(), h)
 } else {
 let a_t = a.to_topods();
 let b_t = b.to_topods();
 let ds_tol = if options.fuzzy_tol > 0.0 { options.fuzzy_tol } else { TOLERANCE_ABS };
 let mut ds = bopds::ds::DS::new_from_topods(&a_t, &b_t, ds_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.configure_fuzzy(options.fuzzy_tol);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let (t, h) = builder.build_with_history()?;
 (t, h)
 }
 } else {
 let a_t = a.to_topods();
 let b_t = b.to_topods();
 let ds_tol = if options.fuzzy_tol > 0.0 { options.fuzzy_tol } else { TOLERANCE_ABS };
 let mut ds = bopds::ds::DS::new_from_topods(&a_t, &b_t, ds_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let mut filler = pave_filler::PaveFiller::new(&mut ds);
 filler.brep = Some(&mut brep);
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.configure_fuzzy(options.fuzzy_tol);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let (t, h) = builder.build_with_history()?;
 (t, h)
 };
 (
 result,
 BooleanExecutionReport {
 input_faces_a,
 input_faces_b,
 used_bvh,
 ..BooleanExecutionReport::default()
 },
 Some(history),
 )
 } else {
 let result = if options.use_bvh {
 if options.fuzzy_tol > 0.0 || options.use_glue {
 let a_t = a.to_topods();
 let b_t = b.to_topods();
 let mut ds = bopds::ds::DS::new_from_topods(&a_t, &b_t, options.fuzzy_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let r = builder.build()?;
 boolean_postprocess_pave_result_topods(op, a, b, r)?
 } else {
 boolean_op(op, &a.to_topods(), &b.to_topods())?
 }
 } else {
 let a_t = a.to_topods();
 let b_t = b.to_topods();
 let ds_tol = if options.fuzzy_tol > 0.0 { options.fuzzy_tol } else { TOLERANCE_ABS };
 let mut ds = bopds::ds::DS::new_from_topods(&a_t, &b_t, ds_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let mut filler = pave_filler::PaveFiller::new(&mut ds);
 filler.brep = Some(&mut brep);
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let r = builder.build()?;
 boolean_postprocess_pave_result_topods(op, a, b, r)?
 };
 (
 result,
 BooleanExecutionReport {
 input_faces_a,
 input_faces_b,
 used_bvh,
 ..BooleanExecutionReport::default()
 },
 None,
 )
 };

 if options.run_healing {
 let mut healing_options = options.healing;
 // If boolean make-connected is enabled, allow healing to use the same
 // connectivity rebuild policy when repair passes stall.
 if options.run_make_connected {
 healing_options.make_connected_prepass_mode = MakeConnectedPrepassMode::IssueDriven;
 healing_options.run_make_connected_on_stall = true;
 healing_options.make_connected_tolerance = options.make_connected_tolerance;
 healing_options.make_connected_max_passes = options.make_connected_max_passes;
 healing_options.make_connected_tolerance_growth =
 options.make_connected_tolerance_growth;
 healing_options.make_connected_tolerance_cap = options.make_connected_tolerance_cap;
 }
 let (healed, heal_report) = analyze_and_heal(&out, healing_options);
 out = healed;
 report.healed = true;
 report.healing_report = Some(heal_report);
 }

 if options.run_make_connected {
 let old_for_mc = rcad_kernel::BRep::from_topods(&out);
 let (connected, connected_report) = run_make_connected_for_boolean_output(
 &old_for_mc,
 history_opt.as_ref(),
 &options,
 &mut report,
 );
 out = connected.to_topods();
 report.made_connected = true;
 report.make_connected_report = Some(connected_report);
 }

 if options.run_simplify {
 let (simplified, simp_report) = simplify_brep_post_ops(&out, options.simplify);
 out = simplified;
 report.simplified = true;
 report.simplify_report = Some(simp_report);
 }

 if options.run_propagate_geom_tolerances {
 let floor = resolved_boolean_fuzzy_tol_for_ds(options.fuzzy_tol);
 let old_out = rcad_kernel::BRep::from_topods(&out);
 out = propagate_tolerances(&old_out, floor, ToleranceFlowDirection::BottomUp).to_topods();
 report.propagated_geom_tolerances = true;
 }

 report.output_faces = face_count_of(&rcad_kernel::BRep::from_topods(&out));
 report.configured_fuzzy_tol = options.fuzzy_tol;
 report.effective_fuzzy_tol = resolved_boolean_fuzzy_tol_for_ds(options.fuzzy_tol);
 report.boolean_history = history_opt.as_ref().cloned();
 if let Some(history) = history_opt {
 report.history_faces = history.len();
 report.history_edges = history.edge_origins.len();
 report.history_vertices = history.vertex_origins.len();
 report.history_shells = history.shell_origins.len();
 report.history_solids = history.solid_origins.len();
 report.persistent_face_labels = persistent_face_labels_from_history(&history);
 report.persistent_edge_labels = persistent_edge_labels_from_history(&history);
 report.persistent_shell_labels = persistent_shell_labels_from_history(&history);
 report.persistent_solid_labels = persistent_solid_labels_from_history(&history);
 }

 Ok((rcad_kernel::BRep::from_topods(&out), report))
}

/// Robust boolean operation with automatic fuzzy-tolerance retries.
///
/// Attempts run in this order:
/// 1. `options.base.fuzzy_tol`
/// 2. each value in `options.fuzzy_retry_ladder`
///
/// The first successful attempt is returned, with retry metadata in
/// [`BooleanExecutionReport`].
pub fn boolean_op_robust(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
 options: BooleanRobustOptions,
) -> Result<(BRep, BooleanExecutionReport), BooleanError> {
 const MAX_RETRY_ESCALATION_ROUNDS: usize = 2;

 let mut pending = std::collections::VecDeque::new();
 pending.push_back((options.base.fuzzy_tol.max(0.0), None, 0usize));
 let mut tried: Vec<(f64, Option<BooleanRetryClass>, usize)> = Vec::new();
 let mut attempt_reports: Vec<BooleanRobustAttemptReport> = Vec::new();
 let mut last_err: Option<BooleanError> = None;

 while let Some((fuzzy, origin_retry_class, retry_round)) = pending.pop_front() {
 if tried.iter().any(|(v, cls, round)| {
 (*v - fuzzy).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP && *cls == origin_retry_class && *round == retry_round
 }) {
 continue;
 }
 tried.push((fuzzy, origin_retry_class, retry_round));

 let mut attempt_options = options.base;
 attempt_options.fuzzy_tol = fuzzy;
 tune_boolean_options_for_retry_class(&mut attempt_options, origin_retry_class, retry_round);
 let attempt_make_connected_scoped_enabled =
 attempt_options.run_make_connected && attempt_options.make_connected_scoped;
 let attempt_scope_seed_mode =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_seed_mode)
 } else {
 None
 };
 let attempt_scope_history_ring_depth =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_history_ring_depth)
 } else {
 None
 };
 let attempt_scope_seed_length =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_seed_length)
 } else {
 None
 };
 let attempt_scope_min_history_edges =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_min_history_edges)
 } else {
 None
 };
 match boolean_op_with_options(op, a, b, attempt_options) {
 Ok((brep, mut report)) => {
 attempt_reports.push(BooleanRobustAttemptReport {
 fuzzy_tol: fuzzy,
 success: true,
 retry_round,
 origin_retry_class,
 make_connected_scoped_enabled: attempt_make_connected_scoped_enabled,
 make_connected_scope_seed_mode: report.make_connected_scope_seed_mode,
 make_connected_scope_history_ring_depth: report
 .make_connected_scope_history_ring_depth,
 make_connected_scope_seed_length: attempt_scope_seed_length,
 make_connected_scope_min_history_edges: attempt_scope_min_history_edges,
 make_connected_scope_seed_source: report.make_connected_scope_seed_source,
 make_connected_scope_history_seed_edge_count: Some(
 report.make_connected_scope_history_seed_edge_count,
 ),
 make_connected_scope_heuristic_seed_edge_count: Some(
 report.make_connected_scope_heuristic_seed_edge_count,
 ),
 make_connected_scope_seed_vertex_count: Some(
 report.make_connected_scope_seed_vertices.len(),
 ),
 make_connected_scope_seed_edge_count: Some(
 report.make_connected_scope_seed_edges.len(),
 ),
 used_glue: attempt_options.use_glue,
 glue_tolerance: attempt_options.glue_tolerance,
 retry_class: None,
 error_message: None,
 output_faces: Some(report.output_faces),
 made_connected: report.made_connected,
 make_connected_scope_fallback_applied: report
 .make_connected_scope_fallback_applied,
 make_connected_scope_fallback_reason: report
 .make_connected_scope_fallback_reason,
 make_connected_scope_seed_edge_coverage: report
 .make_connected_scope_seed_edge_coverage,
 make_connected_scope_seed_face_coverage: report
 .make_connected_scope_seed_face_coverage,
 make_connected_scope_global_fallback_initial_tolerance: report
 .make_connected_scope_global_fallback_initial_tolerance,
 make_connected_scope_global_fallback_max_passes: report
 .make_connected_scope_global_fallback_max_passes,
 });
 report.robust_attempts = attempt_reports;
 report.retry_count = tried.len().saturating_sub(1);
 report.configured_fuzzy_tol = fuzzy;
 report.effective_fuzzy_tol = resolved_boolean_fuzzy_tol_for_ds(fuzzy);
 return Ok((brep, report));
 }
 Err(err) => {
 let retry_class = classify_boolean_retry(&err);
 attempt_reports.push(BooleanRobustAttemptReport {
 fuzzy_tol: fuzzy,
 success: false,
 retry_round,
 origin_retry_class,
 make_connected_scoped_enabled: attempt_make_connected_scoped_enabled,
 make_connected_scope_seed_mode: attempt_scope_seed_mode,
 make_connected_scope_history_ring_depth: attempt_scope_history_ring_depth,
 make_connected_scope_seed_length: attempt_scope_seed_length,
 make_connected_scope_min_history_edges: attempt_scope_min_history_edges,
 make_connected_scope_seed_source: None,
 make_connected_scope_history_seed_edge_count: None,
 make_connected_scope_heuristic_seed_edge_count: None,
 make_connected_scope_seed_vertex_count: None,
 make_connected_scope_seed_edge_count: None,
 used_glue: attempt_options.use_glue,
 glue_tolerance: attempt_options.glue_tolerance,
 retry_class: Some(retry_class),
 error_message: Some(format!("{err:?}")),
 output_faces: None,
 made_connected: false,
 make_connected_scope_fallback_applied: false,
 make_connected_scope_fallback_reason: None,
 make_connected_scope_seed_edge_coverage: None,
 make_connected_scope_seed_face_coverage: None,
 make_connected_scope_global_fallback_initial_tolerance: None,
 make_connected_scope_global_fallback_max_passes: None,
 });
 for candidate in boolean_retry_followup_attempts(
 fuzzy,
 &options.fuzzy_retry_ladder,
 &err,
 options.retry_policy,
 origin_retry_class,
 retry_round,
 MAX_RETRY_ESCALATION_ROUNDS,
 attempt_make_connected_scoped_enabled,
 ) {
 let seen = tried.iter().any(|(v, cls, round)| {
 (*v - candidate.0).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP
 && *cls == candidate.1
 && *round == candidate.2
 }) || pending.iter().any(|(v, cls, round)| {
 (*v - candidate.0).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP
 && *cls == candidate.1
 && *round == candidate.2
 });
 if !seen {
 pending.push_back(candidate);
 }
 }
 last_err = Some(err);
 }
 }
 }

 Err(last_err.unwrap_or(BooleanError::DegenerateResult))
}

/// Run post-operation simplification passes on a BRep.
pub fn simplify_brep_post_ops(brep: &topods::BRep, options: SimplifyOptions) -> (topods::BRep, SimplifyReport) {
 let old = rcad_kernel::BRep::from_topods_with_location(brep, glam::DAffine3::IDENTITY);
 let (result, report) = simplify_brep_post_ops_old(&old, options);
 (result.to_topods(), report)
}

/// Legacy: takes old BRep.
fn simplify_brep_post_ops_old(brep: &BRep, options: SimplifyOptions) -> (BRep, SimplifyReport) {
 fn closure_score(brep: &BRep) -> usize {
 let report = crate::brep_check::validate_solid_closure(brep);
 report
 .issues
 .iter()
 .map(|iss| match iss {
 crate::CheckIssue::SolidNotClosed {
 boundary_edge_count,
 ..
 } => *boundary_edge_count,
 _ => 1,
 })
 .sum()
 }

 let before = brep_check_analyze(brep);
 let mut out = brep.clone();
 let mut report = SimplifyReport {
 issues_before: before.issues.len(),
 ..SimplifyReport::default()
 };

 if options.merge_vertices {
 let (next, merged) = merge_close_vertices(&out, options.merge_tolerance);
 out = next;
 report.vertices_merged = merged;
 }
 if options.recompute_normals {
 let (next, n) = recompute_face_normals(&out);
 out = next;
 report.normals_recomputed = n;
 }
 if options.remove_degenerate_faces {
 let (next, n) = remove_degenerate_faces(&out);
 out = next;
 report.degenerate_faces_removed = n;
 }
 if options.remove_internal_faces {
 let (next, n) = remove_internal_faces(&out);
 out = next;
 report.internal_faces_removed = n;
 }
 if options.fix_wire_orientation {
 let (next, n) = fix_wire_orientation(&out, options.merge_tolerance);
 out = next;
 report.wires_fixed = n;
 }
 if options.unify_same_domain_faces {
 let cur_score = closure_score(&out);
 let (next, n) = unify_same_domain_faces(&out);
 let next_score = closure_score(&next);
 if next_score <= cur_score {
 out = next;
 report.same_domain_face_merges = n;
 }
 }
  // After same-domain unification, run same-domain unification once more to
 // absorb newly adjacent coplanar patches produced by the fuse pass.
 if options.unify_same_domain_faces {
 let cur_score = closure_score(&out);
 let (next, n) = unify_same_domain_faces(&out);
 let next_score = closure_score(&next);
 if next_score <= cur_score {
 out = next;
 report.same_domain_face_merges += n;
 }
 }

 // Kernel-level wire cleanup: collapse consecutive collinear segments so
 // post-boolean faces do not keep fragmented edge chains.
 let collinear_edge_merges = rcad_kernel::merge_collinear_edges_in_wires(
 &mut out,
 options.merge_tolerance.max(tolerance::TOLERANCE_ABS),
 );
 report.wires_fixed += collinear_edge_merges;

 if options.remove_small_edges {
 let cur_score = closure_score(&out);
 let (next, n) = remove_small_edges(&out, options.small_edge_min_length);
 let next_score = closure_score(&next);
 if next_score <= cur_score {
 out = next;
 report.small_edges_removed = n;
 }
 }

 // Final safety net: never return an open solid from simplification if it
 // can be repaired into a closed one with the standard solid fixer.
 if !crate::brep_check::validate_solid_closure(&out).is_clean() {
 let (fixed, _fix_report) =
 fix_solid(&out, options.merge_tolerance.max(tolerance::TOLERANCE_ABS));
 if crate::brep_check::validate_solid_closure(&fixed).is_clean() {
 out = fixed;
 } else {
 let (healed, _heal_report) = heal_comprehensive(&out, &HealingOptions::default());
 if crate::brep_check::validate_solid_closure(&healed).is_clean() {
 out = healed;
 }
 }
 }

 // Face merges (same-domain / orthogonal coplanar) leave `triangles` empty with
 // `mesh_dirty=true`. Callers that use `Tessellator::tessellate(&brep)` without
 // `mesh_brep` would draw only edges and show interior voids ("open box").
 if out
 .solids
 .iter()
 .flat_map(|s| s.shells.iter())
 .flat_map(|sh| sh.faces.iter())
 .any(|f| !f.mesh_is_clean())
 {
 crate::triangulate::mesh_brep(&mut out, &crate::triangulate::TessellationParams::default());
 }

 report.issues_after = brep_check_analyze(&out).issues.len();
 (out, report)
}

/// Boolean + simplification convenience pipeline.
///
/// Mirrors OCCT's `BRepAlgoAPI::SimplifyResult()` which runs
/// `BRepLib_MakeConnected` before simplification to merge coincident
/// vertices and edges, ensuring clean topological connectivity.
pub fn boolean_op_simplified(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
 options: SimplifyOptions,
) -> Result<(BRep, SimplifyReport), BooleanError> {
 let t = boolean_op_topods_simplified(op, a, b, options)?;
 Ok((rcad_kernel::BRep::from_topods(&t.0), t.1))
}

/// Same as boolean_op_simplified but returns topods::BRep.
pub fn boolean_op_topods_simplified(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
 options: SimplifyOptions,
) -> Result<(topods::BRep, SimplifyReport), BooleanError> {
 let raw = boolean_op(op, &a.to_topods(), &b.to_topods())?;
 let (connected, _mc_report) = make_connected_enhanced(
 &raw,
 tolerance::TOLERANCE_ABS,
 3,
 );
 let (simplified, report) = simplify_brep_post_ops(&connected, options);
 Ok((simplified, report))
}

/// Split `target` by one or more `tools` without boolean classification.
///
/// This is a first-stage splitter built on top of [`imprint_shape`]. It keeps
/// target material and iteratively imprints tool boundaries onto the evolving
/// target shape.
pub fn split_shape(target: &BRep, tools: &[BRep]) -> (BRep, SplitterReport) {
 split_shape_with_options(target, tools, SplitterOptions::default())
}

/// Like [`split_shape`] with advanced options.
pub fn split_shape_with_options(
 target: &BRep,
 tools: &[BRep],
 options: SplitterOptions,
) -> (BRep, SplitterReport) {
 let (result, report) = split_brep_internal_with_partial_report(target, tools, options, false);
 match result {
 Ok(brep) => (brep, report),
 Err(_) => unreachable!("unchecked splitter path should not fail"),
 }
}

/// Split `target` by tools and validate each executed step.
///
/// Returns a step-indexed error if an intermediate split result has structural
/// validity issues, excluding `NonManifoldEdge` (which can be expected for
/// split-first intermediate topology).
pub fn split_shape_checked_with_options(
 target: &BRep,
 tools: &[BRep],
 options: SplitterOptions,
) -> Result<(BRep, SplitterReport), SplitterError> {
 let (result, report) = split_brep_internal_with_partial_report(target, tools, options, true);
 result.map(|brep| (brep, report))
}

fn split_brep_internal_with_partial_report(
 target: &BRep,
 tools: &[BRep],
 options: SplitterOptions,
 validate_each_step: bool,
) -> (Result<BRep, SplitterError>, SplitterReport) {
 let mut acc = target.clone();
 let mut report = SplitterReport::default();

 for (step_index, tool) in tools.iter().enumerate() {
 let input_faces = face_count_of(&acc);
 let fuzzy = options.fuzzy_tolerance.max(0.0);
 let skipped_by_broad_phase =
 options.broad_phase_pruning && breps_farther_than_tolerance(&acc, tool, fuzzy);

 if skipped_by_broad_phase {
 report.steps.push(SplitterStepReport {
 step_index,
 input_faces,
 seam_edges: 0,
 output_faces: input_faces,
 healed: false,
 skipped_by_broad_phase: true,
 validation_issue_count: if validate_each_step { Some(0) } else { None },
 validation_first_issue: None,
 });
 continue;
 }

 let mut step = imprint_shape(&acc, tool);
 let seam_edges = step.seam_edges.len();

 if options.heal_after_each_step {
 let mut healing = options.healing;
 align_healing_options_with_boolean_operands(
 &mut healing,
 &acc,
 tool,
 options.fuzzy_tolerance,
 );
 let (healed, _) = analyze_and_heal(&step.brep.to_topods(), healing);
 step.brep = rcad_kernel::BRep::from_topods(&healed);
 }

 let mut validation_issue_count = None;
 let mut validation_first_issue = None;
 let output_faces = face_count_of(&step.brep);
 if validate_each_step {
 let validity = brep_check_analyze(&step.brep);
 let (issue_count, first_issue) =
 splitter_issues_by_level(&validity, options.validation_level);
 validation_issue_count = Some(issue_count);
 validation_first_issue = first_issue.clone();
 if issue_count > 0 {
 report.steps.push(SplitterStepReport {
 step_index,
 input_faces,
 seam_edges,
 output_faces,
 healed: options.heal_after_each_step,
 skipped_by_broad_phase: false,
 validation_issue_count,
 validation_first_issue,
 });
 return (
 Err(SplitterError::StepInvalid {
 step_index,
 issue_count,
 first_issue,
 }),
 report,
 );
 }
 }

 report.total_seam_edges += seam_edges;
 report.steps.push(SplitterStepReport {
 step_index,
 input_faces,
 seam_edges,
 output_faces,
 healed: options.heal_after_each_step,
 skipped_by_broad_phase: false,
 validation_issue_count,
 validation_first_issue,
 });

 acc = step.brep;
 }

 (Ok(acc), report)
}

fn brep_bounds(brep: &BRep) -> Option<(glam::DVec3, glam::DVec3)> {
 let mut it = brep.vertices.iter();
 let first = it.next()?.point;
 let mut min = first;
 let mut max = first;
 for v in it {
 min = min.min(v.point);
 max = max.max(v.point);
 }
 Some((min, max))
}

fn aabb_distance(
 min_a: glam::DVec3,
 max_a: glam::DVec3,
 min_b: glam::DVec3,
 max_b: glam::DVec3,
) -> f64 {
 let dx = if max_a.x < min_b.x {
 min_b.x - max_a.x
 } else if max_b.x < min_a.x {
 min_a.x - max_b.x
 } else {
 0.0
 };
 let dy = if max_a.y < min_b.y {
 min_b.y - max_a.y
 } else if max_b.y < min_a.y {
 min_a.y - max_b.y
 } else {
 0.0
 };
 let dz = if max_a.z < min_b.z {
 min_b.z - max_a.z
 } else if max_b.z < min_a.z {
 min_a.z - max_b.z
 } else {
 0.0
 };
 (dx * dx + dy * dy + dz * dz).sqrt()
}

fn breps_farther_than_tolerance(a: &BRep, b: &BRep, tol: f64) -> bool {
 let Some((min_a, max_a)) = brep_bounds(a) else {
 return false;
 };
 let Some((min_b, max_b)) = brep_bounds(b) else {
 return false;
 };
 aabb_distance(min_a, max_a, min_b, max_b) > tol
}

fn splitter_issues_by_level(
 validity: &CheckResult,
 level: SplitterValidationLevel,
) -> (usize, Option<String>) {
 let filtered: Vec<&CheckIssue> = match level {
 SplitterValidationLevel::Relaxed => validity
 .issues
 .iter()
 .filter(|issue| !matches!(issue, CheckIssue::NonManifoldEdge { .. }))
 .collect(),
 SplitterValidationLevel::Strict => validity.issues.iter().collect(),
 };
 (filtered.len(), filtered.first().map(|it| it.to_string()))
}

/// Split each object by a shared set of tools.
///
/// This is a grouped splitter API similar to object/tool workflows in mature
/// boolean kernels: every input object is split against all tools, and results
/// are returned in object order.
pub fn split_objects_with_tools(
 objects: &[BRep],
 tools: &[BRep],
) -> (Vec<BRep>, SplitterObjectsReport) {
 split_objects_with_tools_options(objects, tools, SplitterOptions::default())
}

/// Like [`split_objects_with_tools`] but with advanced options.
pub fn split_objects_with_tools_options(
 objects: &[BRep],
 tools: &[BRep],
 options: SplitterOptions,
) -> (Vec<BRep>, SplitterObjectsReport) {
 let mut outputs = Vec::with_capacity(objects.len());
 let mut objects_report = Vec::with_capacity(objects.len());

 for (object_index, object) in objects.iter().enumerate() {
 let (split, report) = split_shape_with_options(object, tools, options);
 outputs.push(split);
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: true,
 error: None,
 });
 }

 (
 outputs,
 SplitterObjectsReport {
 objects: objects_report,
 },
 )
}

/// Checked grouped splitter variant.
///
/// Validates each split step for each object and returns the first error.
pub fn split_objects_with_tools_checked_options(
 objects: &[BRep],
 tools: &[BRep],
 options: SplitterOptions,
) -> Result<(Vec<BRep>, SplitterObjectsReport), SplitterError> {
 let mut outputs = Vec::with_capacity(objects.len());
 let mut objects_report = Vec::with_capacity(objects.len());

 for (object_index, object) in objects.iter().enumerate() {
 let (split, report) = split_shape_checked_with_options(object, tools, options)?;
 outputs.push(split);
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: true,
 error: None,
 });
 }

 Ok((
 outputs,
 SplitterObjectsReport {
 objects: objects_report,
 },
 ))
}

/// Checked grouped splitter with per-object failure collection.
///
/// Unlike [`split_objects_with_tools_checked_options`], this function does not
/// fail fast. It records per-object errors in the returned report and keeps
/// processing remaining objects.
pub fn split_objects_with_tools_checked_collect_options(
 objects: &[BRep],
 tools: &[BRep],
 options: SplitterOptions,
) -> (Vec<Option<BRep>>, SplitterObjectsReport) {
 let mut outputs = Vec::with_capacity(objects.len());
 let mut objects_report = Vec::with_capacity(objects.len());

 for (object_index, object) in objects.iter().enumerate() {
 let (result, report) =
 split_brep_internal_with_partial_report(object, tools, options, true);
 match result {
 Ok(split) => {
 outputs.push(Some(split));
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: true,
 error: None,
 });
 }
 Err(err) => {
 outputs.push(None);
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: false,
 error: Some(err),
 });
 }
 }
 }

 (
 outputs,
 SplitterObjectsReport {
 objects: objects_report,
 },
 )
}

/// Like [`boolean_op`] but also returns a [`BooleanHistory`] mapping each result
/// face back to its source in solid A or B.
pub fn boolean_op_with_history(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
 if matches!(op, BooleanOpType::Union) {
 let (t, hist) = bop_occt_union::fuse_with_history(&a.to_topods(), &b.to_topods())?;
 return Ok((rcad_kernel::BRep::from_topods(&t), hist));
 }

 let a_t = a.to_topods();
 let b_t = b.to_topods();
 let mut ds = bopds::ds::DS::new_from_topods(&a_t, &b_t, TOLERANCE_ABS);
 let fuzzy_tol = ds.fuzzy_tol;
 let mut brep = rcad_kernel::topods::BRep::new();
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let (face_refs, ic_edge_map) = {
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.set_run_parallel(false);
 filler.configure_fuzzy(fuzzy_tol);
 filler.set_non_destructive(false);
 filler.configure_glue(false, TOLERANCE_ABS);
 filler.set_use_obb(false);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
 let (t, history) = builder.build_with_history()?;
 Ok((rcad_kernel::BRep::from_topods(&t), history))
}

pub fn boolean_op_par(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
) -> Result<(rcad_kernel::BRep, BooleanHistory), BooleanError> {
 if matches!(op, BooleanOpType::Union) {
 let (t, h) = bop_occt_union::fuse_with_history_par(&a.to_topods(), &b.to_topods())?;
 return Ok((t, h));
 }

 let a_t = a.to_topods();
 let b_t = b.to_topods();
 let mut ds = bopds::ds::DS::new_from_topods(&a_t, &b_t, TOLERANCE_ABS);
 let fuzzy_tol = ds.fuzzy_tol;
 let mut brep = rcad_kernel::topods::BRep::new();
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let (face_refs, ic_edge_map) = {
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.set_run_parallel(true);
 filler.configure_fuzzy(fuzzy_tol);
 filler.set_non_destructive(false);
 filler.configure_glue(false, TOLERANCE_ABS);
 filler.set_use_obb(false);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
 let (t, history) = builder.build_with_history()?;
 Ok((rcad_kernel::BRep::from_topods(&t), history))
}

/// Check if any solid in the BRep has at least one face (deep check across all solids).
fn has_any_face(brep: &BRep) -> bool {
 brep.solids
 .iter()
 .any(|s| s.shells.iter().any(|sh| !sh.faces.is_empty()))
}

/// Build BVHs for both BReps if they have faces; returns None for empty BReps.
pub(crate) fn build_optional_bvhs(a: &BRep, b: &BRep) -> (Option<bvh::Bvh>, Option<bvh::Bvh>) {
 let has_faces_a = a
 .solids
 .first()
 .and_then(|s| s.shells.first())
 .is_some_and(|sh| !sh.faces.is_empty());
 let has_faces_b = b
 .solids
 .first()
 .and_then(|s| s.shells.first())
 .is_some_and(|sh| !sh.faces.is_empty());
 (
 if has_faces_a {
 Some(bvh::Bvh::build(a))
 } else {
 None
 },
 if has_faces_b {
 Some(bvh::Bvh::build(b))
 } else {
 None
 },
 )
}

fn has_faces(brep: &BRep) -> bool {
 brep.solids
 .first()
 .and_then(|s| s.shells.first())
 .is_some_and(|sh| !sh.faces.is_empty())
}

fn make_connected_seed_vertices_from_short_edges(brep: &BRep, seed_length: f64) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 let threshold = seed_length.max(tolerance::TOLERANCE_ABS);
 for e in &brep.edges {
 if e.start >= brep.vertices.len() || e.end >= brep.vertices.len() {
 continue;
 }
 let ps = brep.vertices[e.start].point;
 let pe = brep.vertices[e.end].point;
 if (pe - ps).length() <= threshold {
 out.insert(e.start);
 out.insert(e.end);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_near_duplicates(brep: &BRep, seed_length: f64) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 let threshold = seed_length.max(tolerance::TOLERANCE_ABS);
 let threshold2 = threshold * threshold;
 for i in 0..brep.vertices.len() {
 for j in (i + 1)..brep.vertices.len() {
 let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
 if d2 <= threshold2 {
 out.insert(i);
 out.insert(j);
 }
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_tolerance_tagged_edges(
 brep: &BRep,
 tolerance_threshold: f64,
) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 let threshold = tolerance_threshold.max(tolerance::TOLERANCE_ABS);
 for (ei, e) in brep.edges.iter().enumerate() {
 let edge_tol = brep
 .geom
 .edge_tolerance
 .get(ei)
 .copied()
 .unwrap_or(tolerance::TOLERANCE_ABS);
 if edge_tol >= threshold {
 out.insert(e.start);
 out.insert(e.end);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_multi_pcurve_edges(brep: &BRep) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 for (ei, e) in brep.edges.iter().enumerate() {
 if brep
 .geom
 .edge_pcurves
 .get(ei)
 .map(|pcs| pcs.len() >= 2)
 .unwrap_or(false)
 {
 out.insert(e.start);
 out.insert(e.end);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_topology_seam_candidates(brep: &BRep) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 for ei in rcad_kernel::periodic_seam_edge_indices(brep) {
 if let Some(e) = brep.edges.get(ei) {
 out.insert(e.start);
 out.insert(e.end);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_edges_from_short_edges(brep: &BRep, seed_length: f64) -> Vec<usize> {
 let mut out = Vec::new();
 let threshold = seed_length.max(tolerance::TOLERANCE_ABS);
 for (ei, e) in brep.edges.iter().enumerate() {
 if e.start >= brep.vertices.len() || e.end >= brep.vertices.len() {
 continue;
 }
 let ps = brep.vertices[e.start].point;
 let pe = brep.vertices[e.end].point;
 if (pe - ps).length() <= threshold {
 out.push(ei);
 }
 }
 out
}

fn make_connected_seed_edges_from_near_duplicates(brep: &BRep, seed_length: f64) -> Vec<usize> {
 let dup_vertices: std::collections::HashSet<usize> =
 make_connected_seed_vertices_from_near_duplicates(brep, seed_length)
 .into_iter()
 .collect();
 brep.edges
 .iter()
 .enumerate()
 .filter(|(_, e)| dup_vertices.contains(&e.start) || dup_vertices.contains(&e.end))
 .map(|(ei, _)| ei)
 .collect()
}

fn make_connected_seed_edges_from_tolerance_tagged_edges(
 brep: &BRep,
 tolerance_threshold: f64,
) -> Vec<usize> {
 let threshold = tolerance_threshold.max(tolerance::TOLERANCE_ABS);
 brep.edges
 .iter()
 .enumerate()
 .filter(|(ei, _)| {
 brep.geom
 .edge_tolerance
 .get(*ei)
 .copied()
 .unwrap_or(tolerance::TOLERANCE_ABS)
 >= threshold
 })
 .map(|(ei, _)| ei)
 .collect()
}

fn make_connected_seed_edges_from_multi_pcurve_edges(brep: &BRep) -> Vec<usize> {
 brep.edges
 .iter()
 .enumerate()
 .filter(|(ei, _)| {
 brep.geom
 .edge_pcurves
 .get(*ei)
 .map(|pcs| pcs.len() >= 2)
 .unwrap_or(false)
 })
 .map(|(ei, _)| ei)
 .collect()
}

fn make_connected_seed_edges_from_topology_seam_candidates(brep: &BRep) -> Vec<usize> {
 rcad_kernel::periodic_seam_edge_indices(brep)
}

fn make_connected_seed_edges(
 brep: &BRep,
 seed_length: f64,
 mode: MakeConnectedScopeSeedMode,
) -> Vec<usize> {
 match mode {
 MakeConnectedScopeSeedMode::ShortEdges => {
 make_connected_seed_edges_from_short_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::NearDuplicateVertices => {
 make_connected_seed_edges_from_near_duplicates(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::ToleranceTaggedEdges => {
 make_connected_seed_edges_from_tolerance_tagged_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::MultiPcurveEdges => {
 make_connected_seed_edges_from_multi_pcurve_edges(brep)
 }
 MakeConnectedScopeSeedMode::TopologySeamCandidates => {
 make_connected_seed_edges_from_topology_seam_candidates(brep)
 }
 MakeConnectedScopeSeedMode::Hybrid => {
 let mut set = std::collections::BTreeSet::new();
 for ei in make_connected_seed_edges_from_short_edges(brep, seed_length) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_near_duplicates(brep, seed_length) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_tolerance_tagged_edges(brep, seed_length) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_multi_pcurve_edges(brep) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_topology_seam_candidates(brep) {
 set.insert(ei);
 }
 set.into_iter().collect()
 }
 }
}

pub(crate) fn make_connected_seed_vertices_from_edge_ids(brep: &BRep, edge_ids: &[usize]) -> Vec<usize> {
 let mut set = std::collections::BTreeSet::new();
 for &ei in edge_ids {
 if let Some(e) = brep.edges.get(ei) {
 set.insert(e.start);
 set.insert(e.end);
 }
 }
 set.into_iter().collect()
}

pub(crate) fn select_scoped_seed_edges(
 brep: &BRep,
 history: Option<&BooleanHistory>,
 seed_length: f64,
 mode: MakeConnectedScopeSeedMode,
 history_ring_depth: usize,
 min_history_edges: usize,
) -> (Vec<usize>, usize, usize, MakeConnectedScopeSeedSource) {
 let history_seed_edges_raw = history
 .map(|h| make_connected_seed_edges_from_boolean_history(brep, h))
 .unwrap_or_default();
 // Expand history-derived seeds to configurable ring depth around boolean
 // interface topology while preserving raw-history count semantics for reports.
 let history_seed_edges =
 expand_seed_edges_with_ring_depth(brep, &history_seed_edges_raw, history_ring_depth);
 let heuristic_seed_edges = make_connected_seed_edges(brep, seed_length, mode);

 if history_seed_edges_raw.is_empty() {
 return (
 heuristic_seed_edges.clone(),
 0,
 heuristic_seed_edges.len(),
 MakeConnectedScopeSeedSource::Heuristic,
 );
 }

 if history_seed_edges_raw.len() < min_history_edges {
 let mut set = std::collections::BTreeSet::new();
 for ei in &history_seed_edges {
 set.insert(*ei);
 }
 for ei in &heuristic_seed_edges {
 set.insert(*ei);
 }
 return (
 set.into_iter().collect(),
 history_seed_edges_raw.len(),
 heuristic_seed_edges.len(),
 MakeConnectedScopeSeedSource::HistoryAugmentedHeuristic,
 );
 }

 (
 history_seed_edges.clone(),
 history_seed_edges_raw.len(),
 heuristic_seed_edges.len(),
 MakeConnectedScopeSeedSource::History,
 )
}

fn expand_seed_edges_with_ring_depth(
 brep: &BRep,
 seed_edges: &[usize],
 ring_depth: usize,
) -> Vec<usize> {
 let mut out: std::collections::BTreeSet<usize> = seed_edges.iter().copied().collect();
 if ring_depth == 0 || seed_edges.is_empty() {
 return out.into_iter().collect();
 }

 let mut visited_faces = std::collections::BTreeSet::new();
 let mut frontier = std::collections::BTreeSet::new();
 for &ei in seed_edges {
 for fi in rcad_kernel::edge_adjacent_faces(brep, ei) {
 if visited_faces.insert(fi) {
 frontier.insert(fi);
 }
 }
 }

 for _ in 0..ring_depth {
 if frontier.is_empty() {
 break;
 }
 let current: Vec<usize> = frontier.iter().copied().collect();
 frontier.clear();

 for fi in current {
 for fei in rcad_kernel::face_edges(brep, fi) {
 out.insert(fei);
 for nfi in rcad_kernel::edge_adjacent_faces(brep, fei) {
 if visited_faces.insert(nfi) {
 frontier.insert(nfi);
 }
 }
 }
 }
 }

 out.into_iter().collect()
}

fn make_connected_seed_edges_from_boolean_history(
 brep: &BRep,
 history: &BooleanHistory,
) -> Vec<usize> {
 let mut seed_edges = std::collections::BTreeSet::new();

 // If edge history is available, prefer boundary-like generated/split edges.
 for (ei, origin) in history.edge_origins.iter().enumerate() {
 if ei >= brep.edges.len() {
 break;
 }
 if matches!(
 origin,
 EdgeOrigin::Generated | EdgeOrigin::SplitFromA(_) | EdgeOrigin::SplitFromB(_)
 ) {
 seed_edges.insert(ei);
 }
 }

 // Fallback semantic extraction from face history: edges adjacent to both A and B faces
 // are strong candidates for boolean interface cleanup.
 for ei in 0..brep.edges.len() {
 let adjacent = rcad_kernel::edge_adjacent_faces(brep, ei);
 if adjacent.is_empty() {
 continue;
 }
 let mut has_a = false;
 let mut has_b = false;
 let mut has_generated = false;
 for fi in adjacent {
 if fi >= history.face_origins.len() {
 continue;
 }
 match history.face_origins[fi] {
 FaceOrigin::FromA(_) => has_a = true,
 FaceOrigin::FromB(_) => has_b = true,
 FaceOrigin::Generated => has_generated = true,
 }
 }
 if has_generated || (has_a && has_b) {
 seed_edges.insert(ei);
 }
 }

 seed_edges.into_iter().collect()
}

pub(crate) fn make_connected_seed_edge_labels(brep: &BRep, edge_ids: &[usize]) -> Vec<String> {
 edge_ids
 .iter()
 .map(|&ei| match brep.edges.get(ei) {
 Some(e) => {
 let pa = brep.vertices.get(e.start).map(|v| v.point);
 let pb = brep.vertices.get(e.end).map(|v| v.point);
 match (pa, pb) {
 (Some(a), Some(b)) => {
 let a_label = format!("{:.9},{:.9},{:.9}", a.x, a.y, a.z);
 let b_label = format!("{:.9},{:.9},{:.9}", b.x, b.y, b.z);
 if a_label <= b_label {
 format!("edge.{ei}.{a_label}->{b_label}")
 } else {
 format!("edge.{ei}.{b_label}->{a_label}")
 }
 }
 _ => format!("edge.{ei}.invalid-vertex"),
 }
 }
 None => format!("edge.{ei}.invalid-edge"),
 })
 .collect()
}

pub(crate) fn make_connected_seed_vertices(
 brep: &BRep,
 seed_length: f64,
 mode: MakeConnectedScopeSeedMode,
) -> Vec<usize> {
 match mode {
 MakeConnectedScopeSeedMode::ShortEdges => {
 make_connected_seed_vertices_from_short_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::NearDuplicateVertices => {
 make_connected_seed_vertices_from_near_duplicates(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::ToleranceTaggedEdges => {
 make_connected_seed_vertices_from_tolerance_tagged_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::MultiPcurveEdges => {
 make_connected_seed_vertices_from_multi_pcurve_edges(brep)
 }
 MakeConnectedScopeSeedMode::TopologySeamCandidates => {
 make_connected_seed_vertices_from_topology_seam_candidates(brep)
 }
 MakeConnectedScopeSeedMode::Hybrid => {
 let mut set = std::collections::BTreeSet::new();
 for v in make_connected_seed_vertices_from_short_edges(brep, seed_length) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_near_duplicates(brep, seed_length) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_tolerance_tagged_edges(brep, seed_length) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_multi_pcurve_edges(brep) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_topology_seam_candidates(brep) {
 set.insert(v);
 }
 set.into_iter().collect()
 }
 }
}

/// Create stable per-face labels from boolean history.
pub fn persistent_face_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .face_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 FaceOrigin::FromA(src) => format!("face.{idx}.A.{src}"),
 FaceOrigin::FromB(src) => format!("face.{idx}.B.{src}"),
 FaceOrigin::Generated => format!("face.{idx}.G"),
 })
 .collect()
}

/// Create stable per-edge labels from boolean history.
pub fn persistent_edge_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .edge_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 EdgeOrigin::FromA(src) => format!("edge.{idx}.A.{src}"),
 EdgeOrigin::FromB(src) => format!("edge.{idx}.B.{src}"),
 EdgeOrigin::Generated => format!("edge.{idx}.G"),
 EdgeOrigin::SplitFromA(src) => format!("edge.{idx}.A.split.{src}"),
 EdgeOrigin::SplitFromB(src) => format!("edge.{idx}.B.split.{src}"),
 })
 .collect()
}

/// Create stable per-shell labels from boolean history.
pub fn persistent_shell_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .shell_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 ShellOrigin::FromA => format!("shell.{idx}.A"),
 ShellOrigin::FromB => format!("shell.{idx}.B"),
 ShellOrigin::Generated => format!("shell.{idx}.G"),
 ShellOrigin::Mixed => format!("shell.{idx}.M"),
 })
 .collect()
}

/// Create stable per-solid labels from boolean history.
pub fn persistent_solid_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .solid_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 SolidOrigin::FromA => format!("solid.{idx}.A"),
 SolidOrigin::FromB => format!("solid.{idx}.B"),
 SolidOrigin::Generated => format!("solid.{idx}.G"),
 SolidOrigin::Mixed => format!("solid.{idx}.M"),
 })
 .collect()
}

/// Union two BReps and return both the result and face origin history.
pub fn union_with_history(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
 boolean_op_with_history(BooleanOpType::Union, a, b)
}

/// Intersect two BReps and return both the result and face origin history.
pub fn intersection_with_history(
 a: &BRep,
 b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
 boolean_op_with_history(BooleanOpType::Intersection, a, b)
}

/// Subtract solid B from solid A and return both the result and face origin history.
pub fn difference_with_history(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
 boolean_op_with_history(BooleanOpType::Difference, a, b)
}

/// Run boolean operation followed by structured healing using default options.
pub fn boolean_op_healed(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
) -> Result<(BRep, HealingReport), BooleanError> {
 let raw = boolean_op(op, &a.to_topods(), &b.to_topods())?;
 let mut healing = HealingOptions::default();
 align_healing_options_with_boolean_operands(&mut healing, a, b, 0.0);
 let (healed, report) = analyze_and_heal(&raw, healing);
 Ok((rcad_kernel::BRep::from_topods(&healed), report))
}

/// Run boolean operation followed by structured healing using custom options.
pub fn boolean_op_healed_with_options(
 op: BooleanOpType,
 a: &BRep,
 b: &BRep,
 mut options: HealingOptions,
) -> Result<(BRep, HealingReport), BooleanError> {
 let raw = boolean_op(op, &a.to_topods(), &b.to_topods())?;
 align_healing_options_with_boolean_operands(&mut options, a, b, 0.0);
 let (healed, report) = analyze_and_heal(&raw, options);
 Ok((rcad_kernel::BRep::from_topods(&healed), report))
}

/// Multi-body boolean fuse (union) over a list of solids.
///
/// Delegates to [`general_fuse_with_options`] with [`BooleanOptions::default`]. Each fold step
/// uses [`boolean_op_with_options`], so pairwise [`merge_pairwise_model_tol_into_boolean_options`]
/// runs on every `(accumulator, part)` pair.
pub fn general_fuse(parts: &[BRep]) -> Result<BRep, BooleanError> {
 general_fuse_with_options(parts, BooleanOptions::default())
}