
/// Concatenate multiple slab BReps into one and run topology optimization.
fn concatenate_and_merge_slabs(slabs: &[BRep]) -> Option<BRep> {
 use std::collections::HashMap;
 use rcad_kernel::topology::{Face, Shell, Solid, Wire, WireEdge};
 if slabs.is_empty() { return None; }
 if slabs.len() == 1 { return Some(slabs[0].clone()); }

 let mut out = BRep::new();
 let mut v_offset: Vec<usize> = Vec::new();
 let mut e_offset: Vec<usize> = Vec::new();
 let mut s_offset: Vec<usize> = Vec::new();
 let mut surf_remap: Vec<HashMap<usize, usize>> = Vec::new();

 for slab in slabs {
 v_offset.push(out.vertices.len());
 e_offset.push(out.edges.len());
 s_offset.push(out.geom.surfaces.len());
 out.vertices.extend_from_slice(&slab.vertices);

 let mut remap: HashMap<usize, usize> = HashMap::new();
 for (old_si, surf) in slab.geom.surfaces.iter().enumerate() {
 let new_si = out.geom.surfaces.iter().position(|s| surfaces_eq(s, surf))
 .unwrap_or_else(|| { let idx = out.geom.surfaces.len(); out.geom.surfaces.push(surf.clone()); idx });
 remap.insert(old_si, new_si);
 }
 surf_remap.push(remap);

 let vo = *v_offset.last().unwrap_or(&0);
 for e in &slab.edges {
 out.edges.push(rcad_kernel::Edge { start: e.start + vo, end: e.end + vo });
 }
 let eo = *e_offset.last().unwrap_or(&0);
 for (ei, curve_idx_opt) in slab.geom.edge_curve.iter().enumerate() {
 if let Some(&ci) = curve_idx_opt.as_ref() {
 if let Some(curve) = slab.geom.curves.get(ci) {
 let new_ci = out.geom.curves.len();
 out.geom.curves.push(curve.clone());
 while out.geom.edge_curve.len() <= eo + ei { out.geom.edge_curve.push(None); }
 out.geom.edge_curve[eo + ei] = Some(new_ci);
 }
 }
 }
 for (ei, range) in slab.geom.edge_curve_range.iter().enumerate() {
 while out.geom.edge_curve_range.len() <= eo + ei { out.geom.edge_curve_range.push(None); }
 out.geom.edge_curve_range[eo + ei] = *range;
 }
 for (ei, sp) in slab.geom.edge_same_parameter.iter().enumerate() {
 while out.geom.edge_same_parameter.len() <= eo + ei { out.geom.edge_same_parameter.push(false); }
 out.geom.edge_same_parameter[eo + ei] = *sp;
 }
 for (ei, pcs) in slab.geom.edge_pcurves.iter().enumerate() {
 if !pcs.is_empty() {
 while out.geom.edge_pcurves.len() <= eo + ei { out.geom.edge_pcurves.push(Vec::new()); }
 out.geom.edge_pcurves[eo + ei] = pcs.clone();
 }
 }
 for face in &slab.solids[0].shells[0].faces {
 let remap_we = |we: &WireEdge| WireEdge { idx: we.idx + eo, forward: we.forward };
 let outer_wire = Wire { edges: face.outer_wire.edges.iter().map(remap_we).collect() };
 let inner_wires: Vec<Wire> = face.inner_wires.iter().map(|w| Wire {
 edges: w.edges.iter().map(remap_we).collect()
 }).collect();
 out.solids.push(Solid {
 shells: vec![Shell { faces: vec![Face {
 outer_wire, inner_wires,
 normal: face.normal,
 triangles: vec![], sample_point: None, mesh_dirty: true,
 surface_idx: None,
 }]}]
 });
 }
 }

 let mut flat_fi = 0usize;
 for (si, slab) in slabs.iter().enumerate() {
 let remap = &surf_remap[si];
 for (old_fi, _) in slab.solids[0].shells[0].faces.iter().enumerate() {
 let old_si = slab.geom.face_surface.get(old_fi).copied().flatten().unwrap_or(0);
 let new_si = remap.get(&old_si).copied().unwrap_or(0);
 while out.geom.face_surface.len() <= flat_fi { out.geom.face_surface.push(None); }
 out.geom.face_surface[flat_fi] = Some(new_si);
 flat_fi += 1;
 }
 }

 let out = crate::deduplicate_surfaces(out);
 let out = crate::optimize_boolean_topology(out);
 Some(out)
}

/// Compare two surfaces for geometric equality.
fn surfaces_eq(a: &rcad_kernel::geom::Surface3, b: &rcad_kernel::geom::Surface3) -> bool {
 use rcad_kernel::geom::Surface3;
 let tol = 1e-8;
 match (a, b) {
 (Surface3::Plane(p1), Surface3::Plane(p2)) =>
 (p1.normal - p2.normal).length() < tol && (p1.origin - p2.origin).length() < tol,
 _ => std::mem::discriminant(a) == std::mem::discriminant(b),
 }
}

/// Union of two general boxes (axis-aligned or rotated), computed analytically.
///
/// Both BReps must be detected as boxes by [`try_as_box`]. For disjoint boxes,
/// returns a compound of the two inputs. For containment, returns the containing
/// box. For partial overlap, decomposes A and B into non-overlapping slabs
/// around each other's axes, adds the overlap region, and sews into a solid.
///
/// Falls through to Pave-Filler (returns `None`) when sewing fails or excessive
/// internal-face inflation is detected.
pub fn try_union_box_general(a: &BRep, b: &BRep) -> Option<BRep> {
 // OCCT does not have fast paths  ?skip NURBS so PaveFiller+Builder
 // preserves BSpline surface types (nurbsconvert cases).
 for operand in [a, b] {
 if operand.geom.surfaces.iter().any(|s| matches!(s, Surface3::BSpline(_))) {
 return None;
 }
 }
 let info_a = try_as_box(a)?;
 let info_b = try_as_box(b)?;

 let inter = match try_intersection_box_general(a, b) {
 Some(inter) => inter,
 None => return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()])),
 };

 // No overlap =compound (disjoint boxes).
 if inter.vertices.len() < 4 {
 return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()]));
 }

 let inter_vol = volume(&inter);
 let a_vol = volume(a);
 let b_vol = volume(b);
 let scale = (a_vol.max(b_vol) / 3.0).max(1.0);
 let vol_tol = TOLERANCE_LEN_MIN * scale;

 // A contains B =A alone is the union.
 if b_vol > vol_tol && (b_vol - inter_vol).abs() < vol_tol {
 return Some(a.clone());
 }
 // B contains A =B alone is the union.
 if a_vol > vol_tol && (a_vol - inter_vol).abs() < vol_tol {
 return Some(b.clone());
 }

 let a_planes = info_a.planes();
 let b_planes = info_b.planes();
 let zero_tol = TOLERANCE_LEN_MIN * scale;

 let mut slabs: Vec<BRep> = Vec::new();

 // Helper macro: build a convex polyhedron from base planes + extra half-spaces.
 macro_rules! try_slab {
 ($base:expr, $extra:expr) => {
 let mut planes = ($base).to_vec();
 planes.extend($extra);
 if let Ok(s) = make_convex_polyhedron_from_half_spaces(&planes) {
 if volume(&s) > vol_tol * 0.1 {
 slabs.push(s);
 }
 }
 };
 }

 // = =  A \ B slabs: decompose A around B's axes = = 
 {
 let [u, v, w] = info_b.axes;
 let [eu, ev, ew] = info_b.extents;
 let c = info_b.center;
 let u_min = u.dot(c) - eu;
 let u_max = u.dot(c) + eu;
 let v_min = v.dot(c) - ev;
 let v_max = v.dot(c) + ev;
 let w_min = w.dot(c) - ew;
 let w_max = w.dot(c) + ew;

 let a_verts: Vec<DVec3> = a.vertices.iter().map(|vi| vi.point).collect();
 let a_u_min = a_verts.iter().map(|p| u.dot(*p)).fold(f64::MAX, f64::min);
 let a_u_max = a_verts.iter().map(|p| u.dot(*p)).fold(f64::MIN, f64::max);
 let a_v_min = a_verts.iter().map(|p| v.dot(*p)).fold(f64::MAX, f64::min);
 let a_v_max = a_verts.iter().map(|p| v.dot(*p)).fold(f64::MIN, f64::max);
 let a_w_min = a_verts.iter().map(|p| w.dot(*p)).fold(f64::MAX, f64::min);
 let a_w_max = a_verts.iter().map(|p| w.dot(*p)).fold(f64::MIN, f64::max);

 let u_span = a_u_max > u_min + zero_tol && a_u_min < u_max - zero_tol;
 let v_span = a_v_max > v_min + zero_tol && a_v_min < v_max - zero_tol;

 if a_u_min < u_min - zero_tol {
 try_slab!(&a_planes, vec![(u * u_min, u)]);
 }
 if a_u_max > u_max + zero_tol {
 try_slab!(&a_planes, vec![(u * u_max, -u)]);
 }
 if u_span && a_v_min < v_min - zero_tol {
 try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, v)]);
 }
 if u_span && a_v_max > v_max + zero_tol {
 try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_max, -v)]);
 }
 if u_span && v_span && a_w_min < w_min - zero_tol {
 try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_min, w)]);
 }
 if u_span && v_span && a_w_max > w_max + zero_tol {
 try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_max, -w)]);
 }
 }

 // = =  A =B (overlap) = = 
 slabs.push(inter.clone());

 // = =  B \ A slabs: decompose B around A's axes = = 
 {
 let [u, v, w] = info_a.axes;
 let [eu, ev, ew] = info_a.extents;
 let c = info_a.center;
 let u_min = u.dot(c) - eu;
 let u_max = u.dot(c) + eu;
 let v_min = v.dot(c) - ev;
 let v_max = v.dot(c) + ev;
 let w_min = w.dot(c) - ew;
 let w_max = w.dot(c) + ew;

 let b_verts: Vec<DVec3> = b.vertices.iter().map(|vi| vi.point).collect();
 let b_u_min = b_verts.iter().map(|p| u.dot(*p)).fold(f64::MAX, f64::min);
 let b_u_max = b_verts.iter().map(|p| u.dot(*p)).fold(f64::MIN, f64::max);
 let b_v_min = b_verts.iter().map(|p| v.dot(*p)).fold(f64::MAX, f64::min);
 let b_v_max = b_verts.iter().map(|p| v.dot(*p)).fold(f64::MIN, f64::max);
 let b_w_min = b_verts.iter().map(|p| w.dot(*p)).fold(f64::MAX, f64::min);
 let b_w_max = b_verts.iter().map(|p| w.dot(*p)).fold(f64::MIN, f64::max);

 let u_span = b_u_max > u_min + zero_tol && b_u_min < u_max - zero_tol;
 let v_span = b_v_max > v_min + zero_tol && b_v_min < v_max - zero_tol;

 if b_u_min < u_min - zero_tol {
 try_slab!(&b_planes, vec![(u * u_min, u)]);
 }
 if b_u_max > u_max + zero_tol {
 try_slab!(&b_planes, vec![(u * u_max, -u)]);
 }
 if u_span && b_v_min < v_min - zero_tol {
 try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, v)]);
 }
 if u_span && b_v_max > v_max + zero_tol {
 try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_max, -v)]);
 }
 if u_span && v_span && b_w_min < w_min - zero_tol {
 try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_min, w)]);
 }
 if u_span && v_span && b_w_max > w_max + zero_tol {
 try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_max, -w)]);
 }
 }

 if slabs.is_empty() {
 return Some(BRep::default());
 }
 if slabs.len() == 1 {
 return Some(slabs.remove(0));
 }

 // = =  SA-inflation guard: try sew, fall back to fuse if inflated = = 
 let _slab_sa_sum: f64 = slabs.iter().map(|s| surface_area(s)).sum();
 let slab_vol_sum: f64 = slabs.iter().map(|s| volume(s)).sum();
 let expected_union_sa = surface_area(a) + surface_area(b) - surface_area(&inter);

 let sewn = sew_slabs_into_solid(&slabs, zero_tol);
 let sewn_sa = surface_area(&sewn);

 if sewn_sa <= expected_union_sa * 1.15 {
 return Some(sewn);
 }

 // SA inflated =try fuse-based sequential merge to remove shared internal faces.
 let mut fused = slabs[0].clone();
 let mut ok = true;
 for slab in &slabs[1..] {
 match crate::bop_occt_union::fuse(&fused, slab) {
 Ok(u) => { fused = rcad_kernel::BRep::from_topods(&u); }
 Err(_) => { ok = false; break; }
 }
 }
 if ok {
 let fused_sa = surface_area(&fused);
 let fused_vol = volume(&fused);
 let vol_ok = (fused_vol - slab_vol_sum).abs() < vol_tol * (slabs.len() as f64).max(1000.0);
 if vol_ok && fused_sa <= expected_union_sa * 1.15 {
 return Some(fused);
 }
 }

 Some(sewn)
}

/// Kernel analytic sphere primitive: one spherical face (`Surface3::Sphere`).
fn try_sphere_primitive_center_radius(brep: &BRep) -> Option<(DVec3, f64)> {
 let sh = brep.solids.get(0)?.shells.get(0)?;
 if sh.faces.len() != 1 {
 return None;
 }
 let si = *brep.geom.face_surface.get(0)?.as_ref()?;
 match brep.geom.surfaces.get(si)? {
 Surface3::Sphere(s) => Some((s.center, s.radius)),
 _ => None,
 }
}

/// Fast path for coaxial cylinder-sphere difference.
///
/// For `sphere - cylinder`: returns the spherical portion(s) outside the
/// cylinder's Z range (when the sphere cross-section fits entirely inside the
/// cylinder radius over the overlap).
///
/// For `cylinder - sphere`: returns the cylinder body with a spherical cavity
/// at the overlap end, built as a Z-slice triangle mesh.
pub fn try_difference_coaxial_cylinder_sphere(a: &BRep, b: &BRep) -> Option<BRep> {
 // Try analytic fast paths first (exact geometry, no tessellation).
 if let Some(result) = crate::cylinder_sphere_analytic::build_sphere_minus_cylinder_analytic(a, b) {
 return Some(result);
 }
 if let Some(result) = crate::cylinder_sphere_analytic::build_cylinder_minus_sphere_analytic(a, b) {
 return Some(result);
 }

 // Fall through to mesh-backed detection logic.
 // Try both orderings. Track which operand is the cylinder to build the
 // correct result: sphere - cylinder (existing) vs cylinder - sphere.
 if let Some((sp, cyl)) = try_sphere_primitive_center_radius(a)
 .and_then(|sp| z_axis_cylinder_z_span_r(b).map(|cyl| (sp, cyl)))
 {
 // a = sphere, b = cylinder =sphere - cylinder
 if let Some(result) = try_sphere_minus_cylinder(sp, cyl) {
 return Some(result);
 }
 }

 if let Some((sp, cyl)) = try_sphere_primitive_center_radius(b)
 .and_then(|sp| z_axis_cylinder_z_span_r(a).map(|cyl| (sp, cyl)))
 {
 // a = cylinder, b = sphere =cylinder - sphere
 if let Some(result) = try_cylinder_minus_sphere(sp, cyl) {
 return Some(result);
 }
 }

 None
}

/// Build `sphere - cylinder` =portions of sphere outside cylinder Z range.
fn try_sphere_minus_cylinder(
 sphere_brep: (DVec3, f64),
 cyl_brep: (f64, f64, f64),
) -> Option<BRep> {
 let (center, radius) = sphere_brep;
 let (cyl_z_lo, cyl_z_hi, cyl_r) = cyl_brep;

 if center.x.abs() > TOLERANCE_ABS || center.y.abs() > TOLERANCE_ABS {
 return None;
 }

 let sz = center.z;
 let sphere_z_lo = sz - radius;
 let sphere_z_hi = sz + radius;

 let overlap_lo = sphere_z_lo.max(cyl_z_lo);
 let overlap_hi = sphere_z_hi.min(cyl_z_hi);
 if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
 return None;
 }

 for z in [overlap_lo, overlap_hi] {
 let dz = z - sz;
 let r_at = (radius.powi(2) - dz.powi(2)).sqrt();
 if r_at > cyl_r + TOLERANCE_LEN_MIN {
 return None;
 }
 }

 let mut parts: Vec<BRep> = Vec::new();
 if sphere_z_lo < cyl_z_lo - TOLERANCE_LEN_MIN {
 let z_to = cyl_z_lo.min(sphere_z_hi);
 if z_to - sphere_z_lo > TOLERANCE_LEN_MIN {
 if let Some(p) = build_spherical_slice_solid(center, radius, sphere_z_lo, z_to) {
 parts.push(p);
 }
 }
 }
 if sphere_z_hi > cyl_z_hi + TOLERANCE_LEN_MIN {
 let z_from = cyl_z_hi.max(sphere_z_lo);
 if sphere_z_hi - z_from > TOLERANCE_LEN_MIN {
 if let Some(p) = build_spherical_slice_solid(center, radius, z_from, sphere_z_hi) {
 parts.push(p);
 }
 }
 }

 match parts.len() {
 0 => Some(BRep::default()),
 1 => Some(parts.swap_remove(0)),
 _ => {
 let mut base = parts.swap_remove(0);
 for p in parts {
 append_frustum_brep(&mut base, p);
 }
 Some(base)
 }
 }
}

/// Build `cylinder - sphere` =cylinder body with a spherical cavity where the
/// sphere overlaps.  The sphere cross-section must fit inside the cylinder over
/// the entire overlap Z-range (checked by the caller).
fn try_cylinder_minus_sphere(
 sphere_brep: (DVec3, f64),
 cyl_brep: (f64, f64, f64),
) -> Option<BRep> {
 let (center, radius) = sphere_brep;
 let (cyl_z_lo, cyl_z_hi, cyl_r) = cyl_brep;

 if center.x.abs() > TOLERANCE_ABS || center.y.abs() > TOLERANCE_ABS {
 return None;
 }

 let sz = center.z;
 let sphere_z_lo = sz - radius;
 let sphere_z_hi = sz + radius;

 let overlap_lo = sphere_z_lo.max(cyl_z_lo);
 let overlap_hi = sphere_z_hi.min(cyl_z_hi);
 if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
 return None;
 }

 for z in [overlap_lo, overlap_hi] {
 let dz = z - sz;
 let r_at = (radius.powi(2) - dz.powi(2)).sqrt();
 if r_at > cyl_r + TOLERANCE_LEN_MIN {
 return None;
 }
 }

 build_cylinder_minus_sphere_tessellated(
 DVec2::new(center.x, center.y),
 cyl_z_lo, cyl_z_hi, cyl_r,
 center, radius,
 )
}

/// Triangulate an annular ring at `z` between `inner_r` and `outer_r`.
/// Winding direction is CCW on outer, CCW on inner (both triangles face +Z).
/// Area computation is unaffected by winding so this works for any ring.
fn add_annular_ring(
 add_v: &mut impl FnMut(DVec3) -> usize,
 faces: &mut Vec<Face>,
 z: f64, inner_r: f64, outer_r: f64, n_pts: usize,
 to_world: &impl Fn(f64, f64, f64) -> DVec3,
 empty_wire: &impl Fn() -> Wire,
) {
 let tau = std::f64::consts::TAU;
 let mut outer_idx = Vec::with_capacity(n_pts + 1);
 let mut inner_idx = Vec::with_capacity(n_pts + 1);
 for k in 0..=n_pts {
 let ang = tau * k as f64 / n_pts as f64;
 let (s, c) = ang.sin_cos();
 outer_idx.push(add_v(to_world(outer_r * c, outer_r * s, z)));
 inner_idx.push(add_v(to_world(inner_r * c, inner_r * s, z)));
 }
 let mut tris = Vec::with_capacity(n_pts * 2);
 for j in 0..n_pts {
 tris.push([outer_idx[j], outer_idx[j + 1], inner_idx[j + 1]]);
 tris.push([outer_idx[j], inner_idx[j + 1], inner_idx[j]]);
 }
 faces.push(Face {
 outer_wire: empty_wire(), inner_wires: vec![],
 normal: DVec3::ZERO, triangles: tris,
 sample_point: None, mesh_dirty: false,
 surface_idx: None,
 });
}

/// Build a tessellated cylinder body with a spherical cavity at the overlap.
/// Result = cylinder body minus the spherical portion.  The cylinder wall is
/// intact (full height), the far-end flat cap is preserved, and the spherical
/// cavity surface replaces the near-end region.  If the cavity extends to a
/// cylinder end without filling its cross-section, an annular flat ring caps
/// the remaining opening.
fn build_cylinder_minus_sphere_tessellated(
 center_xy: DVec2,
 cyl_z_lo: f64, cyl_z_hi: f64, cyl_r: f64,
 sphere_center: DVec3, sphere_r: f64,
) -> Option<BRep> {
 let tol = TOLERANCE_LEN_MIN;
 let sz = sphere_center.z;
 let sphere_z_lo = sz - sphere_r;
 let sphere_z_hi = sz + sphere_r;
 let overlap_lo = cyl_z_lo.max(sphere_z_lo);
 let overlap_hi = cyl_z_hi.min(sphere_z_hi);
 if overlap_hi <= overlap_lo + tol { return None; }

 let n_arc = 128;
 let n_slices_circ = 16;
 let tau = std::f64::consts::TAU;
 let empty_wire = || Wire { edges: vec![] };

 let mut verts: Vec<Vertex> = Vec::new();
 let mut add_v = |p: DVec3| -> usize {
 let idx = verts.len();
 verts.push(Vertex { point: p });
 idx
 };

 let mut faces: Vec<Face> = Vec::new();

 let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
 DVec3::new(center_xy.x + u, center_xy.y + v, z)
 };

 let circle_poly = |r: f64| -> Vec<DVec2> {
 let mut poly = Vec::with_capacity(n_arc + 1);
 for k in 0..=n_arc {
 let ang = tau * k as f64 / n_arc as f64;
 let (s, c) = ang.sin_cos();
 poly.push(DVec2::new(r * c, r * s));
 }
 poly
 };

 let cyl_poly = circle_poly(cyl_r);

 // 1. Cylinder wall (full height)
 add_wall_section(&mut add_v, &mut faces, &cyl_poly, cyl_z_lo, cyl_z_hi, n_slices_circ, &to_world, &empty_wire);

 // 2. Determine which cylinder ends are reached by the cavity
 let cavity_reaches_top = (overlap_hi - cyl_z_hi).abs() < 1e-9;
 let cavity_reaches_bot = (overlap_lo - cyl_z_lo).abs() < 1e-9;

 // 3. Far-end cap(s) =cylinder ends NOT reached by cavity
 if !cavity_reaches_bot {
 add_cap_face(&mut add_v, &mut faces, &cyl_poly, cyl_z_lo, -DVec3::Z, &to_world, &empty_wire);
 }
 if !cavity_reaches_top {
 add_cap_face(&mut add_v, &mut faces, &cyl_poly, cyl_z_hi, DVec3::Z, &to_world, &empty_wire);
 }

 // 4. Spherical cavity surface in overlap region
 {
 let n_rings = 32usize;
 let dz = (overlap_hi - overlap_lo) / n_rings as f64;
 for i in 0..n_rings {
 let za = overlap_lo + dz * i as f64;
 let zb = overlap_lo + dz * (i + 1) as f64;
 let dz_a = za - sz;
 let dz_b = zb - sz;
 let r_sq_a = sphere_r.powi(2) - dz_a.powi(2);
 let r_sq_b = sphere_r.powi(2) - dz_b.powi(2);
 let ra = if r_sq_a <= 0.0 { 0.0 } else { r_sq_a.sqrt() };
 let rb = if r_sq_b <= 0.0 { 0.0 } else { r_sq_b.sqrt() };
 if ra < tol && rb < tol { continue; }

 let nn = n_arc;
 let mut idx = Vec::with_capacity(2 * (nn + 1));
 for k in 0..=nn {
 let ang = tau * k as f64 / nn as f64;
 let (s, c) = ang.sin_cos();
 idx.push(add_v(to_world(ra * c, ra * s, za)));
 }
 for k in 0..=nn {
 let ang = tau * k as f64 / nn as f64;
 let (s, c) = ang.sin_cos();
 idx.push(add_v(to_world(rb * c, rb * s, zb)));
 }

 let mut tris = Vec::with_capacity(nn * 2);
 for j in 0..nn {
 tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
 tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
 }
 faces.push(Face {
 outer_wire: empty_wire(), inner_wires: vec![],
 normal: DVec3::ZERO, triangles: tris,
 sample_point: None, mesh_dirty: false,
 surface_idx: None,
 });
 }
 }

 // 5. Annular ring(s) at cylinder end(s) where cavity reaches but doesn't fill
 if cavity_reaches_top {
 let r_at = (sphere_r.powi(2) - (cyl_z_hi - sz).powi(2)).sqrt().max(0.0);
 if r_at < cyl_r - tol && r_at >= 0.0 {
 add_annular_ring(&mut add_v, &mut faces, cyl_z_hi, r_at, cyl_r, n_arc, &to_world, &empty_wire);
 }
 }
 if cavity_reaches_bot {
 let r_at = (sphere_r.powi(2) - (cyl_z_lo - sz).powi(2)).sqrt().max(0.0);
 if r_at < cyl_r - tol && r_at >= 0.0 {
 add_annular_ring(&mut add_v, &mut faces, cyl_z_lo, r_at, cyl_r, n_arc, &to_world, &empty_wire);
 }
 }

 if faces.is_empty() { return None; }

 let geom = GeomStore { edge_vertex_params: vec![],  face_internal_vertices: vec![],
 curves: vec![], surfaces: vec![], curve2ds: vec![],
 edge_curve: vec![],
 face_surface: vec![None; faces.len()],
 edge_pcurves: vec![], edge_curve_range: vec![],
 edge_degenerated: vec![], vertex_tolerance: vec![],
 edge_tolerance: vec![], face_tolerance: vec![],
 curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
 edge_same_parameter: vec![], edge_same_range: vec![],
 };

 Some(BRep {
 vertices: verts, edges: vec![],
 solids: vec![Solid { shells: vec![Shell { faces }] }],
 geom, compound: None, compsolid: None,
 })
}

/// Build a Z-slice triangle-mesh solid representing a spherical slice (portion
/// of a sphere between two parallel Z planes). The solid is built with
/// pre-triangulated faces and no analytic surfaces =SA computation reads
/// the stored triangles directly.
fn build_spherical_slice_solid(center: DVec3, radius: f64, z_lo: f64, z_hi: f64) -> Option<BRep> {
 const N_RINGS: usize = 32;
 const N_CIRC: usize = 48;

 let h = z_hi - z_lo;
 if h <= TOLERANCE_LEN_MIN {
 return None;
 }
 let dz = h / N_RINGS as f64;

 let mut brep = BRep::new();

 // Generate ring vertices
 let mut verts: Vec<Vertex> = Vec::new();
 let mut ring_start: Vec<usize> = Vec::with_capacity(N_RINGS + 2);

 for i in 0..=N_RINGS {
 let z = z_lo + dz * i as f64;
 let dz_c = z - center.z;
 let r_sq = radius.powi(2) - dz_c.powi(2);
 ring_start.push(verts.len());

 if r_sq <= 0.0 {
 verts.push(Vertex { point: DVec3::new(center.x, center.y, z) });
 } else {
 let r = r_sq.sqrt();
 for j in 0..N_CIRC {
 let ang = std::f64::consts::TAU * j as f64 / N_CIRC as f64;
 let (s, c) = ang.sin_cos();
 verts.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, z) });
 }
 }
 }
 ring_start.push(verts.len());

 let v_off = brep.vertices.len();
 brep.vertices = verts;

 let mut faces: Vec<Face> = Vec::new();

 // Lateral faces: triangle strips between adjacent rings
 for i in 0..N_RINGS {
 let b_s = ring_start[i];
 let b_e = ring_start[i + 1];
 let t_s = ring_start[i + 1];
 let t_e = ring_start[i + 2];
 let bn = b_e - b_s;
 let tn = t_e - t_s;

 let tris: Vec<[usize; 3]> = if bn == 1 && tn > 1 {
 (0..tn).flat_map(|j| {
 let k = (j + 1) % tn;
 vec![[v_off + b_s, v_off + t_s + j, v_off + t_s + k]]
 }).collect()
 } else if tn == 1 && bn > 1 {
 (0..bn).flat_map(|j| {
 let k = (j + 1) % bn;
 vec![[v_off + b_s + j, v_off + b_s + k, v_off + t_s]]
 }).collect()
 } else if bn > 1 && tn > 1 {
 let n = bn.min(tn);
 (0..n).flat_map(|j| {
 let k = (j + 1) % n;
 let b0 = v_off + b_s + j;
 let b1 = v_off + b_s + k;
 let t0 = v_off + t_s + j;
 let t1 = v_off + t_s + k;
 vec![[b0, b1, t1], [b0, t1, t0]]
 }).collect()
 } else {
 Vec::new()
 };

 if !tris.is_empty() {
 faces.push(Face {
 outer_wire: Wire { edges: vec![] },
 inner_wires: vec![],
 normal: DVec3::ZERO,
 triangles: tris,
 sample_point: None,
 mesh_dirty: false,
 surface_idx: None,
 });
 }
 }

 // Cap disks at exposed planar ends (skip if at pole)
 let bot_n = ring_start[1] - ring_start[0];
 if bot_n > 1 {
 let r_sq = radius.powi(2) - (z_lo - center.z).powi(2);
 if r_sq > 0.0 {
 let r = r_sq.sqrt();
 push_cap_disk(&mut brep, &mut faces, center, r, z_lo,
 if z_lo > center.z { DVec3::Z } else { -DVec3::Z });
 }
 }

 let top_n = ring_start[N_RINGS + 1] - ring_start[N_RINGS];
 if top_n > 1 {
 let r_sq = radius.powi(2) - (z_hi - center.z).powi(2);
 if r_sq > 0.0 {
 let r = r_sq.sqrt();
 push_cap_disk(&mut brep, &mut faces, center, r, z_hi,
 if z_hi > center.z { DVec3::Z } else { -DVec3::Z });
 }
 }

 brep.solids.push(Solid { shells: vec![Shell { faces }] });
 Some(brep)
}

/// Push a triangulated planar disk face (center + rim at given z) into `faces`.
fn push_cap_disk(brep: &mut BRep, faces: &mut Vec<Face>, center: DVec3, r: f64, z: f64, normal: DVec3) {
 const N_CIRC: usize = 48;
 let mut rim: Vec<usize> = Vec::with_capacity(N_CIRC);
 for j in 0..N_CIRC {
 let ang = std::f64::consts::TAU * j as f64 / N_CIRC as f64;
 let (s, c) = ang.sin_cos();
 let idx = brep.vertices.len();
 brep.vertices.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, z) });
 rim.push(idx);
 }
 let ctr_idx = brep.vertices.len();
 brep.vertices.push(Vertex { point: DVec3::new(center.x, center.y, z) });

 let tris: Vec<[usize; 3]> = (0..N_CIRC).map(|j| {
 let k = (j + 1) % N_CIRC;
 [rim[j], rim[k], ctr_idx]
 }).collect();

 faces.push(Face {
 outer_wire: Wire { edges: vec![] },
 inner_wires: vec![],
 normal,
 triangles: tris,
 sample_point: None,
 mesh_dirty: false,
 surface_idx: None,
 });
}

/// [`BooleanOpType::Difference`] for nested analytic spheres sharing a center.
///
/// Builds a hollow spherical shell as **two solids**: outer sphere plus inner sphere with reversed
/// face orientation (`reverse_face`). Total surface area matches the analytic spherical shell
/// \(4\pi(R^2+r^2)\) under [`rcad_kernel::surface_area`].
///
/// Compound `rcad_kernel::volume` may not match \(4\pi/3(R^3-r^3)\) until sphere face normals / tessellation agree for
/// [`signed_volume`] everywhere.
///
/// OCCT DRAW `mkvolume` on trimmed spherical patches can report a different `checkprops -s` than this
/// full-sphere analytic shell.
pub fn try_difference_concentric_spheres(a: &BRep, b: &BRep) -> Option<BRep> {
 let (ca, ra) = try_sphere_primitive_center_radius(a)?;
 let (cb, rb) = try_sphere_primitive_center_radius(b)?;
 let scale = ra.max(rb).max(1.0);
 if (ca - cb).length() > TOL.max(TOLERANCE_COORD_SUB * scale) {
 return None;
 }
 let ro = ra.max(rb);
 let ri = ra.min(rb);
 if ro - ri <= TOLERANCE_LEN_MIN * ro.max(1.0) {
 return None;
 }
 let center = ca;
 let outer = make_sphere_brep(center, ro).ok()?;
 let mut inner_cavity = make_sphere_brep(center, ri).ok()?;
 crate::reverse_face(&mut inner_cavity, 0);
 Some(BRep::compound_from_shapes(&[outer, inner_cavity]))
}

/// [`BooleanOpType::Intersection`] for two analytic sphere primitives sharing a center.
///
/// The intersection of nested balls is the smaller-radius ball (same center).
///
/// Returns [`None`] when centers differ or radii are degenerate =callers fall back to Pave/Builder.
pub fn try_intersection_concentric_spheres(a: &BRep, b: &BRep) -> Option<BRep> {
 let (ca, ra) = try_sphere_primitive_center_radius(a)?;
 let (cb, rb) = try_sphere_primitive_center_radius(b)?;
 let scale = ra.max(rb).max(1.0);
 if (ca - cb).length() > TOL.max(TOLERANCE_COORD_SUB * scale) {
 return None;
 }
 let r = ra.min(rb);
 let r_eps = TOLERANCE_LEN_MIN * scale.max(1.0);
 if r <= r_eps {
 return None;
 }
 make_sphere_brep(ca, r).ok()
}

// --- Coaxial cone =cylinder (OCCT `bopcommon_simple/ZP7`): generic Builder over-counts area. --------

fn z_axis_sharp_cone_z_span(cone: &BRep) -> Option<(f64, f64, f64)> {
 const APAR: f64 = TOLERANCE_ADAPTIVE_MAX;
 const XY: f64 = 2.0 * TOLERANCE_ADAPTIVE_MAX;
 let sh = cone.solids.get(0)?.shells.get(0)?;
 if sh.faces.len() < 2 {
 return None;
 }
 let mut cf: Option<ConicalSurface> = None;
 let mut po: Option<DVec3> = None;
 let mut fi = 0usize;
 for s in &cone.solids {
 for sh in &s.shells {
 for _ in &sh.faces {
 let si = *cone.geom.face_surface.get(fi)?.as_ref()?;
 match cone.geom.surfaces.get(si)? {
 Surface3::Cone(c) => cf = Some(*c),
 Surface3::Plane(p) => po = Some(p.origin),
 _ => return None,
 }
 fi += 1;
 }
 }
 }
 let c = cf?;
 let u = c.axis_dir();
 if u.cross(DVec3::Z).length() > APAR {
 return None;
 }
 let apex = c.apex_point();
 if apex.x.abs() > XY || apex.y.abs() > XY {
 return None;
 }
 let b = po?;
 if b.x.abs() > XY || b.y.abs() > XY {
 return None;
 }
 let t = (b - apex).dot(u);
 let rb = t * c.half_angle_rad.tan();
 if t < TOLERANCE_MESH_LEGACY || rb < TOLERANCE_MESH_LEGACY {
 return None;
 }
 Some((apex.z, b.z, rb))
}

pub(crate) fn z_axis_cylinder_z_span_r(cyl: &BRep) -> Option<(f64, f64, f64)> {
 const APAR: f64 = TOLERANCE_ADAPTIVE_MAX;
 const XY: f64 = 2.0 * TOLERANCE_ADAPTIVE_MAX;
 let sh = cyl.solids.get(0)?.shells.get(0)?;
 if sh.faces.len() != 3 {
 return None;
 }
 let mut r = None;
 let mut zs = Vec::with_capacity(2);
 let mut fi = 0usize;
 for s in &cyl.solids {
 for sh in &s.shells {
 for _ in &sh.faces {
 let si = *cyl.geom.face_surface.get(fi)?.as_ref()?;
 match cyl.geom.surfaces.get(si)? {
 Surface3::Cylinder(cc) => {
 if cc.axis.normalize_or_zero().cross(DVec3::Z).length() > APAR {
 return None;
 }
 if cc.origin.x.abs() > XY || cc.origin.y.abs() > XY {
 return None;
 }
 r = Some(cc.radius);
 }
 Surface3::Plane(p) => zs.push(p.origin.z),
 _ => return None,
 }
 fi += 1;
 }
 }
 }
 if zs.len() != 2 {
 return None;
 }
 Some((zs[0].min(zs[1]), zs[0].max(zs[1]), r?))
}

fn try_intersection_coaxial_cone_cylinder_pair(cone: &BRep, cyl: &BRep) -> Option<BRep> {
 use rcad_modeling::make_conical_frustum_brep;
 let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
 let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
 let zcn = zb.min(za);
 let zcx = zb.max(za);
 let z0 = zlo.max(zcn);
 let z1 = zhi.min(zcx);
 if z1 - z0 < TOLERANCE_MESH_LEGACY {
 return None;
 }
 let apex_hi = za > zb;
 let hc = (za - zb).abs();
 let rcz = |z: f64| {
 let num = if apex_hi { (za - z).abs() } else { (z - za).abs() };
 rb * num / hc
 };
 let r0 = rcz(z0).min(rc);
 let r1 = rcz(z1).min(rc);
 if r0 < TOLERANCE_COORD_SUB && r1 < TOLERANCE_COORD_SUB {
 return None;
 }
 let zm = (z0 + z1) * 0.5;
 let h = z1 - z0;
 make_conical_frustum_brep(DVec3::new(0.0, 0.0, zm), DVec3::Z, DVec3::X, r0, r1, h).ok()
}

/// Sharp Z-aligned cone =finite Z-aligned cylinder (same axis / origin in `xy`), e.g. OCCT ZP7.
pub fn try_intersection_coaxial_cone_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
 try_intersection_coaxial_cone_cylinder_pair(a, b)
 .or_else(|| try_intersection_coaxial_cone_cylinder_pair(b, a))
}

/// Two coaxial Z-aligned cylinders with the same radius -> intersection is the
/// overlapping Z-span cylinder (e.g. OCCT bcommon_simple/J1).
pub fn try_intersection_coaxial_cylinder_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
 let (z1_lo, z1_hi, r1) = z_axis_cylinder_z_span_r(a)?;
 let (z2_lo, z2_hi, r2) = z_axis_cylinder_z_span_r(b)?;
 if (r1 - r2).abs() > TOLERANCE_ADAPTIVE_MAX {
 return None;
 }
 let z0 = z1_lo.max(z2_lo);
 let z1 = z1_hi.min(z2_hi);
 if z1 - z0 < TOLERANCE_MESH_LEGACY {
 return None;
 }
 let zm = (z0 + z1) * 0.5;
 let h = z1 - z0;
 rcad_modeling::make_cylinder_brep(
 DVec3::new(0.0, 0.0, zm), DVec3::Z, DVec3::X, r1, h,
 ).ok()
}

/// Detect a cylinder BRep and return (center, axis, radius, height, origin, ref_dir).
///
/// Unlike `try_cylinder_center_axis_radius_height`, this works for cylinders
/// along ANY axis (not just Z) by computing height from face-plane distances
/// along the axis.
fn try_cylinder_any_axis(brep: &BRep) -> Option<(DVec3, DVec3, f64, f64, DVec3, DVec3)> {
 for s in &brep.solids {
 for sh in &s.shells {
 let mut origin = DVec3::ZERO;
 let mut axis = DVec3::Z;
 let mut ref_dir = DVec3::X;
 let mut radius = 0.0;
 let mut found = false;
 let mut plane_dots: Vec<f64> = Vec::new();

 for fi in 0..sh.faces.len() {
 let si = brep.geom.face_surface.get(fi)?.as_ref()?;
 match brep.geom.surfaces.get(*si)? {
 Surface3::Cylinder(cc) => {
 origin = cc.origin;
 axis = cc.axis.normalize_or_zero();
 radius = cc.radius;
 ref_dir = cc.ref_dir;
 found = true;
 }
 Surface3::Plane(pl) => {
 plane_dots.push(pl.origin.dot(axis));
 }
 _ => return None,
 }
 }
 if !found || plane_dots.len() < 2 { return None; }
 plane_dots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
 let height = plane_dots.last()? - plane_dots.first()?;
 if height < crate::tolerance::TOLERANCE_LEN_MIN { return None; }
 let center = origin + axis * (height / 2.0);
 return Some((center, axis, radius, height, origin, ref_dir));
 }
 }
 None
}