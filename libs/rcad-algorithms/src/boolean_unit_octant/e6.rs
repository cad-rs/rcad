
/// Generate a quad-grid triangulation between two Z-aligned rings (conical lateral strip).
fn wall_grid(
 add_v: &mut impl FnMut(DVec3) -> usize,
 bot: &[DVec3],
 top: &[DVec3],
 tris: &mut Vec<[usize; 3]>,
) {
 let n = bot.len().min(top.len());
 if n < 3 { return; }
 let mut idx = Vec::with_capacity(2 * n);
 for i in 0..n {
 idx.push(add_v(bot[i]));
 }
 for i in 0..n {
 idx.push(add_v(top[i]));
 }
 for i in 0..n {
 let j = (i + 1) % n;
 let b0 = idx[i];
 let b1 = idx[j];
 let t0 = idx[n + i];
 let t1 = idx[n + j];
 tris.push([b0, b1, t1]);
 tris.push([b0, t1, t0]);
 }
}

/// Triangle fan for a full disk centered at origin in XY plane at given z.
fn disk_tri_fan(
 add_v: &mut impl FnMut(DVec3) -> usize,
 center: DVec3,
 radius: f64,
 n: usize,
 tris: &mut Vec<[usize; 3]>,
) {
 let c_idx = add_v(center);
 let mut ring = Vec::with_capacity(n);
 for i in 0..n {
 let ang = std::f64::consts::TAU * i as f64 / n as f64;
 let (s, c) = ang.sin_cos();
 ring.push(add_v(DVec3::new(radius * c, radius * s, center.z)));
 }
 for i in 0..n {
 let j = (i + 1) % n;
 tris.push([c_idx, ring[i], ring[j]]);
 }
}

/// Triangle fan for an annulus centered at origin in XY plane at given z.
fn annulus_tri_fan(
 add_v: &mut impl FnMut(DVec3) -> usize,
 center: DVec3,
 r_inner: f64,
 r_outer: f64,
 n: usize,
 tris: &mut Vec<[usize; 3]>,
) {
 let mut outer_ring = Vec::with_capacity(n);
 let mut inner_ring = Vec::with_capacity(n);
 for i in 0..n {
 let ang = std::f64::consts::TAU * i as f64 / n as f64;
 let (s, c) = ang.sin_cos();
 outer_ring.push(add_v(DVec3::new(r_outer * c, r_outer * s, center.z)));
 inner_ring.push(add_v(DVec3::new(r_inner * c, r_inner * s, center.z)));
 }
 for i in 0..n {
 let j = (i + 1) % n;
 tris.push([outer_ring[i], outer_ring[j], inner_ring[j]]);
 tris.push([outer_ring[i], inner_ring[j], inner_ring[i]]);
 }
}


// = =  Box= ylinder Difference Fast Path (box =cylinder, Z-axis cylinder) = = = = = = = 

/// A point where the circle intersects a box edge in UV space.
#[derive(Debug, Clone, Copy)]
struct UVEdgePt {
 /// Perimeter parameter t =[0, 8)
 t: f64,
 /// UV coordinates
 u: f64,
 v: f64,
 /// Edge index: 0 = u_min, 1 = u_max, 2 = v_min, 3 = v_max
 edge: usize,
 /// Circle angle  ?= atan2(v =cv, u =cu)
 theta: f64,
}

/// Find which of the 3 box axes is closest to the world Z axis.
pub(crate) fn find_z_axis_index(info: &BoxInfo) -> Option<usize> {
 for i in 0..3 {
 if info.axes[i].dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN {
 return Some(i);
 }
 }
 None
}

/// Extract cylinder parameters from a cylinder primitive.
pub(crate) fn try_cylinder_center_axis_radius_height(brep: &BRep) -> Option<(DVec3, DVec3, f64, f64)> {
 let Some(shell) = brep.solids.first()?.shells.first() else { return None };
 let mut center = DVec3::ZERO;
 let mut axis = DVec3::Z;
 let mut radius = 0.0;
 let mut height = 0.0;
 let mut found = false;
 for fi in 0..shell.faces.len() {
 let Some(Some(si)) = brep.geom.face_surface.get(fi) else { continue };
 let Some(surf) = brep.geom.surfaces.get(*si) else { continue };
 if let Surface3::Cylinder(cc) = surf {
 axis = cc.axis.normalize_or_zero();
 center = cc.origin;
 radius = cc.radius;
 found = true;
 }
 }
 if !found { return None; }
 let mut z_vals = Vec::new();
 for fi in 0..shell.faces.len() {
 let Some(Some(si)) = brep.geom.face_surface.get(fi) else { continue };
 let Some(surf) = brep.geom.surfaces.get(*si) else { continue };
 if let Surface3::Plane(pl) = surf {
 z_vals.push(pl.origin.z);
 }
 }
 if z_vals.len() < 2 { return None; }
 let z_lo = z_vals.iter().min_by(|a,b| a.partial_cmp(b).unwrap()).copied()?;
 let z_hi = z_vals.iter().max_by(|a,b| a.partial_cmp(b).unwrap()).copied()?;
 height = z_hi - z_lo;
 Some((center, axis, radius, height))
}

/// Box perimeter in UV space (CCW from v_min).
/// t =[0,2) =v_min (v== v), u= = u,eu]
/// t =[2,4) =u_max (u=eu),  v= = v,ev]
/// t =[4,6) =v_max (v=ev),  u= eu,= u]
/// t =[6,8) =u_min (u== u), v= ev,= v]
fn box_perimeter_uv(t: f64, eu: f64, ev: f64) -> (f64, f64) {
 let tn = t.rem_euclid(8.0);
 if tn < 2.0 {
 let s = tn / 2.0;
 (-eu + 2.0 * eu * s, -ev)
 } else if tn < 4.0 {
 let s = (tn - 2.0) / 2.0;
 (eu, -ev + 2.0 * ev * s)
 } else if tn < 6.0 {
 let s = (tn - 4.0) / 2.0;
 (eu - 2.0 * eu * s, ev)
 } else {
 let s = (tn - 6.0) / 2.0;
 (-eu, ev - 2.0 * ev * s)
 }
}

/// Map UV coordinate to perimeter t =[0, 8).
fn uv_to_perimeter_t(u: f64, v: f64, eu: f64, ev: f64) -> f64 {
 if v <= -ev + TOLERANCE_LEN_MIN {
 ((u + eu) / (2.0 * eu).max(TOLERANCE_LEN_MIN)) * 2.0
 } else if u >= eu - TOLERANCE_LEN_MIN {
 2.0 + ((v + ev) / (2.0 * ev).max(TOLERANCE_LEN_MIN)) * 2.0
 } else if v >= ev - TOLERANCE_LEN_MIN {
 4.0 + ((eu - u) / (2.0 * eu).max(TOLERANCE_LEN_MIN)) * 2.0
 } else {
 6.0 + ((ev - v) / (2.0 * ev).max(TOLERANCE_LEN_MIN)) * 2.0
 }
}

/// Find circle-box edge intersections in UV space.
fn circle_rect_intersections_uv(cu: f64, cv: f64, r: f64, eu: f64, ev: f64) -> Vec<UVEdgePt> {
 let mut pts = Vec::new();
 let tol = TOLERANCE_LEN_MIN;

 let add_if = |u: f64, v: f64, edge: usize, pts: &mut Vec<UVEdgePt>| {
 if u >= -eu - tol && u <= eu + tol && v >= -ev - tol && v <= ev + tol {
 let u_cl = u.clamp(-eu, eu);
 let v_cl = v.clamp(-ev, ev);
 let t = uv_to_perimeter_t(u_cl, v_cl, eu, ev);
 let theta = (v_cl - cv).atan2(u_cl - cu);
 pts.push(UVEdgePt { t, u: u_cl, v: v_cl, edge, theta });
 }
 };

 let d = -ev - cv; let disc = r * r - d * d;
 if disc >= 0.0 { let off = disc.sqrt(); add_if(cu + off, -ev, 2, &mut pts); if off > tol { add_if(cu - off, -ev, 2, &mut pts); } }
 let d = ev - cv; let disc = r * r - d * d;
 if disc >= 0.0 { let off = disc.sqrt(); add_if(cu + off, ev, 3, &mut pts); if off > tol { add_if(cu - off, ev, 3, &mut pts); } }
 let d = -eu - cu; let disc = r * r - d * d;
 if disc >= 0.0 { let off = disc.sqrt(); add_if(-eu, cv + off, 0, &mut pts); if off > tol { add_if(-eu, cv - off, 0, &mut pts); } }
 let d = eu - cu; let disc = r * r - d * d;
 if disc >= 0.0 { let off = disc.sqrt(); add_if(eu, cv + off, 1, &mut pts); if off > tol { add_if(eu, cv - off, 1, &mut pts); } }

 pts.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
 pts.dedup_by(|a, b| (a.t - b.t).abs() < tol && (a.u - b.u).abs() < tol && (a.v - b.v).abs() < tol);
 pts
}

/// Create a planar BRep face from 4 corner points and a given normal.
fn rect_face_4(corners: [DVec3; 4], normal: DVec3) -> Option<BRep> {
 let mut brep = BRep::default();
 let surface = Surface3::Plane(Plane { origin: corners[0], normal });
 let vs: Vec<usize> = corners.iter().map(|p| make_vertex(&mut brep, *p)).collect();
 let mut wes = Vec::with_capacity(4);
 for i in 0..4 {
 let j = (i + 1) % 4;
 let dir = (corners[j] - corners[i]).normalize();
 let len = (corners[j] - corners[i]).length();
 if len < TOLERANCE_LEN_MIN { return None; }
 let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: corners[i], direction: dir }), 0.0, len, vs[i], vs[j]).ok()?;
 wes.push(WireEdge::new(ei, true));
 }
 let _fi = make_face(&mut brep, surface, make_wire(wes), vec![]).ok()?;
 Some(brep)
}

/// Create a planar face from a polygon with a given normal.
fn planar_face_from_polygon(poly: &[DVec3], normal: DVec3) -> Option<BRep> {
 if poly.len() < 3 { return None; }
 let mut brep = BRep::default();
 let surface = Surface3::Plane(Plane { origin: poly[0], normal });
 let vs: Vec<usize> = poly.iter().map(|p| make_vertex(&mut brep, *p)).collect();
 let mut wes = Vec::with_capacity(poly.len());
 let n = poly.len();
 for i in 0..n {
 let j = (i + 1) % n;
 let dir = (poly[j] - poly[i]).normalize();
 let len = (poly[j] - poly[i]).length();
 if len < TOLERANCE_LEN_MIN { return None; }
 let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: poly[i], direction: dir }), 0.0, len, vs[i], vs[j]).ok()?;
 wes.push(WireEdge::new(ei, true));
 }
 let _fi = make_face(&mut brep, surface, make_wire(wes), vec![]).ok()?;
 Some(brep)
}

/// Create a planar face from an outer polygon and an inner polygon (hole).
/// `outer` must be CCW when viewed along `normal`.
/// `inner` must be CW when viewed along `normal`.
fn planar_face_with_inner_hole(outer: &[DVec3], inner: &[DVec3], normal: DVec3) -> Option<BRep> {
 if outer.len() < 3 || inner.len() < 3 { return None; }
 let mut brep = BRep::default();
 let surface = Surface3::Plane(Plane { origin: outer[0], normal });
 let vs: Vec<usize> = outer.iter().map(|p| make_vertex(&mut brep, *p)).collect();
 let mut outer_wes = Vec::with_capacity(outer.len());
 let n = outer.len();
 for i in 0..n {
 let j = (i + 1) % n;
 let dir = (outer[j] - outer[i]).normalize();
 let len = (outer[j] - outer[i]).length();
 if len < TOLERANCE_LEN_MIN { return None; }
 let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: outer[i], direction: dir }), 0.0, len, vs[i], vs[j]).ok()?;
 outer_wes.push(WireEdge::new(ei, true));
 }
 let inner_vs: Vec<usize> = inner.iter().map(|p| make_vertex(&mut brep, *p)).collect();
 let mut inner_wes = Vec::with_capacity(inner.len());
 let m = inner.len();
 for i in 0..m {
 let j = (i + 1) % m;
 let dir = (inner[j] - inner[i]).normalize();
 let len = (inner[j] - inner[i]).length();
 if len < TOLERANCE_LEN_MIN { return None; }
 let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: inner[i], direction: dir }), 0.0, len, inner_vs[i], inner_vs[j]).ok()?;
 inner_wes.push(WireEdge::new(ei, true));
 }
 let _fi = make_face(&mut brep, surface, make_wire(outer_wes), vec![make_wire(inner_wes)]).ok()?;
 Some(brep)
}

// = =  shared segment types and helpers = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// A merged segment along the box perimeter =either inside or outside the cylinder.
struct MergedSeg {
 t0: f64,
 t1: f64,
 outside: bool,
}

/// Compute merged segments along the box perimeter from intersection points.
fn compute_merged_segments(ints: &[UVEdgePt], cu: f64, cv: f64, cyl_r: f64, eu: f64, ev: f64) -> Vec<MergedSeg> {
 let tol = TOLERANCE_LEN_MIN;
 let r2 = cyl_r * cyl_r;

 let outside_at = |t: f64| -> bool {
 let (u, v) = box_perimeter_uv(t, eu, ev);
 (u - cu).powi(2) + (v - cv).powi(2) > r2 + tol
 };

 if ints.is_empty() {
 return vec![MergedSeg { t0: 0.0, t1: 8.0, outside: outside_at(0.0) }];
 }

 let mut t_vals: Vec<f64> = ints.iter().map(|p| p.t).collect();
 t_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
 t_vals.dedup_by(|a, b| (*a - *b).abs() < tol);

 let mut segs: Vec<MergedSeg> = Vec::new();
 let mut prev_t = 0.0;
 for &t in &t_vals {
 if t <= prev_t + tol { continue; }
 let mid_t = (prev_t + t) / 2.0;
 segs.push(MergedSeg { t0: prev_t, t1: t, outside: outside_at(mid_t) });
 prev_t = t;
 }
 if prev_t < 8.0 - tol {
 let mid_t = (prev_t + 8.0) / 2.0;
 segs.push(MergedSeg { t0: prev_t, t1: 8.0, outside: outside_at(mid_t) });
 }

 let mut merged: Vec<MergedSeg> = Vec::new();
 for s in segs {
 if let Some(last) = merged.last_mut() {
 if last.outside == s.outside {
 last.t1 = s.t1;
 continue;
 }
 }
 merged.push(s);
 }
 merged
}

/// Generate evenly-spaced points on a circular arc, choosing the direction
/// whose midpoint stays inside the box rect.  Returns vertices in world
/// coordinates via `corner(u, v, z)`, including both endpoints.
fn arc_vertices(
 cu: f64, cv: f64, r: f64,
 th_start: f64, th_end: f64,
 eu: f64, ev: f64, z: f64,
 corner: &impl Fn(f64, f64, f64) -> DVec3,
) -> Vec<DVec3> {
 let tol = TOLERANCE_LEN_MIN;
 let inside_box = |th: f64| -> bool {
 let u = cu + r * th.cos();
 let v = cv + r * th.sin();
 u >= -eu - tol && u <= eu + tol && v >= -ev - tol && v <= ev + tol
 };

 let mut cw_len = th_start - th_end;
 if cw_len < 0.0 { cw_len += 2.0 * std::f64::consts::PI; }
 let ccw_len = 2.0 * std::f64::consts::PI - cw_len;

 let cw_mid = th_start - cw_len / 2.0;
 let ccw_mid = th_start + ccw_len / 2.0;
 let cw_in = inside_box(cw_mid);
 let ccw_in = inside_box(ccw_mid);
 let (dtheta, arc_len) = if cw_in && !ccw_in {
 (-cw_len, cw_len)
 } else if ccw_in && !cw_in {
 (ccw_len, ccw_len)
 } else if cw_in && ccw_in {
 // Both midpoints are inside the box. The shorter arc may still go outside
 // the box at non-midpoint positions (e.g., an arc >180 ?whose midpoint
 // happens to be inside but passes through a bulge).  Sample several points
 // along the shorter arc to verify all lie inside; if not, take the longer.
 let shorter_cw = cw_len <= ccw_len;
 let (check_len, check_dir) = if shorter_cw {
 (cw_len, -1.0)  // CW direction = decreasing  ?
 } else {
 (ccw_len, 1.0)  // CCW direction = increasing  ?
 };
 // Check 5 evenly-spaced sample points along the shorter arc.
 let n_samples = 5usize;
 let all_inside = (0..=n_samples).all(|i| {
 let frac = i as f64 / n_samples as f64;
 let th = th_start + check_dir * check_len * frac;
 // Normalize to [0, 2 ?
 let th_n = th.rem_euclid(2.0 * std::f64::consts::PI);
 inside_box(th_n)
 });
 if all_inside {
 // Shorter arc stays inside =use it.
 if shorter_cw { (-cw_len, cw_len) } else { (ccw_len, ccw_len) }
 } else {
 // Shorter arc exits the box =use the longer arc.
 if shorter_cw { (ccw_len, ccw_len) } else { (-cw_len, cw_len) }
 }
 } else if cw_len <= ccw_len {
 (-cw_len, cw_len)
 } else {
 (ccw_len, ccw_len)
 };
 let steps = (arc_len / 0.08).ceil() as usize;
 let steps = steps.max(2).min(200);

 let mut verts = Vec::with_capacity(steps + 1);
 verts.push(corner(cu + r * th_start.cos(), cv + r * th_start.sin(), z));
 for i in 1..steps {
 let th = th_start + dtheta * (i as f64 / steps as f64);
 let pt = corner(cu + r * th.cos(), cv + r * th.sin(), z);
 if (verts.last().unwrap() - pt).length() > tol {
 verts.push(pt);
 }
 }
 let pt_end = corner(cu + r * th_end.cos(), cv + r * th_end.sin(), z);
 if (verts.last().unwrap() - pt_end).length() > tol {
 verts.push(pt_end);
 }
 verts
}

/// Build cap polygon at z level for (box rect =circle).
fn build_cap_polygon(
 merged: &[MergedSeg], cu: f64, cv: f64, cyl_r: f64,
 eu: f64, ev: f64, z: f64,
 corner: &impl Fn(f64, f64, f64) -> DVec3,
) -> Vec<DVec3> {
 let tol = TOLERANCE_LEN_MIN;

 if merged.is_empty() || (merged.len() == 1 && !merged[0].outside) {
 return vec![];
 }

 let mut poly = Vec::new();
 for seg in merged {
 if seg.t1 <= seg.t0 + tol { continue; }
 if seg.outside {
 add_box_perimeter_vertices(seg.t0, seg.t1, eu, ev, z, corner, &mut poly);
 } else {
 let (pu, pv) = box_perimeter_uv(seg.t0, eu, ev);
 let (nu, nv) = box_perimeter_uv(seg.t1, eu, ev);
 let th_prev = (pv - cv).atan2(pu - cu);
 let th_next = (nv - cv).atan2(nu - cu);
 add_circle_arc_vertices(cu, cv, cyl_r, th_prev, th_next, eu, ev, z, corner, &mut poly);
 }
 }

 if poly.len() >= 2 && (poly.last().unwrap() - poly[0]).length() < tol { poly.pop(); }
 poly
}

/// Add box perimeter vertices between two t values.
fn add_box_perimeter_vertices(
 t_start: f64, t_end: f64, eu: f64, ev: f64, z: f64,
 corner: &impl Fn(f64, f64, f64) -> DVec3,
 poly: &mut Vec<DVec3>,
) {
 let tol = TOLERANCE_LEN_MIN;
 if poly.is_empty() || ((t_start - 0.0).abs() < tol || (t_start - 8.0).abs() < tol) {
 let (u0, v0) = box_perimeter_uv(t_start, eu, ev);
 let pt = corner(u0, v0, z);
 if poly.is_empty() || (poly.last().unwrap() - pt).length() > tol {
 poly.push(pt);
 }
 }

 let corner_ts = [0.0, 2.0, 4.0, 6.0];
 let start_norm = if t_start < 0.0 { t_start + 8.0 } else { t_start };
 let end_norm = if t_end <= t_start { t_end + 8.0 } else { t_end };

 // Collect corners within range and sort by normalized t for CCW order.
 let mut corners_in_range: Vec<(f64, f64, f64)> = Vec::new();
 for &ct in &corner_ts {
 let ct_norm = if ct < start_norm { ct + 8.0 } else { ct };
 if ct_norm > start_norm + tol && ct_norm < end_norm - tol {
 let (u, v) = box_perimeter_uv(ct, eu, ev);
 corners_in_range.push((ct_norm, u, v));
 }
 }
 corners_in_range.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
 for (_, u, v) in corners_in_range {
 let pt = corner(u, v, z);
 if (poly.last().unwrap() - pt).length() > tol { poly.push(pt); }
 }

 let (ue, ve) = box_perimeter_uv(t_end, eu, ev);
 let pt = corner(ue, ve, z);
 if (poly.last().unwrap() - pt).length() > tol { poly.push(pt); }
}

/// Add circular arc vertices for the cap polygon, taking direction from
/// `arc_vertices` and skipping the first vertex (already in poly from the
/// preceding box perimeter segment).
fn add_circle_arc_vertices(
 cu: f64, cv: f64, r: f64,
 th_start: f64, th_end: f64,
 eu: f64, ev: f64, z: f64,
 corner: &impl Fn(f64, f64, f64) -> DVec3,
 poly: &mut Vec<DVec3>,
) {
 let verts = arc_vertices(cu, cv, r, th_start, th_end, eu, ev, z, corner);
 // verts[0] is the start point =already in poly from the previous segment,
 // or if poly is empty (first merged segment is a circle arc, not box perimeter)
 // we need to push it here.
 if poly.is_empty() {
 if let Some(&first) = verts.first() {
 poly.push(first);
 }
 }
 for i in 1..verts.len() {
 if (poly.last().unwrap() - verts[i]).length() > TOLERANCE_LEN_MIN {
 poly.push(verts[i]);
 }
 }
}

/// Build side wall pieces trimmed at cylinder intersection parameters.
fn build_trimmed_edge_pieces(
 p_min: f64, p_max: f64,
 p_lo: f64, p_hi: f64,
 z_lo: f64, z_hi: f64, ew: f64,
 corner: &dyn Fn(f64, f64) -> DVec3,
 normal: DVec3,
 pieces: &mut Vec<BRep>,
) {
 let tol = TOLERANCE_LEN_MIN;
 let push = |c: [DVec3; 4], n: DVec3, p: &mut Vec<BRep>| { if let Some(f) = rect_face_4(c, n) { p.push(f); } };

 let strip = |s_min: f64, s_max: f64, pieces: &mut Vec<BRep>| {
 if s_max <= s_min + tol { return; }
 if z_lo > -ew + tol { push([corner(s_min, -ew), corner(s_max, -ew), corner(s_max, z_lo), corner(s_min, z_lo)], normal, pieces); }
 if z_hi < ew - tol { push([corner(s_min, z_hi), corner(s_max, z_hi), corner(s_max, ew), corner(s_min, ew)], normal, pieces); }
 if z_hi > z_lo + tol { push([corner(s_min, z_lo), corner(s_max, z_lo), corner(s_max, z_hi), corner(s_min, z_hi)], normal, pieces); }
 };

 if p_lo <= p_min + tol && p_hi >= p_max - tol { strip(p_min, p_max, pieces); return; }
 if p_lo > p_min + tol { strip(p_min, p_lo, pieces); }
 if p_hi < p_max - tol { strip(p_hi, p_max, pieces); }
}

/// Build quad faces for the cylindrical wall from merged perimeter segments,
/// using the same angular step as the cap polygon for matching vertices.
fn build_cylindrical_wall_from_segs(
 merged: &[MergedSeg],
 cu: f64, cv: f64, cyl_r: f64,
 z_lo: f64, z_hi: f64, eu: f64, ev: f64,
 u_ax: DVec3, v_ax: DVec3, c: DVec3,
 corner: &impl Fn(f64, f64, f64) -> DVec3,
 pieces: &mut Vec<BRep>,
) {
 let _tol = TOLERANCE_LEN_MIN;
 for seg in merged {
 if seg.outside { continue; }
 let (pu, pv) = box_perimeter_uv(seg.t0, eu, ev);
 let (nu, nv) = box_perimeter_uv(seg.t1, eu, ev);
 let th_prev = (pv - cv).atan2(pu - cu);
 let th_next = (nv - cv).atan2(nu - cu);

 let lo_verts = arc_vertices(cu, cv, cyl_r, th_prev, th_next, eu, ev, z_lo, corner);
 let hi_verts = arc_vertices(cu, cv, cyl_r, th_prev, th_next, eu, ev, z_hi, corner);

 let n = lo_verts.len().min(hi_verts.len());
 if n < 2 { continue; }

 for i in 0..n - 1 {
 let b0 = lo_verts[i];
 let b1 = lo_verts[i + 1];
 let t1 = hi_verts[i + 1];
 let t0 = hi_verts[i];

 // Compute outward radial direction for sign convention
 let mid = (b0 + b1 + t1 + t0) / 4.0;
 let u_mid = (mid - c).dot(u_ax);
 let v_mid = (mid - c).dot(v_ax);
 let radial_u = u_mid - cu;
 let radial_v = v_mid - cv;
 let outward = (u_ax * radial_u + v_ax * radial_v).normalize();
 let n_vec = (t0 - b0).cross(b1 - b0).normalize();
 let n_final = if n_vec.dot(outward) > 0.0 { n_vec } else { -n_vec };

 if let Some(f) = rect_face_4([b0, b1, t1, t0], n_final) { pieces.push(f); }
 }
 }
}

/// Main Stage 3 orchestrator: build box =cylinder with partial XY containment.
fn build_box_cylinder_result_partial(
 c: DVec3, u_ax: DVec3, v_ax: DVec3, eu: f64, ev: f64, ew: f64,
 cu: f64, cv: f64, cyl_r: f64, cyl_z_lo: f64, cyl_z_hi: f64,
) -> Option<BRep> {
 let tol = TOLERANCE_LEN_MIN;
 let z_lo = (cyl_z_lo - c.z).max(-ew);
 let z_hi = (cyl_z_hi - c.z).min(ew);
 if z_hi <= z_lo + tol { return None; }

 let corner_z = |u: f64, v: f64, z: f64| -> DVec3 { c + u*u_ax + v*v_ax + z*DVec3::Z };

 let ints = circle_rect_intersections_uv(cu, cv, cyl_r, eu, ev);
 let merged = compute_merged_segments(&ints, cu, cv, cyl_r, eu, ev);
 let cap_bottom = build_cap_polygon(&merged, cu, cv, cyl_r, eu, ev, z_lo, &corner_z);
 let cap_top = build_cap_polygon(&merged, cu, cv, cyl_r, eu, ev, z_hi, &corner_z);

 let mut pieces: Vec<BRep> = Vec::new();

 // Compute per-edge intersection parameters directly (NOT from the global
 // intersection list, which loses corner-shared entries during dedup).
 let add_pt = |v: f64, lo: f64, hi: f64, out: &mut Vec<f64>| {
 if v >= lo - tol && v <= hi + tol { out.push(v.clamp(lo, hi)); }
 };

 // u_min edge (u=-eu): solve (-eu-cu)^2 + (v-cv)^2 = r^2 =v
 let mut ints_u_min: Vec<f64> = Vec::with_capacity(2);
 let d0 = -eu - cu; let disc0 = cyl_r * cyl_r - d0 * d0;
 if disc0 >= 0.0 { let off = disc0.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_min); add_pt(cv+off, -ev, ev, &mut ints_u_min); }
 ints_u_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
 ints_u_min.dedup_by(|a,b| (*a - *b).abs() < tol);

 // u_max edge (u=eu): solve (eu-cu)^2 + (v-cv)^2 = r^2 =v
 let mut ints_u_max: Vec<f64> = Vec::with_capacity(2);
 let d1 = eu - cu; let disc1 = cyl_r * cyl_r - d1 * d1;
 if disc1 >= 0.0 { let off = disc1.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_max); add_pt(cv+off, -ev, ev, &mut ints_u_max); }
 ints_u_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
 ints_u_max.dedup_by(|a,b| (*a - *b).abs() < tol);

 // v_min edge (v=-ev): solve (u-cu)^2 + (-ev-cv)^2 = r^2 =u
 let mut ints_v_min: Vec<f64> = Vec::with_capacity(2);
 let d2 = -ev - cv; let disc2 = cyl_r * cyl_r - d2 * d2;
 if disc2 >= 0.0 { let off = disc2.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_min); add_pt(cu+off, -eu, eu, &mut ints_v_min); }
 ints_v_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
 ints_v_min.dedup_by(|a,b| (*a - *b).abs() < tol);

 // v_max edge (v=ev): solve (u-cu)^2 + (ev-cv)^2 = r^2 =u
 let mut ints_v_max: Vec<f64> = Vec::with_capacity(2);
 let d3 = ev - cv; let disc3 = cyl_r * cyl_r - d3 * d3;
 if disc3 >= 0.0 { let off = disc3.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_max); add_pt(cu+off, -eu, eu, &mut ints_v_max); }
 ints_v_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
 ints_v_max.dedup_by(|a,b| (*a - *b).abs() < tol);

 // Side walls =always use build_trimmed_edge_pieces for proper z-splitting
 // at z_lo/z_hi, matching the cap polygon edge vertices.

 // Helper for building trimmed edge pieces from intersection list.
 // Skips when the intersection covers the full span (both points at
 // corners), since `build_trimmed_edge_pieces` treats p_lo=p_min/p_hi=p_max
 // as "no intersection" and incorrectly keeps the full face.
 let trimmed_or_skip = |ints: &[f64], lo: f64, hi: f64,
 cn: &dyn Fn(f64, f64) -> DVec3, nrm: DVec3,
 pieces: &mut Vec<BRep>| {
 if ints.len() >= 2 && (ints.last().unwrap() - ints.first().unwrap()).abs() >= tol {
 let p_lo = ints[0];
 let p_hi = ints[ints.len() - 1];
 if p_lo > lo + tol || p_hi < hi - tol {
 build_trimmed_edge_pieces(lo, hi, p_lo, p_hi, z_lo, z_hi, ew, cn, nrm, pieces);
 }
 // else: full span =cylinder removes this entire face. Nothing to keep.
 } else if ints.len() == 1 {
 // Single intersection point =check which interval to keep.
 let p = ints[0];
 // Determine which face this is by checking nrm vs u_ax/v_ax.
 let is_u_face = nrm.dot(u_ax).abs() > 0.5;
 let inside: Box<dyn Fn(f64) -> bool> = if is_u_face {
 Box::new(|coord: f64| { (eu - cu).powi(2) + (coord - cv).powi(2) < cyl_r.powi(2) + tol })
 } else {
 Box::new(|coord: f64| { (coord - cu).powi(2) + (ev - cv).powi(2) < cyl_r.powi(2) + tol })
 };
 let mid_lo = (lo + p) * 0.5;
 let mid_hi = (p + hi) * 0.5;
 let ins_lo = inside(mid_lo);
 let ins_hi = inside(mid_hi);
 if !ins_lo && !ins_hi {
 // Both outside =tangent. Keep full face.
 build_trimmed_edge_pieces(lo, hi, lo, hi, z_lo, z_hi, ew, cn, nrm, pieces);
 } else if !ins_lo {
 if p > lo + tol { build_trimmed_edge_pieces(lo, hi, lo, p, z_lo, z_hi, ew, cn, nrm, pieces); }
 } else if !ins_hi {
 if p < hi - tol { build_trimmed_edge_pieces(lo, hi, p, hi, z_lo, z_hi, ew, cn, nrm, pieces); }
 } else {
 // Both inside =shouldn't happen. Fall through to full face.
 build_trimmed_edge_pieces(lo, hi, lo, hi, z_lo, z_hi, ew, cn, nrm, pieces);
 }
 } else {
 build_trimmed_edge_pieces(lo, hi, lo, hi, z_lo, z_hi, ew, cn, nrm, pieces);
 }
 };

 // u_max face (normal = +u_ax, param v)
 {
 let cn = |p: f64, z: f64| -> DVec3 { c + eu*u_ax + p*v_ax + z*DVec3::Z };
 trimmed_or_skip(&ints_u_max, -ev, ev, &cn, u_ax, &mut pieces);
 }
 // u_min face (normal = -u_ax, param v)
 {
 let cn = |p: f64, z: f64| -> DVec3 { c - eu*u_ax + p*v_ax + z*DVec3::Z };
 trimmed_or_skip(&ints_u_min, -ev, ev, &cn, -u_ax, &mut pieces);
 }
 // v_max face (normal = +v_ax, param u)
 {
 let cn = |p: f64, z: f64| -> DVec3 { c + p*u_ax + ev*v_ax + z*DVec3::Z };
 trimmed_or_skip(&ints_v_max, -eu, eu, &cn, v_ax, &mut pieces);
 }
 // v_min face (normal = -v_ax, param u)
 {
 let cn = |p: f64, z: f64| -> DVec3 { c + p*u_ax - ev*v_ax + z*DVec3::Z };
 trimmed_or_skip(&ints_v_min, -eu, eu, &cn, -v_ax, &mut pieces);
 }

 if !cap_bottom.is_empty() {
 if let Some(f) = planar_face_from_polygon(&cap_bottom, DVec3::Z) {
 pieces.push(f);
 }
 }
 if !cap_top.is_empty() {
 if let Some(f) = planar_face_from_polygon(&cap_top, -DVec3::Z) { pieces.push(f); }
 }

 // Cylindrical wall
 build_cylindrical_wall_from_segs(&merged, cu, cv, cyl_r, z_lo, z_hi, eu, ev, u_ax, v_ax, c, &corner_z, &mut pieces);

 if pieces.is_empty() { return None; }

 let sewn = sew_shells(&pieces, tol.max(TOLERANCE_ABS * 100.0));
 Some(sewn.brep)
}

/// Build box =cylinder difference for full XY containment (cylinder inside or
/// tangent to the box rect). Uses inner-wire annular cap faces and a full 360 ?
/// cylindrical wall, with side walls split at the cylinder Z range.
fn build_box_cylinder_full_containment(
 c: DVec3, u_ax: DVec3, v_ax: DVec3, eu: f64, ev: f64, ew: f64,
 cu: f64, cv: f64, cyl_r: f64, cyl_z_lo: f64, cyl_z_hi: f64,
) -> Option<BRep> {
 let tol = TOLERANCE_LEN_MIN;
 let z_lo = (cyl_z_lo - c.z).max(-ew);
 let z_hi = (cyl_z_hi - c.z).min(ew);
 if z_hi <= z_lo + tol { return None; }

 let mut pieces: Vec<BRep> = Vec::new();
 let p = |u: f64, v: f64, z: f64| -> DVec3 { c + u*u_ax + v*v_ax + z*DVec3::Z };

 // Discretize circle (CCW in XY) for inner wires and cylindrical wall.
 let n_circ = 64usize;
 let circle_ccw = |z: f64| -> Vec<DVec3> {
 (0..n_circ).map(|i| {
 let th = 2.0 * std::f64::consts::PI * (i as f64) / (n_circ as f64);
 p(cu + cyl_r * th.cos(), cv + cyl_r * th.sin(), z)
 }).collect()
 };
 let ccw_lo = circle_ccw(z_lo);
 let ccw_hi = circle_ccw(z_hi);

 // = =  1. Side walls split at z_lo/z_hi = = 
 let mut strip = |u: f64, v0: f64, v1: f64, z0: f64, z1: f64, nrm: DVec3| {
 if z1 <= z0 + tol { return; }
 if let Some(f) = rect_face_4([p(u, v0, z0), p(u, v1, z0), p(u, v1, z1), p(u, v0, z1)], nrm) {
 pieces.push(f);
 }
 };
 let mut split_wall = |u: f64, v_min: f64, v_max: f64, nrm: DVec3| {
 if z_lo > -ew + tol { strip(u, v_min, v_max, -ew, z_lo, nrm); }
 strip(u, v_min, v_max, z_lo, z_hi, nrm);
 if z_hi < ew - tol { strip(u, v_min, v_max, z_hi, ew, nrm); }
 };
 split_wall(eu, -ev, ev, u_ax);
 split_wall(-eu, -ev, ev, -u_ax);
 split_wall(ev, -eu, eu, v_ax);
 split_wall(-ev, -eu, eu, -v_ax);

 // = =  2. Annular cap faces = ㈤ ?
 // Bottom region
 if z_lo > -ew + tol {
 // Full bottom face at z=-ew (cylinder doesn't reach bottom).
 let bot = [p(-eu, -ev, -ew), p(-eu, ev, -ew), p(eu, ev, -ew), p(eu, -ev, -ew)];
 if let Some(f) = rect_face_4(bot, -DVec3::Z) { pieces.push(f); }
 // Interior annular cap at z_lo, normal +Z.
 // Outer CCW in XY; inner CW in XY (= reversed CCW).
 let outer = [p(-eu, -ev, z_lo), p(eu, -ev, z_lo), p(eu, ev, z_lo), p(-eu, ev, z_lo)];
 if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_lo.iter().rev().copied().collect::<Vec<_>>(), DVec3::Z) { pieces.push(f); }
 } else {
 // Bottom face IS the annular cap (cylinder goes through bottom).
 // Outer CCW in -Z view (= CW in XY); inner CW in -Z view (= CCW in XY).
 let outer = [p(-eu, -ev, -ew), p(-eu, ev, -ew), p(eu, ev, -ew), p(eu, -ev, -ew)];
 if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_lo, -DVec3::Z) { pieces.push(f); }
 }

 // Top region
 if z_hi < ew - tol {
 // Full top face at z=ew (cylinder doesn't reach top).
 let top = [p(-eu, -ev, ew), p(eu, -ev, ew), p(eu, ev, ew), p(-eu, ev, ew)];
 if let Some(f) = rect_face_4(top, DVec3::Z) { pieces.push(f); }
 // Interior annular cap at z_hi, normal -Z.
 // Outer CCW in -Z view (= CW in XY); inner CW in -Z view (= CCW in XY).
 let outer = [p(-eu, -ev, z_hi), p(-eu, ev, z_hi), p(eu, ev, z_hi), p(eu, -ev, z_hi)];
 if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_hi, -DVec3::Z) { pieces.push(f); }
 } else {
 // Top face IS the annular cap (cylinder goes through top).
 // Outer CCW in +Z view (= CCW in XY); inner CW in +Z view (= CW in XY = reversed CCW).
 let outer = [p(-eu, -ev, ew), p(eu, -ev, ew), p(eu, ev, ew), p(-eu, ev, ew)];
 if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_hi.iter().rev().copied().collect::<Vec<_>>(), DVec3::Z) { pieces.push(f); }
 }

 // = =  3. Cylindrical wall (full 360 ? = = 
 for i in 0..n_circ {
 let b0 = ccw_lo[i];
 let b1 = ccw_lo[(i + 1) % n_circ];
 let t0 = ccw_hi[i];
 let t1 = ccw_hi[(i + 1) % n_circ];
 let mid = (b0 + b1 + t1 + t0) / 4.0;
 let radial = mid - (c + cu*u_ax + cv*v_ax + mid.z*DVec3::Z);
 let inward = (-radial).normalize_or_zero();
 let inward = if inward.length_squared() < 0.5 { DVec3::X } else { inward };
 let n_vec = (t0 - b0).cross(b1 - b0).normalize();
 let n_final = if n_vec.dot(inward) > 0.0 { n_vec } else { -n_vec };
 if let Some(f) = rect_face_4([b0, b1, t1, t0], n_final) { pieces.push(f); }
 }

 // = =  4. Sew = = 
 if pieces.is_empty() { return None; }
 let sewn = sew_shells(&pieces, tol.max(TOLERANCE_ABS * 100.0));
 Some(sewn.brep)
}

/// Compute the closed boundary of `rect - circle` at one Z-level as a sequence of 2D points.
///
/// The boundary is traced clockwise. For the typical case (single closed curve), returns
/// a polygon with `n` sample points. For special cases:
/// - rect fully inside circle =returns empty (no boundary, void removes entire cross-section)
/// - circle fully inside rect =returns full rect perimeter (outer boundary only; the inner
/// circle hole is handled as a separate inner loop by the caller)
/// - circle outside rect (no overlap) =returns full rect perimeter (no void)
fn rect_minus_circle_boundary(
 bmin: DVec2, bmax: DVec2,
 cx: f64, cy: f64, r: f64,
 n: usize,
) -> Vec<DVec2> {
 let tol = TOLERANCE_LEN_MIN;
 if n < 4 { return vec![]; }

 let edges = [
 (DVec2::new(bmin.x, bmin.y), DVec2::new(bmax.x, bmin.y)), // bottom
 (DVec2::new(bmax.x, bmin.y), DVec2::new(bmax.x, bmax.y)), // right
 (DVec2::new(bmax.x, bmax.y), DVec2::new(bmin.x, bmax.y)), // top
 (DVec2::new(bmin.x, bmax.y), DVec2::new(bmin.x, bmin.y)), // left
 ];

 // ---- Step 1: Find circle-rect intersection t-values on each edge ----
 struct Intersection { t: f64, edge: usize, pt: DVec2 }
 let mut ints: Vec<Intersection> = Vec::new();

 for (ei, (p0, p1)) in edges.iter().enumerate() {
 let d = *p1 - *p0;
 let a0 = *p0 - DVec2::new(cx, cy);
 // Quadratic: (d )*t ?+ 2*(a0 )*t + (a0 0 - r ? = 0
 let A = d.dot(d);
 if A < 1e-30 { continue; }
 let B = 2.0 * a0.dot(d);
 let C = a0.dot(a0) - r * r;
 let disc = B * B - 4.0 * A * C;
 if disc < 0.0 { continue; }
 let sqrt_disc = disc.sqrt();
 for t in [(-B - sqrt_disc) / (2.0 * A), (-B + sqrt_disc) / (2.0 * A)] {
 if t >= -tol && t <= 1.0 + tol {
 let tc = t.clamp(0.0, 1.0);
 let pt = *p0 + d * tc;
 ints.push(Intersection { t: tc, edge: ei, pt });
 }
 }
 }

 // Deduplicate near-identical intersections (same edge and t)
 ints.sort_by(|a, b| a.edge.cmp(&b.edge).then(a.t.partial_cmp(&b.t).unwrap()));
 ints.dedup_by(|a, b| a.edge == b.edge && (a.t - b.t).abs() < tol);

 // Deduplicate by spatial proximity =a corner may appear as an intersection
 // on two adjacent edges when the circle passes exactly through the corner.
 // Keeping both creates zero-length perimeter segments that produce full-wrap
 // rect loops and self-intersecting polygons.
 ints.sort_by(|a, b| a.pt.x.partial_cmp(&b.pt.x).unwrap()
 .then(a.pt.y.partial_cmp(&b.pt.y).unwrap()));
 ints.dedup_by(|a, b| (a.pt - b.pt).length_squared() < tol * tol);

 // ---- Step 2: Handle no-intersection cases ----
 if ints.is_empty() {
 // Check if the rect center is inside the circle
 let rect_center = (bmin + bmax) * 0.5;
 let inside = (rect_center - DVec2::new(cx, cy)).length_squared() <= r * r + tol;
 if inside {
 // Rect is entirely inside the circle =empty cross-section
 return vec![];
 }
 // No overlap =full rect perimeter
 let mut result = Vec::with_capacity(n);
 for i in 0..n {
 let t = i as f64 / n as f64;
 let total_perim = 2.0 * ((bmax.x - bmin.x) + (bmax.y - bmin.y));
 let t_abs = t * total_perim;
 result.push(rect_perimeter_point(bmin, bmax, t_abs));
 }
 return result;
 }

 // ---- Step 3: Sort intersections along clockwise perimeter ----
 // Clockwise perimeter parameterization: edge 0: t= 0,1), edge 1: t= 1,2), edge 2: t= 2,3), edge 3: t= 3,4)
 let perim_pos = |ei: usize, t: f64| -> f64 { ei as f64 + t };
 ints.sort_by(|a, b| perim_pos(a.edge, a.t).partial_cmp(&perim_pos(b.edge, b.t)).unwrap());

 let total_perim = 2.0 * ((bmax.x - bmin.x) + (bmax.y - bmin.y));
 let n_per_edge = n / 4;
 let mut result = Vec::new();
 result.reserve(n);
 let tau = std::f64::consts::TAU;

 // Test if a 2D point is inside the rectangle
 let point_in_rect = |p: DVec2| -> bool {
 p.x >= bmin.x - tol && p.x <= bmax.x + tol
 && p.y >= bmin.y - tol && p.y <= bmax.y + tol
 };

 // ---- Step 4: Trace boundary ----
 // Walk clockwise along the rect perimeter. Between consecutive intersections,
 // the rect segment is either outside the circle (=keep rect points) or
 // inside the circle (=replace with circle arc).
 let m = ints.len();
 for i in 0..m {
 let j = (i + 1) % m;
 let ei = &ints[i];
 let ej = &ints[j];

 // Compute midpoint of the rect segment between these intersections
 let pi_pos = perim_pos(ei.edge, ei.t);
 let pj_pos = perim_pos(ej.edge, ej.t);
 let pm = if pj_pos > pi_pos {
 (pi_pos + pj_pos) * 0.5
 } else {
 // Wraps around (across the 0/4 boundary)
 let wrapped = (pi_pos + pj_pos + 4.0) * 0.5;
 if wrapped >= 4.0 { wrapped - 4.0 } else { wrapped }
 };
 let mid_pt = rect_perimeter_point(bmin, bmax, pm * total_perim / 4.0);
 let mid_inside = (mid_pt - DVec2::new(cx, cy)).length_squared() <= r * r + tol;

 if mid_inside {
 // Rect segment is inside the circle =follow circle arc from ei to ej
 // The arc must stay inside the rect. Test both possible arcs and pick
 // the one whose midpoint is inside the rect.
 let v1 = ei.pt - DVec2::new(cx, cy);
 let v2 = ej.pt - DVec2::new(cx, cy);
 let a1 = f64::atan2(v1.y, v1.x);
 let a2 = f64::atan2(v2.y, v2.x);

 // Positive (CCW) delta from a1 to a2, in [0,  ?
 let da_ccw = (a2 - a1).rem_euclid(tau);

 // Midpoint of CCW arc
 let mid_ccw = a1 + da_ccw * 0.5;
 let mid_ccw_pt = DVec2::new(cx + r * mid_ccw.cos(), cy + r * mid_ccw.sin());

 // Midpoint of CW arc (negative sweep)
 let mid_cw = a1 + (da_ccw - tau) * 0.5;
 let _mid_cw_pt = DVec2::new(cx + r * mid_cw.cos(), cy + r * mid_cw.sin());

 // Pick the arc whose midpoint is inside the rect
 let sweep = if point_in_rect(mid_ccw_pt) {
 da_ccw
 } else {
 da_ccw - tau // negative =clockwise
 };

 // Number of arc sample points (at least 4, scale with arc size)
 let arc_pts = (n_per_edge as f64 * sweep.abs() / std::f64::consts::PI).ceil().max(2.0) as usize;
 let arc_pts = arc_pts.min(64);

 // Add points along the arc (include start ei, exclude end ej)
 if i == 0 { result.push(ei.pt); }
 for k in 1..arc_pts {
 let frac = k as f64 / arc_pts as f64;
 let ang = a1 + sweep * frac;
 let (s, c) = ang.sin_cos();
 result.push(DVec2::new(cx + r * c, cy + r * s));
 }
 } else {
 // Rect segment is outside the circle =follow rect perimeter from ei to ej
 // Walk along the clockwise perimeter from position pi to pj
 let p_start = pi_pos * (total_perim / 4.0);
 let p_end = pj_pos * (total_perim / 4.0);
 let (t_start, t_end) = if p_end > p_start {
 (p_start, p_end)
 } else {
 (p_start, p_end + total_perim)
 };
 let seg_len = t_end - t_start;
 let n_pts = (n_per_edge as f64 * seg_len / total_perim).ceil().max(2.0) as usize;
 let n_pts = n_pts.min(64);

 if i == 0 { result.push(ei.pt); }
 for k in 1..n_pts {
 let frac = k as f64 / n_pts as f64;
 let t_abs = t_start + frac * seg_len;
 result.push(rect_perimeter_point(bmin, bmax, t_abs));
 }
 }
 }

 if result.len() < 3 { return vec![]; }
 result
}