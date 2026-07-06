
/// Add circle arc vertices (CCW) from `a1` to `a2` with CCW angular delta `da_ccw`.
/// Used by `build_circle_union_rect_polygon`.
fn add_circle_arc_pts(
 pts: &mut Vec<DVec2>, cu: f64, cv: f64, r: f64,
 a1: f64, da_ccw: f64, n_arc: usize, tol: f64,
) {
 for k in 0..=n_arc {
 let frac = k as f64 / n_arc as f64;
 let ang = a1 + da_ccw * frac;
 let (s, c) = ang.sin_cos();
 let p = DVec2::new(cu + r * c, cv + r * s);
 if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
 pts.push(p);
 }
 }
}

/// Add rect perimeter vertices along CCW direction from `t_start` to `t_end`.
/// Adds intermediate corner vertices (at integer t values) and the endpoint.
fn add_rect_perimeter_pts(
 pts: &mut Vec<DVec2>, t_start: f64, t_end: f64, eu: f64, ev: f64, tol: f64,
) {
 let ts = t_start.rem_euclid(8.0);
 let te = t_end.rem_euclid(8.0);

 if (ts - te).abs() < tol {
 let (u, v) = box_perimeter_uv(ts, eu, ev);
 let p = DVec2::new(u, v);
 if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
 pts.push(p);
 }
 return;
 }

 // Walk CCW: if ts < te walk forward; if ts > te wrap through 8.0
 let range: Vec<usize> = if ts < te {
 ((ts.ceil() as usize)..=(te.floor() as usize)).collect()
 } else {
 let mut r: Vec<usize> = ((ts.ceil() as usize)..8).collect();
 r.extend(0..=(te.floor() as usize));
 r
 };
 // Walk CW along rect perimeter.
 // For non-wrapping (ts < te): keep integer points between ts and te.
 // For wrapping (ts > te): the range already separates into two portions
 // ([ceil(ts), 8) and [0, floor(te)]), so accept points in either portion
 // that are strictly between ts and te (use `||` since the ranges are disjoint).
 for &tc in &range {
 let tt = tc as f64;
 if if ts < te {
 tt > ts + tol && tt < te - tol
 } else {
 tt > ts + tol || tt < te - tol
 } {
 let (u, v) = box_perimeter_uv(tt, eu, ev);
 let p = DVec2::new(u, v);
 if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
 pts.push(p);
 }
 }
 }

 // Endpoint
 let (u, v) = box_perimeter_uv(te, eu, ev);
 let p = DVec2::new(u, v);
 if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
 pts.push(p);
 }
}

/// Compute the 2D boundary polygon for `circle =rect` (CCW closed polygon).
///
/// Returns empty `Vec` if the shapes are disjoint.
fn build_circle_union_rect_polygon(
 cu: f64, cv: f64, r: f64,
 eu: f64, ev: f64,
) -> Vec<DVec2> {
 let tol = TOLERANCE_LEN_MIN;
 let tau = std::f64::consts::TAU;
 let n_arc = 128usize;
 let bmin = DVec2::new(-eu, -ev);
 let bmax = DVec2::new(eu, ev);

 let point_in_rect = |p: DVec2| -> bool {
 p.x >= bmin.x - tol && p.x <= bmax.x + tol
 && p.y >= bmin.y - tol && p.y <= bmax.y + tol
 };

 let raw_ints = circle_rect_intersections_uv(cu, cv, r, eu, ev);

 if raw_ints.is_empty() {
 // No intersections =check containment
 let corners = [
 DVec2::new(bmin.x, bmin.y),
 DVec2::new(bmax.x, bmin.y),
 DVec2::new(bmax.x, bmax.y),
 DVec2::new(bmin.x, bmax.y),
 ];
 let all_inside_circle = corners.iter().all(|c| {
 (c.x - cu).powi(2) + (c.y - cv).powi(2) <= (r + tol).powi(2)
 });
 if all_inside_circle {
 // Circle fully contains rect =return full circle
 let mut poly = Vec::with_capacity(n_arc + 1);
 for k in 0..=n_arc {
 let ang = tau * k as f64 / n_arc as f64;
 let (s, c) = ang.sin_cos();
 poly.push(DVec2::new(cu + r * c, cv + r * s));
 }
 return poly;
 }
 let center_in_rect = point_in_rect(DVec2::new(cu, cv));
 let cyl_in_rect = cu - r >= bmin.x - tol && cu + r <= bmax.x + tol
 && cv - r >= bmin.y - tol && cv + r <= bmax.y + tol;
 if center_in_rect && cyl_in_rect {
 // Rect fully contains circle =return rect perimeter
 return vec![
 DVec2::new(bmin.x, bmin.y),
 DVec2::new(bmax.x, bmin.y),
 DVec2::new(bmax.x, bmax.y),
 DVec2::new(bmin.x, bmax.y),
 ];
 }
 return Vec::new();
 }

 // Sort intersections by CCW angle around circle center
 let mut sorted: Vec<&UVEdgePt> = raw_ints.iter().collect();
 sorted.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());

 let m = sorted.len();
 if m < 2 { return Vec::new(); }

 let mut pts: Vec<DVec2> = Vec::new();

 for i in 0..m {
 let j = (i + 1) % m;
 let a1 = sorted[i].theta;
 let a2 = sorted[j].theta;
 let da_ccw = (a2 - a1).rem_euclid(tau);
 if da_ccw < 1e-12 { continue; } // zero-length arc, skip

 let mid_ang = a1 + da_ccw * 0.5;
 let mid_pt = DVec2::new(cu + r * mid_ang.cos(), cv + r * mid_ang.sin());

 if point_in_rect(mid_pt) {
 // Arc midpoint INSIDE rect =rect perimeter is the boundary
 add_rect_perimeter_pts(&mut pts, sorted[i].t, sorted[j].t, eu, ev, tol);
 } else {
 // Arc midpoint OUTSIDE rect =circle arc is the boundary
 add_circle_arc_pts(&mut pts, cu, cv, r, a1, da_ccw, n_arc, tol);
 }
 }

 // Close polygon: remove last point if it duplicates the first
 if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length_squared() < tol * tol {
 pts.pop();
 }

 pts
}

/// Compute the 2D boundary polygon for `circle =rect` (CCW closed polygon).
///
/// Returns empty `Vec` if the shapes are disjoint.
fn build_circle_intersect_rect_polygon(
 cu: f64, cv: f64, r: f64,
 eu: f64, ev: f64,
) -> Vec<DVec2> {
 let tol = TOLERANCE_LEN_MIN;
 let tau = std::f64::consts::TAU;
 let n_arc = 128usize;
 let bmin = DVec2::new(-eu, -ev);
 let bmax = DVec2::new(eu, ev);

 let point_in_rect = |p: DVec2| -> bool {
 p.x >= bmin.x - tol && p.x <= bmax.x + tol
 && p.y >= bmin.y - tol && p.y <= bmax.y + tol
 };

 let raw_ints = circle_rect_intersections_uv(cu, cv, r, eu, ev);

 if raw_ints.is_empty() {
 // No intersections =check containment
 let corners = [
 DVec2::new(bmin.x, bmin.y),
 DVec2::new(bmax.x, bmin.y),
 DVec2::new(bmax.x, bmax.y),
 DVec2::new(bmin.x, bmax.y),
 ];
 let all_inside_circle = corners.iter().all(|c| {
 (c.x - cu).powi(2) + (c.y - cv).powi(2) <= (r + tol).powi(2)
 });
 if all_inside_circle {
 // Circle fully contains rect =intersection = rect perimeter
 return vec![
 DVec2::new(bmin.x, bmin.y),
 DVec2::new(bmax.x, bmin.y),
 DVec2::new(bmax.x, bmax.y),
 DVec2::new(bmin.x, bmax.y),
 ];
 }
 let center_in_rect = point_in_rect(DVec2::new(cu, cv));
 let cyl_in_rect = cu - r >= bmin.x - tol && cu + r <= bmax.x + tol
 && cv - r >= bmin.y - tol && cv + r <= bmax.y + tol;
 if center_in_rect && cyl_in_rect {
 // Rect fully contains circle =intersection = circle polygon
 let mut poly = Vec::with_capacity(n_arc + 1);
 for k in 0..=n_arc {
 let ang = tau * k as f64 / n_arc as f64;
 let (s, c) = ang.sin_cos();
 poly.push(DVec2::new(cu + r * c, cv + r * s));
 }
 return poly;
 }
 return Vec::new();
 }

 // Sort intersections by CCW angle around circle center
 let mut sorted: Vec<&UVEdgePt> = raw_ints.iter().collect();
 sorted.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());

 let m = sorted.len();
 if m < 2 { return Vec::new(); }

 let mut pts: Vec<DVec2> = Vec::new();

 for i in 0..m {
 let j = (i + 1) % m;
 let a1 = sorted[i].theta;
 let a2 = sorted[j].theta;
 let da_ccw = (a2 - a1).rem_euclid(tau);
 if da_ccw < 1e-12 { continue; } // zero-length arc, skip

 let mid_ang = a1 + da_ccw * 0.5;
 let mid_pt = DVec2::new(cu + r * mid_ang.cos(), cv + r * mid_ang.sin());

 if point_in_rect(mid_pt) {
 // Arc midpoint INSIDE rect =circle arc is part of intersection boundary
 add_circle_arc_pts(&mut pts, cu, cv, r, a1, da_ccw, n_arc, tol);
 } else {
 // Arc midpoint OUTSIDE rect =rect perimeter is part of intersection boundary
 add_rect_perimeter_pts(&mut pts, sorted[i].t, sorted[j].t, eu, ev, tol);
 }
 }

 // Close polygon: remove last point if it duplicates the first
 if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length_squared() < tol * tol {
 pts.pop();
 }

 pts
}

/// Compute the 2D boundary polygon for `rect - circle` (CCW closed polygon).
///
/// Returns empty `Vec` if the result is empty (rect entirely inside circle).
/// This is the inverse of `build_circle_intersect_rect_polygon`: where the
/// intersection polygon uses circle arcs (arc midpoint inside rect), this
/// polygon uses rect perimeter, and vice versa.
fn build_rect_minus_circle_polygon(
 cu: f64, cv: f64, r: f64,
 eu: f64, ev: f64,
) -> Vec<DVec2> {
 let tol = TOLERANCE_LEN_MIN;
 let tau = std::f64::consts::TAU;
 let n_arc = 128usize;
 let bmin = DVec2::new(-eu, -ev);
 let bmax = DVec2::new(eu, ev);

 let point_in_rect = |p: DVec2| -> bool {
 p.x >= bmin.x - tol && p.x <= bmax.x + tol
 && p.y >= bmin.y - tol && p.y <= bmax.y + tol
 };

 let raw_ints = circle_rect_intersections_uv(cu, cv, r, eu, ev);

 if raw_ints.is_empty() {
 // No intersections =check containment
 let corners = [
 DVec2::new(bmin.x, bmin.y),
 DVec2::new(bmax.x, bmin.y),
 DVec2::new(bmax.x, bmax.y),
 DVec2::new(bmin.x, bmax.y),
 ];
 let all_inside_circle = corners.iter().all(|c| {
 (c.x - cu).powi(2) + (c.y - cv).powi(2) <= (r + tol).powi(2)
 });
 if all_inside_circle {
 // Circle fully contains rect =rect - circle = empty
 return Vec::new();
 }
 let center_in_rect = point_in_rect(DVec2::new(cu, cv));
 let cyl_in_rect = cu - r >= bmin.x - tol && cu + r <= bmax.x + tol
 && cv - r >= bmin.y - tol && cv + r <= bmax.y + tol;
 if center_in_rect && cyl_in_rect {
 // Rect fully contains circle =rect - circle = rect perimeter
 return vec![
 DVec2::new(bmin.x, bmin.y),
 DVec2::new(bmax.x, bmin.y),
 DVec2::new(bmax.x, bmax.y),
 DVec2::new(bmin.x, bmax.y),
 ];
 }
 // Disjoint =rect - circle = full rect
 return vec![
 DVec2::new(bmin.x, bmin.y),
 DVec2::new(bmax.x, bmin.y),
 DVec2::new(bmax.x, bmax.y),
 DVec2::new(bmin.x, bmax.y),
 ];
 }

 // Sort intersections by CCW angle around circle center
 let mut sorted: Vec<&UVEdgePt> = raw_ints.iter().collect();
 sorted.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());

 let m = sorted.len();
 if m < 2 { return Vec::new(); }

 let mut pts: Vec<DVec2> = Vec::new();

 for i in 0..m {
 let j = (i + 1) % m;
 let a1 = sorted[i].theta;
 let a2 = sorted[j].theta;
 let da_ccw = (a2 - a1).rem_euclid(tau);
 if da_ccw < 1e-12 { continue; }

 let mid_ang = a1 + da_ccw * 0.5;
 let mid_pt = DVec2::new(cu + r * mid_ang.cos(), cv + r * mid_ang.sin());

 // INVERTED from build_circle_intersect_rect_polygon:
 if point_in_rect(mid_pt) {
 // Arc midpoint INSIDE rect =rect perimeter is boundary
 add_rect_perimeter_pts(&mut pts, sorted[i].t, sorted[j].t, eu, ev, tol);
 } else {
 // Arc midpoint OUTSIDE rect =circle arc cuts away the rect.
 // Use the SHORT arc between a1 and a2 (whichever of CCW and CW
 // is shorter) to avoid overlapping complement arcs when the
 // circle is larger than the rect (all midpoints outside rect).
 if da_ccw < tau - da_ccw {
 add_circle_arc_pts(&mut pts, cu, cv, r, a1, da_ccw, n_arc, tol);
 } else {
 add_circle_arc_pts(&mut pts, cu, cv, r, a1, -(tau - da_ccw), n_arc, tol);
 }
 }
 }

 // Close polygon
 if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length_squared() < tol * tol {
 pts.pop();
 }

 pts
}

/// Build a closed circle polygon with `n_arc` segments.
fn build_circle_polygon(cu: f64, cv: f64, r: f64) -> Vec<DVec2> {
 let n_arc = 128usize;
 let tau = std::f64::consts::TAU;
 let mut poly = Vec::with_capacity(n_arc + 1);
 for k in 0..=n_arc {
 let ang = tau * k as f64 / n_arc as f64;
 let (s, c) = ang.sin_cos();
 poly.push(DVec2::new(cu + r * c, cv + r * s));
 }
 poly
}

/// Add a wall section from `z0` to `z1` using polygon `pts` (shared vertex pool).
fn add_wall_section(
 add_v: &mut impl FnMut(DVec3) -> usize,
 faces: &mut Vec<Face>,
 pts: &[DVec2],
 z0: f64, z1: f64, n_slices: usize,
 to_world: &impl Fn(f64, f64, f64) -> DVec3,
 empty_wire: &impl Fn() -> Wire,
) {
 let dz = (z1 - z0) / n_slices as f64;
 let n = pts.len();
 if n < 3 { return; }

 for i in 0..n_slices {
 let za = z0 + dz * i as f64;
 let zb = z0 + dz * (i + 1) as f64;
 let mut idx = Vec::with_capacity(2 * n);
 for p in pts { idx.push(add_v(to_world(p.x, p.y, za))); }
 for p in pts { idx.push(add_v(to_world(p.x, p.y, zb))); }
 let mut tris = Vec::with_capacity(n * 2);
 for j in 0..n {
 let k = (j + 1) % n;
 tris.push([idx[j], idx[k], idx[n + k]]);
 tris.push([idx[j], idx[n + k], idx[n + j]]);
 }
 faces.push(Face {
 outer_wire: empty_wire(), inner_wires: vec![],
 normal: DVec3::ZERO, triangles: tris,
 sample_point: None, mesh_dirty: false,
 surface_idx: None,
 });
 }
}

/// Add a triangulated cap face at Z level `z` using polygon `pts`.
fn add_cap_face(
 add_v: &mut impl FnMut(DVec3) -> usize,
 faces: &mut Vec<Face>,
 pts: &[DVec2],
 z: f64,
 normal: DVec3,
 to_world: &impl Fn(f64, f64, f64) -> DVec3,
 empty_wire: &impl Fn() -> Wire,
) {
 let poly: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z)).collect();
 let tris = crate::triangulate::triangulate_polygon(&poly, normal);
 if tris.is_empty() { return; }
 let mut remapped = Vec::with_capacity(tris.len());
 let local: Vec<usize> = poly.iter().map(|p| add_v(*p)).collect();
 for t in &tris { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
 faces.push(Face {
 outer_wire: empty_wire(), inner_wires: vec![],
 normal: DVec3::ZERO, triangles: remapped,
 sample_point: None, mesh_dirty: false,
 surface_idx: None,
 });
}

/// Add a ring face at the interface between the union-polygon wall and the circle-polygon wall.
///
/// At `z=box_z_hi`, the cross-section transitions from `circle =rect` to just `circle`.
/// The box top face outside the cylinder adds surface area.  Triangulates the union polygon
/// and keeps only triangles whose UV centroid lies outside the circle.
fn add_interface_face(
 add_v: &mut impl FnMut(DVec3) -> usize,
 faces: &mut Vec<Face>,
 pts: &[DVec2],
 z: f64,
 normal: DVec3,
 circle_center_uv: DVec2,
 circle_r: f64,
 to_world: &impl Fn(f64, f64, f64) -> DVec3,
 empty_wire: &impl Fn() -> Wire,
) {
 if pts.len() < 3 { return; }
 let poly3d: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z)).collect();
 let tris = crate::triangulate::triangulate_polygon(&poly3d, normal);
 if tris.is_empty() { return; }
 let r2 = circle_r * circle_r + 1e-12;
 let mut kept: Vec<[usize; 3]> = Vec::new();
 for t in &tris {
 let c = (pts[t[0]] + pts[t[1]] + pts[t[2]]) / 3.0;
 if (c - circle_center_uv).length_squared() > r2 {
 kept.push(*t);
 }
 }
 if kept.is_empty() { return; }
 let local: Vec<usize> = poly3d.iter().map(|p| add_v(*p)).collect();
 let mut remapped = Vec::with_capacity(kept.len());
 for t in &kept { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
 faces.push(Face {
 outer_wire: empty_wire(), inner_wires: vec![],
 normal: DVec3::ZERO, triangles: remapped,
 sample_point: None, mesh_dirty: false,
 surface_idx: None,
 });
}

/// Build a tessellated BRep for `cylinder =box` via Z-slice tessellation.
///
/// Builds three sections: below-box (circle), overlap (circle =rect), above-box (circle).
/// The full cylinder Z range is [cyl_z_lo, cyl_z_hi]; the box occupies [box_z_lo, box_z_hi].
fn build_cylinder_box_union_tessellated(
 bc: DVec3,
 u_ax: DVec3,
 v_ax: DVec3,
 cu: f64,
 cv: f64,
 r: f64,
 eu: f64,
 ev: f64,
 cyl_z_lo: f64,
 cyl_z_hi: f64,
 box_z_lo: f64,
 box_z_hi: f64,
) -> Option<BRep> {
 let tol = TOLERANCE_LEN_MIN;
 if cyl_z_hi <= cyl_z_lo + tol { return None; }
 if r < tol { return None; }

 let n_slices = 16usize;
 let n_slices_circ = 8usize; // fewer for plain cylinder sections
 let empty_wire = || Wire { edges: vec![] };

 let mut verts: Vec<Vertex> = Vec::new();
 let mut add_v = |p: DVec3| -> usize {
 for (i, v) in verts.iter().enumerate() {
 if (v.point - p).length() < 1e-12 { return i; }
 }
 let idx = verts.len();
 verts.push(Vertex { point: p });
 idx
 };

 let mut faces: Vec<Face> = Vec::new();

 let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
 bc + u_ax * u + v_ax * v + DVec3::new(0.0, 0.0, z)
 };

 // Pre-compute polygons
 let union_poly = build_circle_union_rect_polygon(cu, cv, r, eu, ev);
 let circle_poly = build_circle_polygon(cu, cv, r);

 if union_poly.len() < 3 || circle_poly.len() < 3 {
 return None;
 }

 // Section 1: below box (circle polygon)
 if box_z_lo > cyl_z_lo + tol {
 add_wall_section(
 &mut add_v, &mut faces,
 &circle_poly, cyl_z_lo, box_z_lo, n_slices_circ,
 &to_world, &empty_wire,
 );
 }

 // Section 2: overlap (union polygon)
 if box_z_hi > box_z_lo + tol {
 add_wall_section(
 &mut add_v, &mut faces,
 &union_poly, box_z_lo, box_z_hi, n_slices,
 &to_world, &empty_wire,
 );
 }

 // Section 3: above box (circle polygon)
 if cyl_z_hi > box_z_hi + tol {
 add_wall_section(
 &mut add_v, &mut faces,
 &circle_poly, box_z_hi, cyl_z_hi, n_slices_circ,
 &to_world, &empty_wire,
 );
 }

 // Interface face at box top: the ring between union_poly and circle_poly (box top outside cylinder)
 if cyl_z_hi > box_z_hi + tol && union_poly.len() >= 3 {
 add_interface_face(
 &mut add_v, &mut faces,
 &union_poly, box_z_hi, DVec3::Z,
 DVec2::new(cu, cv), r,
 &to_world, &empty_wire,
 );
 }

 // Bottom cap: use union polygon if box reaches bottom, circle otherwise
 if box_z_lo <= cyl_z_lo + tol {
 add_cap_face(
 &mut add_v, &mut faces,
 &union_poly, cyl_z_lo, -DVec3::Z,
 &to_world, &empty_wire,
 );
 } else {
 add_cap_face(
 &mut add_v, &mut faces,
 &circle_poly, cyl_z_lo, -DVec3::Z,
 &to_world, &empty_wire,
 );
 }

 // Top cap: use union polygon if box reaches top, circle otherwise
 if box_z_hi >= cyl_z_hi - tol {
 add_cap_face(
 &mut add_v, &mut faces,
 &union_poly, cyl_z_hi, DVec3::Z,
 &to_world, &empty_wire,
 );
 } else {
 add_cap_face(
 &mut add_v, &mut faces,
 &circle_poly, cyl_z_hi, DVec3::Z,
 &to_world, &empty_wire,
 );
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

// ── Cone-box union fast path ─────────────────────────────────────

/// Remap a closed polygon to `n` equally-spaced arc-length points,
/// starting from the point closest to `ref_pt`.
fn remap_polygon_arclength(poly: &[DVec2], n: usize, ref_pt: DVec2) -> Vec<DVec2> {
 let tol = TOLERANCE_LEN_MIN;
 if poly.len() < 3 || n < 3 { return poly.to_vec(); }

 // Rotate to start from the point closest to ref_pt
 let mut best_idx = 0;
 let mut best_dist = (poly[0] - ref_pt).length_squared();
 for (i, p) in poly.iter().enumerate() {
 let d = (*p - ref_pt).length_squared();
 if d < best_dist { best_dist = d; best_idx = i; }
 }
 let m = poly.len();
 let mut aligned = poly[best_idx..].to_vec();
 aligned.extend_from_slice(&poly[..best_idx]);

 // Compute cumulative arc length
 let mut arc_len = vec![0.0_f64; m + 1];
 for i in 1..=m {
 let j = i % m;
 let k = (i - 1) % m;
 arc_len[i] = arc_len[i - 1] + aligned[k].distance(aligned[j]);
 }
 let total = arc_len[m];
 if total <= tol { return poly.to_vec(); }

 // Resample to n equally-spaced points
 let mut result = Vec::with_capacity(n);
 let mut src_idx = 0;
 for i in 0..n {
 let target = total * i as f64 / n as f64;
 while src_idx < m && arc_len[src_idx + 1] < target {
 src_idx += 1;
 }
 let t0 = arc_len[src_idx];
 let t1 = arc_len[(src_idx + 1) % (m + 1)];
 if (t1 - t0).abs() < 1e-15 {
 result.push(aligned[src_idx % m]);
 } else {
 let frac = (target - t0) / (t1 - t0);
 let a = src_idx % m;
 let b = (src_idx + 1) % m;
 result.push(aligned[a].lerp(aligned[b], frac));
 }
 }
 result
}

/// Build a tessellated BRep for `cone =box` via Z-slice tessellation.
///
/// The cone is a Z-aligned conical frustum with center at `(cx, cy)` in XY,
/// extending from Z `cz_lo` to `cz_hi`, with bottom radius `cr_lo` and top
/// radius `cr_hi`. The box is axis-aligned `[bmin, bmax]`.
///
/// Builds three sections: below-box (circle), overlap (circle =rect), above-box (circle).
fn build_cone_box_union_tessellated(
 bmin: DVec3, bmax: DVec3,
 cx: f64, cy: f64,
 cz_lo: f64, cz_hi: f64,
 cr_lo: f64, cr_hi: f64,
) -> Option<BRep> {
 let tol = TOLERANCE_LEN_MIN;
 if cz_hi <= cz_lo + tol { return None; }
 if cr_lo < tol && cr_hi < tol { return None; }

 let box_z_lo = bmin.z;
 let box_z_hi = bmax.z;
 let box_center = (bmin + bmax) * 0.5;
 let eu = (bmax.x - bmin.x) * 0.5;
 let ev = (bmax.y - bmin.y) * 0.5;
 let cu = cx - box_center.x;
 let cv = cy - box_center.y;

 let n_slices = 64usize;
 let n_slices_circ = 16usize;
 let n_boundary = 256usize;
 let n_arc = 128usize;
 let tau = std::f64::consts::TAU;
 let empty_wire = || Wire { edges: vec![] };

 let mut verts: Vec<Vertex> = Vec::new();
 let mut add_v = |p: DVec3| -> usize {
 for (i, v) in verts.iter().enumerate() {
 if (v.point - p).length() < 1e-12 { return i; }
 }
 let idx = verts.len();
 verts.push(Vertex { point: p });
 idx
 };

 let mut faces: Vec<Face> = Vec::new();

 let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
 DVec3::new(box_center.x + u, box_center.y + v, z)
 };

 let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);

 // ---- Section 1: Below box (circle cross-section) ----
 if box_z_lo > cz_lo + tol {
 let z0 = cz_lo;
 let z1 = box_z_lo;
 let n = n_slices_circ;
 let dz = (z1 - z0) / n as f64;
 for i in 0..n {
 let za = z0 + dz * i as f64;
 let zb = z0 + dz * (i + 1) as f64;
 let ra = (cr_lo + dr_dz * (za - cz_lo)).max(0.0);
 let rb = (cr_lo + dr_dz * (zb - cz_lo)).max(0.0);
 if ra < tol && rb < tol { continue; }

 let nn = n_arc;
 let mut idx = Vec::with_capacity(2 * (nn + 1));
 for k in 0..=nn {
 let ang = tau * k as f64 / nn as f64;
 let (s, c) = ang.sin_cos();
 idx.push(add_v(to_world(cu + ra * c, cv + ra * s, za)));
 }
 for k in 0..=nn {
 let ang = tau * k as f64 / nn as f64;
 let (s, c) = ang.sin_cos();
 idx.push(add_v(to_world(cu + rb * c, cv + rb * s, zb)));
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

 // ---- Section 2: Overlap (circle =rect cross-section) ----
 if box_z_hi > box_z_lo + tol {
 let z0 = box_z_lo;
 let z1 = box_z_hi;
 let n = n_slices;
 let dz = (z1 - z0) / n as f64;
 let ref_pt = DVec2::new(-eu, -ev);

 // Pre-compute and remap all slice polygons
 let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n + 1);
 for i in 0..=n {
 let z = z0 + dz * i as f64;
 let r = (cr_lo + dr_dz * (z - cz_lo)).max(0.0);
 if r < tol {
 slices.push(vec![]);
 } else {
 let poly = build_circle_union_rect_polygon(cu, cv, r, eu, ev);
 if poly.len() >= 3 {
 slices.push(remap_polygon_arclength(&poly, n_boundary, ref_pt));
 } else {
 slices.push(vec![]);
 }
 }
 }

 // Build wall faces between adjacent remapped slices
 for i in 0..n {
 let bot = &slices[i];
 let top = &slices[i + 1];
 if bot.len() < 3 || top.len() < 3 { continue; }

 let n_pts = bot.len().min(top.len());
 let z_bot = z0 + dz * i as f64;
 let z_top = z0 + dz * (i + 1) as f64;

 let mut idx = Vec::with_capacity(2 * n_pts);
 for p in bot.iter() { idx.push(add_v(to_world(p.x, p.y, z_bot))); }
 for p in top.iter() { idx.push(add_v(to_world(p.x, p.y, z_top))); }

 let mut tris = Vec::with_capacity(n_pts * 2);
 for j in 0..n_pts {
 let k = (j + 1) % n_pts;
 tris.push([idx[j], idx[k], idx[n_pts + k]]);
 tris.push([idx[j], idx[n_pts + k], idx[n_pts + j]]);
 }
 faces.push(Face {
 outer_wire: empty_wire(), inner_wires: vec![],
 normal: DVec3::ZERO, triangles: tris,
 sample_point: None, mesh_dirty: false,
 surface_idx: None,
 });
 }
 }

 // ---- Section 3: Above box (circle cross-section) ----
 if cz_hi > box_z_hi + tol {
 let z0 = box_z_hi;
 let z1 = cz_hi;
 let n = n_slices_circ;
 let dz = (z1 - z0) / n as f64;
 for i in 0..n {
 let za = z0 + dz * i as f64;
 let zb = z0 + dz * (i + 1) as f64;
 let ra = (cr_lo + dr_dz * (za - cz_lo)).max(0.0);
 let rb = (cr_lo + dr_dz * (zb - cz_lo)).max(0.0);
 if ra < tol && rb < tol { continue; }

 let nn = n_arc;
 let mut idx = Vec::with_capacity(2 * (nn + 1));
 for k in 0..=nn {
 let ang = tau * k as f64 / nn as f64;
 let (s, c) = ang.sin_cos();
 idx.push(add_v(to_world(cu + ra * c, cv + ra * s, za)));
 }
 for k in 0..=nn {
 let ang = tau * k as f64 / nn as f64;
 let (s, c) = ang.sin_cos();
 idx.push(add_v(to_world(cu + rb * c, cv + rb * s, zb)));
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

 // ---- Interface face at box top (ring between union and circle) ----
 if cz_hi > box_z_hi + tol {
 let r_at = (cr_lo + dr_dz * (box_z_hi - cz_lo)).max(0.0);
 let union_poly = build_circle_union_rect_polygon(cu, cv, r_at, eu, ev);
 if union_poly.len() >= 3 {
 add_interface_face(
 &mut add_v, &mut faces,
 &union_poly, box_z_hi, DVec3::Z,
 DVec2::new(cu, cv), r_at,
 &to_world, &empty_wire,
 );
 }
 }

 // ---- Bottom cap ----
 let r_bottom = (cr_lo + dr_dz * (cz_lo - cz_lo)).max(0.0);
 if box_z_lo <= cz_lo + tol {
 let poly = build_circle_union_rect_polygon(cu, cv, r_bottom, eu, ev);
 if poly.len() >= 3 {
 add_cap_face(&mut add_v, &mut faces, &poly, cz_lo, -DVec3::Z, &to_world, &empty_wire);
 }
 } else if r_bottom > tol {
 let poly = build_circle_polygon(cu, cv, r_bottom);
 if poly.len() >= 3 {
 add_cap_face(&mut add_v, &mut faces, &poly, cz_lo, -DVec3::Z, &to_world, &empty_wire);
 }
 }

 // ---- Top cap ----
 let r_top = (cr_lo + dr_dz * (cz_hi - cz_lo)).max(0.0);
 if box_z_hi >= cz_hi - tol {
 let poly = build_circle_union_rect_polygon(cu, cv, r_top, eu, ev);
 if poly.len() >= 3 {
 add_cap_face(&mut add_v, &mut faces, &poly, cz_hi, DVec3::Z, &to_world, &empty_wire);
 }
 } else if r_top > tol {
 let poly = build_circle_polygon(cu, cv, r_top);
 if poly.len() >= 3 {
 add_cap_face(&mut add_v, &mut faces, &poly, cz_hi, DVec3::Z, &to_world, &empty_wire);
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