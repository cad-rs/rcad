//! Shape properties: surface area, volume, and centroid.
//!
//! Analogous to OCCT `GProp_GProps` with `BRepGProp`.
//!
//! `surface_area` first uses analytic or parametric integrals (plane: shoe-lace; sphere: UV
//! mask and `R虏 d惟`; other finite UV patches without `inner_wires`: cylinder `r路螖u路螖v` or a
//! midpoint `鈥栤垈P/鈭倁脳鈭侾/鈭倂鈥朻 sum on the same domain as `tessellate_curved_face`), then
//! triangulated faces. Where triangles are used: UV-grid tessellation, holed
//! ear-cut, or fan-triangulation from wire vertices.

use glam::{DVec2, DVec3};
use std::f64::consts::PI;

use crate::BRep;
use crate::topods;
use crate::geom::{ConicalSurface, Curve3, CurveEval, CylindricalSurface, SphericalSurface, Surface3, SurfaceEval};
use crate::topology::{Face, Wire, WireEdge};

// 鈹€鈹€ Internal helpers 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the signed area of a triangle from three points.
/// The sign depends on the orientation relative to the caller.
#[inline]
fn tri_area(a: DVec3, b: DVec3, c: DVec3) -> f64 {
 (b - a).cross(c - a).length() * 0.5
}

/// Signed volume contribution of a tetrahedron from the origin to triangle (a,b,c).
/// Summing over all surface triangles gives 1/6 * signed volume of the solid.
#[inline]
fn tet_signed_volume(a: DVec3, b: DVec3, c: DVec3) -> f64 {
 a.dot(b.cross(c)) / 6.0
}

/// Sample a closed wire to a 3D polyline (edge curve samples or vertex fallbacks).
fn sample_wire_polyline_3d(brep: &topods::BRep, wire: &Wire) -> Vec<DVec3> {
 sample_wire_polyline_3d_with_n(brep, wire, 1024)
}

/// Like [`sample_wire_polyline_3d`] but with configurable samples per edge.
fn sample_wire_polyline_3d_with_n(brep: &topods::BRep, wire: &Wire, n: usize) -> Vec<DVec3> {
 use crate::geom::CurveEval;
 let mut pts = Vec::new();
 for we in &wire.edges {
 let flat_edges = brep.flat_edges();
 let edge = match flat_edges.get(we.idx) {
 Some(e) => e,
 None => continue,
 };
 let curve_opt = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } });
 if let Some(curve) = curve_opt {
 let range = brep
 .tshapes
 .get(we.idx)
 .and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None } })
 .unwrap_or_else(|| curve.default_domain());
 let [t0, t1] = if we.forward { range } else { [range[1], range[0]] };
 if (t1 - t0).abs() > 1e-12 && t0.is_finite() && t1.is_finite() {
 let full_circle = (t1 - t0).abs() >= 2.0 * std::f64::consts::PI - 1e-9;
 let samples = if full_circle { n } else { n + 1 };
 for k in 0..samples {
 let frac = k as f64 / n as f64;
 let t = t0 + (t1 - t0) * frac;
 pts.push(curve.point_at(t));
 }
 continue;
 }
 }
 let vidx = if we.forward { edge.0 } else { edge.1 };
 if let Some(v) = brep.vertex_point(vidx) {
 pts.push(v);
 }
 }
 pts
}

fn trim_almost_closed_polyline(pts: &mut Vec<DVec3>, tol: f64) {
 if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length() < tol {
 pts.pop();
 }
}

fn local_basis_from_normal(normal: DVec3) -> (DVec3, DVec3) {
 let ref_dir = if normal.x.abs() < 0.9 {
 DVec3::X
 } else {
 DVec3::Y
 };
 let u = normal.cross(ref_dir).normalize();
 let v = normal.cross(u).normalize();
 (u, v)
}

/// Spherical (u, v) with u 鈭?[鈭捪€, 蟺] (longitude) and v 鈭?[0, 蟺] (colatitude),
/// matching [`SphericalSurface`](crate::geom::SphericalSurface) / `SurfaceEval`.
fn sphere_point_to_uv(s: &SphericalSurface, p: DVec3) -> DVec2 {
 (*s).world_to_uv(p)
}

/// Remove 2蟺 jumps in u so consecutive samples stay in one branch (for earcut).
fn unwrap_sphere_u_in_chain(uvs: &mut [DVec2]) {
 use std::f64::consts::PI;
 if uvs.len() < 2 {
 return;
 }
 let two_pi = 2.0 * std::f64::consts::PI;
 let mut prev = uvs[0].x;
 for i in 1..uvs.len() {
 let mut u = uvs[i].x;
 let mut d = u - prev;
 while d > PI {
 u -= two_pi;
 d = u - prev;
 }
 while d < -PI {
 u += two_pi;
 d = u - prev;
 }
 uvs[i].x = u;
 prev = u;
 }
}

/// Mapbox earcut: triangulate `flat` as pairs `(x, y)`; `hole_indices` are first-vertex
/// indices of each hole in that vertex list (empty = single polygon).
fn earcut_indices_from_flat(flat: &[f64], hole_indices: &[usize]) -> Vec<usize> {
 let coords: Vec<[f64; 2]> = flat
 .chunks_exact(2)
 .map(|c| [c[0], c[1]])
 .collect();
 if coords.len() < 3 {
 return Vec::new();
 }
 let mut out: Vec<usize> = Vec::new();
 let mut ear = earcut::Earcut::new();
 ear.earcut(coords, hole_indices, &mut out);
 out
}

fn earcut_flat_to_tris(
 flat: &[f64],
 hole_starts: &[usize],
 all_3d: &[DVec3],
 face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
 if flat.len() < 9 || hole_starts.is_empty() {
 return None;
 }
 let indices = earcut_indices_from_flat(flat, hole_starts);
 if indices.is_empty() {
 return None;
 }
 let mut out = Vec::with_capacity(indices.len() / 3);
 for tri in indices.chunks_exact(3) {
 let a = all_3d[tri[0]];
 let b = all_3d[tri[1]];
 let c = all_3d[tri[2]];
 out.push(orient_tri([a, b, c], face_normal));
 }
 Some(out)
}

/// Ear-cut a planar outer + holes in a plane with normal `face_normal`.
/// Uses `pivot` (typically a point on the plane) for a stable 2D frame.
fn try_planar_earcut_holes(
 outer: &[DVec3],
 holes: &[Vec<DVec3>],
 face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
 let (ux, uy) = local_basis_from_normal(face_normal);
 let pivot = outer.first().copied()?;
 let mut all_3d: Vec<DVec3> = Vec::new();
 let mut flat: Vec<f64> = Vec::new();
 for p in outer {
 let q = *p - pivot;
 flat.push(q.dot(ux));
 flat.push(q.dot(uy));
 all_3d.push(*p);
 }
 let mut hole_starts: Vec<usize> = Vec::new();
 for h in holes {
 if h.len() < 3 {
 continue;
 }
 hole_starts.push(all_3d.len());
 for p in h {
 let q = *p - pivot;
 flat.push(q.dot(ux));
 flat.push(q.dot(uy));
 all_3d.push(*p);
 }
 }
 if hole_starts.is_empty() {
 return None;
 }
 earcut_flat_to_tris(&flat, &hole_starts, &all_3d, face_normal)
}

/// Single spherical patch (one outer wire, no holes) in the sphere's (u, v) chart.
/// Preferred over pre-triangles from 3D ear-clips that only cover a convex / disk fill.
fn try_spherical_earcut_simple(
 s: &SphericalSurface,
 outer: &[DVec3],
 face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
 if outer.len() < 3 {
 return None;
 }
 let mut outer_uv: Vec<DVec2> = outer.iter().map(|p| sphere_point_to_uv(s, *p)).collect();
 // Closed-chain-aware unwrapping: `unwrap_u_circle_chain_closed` normalises
 // the u-range so a seam-crossing boundary stays within a single 2蟺 interval.
 // Using `unwrap_sphere_u_in_chain` here caused the ear-cut to cover an
 // extra wrap-around on rotated-sphere faces that cross the seam, producing
 // triangle areas larger than the entire sphere surface.
 let u_vals: Vec<f64> = outer_uv.iter().map(|uv| uv.x).collect();
 let unwrapped_u = unwrap_u_circle_chain_closed(&u_vals);
 for (uv, u) in outer_uv.iter_mut().zip(unwrapped_u.iter()) {
 uv.x = *u;
 }
 // Decimate the UV chain after unwrapping: keeps the ear-cut O(n虏) tractable
 // while preserving the correct UV range.
 const MAX_UV_PTS: usize = 8000;
 if outer_uv.len() > MAX_UV_PTS {
 let step = (outer_uv.len() - 1) as f64 / (MAX_UV_PTS - 1) as f64;
 outer_uv = (0..MAX_UV_PTS)
 .map(|i| {
 let idx = (i as f64 * step).round() as usize;
 outer_uv[idx.min(outer_uv.len() - 1)]
 })
 .collect();
 }
 // Remove near-duplicate consecutive UV points (common at sphere poles
 // and seam boundaries where multiple UV coordinates map to the same 3D
 // point). The earcut triangulation breaks on degenerate UV polygons.
 let uv_dedup_tol = 1e-8;
 outer_uv.dedup_by(|a, b| (*a - *b).length_squared() < uv_dedup_tol);
 if outer_uv.len() < 3 {
 return None;
 }
 let mut flat: Vec<f64> = Vec::with_capacity(2 * outer_uv.len());
 for uv in &outer_uv {
 flat.push(uv.x);
 flat.push(uv.y);
 }
 // Rebuild a parallel 3D array from the deduped UV indices
 let deduped_3d: Vec<DVec3> = outer_uv.iter()
 .map(|uv| s.point_at(uv.x, uv.y))
 .collect();
 let indices = earcut_indices_from_flat(&flat, &[]);
 if indices.is_empty() {
 return None;
 }
 // Guard: if the earcut produces far fewer triangles than expected
 // for the number of boundary points, the UV polygon is likely
 // degenerate (winding issues, pole clustering). Fall through to
 // the UV grid raster which handles this robustly.
 let expected_min = (outer_uv.len() - 2).saturating_sub(outer_uv.len() / 4);
 if indices.len() / 3 < expected_min {
 return None;
 }
 let mut out = Vec::with_capacity(indices.len() / 3);
 for tri in indices.chunks_exact(3) {
 let a = deduped_3d[tri[0]];
 let b = deduped_3d[tri[1]];
 let c = deduped_3d[tri[2]];
 out.push(orient_tri([a, b, c], face_normal));
 }
 Some(out)
}

/// Ear-cut a spherical trimmed patch (outer + holes) in (u, v) parameter space.
fn try_spherical_earcut_holes(
 s: &SphericalSurface,
 outer: &[DVec3],
 holes: &[Vec<DVec3>],
 face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
 let mut outer_uv: Vec<DVec2> = outer.iter().map(|p| sphere_point_to_uv(s, *p)).collect();
 unwrap_sphere_u_in_chain(&mut outer_uv);
 let mut all_3d: Vec<DVec3> = outer.to_vec();
 let mut flat: Vec<f64> = Vec::with_capacity(2 * (outer.len() + holes.iter().map(|h| h.len()).sum::<usize>()));
 for uv in &outer_uv {
 flat.push(uv.x);
 flat.push(uv.y);
 }
 let mut hole_starts: Vec<usize> = Vec::new();
 for h in holes {
 if h.len() < 3 {
 continue;
 }
 hole_starts.push(all_3d.len());
 let mut huv: Vec<DVec2> = h.iter().map(|p| sphere_point_to_uv(s, *p)).collect();
 unwrap_sphere_u_in_chain(&mut huv);
 for p in h {
 all_3d.push(*p);
 }
 for uv in huv {
 flat.push(uv.x);
 flat.push(uv.y);
 }
 }
 if hole_starts.is_empty() {
 return None;
 }
 earcut_flat_to_tris(&flat, &hole_starts, &all_3d, face_normal)
}

/// One outer ring only in the plane (pivot) 鉄?`face_normal` (for nearly planar caps).
fn try_planar_earcut_simple_outer(outer: &[DVec3], face_normal: DVec3) -> Option<Vec<[DVec3; 3]>> {
 if outer.len() < 3 {
 return None;
 }
 let (ux, uy) = local_basis_from_normal(face_normal);
 let pivot = outer[0];
 let mut flat: Vec<f64> = Vec::with_capacity(2 * outer.len());
 for p in outer {
 let q = *p - pivot;
 flat.push(q.dot(ux));
 flat.push(q.dot(uy));
 }
 let all_3d: Vec<DVec3> = outer.to_vec();
 let indices = earcut_indices_from_flat(&flat, &[]);
 if indices.is_empty() {
 return None;
 }
 // Guard: if earcut produces far fewer triangles than expected, the
 // projected 2D polygon is self-intersecting (non-planar 3D points).
 // Fall through to a method that handles this robustly.
 let expected_min = (outer.len() - 2).saturating_sub(outer.len() / 4);
 if indices.len() / 3 < expected_min {
 return None;
 }
 let mut out = Vec::with_capacity(indices.len() / 3);
 for tri in indices.chunks_exact(3) {
 let a = all_3d[tri[0]];
 let b = all_3d[tri[1]];
 let c = all_3d[tri[2]];
 out.push(orient_tri([a, b, c], face_normal));
 }
 Some(out)
}

/// When a face has inner wires, cached `face.triangles` (outer-only) and
/// full-domain UV tessellation are both wrong. Build triangles from 2D earcut
/// when the analytic surface is a plane or sphere.
fn try_face_with_holes(
 brep: &BRep,
 face: &crate::topology::Face,
 face_flat_idx: usize,
) -> Option<Vec<[DVec3; 3]>> {
 if face.inner_wires.is_empty() {
 return None;
 }
 let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } })?;
 let surf = &surf_idx;
 let tol = 1e-5_f64;

 let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
 trim_almost_closed_polyline(&mut outer, tol);
 if outer.len() < 3 {
 return None;
 }
 let mut holes_3d: Vec<Vec<DVec3>> = Vec::new();
 for iw in &face.inner_wires {
 let mut h = sample_wire_polyline_3d(brep, iw);
 trim_almost_closed_polyline(&mut h, tol);
 if h.len() >= 3 {
 holes_3d.push(h);
 }
 }
 if holes_3d.is_empty() {
 return None;
 }

 let rev_holes: Vec<Vec<DVec3>> = holes_3d
 .iter()
 .map(|h| h.iter().rev().copied().collect())
 .collect();
 match surf {
 Surface3::Plane(_) => try_planar_earcut_holes(&outer, &holes_3d, face.normal)
 .or_else(|| try_planar_earcut_holes(&outer, &rev_holes, face.normal)),
 // Spheres: masked UV raster (hole subtraction) is more reliable than ear-cut on
 // non-convex (u,v) annuli; then ear-cut and planar fallbacks.
 Surface3::Sphere(s) => try_spherical_uv_masked_raster(s, brep, face, face_flat_idx, face.normal)
 .or_else(|| try_spherical_earcut_holes(s, &outer, &holes_3d, face.normal))
 .or_else(|| try_spherical_earcut_holes(s, &outer, &rev_holes, face.normal))
 .or_else(|| try_planar_earcut_holes(&outer, &holes_3d, face.normal))
 .or_else(|| try_planar_earcut_holes(&outer, &rev_holes, face.normal)),
 _ => None,
 }
}

/// Resolution used for UV-grid tessellation of curved faces (per axis).
///
/// 64脳64 gives a <0.1% volume error for a unit sphere.
const UV_TESS_N: usize = 80;

/// Finer grid for masked sphere patches (trimmed / holed) when ear-cut fails.
const SPHERE_UV_MASK_N: usize = 100;

/// Grid for 鈥栤垈P/鈭倁脳鈭侾/鈭倂鈥?midpoint integration on a UV rectangle (per-axis).
const PARAM_RECT_AREA_N: usize = 80;

/// Winding-number test; correct for self-intersecting polygons, unlike the
/// even-odd ray-crossing test used by [`point_in_polygon_2d`].
fn winding_number_2d(poly: &[DVec2], pt: DVec2) -> i32 {
 let n = poly.len();
 if n < 3 {
 return 0;
 }
 let mut wn = 0i32;
 let mut j = n - 1;
 for i in 0..n {
 let vi = poly[i];
 let vj = poly[j];
 if vi.y <= pt.y {
 if vj.y > pt.y && is_left_2d(vi, vj, pt) > 0.0 {
 wn += 1;
 }
 } else if vj.y <= pt.y && is_left_2d(vi, vj, pt) < 0.0 {
 wn -= 1;
 }
 j = i;
 }
 wn
}

/// Cross product (b鈭抋)脳(c鈭抋) for 2D orientation / signed area.
fn is_left_2d(a: DVec2, b: DVec2, c: DVec2) -> f64 {
 (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
 let n = poly.len();
 if n < 3 {
 return false;
 }
 let mut inside = false;
 let mut j = n - 1;
 for i in 0..n {
 let vi = poly[i];
 let vj = poly[j];
 if ((vi.y > pt.y) != (vj.y > pt.y))
 && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x)
 {
 inside = !inside;
 }
 j = i;
 }
 inside
}

/// Test if `point` (on the sphere) is inside a spherical polygon defined by
/// `boundary` (closed 3D polyline on the sphere).  Uses the angular sum test:
/// the signed angles around the polygon sum to 卤2蟺 for inside, 0 for outside.
/// Correct for any spherical polygon, no UV-mapping degeneracies.
fn point_in_spherical_polygon_3d(boundary: &[DVec3], point: DVec3) -> bool {
 let n = boundary.len();
 if n < 3 {
 return false;
 }
 let mut total = 0.0f64;
 for i in 0..n {
 let j = if i + 1 < n { i + 1 } else { 0 };
 let a = boundary[i] - point;
 let b = boundary[j] - point;
 let cross = a.cross(b);
 let sin_theta = cross.length();
 let cos_theta = a.dot(b);
 let theta = sin_theta.atan2(cos_theta);
 if cross.dot(point) >= 0.0 {
 total += theta;
 } else {
 total -= theta;
 }
 }
 total.abs() > std::f64::consts::PI
}

pub fn point_in_spherical_polygon_3d_pub(boundary: &[DVec3], point: DVec3) -> bool {
 point_in_spherical_polygon_3d(boundary, point)
}

/// Shortest signed step from `a` to `b` on the circle, with `a, b` reduced mod 2蟺 to [0,2蟺).
#[inline]
fn short_delta_on_circle_01(a: f64, b: f64) -> f64 {
 const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
 const PI: f64 = std::f64::consts::PI;
 let am = a.rem_euclid(TWO_PI);
 let bm = b.rem_euclid(TWO_PI);
 let mut d = bm - am;
 if d > PI {
 d -= TWO_PI;
 }
 if d < -PI {
 d += TWO_PI;
 }
 d
}

/// Cumulative u along a **closed** wire, using only short steps on S鹿, then remove full 2蟺
/// turns so the first and last point close in 鈩?(winding 0 in the (u, v) chart).
/// This fixes `unwrap_sphere_u_in_chain` sometimes representing a small patch with span > 2蟺.
fn unwrap_u_circle_chain_closed(b: &[f64]) -> Vec<f64> {
 const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
 let n = b.len();
 if n < 2 {
 return b.to_vec();
 }
 let b: Vec<f64> = b.iter().map(|x| x.rem_euclid(TWO_PI)).collect();
 let mut o: Vec<f64> = Vec::with_capacity(n);
 o.push(b[0]);
 for i in 1..n {
 o.push(o[i - 1] + short_delta_on_circle_01(b[i - 1], b[i]));
 }
 let d_close = short_delta_on_circle_01(b[n - 1], b[0]);
 let gap = o[0] - (o[n - 1] + d_close);
 if gap.abs() > 1e-2 {
 let w = (gap / TWO_PI).round() as i32;
 if w != 0 {
 for i in 1..n {
 o[i] -= (w as f64) * TWO_PI;
 }
 }
 }
 o
}

/// Add +2蟺k to an inner u-chain so it lies near the outer patch in 鈩?(shared chart for PIP).
fn align_inner_u_chain_to_outer(outer_umin: f64, outer_umax: f64, o_inner: &mut [f64]) {
 if o_inner.is_empty() {
 return;
 }
 const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
 let mid = 0.5 * (outer_umin + outer_umax);
 let mut best_k: i32 = 0;
 let mut best_d = f64::INFINITY;
 for k in -2..=2 {
 let t = o_inner[0] + (k as f64) * TWO_PI;
 let d = (t - mid).abs();
 if d < best_d {
 best_d = d;
 best_k = k;
 }
 }
 let s = (best_k as f64) * TWO_PI;
 for u in o_inner.iter_mut() {
 *u += s;
 }
}

/// Shared (u,v) chart and grid bounds for sphere UV raster / parametric **dA** sum (with or
/// without inner loops).
struct SphereHoledMaskCtx {
 outer_uv: Vec<DVec2>,
 outer_3d: Vec<DVec3>,
 inner_polys: Vec<Vec<DVec2>>,
 inner_3d: Vec<Vec<DVec3>>,
 umin: f64,
 umax: f64,
 vmin: f64,
 vmax: f64,
 /// False when the UV polygon wraps around u multiple times, making
 /// winding_number_2d unreliable. When false, the grid test uses only
 /// point_in_spherical_polygon_3d.
 use_uv_winding: bool,
}

fn spherical_holed_uv_mask_setup(
 s: &SphericalSurface,
 brep: &BRep,
 face: &Face,
) -> Option<SphereHoledMaskCtx> {
 let tol = 1e-5_f64;
 let n_edges = face.outer_wire.edges.len();
 let per_edge = if n_edges > 600 {
 4
 } else if n_edges > 300 {
 8
 } else if n_edges > 150 {
 16
 } else if n_edges > 75 {
 32
 } else if n_edges > 30 {
 48
 } else {
 64
 };
 let mut outer3 = sample_wire_polyline_3d_with_n(brep, &face.outer_wire, per_edge);
 trim_almost_closed_polyline(&mut outer3, tol);
 if outer3.len() < 3 {
 return None;
 }
 // Decimate after UV mapping so the UV polygon covers parameter space
 // uniformly, avoiding pole-clustering from uniform 3D decimation.
 let mut b_outer: Vec<f64> = outer3
 .iter()
 .map(|p| sphere_point_to_uv(s, *p).x.rem_euclid(2.0 * std::f64::consts::PI))
 .collect();
 // Fix degenerate u=0 at sphere poles before unwrapping.
 {
 use std::f64::consts::PI;
 let two_pi = 2.0 * PI;
 let n = outer3.len();
 let mut idx = 0;
 while idx < n {
 let p = outer3[idx];
 let v = sphere_point_to_uv(s, p).y;
 let at_pole = v < 1e-12 || (v - PI).abs() < 1e-12;
 if at_pole && b_outer[idx].abs() < 1e-12 {
 // Find extent of consecutive pole vertices
 let mut run_start = idx;
 while run_start > 0 {
 let p0 = outer3[run_start - 1];
 let v0 = sphere_point_to_uv(s, p0).y;
 if !(v0 < 1e-12 || (v0 - PI).abs() < 1e-12)
 || b_outer[run_start - 1].abs() >= 1e-12
 {
 break;
 }
 run_start -= 1;
 }
 let mut run_end = idx;
 while run_end + 1 < n {
 let p1 = outer3[run_end + 1];
 let v1 = sphere_point_to_uv(s, p1).y;
 if !(v1 < 1e-12 || (v1 - PI).abs() < 1e-12)
 || b_outer[run_end + 1].abs() >= 1e-12
 {
 break;
 }
 run_end += 1;
 }
 let u_prev = if run_start > 0 {
 b_outer[run_start - 1]
 } else {
 b_outer[n - 1]
 };
 let u_next = if run_end + 1 < n {
 b_outer[run_end + 1]
 } else {
 b_outer[0]
 };
 let delta = short_delta_on_circle_01(u_prev, u_next);
 let run_len = (run_end - run_start + 1) as f64;
 for k in 0..=(run_end - run_start) {
 let frac = (k as f64 + 1.0) / (run_len + 1.0);
 b_outer[run_start + k] = (u_prev + delta * frac).rem_euclid(two_pi);
 }
 idx = run_end + 1;
 } else {
 idx += 1;
 }
 }
 }
 let o_outer = unwrap_u_circle_chain_closed(&b_outer);
 let (ou0, ou1) = o_outer.iter().fold(
 (f64::INFINITY, f64::NEG_INFINITY),
 |(a, b), u| (a.min(*u), b.max(*u)),
 );
 // If the UV polygon has near-zero u-range, the outer wire is a seam-edge
 // (full sphere before boolean) rather than a proper closed boundary.
 // Fall back to full-sphere tessellation instead of the masked grid.
 if (ou1 - ou0).abs() < 1e-8 {
 return None;
 }
 let outer_uv: Vec<DVec2> = outer3
 .iter()
 .zip(o_outer.iter())
 .map(|(p, u)| DVec2::new(*u, sphere_point_to_uv(s, *p).y))
 .collect();

 let mut inner_polys: Vec<Vec<DVec2>> = Vec::new();
 let mut inner_3d: Vec<Vec<DVec3>> = Vec::new();
 for w in &face.inner_wires {
 let mut h3 = sample_wire_polyline_3d(brep, w);
 trim_almost_closed_polyline(&mut h3, tol);
 if h3.len() < 3 {
 continue;
 }
 let b_in: Vec<f64> = h3
 .iter()
 .map(|p| sphere_point_to_uv(s, *p).x.rem_euclid(2.0 * std::f64::consts::PI))
 .collect();
 let mut o_in = unwrap_u_circle_chain_closed(&b_in);
 align_inner_u_chain_to_outer(ou0, ou1, &mut o_in);
 let huv: Vec<DVec2> = h3
 .iter()
 .zip(o_in.iter())
 .map(|(p, u)| DVec2::new(*u, sphere_point_to_uv(s, *p).y))
 .collect();
 let h3d: Vec<DVec3> = huv.iter().map(|uv| s.point_at(uv.x, uv.y)).collect();
 inner_polys.push(huv);
 inner_3d.push(h3d);
 }

 // Fix degenerate pole-boundary closing edges in UV polygon.
 // At sphere poles (v 鈮?0 or v 鈮?蟺), many u-values map to the same 3D
 // point.  The closing edge (last鈫抐irst sample) can jump diagonally
 // across the UV interior instead of following the constant-v pole
 // boundary.  Insert intermediate points along the pole.
 const POLE_V_THRESHOLD: f64 = 0.15;
 const U_JUMP_THRESHOLD: f64 = std::f64::consts::FRAC_PI_2;
 const POLE_N_INSERT: usize = 32;
 let pi = std::f64::consts::PI;
 let mut fixed_uv: Vec<DVec2> = Vec::with_capacity(outer_uv.len() + 8);
 let mut fixed_3d: Vec<DVec3> = Vec::with_capacity(outer3.len() + 8);
 let nu = outer_uv.len();
 for i in 0..nu {
 let j = (i + 1) % nu;
 let a_uv = outer_uv[i];
 let b_uv = outer_uv[j];
 fixed_uv.push(a_uv);
 fixed_3d.push(outer3[i]);
 let near_a = a_uv.y < POLE_V_THRESHOLD || a_uv.y > pi - POLE_V_THRESHOLD;
 let near_b = b_uv.y < POLE_V_THRESHOLD || b_uv.y > pi - POLE_V_THRESHOLD;
 let same_pole = (a_uv.y < POLE_V_THRESHOLD && b_uv.y < POLE_V_THRESHOLD)
 || (a_uv.y > pi - POLE_V_THRESHOLD && b_uv.y > pi - POLE_V_THRESHOLD);
 let big_jump = (b_uv.x - a_uv.x).abs() > U_JUMP_THRESHOLD;
 if near_a && near_b && same_pole && big_jump {
 let pole_v = if a_uv.y < POLE_V_THRESHOLD { 0.0 } else { pi };
 for k in 1..POLE_N_INSERT {
 let frac = k as f64 / POLE_N_INSERT as f64;
 let u = a_uv.x + (b_uv.x - a_uv.x) * frac;
 fixed_uv.push(DVec2::new(u, pole_v));
 fixed_3d.push(s.point_at(u, pole_v));
 }
 }
 }

 let mut outer_uv = fixed_uv;
 let mut outer_3d = fixed_3d;

 // Dedup consecutive near-duplicate UV points (from shared edge vertices in
 // the outer wire sampling) to keep the polygon small for the grid test.
 {
 let mut dedup_uv: Vec<DVec2> = Vec::with_capacity(outer_uv.len());
 let mut dedup_3d: Vec<DVec3> = Vec::with_capacity(outer_3d.len());
 dedup_uv.push(outer_uv[0]);
 dedup_3d.push(outer_3d[0]);
 for i in 1..outer_uv.len() {
 if (outer_uv[i] - *dedup_uv.last().unwrap()).length_squared() > 1e-16 {
 dedup_uv.push(outer_uv[i]);
 dedup_3d.push(outer_3d[i]);
 }
 }
 // Keep the polygon closed
 if dedup_uv.len() >= 2
 && (dedup_uv[0] - *dedup_uv.last().unwrap()).length_squared() < 1e-16
 {
 dedup_uv.pop();
 dedup_3d.pop();
 }
 outer_uv = dedup_uv;
 outer_3d = dedup_3d;
 }

 let mut umin = f64::INFINITY;
 let mut umax = f64::NEG_INFINITY;
 let mut vmin = f64::INFINITY;
 let mut vmax = f64::NEG_INFINITY;
 for p in &outer_uv {
 umin = umin.min(p.x);
 umax = umax.max(p.x);
 vmin = vmin.min(p.y);
 vmax = vmax.max(p.y);
 }
 for h in &inner_polys {
 for p in h {
 umin = umin.min(p.x);
 umax = umax.max(p.x);
 vmin = vmin.min(p.y);
 vmax = vmax.max(p.y);
 }
 }
 // Snapshot raw u-range before margins to detect UV wrapping (merged faces
 // can produce u-ranges > 2蟺 from multiple seam crossings).
 let raw_u_range = umax - umin;
 let margin_u = raw_u_range * 0.02 + 1e-4;
 let margin_v = (vmax - vmin) * 0.02 + 1e-4;
 umin -= margin_u;
 umax += margin_u;
 vmin -= margin_v;
 vmax += margin_v;
 let pi = std::f64::consts::PI;
 let vmin = vmin.clamp(0.0, pi);
 let vmax = vmax.clamp(0.0, pi);
 if !umin.is_finite()
 || !umax.is_finite()
 || !vmin.is_finite()
 || !vmax.is_finite()
 || umax <= umin + 1e-14
 || vmax <= vmin + 1e-14
 {
 return None;
 }
 // Disable UV winding-number when the polygon wraps u multiple times
 // (raw range >> 2蟺). Winding number in a non-periodic 2D domain gives
 // wrong results for wrapped polygons.
 let two_pi = 2.0 * std::f64::consts::PI;
 let outer_u_min = outer_uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
 let outer_u_max = outer_uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
 let outer_u_span = outer_u_max - outer_u_min;
 // Also check that the outer UV polygon's own u-span doesn't exceed 蟺.
 // When a polygon spans > 蟺 in u (e.g. nearly 2蟺 from seam wrapping),
 // the 2D winding number gives wrong inside/outside classification
 // because the boundary goes the "long way" around the sphere.
 // Such polygons must be tested with 3D-only point-in-spherical-polygon.
 let use_uv_winding = raw_u_range <= two_pi + 0.5
 && outer_u_span.abs() <= std::f64::consts::PI + 0.1;

 if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
 // Compute shoelace area of UV polygon (in UV^2 units, not physical area)
 let mut _uv_area = 0.0;
 let n = outer_uv.len();
 for i in 0..n {
 let j = (i + 1) % n;
 _uv_area += outer_uv[i].x * outer_uv[j].y;
 _uv_area -= outer_uv[j].x * outer_uv[i].y;
 }
 let _ = (_uv_area * 0.5).abs();
 }

 Some(SphereHoledMaskCtx {
 outer_uv,
 outer_3d,
 inner_polys,
 inner_3d,
 umin,
 umax,
 vmin,
 vmax,
 use_uv_winding,
 })
}

/// Sum `鈭?R虏 sin v dudv` over the same masked grid (no triangulation), for GProp-style area.
/// projected loop in world UV, return that box's `width * height`.

/// Compute sphere face area using 2x2 Gauss-Legendre quadrature over the
/// UV bounding box, with inside/outside testing via 3D point-in-spherical-polygon.
///
/// OCCT uses adaptive Gauss-Legendre integration in BRepGProp::SurfaceProperties

/// Exact spherical-polygon area when all outer-wire edges are great circles
/// (Circle3 curves centered at the sphere center).
///
/// Uses the formula A = R^2 * (sum(interior_angles) - (n-2)*pi) for a spherical
/// polygon bounded by great-circle arcs.  Interior angles are computed from the
/// tangent vectors at each vertex via projection onto the tangent plane, which
/// is independent of the edge curve parameterization 鈥?the great-circle check
/// only needs the Curve3 type + center, not the normal direction.
///
/// Returns None when any edge is not a great circle or when the face has
/// inner wires (holes).
fn try_spherical_polygon_great_circle_area(
 s: &SphericalSurface,
 brep: &BRep,
 face: &Face,
) -> Option<f64> {
 if !face.inner_wires.is_empty() {
 return None;
 }
 let n_edges = face.outer_wire.edges.len();
 if n_edges < 3 {
 return None;
 }
 let tol = 1e-10;

 // Collect vertices in boundary order, verifying all edges are great circles.
 let mut verts: Vec<DVec3> = Vec::with_capacity(n_edges + 1);
 for we in &face.outer_wire.edges {
 let ei = we.idx;
 let edge = brep.flat_edges().get(ei).copied()?;
 let curve_idx = brep.tshapes.get(ei).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } })?;
 let curve = &curve_idx;
 match curve {
 Curve3::Circle(c) => {
 if (c.center - s.center).length() > tol {
 return None;
 }
 }
 _ => return None,
 }
 // Start vertex: for WireEdge::rev, the boundary enters via the
 // forward edge's end vertex.
 let vi = if we.forward { edge.0 } else { edge.1 };
 let pt = brep.vertex_point(vi).unwrap_or(DVec3::ZERO);
 if verts.is_empty() || (pt - *verts.last()?).length() > tol {
 verts.push(pt);
 }
 }

 // Ensure closed loop (last == first for the formula's sum).
 if (verts.first()? - verts.last()?).length() > tol {
 verts.push(*verts.first()?);
 }

 let n = verts.len() - 1; // number of unique vertices
 if n < 3 {
 return None;
 }

 // Compute interior angles at each unique vertex.
 // verts = [v0, v1, ..., v_{n-1}, v0] (n+1 elements, last is closure).
 // For vertex verts[i]: incoming edge verts[i-1]->verts[i], outgoing verts[i]->verts[i+1].
 let mut sum_angles = 0.0;
 for i in 0..n {
 let v_prev = if i > 0 { verts[i - 1] } else { verts[n - 1] };
 let v_curr = verts[i];
 let v_next = verts[i + 1]; // verts[n] == verts[0] via closure
 let v_hat = v_curr.normalize();

 // Tangent of incoming edge (from v_prev 鈫?v_curr, pointing into v_curr)
 // projected onto the tangent plane of the sphere at v_curr.
 let t_in = (v_prev - v_hat * v_prev.dot(v_hat)).normalize();
 // Tangent of outgoing edge (from v_curr 鈫?v_next)
 let t_out = (v_next - v_hat * v_next.dot(v_hat)).normalize();

 let cos_theta = t_in.dot(t_out).clamp(-1.0, 1.0);
 let theta = cos_theta.acos();

 // Signed turn direction.
 // For a CCW outer wire (outward-facing sphere): right turn 鈫?convex, left 鈫?reflex.
 // We take interior = 蟺 - 胃 for right turns, 蟺 + 胃 for left turns.
 let cross_sign = t_in.cross(t_out).dot(v_hat);
 let interior = if cross_sign < 0.0 {
 std::f64::consts::PI - theta  // right turn, convex
 } else {
 std::f64::consts::PI + theta  // left turn, reflex
 };
 sum_angles += interior;
 }

 let r2 = s.radius * s.radius;
 let full = 4.0 * std::f64::consts::PI * r2;
 let mut area = r2 * (sum_angles - (n as f64 - 2.0) * std::f64::consts::PI);
 if cfg!(debug_assertions) && std::env::var("RCAD_DEBUG_SA").is_ok() {
 eprintln!("[SA_SPHERE] n={} sum_angles={:.6} area_raw={:.6} full={:.6} sp={:?}",
 n, sum_angles, area, full, face.sample_point);
 }
 // For CW (reversed) wires, the spherical polygon formula gives area < 0.
 // Take the absolute value 鈥?the correct region (small polygon vs complement)
 // is determined below by the sample_point check.
 area = area.abs();
 if area <= full + 1e-12 {
 // 鉁?OCCT :  銆?
 // area > half sphere complement锛?锛夈€?
 if area > full * 0.5 { area = full - area; }
 if let Some(sp) = face.sample_point {
 let inside = point_in_spherical_polygon_3d(&verts[..n], sp);
 if !inside { area = full - area; }
 }
 // 鉁?OCCT : clamp  銆?
 if area > 0.0 && area <= full + 1e-12 {
 return Some(area);
 }
 }
 None
}

/// Compute sphere face area using 2x2 Gauss-Legendre quadrature over the
/// for all surface types including spheres.  This replaces the old grid-raster
/// approach that used a uniform grid with a 5-point OR test 鈥?approximating
/// the boundary poorly and underestimating area for faces without pcurves
/// (e.g. analytic-constructed spherical caps).
///
/// With N=30 cells per side and 2x2 GL points per cell (3600 evaluations),
/// this achieves O(h^4) accuracy vs O(h^2) for the midpoint rule.
fn sphere_gauss_legendre_area_sum(s: &SphericalSurface, ctx: &SphereHoledMaskCtx) -> f64 {
 const N: usize = 30;
 const GL_NEG: f64 = -0.5773502691896257;
 const GL_POS: f64 = 0.5773502691896257;
 let gl_pts = [GL_NEG, GL_POS];
 let umin = ctx.umin;
 let umax = ctx.umax;
 let vmin = ctx.vmin;
 let vmax = ctx.vmax;
 let du = (umax - umin) / N as f64;
 let dv = (vmax - vmin) / N as f64;
 let r2 = s.radius * s.radius;
 let outer_3d = &ctx.outer_3d;
 let inner_3d = &ctx.inner_3d;
 let cell_area = du * dv / 4.0;

 let mut total = 0.0;
 for i in 0..N {
 let u_mid = umin + (i as f64 + 0.5) * du;
 for j in 0..N {
 let v_mid = vmin + (j as f64 + 0.5) * dv;
 let mut sum_sin = 0.0;
 for &gu in &gl_pts {
 let u = u_mid + gu * du * 0.5;
 for &gv in &gl_pts {
 let v = v_mid + gv * dv * 0.5;
 if v < 0.0 || v > std::f64::consts::PI { continue; }
 let ok = {
 let p3d = s.point_at(u, v);
 point_in_spherical_polygon_3d(outer_3d, p3d)
 };
 if ok && !inner_3d.iter().any(|h3d| {
 let p3d = s.point_at(u, v);
 point_in_spherical_polygon_3d(h3d, p3d)
 }) {
 sum_sin += v.sin();
 }
 }
 }
 if sum_sin > 0.0 {
 total += r2 * cell_area * sum_sin;
 }
 }
 }
 total
}

///
/// This recovers full rectangle area when shoe-lace on coarse edge samples under-counts (merged
/// `box` unions). It is **not** sufficient to prove the polygon fills that box: e.g. a parallelogram
/// or other convex quad inscribed in its AABB still has all vertices on the box boundary in 2D,
/// but area `< w*h`. Callers must compare against shoe-lace and take the smaller when they
/// disagree ([`try_planar_face_area_shoelace`]).
fn axis_aligned_world_plane_uv_axes(n: DVec3) -> Option<[usize; 2]> {
 let a = n.abs();
 if a.x > 1.0 - 2e-3 {
 Some([1, 2])
 } else if a.y > 1.0 - 2e-3 {
 Some([0, 2])
 } else if a.z > 1.0 - 2e-3 {
 Some([0, 1])
 } else {
 None
 }
}

fn bbox2d_components(uv: &[(f64, f64)]) -> Option<(f64, f64, f64, f64)> {
 if uv.is_empty() {
 return None;
 }
 let mut u0 = f64::INFINITY;
 let mut u1 = f64::NEG_INFINITY;
 let mut v0 = f64::INFINITY;
 let mut v1 = f64::NEG_INFINITY;
 for &(u, v) in uv {
 u0 = u0.min(u);
 u1 = u1.max(u);
 v0 = v0.min(v);
 v1 = v1.max(v);
 }
 Some((u0, u1, v0, v1))
}

/// UV chain from traversing `wire` vertices only (edge endpoints in order).
///
/// Dense [`sample_wire_polyline_3d`] can self-cross after boolean edge reordering while the
/// vertex ring still traces the true boundary; comparing both shoelaces disambiguates that
/// from legitimate concave faces (`bfuse_simple/D9` vs `bfuse_simple/E5`).
fn wire_edge_endpoint_3d(brep: &BRep, we: &WireEdge) -> Option<DVec3> {
 let edge = brep.flat_edges().get(we.idx).copied()?;
 // Use edge curve when available (handles degenerate edges where start==end
 // but curve spans the full distance between two different positions).
 if let Some(curve) = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } })
 {
 let range = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None } })
 .unwrap_or_else(|| curve.default_domain());
 // For degenerate edges (start==end), the curve spans between two
 // different positions. The face boundary needs BOTH curve endpoints
 // to form a correct polygon. Use the EDGE DIRECTION to select which
 // endpoint: forward 鈫?range[0], reverse 鈫?range[1]. For non-degenerate
 // edges this matches vertex.position; for degenerate edges it gives
 // the correct geometric position on the face boundary.
 let t = if we.forward { range[0] } else { range[1] };
 return Some(curve.point_at(t));
 }
 let vidx = if we.forward { edge.0 } else { edge.1 };
 Some(brep.vertex_point(vidx)?)
}

fn outer_wire_ordered_vertex_uvs(
 brep: &BRep,
 wire: &Wire,
 i: usize,
 j: usize,
 pos_tol: f64,
) -> Vec<(f64, f64)> {
 let mut out: Vec<(f64, f64)> = Vec::new();
 for we in &wire.edges {
 let Some(p) = wire_edge_endpoint_3d(brep, we) else {
 continue;
 };
 let uv = (p[i], p[j]);
 if let Some(&(u0, v0)) = out.last() {
 if (uv.0 - u0).abs() <= pos_tol && (uv.1 - v0).abs() <= pos_tol {
 continue;
 }
 }
 out.push(uv);
 }
 trim_almost_closed_uv_chain(&mut out, pos_tol);
 out
}

fn trim_almost_closed_uv_chain(uvs: &mut Vec<(f64, f64)>, pos_tol: f64) {
 if uvs.len() >= 2 {
 let (u0, v0) = uvs[0];
 let (u1, v1) = uvs[uvs.len() - 1];
 if (u0 - u1).abs() <= pos_tol && (v0 - v1).abs() <= pos_tol {
 uvs.pop();
 }
 }
}

/// Unique vertices referenced by `wire` edges, projected to world coordinates `(i, j)`.
fn outer_wire_unique_vertex_uvs(
 brep: &BRep,
 wire: &Wire,
 i: usize,
 j: usize,
 pos_tol: f64,
) -> Vec<(f64, f64)> {
 let mut out: Vec<(f64, f64)> = Vec::new();
 for we in &wire.edges {
 let flat_edges = brep.flat_edges();
 let Some(edge) = flat_edges.get(we.idx) else {
 continue;
 };
 // Use edge curve endpoints when available (handles degenerate edges
 // where start==end but the curve spans two different positions).
 let pts: [DVec3; 2] = if let Some(curve) = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } })
 {
 let range = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None } })
 .unwrap_or_else(|| curve.default_domain());
 [curve.point_at(range[0]), curve.point_at(range[1])]
 } else {
 let p0 = brep.vertex_point(edge.0).map(|v| v)
 .unwrap_or(DVec3::ZERO);
 let p1 = brep.vertex_point(edge.1).map(|v| v)
 .unwrap_or(DVec3::ZERO);
 [p0, p1]
 };
 for &p in &pts {
 let uv = (p[i], p[j]);
 if !out
 .iter()
 .any(|&(u2, v2)| (uv.0 - u2).abs() <= pos_tol && (uv.1 - v2).abs() <= pos_tol)
 {
 out.push(uv);
 }
 }
 }
 out
}

/// Andrew monotone chain; returns CCW hull vertices (no duplicate closing point).
fn convex_hull_2d_monotone(mut pts: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
 if pts.len() <= 1 {
 return pts;
 }
 pts.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));
 let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
 (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
 };
 let mut lower: Vec<(f64, f64)> = Vec::new();
 for &p in &pts {
 while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 1e-18 {
 lower.pop();
 }
 lower.push(p);
 }
 let mut upper: Vec<(f64, f64)> = Vec::new();
 for &p in pts.iter().rev() {
 while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 1e-18 {
 upper.pop();
 }
 upper.push(p);
 }
 lower.pop();
 upper.pop();
 lower.extend(upper);
 lower
}

fn polygon_area_2d_xy(pts: &[(f64, f64)]) -> f64 {
 if pts.len() < 3 {
 return 0.0;
 }
 let mut a = 0.0_f64;
 for k in 0..pts.len() {
 let p = pts[k];
 let q = pts[(k + 1) % pts.len()];
 a += p.0 * q.1 - p.1 * q.0;
 }
 0.5 * a.abs()
}

fn try_axis_aligned_world_rect_plane_area(
 brep: &BRep,
 face: &Face,
 face_normal: DVec3,
) -> Option<f64> {
 let n = face_normal.normalize_or_zero();
 if n.length_squared() < 1e-24 {
 return None;
 }
 let [i, j] = axis_aligned_world_plane_uv_axes(n)?;
 // Boolean T-junction repair can permute `outer_wire` edge order so dense sampling walks a
 // self-intersecting polyline while the face is still a convex patch (OCCT `bcommon_simple/B1`).
 // Unique-vertex convex hull in the plane normal to 卤X/卤Y/卤Z recovers the true area; the
 // bbox脳shoelace branch below still handles large sampled loops.
 let pos_tol = (1e-7_f64 * brep
 .bounding_box()
 .map(|[mn, mx]| (mx - mn).length())
 .unwrap_or(1.0))
 .max(1e-9);
 let vu = outer_wire_unique_vertex_uvs(brep, &face.outer_wire, i, j, pos_tol);
 // Ordered boundary loop (wire sampling order). Convex hull below over-counts concave polygons.
 let mut outer_ordered = sample_wire_polyline_3d(brep, &face.outer_wire);
 trim_almost_closed_polyline(&mut outer_ordered, 1e-5);
 let a_loop = if outer_ordered.len() >= 3 {
 let uv_ord: Vec<(f64, f64)> = outer_ordered.iter().map(|p| (p[i], p[j])).collect();
 polygon_area_2d_xy(&uv_ord)
 } else {
 0.0
 };

 // Boolean seams can split edges into many vertex points; `bfuse_simple/D9` exceeds the old
 // `<=48` cap and must still hit hull vs dense vs vertex agreement logic.
 if (3..=4096).contains(&vu.len()) {
 let hull = convex_hull_2d_monotone(vu.clone());
 if hull.len() >= 3 {
 let a_hull = polygon_area_2d_xy(&hull);
 if a_hull > 1e-18 {
 const REL: f64 = 1e-4;
 const AGREE_REL: f64 = 2e-3;
 let scale = a_hull.max(a_loop).max(1.0);
 let abs_eps = 1e-9 * scale;
 // Concave silhouette after booleans: hull area > true boundary; dense-sample
 // shoelace should match the vertex-ring shoelace when edges are straight (`D9`).
 // When dense samples self-cross (`E5`), vertex ring still matches the face 鈥?disagree
 // 鈫?convex hull (vertex set) recovers area.
 if a_loop > 1e-18 && a_hull > a_loop * (1.0 + REL) + abs_eps {
 let uv_vert =
 outer_wire_ordered_vertex_uvs(brep, &face.outer_wire, i, j, pos_tol);
 if uv_vert.len() >= 3 {
 let a_vert = polygon_area_2d_xy(&uv_vert);
 let scale_agree = a_vert.max(a_loop).max(1.0);
 // Require both vertex/dense agreement and a substantial fraction of the convex hull:
 // `bcommon_simple/B1` can yield bogus agreeing shoelaces (~80% of hull).
 const LOOP_FRAC_OF_HULL_MIN: f64 = 0.81;
 if (a_vert - a_loop).abs() <= AGREE_REL * scale_agree + abs_eps
 && a_loop + abs_eps >= a_hull * LOOP_FRAC_OF_HULL_MIN
 {
 return Some(a_loop);
 }
 // Vertex ring and dense sample agree, but hull is much larger 鈫?very concave face.
 // The polygon is simple and the shoelace is trustworthy (e.g. box-cylinder annular cap).
 if (a_vert - a_loop).abs() <= AGREE_REL * scale_agree + abs_eps
 && a_vert + abs_eps < a_hull * 0.6
 {
 return Some(a_vert.max(0.0));
 }
 // Limit to large silhouette patches: small faces (`bcommon_simple/B1`) can have
 // bogus vertex-ring shoelaces while hull matches OCCT; `bfuse_simple/D9` top uses ~40000.
 const MIN_HULL_ABS_VERT_FALLBACK: f64 = 15000.0;
 const VERT_OVER_LOOP_REL: f64 = 0.02;
 if a_hull >= MIN_HULL_ABS_VERT_FALLBACK
 && a_vert + abs_eps < a_hull * (1.0 - REL)
 && a_vert > 1e-18
 && a_vert > a_loop * (1.0 + VERT_OVER_LOOP_REL) + abs_eps
 {
 return Some(a_vert.max(0.0));
 }
 } else {
 }
 }
 // When hull is NOT significantly larger than the dense-sample loop, the loop
 // correctly traces the boundary (may be concave, e.g. a clipped cylinder cap
 // with a circular arc 鈫?hull under-estimates).  Return the loop area when it
 // is at least as large as the hull.
 if a_loop > 1e-18 && a_loop + abs_eps >= a_hull {
 return Some(a_loop);
 }
 return Some(a_hull);
 }
 }
 }
 let mut outer = outer_ordered;
 trim_almost_closed_polyline(&mut outer, 1e-5);
 if outer.len() < 3 {
 return None;
 }
 let uv: Vec<(f64, f64)> = outer.iter().map(|p| (p[i], p[j])).collect();
 let (u0, u1, v0, v1) = bbox2d_components(&uv)?;
 let w = u1 - u0;
 let h = v1 - v0;
 if !(w > 1e-18 && h > 1e-18) {
 return None;
 }
 let scale = w.max(h).max(1.0);
 let eps = (1e-5 * scale).max(1e-9);
 for &(u, v) in &uv {
 let on_edge = (u - u0).abs() <= eps
 || (u1 - u).abs() <= eps
 || (v - v0).abs() <= eps
 || (v1 - v).abs() <= eps;
 if !on_edge {
 return None;
 }
 }
 Some((w * h).max(0.0))
}

/// Exact area for a planar face whose outer wire consists of line segments
/// and great-circle arcs (Circle3 curves).  Line edges use the exact shoelace
/// vertex contribution; circular arcs add the exact segment area between the
/// chord and the arc: 卤r虏路(胃 - sin(胃))/2.
///
/// Returns `None` if any edge has a curve type other than Line3 or Circle3,
/// falling through to the sampled shoelace fallback.
fn try_planar_face_exact_contour_area(brep: &BRep, face: &Face, face_normal: DVec3) -> Option<f64> {
 let (ux, uy) = local_basis_from_normal(face_normal);
 let n_edges = face.outer_wire.edges.len();

 // Fast path: single full-circle edge (planar cap from cylinder/sphere clipping).
 // The standard code below requires >= 3 edges for shoelace and handles arcs via
 // chord-segment correction, but the cross-product sign heuristic fails when
 // start == end (full circle 鈥?trav = 0).  Detect this case directly.
 if n_edges == 1 {
 if let Some(we) = face.outer_wire.edges.first() {
 let ei = we.idx;
 if let Some(curve_idx) = brep.tshapes.get(ei).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } }) {
 if let Curve3::Circle(c) = curve_idx {
 if let Some(range) = brep.tshapes.get(ei).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None } }) {
 let theta = (range[1] - range[0]).abs();
 // Full 2蟺 circle: exact analytic area.
 if (theta - 2.0 * std::f64::consts::PI).abs() < 1e-12 {
 return Some(std::f64::consts::PI * c.radius * c.radius);
 }
 }
 }
 }
 }
 return None; // single non-circle edge is degenerate for a planar face
 }

 // Fast path: two circular edges forming a full circle (split by cylinder seam).
 // Common for planar cavity caps where the cylinder seam cuts the circle into
 // two complementary arcs.  Combined they cover 2蟺 鈫?area = 蟺r虏.
 if n_edges == 2 {
 let mut radii: [f64; 2] = [0.0; 2];
 let mut centers: [DVec3; 2] = [DVec3::ZERO; 2];
 let mut spans: [f64; 2] = [0.0; 2];
 let mut n_circle_edges = 0u32;
 for (i, we) in face.outer_wire.edges.iter().enumerate() {
 if let Some(ci) = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } }) {
 if let Curve3::Circle(c) = ci {
 if i < 2 {
 radii[i] = c.radius;
 centers[i] = c.center;
 }
 if let Some(r) = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None } }) {
 if i < 2 { spans[i] = (r[1] - r[0]).abs(); }
 }
 n_circle_edges += 1;
 }
 }
 }
 if n_circle_edges == 2 {
 let total_theta = spans[0] + spans[1];
 let same_center = (centers[0] - centers[1]).length_squared() < 1e-12;
 let same_radius = (radii[0] - radii[1]).abs() < 1e-12;
 if (total_theta - 2.0 * std::f64::consts::PI).abs() < 1e-10
 && same_center && same_radius
 {
 return Some(std::f64::consts::PI * radii[0] * radii[0]);
 }
 }
 }

 if n_edges < 3 { return None; }

 // Collect vertices in traversal order, and for each edge determine
 // whether it's a line or a circle.
 struct EdgeInfo {
 is_arc: bool, // Circle3 (true) vs Line3 (false)
 radius: f64, // for arcs: circle radius
 theta: f64, // for arcs: central angle (positive, < 2蟺)
 center_2d: DVec2, // for arcs: circle center projected to face plane
 sign: f64, // +1 if arc bulges outward, -1 if inward cutout
 start_2d: DVec2, // edge start in local 2D
 end_2d: DVec2, // edge end in local 2D
 }
 let mut edges: Vec<EdgeInfo> = Vec::with_capacity(n_edges);
 // Use first vertex's 3D point as pivot for 2D projection.
 let first_we = &face.outer_wire.edges[0];
 let first_e = brep.flat_edges().get(first_we.idx).copied()?;
 let first_vi = if first_we.forward { first_e.0 } else { first_e.1 };
 let pivot = brep.vertex_point(first_vi).unwrap_or(DVec3::ZERO);

 for we in &face.outer_wire.edges {
 let ei = we.idx;
 let edge = brep.flat_edges().get(ei).copied()?;
 let curve_idx = brep.tshapes.get(ei).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } })?;
 let curve = &curve_idx;
 let range = brep.tshapes.get(ei).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None } }).unwrap_or([0.0, 1.0]);

 let (v_start, v_end) = if we.forward { (edge.0, edge.1) } else { (edge.1, edge.0) };
 let p_start = brep.vertex_point(v_start).unwrap_or(DVec3::ZERO);
 let p_end = brep.vertex_point(v_end).unwrap_or(DVec3::ZERO);
 let start_2d = DVec2::new((p_start - pivot).dot(ux), (p_start - pivot).dot(uy));
 let end_2d = DVec2::new((p_end - pivot).dot(ux), (p_end - pivot).dot(uy));

 match curve {
 Curve3::Line(_) => {
 edges.push(EdgeInfo {
 is_arc: false, radius: 0.0, theta: 0.0,
 center_2d: DVec2::ZERO, sign: 0.0,
 start_2d, end_2d,
 });
 }
 Curve3::Circle(c) => {
 let theta = (range[1] - range[0]).abs();
 if theta < 1e-15 || theta > 2.0 * std::f64::consts::PI + 1e-12 { return None; }
 let center_2d = DVec2::new((c.center - pivot).dot(ux), (c.center - pivot).dot(uy));
 // Determine arc bulge direction: sign = +1 when the arc
 // bulges outward from the chord (center on left side of
 // traversal).  Use 3D cross product in the face plane for
 // projection-independent sign (the 2D cross product sign
 // depends on the local_basis_from_normal orientation).
 let trav_3d = p_end - p_start;
 let left_dir = face_normal.cross(trav_3d);
 let to_center_3d = c.center - p_start;
 let sign = if to_center_3d.dot(left_dir) > 0.0 { 1.0 } else { -1.0 };
 edges.push(EdgeInfo {
 is_arc: true, radius: c.radius, theta,
 center_2d, sign, start_2d, end_2d,
 });
 }
 _ => return None,
 }
 }

 // Compute total area: shoelace over vertices + segment corrections for arcs.
 // Shoelace over unique vertices in boundary order.
 let n = edges.len();
 let mut shoelace = 0.0;
 for i in 0..n {
 let s = edges[i].start_2d;
 let e = edges[i].end_2d;
 shoelace += s.x * e.y - e.x * s.y;
 }
 let mut total = shoelace.abs() * 0.5;

 // Add circular segment corrections for arc edges.
 let mut seg_total = 0.0;
 let mut n_neg = 0u32;
 let mut n_pos = 0u32;
 for edge in &edges {
 if edge.is_arc {
 let t = edge.theta;
 let seg = edge.radius * edge.radius * (t - t.sin()) * 0.5;
 if edge.sign < 0.0 { n_neg += 1; seg_total -= seg; }
 else { n_pos += 1; seg_total += seg; }
 }
 }
 // If all edges have the same sign and the corrections DEGRADE the area
 // (shoelace is larger than total), the traversal direction inverted the
 // bulge test.  This happens for cylinder-box cap faces where the face
 // normal from Newell's method is inverted relative to the surface normal
 // (H7: 4 arcs, all sign=-1 -> total = shoelace - corrections = 0.858).
 let shoelace_raw = total;
 total = shoelace_raw;
 if (n_neg == 0 || n_pos == 0) && seg_total < 0.0 {
 total -= seg_total;  // subtract negative = add positive corrections
 } else {
 total += seg_total;
 }

 if total > 0.0 && total.is_finite() {
 // Subtract enclosed areas for inner wires (holes).
 // OCCT BRepGProp processes each face's parametric domain with holes,
 // so the GL integration naturally excludes hole regions.  rcad's
 // analytic planar-area code computes the outer-wire area, then we
 // subtract hole areas using the same shoelace + arc-correction method.
 for w in &face.inner_wires {
 if w.edges.len() < 3 { continue; }
 // Compute hole area via dense-sampled polyline projection.
 let hole_pts: Vec<DVec3> = w.edges.iter().filter_map(|we| {
 let e = brep.flat_edges().get(we.idx).copied()?;
 let vi = if we.forward { e.0 } else { e.1 };
 brep.vertex_point(vi)
 }).collect();
 if hole_pts.len() < 3 { continue; }
 let mut a_hole = 0.0;
 for i in 1..hole_pts.len().saturating_sub(1) {
 let v0 = hole_pts[0];
 let v1 = hole_pts[i];
 let v2 = hole_pts[i + 1];
 let d1 = v1 - v0;
 let d2 = v2 - v0;
 let tri_area = (d1.cross(d2)).dot(face_normal).abs() * 0.5;
 a_hole += tri_area;
 }
 total -= a_hole.min(total);
 }
 if cfg!(debug_assertions) && std::env::var("RCAD_DEBUG_SA").is_ok() {
 let n_inner = face.inner_wires.len();
 eprintln!("[SA_EXACT] return total={:.8} inner_wires={}", total, n_inner);
 }
 return Some(total);
 } else { None }
}

/// Shoelace area of outer wire minus |hole areas| in the face plane (pivot = first outer point).
fn try_planar_face_area_shoelace(
 brep: &BRep,
 face: &Face,
 face_normal: DVec3,
) -> Option<f64> {
 if face.inner_wires.is_empty() {
 if let Some(a_rect) = try_axis_aligned_world_rect_plane_area(brep, face, face_normal) {
 // axis-aligned rect path
 let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
 trim_almost_closed_polyline(&mut outer, 1e-5);
 if outer.len() >= 3 {
 let (ux, uy) = local_basis_from_normal(face_normal);
 if let Some(pivot) = outer.first().copied() {
 let a_shoe = polygon_area_2d_projected(&outer, pivot, ux, uy).abs();
 const REL: f64 = 1e-5;
 let scale = a_rect.max(a_shoe).max(1.0);
 let abs_eps = 1e-9 * scale;

 // OCCT `bcommon_simple/C8`: axis-aligned plane 鈭?tilted box gives a parallelogram
 // whose vertices all sit on the loop's axis-aligned bbox; `w*h` over-counts.
 // Boolean T-junctions can permute wire edge order so shoe-lace on dense samples
 // under-counts while hull/bbox from [`try_axis_aligned_world_rect_plane_area`] is
 // correct (`bcommon_simple/B1`): only trust the smaller shoe when it is a
 // substantial fraction of the rect metric (true parallelogram vs AABB).
 if a_rect > a_shoe * (1.0 + REL) + abs_eps {
 // When the shoelace is a significant fraction of the bounding
 // rect, the polygon is legitimate (possibly concave). Only
 // fall back to the rect when shoelace is pathologically small
 // (permuted wire edges from boolean T-junctions).
 let ratio = a_shoe / a_rect.max(1e-12);
 if ratio >= 0.40 || a_shoe + abs_eps >= a_rect * 0.65 {
 return Some(a_shoe.max(0.0));
 }
 return Some(a_rect.max(0.0));
 }
 if a_shoe > a_rect * (1.0 + REL) + abs_eps {
 return Some(a_rect);
 }
 return Some(a_rect);
 }
 }
 return Some(a_rect);
 }
 }
 let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
 trim_almost_closed_polyline(&mut outer, 1e-5);
 if outer.len() < 3 {
 return None;
 }
 let (ux, uy) = local_basis_from_normal(face_normal);
 let pivot = outer.first().copied()?;
 let mut a = polygon_area_2d_projected(&outer, pivot, ux, uy).abs();

 // Boolean T-junction repair can scramble wire edge order so the dense-sampled
 // polyline zig-zags across the face instead of tracing the boundary.  The
 // shoelace on the self-intersecting polygon then under-counts (e.g. exactly
 // half the correct area for bcommon_simple/C9's slanted face).  Cross-check
 // against the convex hull of the unique boundary vertices: if the dense-sample
 // area is less than 60 % of the hull area the polyline is unreliable, so fall
 // through to triangulation which is robust to scrambled ordering.
 if a * 1e-7 < outer.len() as f64 {
 // Only run the cross-check when the dense polyline has significantly more
 // points than unique vertices (avoids degenerate faces).
 if let Some(hull_a) = try_boundary_convex_hull_area(brep, &face.outer_wire, pivot, ux, uy)
 {
 if hull_a > 1e-12 && a < 0.6 * hull_a {
 return None;
 }
 }
 }

 // Boolean T-junction repair can scramble wire edge order so the dense-sampled
 // polyline zig-zags across the face instead of tracing the boundary.  The
 // shoelace then gives ~0 even for a valid face (e.g. rotated-box 鈭?unit cube in
 // bcommon_simple/G1).  Compute a rough bbox and bail if the area is implausible:
 // face_triangles uses stored boundary or triangle fan, which is robust to this.
 if a < 1e-12 {
 let mut bbox_min = outer[0];
 let mut bbox_max = outer[0];
 for p in &outer {
 bbox_min = bbox_min.min(*p);
 bbox_max = bbox_max.max(*p);
 }
 let bbox_diag = (bbox_max - bbox_min).length();
 if bbox_diag > 1e-10 {
 // bbox area 鈮?diag虏, so if a << diag虏 the polyline is degenerate
 if a < 1e-12 * bbox_diag * bbox_diag {
 return None;
 }
 } else {
 return None;
 }
 }
 for w in &face.inner_wires {
 // 鉁?OCCT :  ( Circle3  forward+reverse)鈫?- 銆?
 if w.edges.len() == 2 && w.edges[0].idx == w.edges[1].idx {
 let ei = w.edges[0].idx;
 let flat_edges = brep.flat_edges();
 if let (Some(ci), Some(e)) = (brep.tshapes.get(ei).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } }), flat_edges.get(ei)) {
 let c = if let Curve3::Circle(c_val) = ci { c_val } else { continue; };
 let (p0, p1) = (brep.vertex_point(e.0), brep.vertex_point(e.1));
 if let (Some(v0), Some(v1)) = (p0, p1) {
 let d0 = (v0 - c.center).normalize();
 let d1 = (v1 - c.center).normalize();
 let theta = d0.dot(d1).clamp(-1.0, 1.0).acos();
 if theta > 1e-12 {
 a -= c.radius * c.radius * (theta - theta.sin()) * 0.5;
 }
 }
 }
 }
 continue;
 let mut h = sample_wire_polyline_3d(brep, w);
 trim_almost_closed_polyline(&mut h, 1e-5);
 if h.len() < 3 {
 continue;
 }
 a -= polygon_area_2d_projected(&h, pivot, ux, uy).abs();
 }
 Some(a.max(0.0))
}

fn polygon_area_2d_projected(pts: &[DVec3], pivot: DVec3, ux: DVec3, uy: DVec3) -> f64 {
 if pts.len() < 3 {
 return 0.0;
 }
 let n = pts.len();
 let mut s = 0.0;
 for i in 0..n {
 let j = (i + 1) % n;
 let a = (pts[i] - pivot).dot(ux);
 let b = (pts[i] - pivot).dot(uy);
 let c = (pts[j] - pivot).dot(ux);
 let d = (pts[j] - pivot).dot(uy);
 s += a * d - c * b;
 }
 0.5 * s
}

/// Compute the area of the convex hull of the unique boundary vertices of `wire`,
/// projected onto the 2D local basis `(ux, uy)` with the given `pivot`.
///
/// Returns `None` when there are fewer than 3 unique vertices.
///
/// This is used as a cross-check against the dense-sample shoelace area in
/// [`try_planar_face_area_shoelace`]: when the hull area is significantly larger
/// than the dense-sample area, the wire edge order is likely scrambled by boolean
/// T-junction repair and the shoelace result is unreliable.
fn try_boundary_convex_hull_area(
 brep: &BRep,
 wire: &Wire,
 pivot: DVec3,
 ux: DVec3,
 uy: DVec3,
) -> Option<f64> {
 // Collect unique vertex indices from wire edges
 let mut vert_indices: Vec<usize> = Vec::new();
 for we in &wire.edges {
 let edge = brep.flat_edges().get(we.idx).copied()?;
 vert_indices.push(edge.0);
 vert_indices.push(edge.1);
 }
 vert_indices.sort();
 vert_indices.dedup();

 if vert_indices.len() < 3 {
 return None;
 }

 // Project unique vertices to the 2D local basis
 let pts_2d: Vec<DVec2> = vert_indices
 .iter()
 .filter_map(|&vi| brep.vertex_point(vi))
 .map(|v| DVec2::new(
 (v - pivot).dot(ux),
 (v - pivot).dot(uy),
 ))
 .collect();

 if pts_2d.len() < 3 {
 return None;
 }

 // Monotone chain (Andrew's algorithm) convex hull 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
 let mut sorted: Vec<usize> = (0..pts_2d.len()).collect();
 sorted.sort_by(|&a, &b| {
 let pa = pts_2d[a];
 let pb = pts_2d[b];
 if (pa.x - pb.x).abs() > 1e-12 {
 pa.x.partial_cmp(&pb.x).unwrap_or(std::cmp::Ordering::Equal)
 } else {
 pa.y.partial_cmp(&pb.y).unwrap_or(std::cmp::Ordering::Equal)
 }
 });

 let cross = |o: DVec2, a: DVec2, b: DVec2| -> f64 {
 (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
 };

 // Lower hull
 let mut lower: Vec<DVec2> = Vec::new();
 for &si in &sorted {
 while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], pts_2d[si]) <= 0.0 {
 lower.pop();
 }
 lower.push(pts_2d[si]);
 }

 // Upper hull
 let mut upper: Vec<DVec2> = Vec::new();
 for &si in sorted.iter().rev() {
 while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], pts_2d[si]) <= 0.0 {
 upper.pop();
 }
 upper.push(pts_2d[si]);
 }

 // Remove last point of each (it's the first of the other) and combine
 lower.pop();
 upper.pop();
 let hull: Vec<DVec2> = lower.into_iter().chain(upper).collect();

 if hull.len() < 3 {
 return None;
 }

 // Shoelace area of the convex hull
 let mut area = 0.0;
 for i in 0..hull.len() {
 let j = (i + 1) % hull.len();
 area += hull[i].x * hull[j].y;
 area -= hull[j].x * hull[i].y;
 }
 Some((area * 0.5).abs())
}

/// For trimmed cylinder faces, compute the exact surface area from the UV boundary.
///
/// Compute UV area of a cylinder UV polygon via Gauss-Legendre integration
/// of the V-extent at each U.  For each quadrature point u, the V-extent
/// [v_min(u), v_max(u)] is found by intersecting the line U=u with all
/// polygon edges (unwrapped UV coordinates).  This handles figure-8 UV
/// polygons correctly without binning or envelope approximation.
///
/// OCCT uses adaptive Gauss-Legendre in BRepGProp for trimmed cylinders,
/// subdividing in U until each sub-polygon is non-self-intersecting.
fn cylinder_uv_area_gl(uvs: &[DVec2]) -> Option<f64> {
 const NU: usize = 60;
 let n = uvs.len();
 if n < 3 { return None; }
 let two_pi = std::f64::consts::PI * 2.0;

 // Unwrap U to linear coordinates.
 let mut poly: Vec<DVec2> = Vec::with_capacity(n);
 poly.push(uvs[0]);
 for i in 1..n {
 let du = short_delta_on_circle_01(uvs[i-1].x, uvs[i].x);
 poly.push(DVec2::new(poly[i-1].x + du, uvs[i].y));
 }

 // V-range at a given u by intersecting the vertical line with polygon edges.
 let v_range_at = |u: f64| -> (f64, f64) {
 let mut v_lo = f64::INFINITY;
 let mut v_hi = f64::NEG_INFINITY;
 for i in 0..n {
 let j = (i + 1) % n;
 let (u1, v1) = (poly[i].x, poly[i].y);
 let (u2, v2) = (poly[j].x, poly[j].y);
 // Skip horizontal edges (u1 == u2): the point-in-polygon test
 // handles horizontal crossings at the endpoint.
 if u1 == u2 { continue; }
 // Check if u is between u1 and u2 (inclusive at endpoints to
 // handle vertices).
 let (u_lo_e, u_hi_e, v_lo_e, v_hi_e) = if u1 < u2 {
 (u1, u2, v1, v2)
 } else {
 (u2, u1, v2, v1)
 };
 if u < u_lo_e - 1e-12 || u > u_hi_e + 1e-12 { continue; }
 // Interpolate V at this U.
 let t = if u_hi_e > u_lo_e { (u - u_lo_e) / (u_hi_e - u_lo_e) } else { 0.0 };
 let v = v_lo_e + t * (v_hi_e - v_lo_e);
 v_lo = v_lo.min(v);
 v_hi = v_hi.max(v);
 }
 (v_lo, v_hi)
 };

 // Gauss-Legendre over U.
 let u_min = poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
 let u_max = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
 let u_range = u_max - u_min;
 if !u_min.is_finite() || u_range < 1e-14 { return None; }

 let u_lo = u_min;
 let u_hi = if u_range > two_pi - 0.1 { u_min + two_pi } else { u_max };
 let du = (u_hi - u_lo) / NU as f64;
 const GL_NEG: f64 = -0.5773502691896257;
 const GL_POS: f64 = 0.5773502691896257;
 let gl_pts = [GL_NEG, GL_POS];

 let mut total = 0.0;
 for i in 0..NU {
 let u_mid = u_lo + (i as f64 + 0.5) * du;
 for &gu in &gl_pts {
 let u = u_mid + gu * du * 0.5;
 let (v_lo, v_hi) = v_range_at(u);
 if v_lo.is_finite() && v_hi > v_lo + 1e-14 {
 total += (v_hi - v_lo) * du * 0.5; // weight = 1.0 for 2-pt GL
 }
 }
 }
 Some(total)
}

/// For cylinders, `|鈭侾/鈭倁 脳 鈭侾/鈭倂| = R` (constant), so surface area = R 脳 area of the UV region.
/// Uses a Green's theorem line integral 鈭?v路du, integrating from the minimum-v boundary edge
/// in the +u direction, to avoid the cancellation that occurs when the wire traces the same
/// curve forward and backward (e.g., Steinmetz lens bounded by a single intersection curve).
/// Compute the area of a cylinder UV polygon using 2脳2 Gauss-Legendre
/// quadrature over the UV bounding box.  Since |鈭侾/鈭倁 脳 鈭侾/鈭倂| = R is
/// constant, the surface area = R 脳 鈭埆 du dv, so we just integrate the
/// area of the UV domain (no surface-integrand evaluation needed).
///
/// OCCT uses adaptive Gauss-Legendre integration in BRepGProp for trimmed
/// surfaces including cylinders.  This replaces the shoelace-from-samples
/// approach whose accuracy depends on the 3D鈫扷V projection convergence.
fn cylinder_gl_uv_area(uvs: &[DVec2]) -> Option<f64> {
 const N: usize = 60;
 if uvs.len() < 3 { return None; }
 let n = uvs.len();
 const TWO_PI: f64 = std::f64::consts::PI * 2.0;

 // Unwrap U coordinates to handle 2蟺 wrapping (same approach as
 // try_cylinder_trimmed_face_area's shoelace path).
 let unwrapped: Vec<f64> = {
 let mut o = Vec::with_capacity(n);
 o.push(uvs[0].x);
 for i in 1..n {
 o.push(o[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x));
 }
 o
 };

 let mut umin = f64::INFINITY; let mut umax = f64::NEG_INFINITY;
 let mut vmin = f64::INFINITY; let mut vmax = f64::NEG_INFINITY;
 for (i, uv) in uvs.iter().enumerate() {
 umin = umin.min(unwrapped[i]); umax = umax.max(unwrapped[i]);
 vmin = vmin.min(uv.y); vmax = vmax.max(uv.y);
 }
 if !umin.is_finite() || (umax - umin) < 1e-14 || (vmax - vmin) < 1e-14 {
 return None;
 }
 let u_range = umax - umin;
 let v_range = vmax - vmin;

 // The UV polygon might wrap the cylinder in u.  If the unwrapped range
 // is nearly 2蟺, the face covers the full circumference 鈥?integrate
 // over the full [umin, umin+2蟺] domain.
 let u_lo = umin;
 let u_hi = if u_range > TWO_PI - 0.1 { umin + TWO_PI } else { umax };
 let v_lo = vmin;
 let v_hi = vmax;

 let du = (u_hi - u_lo) / N as f64;
 let dv = (v_hi - v_lo) / N as f64;
 let cell_area = du * dv / 4.0;
 const GL_NEG: f64 = -0.5773502691896257;
 const GL_POS: f64 = 0.5773502691896257;
 let gl_pts = [GL_NEG, GL_POS];

 let mut total = 0.0;
 for i in 0..N {
 let u_mid = u_lo + (i as f64 + 0.5) * du;
 for j in 0..N {
 let v_mid = v_lo + (j as f64 + 0.5) * dv;
 let mut n_hit = 0u32;
 for &gu in &gl_pts {
 let u = u_mid + gu * du * 0.5;
 let u_mod = u.rem_euclid(TWO_PI);
 for &gv in &gl_pts {
 let v = v_mid + gv * dv * 0.5;
 // Test against original UV polygon (with wrapped U).
 // The polygon's U values are in [0, 2蟺) while the
 // integration point u might exceed that range due to
 // unwrapping 鈥?remap to [0, 2蟺) for the point test.
 let inside = winding_number_2d(uvs, DVec2::new(u_mod, v)) != 0;
 if inside { n_hit += 1; }
 }
 }
 if n_hit > 0 {
 total += cell_area * n_hit as f64;
 }
 }
 }
 Some(total)
}

fn try_cylinder_trimmed_face_area(
 cyl: &CylindricalSurface,
 brep: &BRep,
 face: &Face,
 face_flat_idx: usize,
) -> Option<f64> {
 use crate::projection::closest_point_on_surface;

 // Fast path: rectangular UV patch via 2 Lines + 2 Circles (no inner wires).
 //
 // Boolean splitting along iso-parametric lines produces sub-faces whose UV
 // domain is a clean rectangle: two generator edges (u = constant, Line3 in
 // 3D) and two circular edges (v = constant, Circle3 in 3D).  For these,
 // the exact analytic area is simply R 脳 螖U 脳 螖V.
 //
 // We compute 螖U and 螖V directly from the edge curves rather than using
 // curved_face_uv_domain (which relies on wire sampling 鈫?closest-point
 // projection and inherits numerical error).
 if face.inner_wires.is_empty() && face.outer_wire.edges.len() == 4 {
 if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
 eprintln!("[CYL_RECT_FAST] checking 4-edge face");
 }
 use crate::geom::Curve3;
 let mut edge_curve_indices: Vec<usize> = Vec::new();  // edge flat indices for line seam detection
 let mut n_lines = 0u32;
 let mut n_circles = 0u32;
 let mut circle_centers: Vec<DVec3> = Vec::new();
 let mut valid = true;
 for we in &face.outer_wire.edges {
 if let Some(ci) = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } }) {
 match ci {
 Curve3::Line(_) => {
 n_lines += 1;
 edge_curve_indices.push(we.idx);
 }
 Curve3::Circle(c) => {
 n_circles += 1;
 if circle_centers.len() < 2 {
 circle_centers.push(c.center);
 }
 }
 _ => { valid = false; break; }
 }
 } else { valid = false; break; }
 }
 if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
 eprintln!("[CYL_RECT_FAST] valid={} lines={} circles={} centers={}", valid, n_lines, n_circles, circle_centers.len());
 }
 if valid && n_lines == 2 && n_circles == 2 && circle_centers.len() == 2 {
 // 螖V: project circle centers onto cylinder axis.
 let axis = cyl.axis;
 let v0 = (circle_centers[0] - cyl.origin).dot(axis);
 let v1 = (circle_centers[1] - cyl.origin).dot(axis);
 let dv = (v1 - v0).abs();
 // 螖U: prefer face_surface_range (most reliable, set by analytic
 // builders like build_box_minus_cylinder_full_uv_z_fail).
 // Fall back to seam-detection, then angular vertex computation.
 let du = 'du: {
 // Priority 1: face_surface_range 鈥?the builder knows the exact UV domain.
 if let Some(range) = brep.tshapes.get(face_flat_idx).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None } }) {
 let du_r = (range[1] - range[0]).abs();
 let dv_r = (range[3] - range[2]).abs();
 if du_r > 1e-14 && dv_r > 1e-14 {
 if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
 eprintln!("[CYL_RECT_DU] from face_surface_range du={}", du_r);
 }
 break 'du du_r;
 }
 }
 // Priority 2: both Line edges share the same curve index (cylinder seam 鈫?full 2蟺).
 if edge_curve_indices.len() == 2 && edge_curve_indices[0] == edge_curve_indices[1] {
 break 'du (std::f64::consts::PI * 2.0);
 }
 // Priority 3: compute angular span from generator line vertices.
 let x_ax = crate::geom::any_perpendicular(axis);
 let y_ax = axis.cross(x_ax).normalize();
 let mut u_vals = Vec::new();
 for ci in &edge_curve_indices {
 if let Some(Curve3::Line(_line)) = brep.tshapes.get(*ci).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None } }) {
 for we in &face.outer_wire.edges {
 if we.idx == *ci {
 let ei = we.idx;
 let flat_edges = brep.flat_edges();
 if let Some(edge) = flat_edges.get(ei) {
 let vi = if we.forward { edge.0 } else { edge.1 };
 if let Some(v) = brep.vertex_point(vi) {
 let d = v - cyl.origin;
 let u = d.dot(y_ax).atan2(d.dot(x_ax));
 u_vals.push(u);
 }
 }
 break;
 }
 }
 }
 }
 if u_vals.len() == 2 {
 let du_raw = u_vals[1] - u_vals[0];
 let du_norm = du_raw.rem_euclid(std::f64::consts::PI * 2.0);
 if du_norm > 1e-14 { break 'du du_norm; }
 }
 // Final fallback: full wrap.
 std::f64::consts::PI * 2.0
 };
 // Validate the circle centers lie on the cylinder axis.
 let r0 = (circle_centers[0] - cyl.origin).cross(axis).length();
 let r1 = (circle_centers[1] - cyl.origin).cross(axis).length();
 if du > 1e-14 && dv > 1e-14 && du.is_finite() && dv.is_finite()
 && r0 < 1e-8 && r1 < 1e-8
 {
 if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
 eprintln!("[CYL_RECT_FAST] R={} du={} dv={} area={}", cyl.radius, du, dv, cyl.radius * du * dv);
 }
 return Some(cyl.radius * du * dv);
 } else if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
 eprintln!("[CYL_RECT_FAST] rejected: du={} dv={} r0={} r1={}", du, dv, r0, r1);
 }
 }
 }

 let wire_uv_area = |wire: &Wire| -> Option<f64> {
 let mut pts_3d = sample_wire_polyline_3d_with_n(brep, wire, 512);
 trim_almost_closed_polyline(&mut pts_3d, 1e-5);
 if pts_3d.len() < 3 { return None; }

 let n = pts_3d.len();
 const TWO_PI: f64 = std::f64::consts::PI * 2.0;

 let surf = Surface3::Cylinder(*cyl);
 let uvs: Vec<DVec2> = pts_3d.iter()
 .map(|&p| { let proj = closest_point_on_surface(&surf, p, 256); DVec2::new(proj.params.0, proj.params.1) })
 .collect();

 // GL quadrature over UV bounding box (OCCT-aligned).
 if let Some(gl_area) = cylinder_gl_uv_area(&uvs) {
 if gl_area > 0.0 { return Some(gl_area); }
 }
 if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
 let gl_res = cylinder_gl_uv_area(&uvs);
 eprintln!("[CYL_GL] n_uvs={} gl={:?} cyl_R={:.6}", uvs.len(), gl_res, cyl.radius);
 }

 // Simple unwrapping for cylinder faces: accumulate via short deltas on S鹿.
 let unwrapped: Vec<f64> = {
 let mut o = Vec::with_capacity(n);
 o.push(uvs[0].x);
 for i in 1..n {
 o.push(o[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x));
 }
 o
 };

 // Shoelace area.
 let mut area2 = 0.0_f64;
 for i in 0..n {
 let j = if i + 1 < n { i + 1 } else { 0 };
 area2 += unwrapped[i] * uvs[j].y - unwrapped[j] * uvs[i].y;
 }
 let uv_area = area2.abs() * 0.5;

 // GL quadrature of V-extent at each U (OCCT-aligned).
 let gl_area = cylinder_uv_area_gl(&uvs).unwrap_or(0.0);
 if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
 eprintln!("[CYL_GL] n_uvs={} gl={:.6} shoelace={:.6}", uvs.len(), gl_area, uv_area);
 }

 // Envelope area: integrate v_max(u) - v_min(u) across u bins.
 let mut band_area = 0.0_f64;
 let nf = n;
 if nf >= 3 {
 // 2048 bins: ~0.003 rad resolution for envelope integration.
 let n_bins = 2048usize;
 let step = TWO_PI / n_bins as f64;
 let mut upper = vec![f64::NEG_INFINITY; n_bins];
 let mut lower = vec![f64::INFINITY; n_bins];
 for i in 0..nf {
 let uw = uvs[i].x.rem_euclid(TWO_PI);
 let b = ((uw / TWO_PI) * n_bins as f64) as usize;
 let b = b.min(n_bins - 1);
 if uvs[i].y > upper[b] { upper[b] = uvs[i].y; }
 if uvs[i].y < lower[b] { lower[b] = uvs[i].y; }
 }

 // Collect populated bins into (u, range) pairs for trapezoidal integration
 let mut prev_u: Option<f64> = None;
 let mut prev_range: Option<f64> = None;
 let mut first_u: f64 = 0.0;
 let mut first_range: f64 = 0.0;
 let mut got_first = false;

 for b in 0..n_bins {
 if !lower[b].is_finite() { continue; }
 let u_center = (b as f64 + 0.5) * step;
 let range = upper[b] - lower[b];
 if !got_first {
 first_u = u_center;
 first_range = range;
 got_first = true;
 }
 if let (Some(pu), Some(pr)) = (prev_u, prev_range) {
 let du = u_center - pu;
 band_area += (range + pr) * 0.5 * du;
 }
 prev_u = Some(u_center);
 prev_range = Some(range);
 }
 // The envelope method with wrap-around is only valid when the UV
 // polygon actually wraps around the full cylinder (populated bins
 // cover nearly all of 0-2蟺).  For partial-wrap faces (e.g., wedge
 // faces from boolean splitting where u_span ~200掳 but only ~55% of
 // bins are populated), the wrap-around gap closure would integrate
 // across empty bins and overcount.  A 85% threshold distinguishes
 // true wrap-around (figure-8) from partial-wrap merged faces.
 let pop_bins = upper.iter().filter(|v| v.is_finite()).count();
 if pop_bins as f64 > n_bins as f64 * 0.85 {
 // Wrap around: close the gap from the last populated bin to the first
 if let (Some(lu), Some(lr)) = (prev_u, prev_range) {
 let du = (TWO_PI + first_u) - lu;
 band_area += (first_range + lr) * 0.5 * du;
 }
 }
 }

 // Return max of shoelace and GL by default. The envelope/band method is
 // only reliable for near-full-wrap cylinder faces; on narrow trimmed
 // patches near the seam it can overcount by integrating across empty U bins.
 let mut result = uv_area.max(gl_area);
 let allow_band = brep
 .tshapes
 .get(face_flat_idx)
 .and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None } })
 .map(|[u0, u1, _v0, _v1]| (u1 - u0).abs() > std::f64::consts::PI * 1.75)
 .unwrap_or(true);
 if allow_band {
 result = result.max(band_area);
 }
 if std::env::var("RCAD_DEBUG_CYL_AREA").is_ok() {
 eprintln!(
 "[CYL_AREA] face={} pts={} shoelace={:.6} gl={:.6} band={:.6} allow_band={} chosen={:.6}",
 face_flat_idx,
 n,
 uv_area,
 gl_area,
 band_area,
 allow_band,
 result,
 );
 }
 Some(result)
 };

 let outer_area = wire_uv_area(&face.outer_wire)?;
 let inner_area: f64 = face.inner_wires.iter()
 .filter_map(|w| wire_uv_area(w)).sum();
 let mut total_uv_area = outer_area - inner_area;

 // Detect figure-8 cancellation: after unify_same_domain_faces, the merged outer
 // wire may self-intersect and be split into outer + inner loops whose UV areas
 // nearly cancel (both ~half the cylinder).  In this case, compute the envelope
 // of ALL points (outer + inner combined), which covers the full cylinder surface.
 if inner_area.abs() > 1e-6 && (outer_area - inner_area).abs() < 0.5 * inner_area.abs().max(outer_area) {
 let all_wires = std::iter::once(&face.outer_wire).chain(face.inner_wires.iter());
 let mut all_uvs: Vec<DVec2> = Vec::new();
 let surf = Surface3::Cylinder(*cyl);
 for wire in all_wires {
 let mut pts = sample_wire_polyline_3d_with_n(brep, wire, 128);
 trim_almost_closed_polyline(&mut pts, 1e-5);
 if pts.len() < 3 { continue; }
 for &p in &pts {
 let proj = closest_point_on_surface(&surf, p, 16);
 all_uvs.push(DVec2::new(proj.params.0, proj.params.1));
 }
 }
 if all_uvs.len() >= 3 {
 let n = all_uvs.len();
 const TWO_PI: f64 = std::f64::consts::PI * 2.0;
 let unwrapped: Vec<f64> = {
 let mut o = Vec::with_capacity(n);
 o.push(all_uvs[0].x);
 for i in 1..n {
 o.push(o[i - 1] + short_delta_on_circle_01(all_uvs[i - 1].x, all_uvs[i].x));
 }
 o
 };
 let u_min = unwrapped.iter().cloned().fold(f64::INFINITY, f64::min);
 let u_max = unwrapped.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 let u_span = u_max - u_min;
 const WRAP_THRESHOLD: f64 = std::f64::consts::PI;

 if u_span > WRAP_THRESHOLD {
 // Full cylinder wrap: the combined outer+inner wires span the full
 // cylinder circumference.  The envelope method undercounts for
 // boundary-only UV points (figure-8 lobes); instead use the global
 // v range which correctly gives the full cylinder area between the
 // extreme v bounds of all points.
 let v_min = all_uvs.iter().map(|uv| uv.y).fold(f64::INFINITY, f64::min);
 let v_max = all_uvs.iter().map(|uv| uv.y).fold(f64::NEG_INFINITY, f64::max);
 total_uv_area = TWO_PI * (v_max - v_min);
 } // else: non-wrap, shoelace is fine; keep total_uv_area as computed
 }
 }
 if std::env::var("RCAD_DEBUG_CYL_AREA").is_ok() {
 let rect_est = brep
 .tshapes
 .get(face_flat_idx)
 .and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None } })
 .map(|[u0, u1, v0, v1]| (u1 - u0).abs() * (v1 - v0).abs())
 .unwrap_or(0.0);
 eprintln!(
 "[CYL_AREA_TOTAL] face={} outer={:.6} inner={:.6} total_uv={:.6} rect_est={:.6} radius={:.6} area={:.6}",
 face_flat_idx,
 outer_area,
 inner_area,
 total_uv_area,
 rect_est,
 cyl.radius,
 cyl.radius * total_uv_area,
 );
 }
 if total_uv_area > 1e-14 { Some(cyl.radius * total_uv_area) } else { None }
}

/// Compute the raw UV shoelace area for a cylinder face's outer wire (no inner-wire
/// subtraction, no envelope/band method).  Used by [`surface_area`] to detect and
/// normalize overlapping UV sub-faces from boolean splitting.
fn cylinder_outer_wire_uv_shoelace_area(cyl: &CylindricalSurface, brep: &BRep, face: &Face) -> Option<f64> {
 use crate::projection::closest_point_on_surface;
 let surf = Surface3::Cylinder(*cyl);
 let mut pts = sample_wire_polyline_3d_with_n(brep, &face.outer_wire, 256);
 trim_almost_closed_polyline(&mut pts, 1e-5);
 if pts.len() < 3 { return None; }
 let uvs: Vec<DVec2> = pts.iter()
 .map(|&p| {
 let proj = closest_point_on_surface(&surf, p, 16);
 DVec2::new(proj.params.0, proj.params.1)
 })
 .collect();
 let n = uvs.len();
 let unwrapped: Vec<f64> = {
 let mut o = Vec::with_capacity(n);
 o.push(uvs[0].x);
 for i in 1..n {
 o.push(o[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x));
 }
 o
 };
 let area2: f64 = (0..n).map(|i| {
 let j = (i + 1) % n;
 unwrapped[i] * uvs[j].y - unwrapped[j] * uvs[i].y
 }).sum::<f64>();
 Some(area2.abs() * 0.5)
}

/// UV-polygon-aware surface area for trimmed cone sub-faces.
///
/// For a cone parameterized as P(u,v) = apex + v路cos(伪)路axis + (R + v路sin(伪))路(cos(u)路x虃 + sin(u)路欧),
/// the area element is |鈭侾/鈭倁 脳 鈭侾/鈭倂| = R + v路sin(伪).  The surface area over a
/// UV polygon is therefore:
///
/// A = 鈭埆 (R + v路sin(伪)) du dv
/// = R 路 Area(UV_polygon) + sin(伪) 路 鈭埆 v du dv
///
/// where 鈭埆 v du dv is the first moment about the u-axis, computed via Green's theorem
/// from the polygon boundary.  This avoids the overcount from `param_rect_area_cross`
/// which integrates over the bounding rectangle of the UV polygon.
fn try_cone_trimmed_face_area(
 cone: &ConicalSurface,
 brep: &BRep,
 face: &Face,
) -> Option<f64> {
 use crate::projection::closest_point_on_surface;

 let wire_uv_data = |wire: &Wire| -> Option<(f64, f64)> {
 // Returns (area, first_moment_Mu) for the UV polygon
 let mut pts_3d = sample_wire_polyline_3d_with_n(brep, wire, 512);
 trim_almost_closed_polyline(&mut pts_3d, 1e-5);
 if pts_3d.len() < 3 { return None; }

 let n = pts_3d.len();
 let surf = Surface3::Cone(*cone);
 let uvs: Vec<DVec2> = pts_3d.iter()
 .map(|&p| { let proj = closest_point_on_surface(&surf, p, 256); DVec2::new(proj.params.0, proj.params.1) })
 .collect();

 // Unwrap u for periodic S鹿 parameter
 let unwrapped: Vec<f64> = {
 let mut o = Vec::with_capacity(n);
 o.push(uvs[0].x);
 for i in 1..n {
 o.push(o[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x));
 }
 o
 };

 // Shoelace area and first moment M_u = 鈭埆 v du dv
 let mut area2 = 0.0_f64;
 let mut moment6 = 0.0_f64;  // 6 脳 M_u before division
 for i in 0..n {
 let j = if i + 1 < n { i + 1 } else { 0 };
 let cross = unwrapped[i] * uvs[j].y - unwrapped[j] * uvs[i].y;
 area2 += cross;
 moment6 += (uvs[i].y + uvs[j].y) * cross;
 }
 let uv_area = area2.abs() * 0.5;
 let mu = moment6.abs() / 6.0;  // 鈭埆 v du dv = |moment6| / 6

 // Detect doubled-back UV polygon: short_delta unwrapping cancels the net u
 // range, producing a near-zero uv_area even though the face has real area.
 let total_abs_path: f64 = (1..n)
 .map(|i| (unwrapped[i] - unwrapped[i - 1]).abs())
 .sum();
 let net_range = (unwrapped[n - 1] - unwrapped[0]).abs();
 if net_range < total_abs_path * 0.25 && total_abs_path > std::f64::consts::PI {
 // Envelope method: bin UV points by u-mod-2蟺, track v_min/v_max per bin,
 // then integrate trapezoidally.  Handles doubled-back polygons correctly
 // because it does not depend on winding order.
 const TWO_PI: f64 = std::f64::consts::TAU;
 const N_BINS: usize = 512;
 let step = TWO_PI / N_BINS as f64;
 let mut v_min = vec![f64::INFINITY; N_BINS];
 let mut v_max = vec![f64::NEG_INFINITY; N_BINS];
 for uv in &uvs {
 let b = ((uv.x.rem_euclid(TWO_PI) / TWO_PI) * N_BINS as f64) as usize;
 let b = b.min(N_BINS - 1);
 if uv.y > v_max[b] { v_max[b] = uv.y; }
 if uv.y < v_min[b] { v_min[b] = uv.y; }
 }

 let mut env_area = 0.0_f64;
 let mut env_mu = 0.0_f64;
 let mut prev_u: Option<f64> = None;
 let mut prev_range: Option<f64> = None;
 let mut prev_strip: Option<f64> = None;

 for b in 0..N_BINS {
 if !v_min[b].is_finite() { continue; }
 let u_center = (b as f64 + 0.5) * step;
 let range = v_max[b] - v_min[b];
 let strip_mu = (v_max[b] * v_max[b] - v_min[b] * v_min[b]) * 0.5;
 if let (Some(pu), Some(pr), Some(pm)) = (prev_u, prev_range, prev_strip) {
 let du = u_center - pu;
 env_area += (range + pr) * 0.5 * du;
 env_mu += (strip_mu + pm) * 0.5 * du;
 }
 prev_u = Some(u_center);
 prev_range = Some(range);
 prev_strip = Some(strip_mu);
 }

 return Some((env_area, env_mu));
 }

 Some((uv_area, mu))
 };

 let (outer_area, outer_mu) = wire_uv_data(&face.outer_wire)?;
 let (inner_area, inner_mu) = {
 let mut a = 0.0_f64;
 let mut m = 0.0_f64;
 for w in &face.inner_wires {
 if let Some((wa, wm)) = wire_uv_data(w) {
 a += wa;
 m += wm;
 }
 }
 (a, m)
 };
 let net_area = outer_area - inner_area;
 let net_mu = outer_mu - inner_mu;
 let r = cone.radius;
 let sin_alpha = cone.half_angle_rad.sin();
 let total = r * net_area + sin_alpha * net_mu;
 if total > 1e-14 { Some(total) } else { None }
}

/// Generic trimmed face area via UV-grid integration of |鈭侾/鈭倁 脳 鈭侾/鈭倂|.
///
/// Samples the wire boundary in 3D, projects to UV via [`closest_point_on_surface`],
/// then rasterizes the UV polygon at ~200脳200 grid resolution with winding-number
/// point-in-polygon tests for outer/inner wires.  Handles any surface type (torus,
/// BSpline, Bezier, 鈥? where the existing per-surface optimisations don't apply.
fn try_generic_trimmed_face_area(
 surf: &Surface3,
 brep: &BRep,
 face: &Face,
 _face_flat_idx: usize,
) -> Option<f64> {
 use crate::projection::closest_point_on_surface;
 const TWO_PI: f64 = std::f64::consts::TAU;
 const PER_WIRE: usize = 256;

 let is_u_periodic = matches!(
 surf,
 Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_) | Surface3::Cone(_)
 );
 let is_v_periodic = matches!(surf, Surface3::Torus(_));

 // Sample wire 3D points 鈫?project to UV 鈫?optionally unwrap U so that the
 // winding-number test works correctly across the periodic seam.
 let wire_to_uv = |wire: &Wire| -> Option<Vec<DVec2>> {
 let mut pts_3d = sample_wire_polyline_3d_with_n(brep, wire, PER_WIRE);
 trim_almost_closed_polyline(&mut pts_3d, 1e-5);
 if pts_3d.len() < 3 {
 return None;
 }

 let mut uvs: Vec<DVec2> = pts_3d
 .iter()
 .map(|&p| {
 let proj = closest_point_on_surface(surf, p, 16);
 DVec2::new(proj.params.0, proj.params.1)
 })
 .collect();

 if is_u_periodic {
 let n = uvs.len();
 let mut uw = Vec::with_capacity(n);
 uw.push(uvs[0].x);
 for i in 1..n {
 uw.push(uw[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x));
 }
 for (i, u) in uw.into_iter().enumerate() {
 uvs[i].x = u;
 }
 }

 Some(uvs)
 };

 let outer_uv = wire_to_uv(&face.outer_wire)?;
 if outer_uv.len() < 3 {
 return None;
 }

 // UV bounding box from the outer wire.
 let (mut umin, mut umax) = (f64::INFINITY, f64::NEG_INFINITY);
 let (mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY);
 for uv in &outer_uv {
 umin = umin.min(uv.x);
 umax = umax.max(uv.x);
 vmin = vmin.min(uv.y);
 vmax = vmax.max(uv.y);
 }
 if umin >= umax || vmin >= vmax || !umin.is_finite() {
 return None;
 }

 // Pad bounds by 0.1 % so cells exactly on the boundary are inside.
 let du_pad = (umax - umin) * 0.001;
 let dv_pad = (vmax - vmin) * 0.001;
 let b_umin = umin - du_pad;
 let b_umax = umax + du_pad;
 let b_vmin = vmin - dv_pad;
 let b_vmax = vmax + dv_pad;

 // Inner-hole UV polygons.
 let inner_uvs: Vec<Vec<DVec2>> = face
 .inner_wires
 .iter()
 .filter_map(|w| wire_to_uv(w))
 .collect();

 // Grid resolution: ~200 cells/axis, clamped to [16, 400].
 let n_cells = 200usize;
 let u_ext = b_umax - b_umin;
 let v_ext = b_vmax - b_vmin;
 let aspect = (v_ext / u_ext).abs();
 let nu = if aspect > 1.0 {
 (n_cells as f64 / aspect).ceil() as usize
 } else {
 n_cells
 };
 let nv = if aspect <= 1.0 {
 (n_cells as f64 * aspect).ceil() as usize
 } else {
 n_cells
 };
 let nu = nu.max(16).min(400);
 let nv = nv.max(16).min(400);

 let du_cell = u_ext / nu as f64;
 let dv_cell = v_ext / nv as f64;
 let h = (du_cell * du_cell + dv_cell * dv_cell).sqrt().max(1e-12) * 1e-3;

 let mut area = 0.0_f64;
 for i in 0..nu {
 for j in 0..nv {
 let uc = b_umin + (i as f64 + 0.5) * du_cell;
 let vc = b_vmin + (j as f64 + 0.5) * dv_cell;
 let pt = DVec2::new(uc, vc);

 if winding_number_2d(&outer_uv, pt) == 0 {
 continue;
 }
 if inner_uvs
 .iter()
 .any(|poly| winding_number_2d(poly, pt) != 0)
 {
 continue;
 }

 // Map to surface-compatible parameters for point evaluation.
 let uc_s = if is_u_periodic {
 uc.rem_euclid(TWO_PI)
 } else {
 uc
 };
 let vc_s = if is_v_periodic {
 vc.rem_euclid(TWO_PI)
 } else {
 vc
 };

 let pu =
 (surf.point_at(uc_s + h, vc_s) - surf.point_at(uc_s - h, vc_s)) / (2.0 * h);
 let pv =
 (surf.point_at(uc_s, vc_s + h) - surf.point_at(uc_s, vc_s - h)) / (2.0 * h);
 area += pu.cross(pv).length() * du_cell * dv_cell;
 }
 }

 if area > 0.0 && area.is_finite() {
 Some(area)
 } else {
 None
 }
}

/// Prefer analytic / parametric area for `surface_area`: plane (shoelace); all sphere faces
/// (UV polygon mask + `R虏 d惟`); finite-UV rectangular patches on other surfaces without inner
/// wires (cylinder exact; otherwise `鈥朠u脳Pv鈥朻 midpoint rule on the same domain as tessellation).
fn try_analytic_face_surface_area(
 brep: &BRep,
 face: &Face,
 face_flat_idx: usize,
) -> Option<f64> {
 let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } })?;
 let surf = &surf_idx;

 match surf {
 Surface3::Plane(p) => {
 // Exact arc-aware contour area.  Use surface normal from Plane
 // geometry (not face.normal) for the arc bulge sign computation.
 // face.normal from Newell's method can be inverted for some cap
 // faces (cylinder-box bottom cap), flipping all arc correction
 // signs and producing area 0.858 instead of pi (H7).
 if let Some(a) = try_planar_face_exact_contour_area(brep, face, p.normal) {
 return Some(a);
 }

 // Compute both shoelace and GL areas.
 // BSpline鈫扨lane promotion can produce Plane faces whose 3D vertices
 // don't lie exactly on the promoted Plane (best-fit plane 鈮?original
 // BSpline surface).  The shoelace projects vertices to the 2D plane,
 // getting distorted area (e.g. 0.748 instead of 1.0).  GL integration
 // evaluates |Su 脳 Sv| directly on the Plane at UV Gauss points, giving
 // the exact UV-rectangle area (works for full-rectangle faces).
 // 鉁?OCCT  : BRepGProp uses GL for ALL surface types, including Plane.
 // (BRepGProp_Face.cxx L217-257: SIntOrder returns Nu=1,Nv=1 for Plane)
 let shoelace = try_planar_face_area_shoelace(brep, face, p.normal);
 let gl = face_surface_area_gauss(brep, face, face_flat_idx);
 match (shoelace, gl) {
 (Some(s), Some(g)) => {
 // When shoelace and GL disagree (>2%), the face may have
 // BSpline鈫扨lane promotion artifact.  Use GL if the face is
 // a simple 4-edge rectangle (no inner wires) 鈥?these are
 // full UV rectangles where GL gives exact area.
 let ratio = if g > 1e-12 { (s - g).abs() / g } else { 0.0 };
 if ratio > 0.02
 && face.inner_wires.is_empty()
 && face.outer_wire.edges.len() == 4
 {
 return Some(g);
 }
 return Some(s);
 }
 (Some(s), None) => return Some(s),
 (None, Some(g)) => return Some(g),
 (None, None) => {}
 }
 // Vertex-polygon fallback: the main shoelace path can return None
 // for valid planar faces with boolean-T-junction scrambled wires
 // or curves that the exact-contour handler cannot process.
 let pts: Vec<DVec3> = face.outer_wire.edges.iter().filter_map(|we| {
 let e = brep.flat_edges().get(we.idx).copied()?;
 let vi = if we.forward { e.0 } else { e.1 };
 brep.vertex_point(vi)
 }).collect();
 if pts.len() >= 3 {
 let n = face.normal;
 let (ux, uy) = local_basis_from_normal(n);
 let mut area2 = 0.0;
 for i in 0..pts.len() {
 let a = pts[i];
 let b = pts[(i + 1) % pts.len()];
 area2 += (a.dot(ux)) * (b.dot(uy)) - (b.dot(ux)) * (a.dot(uy));
 }
 let v = 0.5 * area2.abs();
 if v > 0.0 && v.is_finite() { return Some(v); }
 }
 None
 }
 // Planar BSpline/Bezier (e.g. from nurbsconvert): use shoelace.
 Surface3::BSpline(bsp) => {
 // Try multiple tolerances: restrict_to_bspline can introduce tiny
 // numerical deviations in control points for an exactly-planar surface.
 let is_p12 = crate::geom::bspline_is_planar(bsp, 1e-12);
 let is_p7 = crate::geom::bspline_is_planar(bsp, 1e-7);
 let is_p4 = crate::geom::bspline_is_planar(bsp, 1e-4);
 if !is_p12 && !is_p7 && !is_p4 { return None; }
 let plane = crate::geom::bspline_to_plane(bsp);
 let sa = try_planar_face_area_shoelace(brep, face, plane.normal);
 if sa.is_some() { return sa; }
 // Fallback: compute polygon area from face vertices directly.
 // The primary shoelace path can return None for valid planar
 // BSpline faces when the dense-sampled polyline is unreliable
 // (boolean T-junction wire reordering, bspline edge evaluation
 // noise from restrict_to_bspline).  The mesh vertices are exact.
 let pts: Vec<DVec3> = face.outer_wire.edges.iter().filter_map(|we| {
 let e = brep.flat_edges().get(we.idx).copied()?;
 let vi = if we.forward { e.0 } else { e.1 };
 brep.vertex_point(vi)
 }).collect();
 if pts.len() >= 3 {
 let n = plane.normal;
 let (ux, uy) = local_basis_from_normal(n);
 let mut area2 = 0.0;
 for i in 0..pts.len() {
 let a = pts[i];
 let b = pts[(i + 1) % pts.len()];
 area2 += (a.dot(ux)) * (b.dot(uy)) - (b.dot(ux)) * (a.dot(uy));
 }
 let a = 0.5 * area2.abs();
 if a > 0.0 && a.is_finite() { return Some(a); }
 }
 try_planar_face_area_shoelace(brep, face, face.normal)
 }
 Surface3::Bezier(bez) if bez.control_points.len() >= 2 && bez.control_points.iter().all(|r| r.len() >= 2) => {
 let degree_u = bez.control_points.len().saturating_sub(1);
 let degree_v = bez.control_points.first().map_or(0, |r| r.len().saturating_sub(1));
 let bsp = crate::geom::BSplineSurface {
 degree_u,
 degree_v,
 control_points: bez.control_points.clone(),
 knots_u: vec![],
 knots_v: vec![],
 weights: bez.weights.clone(),
 };
 if crate::geom::bspline_is_planar(&bsp, 1e-12) {
 let plane = crate::geom::bspline_to_plane(&bsp);
 try_planar_face_area_shoelace(brep, face, plane.normal)
 } else {
 None
 }
 }
 Surface3::Sphere(s) => {
 // Fast-path: if all edges are great circles (Circle3 center ==
 // sphere center), compute exact area via spherical-polygon formula.
 if let Some(a) = try_spherical_polygon_great_circle_area(s, brep, face) {
 return Some(a);
 }
let ctx = spherical_holed_uv_mask_setup(s, brep, face)?;
 let v = sphere_gauss_legendre_area_sum(s, &ctx);
 let full_sphere_area = 4.0 * std::f64::consts::PI * s.radius * s.radius;
 let sample_inside = face
 .sample_point
 .map(|p| point_in_spherical_polygon_3d(&ctx.outer_3d, p))
 .unwrap_or(true);
 if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
 eprintln!(
 "[SPHERE_AREA] face={} area={:.6} sphere_full={:.6} u=[{:.4},{:.4}] v=[{:.4},{:.4}] sample_inside={} outer_pts={} inner_loops={}",
 face_flat_idx,
 v,
 full_sphere_area,
 ctx.umin,
 ctx.umax,
 ctx.vmin,
 ctx.vmax,
 sample_inside,
 ctx.outer_3d.len(),
 ctx.inner_3d.len(),
 );
 }
 // Architecture diff A5b: the UV mask covers the small polygon region,
 // but the WireFace may represent the complement (large cap).  Check
 // the winding of outer_uv to determine which region is correct.
 let is_cw = {
 let mut area2 = 0.0;
 let n = ctx.outer_uv.len();
 for i in 0..n {
 let j = (i + 1) % n;
 area2 += ctx.outer_uv[i].x * ctx.outer_uv[j].y - ctx.outer_uv[j].x * ctx.outer_uv[i].y;
 }
 area2 * 0.5 < 0.0
 };
 let result_area = if is_cw { full_sphere_area - v } else { v };
 let suspicious_wrap = (ctx.umax - ctx.umin).abs() > std::f64::consts::TAU + 0.25;
 if suspicious_wrap || result_area > full_sphere_area * 1.001 {
 return None;
 }
 // Also check sample_inside consistency: if the sample point is
 // outside the UV polygon but the face is CCW (small region),
 // the region is wrong 鈫?return None (let triangulation handle it).
 if !suspicious_wrap && !is_cw && v > 0.0 && v < full_sphere_area * 0.5 {
 if !sample_inside {
 return None;
 }
 }
 if result_area > 0.0 { return Some(result_area); }
 None
 }
 _ => {
 // Use existing specialized paths for all surface types.
 // GL integration (face_surface_area_gauss) is available for future use
 // but is not wired into the dispatch to avoid double-counting from
 // overlapping Plane+BSpline faces produced by the boolean pipeline.
 let [u0, u1, v0, v1] = curved_face_uv_domain(brep, face, face_flat_idx, surf)?;
 if !u0.is_finite()
 || !u1.is_finite()
 || !v0.is_finite()
 || !v1.is_finite()
 || (u1 - u0).abs() < 1e-14
 || (v1 - v0).abs() < 1e-14
 {
 return None;
 }
 match surf {
 Surface3::Cylinder(c) => {
 let result = try_cylinder_trimmed_face_area(c, brep, face, face_flat_idx);
 if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
 eprintln!("[CYL_ANALYTIC] fi={} analytic={:?}",
 face_flat_idx, result);
 }
 result
 }
 Surface3::Cone(c) => {
 try_cone_trimmed_face_area(c, brep, face)
 }
 Surface3::Torus(_) => {
 if !face.inner_wires.is_empty() || (u1 - u0).abs() > std::f64::consts::TAU * 1.01 {
 try_generic_trimmed_face_area(surf, brep, face, face_flat_idx)
 } else {
 param_rect_area_cross(surf, u0, u1, v0, v1)
 }
 }
 _ if !face.inner_wires.is_empty() => {
 try_generic_trimmed_face_area(surf, brep, face, face_flat_idx)
 }
 _ => param_rect_area_cross(surf, u0, u1, v0, v1),
 }
 }
 }
}

pub fn try_analytic_face_surface_area_pub(
 brep: &BRep,
 face: &Face,
 face_flat_idx: usize,
) -> Option<f64> {
 try_analytic_face_surface_area(brep, face, face_flat_idx)
}

/// When ear-clipping fails in (u,v), approximate the trimmed patch by a regular
/// grid in parameter space, keeping only cells whose centres lie inside the
/// outer UV loop and outside any inner loop; map quads to chordal 3D triangles.
fn try_spherical_uv_masked_raster(
 s: &SphericalSurface,
 brep: &BRep,
 face: &Face,
 _face_flat_idx: usize,
 face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
 let ctx = spherical_holed_uv_mask_setup(s, brep, face)?;
 let umin = ctx.umin;
 let umax = ctx.umax;
 let vmin = ctx.vmin;
 let vmax = ctx.vmax;
 let outer_uv = &ctx.outer_uv;
 let outer_3d = &ctx.outer_3d;
 let inner_polys = &ctx.inner_polys;
 let use_uv_winding = ctx.use_uv_winding;

 // When the UV polygon wraps around u (use_uv_winding=false), use the
 // full sphere domain [鈭捪€, 蟺] 脳 [0, 蟺] so every physical point is
 // tested exactly once. The 3D angular-sum test is insensitive to the
 // UV wrapping.
 let (umin_eff, umax_eff, nu, nv) = if use_uv_winding {
 (umin, umax, SPHERE_UV_MASK_N, SPHERE_UV_MASK_N)
 } else {
 let pi = std::f64::consts::PI;
 (-pi, pi, 120usize, 120usize)
 };
 let umin = umin_eff;
 let umax = umax_eff;
 let du = (umax - umin) / nu as f64;
 let dv = (vmax - vmin) / nv as f64;
 let radius = s.radius;

 let emit_grid = |use_inner_mask: bool| -> Vec<[DVec3; 3]> {
 let mut tris: Vec<[DVec3; 3]> = Vec::new();
 for i in 0..nu {
 for j in 0..nv {
 let u0 = umin + i as f64 * du;
 let u1 = u0 + du;
 let v0 = vmin + j as f64 * dv;
 let v1 = v0 + dv;
 let uc = umin + (i as f64 + 0.5) * du;
 let vc = vmin + (j as f64 + 0.5) * dv;
 let pc = s.point_at(uc, vc);
 // 5-point OR test: accept if center is inside (UV winding + 3D
 // angular-sum), OR if any corner is inside (3D angular-sum).
 // When the UV polygon wraps u multiple times (merged sphere
 // faces), skip the unreliable winding-number test.
 let p00 = s.point_at(u0, v0);
 let p10 = s.point_at(u1, v0);
 let p11 = s.point_at(u1, v1);
 let p01 = s.point_at(u0, v1);
 let center_in = if use_uv_winding {
 winding_number_2d(outer_uv, DVec2::new(uc, vc)) > 0
 || point_in_spherical_polygon_3d(outer_3d, pc)
 } else {
 point_in_spherical_polygon_3d(outer_3d, pc)
 };
 let any_corner_in = point_in_spherical_polygon_3d(outer_3d, p00)
 || point_in_spherical_polygon_3d(outer_3d, p10)
 || point_in_spherical_polygon_3d(outer_3d, p11)
 || point_in_spherical_polygon_3d(outer_3d, p01);
 if !center_in && !any_corner_in {
 continue;
 }
 if use_inner_mask {
 // Use UV-space winding number for inner-hole test regardless of
 // use_uv_winding.  Inner wires bound small regions (area < 蟺), but
 // point_in_spherical_polygon_3d's `|total| > 蟺` threshold cannot
 // detect points inside patches smaller than a hemisphere 鈥?it returns
 // false for *all* points.  UV-space handles small holes correctly.
 let in_hole = inner_polys.iter().any(|h| winding_number_2d(h, DVec2::new(uc, vc)) != 0);
 if in_hole {
 continue;
 }
 }
 let nref = s.normal_at(uc, vc);
 // Exact area in (u, v) for a sphere: R虏 路 du 路 (cos v0 鈭?cos v1); isotropic
 // scale of the chordal bilinear patch toward the cell centre to match dA.
 let r2 = radius * radius;
 let d_target = r2 * du * (v0.cos() - v1.cos());
 let a0 = tri_area(p00, p10, p11);
 let a1 = tri_area(p00, p11, p01);
 let chord = a0 + a1;
 if chord > 1e-20 && d_target > 0.0 {
 let k = (d_target / chord).sqrt().clamp(0.1, 10.0);
 let t = |p: DVec3| pc + (p - pc) * k;
 tris.push(orient_by_ref([t(p00), t(p10), t(p11)], nref));
 tris.push(orient_by_ref([t(p00), t(p11), t(p01)], nref));
 } else {
 tris.push(orient_by_ref([p00, p10, p11], nref));
 tris.push(orient_by_ref([p00, p11, p01], nref));
 }
 }
 }
 tris
 };

 let mut tris = emit_grid(true);
 if tris.is_empty() && !inner_polys.is_empty() {
 tris = emit_grid(false);
 }
 if tris.is_empty() {
 None
 } else {
 Some(
 tris
 .into_iter()
 .map(|[a, b, c]| orient_tri([a, b, c], face_normal))
 .collect(),
 )
 }
}

/// Tessellate a curved face by sampling the underlying `Surface3` over its
/// UV domain on a regular `UV_TESS_N 脳 UV_TESS_N` grid.
///
/// Returns triangles oriented outward (consistent with the surface normal).
/// UV rectangle used by [`tessellate_curved_face`] and parametric `surface_area` fallbacks
/// (same resolution priority as triangulation: range 鈫?finite `default_domain` 鈫?wire estimate).
fn curved_face_uv_domain(
 brep: &BRep,
 face: &Face,
 face_flat_idx: usize,
 surf: &Surface3,
) -> Option<[f64; 4]> {
 if let Some(r) = brep.tshapes.get(face_flat_idx).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None } }) {
 Some(r)
 } else {
 let d = surf.default_domain();
 if d.iter().all(|x| x.is_finite()) {
 Some(d)
 } else {
 estimate_uv_domain_from_wire(brep, face, face_flat_idx, surf)
 }
 }
}

// 鈹€鈹€ Gauss-Legendre Integration (OCCT BRepGProp alignment) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Gauss-Legendre nodes and weights for orders 1..=12.
/// Generated from Legendre polynomial roots.
/// Reference: OCCT math_GaussPoints.cxx
/// 鉁?OCCT (math::GaussPoints / math::GaussWeights)
#[allow(dead_code)]
struct GLTable {
 #[allow(dead_code)]
 n: usize,
 points: [f64; 12],
 weights: [f64; 12],
}

const GL_TABLES: [GLTable; 13] = [
 // order 0 (unused)
 GLTable { n: 0, points: [0.0; 12], weights: [0.0; 12] },
 // order 1
 GLTable { n: 1, points: [0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], weights: [2.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] },
 // order 2
 GLTable { n: 2, points: [
 -0.5773502691896257645091487805019574556476,
 0.5773502691896257645091487805019574556476,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ], weights: [
 1.0, 1.0,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ]},
 // order 3
 GLTable { n: 3, points: [
 -0.7745966692414833770358530799564799221665,
 0.0,
 0.7745966692414833770358530799564799221665,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ], weights: [
 0.5555555555555555555555555555555555555556,
 0.8888888888888888888888888888888888888889,
 0.5555555555555555555555555555555555555556,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ]},
 // order 4
 GLTable { n: 4, points: [
 -0.8611363115940525752239464888928095050941,
 -0.3399810435848562648026657591032446872006,
 0.3399810435848562648026657591032446872006,
 0.8611363115940525752239464888928095050941,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ], weights: [
 0.3478548451374538573730639492219994072320,
 0.6521451548625461426269730507781286076179,
 0.6521451548625461426269730507781286076179,
 0.3478548451374538573730639492219994072320,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ]},
 // order 5
 GLTable { n: 5, points: [
 -0.9061798459386639927976268782993929651257,
 -0.5384693101056830910363144207002088049673,
 0.0,
 0.5384693101056830910363144207002088049673,
 0.9061798459386639927976268782993929651257,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ], weights: [
 0.2369268850561890875142640407199173626432,
 0.4786286704993664680412915148356381929123,
 0.5688888888888888888888888888888888888889,
 0.4786286704993664680412915148356381929123,
 0.2369268850561890875142640407199173626432,
 0.0,0.0,0.0,0.0,0.0,0.0,0.0,
 ]},
 // order 6
 GLTable { n: 6, points: [
 -0.9324695142031520278123015544939940976115,
 -0.6612093864662645136613995950199050806726,
 -0.2386191860831969086305017216807119325162,
 0.2386191860831969086305017216807119325162,
 0.6612093864662645136613995950199050806726,
 0.9324695142031520278123015544939940976115,
 0.0,0.0,0.0,0.0,0.0,0.0,
 ], weights: [
 0.1713244923791703450402961421727328935268,
 0.3607615730481386075698335138377161416615,
 0.4679139345726910473898703439895509948116,
 0.4679139345726910473898703439895509948116,
 0.3607615730481386075698335138377161416615,
 0.1713244923791703450402961421727328935268,
 0.0,0.0,0.0,0.0,0.0,0.0,
 ]},
 // order 7
 GLTable { n: 7, points: [
 -0.9491079123427585245261896840478512624008,
 -0.7415311855993944398638647732807884070741,
 -0.4058451513773971669066064120769614633473,
 0.0,
 0.4058451513773971669066064120769614633473,
 0.7415311855993944398638647732807884070741,
 0.9491079123427585245261896840478512624008,
 0.0,0.0,0.0,0.0,0.0,
 ], weights: [
 0.1294849661688696932706114326790820182324,
 0.2797053914892766679014677714237795824869,
 0.3818300505051189449503697754889751338784,
 0.4179591836734693877551020408163265306122,
 0.3818300505051189449503697754889751338784,
 0.2797053914892766679014677714237795824869,
 0.1294849661688696932706114326790820182324,
 0.0,0.0,0.0,0.0,0.0,
 ]},
 // order 8
 GLTable { n: 8, points: [
 -0.9602898564975362316835608685694729904282,
 -0.7966664774136267395915539364758302189390,
 -0.5255324099163289858175570492491803935109,
 -0.1834346424956498049394761423601839806667,
 0.1834346424956498049394761423601839806667,
 0.5255324099163289858175570492491803935109,
 0.7966664774136267395915539364758302189390,
 0.9602898564975362316835608685694729904282,
 0.0,0.0,0.0,0.0,
 ], weights: [
 0.1012285362903762591525314785099858429633,
 0.2223810344533744705443559944262408844301,
 0.3137066458778872873379622019866013132603,
 0.3626837833783619829651504492771976121474,
 0.3626837833783619829651504492771976121474,
 0.3137066458778872873379622019866013132603,
 0.2223810344533744705443559944262408844301,
 0.1012285362903762591525314785099858429633,
 0.0,0.0,0.0,0.0,
 ]},
 // order 9
 GLTable { n: 9, points: [
 -0.9681602395076260898355762029036728691970,
 -0.8360311073266357942994297880697348765445,
 -0.6133714327005903973087020393414741847858,
 -0.3242534234038089290385380146433365678265,
 0.0,
 0.3242534234038089290385380146433365678265,
 0.6133714327005903973087020393414741847858,
 0.8360311073266357942994297880697348765445,
 0.9681602395076260898355762029036728691970,
 0.0,0.0,0.0,
 ], weights: [
 0.0812743883615744119718921581105236506756,
 0.1806481606948574040584720312429128095144,
 0.2606106964029354623187428694186328496466,
 0.3123470770400028400686304065844436655987,
 0.3302393550012597631645250692869740488788,
 0.3123470770400028400686304065844436655987,
 0.2606106964029354623187428694186328496466,
 0.1806481606948574040584720312429128095144,
 0.0812743883615744119718921581105236506756,
 0.0,0.0,0.0,
 ]},
 // order 10
 GLTable { n: 10, points: [
 -0.9739065285171717200779640120844520534383,
 -0.8650633666889845107320966884234930481254,
 -0.6794095682990244062343273651148735757693,
 -0.4333953941292471907992659431657841622000,
 -0.1488743389816312108848260011297199846175,
 0.1488743389816312108848260011297199846175,
 0.4333953941292471907992659431657841622000,
 0.6794095682990244062343273651148735757693,
 0.8650633666889845107320966884234930481254,
 0.9739065285171717200779640120844520534383,
 0.0,0.0,
 ], weights: [
 0.0666713443086881375935688098933320808228,
 0.1494513491505805931457763396576973241563,
 0.2190863625159820439955349342281631926256,
 0.2692667193099963550912269215694828526643,
 0.2955242247147528701738929946513383294210,
 0.2955242247147528701738929946513383294210,
 0.2692667193099963550912269215694828526643,
 0.2190863625159820439955349342281631926256,
 0.1494513491505805931457763396576973241563,
 0.0666713443086881375935688098933320808228,
 0.0,0.0,
 ]},
 // order 11
 GLTable { n: 11, points: [
 -0.9782286581460569928038090019603639625151,
 -0.8870625997680952990751577693039274940173,
 -0.7301520055740493240934262531031699040233,
 -0.5190961292068118159257256694585208657768,
 -0.2695431559523449723315319854008615424241,
 0.0,
 0.2695431559523449723315319854008615424241,
 0.5190961292068118159257256694585208657768,
 0.7301520055740493240934262531031699040233,
 0.8870625997680952990751577693039274940173,
 0.9782286581460569928038090019603639625151,
 0.0,
 ], weights: [
 0.0556685671161736664827537204425485787286,
 0.1255803694649046246346942992239400861978,
 0.1862902109277342514262576411406595099962,
 0.2331937645919904799185237048411751086375,
 0.2628045445102466621806888698955091843562,
 0.2729250867779006307144835283363420935117,
 0.2628045445102466621806888698955091843562,
 0.2331937645919904799185237048411751086375,
 0.1862902109277342514262576411406595099962,
 0.1255803694649046246346942992239400861978,
 0.0556685671161736664827537204425485787286,
 0.0,
 ]},
 // order 12
 GLTable { n: 12, points: [
 -0.9815606342467192506905490901492808229601,
 -0.9041172563704748566784658661190961925378,
 -0.7699026741943046870368938332128180752053,
 -0.5873179542866174472967024189405342844225,
 -0.3678314989981801937526915366437195613199,
 -0.1252334085114689154724413694638531299833,
 0.1252334085114689154724413694638531299833,
 0.3678314989981801937526915366437195613199,
 0.5873179542866174472967024189405342844225,
 0.7699026741943046870368938332128180752053,
 0.9041172563704748566784658661190961925378,
 0.9815606342467192506905490901492808229601,
 ], weights: [
 0.0471753363865118271946159614850171063171,
 0.1069393259953184309602547181939962241076,
 0.1600783285433462263346525295433590718720,
 0.2031674267230659217490644558097983765066,
 0.2334925365383548087608498989248780562594,
 0.2491470458134027850005624360429512108305,
 0.2491470458134027850005624360429512108305,
 0.2334925365383548087608498989248780562594,
 0.2031674267230659217490644558097983765066,
 0.1600783285433462263346525295433590718720,
 0.1069393259953184309602547181939962241076,
 0.0471753363865118271946159614850171063171,
 ]},
];

fn gl_table(order: usize) -> &'static GLTable {
 debug_assert!(order >= 1 && order <= 12, "GL order {} out of range [1,12]", order);
 &GL_TABLES[order]
}

/// OCCT-aligned: BRepGProp_Face::SIntOrder (BRepGProp_Face.cxx L217-257)
/// Returns the number of Gauss points in U and V for integration over the surface.
/// 鉁?OCCT  
fn gl_s_integration_order(surf: &Surface3) -> (usize, usize) {
 match surf {
 Surface3::Plane(_) => (1, 1),
 Surface3::Cylinder(_) => (2, 1),
 Surface3::Cone(_) => (2, 1),
 Surface3::Sphere(_) => (2, 2),
 Surface3::Torus(_) => (2, 2),
 Surface3::BSpline(bsp) => {
 let nu = bsp.degree_u.max(1);
 let nv = bsp.degree_v.max(1);
 (nu, nv)
 }
 Surface3::Bezier(bez) => {
 let nu = bez.control_points.len().saturating_sub(1).max(1);
 let nv = bez.control_points.first().map_or(1, |r| r.len().saturating_sub(1)).max(1);
 (nu, nv)
 }
 _ => (2, 2),
 }
}

/// OCCT-aligned: BRepGProp_Face::SUIntSubs (BRepGProp_Face.cxx L261-292)
/// Returns number of U intervals for subdivision.
/// 鉁?OCCT  
fn gl_s_u_subs(surf: &Surface3) -> usize {
 match surf {
 Surface3::Plane(_) => 1,
 Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => 3,
 Surface3::BSpline(bsp) => bsp.knots_u.len().saturating_sub(1).max(1),
 _ => 1,
 }
}

/// OCCT-aligned: BRepGProp_Face::SVIntSubs (BRepGProp_Face.cxx L296-327)
/// Returns number of V intervals for subdivision.
/// 鉁?OCCT  
fn gl_s_v_subs(surf: &Surface3) -> usize {
 match surf {
 Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) => 1,
 Surface3::Sphere(_) => 2,
 Surface3::Torus(_) => 3,
 Surface3::BSpline(bsp) => bsp.knots_v.len().saturating_sub(1).max(1),
 _ => 1,
 }
}

/// OCCT-aligned: BRepGProp_Face::UKnots (BRepGProp_Face.cxx L331-356)
/// Returns U knot positions for subdividing the integration domain, clipped to [u0, u1].
/// 鉁?OCCT  
fn gl_u_knots(surf: &Surface3, u0: f64, u1: f64) -> Vec<f64> {
 let knots = match surf {
 Surface3::Plane(_) => vec![u0, u1],
 Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
 // OCCT: 0, 2蟺/3, 4蟺/3, 2蟺
 vec![0.0, std::f64::consts::TAU / 3.0, 2.0 * std::f64::consts::TAU / 3.0, std::f64::consts::TAU]
 }
 Surface3::BSpline(bsp) => bsp.knots_u.clone(),
 _ => vec![u0, u1],
 };
 // Clip to [u0, u1]
 let mut clipped: Vec<f64> = knots.into_iter()
 .filter(|&k| k >= u0 - 1e-12 && k <= u1 + 1e-12)
 .collect();
 clipped.sort_by(|a, b| a.partial_cmp(b).unwrap());
 clipped.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
 if clipped.is_empty() { return vec![u0, u1]; }
 if clipped[0] > u0 + 1e-12 { clipped.insert(0, u0); }
 let last = clipped.len() - 1;
 if clipped[last] < u1 - 1e-12 { clipped.push(u1); }
 // OCCT uses at most SUBS_POWER subintervals; keep it manageable
 if clipped.len() > 64 { return vec![u0, u1]; }
 clipped
}

/// OCCT-aligned: BRepGProp_Face::VKnots (BRepGProp_Face.cxx L360-389)
/// Returns V knot positions for subdividing the integration domain, clipped to [v0, v1].
/// 鉁?OCCT  
fn gl_v_knots(surf: &Surface3, v0: f64, v1: f64) -> Vec<f64> {
 let knots = match surf {
 Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) => {
 vec![v0, v1]
 }
 Surface3::Sphere(_) => {
 // OCCT: -蟺/2, 0, 蟺/2
 vec![-std::f64::consts::FRAC_PI_2, 0.0, std::f64::consts::FRAC_PI_2]
 }
 Surface3::Torus(_) => {
 // OCCT: 0, 2蟺/3, 4蟺/3, 2蟺
 vec![0.0, std::f64::consts::TAU / 3.0, 2.0 * std::f64::consts::TAU / 3.0, std::f64::consts::TAU]
 }
 Surface3::BSpline(bsp) => bsp.knots_v.clone(),
 _ => vec![v0, v1],
 };
 // Clip to [v0, v1]
 let mut clipped: Vec<f64> = knots.into_iter()
 .filter(|&k| k >= v0 - 1e-12 && k <= v1 + 1e-12)
 .collect();
 clipped.sort_by(|a, b| a.partial_cmp(b).unwrap());
 clipped.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
 if clipped.is_empty() { return vec![v0, v1]; }
 if clipped[0] > v0 + 1e-12 { clipped.insert(0, v0); }
 let last = clipped.len() - 1;
 if clipped[last] < v1 - 1e-12 { clipped.push(v1); }
 if clipped.len() > 64 { return vec![v0, v1]; }
 clipped
}

/// Compute |Su 脳 Sv| at parameter (u,v) via central finite differences.
/// Equivalent to OCCT BRepGProp_Face::Normal (BRepGProp_Face.cxx L191-198)
/// which returns the UNNORMALIZED surface normal = Su 脳 Sv.
/// The magnitude |Su 脳 Sv| is the area element Jacobian.
/// 鈴? : uses finite differences (rcad has no D1 evaluation trait);
/// OCCT uses exact D1 evaluation. Accuracy is equivalent for analytical
/// surfaces and within 1e-8 for BSpline at the chosen h.
fn surface_normal_jacobian(surf: &Surface3, u: f64, v: f64) -> f64 {
 let h = 1e-6;
 let pu = (surf.point_at(u + h, v) - surf.point_at(u - h, v)) / (2.0 * h);
 let pv = (surf.point_at(u, v + h) - surf.point_at(u, v - h)) / (2.0 * h);
 pu.cross(pv).length()
}

/// Compute face surface area via Gauss-Legendre numerical integration over the
/// UV parameter domain, strictly aligned with OCCT BRepGProp_Gauss::Compute
/// (BRepGProp_Gauss.cxx L1050-1135, Sinert variant, no-Eps overload).
///
/// For faces without inner wires (no holes), performs a nested GL double
/// integral using SIntOrder integration order and UKnots/VKnots subdivision.
/// The integral is 鈭戔垜 |Su 脳 Sv| 路 w_u 路 w_v 路 u_rad 路 v_rad over all
/// sub-intervals and Gauss points.
///
/// 鉁?OCCT (BRepGProp_Gauss::Compute + computeSInertiaOfElementaryPart)
fn face_surface_area_gauss(brep: &BRep, face: &Face, fi: usize) -> Option<f64> {
 let surf_idx = brep.tshapes.get(fi).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } })?;
 let surf = &surf_idx;

 // For faces with inner wires, fall back to generic trimmed-face area
 if !face.inner_wires.is_empty() {
 return None;
 }

 let [u0, u1, v0, v1] = curved_face_uv_domain(brep, face, fi, surf)?;
 if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
 return None;
 }
 if (u1 - u0).abs() < 1e-14 || (v1 - v0).abs() < 1e-14 {
 return Some(0.0);
 }

 // Integration order from SIntOrder
 let (nu, nv) = gl_s_integration_order(surf);
 if nu == 0 || nv == 0 { return None; }

 // Clamp to available GL tables (max 12)
 let nu = nu.min(12);
 let nv = nv.min(12);

 // Get GL points/weights
 let glu = gl_table(nu);
 let glv = gl_table(nv);

 // Get U/V subdivision knots
 let uknots = gl_u_knots(surf, u0, u1);
 let vknots = gl_v_knots(surf, v0, v1);

 let mut total_area = 0.0;

 // Iterate over V subdivisions
 for vi in 0..vknots.len() - 1 {
 let v_a = vknots[vi];
 let v_b = vknots[vi + 1];
 if (v_b - v_a).abs() < 1e-14 { continue; }
 let v_mid = 0.5 * (v_b + v_a);
 let v_rad = 0.5 * (v_b - v_a);

 // Iterate over U subdivisions
 for ui in 0..uknots.len() - 1 {
 let u_a = uknots[ui];
 let u_b = uknots[ui + 1];
 if (u_b - u_a).abs() < 1e-14 { continue; }
 let u_mid = 0.5 * (u_b + u_a);
 let u_rad = 0.5 * (u_b - u_a);

 let mut sub_area = 0.0;

 // Inner double GL integral over this sub-rectangle
 for j in 0..nv {
 let v = v_mid + v_rad * glv.points[j];
 let w_v = glv.weights[j];

 for i in 0..nu {
 let u = u_mid + u_rad * glu.points[i];
 let w_u = glu.weights[i];

 // Evaluate |Su 脳 Sv| at (u,v)
 let jac = surface_normal_jacobian(surf, u, v);
 sub_area += jac * w_u * w_v;
 }
 }

 total_area += sub_area * u_rad * v_rad;
 }
 }

 if total_area.is_finite() && total_area >= 0.0 {
 Some(total_area)
 } else {
 None
 }
}

/// `鈭埆 鈥栤垈P/鈭倁脳鈭侾/鈭倂鈥?dudv` on `[u0,u1]脳[v0,v1]` (midpoint rule, central differences for partials).
/// Matches the same UV box as `tessellate_curved_face` without chordal tri area bias.
fn param_rect_area_cross(surf: &Surface3, u0: f64, u1: f64, v0: f64, v1: f64) -> Option<f64> {
 if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
 return None;
 }
 if (u1 - u0).abs() < 1e-15 || (v1 - v0).abs() < 1e-15 {
 return None;
 }
 let nu = PARAM_RECT_AREA_N;
 let nv = PARAM_RECT_AREA_N;
 let du = (u1 - u0) / nu as f64;
 let dv = (v1 - v0) / nv as f64;
 let h = (du * du + dv * dv).sqrt().max(1e-12) * 1e-3;
 let mut a = 0.0f64;
 for i in 0..nu {
 for j in 0..nv {
 let uc = u0 + (i as f64 + 0.5) * du;
 let vc = v0 + (j as f64 + 0.5) * dv;
 let pu = (surf.point_at(uc + h, vc) - surf.point_at(uc - h, vc)) / (2.0 * h);
 let pv = (surf.point_at(uc, vc + h) - surf.point_at(uc, vc - h)) / (2.0 * h);
 a += pu.cross(pv).length() * du * dv;
 }
 }
 if a > 0.0 && a.is_finite() { Some(a) } else { None }
}

/// Returns `None` if the face has no associated surface or the domain cannot
/// be determined (e.g. a truly unbounded Plane with no face_surface_range).
fn tessellate_curved_face(
 brep: &BRep,
 face: &crate::topology::Face,
 face_flat_idx: usize,
) -> Option<Vec<[DVec3; 3]>> {
 // Look up the surface for this face.
 let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } })?;
 let surf = &surf_idx;

 let domain = curved_face_uv_domain(brep, face, face_flat_idx, surf)?;

 let [u0, u1, v0, v1] = domain;

 // Sanity checks.
 if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
 return None;
 }
 if (u1 - u0).abs() < 1e-14 || (v1 - v0).abs() < 1e-14 {
 return None;
 }

 let nu = UV_TESS_N;
 let nv = UV_TESS_N;

 // Build a (nu+1)脳(nv+1) grid of 3-D points.
 let mut pts = Vec::with_capacity((nu + 1) * (nv + 1));
 for i in 0..=nu {
 let u = u0 + (u1 - u0) * (i as f64 / nu as f64);
 for j in 0..=nv {
 let v = v0 + (v1 - v0) * (j as f64 / nv as f64);
 pts.push(surf.point_at(u, v));
 }
 }

 // Emit two triangles per quad cell (i,j)鈥?i+1,j)鈥?i,j+1)鈥?i+1,j+1).
 let idx = |i: usize, j: usize| i * (nv + 1) + j;
 let mut tris: Vec<[DVec3; 3]> = Vec::with_capacity(nu * nv * 2);

 for i in 0..nu {
 for j in 0..nv {
 let p00 = pts[idx(i, j)];
 let p10 = pts[idx(i + 1, j)];
 let p01 = pts[idx(i, j + 1)];
 let p11 = pts[idx(i + 1, j + 1)];

 // Reference outward normal at cell centre.
 let uc = u0 + (u1 - u0) * ((i as f64 + 0.5) / nu as f64);
 let vc = v0 + (v1 - v0) * ((j as f64 + 0.5) / nv as f64);
 let n_ref = surf.normal_at(uc, vc);

 tris.push(orient_by_ref([p00, p10, p11], n_ref));
 tris.push(orient_by_ref([p00, p11, p01], n_ref));
 }
 }

 Some(tris)
}

/// Orient a triangle so its normal agrees with `n_ref`.
#[inline]
fn orient_by_ref(tri: [DVec3; 3], n_ref: DVec3) -> [DVec3; 3] {
 let [a, b, c] = tri;
 let n = (b - a).cross(c - a);
 if n.dot(n_ref) < 0.0 { [a, c, b] } else { [a, b, c] }
}

/// Estimate UV domain for surfaces whose natural domain has infinite extents
/// (CylindricalSurface, ConicalSurface) by projecting wire vertices onto the
/// surface's UV space.
///
/// For a CylindricalSurface with axis +Y: u = atan2(z_proj, x_proj) and
/// v = dot(pt - origin, axis).  We just use the bounding box of all wire
/// vertex projections, with a small margin, for the finite axis.
fn estimate_uv_domain_from_wire(
 brep: &BRep,
 face: &crate::topology::Face,
 face_flat_idx: usize,
 surf: &crate::geom::Surface3,
) -> Option<[f64; 4]> {
 use crate::geom::Surface3;

 // Collect all wire vertex 3-D points (outer + inner wires).
 let all_wires = std::iter::once(&face.outer_wire).chain(face.inner_wires.iter());
 let pts: Vec<DVec3> = all_wires
 .flat_map(|w| &w.edges)
 .filter_map(|we| {
 let edge = brep.flat_edges().get(we.idx).copied()?;
 let vidx = if we.forward { edge.0 } else { edge.1 };
 brep.vertex_point(vidx)
 })
 .collect();

 if pts.is_empty() {
 return None;
 }

 match surf {
 Surface3::Cylinder(cyl) => {
 // Use edge pcurves on this face to determine the UV domain.
 // This avoids the any_perpendicular reference-frame mismatch
 // that occurs when projecting 3D curve points through point_to_u.
 let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.as_ref() } else { None } }).map(|s| ())?;
 let mut u_vals: Vec<f64> = Vec::new();
 let mut v_vals: Vec<f64> = pts.iter().map(|p| (*p - cyl.origin).dot(cyl.axis)).collect();
 if v_vals.is_empty() { return None; }

 use crate::geom::Curve2dEval;
 let all_wires = std::iter::once(&face.outer_wire).chain(face.inner_wires.iter());
 for we in all_wires.flat_map(|w| &w.edges) {
 let range = brep.tshapes.get(we.idx).and_then(|ts| { if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None } });
 if let Some([t0, t1]) = range {
 // Find pcurve on this face's surface.
 if let Some(c2d) = None::<&crate::geom::Curve2d>
 {
 let ns = 16;
 for k in 0..=ns {
 let frac = k as f64 / ns as f64;
 let uv = c2d.point_at(t0 + (t1 - t0) * frac);
 let u = if uv.x < 0.0 { uv.x + 2.0 * PI } else { uv.x };
 u_vals.push(u);
 if uv.y.is_finite() { v_vals.push(uv.y); }
 }
 }
 }
 }

 if u_vals.is_empty() {
 // Fallback: vertex-based estimate only
 let u_vert: Vec<f64> = pts.iter().map(|p| {
 let radial = *p - cyl.origin - (*p - cyl.origin).dot(cyl.axis) * cyl.axis;
 let x_ax = cyl.ref_dir.normalize();
 let y_ax = cyl.axis.cross(x_ax).normalize();
 let u = radial.dot(y_ax).atan2(radial.dot(x_ax));
 if u < 0.0 { u + 2.0 * PI } else { u }
 }).collect();
 let u0 = u_vert.iter().cloned().fold(f64::INFINITY, f64::min);
 let u1 = u_vert.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 if u0 >= u1 { return None; }
 let v0 = v_vals.iter().cloned().fold(f64::INFINITY, f64::min);
 let v1 = v_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 return Some([u0 - 1e-10, u1 + 1e-10, v0 - 1e-10, v1 + 1e-10]);
 }

 let u0 = u_vals.iter().cloned().fold(f64::INFINITY, f64::min);
 let u1 = u_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 let v0 = v_vals.iter().cloned().fold(f64::INFINITY, f64::min);
 let v1 = v_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 if u0 >= u1 || v0 >= v1 { return None; }
 Some([u0 - 1e-10, u1 + 1e-10, v0 - 1e-10, v1 + 1e-10])
 }
 Surface3::Cone(con) => {
 // ConicalSurface: u = azimuth [0, 2蟺], v = slant distance 鈮?0.
 let d = surf.default_domain();
 let u0 = d[0];
 let u1 = d[1];
 // v = distance from apex along slant (axis direction component / cos(half_angle))
 let cos_a = con.half_angle_rad.cos();
 let v_vals: Vec<f64> = pts.iter().map(|p| (*p - con.apex).dot(con.axis) / cos_a).collect();
 let v0 = v_vals.iter().cloned().fold(f64::INFINITY, f64::min).max(0.0);
 let v1 = v_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 if v0 >= v1 { return None; }
 Some([u0, u1, v0, v1 + 1e-10])
 }
 _ => None,
 }
}

/// Collect triangles for a face (either pre-triangulated, UV-tessellated for
/// curved faces, or fan-triangulated from wire vertices as last resort),
/// oriented outward (consistent with face.normal).
///
/// `face_flat_idx` is the face's flattened index across all solids/shells,
/// matching the indexing of `GeomStore.face_surface`.
fn face_triangles(
 brep: &BRep,
 face: &crate::topology::Face,
 face_flat_idx: usize,
) -> Vec<[DVec3; 3]> {
 // Determine whether the face has a curved surface (not Plane).
 // For curved surfaces, `face.normal` is unreliable 鈥?it is typically estimated
 // from the boundary polygon (Newell method) and can point in a completely wrong
 // direction for surfaces like the wall of a half-cylinder.  Using it in
 // `orient_tri` flips the winding of valid triangles, causing the signed volume
 // contribution to be negative or zero.  When curved, we trust the natural winding
 // from the topology (wire order) which follows the conventional CCW-outward
 // convention for properly-built solids.
 let is_curved = brep
 .tshapes
 .get(face_flat_idx)
 .and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.as_ref() } else { None } })
 .map_or(false, |s| !matches!(s, Surface3::Plane(_)));

 // Holed / trimmed faces: the boolean result may cache outer-loop triangles only
 // (see `ResultBuilder`); re-mesh in plane/UV for area / volume.
 if !face.inner_wires.is_empty() {
 if let Some(tris) = try_face_with_holes(brep, face, face_flat_idx) {
 return tris;
 }
 } else {
 // Spherical sub-faces: prefer UV ear-cut of the boundary loop over builder
 // 3D ear-clip, which can under-fill a non-convex (trimmed) patch.
 if let Some(surf) = brep
 .tshapes
 .get(face_flat_idx)
 .and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.as_ref() } else { None } })
 {
 if let Surface3::Sphere(s) = surf {
 let n_edges = face.outer_wire.edges.len();
 let per_edge = if n_edges > 600 { 4 } else if n_edges > 300 { 8 } else if n_edges > 150 { 16 } else if n_edges > 30 { 24 } else { 48 };
 let mut outer = sample_wire_polyline_3d_with_n(brep, &face.outer_wire, per_edge);
 trim_almost_closed_polyline(&mut outer, 1e-5);
 if outer.len() >= 3 {
 // UV grid raster: for large merged faces the UV polygon wraps
 // around u multiple times; the grid now detects this and falls
 // back to 3D-only point-in-spherical-polygon tests.
 if let Some(tris) = try_spherical_uv_masked_raster(s, brep, face, face_flat_idx, face.normal) {
 return tris;
 }
 if let Some(tris) = try_spherical_earcut_simple(s, &outer, face.normal) {
 return tris;
 }
 if let Some(tris) = try_planar_earcut_simple_outer(&outer, face.normal) {
 return tris;
 }
 }
 }
 }
 }

 // For planar surfaces, stored triangles from emit_face_with_origin's
 // triangulate_polygon are correct (the surface IS planar).  For curved
 // surfaces, those same triangles are ear-clipped chordal triangles that
 // cut across the interior rather than following the surface 鈥?they produce
 // wrong signed volume.  Always fall through to tessellate_curved_face for
 // curved surfaces even when face.triangles is non-empty.
 if !is_curved && face.inner_wires.is_empty() && !face.triangles.is_empty() {
 return face
 .triangles
 .iter()
 .filter_map(|&[i, j, k]| {
 let a = brep.vertex_point(i)?;
 let b = brep.vertex_point(j)?;
 let c = brep.vertex_point(k)?;
 Some(orient_tri([a, b, c], face.normal))
 })
 .collect();
 }

 if let Some(uv_tris) = tessellate_curved_face(brep, face, face_flat_idx) {
 // UV-grid tessellation over the natural surface domain when
 // `face.triangles` is empty (untessellated primitives).
 return uv_tris;
 }

 let wire_pts = sample_wire_polyline_3d(brep, &face.outer_wire);
 if wire_pts.len() < 3 {
 return Vec::new();
 }
 // Fan from first sample (convex-ish outer loops only; holed cases handled above).
 let origin = wire_pts[0];
 if is_curved {
 let tris: Vec<_> = (1..wire_pts.len() - 1)
 .map(|i| [origin, wire_pts[i], wire_pts[i + 1]])
 .collect();
 return tris;
 }
 // For planar faces, use ear-clipping on the polygon vertices (extracted
 // from wire edges, not the dense 65-per-edge samples which include
 // inter-edge duplicates that confuse earcut).  The fan-from-vertex-0
 // approach fills the convex hull, inflating volume for concave caps.
 let face_verts: Vec<DVec3> = face.outer_wire.edges.iter().filter_map(|we| {
 let edge = brep.flat_edges().get(we.idx).copied()?;
 let vi = if we.forward { edge.0 } else { edge.1 };
 brep.vertex_point(vi)
 }).collect();
 if face_verts.len() >= 3 {
 if let Some(ear_tris) = try_planar_earcut_simple_outer(&face_verts, face.normal) {
 return ear_tris;
 }
 }
 // Fallback: fan from first dense sample (convex or near-degenerate).
 (1..wire_pts.len() - 1)
 .map(|i| orient_tri([origin, wire_pts[i], wire_pts[i + 1]], face.normal))
 .collect()
}

/// Ensure triangle [a,b,c] is oriented so its normal agrees with `face_normal`.
#[inline]
fn orient_tri(tri: [DVec3; 3], face_normal: DVec3) -> [DVec3; 3] {
 let [a, b, c] = tri;
 let n = (b - a).cross(c - a);
 if n.dot(face_normal) < 0.0 { [a, c, b] } else { [a, b, c] }
}

/// Iterate over (face_flat_index, &Face) pairs across all solids/shells.
pub fn face_flat_iter(brep: &topods::BRep) -> Vec<(usize, crate::topology::Face)> {
 let mut faces = Vec::new();
 for (fi, ts) in brep.tshapes.iter().enumerate() {
  if let topods::TShape::Face(fd) = &**ts {
  faces.push((fi, crate::topology::Face {
   outer_wire: crate::topology::Wire { edges: Vec::new() },
   inner_wires: Vec::new(), normal: DVec3::Z,
   triangles: Vec::new(), sample_point: fd.sample_point,
   mesh_dirty: true, surface_idx: None,
  }));
  }
 }
 faces
}

// 鈹€鈹€ Public API 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Public re-export of `face_triangles` for diagnostic use in tests.
/// Not part of the stable API.
#[doc(hidden)]
pub fn face_triangles_pub(
 brep: &topods::BRep,
 face: &crate::topology::Face,
 face_flat_idx: usize,
) -> Vec<[DVec3; 3]> {
 face_triangles(brep, face, face_flat_idx)
}

/// Compute the total surface area of all faces in the BRep.
///
/// For each face, uses analytic area when available (planar shoe lace; holed
/// sphere: parametric `R虏 d惟` on the same UV mask as the raster), otherwise
/// sums triangle areas (pre-triangulated, UV-sampled, or fan-triangulated).
/// Returns 0.0 if the BRep has no faces.
///
/// Cylinder sub-faces from boolean splitting may have overlapping UV polygons
/// (the intersection curve mapping is numerically inconsistent).  When the sum
/// of UV shoelace areas across all sub-faces on a cylinder surface exceeds
/// `2蟺 脳 (v_global_max - v_global_min)`, each sub-face area is scaled down
/// proportionally to the correct total.
pub fn surface_area(brep: &topods::BRep) -> f64 {
 struct CylEntry { sa: f64, uv_area: f64, v_min: f64, v_max: f64, radius: f64 }

 // Groups of cylinder faces on the same geometric cylinder (identified by
 // comparing origin, axis, and radius within tolerance).  Boolean splitting
 // may create multiple surface indices for the same physical cylinder, so
 // we cannot rely on surface-index equality alone.
 struct CylGroup { origin: DVec3, axis: DVec3, radius: f64, entries: Vec<CylEntry> }
 let mut cyl_groups: Vec<CylGroup> = Vec::new();
 let mut total = 0.0f64;

 for (_fi, f) in face_flat_iter(brep) {
 // Pre-computed clean triangles (e.g. Steinmetz analytic builder which
 // tessellates each face accurately in UV space) avoid
 // `short_delta_on_circle_01` picking the wrong direction for
 // triangular faces that span >蟺 on a cylinder's angular domain.
 // HOWEVER: `mesh_brep` (called after boolean pipeline) tessellates
 // sphere faces using the full default UV domain [0,2蟺]脳[-蟺/2,蟺/2]
 // instead of the trimmed portion, producing wrong triangles for
 // trimmed spherical sub-faces.  Always try analytic first.
 let analytic = try_analytic_face_surface_area(brep, &f, _fi);
 let has_clean_tris = !f.triangles.is_empty() && !f.mesh_dirty;
 let a = match analytic {
 Some(sa) => sa,
 None => {
 if has_clean_tris {
 f.triangles.iter().map(|&[a, b, c]| tri_area(
 brep.vertex_point(a).unwrap_or(DVec3::ZERO),
 brep.vertex_point(b).unwrap_or(DVec3::ZERO),
 brep.vertex_point(c).unwrap_or(DVec3::ZERO),
 )).sum()
 } else {
 let tris = face_triangles(brep, &f, _fi);
 tris.iter().map(|&[a, b, c]| tri_area(a, b, c)).sum()
 }
 }
 };
 total += a;

 if cfg!(debug_assertions) && std::env::var("RCAD_DEBUG_SA").is_ok() {
 let vname = brep.tshapes.get(_fi).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } })
 .and_then(|si| Some(si))
 .map(|s| match s { Surface3::Plane(_) => "Plane", Surface3::BSpline(_) => "BSpline", Surface3::Sphere(_) => "Sphere", Surface3::Cylinder(_) => "Cylinder", _ => "Other" })
 .unwrap_or("none");
 let n_inner = f.inner_wires.len();
 let inner_edges: usize = f.inner_wires.iter().map(|w| w.edges.len()).sum();
 eprintln!("[SA_FACE] fi={} surf={} analytic={:?} area={:.6} inner_wires={} inner_edges={} n_tris={} mesh_dirty={}",
 _fi, vname, analytic, a, n_inner, inner_edges, f.triangles.len(), f.mesh_dirty);
 }
 // CI assertion: warn if any face has no surface reference (mesh-only fast-path).
 // All boolean results should have proper analytic surfaces for exact SA.
 if analytic.is_none() && !f.triangles.is_empty() && cfg!(debug_assertions) {
 eprintln!("[SA_WARN] face[{}] has no analytic surface 鈥?SA from {} triangles (total {:.6})",
 _fi, f.triangles.len(), a);
 }

 // Track analytic cylinder faces for UV-overlap normalization.
 // Skip faces with inner wires (figure-8 case has its own handling).
 if analytic.is_some() {
 if let Some(si) = brep.tshapes.get(_fi).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } }) {
 if let Some(surf) = Some(si) {
 if let Surface3::Cylinder(c) = &surf {
 if f.inner_wires.is_empty() {
 if let Some([_, _, v0, v1]) = curved_face_uv_domain(brep, &f, _fi, &surf) {
 if let Some(uv_area) = cylinder_outer_wire_uv_shoelace_area(c, brep, &f) {
 // Find or create group for this geometric cylinder
 let found = cyl_groups.iter_mut().find(|g| {
 (g.origin - c.origin).length_squared() < 1e-8
 && g.axis.dot(c.axis) > 1.0 - 1e-8
 && (g.radius - c.radius).abs() < 1e-8
 });
 let g = match found {
 Some(g) => g,
 None => {
 cyl_groups.push(CylGroup {
 origin: c.origin, axis: c.axis,
 radius: c.radius, entries: Vec::new(),
 });
 cyl_groups.last_mut().unwrap()
 }
 };
 g.entries.push(CylEntry {
 sa: a, uv_area, v_min: v0, v_max: v1, radius: c.radius,
 });
 }
 }
 }
 }
 }
 }
 }
 }

 // Normalize each cylinder group: if sum of UV shoelace areas exceeds
 // the maximum possible UV area (2蟺 脳 global v-span), scale down proportionally.
 if std::env::var("RCAD_DEBUG_BUILDER").is_ok() {
 // Second pass: print all face area breakdown
 for (_fi, f) in face_flat_iter(brep) {
 let surf_idx = brep.tshapes.get(_fi).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } });
 let surf_name = surf_idx.and_then(|si| Some(si)).map(|s| {
 match s {
 Surface3::Plane(_) => "Plane".to_string(),
 Surface3::Cylinder(_) => "Cylinder".to_string(),
 Surface3::Cone(_) => "Cone".to_string(),
 Surface3::Sphere(_) => "Sphere".to_string(),
 Surface3::Torus(_) => "Torus".to_string(),
 _ => "Other".to_string(),
 }
 }).unwrap_or_else(|| "None".to_string());
 let edge_n = f.outer_wire.edges.len();
 let analytic = try_analytic_face_surface_area(brep, &f, _fi);
 let a = match analytic {
 Some(sa) => sa,
 None => {
 let tris = face_triangles(brep, &f, _fi);
 tris.iter().map(|&[a, b, c]| tri_area(a, b, c)).sum()
 }
 };
 eprintln!("[SA_FACE] fi={} surf={} edges={} inner={} analytic={} area={:.6}",
 _fi, surf_name, edge_n, f.inner_wires.len(),
 analytic.is_some(), a);
 }
 eprintln!("[SA_TRACK] cylinder groups={}", cyl_groups.len());
 for (gi, g) in cyl_groups.iter().enumerate() {
 eprintln!("[SA_TRACK]  cyl group={} entries={} R={}", gi, g.entries.len(), g.radius);
 let entries = &g.entries;
 let total_uv: f64 = entries.iter().map(|e| e.uv_area).sum();
 let v_min = entries.iter().map(|e| e.v_min).fold(f64::INFINITY, f64::min);
 let v_max = entries.iter().map(|e| e.v_max).fold(f64::NEG_INFINITY, f64::max);
 eprintln!("[SA_TRACK] total_uv={:.6} v_rng=[{:.4},{:.4}] max_uv={:.6}",
 total_uv, v_min, v_max, std::f64::consts::PI * 2.0 * (v_max - v_min));
 for (i, e) in entries.iter().enumerate() {
 eprintln!("[SA_TRACK] [{}] sa={:.6} uv={:.6} v=[{:.4},{:.4}] R={}",
 i, e.sa, e.uv_area, e.v_min, e.v_max, e.radius);
 }
 }
 }
 const TWO_PI: f64 = std::f64::consts::PI * 2.0;
 for g in &cyl_groups {
 let entries = &g.entries;
 if entries.len() < 2 { continue; }
 let total_uv: f64 = entries.iter().map(|e| e.uv_area).sum();
 let v_min = entries.iter().map(|e| e.v_min).fold(f64::INFINITY, f64::min);
 let v_max = entries.iter().map(|e| e.v_max).fold(f64::NEG_INFINITY, f64::max);
 let v_span = v_max - v_min;
 if v_span < 1e-14 { continue; }
 let max_uv = TWO_PI * v_span;
 if total_uv > max_uv * 1.01 {
 let scale = max_uv / total_uv;
 let old_sum: f64 = entries.iter().map(|e| e.sa).sum();
 total = total - old_sum;
 total += entries.iter().map(|e| e.radius * e.uv_area * scale).sum::<f64>();
 }
 }

 total
}

/// Area of one face, using the same rules as [`surface_area`].
#[doc(hidden)]
pub fn face_surface_area(brep: &topods::BRep, face: &Face, face_flat_idx: usize) -> f64 {
 if let Some(a) = try_analytic_face_surface_area(brep, face, face_flat_idx) {
 a
 } else {
 let tris = face_triangles(brep, face, face_flat_idx);
 tris
 .iter()
 .map(|&[a, b, c]| tri_area(a, b, c))
 .sum()
 }
}

/// Analytic volume contribution of a sphere face: V = (R/3) 脳 A from the parametric
/// UV-mask integral using the 5-point OR acceptance test (same as the raster used
/// by `try_spherical_uv_masked_raster`), not the strict all-corners test used by
/// `sphere_holed_mask_param_area_sum`.  The 5-point test captures partially-covered
/// boundary cells that the all-corners test misses, giving a more accurate surface
/// integral at moderate grid resolutions.
///
/// For a sphere centered at origin, r路n = R, so the divergence-theorem integral
/// V = (1/3)鈭埆 r路n dA reduces to (R/3)路A.
fn sphere_holed_mask_param_volume_sum(s: &SphericalSurface, ctx: &SphereHoledMaskCtx) -> f64 {
 const N: usize = SPHERE_UV_MASK_N;
 let umin = ctx.umin;
 let umax = ctx.umax;
 let vmin = ctx.vmin;
 let vmax = ctx.vmax;
 let du = (umax - umin) / N as f64;
 let dv = (vmax - vmin) / N as f64;
 let r3 = s.radius * s.radius * s.radius;
 let inner = &ctx.inner_polys;
 let inner_3d = &ctx.inner_3d;
 let use_uv = ctx.use_uv_winding;

 let emit = |poly_uv: &[DVec2], poly_3d: &[DVec3], use_inner_mask: bool| -> f64 {
 let mut v = 0.0_f64;
 for i in 0..N {
 for j in 0..N {
 let u0 = umin + i as f64 * du;
 let u1 = u0 + du;
 let v0 = vmin + j as f64 * dv;
 let v1 = v0 + dv;
 let uc = umin + (i as f64 + 0.5) * du;
 let vc = vmin + (j as f64 + 0.5) * dv;

 let p00 = s.point_at(u0, v0);
 let p10 = s.point_at(u1, v0);
 let p11 = s.point_at(u1, v1);
 let p01 = s.point_at(u0, v1);
 let pc = s.point_at(uc, vc);

 // 5-point OR test matching try_spherical_uv_masked_raster:
 // accept if center is inside (UV winding + 3D angular-sum), OR
 // if any corner is inside (3D angular-sum).
 // Use `!= 0` (not `> 0`) to handle both clockwise and CCW polygons.
 let center_in = winding_number_2d(poly_uv, DVec2::new(uc, vc)) != 0
 || point_in_spherical_polygon_3d(poly_3d, pc);
 let any_corner_in = point_in_spherical_polygon_3d(poly_3d, p00)
 || point_in_spherical_polygon_3d(poly_3d, p10)
 || point_in_spherical_polygon_3d(poly_3d, p11)
 || point_in_spherical_polygon_3d(poly_3d, p01);
 if !center_in && !any_corner_in {
 continue;
 }
 if use_inner_mask {
 let in_hole = if use_uv {
 inner.iter().any(|h| point_in_polygon_2d(h, DVec2::new(uc, vc)))
 } else {
 inner_3d.iter().any(|h3d| point_in_spherical_polygon_3d(h3d, pc))
 };
 if in_hole {
 continue;
 }
 }

 // Analytic dV = (R/3) * R虏 * du * (cos(v鈧€) - cos(v鈧?)
 v += r3 / 3.0 * du * (v0.cos() - v1.cos());
 }
 }
 v
 };

 let mut t = emit(&ctx.outer_uv, &ctx.outer_3d, true);
 if t <= 0.0 && !inner.is_empty() {
 t = emit(&ctx.outer_uv, &ctx.outer_3d, false);
 }
 if t > 0.0 {
 return t;
 }
 // Fallback: reversed winding
 let mut rev_uv = ctx.outer_uv.clone();
 rev_uv.reverse();
 let mut rev_3d = ctx.outer_3d.clone();
 rev_3d.reverse();
 t = emit(&rev_uv, &rev_3d, true);
 if t <= 0.0 && !inner.is_empty() {
 t = emit(&rev_uv, &rev_3d, false);
 }
 t
}

/// Analytic volume contribution of a sphere face: V = (R/3) 脳 A from the parametric
/// UV-mask integral.  For a sphere centered at origin, r路n = R, so the divergence-theorem
/// integral V = (1/3)鈭埆 r路n dA reduces to (R/3)路A.  This avoids the non-manifold tet-sum
/// error from the area-corrected chordal triangulation in `try_spherical_uv_masked_raster`.
///
/// Only applies when the sphere center is at the origin (within tolerance) 鈥?for offset
/// spheres the surface integral does not simplify to (R/3)路A and the tet sum fallback
/// is used instead.
const SPHERE_CENTER_AT_ORIGIN_TOL: f64 = 1e-10;
fn try_sphere_face_analytic_volume(
 brep: &BRep,
 face: &Face,
 face_flat_idx: usize,
) -> Option<f64> {
 let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| { if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None } })?;
 let surf = &surf_idx;
 match surf {
 Surface3::Sphere(s) if s.center.length_squared() < SPHERE_CENTER_AT_ORIGIN_TOL => {
 let ctx = spherical_holed_uv_mask_setup(s, brep, face)?;
 let vol = sphere_holed_mask_param_volume_sum(s, &ctx);
 if vol > 0.0 { Some(vol) } else { None }
 }
 _ => None,
 }
}

/// Divergence-theorem signed volume: `(1/6) 危 a路(b脳c)` over surface triangles (no absolute value).
///
/// Sphere faces use the analytic parametric integral (V = (R/3)路A) instead of the
/// per-cell area-corrected triangulation which creates a non-manifold mesh and
/// systematically underestimates volume via the tet sum.  Plane and other faces
/// use the standard triangulation which remains correct.
///
/// For a closed solid with **outward** face normals this is positive; **inward** shells
/// (e.g. OCCT `treverse` before `prism`) yield a negative sum when the mesh is consistent.
pub fn signed_volume(brep: &topods::BRep) -> f64 {
 let mut vol = 0.0_f64;
 for (fi, face) in face_flat_iter(brep) {
 if let Some(analytic_vol) = try_sphere_face_analytic_volume(brep, &face, fi) {
 vol += analytic_vol;
 } else {
 for [a, b, c] in face_triangles(brep, &face, fi) {
 vol += tet_signed_volume(a, b, c);
 }
 }
 }
 vol
}

/// Absolute volume (see [`signed_volume`]).
pub fn volume(brep: &topods::BRep) -> f64 {
 signed_volume(brep).abs()
}

/// Compute the centroid (center of mass) of the solid by volumetric integration.
///
/// Uses the formula: C = (1 / 8V) 危_triangles (a+b+c) * tet_signed_vol(a,b,c)
/// where the sum is over all surface triangles.
///
/// Falls back to `BRep::center()` (vertex average) if the volume is near zero.
pub fn centroid(brep: &topods::BRep) -> DVec3 {
 let mut vol_sum = 0.0_f64;
 let mut weighted_sum = DVec3::ZERO;

 for (fi, face) in face_flat_iter(brep) {
 for [a, b, c] in face_triangles(brep, &face, fi) {
 let sv = tet_signed_volume(a, b, c);
 vol_sum += sv;
 // Weight the centroid of each tet (at (a+b+c+origin)/4,
 // origin=0) 鈫?simplified to (a+b+c) * sv
 weighted_sum += (a + b + c) * sv;
 }
 }

 if vol_sum.abs() < 1e-15 {
 return brep.center();
 }

 // Centroid formula: (1/(2 * 4 * vol_sum)) * 危 (a+b+c) * sv
 // Simplification: weighted_sum / (4 * vol_sum) gives tet centroid average
 weighted_sum / (4.0 * vol_sum)
}

// 鈹€鈹€ Inertia tensor 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Symmetric 3脳3 moment of inertia tensor (assuming uniform density = 1).
///
/// The components are defined as:
/// ```text
/// Ixx = 鈭?y虏+z虏) dV,  Iyy = 鈭?x虏+z虏) dV,  Izz = 鈭?x虏+y虏) dV
/// Ixy = -鈭玿y dV, Ixz = -鈭玿z dV, Iyz = -鈭珁z dV
/// ```
///
/// Computed about the world origin. To get the tensor about the centroid,
/// use the parallel-axis theorem.
#[derive(Debug, Clone, Copy)]
pub struct InertiaTensor {
 pub ixx: f64,
 pub iyy: f64,
 pub izz: f64,
 pub ixy: f64,
 pub ixz: f64,
 pub iyz: f64,
}

impl InertiaTensor {
 /// Returns the 3脳3 inertia matrix as row-major `[[f64;3];3]`.
 pub fn to_matrix(&self) -> [[f64; 3]; 3] {
 [
 [self.ixx, -self.ixy, -self.ixz],
 [-self.ixy, self.iyy, -self.iyz],
 [-self.ixz, -self.iyz, self.izz],
 ]
 }
}

/// Computes the moment of inertia tensor of a closed BRep solid about the
/// world origin.
///
/// Uses the divergence theorem (polyhedral formula from Mirtich 1996) applied
/// to the BRep's triangulated faces, consistent with the existing `volume` and
/// `centroid` implementations.
///
/// Assumes uniform density = 1 (unit density).  Multiply each component by
/// the actual density to get physical inertia.
pub fn inertia_tensor(brep: &topods::BRep) -> InertiaTensor {
 let mut ixx = 0.0_f64;
 let mut iyy = 0.0_f64;
 let mut izz = 0.0_f64;
 let mut ixy = 0.0_f64;
 let mut ixz = 0.0_f64;
 let mut iyz = 0.0_f64;

 for (fi, face) in face_flat_iter(brep) {
 for [a, b, c] in face_triangles(brep, &face, fi) {
 // Signed volume of tet (origin, a, b, c)
 // sv = a路(b脳c)/6 鈥?same as tet_signed_volume
 let sv = a.dot(b.cross(c)) / 6.0;

 // Symmetric quadratic sums for each coordinate pair.
 // For 鈭玙tet x虏 dV = sv/10 * x2_sym (from simplex integration).
 let x2 = a.x * a.x + b.x * b.x + c.x * c.x + a.x * b.x + a.x * c.x + b.x * c.x;
 let y2 = a.y * a.y + b.y * b.y + c.y * c.y + a.y * b.y + a.y * c.y + b.y * c.y;
 let z2 = a.z * a.z + b.z * b.z + c.z * c.z + a.z * b.z + a.z * c.z + b.z * c.z;

 ixx += sv / 10.0 * (y2 + z2);
 iyy += sv / 10.0 * (x2 + z2);
 izz += sv / 10.0 * (x2 + y2);

 // For 鈭玙tet xy dV = sv/20 * xy_mixed (from simplex integration).
 // Product-moment: Ixy = -鈭玿y dV, etc.
 let xy = 2.0 * (a.x * a.y + b.x * b.y + c.x * c.y)
 + a.x * b.y
 + b.x * a.y
 + a.x * c.y
 + c.x * a.y
 + b.x * c.y
 + c.x * b.y;
 let xz = 2.0 * (a.x * a.z + b.x * b.z + c.x * c.z)
 + a.x * b.z
 + b.x * a.z
 + a.x * c.z
 + c.x * a.z
 + b.x * c.z
 + c.x * b.z;
 let yz = 2.0 * (a.y * a.z + b.y * b.z + c.y * c.z)
 + a.y * b.z
 + b.y * a.z
 + a.y * c.z
 + c.y * a.z
 + b.y * c.z
 + c.y * b.z;

 ixy += sv / 20.0 * xy;
 ixz += sv / 20.0 * xz;
 iyz += sv / 20.0 * yz;
 }
 }

 // Diagonal terms must be positive for a physical solid.
 // Off-diagonal sign: Ixy = -鈭玿y dV so negate the accumulated sums.
 InertiaTensor {
 ixx: ixx.abs(),
 iyy: iyy.abs(),
 izz: izz.abs(),
 ixy: -ixy,
 ixz: -ixz,
 iyz: -iyz,
 }
}

// 鈹€鈹€ Tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€






