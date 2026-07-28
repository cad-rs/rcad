//! Triangulation helpers for GProp surface/volume computation.
//!
//! OCCT: internal helpers used by BRepGProp. Contains triangle decomposition,
//! UV unwrapping, ear-cutting, and winding-number utilities.

use glam::{DVec2, DVec3};
use std::f64::consts::PI;

use crate::BRep;
use crate::geom::{Curve3, CurveEval, SphericalSurface, Surface3, SurfaceEval};
use crate::topo::topods;
use crate::topo::topology::{Face, Wire, WireEdge};

/// Compute the signed area of a triangle from three points.
#[inline]
pub fn tri_area(a: DVec3, b: DVec3, c: DVec3) -> f64 {
    (b - a).cross(c - a).length() * 0.5
}

/// Signed volume contribution of a tetrahedron from the origin.
#[inline]
pub fn tet_signed_volume(a: DVec3, b: DVec3, c: DVec3) -> f64 {
    a.dot(b.cross(c)) / 6.0
}

/// Sample a closed wire to a 3D polyline.
pub fn sample_wire_polyline_3d(brep: &topods::BRep, wire: &Wire) -> Vec<DVec3> {
    sample_wire_polyline_3d_with_n(brep, wire, 1024)
}

/// Like sample_wire_polyline_3d with configurable samples per edge.
pub fn sample_wire_polyline_3d_with_n(brep: &topods::BRep, wire: &Wire, n: usize) -> Vec<DVec3> {
    let mut pts = Vec::new();
    for we in &wire.edges {
        let flat_edges = brep.flat_edges();
        let edge = match flat_edges.get(we.idx) {
            Some(e) => e,
            None => continue,
        };
        let curve_opt = brep.tshapes.get(we.idx).and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts {
                ed.curve.as_ref()
            } else { None }
        });
        if let Some(curve) = curve_opt {
            let range = brep.tshapes.get(we.idx).and_then(|ts| {
                if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None }
            }).unwrap_or_else(|| curve.default_domain());
            let [t0, t1] = if we.forward { range } else { [range[1], range[0]] };
            if (t1 - t0).abs() > 1e-12 && t0.is_finite() && t1.is_finite() {
                let full_circle = (t1 - t0).abs() >= 2.0 * PI - 1e-9;
                let samples = if full_circle { n } else { n + 1 };
                for k in 0..samples {
                    let frac = k as f64 / n as f64;
                    pts.push(curve.point_at(t0 + (t1 - t0) * frac));
                }
                continue;
            }
        }
        let vidx = if we.forward { edge.0 } else { edge.1 };
        if let Some(v) = brep.vertex_point(vidx) { pts.push(v); }
    }
    pts
}

pub fn trim_almost_closed_polyline(pts: &mut Vec<DVec3>, tol: f64) {
    if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length() < tol {
        pts.pop();
    }
}

pub fn local_basis_from_normal(normal: DVec3) -> (DVec3, DVec3) {
    let ref_dir = if normal.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let u = normal.cross(ref_dir).normalize();
    let v = normal.cross(u).normalize();
    (u, v)
}

pub fn sphere_point_to_uv(s: &SphericalSurface, p: DVec3) -> DVec2 {
    (*s).world_to_uv(p)
}

pub fn unwrap_sphere_u_in_chain(uvs: &mut [DVec2]) {
    if uvs.len() < 2 { return; }
    let two_pi = 2.0 * PI;
    let mut prev = uvs[0].x;
    for i in 1..uvs.len() {
        let mut u = uvs[i].x;
        let mut d = u - prev;
        while d > PI { u -= two_pi; d = u - prev; }
        while d < -PI { u += two_pi; d = u - prev; }
        uvs[i].x = u;
        prev = u;
    }
}

// ── Earcut triangulation helpers ─────────────────────────────────────────

pub fn earcut_indices_from_flat(flat: &[f64], hole_indices: &[usize]) -> Vec<usize> {
    let coords: Vec<[f64; 2]> = flat.chunks_exact(2).map(|c| [c[0], c[1]]).collect();
    if coords.len() < 3 { return Vec::new(); }
    let mut out = Vec::new();
    let mut ear = earcut::Earcut::new();
    ear.earcut(coords, hole_indices, &mut out);
    out
}

pub fn earcut_flat_to_tris(
    flat: &[f64], hole_starts: &[usize], all_3d: &[DVec3], face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
    if flat.len() < 9 || hole_starts.is_empty() { return None; }
    let indices = earcut_indices_from_flat(flat, hole_starts);
    if indices.is_empty() { return None; }
    let mut out = Vec::with_capacity(indices.len() / 3);
    for tri in indices.chunks_exact(3) {
        let a = all_3d[tri[0]]; let b = all_3d[tri[1]]; let c = all_3d[tri[2]];
        out.push(orient_tri([a, b, c], face_normal));
    }
    Some(out)
}

pub fn try_planar_earcut_holes(
    outer: &[DVec3], holes: &[Vec<DVec3>], face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
    let (ux, uy) = local_basis_from_normal(face_normal);
    let pivot = outer.first().copied()?;
    let mut all_3d = Vec::new();
    let mut flat = Vec::new();
    for p in outer { let q = *p - pivot; flat.push(q.dot(ux)); flat.push(q.dot(uy)); all_3d.push(*p); }
    let mut hole_starts = Vec::new();
    for h in holes {
        if h.len() < 3 { continue; }
        hole_starts.push(all_3d.len());
        for p in h { let q = *p - pivot; flat.push(q.dot(ux)); flat.push(q.dot(uy)); all_3d.push(*p); }
    }
    if hole_starts.is_empty() { return None; }
    earcut_flat_to_tris(&flat, &hole_starts, &all_3d, face_normal)
}

pub fn try_planar_earcut_simple_outer(outer: &[DVec3], face_normal: DVec3) -> Option<Vec<[DVec3; 3]>> {
    if outer.len() < 3 { return None; }
    let (ux, uy) = local_basis_from_normal(face_normal);
    let pivot = outer[0];
    let mut flat = Vec::with_capacity(2 * outer.len());
    for p in outer { let q = *p - pivot; flat.push(q.dot(ux)); flat.push(q.dot(uy)); }
    let all_3d = outer.to_vec();
    let indices = earcut_indices_from_flat(&flat, &[]);
    if indices.is_empty() { return None; }
    let expected_min = (outer.len() - 2).saturating_sub(outer.len() / 4);
    if indices.len() / 3 < expected_min { return None; }
    let mut out = Vec::with_capacity(indices.len() / 3);
    for tri in indices.chunks_exact(3) {
        out.push(orient_tri([all_3d[tri[0]], all_3d[tri[1]], all_3d[tri[2]]], face_normal));
    }
    Some(out)
}

// ── Spherical earcut helpers ─────────────────────────────────────────────

pub fn try_spherical_earcut_simple(
    s: &SphericalSurface, outer: &[DVec3], face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
    if outer.len() < 3 { return None; }
    let mut outer_uv: Vec<DVec2> = outer.iter().map(|p| sphere_point_to_uv(s, *p)).collect();
    let u_vals: Vec<f64> = outer_uv.iter().map(|uv| uv.x).collect();
    let unwrapped_u = unwrap_u_circle_chain_closed(&u_vals);
    for (uv, u) in outer_uv.iter_mut().zip(unwrapped_u.iter()) { uv.x = *u; }
    const MAX_UV_PTS: usize = 8000;
    if outer_uv.len() > MAX_UV_PTS {
        let step = (outer_uv.len() - 1) as f64 / (MAX_UV_PTS - 1) as f64;
        outer_uv = (0..MAX_UV_PTS).map(|i| {
            let idx = (i as f64 * step).round() as usize;
            outer_uv[idx.min(outer_uv.len() - 1)]
        }).collect();
    }
    outer_uv.dedup_by(|a, b| (*a - *b).length_squared() < 1e-8);
    if outer_uv.len() < 3 { return None; }
    let mut flat = Vec::with_capacity(2 * outer_uv.len());
    for uv in &outer_uv { flat.push(uv.x); flat.push(uv.y); }
    let deduped_3d: Vec<DVec3> = outer_uv.iter().map(|uv| s.point_at(uv.x, uv.y)).collect();
    let indices = earcut_indices_from_flat(&flat, &[]);
    if indices.is_empty() { return None; }
    let expected_min = (outer_uv.len() - 2).saturating_sub(outer_uv.len() / 4);
    if indices.len() / 3 < expected_min { return None; }
    let mut out = Vec::with_capacity(indices.len() / 3);
    for tri in indices.chunks_exact(3) {
        out.push(orient_tri([deduped_3d[tri[0]], deduped_3d[tri[1]], deduped_3d[tri[2]]], face_normal));
    }
    Some(out)
}

pub fn try_spherical_earcut_holes(
    s: &SphericalSurface, outer: &[DVec3], holes: &[Vec<DVec3>], face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
    let mut outer_uv: Vec<DVec2> = outer.iter().map(|p| sphere_point_to_uv(s, *p)).collect();
    unwrap_sphere_u_in_chain(&mut outer_uv);
    let mut all_3d = outer.to_vec();
    let mut flat = Vec::with_capacity(2 * (outer.len() + holes.iter().map(|h| h.len()).sum::<usize>()));
    for uv in &outer_uv { flat.push(uv.x); flat.push(uv.y); }
    let mut hole_starts = Vec::new();
    for h in holes {
        if h.len() < 3 { continue; }
        hole_starts.push(all_3d.len());
        let mut huv = h.iter().map(|p| sphere_point_to_uv(s, *p)).collect::<Vec<_>>();
        unwrap_sphere_u_in_chain(&mut huv);
        for p in h { all_3d.push(*p); }
        for uv in huv { flat.push(uv.x); flat.push(uv.y); }
    }
    if hole_starts.is_empty() { return None; }
    earcut_flat_to_tris(&flat, &hole_starts, &all_3d, face_normal)
}

pub fn try_face_with_holes(
    brep: &BRep, face: &Face, face_flat_idx: usize,
) -> Option<Vec<[DVec3; 3]>> {
    if face.inner_wires.is_empty() { return None; }
    let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
    })?;
    let surf = &surf_idx;
    let tol = 1e-5;
    let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
    trim_almost_closed_polyline(&mut outer, tol);
    if outer.len() < 3 { return None; }
    let mut holes_3d = Vec::new();
    for iw in &face.inner_wires {
        let mut h = sample_wire_polyline_3d(brep, iw);
        trim_almost_closed_polyline(&mut h, tol);
        if h.len() >= 3 { holes_3d.push(h); }
    }
    if holes_3d.is_empty() { return None; }
    let rev_holes: Vec<Vec<DVec3>> = holes_3d.iter().map(|h| h.iter().rev().copied().collect()).collect();
    match surf {
        Surface3::Plane(_) =>
            try_planar_earcut_holes(&outer, &holes_3d, face.normal)
            .or_else(|| try_planar_earcut_holes(&outer, &rev_holes, face.normal)),
        Surface3::Sphere(_s) => {
            try_spherical_earcut_holes(_s, &outer, &holes_3d, face.normal)
                .or_else(|| try_spherical_earcut_holes(_s, &outer, &rev_holes, face.normal))
                .or_else(|| try_planar_earcut_holes(&outer, &holes_3d, face.normal))
                .or_else(|| try_planar_earcut_holes(&outer, &rev_holes, face.normal))
        }
        _ => None,
    }
}

// ── 2D geometry predicates ───────────────────────────────────────────────

pub fn winding_number_2d(poly: &[DVec2], pt: DVec2) -> i32 {
    let n = poly.len();
    if n < 3 { return 0; }
    let mut wn = 0i32;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i]; let vj = poly[j];
        if vi.y <= pt.y {
            if vj.y > pt.y && is_left_2d(vi, vj, pt) > 0.0 { wn += 1; }
        } else if vj.y <= pt.y && is_left_2d(vi, vj, pt) < 0.0 { wn -= 1; }
        j = i;
    }
    wn
}

pub fn is_left_2d(a: DVec2, b: DVec2, c: DVec2) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

pub fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
    let n = poly.len();
    if n < 3 { return false; }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i]; let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x) { inside = !inside; }
        j = i;
    }
    inside
}

pub fn point_in_spherical_polygon_3d(boundary: &[DVec3], point: DVec3) -> bool {
    let n = boundary.len();
    if n < 3 { return false; }
    let mut total = 0.0;
    for i in 0..n {
        let j = if i + 1 < n { i + 1 } else { 0 };
        let a = boundary[i] - point;
        let b = boundary[j] - point;
        let cross = a.cross(b);
        let theta = cross.length().atan2(a.dot(b));
        if cross.dot(point) >= 0.0 { total += theta; } else { total -= theta; }
    }
    total.abs() > PI
}

pub fn point_in_spherical_polygon_3d_pub(boundary: &[DVec3], point: DVec3) -> bool {
    point_in_spherical_polygon_3d(boundary, point)
}

pub fn short_delta_on_circle_01(a: f64, b: f64) -> f64 {
    const TWO_PI: f64 = 2.0 * PI;
    let am = a.rem_euclid(TWO_PI); let bm = b.rem_euclid(TWO_PI);
    let mut d = bm - am;
    if d > PI { d -= TWO_PI; } if d < -PI { d += TWO_PI; }
    d
}

pub fn unwrap_u_circle_chain_closed(b: &[f64]) -> Vec<f64> {
    const TWO_PI: f64 = 2.0 * PI;
    let n = b.len();
    if n < 2 { return b.to_vec(); }
    let b: Vec<f64> = b.iter().map(|x| x.rem_euclid(TWO_PI)).collect();
    let mut o = Vec::with_capacity(n);
    o.push(b[0]);
    for i in 1..n { o.push(o[i - 1] + short_delta_on_circle_01(b[i - 1], b[i])); }
    let d_close = short_delta_on_circle_01(b[n - 1], b[0]);
    let gap = o[0] - (o[n - 1] + d_close);
    if gap.abs() > 1e-2 {
        let w = (gap / TWO_PI).round() as i32;
        if w != 0 { for i in 1..n { o[i] -= (w as f64) * TWO_PI; } }
    }
    o
}

pub fn align_inner_u_chain_to_outer(outer_umin: f64, outer_umax: f64, o_inner: &mut [f64]) {
    if o_inner.is_empty() { return; }
    const TWO_PI: f64 = 2.0 * PI;
    let mid = 0.5 * (outer_umin + outer_umax);
    let mut best_k = 0i32; let mut best_d = f64::INFINITY;
    for k in -2..=2 {
        let t = o_inner[0] + (k as f64) * TWO_PI;
        let d = (t - mid).abs();
        if d < best_d { best_d = d; best_k = k; }
    }
    let s = (best_k as f64) * TWO_PI;
    for u in o_inner.iter_mut() { *u += s; }
}

// ── Sphere holed mask context ────────────────────────────────────────────

pub struct SphereHoledMaskCtx {
    pub outer_uv: Vec<DVec2>,
    pub outer_3d: Vec<DVec3>,
    pub inner_polys: Vec<Vec<DVec2>>,
    pub inner_3d: Vec<Vec<DVec3>>,
    pub umin: f64, pub umax: f64, pub vmin: f64, pub vmax: f64,
    pub use_uv_winding: bool,
}

pub fn spherical_holed_uv_mask_setup(
    s: &SphericalSurface, brep: &BRep, face: &Face,
) -> Option<SphereHoledMaskCtx> {
    let tol = 1e-5;
    let n_edges = face.outer_wire.edges.len();
    let per_edge = if n_edges > 600 { 4 } else if n_edges > 300 { 8 }
        else if n_edges > 150 { 16 } else if n_edges > 75 { 32 }
        else if n_edges > 30 { 48 } else { 64 };
    let mut outer3 = sample_wire_polyline_3d_with_n(brep, &face.outer_wire, per_edge);
    trim_almost_closed_polyline(&mut outer3, tol);
    if outer3.len() < 3 { return None; }
    let mut b_outer: Vec<f64> = outer3.iter().map(|p| {
        sphere_point_to_uv(s, *p).x.rem_euclid(2.0 * PI)
    }).collect();
    // Fix degenerate u=0 at sphere poles
    {
        let two_pi = 2.0 * PI;
        let n = outer3.len();
        let mut idx = 0;
        while idx < n {
            let p = outer3[idx];
            let v = sphere_point_to_uv(s, p).y;
            let at_pole = v < 1e-12 || (v - PI).abs() < 1e-12;
            if at_pole && b_outer[idx].abs() < 1e-12 {
                let mut run_start = idx;
                while run_start > 0 {
                    let p0 = outer3[run_start - 1];
                    let v0 = sphere_point_to_uv(s, p0).y;
                    if !(v0 < 1e-12 || (v0 - PI).abs() < 1e-12) || b_outer[run_start - 1].abs() >= 1e-12 { break; }
                    run_start -= 1;
                }
                let mut run_end = idx;
                while run_end + 1 < n {
                    let p1 = outer3[run_end + 1];
                    let v1 = sphere_point_to_uv(s, p1).y;
                    if !(v1 < 1e-12 || (v1 - PI).abs() < 1e-12) || b_outer[run_end + 1].abs() >= 1e-12 { break; }
                    run_end += 1;
                }
                let u_prev = if run_start > 0 { b_outer[run_start - 1] } else { b_outer[n - 1] };
                let u_next = if run_end + 1 < n { b_outer[run_end + 1] } else { b_outer[0] };
                let delta = short_delta_on_circle_01(u_prev, u_next);
                let run_len = (run_end - run_start + 1) as f64;
                for k in 0..=(run_end - run_start) {
                    let frac = (k as f64 + 1.0) / (run_len + 1.0);
                    b_outer[run_start + k] = (u_prev + delta * frac).rem_euclid(two_pi);
                }
                idx = run_end + 1;
            } else { idx += 1; }
        }
    }
    let o_outer = unwrap_u_circle_chain_closed(&b_outer);
    let (ou0, ou1) = o_outer.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), u| (a.min(*u), b.max(*u)));
    if (ou1 - ou0).abs() < 1e-8 { return None; }
    let outer_uv: Vec<DVec2> = outer3.iter().zip(o_outer.iter())
        .map(|(p, u)| DVec2::new(*u, sphere_point_to_uv(s, *p).y)).collect();

    let mut inner_polys = Vec::new();
    let mut inner_3d = Vec::new();
    for w in &face.inner_wires {
        let mut h3 = sample_wire_polyline_3d(brep, w);
        trim_almost_closed_polyline(&mut h3, tol);
        if h3.len() < 3 { continue; }
        let b_in: Vec<f64> = h3.iter().map(|p| sphere_point_to_uv(s, *p).x.rem_euclid(2.0 * PI)).collect();
        let mut o_in = unwrap_u_circle_chain_closed(&b_in);
        align_inner_u_chain_to_outer(ou0, ou1, &mut o_in);
        let huv: Vec<DVec2> = h3.iter().zip(o_in.iter())
            .map(|(p, u)| DVec2::new(*u, sphere_point_to_uv(s, *p).y)).collect();
        let h3d = huv.iter().map(|uv| s.point_at(uv.x, uv.y)).collect();
        inner_polys.push(huv); inner_3d.push(h3d);
    }

    // Pole boundary fix
    const POLE_V_THRESHOLD: f64 = 0.15;
    const U_JUMP_THRESHOLD: f64 = PI / 2.0;
    const POLE_N_INSERT: usize = 32;
    let mut fixed_uv = Vec::with_capacity(outer_uv.len() + 8);
    let mut fixed_3d = Vec::with_capacity(outer3.len() + 8);
    let nu = outer_uv.len();
    for i in 0..nu {
        let j = (i + 1) % nu;
        let a_uv = outer_uv[i]; let b_uv = outer_uv[j];
        fixed_uv.push(a_uv); fixed_3d.push(outer3[i]);
        let near_a = a_uv.y < POLE_V_THRESHOLD || a_uv.y > PI - POLE_V_THRESHOLD;
        let near_b = b_uv.y < POLE_V_THRESHOLD || b_uv.y > PI - POLE_V_THRESHOLD;
        let same_pole = (a_uv.y < POLE_V_THRESHOLD && b_uv.y < POLE_V_THRESHOLD)
            || (a_uv.y > PI - POLE_V_THRESHOLD && b_uv.y > PI - POLE_V_THRESHOLD);
        let big_jump = (b_uv.x - a_uv.x).abs() > U_JUMP_THRESHOLD;
        if near_a && near_b && same_pole && big_jump {
            let pole_v = if a_uv.y < POLE_V_THRESHOLD { 0.0 } else { PI };
            for k in 1..POLE_N_INSERT {
                let frac = k as f64 / POLE_N_INSERT as f64;
                let u = a_uv.x + (b_uv.x - a_uv.x) * frac;
                fixed_uv.push(DVec2::new(u, pole_v));
                fixed_3d.push(s.point_at(u, pole_v));
            }
        }
    }
    let mut outer_uv = fixed_uv; let mut outer_3d = fixed_3d;

    // Dedup consecutive near-duplicate UV
    {
        let mut dedup_uv = Vec::with_capacity(outer_uv.len());
        let mut dedup_3d = Vec::with_capacity(outer_3d.len());
        dedup_uv.push(outer_uv[0]); dedup_3d.push(outer_3d[0]);
        for i in 1..outer_uv.len() {
            if (outer_uv[i] - *dedup_uv.last().unwrap()).length_squared() > 1e-16 {
                dedup_uv.push(outer_uv[i]); dedup_3d.push(outer_3d[i]);
            }
        }
        if dedup_uv.len() >= 2 && (dedup_uv[0] - *dedup_uv.last().unwrap()).length_squared() < 1e-16 {
            dedup_uv.pop(); dedup_3d.pop();
        }
        outer_uv = dedup_uv; outer_3d = dedup_3d;
    }

    let mut umin = f64::INFINITY; let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY; let mut vmax = f64::NEG_INFINITY;
    for p in &outer_uv { umin = umin.min(p.x); umax = umax.max(p.x); vmin = vmin.min(p.y); vmax = vmax.max(p.y); }
    for h in &inner_polys { for p in h { umin = umin.min(p.x); umax = umax.max(p.x); vmin = vmin.min(p.y); vmax = vmax.max(p.y); } }
    let raw_u_range = umax - umin;
    let margin_u = raw_u_range * 0.02 + 1e-4;
    let margin_v = (vmax - vmin) * 0.02 + 1e-4;
    umin -= margin_u; umax += margin_u; vmin -= margin_v; vmax += margin_v;
    let vmin = vmin.clamp(0.0, PI); let vmax = vmax.clamp(0.0, PI);

    if !umin.is_finite() || !umax.is_finite() || !vmin.is_finite() || !vmax.is_finite()
        || umax <= umin + 1e-14 || vmax <= vmin + 1e-14 { return None; }

    let two_pi = 2.0 * PI;
    let outer_u_min = outer_uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let outer_u_max = outer_uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let outer_u_span = outer_u_max - outer_u_min;
    let use_uv_winding = raw_u_range <= two_pi + 0.5 && outer_u_span.abs() <= PI + 0.1;

    Some(SphereHoledMaskCtx { outer_uv, outer_3d, inner_polys, inner_3d, umin, umax, vmin, vmax, use_uv_winding })
}

// ── Convex hull & polygon area ──────────────────────────────────────────

pub fn convex_hull_2d_monotone(mut pts: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    if pts.len() <= 1 { return pts; }
    pts.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut lower = Vec::new();
    for &p in &pts {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 1e-18 { lower.pop(); }
        lower.push(p);
    }
    let mut upper = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 1e-18 { upper.pop(); }
        upper.push(p);
    }
    lower.pop(); upper.pop();
    lower.extend(upper); lower
}

pub fn polygon_area_2d_xy(pts: &[(f64, f64)]) -> f64 {
    if pts.len() < 3 { return 0.0; }
    let mut a = 0.0;
    for k in 0..pts.len() {
        let p = pts[k]; let q = pts[(k + 1) % pts.len()];
        a += p.0 * q.1 - p.1 * q.0;
    }
    0.5 * a.abs()
}

pub fn polygon_area_2d_projected(pts: &[DVec3], pivot: DVec3, ux: DVec3, uy: DVec3) -> f64 {
    if pts.len() < 3 { return 0.0; }
    let mut s = 0.0;
    for i in 0..pts.len() { let j = (i + 1) % pts.len();
        let a = (pts[i] - pivot).dot(ux); let b = (pts[i] - pivot).dot(uy);
        let c = (pts[j] - pivot).dot(ux); let d = (pts[j] - pivot).dot(uy);
        s += a * d - c * b;
    }
    0.5 * s
}

pub fn try_boundary_convex_hull_area(
    brep: &BRep, wire: &Wire, pivot: DVec3, ux: DVec3, uy: DVec3,
) -> Option<f64> {
    let mut vert_indices = Vec::new();
    for we in &wire.edges { let edge = brep.flat_edges().get(we.idx).copied()?;
        vert_indices.push(edge.0); vert_indices.push(edge.1); }
    vert_indices.sort(); vert_indices.dedup();
    if vert_indices.len() < 3 { return None; }
    let pts_2d: Vec<DVec2> = vert_indices.iter()
        .filter_map(|&vi| brep.vertex_point(vi))
        .map(|v| DVec2::new((v - pivot).dot(ux), (v - pivot).dot(uy))).collect();
    if pts_2d.len() < 3 { return None; }

    let mut sorted: Vec<usize> = (0..pts_2d.len()).collect();
    sorted.sort_by(|&a, &b| {
        let pa = pts_2d[a]; let pb = pts_2d[b];
        if (pa.x - pb.x).abs() > 1e-12 { pa.x.partial_cmp(&pb.x).unwrap_or(std::cmp::Ordering::Equal) }
        else { pa.y.partial_cmp(&pb.y).unwrap_or(std::cmp::Ordering::Equal) }
    });
    let cross = |o: DVec2, a: DVec2, b: DVec2| (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    let mut lower = Vec::new();
    for &si in &sorted {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], pts_2d[si]) <= 0.0 { lower.pop(); }
        lower.push(pts_2d[si]);
    }
    let mut upper = Vec::new();
    for &si in sorted.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], pts_2d[si]) <= 0.0 { upper.pop(); }
        upper.push(pts_2d[si]);
    }
    lower.pop(); upper.pop();
    let hull: Vec<DVec2> = lower.into_iter().chain(upper).collect();
    if hull.len() < 3 { return None; }
    let mut area = 0.0;
    for i in 0..hull.len() {
        let j = (i + 1) % hull.len();
        area += hull[i].x * hull[j].y; area -= hull[j].x * hull[i].y;
    }
    Some((area * 0.5).abs())
}

// ── Face triangulation & orientation ─────────────────────────────────────

pub fn orient_by_ref(tri: [DVec3; 3], n_ref: DVec3) -> [DVec3; 3] {
    let n = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
    if n.dot(n_ref) < 0.0 { [tri[0], tri[2], tri[1]] } else { tri }
}

pub fn orient_tri(tri: [DVec3; 3], face_normal: DVec3) -> [DVec3; 3] {
    let n = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
    if n.dot(face_normal) < 0.0 { [tri[0], tri[2], tri[1]] } else { tri }
}

pub fn tessellate_curved_face(
    brep: &BRep, face: &Face, face_flat_idx: usize,
) -> Option<Vec<[DVec3; 3]>> {
    let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
    })?;

    // For sphere faces with inner wires, try spherical earcut + planar fallback
    if let Surface3::Sphere(s) = &surf_idx {
        if !face.inner_wires.is_empty() {
            // Spherical earcut is the most reliable for holed sphere faces
            let outer_pts = sample_wire_polyline_3d(brep, &face.outer_wire);
            let mut outer_cln = outer_pts.clone();
            trim_almost_closed_polyline(&mut outer_cln, 1e-5);
            if outer_cln.len() >= 3 {
                if let Some(tris) = try_spherical_earcut_simple(s, &outer_cln, face.normal) {
                    return Some(tris);
                }
            }
        }
        // Try great-circle analytic earcut
        if let Some(tris) = try_spherical_earcut_simple(s, &sample_wire_polyline_3d(brep, &face.outer_wire), face.normal) {
            return Some(tris);
        }
    }

    // Generic: UV grid tessellation for curved surfaces
    let [u0, u1, v0, v1] = estimate_uv_domain_from_wire(brep, face, &surf_idx)?;
    const N: usize = 40;
    let du = (u1 - u0) / N as f64;
    let dv = (v1 - v0) / N as f64;
    let mut tris = Vec::new();
    for i in 0..N {
        for j in 0..N {
            let u_lo = u0 + i as f64 * du; let u_hi = u_lo + du;
            let v_lo = v0 + j as f64 * dv; let v_hi = v_lo + dv;
            let p00 = surf_idx.point_at(u_lo, v_lo);
            let p10 = surf_idx.point_at(u_hi, v_lo);
            let p01 = surf_idx.point_at(u_lo, v_hi);
            let p11 = surf_idx.point_at(u_hi, v_hi);
            tris.push(orient_tri([p00, p10, p11], face.normal));
            tris.push(orient_tri([p00, p11, p01], face.normal));
        }
    }
    Some(tris)
}

pub fn estimate_uv_domain_from_wire(
    brep: &BRep, face: &Face, surf: &Surface3,
) -> Option<[f64; 4]> {
    // Prefer stored uv_domain
    if let Some(face_ts) = brep.tshapes.iter().find(|ts| matches!(ts.as_ref(), topods::TShape::Face(_))) {
        if let topods::TShape::Face(fd) = &**face_ts {
            if let Some(domain) = fd.uv_domain {
                return Some(domain);
            }
        }
    }
    let pts = sample_wire_polyline_3d(brep, &face.outer_wire);
    if pts.is_empty() { return Some(surf.default_domain()); }
    use crate::geom::SurfaceEval;
    let proj: Vec<DVec2> = pts.iter().map(|&p| {
        let r = crate::math::geom_api::project::closest_point_on_surface(surf, p, 64);
        DVec2::new(r.params.0, r.params.1)
    }).collect();
    let mut u0 = f64::INFINITY; let mut u1 = f64::NEG_INFINITY;
    let mut v0 = f64::INFINITY; let mut v1 = f64::NEG_INFINITY;
    for p in &proj { u0 = u0.min(p.x); u1 = u1.max(p.x); v0 = v0.min(p.y); v1 = v1.max(p.y); }
    let margin_u = (u1 - u0) * 0.02 + 1e-4;
    let margin_v = (v1 - v0) * 0.02 + 1e-4;
    Some([u0 - margin_u, u1 + margin_u, v0 - margin_v, v1 + margin_v])
}

pub fn face_flat_iter(brep: &topods::BRep) -> Vec<(usize, Face)> {
    let mut faces = Vec::new();
    for (ti, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Face(fd) = &**ts {
            let shape_to_edge = |sh: &crate::topo::topo_shape::Shape| -> WireEdge {
                WireEdge { idx: sh.index, forward: sh.orientation.is_forward() }
            };
            let outer_wire = {
                let wi = fd.outer_wire.index;
                if let topods::TShape::Wire(wd) = &*brep.tshapes[wi] {
                    Wire { edges: wd.edges.iter().map(shape_to_edge).collect() }
                } else { continue; }
            };
            let inner_wires: Vec<Wire> = fd.inner_wires.iter().filter_map(|sh| {
                let wi = sh.index;
                if let topods::TShape::Wire(wd) = &*brep.tshapes[wi] {
                    Some(Wire { edges: wd.edges.iter().map(shape_to_edge).collect() })
                } else { None }
            }).collect();
            let normal = fd.surface.as_ref()
                .map(|s| crate::geom::SurfaceEval::normal_at(s, 0.0, 0.0))
                .unwrap_or_default();
            let face = Face {
                outer_wire,
                inner_wires,
                normal,
                triangles: Vec::new(),
                sample_point: fd.sample_point,
                mesh_dirty: true,
                surface_idx: None,
            };
            faces.push((ti, face));
        }
    }
    faces
}

pub fn face_triangles(
    brep: &BRep, face: &Face, face_flat_idx: usize,
) -> Vec<[DVec3; 3]> {
    // Try holes first (earcut-based)
    if let Some(tris) = try_face_with_holes(brep, face, face_flat_idx) { return tris; }
    // Try planar earcut for simple outer
    let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
    trim_almost_closed_polyline(&mut outer, 1e-5);
    if outer.len() >= 3 {
        if let Some(tris) = try_planar_earcut_simple_outer(&outer, face.normal) { return tris; }
    }
    // Try UV grid tessellation
    if let Some(tris) = tessellate_curved_face(brep, face, face_flat_idx) { return tris; }
    // Fan triangulation fallback
    if outer.len() >= 3 {
        let mut tris = Vec::new();
        let ctr = outer.iter().sum::<DVec3>() / outer.len() as f64;
        for i in 0..outer.len() {
            let j = (i + 1) % outer.len();
            tris.push(orient_tri([outer[i], outer[j], ctr], face.normal));
        }
        return tris;
    }
    Vec::new()
}

pub fn face_triangles_pub(brep: &BRep, face_flat_idx: usize) -> Vec<[DVec3; 3]> {
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        if *fi == face_flat_idx { return face_triangles(brep, face, *fi); }
    }
    Vec::new()
}
