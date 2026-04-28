//! Shape properties: surface area, volume, and centroid.
//!
//! Analogous to OCCT `GProp_GProps` with `BRepGProp`.
//!
//! `surface_area` first uses analytic or parametric integrals (plane: shoe-lace; sphere: UV
//! mask and `R² dΩ`; other finite UV patches without `inner_wires`: cylinder `r·Δu·Δv` or a
//! midpoint `‖∂P/∂u×∂P/∂v‖` sum on the same domain as `tessellate_curved_face`), then
//! triangulated faces. Where triangles are used: UV-grid tessellation, holed
//! ear-cut, or fan-triangulation from wire vertices.

use glam::{DVec2, DVec3};

use crate::BRep;
use crate::geom::{SphericalSurface, Surface3, SurfaceEval};
use crate::topology::{Face, Wire};

// ── Internal helpers ──────────────────────────────────────────────────────────

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
fn sample_wire_polyline_3d(brep: &BRep, wire: &Wire) -> Vec<DVec3> {
    const EDGE_SAMPLE_N: usize = 64;
    use crate::geom::CurveEval;
    let mut pts = Vec::new();
    for we in &wire.edges {
        let edge = match brep.edges.get(we.idx) {
            Some(e) => e,
            None => continue,
        };
        let curve_opt = brep.geom.edge_curve.get(we.idx).and_then(|o| *o)
            .and_then(|ci| brep.geom.curves.get(ci));
        if let Some(curve) = curve_opt {
            let range = brep
                .geom
                .edge_curve_range
                .get(we.idx)
                .and_then(|o| *o)
                .unwrap_or_else(|| curve.default_domain());
            let [t0, t1] = if we.forward { range } else { [range[1], range[0]] };
            if (t1 - t0).abs() > 1e-12 && t0.is_finite() && t1.is_finite() {
                let n = EDGE_SAMPLE_N;
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
        let vidx = if we.forward { edge.start } else { edge.end };
        if let Some(v) = brep.vertices.get(vidx) {
            pts.push(v.point);
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

/// Spherical (u, v) with u ∈ [−π, π] (longitude) and v ∈ [0, π] (colatitude),
/// matching [`SphericalSurface`](crate::geom::SphericalSurface) / `SurfaceEval`.
fn sphere_point_to_uv(s: &SphericalSurface, p: DVec3) -> DVec2 {
    (*s).world_to_uv(p)
}

/// Remove 2π jumps in u so consecutive samples stay in one branch (for earcut).
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
    unwrap_sphere_u_in_chain(&mut outer_uv);
    let mut flat: Vec<f64> = Vec::with_capacity(2 * outer.len());
    for uv in &outer_uv {
        flat.push(uv.x);
        flat.push(uv.y);
    }
    let all_3d: Vec<DVec3> = outer.to_vec();
    let indices = earcut_indices_from_flat(&flat, &[]);
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

/// One outer ring only in the plane (pivot) ⟂ `face_normal` (for nearly planar caps).
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
    let surf_idx = brep.geom.face_surface.get(face_flat_idx)?.as_ref().copied()?;
    let surf = brep.geom.surfaces.get(surf_idx)?;
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
/// 64×64 gives a <0.1% volume error for a unit sphere.
const UV_TESS_N: usize = 64;

/// Finer grid for masked sphere patches (trimmed / holed) when ear-cut fails.
const SPHERE_UV_MASK_N: usize = 160;

/// Grid for ‖∂P/∂u×∂P/∂v‖ midpoint integration on a UV rectangle (per-axis).
const PARAM_RECT_AREA_N: usize = 64;

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

/// Shortest signed step from `a` to `b` on the circle, with `a, b` reduced mod 2π to [0,2π).
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

/// Cumulative u along a **closed** wire, using only short steps on S¹, then remove full 2π
/// turns so the first and last point close in ℝ (winding 0 in the (u, v) chart).
/// This fixes `unwrap_sphere_u_in_chain` sometimes representing a small patch with span > 2π.
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

/// Add +2πk to an inner u-chain so it lies near the outer patch in ℝ (shared chart for PIP).
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
    inner_polys: Vec<Vec<DVec2>>,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
}

fn spherical_holed_uv_mask_setup(
    s: &SphericalSurface,
    brep: &BRep,
    face: &Face,
) -> Option<SphereHoledMaskCtx> {
    let tol = 1e-5_f64;
    let mut outer3 = sample_wire_polyline_3d(brep, &face.outer_wire);
    trim_almost_closed_polyline(&mut outer3, tol);
    if outer3.len() < 3 {
        return None;
    }
    let b_outer: Vec<f64> = outer3
        .iter()
        .map(|p| sphere_point_to_uv(s, *p).x.rem_euclid(2.0 * std::f64::consts::PI))
        .collect();
    let o_outer = unwrap_u_circle_chain_closed(&b_outer);
    let (ou0, ou1) = o_outer.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(a, b), u| (a.min(*u), b.max(*u)),
    );
    let outer_uv: Vec<DVec2> = outer3
        .iter()
        .zip(o_outer.iter())
        .map(|(p, u)| DVec2::new(*u, sphere_point_to_uv(s, *p).y))
        .collect();

    let mut inner_polys: Vec<Vec<DVec2>> = Vec::new();
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
        inner_polys.push(huv);
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
    let margin_u = (umax - umin) * 0.02 + 1e-4;
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
    Some(SphereHoledMaskCtx {
        outer_uv,
        inner_polys,
        umin,
        umax,
        vmin,
        vmax,
    })
}

/// Sum `∫ R² sin v dudv` over the same masked grid (no triangulation), for GProp-style area.
fn sphere_holed_mask_param_area_sum(s: &SphericalSurface, ctx: &SphereHoledMaskCtx) -> f64 {
    let nu = SPHERE_UV_MASK_N;
    let nv = SPHERE_UV_MASK_N;
    let umin = ctx.umin;
    let umax = ctx.umax;
    let vmin = ctx.vmin;
    let vmax = ctx.vmax;
    let du = (umax - umin) / nu as f64;
    let dv = (vmax - vmin) / nv as f64;
    let r2 = s.radius * s.radius;
    let inner = &ctx.inner_polys;
    let emit = |outer_poly: &[DVec2], use_inner: bool| -> f64 {
        let mut a = 0.0f64;
        for i in 0..nu {
            for j in 0..nv {
                let u0 = umin + i as f64 * du;
                let u1 = u0 + du;
                let v0 = vmin + j as f64 * dv;
                let v1 = v0 + dv;
                let corners = [
                    DVec2::new(u0, v0),
                    DVec2::new(u1, v0),
                    DVec2::new(u1, v1),
                    DVec2::new(u0, v1),
                ];
                if !corners
                    .iter()
                    .all(|q| point_in_polygon_2d(outer_poly, *q))
                {
                    continue;
                }
                if use_inner
                    && inner.iter().any(|h| {
                        corners
                            .iter()
                            .any(|q| point_in_polygon_2d(h, *q))
                    })
                {
                    continue;
                }
                a += r2 * du * (v0.cos() - v1.cos());
            }
        }
        a
    };
    let mut t = emit(&ctx.outer_uv, true);
    if t <= 0.0 && !inner.is_empty() {
        t = emit(&ctx.outer_uv, false);
    }
    if t > 0.0 {
        return t;
    }
    let mut rev = ctx.outer_uv.clone();
    rev.reverse();
    t = emit(&rev, true);
    if t <= 0.0 && !inner.is_empty() {
        t = emit(&rev, false);
    }
    t
}

/// When the face normal is world-axis-aligned and the outer loop is exactly the boundary of its
/// own axis-aligned bounding box in world (e.g. each side of a merged `box` union), dense edge
/// sampling can make simple shoe-lace on 3D samples lose area. If every sample lies on one of the
/// four sides of that box, use `width * height` in world UV. Fails for L-shapes (re-entrant
/// points inside the bbox), circles, etc. — then we fall back to shoe-lace below.
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
    let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
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

/// Shoelace area of outer wire minus |hole areas| in the face plane (pivot = first outer point).
fn try_planar_face_area_shoelace(
    brep: &BRep,
    face: &Face,
    face_normal: DVec3,
) -> Option<f64> {
    if face.inner_wires.is_empty() {
        if let Some(a) = try_axis_aligned_world_rect_plane_area(brep, face, face_normal) {
            return Some(a);
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
    for w in &face.inner_wires {
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

/// Prefer analytic / parametric area for `surface_area`: plane (shoelace); all sphere faces
/// (UV polygon mask + `R² dΩ`); finite-UV rectangular patches on other surfaces without inner
/// wires (cylinder exact; otherwise `‖Pu×Pv‖` midpoint rule on the same domain as tessellation).
fn try_analytic_face_surface_area(
    brep: &BRep,
    face: &Face,
    face_flat_idx: usize,
) -> Option<f64> {
    let surf_idx = brep.geom.face_surface.get(face_flat_idx).copied().flatten()?;
    let surf = brep.geom.surfaces.get(surf_idx)?;
    match surf {
        Surface3::Plane(_) => try_planar_face_area_shoelace(brep, face, face.normal),
        Surface3::Sphere(s) => {
            let ctx = spherical_holed_uv_mask_setup(s, brep, face)?;
            let v = sphere_holed_mask_param_area_sum(s, &ctx);
            if v > 0.0 { Some(v) } else { None }
        }
        _ if !face.inner_wires.is_empty() => None,
        _ => {
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
                Surface3::Cylinder(c) => Some(
                    c.radius * (u1 - u0).abs() * (v1 - v0).abs(),
                ),
                _ => param_rect_area_cross(surf, u0, u1, v0, v1),
            }
        }
    }
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
    let inner_polys = &ctx.inner_polys;

    let nu = SPHERE_UV_MASK_N;
    let nv = SPHERE_UV_MASK_N;
    let du = (umax - umin) / nu as f64;
    let dv = (vmax - vmin) / nv as f64;

    let emit_grid = |outer_poly: &[DVec2], use_inner_mask: bool| -> Vec<[DVec3; 3]> {
        let mut tris: Vec<[DVec3; 3]> = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let u0 = umin + i as f64 * du;
                let u1 = u0 + du;
                let v0 = vmin + j as f64 * dv;
                let v1 = v0 + dv;
                let corners = [
                    DVec2::new(u0, v0),
                    DVec2::new(u1, v0),
                    DVec2::new(u1, v1),
                    DVec2::new(u0, v1),
                ];
                // All corners in the outer UV loop; any corner in a hole excludes the cell
                // (tighter than centre-only, reduces PIP over-count at the trim).
                if !corners
                    .iter()
                    .all(|q| point_in_polygon_2d(outer_poly, *q))
                {
                    continue;
                }
                if use_inner_mask
                    && inner_polys.iter().any(|h| {
                        corners
                            .iter()
                            .any(|q| point_in_polygon_2d(h, *q))
                    })
                {
                    continue;
                }
                let uc = umin + (i as f64 + 0.5) * du;
                let vc = vmin + (j as f64 + 0.5) * dv;
                let p00 = s.point_at(u0, v0);
                let p10 = s.point_at(u1, v0);
                let p11 = s.point_at(u1, v1);
                let p01 = s.point_at(u0, v1);
                let pc = s.point_at(uc, vc);
                let nref = s.normal_at(uc, vc);
                // Exact area in (u, v) for a sphere: R² · du · (cos v0 − cos v1); isotropic
                // scale of the chordal bilinear patch toward the cell centre to match dA.
                let r2 = s.radius * s.radius;
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

    let mut tris = emit_grid(outer_uv, true);
    if tris.is_empty() && !inner_polys.is_empty() {
        tris = emit_grid(outer_uv, false);
    }
    if tris.is_empty() {
        let mut ou = outer_uv.clone();
        ou.reverse();
        tris = emit_grid(&ou, true);
        if tris.is_empty() && !inner_polys.is_empty() {
            tris = emit_grid(&ou, false);
        }
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
/// UV domain on a regular `UV_TESS_N × UV_TESS_N` grid.
///
/// Returns triangles oriented outward (consistent with the surface normal).
/// UV rectangle used by [`tessellate_curved_face`] and parametric `surface_area` fallbacks
/// (same resolution priority as triangulation: range → finite `default_domain` → wire estimate).
fn curved_face_uv_domain(
    brep: &BRep,
    face: &Face,
    face_flat_idx: usize,
    surf: &Surface3,
) -> Option<[f64; 4]> {
    if let Some(Some(r)) = brep.geom.face_surface_range.get(face_flat_idx) {
        Some(*r)
    } else {
        let d = surf.default_domain();
        if d.iter().all(|x| x.is_finite()) {
            Some(d)
        } else {
            estimate_uv_domain_from_wire(brep, face, surf)
        }
    }
}

/// `∫∫ ‖∂P/∂u×∂P/∂v‖ dudv` on `[u0,u1]×[v0,v1]` (midpoint rule, central differences for partials).
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
    let surf_idx = brep.geom.face_surface.get(face_flat_idx)?.as_ref().copied()?;
    let surf = brep.geom.surfaces.get(surf_idx)?;

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

    // Build a (nu+1)×(nv+1) grid of 3-D points.
    let mut pts = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        let u = u0 + (u1 - u0) * (i as f64 / nu as f64);
        for j in 0..=nv {
            let v = v0 + (v1 - v0) * (j as f64 / nv as f64);
            pts.push(surf.point_at(u, v));
        }
    }

    // Emit two triangles per quad cell (i,j)–(i+1,j)–(i,j+1)–(i+1,j+1).
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
    surf: &crate::geom::Surface3,
) -> Option<[f64; 4]> {
    use crate::geom::Surface3;

    // Collect all wire vertex 3-D points (outer + inner wires).
    let all_wires = std::iter::once(&face.outer_wire).chain(face.inner_wires.iter());
    let pts: Vec<DVec3> = all_wires
        .flat_map(|w| &w.edges)
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let vidx = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vidx).map(|v| v.point)
        })
        .collect();

    if pts.is_empty() {
        return None;
    }

    match surf {
        Surface3::Cylinder(cyl) => {
            // CylindricalSurface: u = azimuth [0, 2π], v = height along axis.
            let d = surf.default_domain(); // [0, 2π, -inf, inf]
            let u0 = d[0];
            let u1 = d[1];
            let v_vals: Vec<f64> = pts.iter().map(|p| (*p - cyl.origin).dot(cyl.axis)).collect();
            let v0 = v_vals.iter().cloned().fold(f64::INFINITY, f64::min);
            let v1 = v_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if v0 >= v1 { return None; }
            Some([u0, u1, v0 - 1e-10, v1 + 1e-10])
        }
        Surface3::Cone(con) => {
            // ConicalSurface: u = azimuth [0, 2π], v = slant distance ≥ 0.
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
    // Holed / trimmed faces: the boolean result may cache outer-loop triangles only
    // (see `ResultBuilder`); re-mesh in plane/UV for area / volume.
    if !face.inner_wires.is_empty() {
        if let Some(tris) = try_face_with_holes(brep, face, face_flat_idx) {
            return tris;
        }
    } else {
        // Spherical sub-faces: prefer UV ear-cut of the boundary loop over builder
        // 3D ear-clip, which can under-fill a non-convex (trimmed) patch.
        if let Some(sidx) = brep
            .geom
            .face_surface
            .get(face_flat_idx)
            .and_then(|o| *o)
        {
            if let Some(Surface3::Sphere(s)) = brep.geom.surfaces.get(sidx) {
                let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
                trim_almost_closed_polyline(&mut outer, 1e-5);
                if outer.len() >= 3 {
                    if let Some(tris) = try_spherical_earcut_simple(s, &outer, face.normal)
                        .or_else(|| try_planar_earcut_simple_outer(&outer, face.normal))
                        .or_else(|| {
                            try_spherical_uv_masked_raster(
                                s,
                                brep,
                                face,
                                face_flat_idx,
                                face.normal,
                            )
                        })
                    {
                        return tris;
                    }
                }
            }
        }
    }

    if face.inner_wires.is_empty() && !face.triangles.is_empty() {
        return face
            .triangles
            .iter()
            .filter_map(|&[i, j, k]| {
                let a = brep.vertices.get(i)?.point;
                let b = brep.vertices.get(j)?.point;
                let c = brep.vertices.get(k)?.point;
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
fn face_flat_iter(brep: &BRep) -> impl Iterator<Item = (usize, &crate::topology::Face)> {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .scan(0usize, |idx, face| {
            let i = *idx;
            *idx += 1;
            Some((i, face))
        })
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Public re-export of `face_triangles` for diagnostic use in tests.
/// Not part of the stable API.
#[doc(hidden)]
pub fn face_triangles_pub(
    brep: &BRep,
    face: &crate::topology::Face,
    face_flat_idx: usize,
) -> Vec<[DVec3; 3]> {
    face_triangles(brep, face, face_flat_idx)
}

/// Compute the total surface area of all faces in the BRep.
///
/// For each face, uses analytic area when available (planar shoe lace; holed
/// sphere: parametric `R² dΩ` on the same UV mask as the raster), otherwise
/// sums triangle areas (pre-triangulated, UV-sampled, or fan-triangulated).
/// Returns 0.0 if the BRep has no faces.
pub fn surface_area(brep: &BRep) -> f64 {
    let mut total = 0.0f64;
    for (fi, f) in face_flat_iter(brep) {
        total += face_surface_area(brep, f, fi);
    }
    total
}

/// Area of one face, using the same rules as [`surface_area`].
#[doc(hidden)]
pub fn face_surface_area(brep: &BRep, face: &Face, face_flat_idx: usize) -> f64 {
    if let Some(a) = try_analytic_face_surface_area(brep, face, face_flat_idx) {
        a
    } else {
        face_triangles(brep, face, face_flat_idx)
            .iter()
            .map(|&[a, b, c]| tri_area(a, b, c))
            .sum()
    }
}

/// Compute the signed volume of the closed BRep solid.
///
/// Uses the divergence theorem: V = (1/6) Σ_triangles a·(b×c).
/// Works correctly for a closed, consistently-oriented mesh.
/// Returns 0.0 for open shells or empty BReps.
pub fn volume(brep: &BRep) -> f64 {
    face_flat_iter(brep)
        .flat_map(|(fi, f)| face_triangles(brep, f, fi))
        .map(|[a, b, c]| tet_signed_volume(a, b, c))
        .sum::<f64>()
        .abs()
}

/// Compute the centroid (center of mass) of the solid by volumetric integration.
///
/// Uses the formula: C = (1 / 8V) Σ_triangles (a+b+c) * tet_signed_vol(a,b,c)
/// where the sum is over all surface triangles.
///
/// Falls back to `BRep::center()` (vertex average) if the volume is near zero.
pub fn centroid(brep: &BRep) -> DVec3 {
    let mut vol_sum = 0.0_f64;
    let mut weighted_sum = DVec3::ZERO;

    for (fi, face) in face_flat_iter(brep) {
        for [a, b, c] in face_triangles(brep, face, fi) {
            let sv = tet_signed_volume(a, b, c);
            vol_sum += sv;
            // Weight the centroid of each tet (at (a+b+c+origin)/4,
            // origin=0) → simplified to (a+b+c) * sv
            weighted_sum += (a + b + c) * sv;
        }
    }

    if vol_sum.abs() < 1e-15 {
        return brep.center();
    }

    // Centroid formula: (1/(2 * 4 * vol_sum)) * Σ (a+b+c) * sv
    // Simplification: weighted_sum / (4 * vol_sum) gives tet centroid average
    weighted_sum / (4.0 * vol_sum)
}

// ── Inertia tensor ────────────────────────────────────────────────────────────

/// Symmetric 3×3 moment of inertia tensor (assuming uniform density = 1).
///
/// The components are defined as:
/// ```text
/// Ixx = ∫(y²+z²) dV,  Iyy = ∫(x²+z²) dV,  Izz = ∫(x²+y²) dV
/// Ixy = -∫xy dV,       Ixz = -∫xz dV,       Iyz = -∫yz dV
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
    /// Returns the 3×3 inertia matrix as row-major `[[f64;3];3]`.
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
pub fn inertia_tensor(brep: &BRep) -> InertiaTensor {
    let mut ixx = 0.0_f64;
    let mut iyy = 0.0_f64;
    let mut izz = 0.0_f64;
    let mut ixy = 0.0_f64;
    let mut ixz = 0.0_f64;
    let mut iyz = 0.0_f64;

    for (fi, face) in face_flat_iter(brep) {
        for [a, b, c] in face_triangles(brep, face, fi) {
            // Signed volume of tet (origin, a, b, c)
            // sv = a·(b×c)/6 — same as tet_signed_volume
            let sv = a.dot(b.cross(c)) / 6.0;

            // Symmetric quadratic sums for each coordinate pair.
            // For ∫_tet x² dV = sv/10 * x2_sym (from simplex integration).
            let x2 = a.x * a.x + b.x * b.x + c.x * c.x + a.x * b.x + a.x * c.x + b.x * c.x;
            let y2 = a.y * a.y + b.y * b.y + c.y * c.y + a.y * b.y + a.y * c.y + b.y * c.y;
            let z2 = a.z * a.z + b.z * b.z + c.z * c.z + a.z * b.z + a.z * c.z + b.z * c.z;

            ixx += sv / 10.0 * (y2 + z2);
            iyy += sv / 10.0 * (x2 + z2);
            izz += sv / 10.0 * (x2 + y2);

            // For ∫_tet xy dV = sv/20 * xy_mixed (from simplex integration).
            // Product-moment: Ixy = -∫xy dV, etc.
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
    // Off-diagonal sign: Ixy = -∫xy dV so negate the accumulated sums.
    InertiaTensor {
        ixx: ixx.abs(),
        iyy: iyy.abs(),
        izz: izz.abs(),
        ixy: -ixy,
        ixz: -ixz,
        iyz: -iyz,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PrimitiveSolid;

    const EPS: f64 = 1e-6;

    #[test]
    fn unit_box_surface_area() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let area = surface_area(&brep);
        assert!(
            (area - 6.0).abs() < EPS,
            "unit box surface area should be 6, got {area}"
        );
    }

    #[test]
    fn unit_box_volume() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let vol = volume(&brep);
        assert!(
            (vol - 1.0).abs() < EPS,
            "unit box volume should be 1, got {vol}"
        );
    }

    #[test]
    fn box_2x3x4_volume() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        let vol = volume(&brep);
        assert!(
            (vol - 24.0).abs() < EPS,
            "2×3×4 box volume should be 24, got {vol}"
        );
    }

    #[test]
    fn box_2x3x4_surface_area() {
        // SA = 2*(2*3 + 3*4 + 2*4) = 2*(6+12+8) = 52
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        let area = surface_area(&brep);
        assert!(
            (area - 52.0).abs() < EPS,
            "2×3×4 box SA should be 52, got {area}"
        );
    }

    /// Each face is a rectangle: two 2×3, two 3×4, two 2×4 — exercises per-face analytic plane area.
    #[test]
    fn box_non_cuboid_per_face_rect_areas() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        let mut areas = Vec::new();
        let mut i = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    areas.push(face_surface_area(&brep, face, i));
                    i += 1;
                }
            }
        }
        assert_eq!(areas.len(), 6, "box should have 6 faces");
        areas.sort_by(|a, b| a.total_cmp(b));
        for (a, e) in areas.iter().zip([6.0_f64, 6.0, 8.0, 8.0, 12.0, 12.0].iter()) {
            assert!((a - e).abs() < 1e-3, "face area {a} expected {e}");
        }
    }

    #[test]
    fn box_per_face_area_sum_matches_surface_area() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        let mut sum = 0.0;
        let mut i = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    sum += face_surface_area(&brep, face, i);
                    i += 1;
                }
            }
        }
        let tot = surface_area(&brep);
        assert!((sum - tot).abs() < 1e-4, "sum of face areas {sum} vs surface_area {tot}");
    }

    #[test]
    fn unit_box_centroid() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let c = centroid(&brep);
        // unit box: centroid at (0.5, 0.5, 0.5)
        assert!(
            (c - DVec3::splat(0.5)).length() < 1e-4,
            "centroid should be (0.5,0.5,0.5), got {c}"
        );
    }

    #[test]
    fn unit_box_inertia_tensor_diagonal_equal() {
        // Unit box [0,1]^3 about the world origin:
        // Ixx = ∫(y²+z²)dV = (1/3 + 1/3) = 2/3
        // By symmetry, Iyy = Izz = 2/3
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let it = inertia_tensor(&brep);
        let expected = 2.0 / 3.0;
        let tol = 1e-4;
        assert!(
            (it.ixx - expected).abs() < tol,
            "Ixx = {} expected {}",
            it.ixx,
            expected
        );
        assert!(
            (it.iyy - expected).abs() < tol,
            "Iyy = {} expected {}",
            it.iyy,
            expected
        );
        assert!(
            (it.izz - expected).abs() < tol,
            "Izz = {} expected {}",
            it.izz,
            expected
        );
    }

    #[test]
    fn box_2x1x1_inertia_tensor() {
        // Box [0,2]×[0,1]×[0,1] about origin:
        // Ixx = ∫(y²+z²)dV = V*(1/3+1/3) = 2*(2/3) = 4/3
        // Iyy = ∫(x²+z²)dV = V*(4/3÷2 + 1/3) = 2*(2/3+1/3) = 2*(1) = wait:
        //   ∫₀²∫₀¹∫₀¹ (x²+z²) dx dy dz  but order matters since box is [0,2]x[0,1]x[0,1]
        //   = 1*1*(∫₀² x² dx) + 1*2*(∫₀¹ z² dz) = (8/3) + 2*(1/3) = 8/3+2/3 = 10/3
        // Izz = ∫(x²+y²)dV = (8/3) + 2*(1/3) = 10/3
        // Ixx = 2*(1/3) + 2*(1/3) = 4/3
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 1.0,
            depth: 1.0,
        });
        let it = inertia_tensor(&brep);
        let tol = 1e-3;
        let expected_ixx = 4.0 / 3.0;
        let expected_iyy = 10.0 / 3.0;
        let expected_izz = 10.0 / 3.0;
        assert!(
            (it.ixx - expected_ixx).abs() < tol,
            "Ixx = {} expected {}",
            it.ixx,
            expected_ixx
        );
        assert!(
            (it.iyy - expected_iyy).abs() < tol,
            "Iyy = {} expected {}",
            it.iyy,
            expected_iyy
        );
        assert!(
            (it.izz - expected_izz).abs() < tol,
            "Izz = {} expected {}",
            it.izz,
            expected_izz
        );
    }

    // ── Curved primitive tests (UV tessellation path) ─────────────────────────

    #[test]
    fn unit_sphere_volume() {
        // V = (4/3)π r³ = 4.18879...  for r=1
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let vol = volume(&brep);
        let expected = 4.0 / 3.0 * std::f64::consts::PI;
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 5e-3,
            "unit sphere volume: got {vol:.6}, expected {expected:.6}, rel_err={rel_err:.4}"
        );
    }

    #[test]
    fn sphere_r2_volume() {
        // V = (4/3)π·8 = 33.5103...  for r=2
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });
        let vol = volume(&brep);
        let expected = 4.0 / 3.0 * std::f64::consts::PI * 8.0;
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 5e-3,
            "r=2 sphere volume: got {vol:.6}, expected {expected:.6}, rel_err={rel_err:.4}"
        );
    }

    #[test]
    fn unit_cylinder_volume() {
        // V = π r² h = π for r=1, h=1
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 1.0,
        });
        let vol = volume(&brep);
        let expected = std::f64::consts::PI;
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 5e-3,
            "unit cylinder volume: got {vol:.6}, expected {expected:.6}, rel_err={rel_err:.4}"
        );
    }
}
