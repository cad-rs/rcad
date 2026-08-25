//! OCCT BRepGProp::SurfaceProperties: surface area computation.
//!
//! Split from the original properties.rs. Uses analytic per-type paths
//! (plane shoe-lace, cylinder R×ΔU×ΔV, sphere great-circle, etc.)
//! and falls back to Gauss-Legendre integration or triangulation.

use glam::{DVec2, DVec3};
use std::f64::consts::PI;

use crate::BRep;
use crate::geom::{
    ConicalSurface, Curve3, CurveEval, CylindricalSurface, SphericalSurface, Surface3, SurfaceEval,
};
use crate::topo::topods;
use crate::topo::topology::{Face, Wire, WireEdge};
use crate::base::gprop::tri::*;
use crate::base::geom_api::project::closest_point_on_surface;

/// Total surface area of a BRep.
pub fn surface_area(brep: &topods::BRep) -> f64 {
    let mut total = 0.0;
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        total += face_surface_area(brep, face, *fi);
    }
    total
}

/// Surface area of a single face.
pub fn face_surface_area(brep: &BRep, face: &Face, face_flat_idx: usize) -> f64 {
    try_analytic_face_surface_area(brep, face, face_flat_idx)
        .unwrap_or_else(|| {
            let tris = face_triangles_pub(brep, face_flat_idx);
            tris.iter().map(|[a, b, c]| tri_area(*a, *b, *c)).sum()
        })
}

pub fn try_analytic_face_surface_area_pub(brep: &BRep, face: &Face, face_flat_idx: usize) -> Option<f64> {
    try_analytic_face_surface_area(brep, face, face_flat_idx)
}

fn try_analytic_face_surface_area(brep: &BRep, face: &Face, face_flat_idx: usize) -> Option<f64> {
    // OCCT BRepGProp_Domain::Next (BRepGProp_Domain.cxx L27-38) skips
    // INTERNAL/EXTERNAL edges in the boundary integral.  The analytic paths
    // below integrate whole wires (shoelace over sampled wire polylines) and
    // cannot represent a skipped single edge, so a face carrying internal
    // edges falls back to the per-edge Green integral
    // (base::gprop::surface::face_surface_area = OCCT BRepGProp_Gauss::Compute,
    // which skips internal edges).
    let has_internal = face.outer_wire.edges.iter().any(|we| we.internal)
        || face.inner_wires.iter().any(|w| w.edges.iter().any(|we| we.internal));
    if has_internal {
        return Some(crate::base::gprop::surface::face_surface_area(
            brep, face, face_flat_idx,
        ));
    }
    let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
    })?;
    let surf = &surf_idx;

    match surf {
        Surface3::Plane(_) => {
            let r = try_planar_face_area_shoelace(brep, face, face.normal)
                .or_else(|| try_planar_face_exact_contour_area(brep, face, face.normal));
            if std::env::var("RCAD_AREA_DEBUG").is_ok() {
                eprintln!("[AREA] plane face n_wires={} outer_edges={} inner={} shoelace_or_exact={:?}",
                    face.inner_wires.len() + 1, face.outer_wire.edges.len(),
                    face.inner_wires.iter().map(|w| w.edges.len()).collect::<Vec<_>>(),
                    r);
            }
            r
        }

        Surface3::Cylinder(cyl) =>
            try_cylinder_trimmed_face_area(cyl, brep, face, face_flat_idx),

        Surface3::Cone(cone) =>
            try_cone_trimmed_face_area(cone, brep, face),

        Surface3::Sphere(s) =>
            try_spherical_polygon_great_circle_area(s, brep, face)
            .or_else(|| face_surface_area_gauss(brep, face, face_flat_idx)),

        Surface3::Torus(_) | Surface3::BSpline(_) | Surface3::Ellipsoid(_) =>
            face_surface_area_gauss(brep, face, face_flat_idx)
            .or_else(|| param_rect_area_cross(surf, 0.0, 1.0, 0.0, 1.0)),

        _ => {
            let [u0, u1, v0, v1] = surf.default_domain();
            param_rect_area_cross(surf, u0, u1, v0, v1)
        }
    }
}

// ══════════════════════════════════════════════════════════════════════════
// GL integration table
// ══════════════════════════════════════════════════════════════════════════

struct GLTable {
    points: usize,
    x: Vec<f64>,
    w: Vec<f64>,
}

fn gl_table(order: usize) -> &'static GLTable {
    use std::sync::OnceLock;
    static TABLES: OnceLock<Vec<GLTable>> = OnceLock::new();
    let tables = TABLES.get_or_init(|| {
        vec![
            GLTable { points: 2, x: vec![-0.5773502691896257, 0.5773502691896257], w: vec![1.0, 1.0] },
            GLTable { points: 3, x: vec![-0.7745966692414834, 0.0, 0.7745966692414834], w: vec![0.5555555555555556, 0.8888888888888888, 0.5555555555555556] },
            GLTable { points: 4, x: vec![-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526], w: vec![0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538] },
            GLTable { points: 5, x: vec![-0.906179845938664, -0.5384693101056831, 0.0, 0.5384693101056831, 0.906179845938664], w: vec![0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891] },
            GLTable { points: 7, x: vec![-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585], w: vec![0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697] },
        ]
    });
    let idx = tables.iter().position(|t| t.points == order).unwrap_or(0);
    &tables[idx]
}

fn gl_s_integration_order(surf: &Surface3) -> (usize, usize) {
    match surf {
        Surface3::Plane(_) => (2, 2),
        Surface3::Cylinder(_) => (4, 2),
        Surface3::Cone(_) => (4, 2),
        Surface3::Sphere(_) => (7, 7),
        Surface3::Torus(_) => (7, 5),
        _ => (7, 7),
    }
}

fn gl_s_u_subs(surf: &Surface3) -> usize {
    match surf {
        Surface3::BSpline(b) => b.knots_u.len().max(8),
        _ => 8,
    }
}

fn gl_s_v_subs(surf: &Surface3) -> usize {
    match surf {
        Surface3::BSpline(b) => b.knots_v.len().max(8),
        _ => 8,
    }
}

fn gl_u_knots(surf: &Surface3, u0: f64, u1: f64) -> Vec<f64> {
    match surf {
        Surface3::BSpline(b) => {
            let mut k = b.knots_u.clone();
            if k.first().map_or(true, |&k0| (k0 - u0).abs() > 1e-12) { k.insert(0, u0); }
            if k.last().map_or(true, |&k1| (k1 - u1).abs() > 1e-12) { k.push(u1); }
            k
        }
        _ => vec![u0, u1],
    }
}

fn gl_v_knots(surf: &Surface3, v0: f64, v1: f64) -> Vec<f64> {
    match surf {
        Surface3::BSpline(b) => {
            let mut k = b.knots_v.clone();
            if k.first().map_or(true, |&k0| (k0 - v0).abs() > 1e-12) { k.insert(0, v0); }
            if k.last().map_or(true, |&k1| (k1 - v1).abs() > 1e-12) { k.push(v1); }
            k
        }
        _ => vec![v0, v1],
    }
}

fn surface_normal_jacobian(surf: &Surface3, u: f64, v: f64) -> f64 {
    let (_p, du, dv) = surf.derivatives(u, v);
    du.cross(dv).length()
}

fn face_surface_area_gauss(brep: &BRep, face: &Face, fi: usize) -> Option<f64> {
    let surf_idx = brep.tshapes.get(fi).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
    })?;
    let surf = &surf_idx;
    let [u0, u1, v0, v1] = curved_face_uv_domain(brep, face, surf)?;
    let (order_u, order_v) = gl_s_integration_order(surf);
    let n_u = gl_s_u_subs(surf);
    let n_v = gl_s_v_subs(surf);
    let tab = gl_table(order_u);
    let tab_v = gl_table(order_v);

    let u_knots = gl_u_knots(surf, u0, u1);
    let v_knots = gl_v_knots(surf, v0, v1);

    let mut total = 0.0;
    for ui in 0..u_knots.len() - 1 {
        let u_lo = u_knots[ui]; let u_hi = u_knots[ui + 1];
        let u_mid = (u_lo + u_hi) * 0.5; let u_half = (u_hi - u_lo) * 0.5;
        for vi in 0..v_knots.len() - 1 {
            let v_lo = v_knots[vi]; let v_hi = v_knots[vi + 1];
            let v_mid = (v_lo + v_hi) * 0.5; let v_half = (v_hi - v_lo) * 0.5;

            let mut sum = 0.0;
            for gi in 0..order_u {
                let u = u_mid + u_half * tab.x[gi];
                for gj in 0..order_v {
                    let v = v_mid + v_half * tab_v.x[gj];
                    let jac = surface_normal_jacobian(surf, u, v);
                    sum += tab.w[gi] * tab_v.w[gj] * jac;
                }
            }
            total += sum * u_half * v_half;
        }
    }
    Some(total)
}

fn param_rect_area_cross(surf: &Surface3, u0: f64, u1: f64, v0: f64, v1: f64) -> Option<f64> {
    if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite()
        || (u1 - u0).abs() < 1e-14 || (v1 - v0).abs() < 1e-14 {
        return None;
    }
    const N: usize = 32;
    let du = (u1 - u0) / N as f64;
    let dv = (v1 - v0) / N as f64;
    let mut area = 0.0;
    for i in 0..N {
        let u = u0 + (i as f64 + 0.5) * du;
        for j in 0..N {
            let v = v0 + (j as f64 + 0.5) * dv;
            area += surface_normal_jacobian(surf, u, v) * du * dv;
        }
    }
    Some(area)
}

fn curved_face_uv_domain(brep: &BRep, face: &Face, surf: &Surface3) -> Option<[f64; 4]> {
    if let Some(face_ts) = brep.tshapes.get(brep.tshapes.iter().position(|ts| {
        matches!(**ts, topods::TShape::Face(_))
    })?) {
        if let topods::TShape::Face(fd) = &**face_ts {
            if let Some(domain) = fd.uv_domain { return Some(domain); }
        }
    }
    let pts = sample_wire_polyline_3d(brep, &face.outer_wire);
    if pts.is_empty() { return Some(surf.default_domain()); }
    let proj: Vec<DVec2> = pts.iter().map(|&p| {
        let r = closest_point_on_surface(surf, p, 64);
        DVec2::new(r.params.0, r.params.1)
    }).collect();
    let mut u0 = f64::INFINITY; let mut u1 = f64::NEG_INFINITY;
    let mut v0 = f64::INFINITY; let mut v1 = f64::NEG_INFINITY;
    for p in &proj { u0 = u0.min(p.x); u1 = u1.max(p.x); v0 = v0.min(p.y); v1 = v1.max(p.y); }
    if !u0.is_finite() { return Some(surf.default_domain()); }
    let margin_u = (u1 - u0) * 0.02 + 1e-4;
    let margin_v = (v1 - v0) * 0.02 + 1e-4;
    Some([u0 - margin_u, u1 + margin_u, v0 - margin_v, v1 + margin_v])
}

// ══════════════════════════════════════════════════════════════════════════
// Planar face area
// ══════════════════════════════════════════════════════════════════════════

fn axis_aligned_world_plane_uv_axes(n: DVec3) -> Option<[usize; 2]> {
    let a = n.abs();
    if a.x > 1.0 - 2e-3 { Some([1, 2]) } else if a.y > 1.0 - 2e-3 { Some([0, 2]) }
    else if a.z > 1.0 - 2e-3 { Some([0, 1]) } else { None }
}

fn bbox2d_components(uv: &[(f64, f64)]) -> Option<(f64, f64, f64, f64)> {
    if uv.is_empty() { return None; }
    let mut u0 = f64::INFINITY; let mut u1 = f64::NEG_INFINITY;
    let mut v0 = f64::INFINITY; let mut v1 = f64::NEG_INFINITY;
    for &(u, v) in uv { u0 = u0.min(u); u1 = u1.max(u); v0 = v0.min(v); v1 = v1.max(v); }
    Some((u0, u1, v0, v1))
}

fn try_axis_aligned_world_rect_plane_area(brep: &BRep, face: &Face, face_normal: DVec3) -> Option<f64> {
    let n = face_normal.normalize_or_zero();
    if n.length_squared() < 1e-24 { return None; }
    let [i, j] = axis_aligned_world_plane_uv_axes(n)?;
    let pos_tol = (1e-7 * brep.bounding_box().map(|[mn, mx]| (mx - mn).length()).unwrap_or(1.0)).max(1e-9);
    let vu = outer_wire_unique_vertex_uvs(brep, &face.outer_wire, i, j, pos_tol);
    let mut outer_ordered = sample_wire_polyline_3d(brep, &face.outer_wire);
    trim_almost_closed_polyline(&mut outer_ordered, 1e-5);
    let a_loop = if outer_ordered.len() >= 3 {
        let uv_ord: Vec<(f64, f64)> = outer_ordered.iter().map(|p| (p[i], p[j])).collect();
        polygon_area_2d_xy(&uv_ord)
    } else { 0.0 };
    if (3..=4096).contains(&vu.len()) {
        let hull = convex_hull_2d_monotone(vu.clone());
        if hull.len() >= 3 {
            let a_hull = polygon_area_2d_xy(&hull);
            if a_hull > 1e-18 {
                const REL: f64 = 1e-4;
                const AGREE_REL: f64 = 2e-3;
                let scale = a_hull.max(a_loop).max(1.0);
                let abs_eps = 1e-9 * scale;
                if a_loop > 1e-18 && a_hull > a_loop * (1.0 + REL) + abs_eps {
                    let uv_vert = outer_wire_ordered_vertex_uvs(brep, &face.outer_wire, i, j, pos_tol);
                    if uv_vert.len() >= 3 {
                        let a_vert = polygon_area_2d_xy(&uv_vert);
                        let scale_agree = a_vert.max(a_loop).max(1.0);
                        const LOOP_FRAC_OF_HULL_MIN: f64 = 0.81;
                        if (a_vert - a_loop).abs() <= AGREE_REL * scale_agree + abs_eps && a_loop + abs_eps >= a_hull * LOOP_FRAC_OF_HULL_MIN { return Some(a_loop); }
                        if (a_vert - a_loop).abs() <= AGREE_REL * scale_agree + abs_eps && a_vert + abs_eps < a_hull * 0.6 { return Some(a_vert.max(0.0)); }
                        const MIN_HULL_ABS_VERT_FALLBACK: f64 = 15000.0;
                        const VERT_OVER_LOOP_REL: f64 = 0.02;
                        if a_hull >= MIN_HULL_ABS_VERT_FALLBACK && a_vert + abs_eps < a_hull * (1.0 - REL) && a_vert > 1e-18 && a_vert > a_loop * (1.0 + VERT_OVER_LOOP_REL) + abs_eps { return Some(a_vert.max(0.0)); }
                    }
                }
                if a_loop > 1e-18 && a_loop + abs_eps >= a_hull { return Some(a_loop); }
                return Some(a_hull);
            }
        }
    }
    let mut outer = outer_ordered;
    trim_almost_closed_polyline(&mut outer, 1e-5);
    if outer.len() < 3 { return None; }
    let uv: Vec<(f64, f64)> = outer.iter().map(|p| (p[i], p[j])).collect();
    let (u0, u1, v0, v1) = bbox2d_components(&uv)?;
    let w = u1 - u0; let h = v1 - v0;
    if !(w > 1e-18 && h > 1e-18) { return None; }
    let scale = w.max(h).max(1.0);
    let eps = (1e-5 * scale).max(1e-9);
    if uv.iter().all(|&(u, v)| (u - u0).abs() <= eps || (u1 - u).abs() <= eps || (v - v0).abs() <= eps || (v1 - v).abs() <= eps) {
        Some((w * h).max(0.0))
    } else { None }
}

fn wire_edge_endpoint_3d(brep: &BRep, we: &WireEdge) -> Option<DVec3> {
    let edge = brep.flat_edges().get(we.idx).copied()?;
    if let Some(curve) = brep.tshapes.get(we.idx).and_then(|ts| {
        if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
    }) {
        let range = brep.tshapes.get(we.idx).and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None }
        }).unwrap_or_else(|| curve.default_domain());
        let t = if we.forward { range[0] } else { range[1] };
        return Some(curve.point_at(t));
    }
    let vidx = if we.forward { edge.0 } else { edge.1 };
    Some(brep.vertex_point(vidx)?)
}

fn outer_wire_ordered_vertex_uvs(brep: &BRep, wire: &Wire, i: usize, j: usize, pos_tol: f64) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::new();
    for we in &wire.edges {
        let Some(p) = wire_edge_endpoint_3d(brep, we) else { continue; };
        let uv = (p[i], p[j]);
        if let Some(&prev) = out.last() { if (uv.0 - prev.0).abs() <= pos_tol && (uv.1 - prev.1).abs() <= pos_tol { continue; } }
        out.push(uv);
    }
    trim_almost_closed_uv_chain(&mut out, pos_tol);
    out
}

fn trim_almost_closed_uv_chain(uvs: &mut Vec<(f64, f64)>, pos_tol: f64) {
    if uvs.len() >= 2 {
        let (u0, v0) = uvs[0]; let (u1, v1) = uvs[uvs.len() - 1];
        if (u0 - u1).abs() <= pos_tol && (v0 - v1).abs() <= pos_tol { uvs.pop(); }
    }
}

fn outer_wire_unique_vertex_uvs(brep: &BRep, wire: &Wire, i: usize, j: usize, pos_tol: f64) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::new();
    for we in &wire.edges {
        let flat_edges = brep.flat_edges();
        let Some(edge) = flat_edges.get(we.idx) else { continue; };
        let pts: [DVec3; 2] = if let Some(curve) = brep.tshapes.get(we.idx).and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
        }) {
            let range = brep.tshapes.get(we.idx).and_then(|ts| {
                if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None }
            }).unwrap_or_else(|| curve.default_domain());
            [curve.point_at(range[0]), curve.point_at(range[1])]
        } else {
            let p0 = brep.vertex_point(edge.0).unwrap_or(DVec3::ZERO);
            let p1 = brep.vertex_point(edge.1).unwrap_or(DVec3::ZERO);
            [p0, p1]
        };
        for &p in &pts {
            let uv = (p[i], p[j]);
            let already = out.iter().any(|&(u2, v2)| (uv.0 - u2).abs() <= pos_tol && (uv.1 - v2).abs() <= pos_tol);
            if !already { out.push(uv); }
        }
    }
    out
}

fn try_planar_face_exact_contour_area(brep: &BRep, face: &Face, face_normal: DVec3) -> Option<f64> {
    let (ux, uy) = local_basis_from_normal(face_normal);
    let n_edges = face.outer_wire.edges.len();
    // Single full circle
    if n_edges == 1 {
        if let Some(we) = face.outer_wire.edges.first() {
            if let Some(curve_idx) = brep.tshapes.get(we.idx).and_then(|ts| {
                if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
            }) {
                if let Curve3::Circle(c) = curve_idx {
                    if let Some(range) = brep.tshapes.get(we.idx).and_then(|ts| {
                        if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None }
                    }) {
                        let theta = (range[1] - range[0]).abs();
                        if (theta - 2.0 * PI).abs() < 1e-12 { return Some(PI * c.radius * c.radius); }
                    }
                }
            }
        }
        return None;
    }
    // Two circles forming a full circle (cylinder seam split)
    if n_edges == 2 {
        let mut radii = [0.0f64; 2]; let mut centers = [DVec3::ZERO; 2]; let mut spans = [0.0f64; 2]; let mut n_circle = 0u32;
        for (i, we) in face.outer_wire.edges.iter().enumerate() {
            if let Some(ci) = brep.tshapes.get(we.idx).and_then(|ts| {
                if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
            }) {
                if let Curve3::Circle(c) = ci {
                    if i < 2 { radii[i] = c.radius; centers[i] = c.center; }
                    if let Some(r) = brep.tshapes.get(we.idx).and_then(|ts| {
                        if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None }
                    }) { if i < 2 { spans[i] = (r[1] - r[0]).abs(); } }
                    n_circle += 1;
                }
            }
        }
        if n_circle == 2 {
            let total_theta = spans[0] + spans[1];
            if (total_theta - 2.0 * PI).abs() < 1e-10 && (centers[0] - centers[1]).length_squared() < 1e-12 && (radii[0] - radii[1]).abs() < 1e-12 {
                return Some(PI * radii[0] * radii[0]);
            }
        }
    }
    if n_edges < 3 { return None; }

    let mut edges = Vec::with_capacity(n_edges);
    let first_we = &face.outer_wire.edges[0];
    let first_e = brep.flat_edges().get(first_we.idx).copied()?;
    let first_vi = if first_we.forward { first_e.0 } else { first_e.1 };
    let pivot = brep.vertex_point(first_vi).unwrap_or(DVec3::ZERO);

    for we in &face.outer_wire.edges {
        let ei = we.idx;
        let edge = brep.flat_edges().get(ei).copied()?;
        let curve_idx = brep.tshapes.get(ei).and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
        })?;
        let range = brep.tshapes.get(ei).and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { Some(ed.range) } else { None }
        }).unwrap_or([0.0, 1.0]);
        let (v_start, v_end) = if we.forward { (edge.0, edge.1) } else { (edge.1, edge.0) };
        let p_start = brep.vertex_point(v_start).unwrap_or(DVec3::ZERO);
        let p_end = brep.vertex_point(v_end).unwrap_or(DVec3::ZERO);
        let start_2d = DVec2::new((p_start - pivot).dot(ux), (p_start - pivot).dot(uy));
        let end_2d = DVec2::new((p_end - pivot).dot(ux), (p_end - pivot).dot(uy));

        match curve_idx {
            Curve3::Line(_) => edges.push(EdgeArcInfo { is_arc: false, radius: 0.0, theta: 0.0, center_2d: DVec2::ZERO, sign: 0.0, start_2d, end_2d }),
            Curve3::Circle(c) => {
                let theta = (range[1] - range[0]).abs();
                if theta < 1e-15 || theta > 2.0 * PI + 1e-12 { return None; }
                let center_2d = DVec2::new((c.center - pivot).dot(ux), (c.center - pivot).dot(uy));
                let trav_3d = p_end - p_start;
                let left_dir = face_normal.cross(trav_3d);
                let sign = if (c.center - p_start).dot(left_dir) > 0.0 { 1.0 } else { -1.0 };
                edges.push(EdgeArcInfo { is_arc: true, radius: c.radius, theta, center_2d, sign, start_2d, end_2d });
            }
            _ => return None,
        }
    }

    let n = edges.len();
    let mut shoelace = 0.0;
    for i in 0..n { let s = edges[i].start_2d; let e = edges[i].end_2d; shoelace += s.x * e.y - e.x * s.y; }
    let mut total = shoelace.abs() * 0.5;
    let shoelace_raw = total;
    let mut seg_total = 0.0; let mut n_neg = 0u32; let mut n_pos = 0u32;
    for edge in &edges {
        if edge.is_arc {
            let t = edge.theta;
            let seg = edge.radius * edge.radius * (t - t.sin()) * 0.5;
            if edge.sign < 0.0 { n_neg += 1; seg_total -= seg; } else { n_pos += 1; seg_total += seg; }
        }
    }
    total = shoelace_raw;
    if (n_neg == 0 || n_pos == 0) && seg_total < 0.0 { total -= seg_total; } else { total += seg_total; }

    if total > 0.0 && total.is_finite() {
        for w in &face.inner_wires {
            if w.edges.len() < 3 { continue; }
            let hole_pts: Vec<DVec3> = w.edges.iter().filter_map(|we| {
                let e = brep.flat_edges().get(we.idx).copied()?;
                let vi = if we.forward { e.0 } else { e.1 };
                brep.vertex_point(vi)
            }).collect();
            if hole_pts.len() < 3 { continue; }
            let mut a_hole = 0.0;
            for i in 1..hole_pts.len().saturating_sub(1) {
                let d1 = hole_pts[i] - hole_pts[0];
                let d2 = hole_pts[i + 1] - hole_pts[0];
                a_hole += d1.cross(d2).dot(face_normal).abs() * 0.5;
            }
            total -= a_hole.min(total);
        }
        return Some(total);
    }
    None
}

struct EdgeArcInfo {
    is_arc: bool, radius: f64, theta: f64, center_2d: DVec2, sign: f64, start_2d: DVec2, end_2d: DVec2,
}

fn try_planar_face_area_shoelace(brep: &BRep, face: &Face, face_normal: DVec3) -> Option<f64> {
    if face.inner_wires.is_empty() {
        if let Some(a_rect) = try_axis_aligned_world_rect_plane_area(brep, face, face_normal) {
            let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
            trim_almost_closed_polyline(&mut outer, 1e-5);
            if outer.len() >= 3 {
                let (ux, uy) = local_basis_from_normal(face_normal);
                if let Some(pivot) = outer.first().copied() {
                    let a_shoe = polygon_area_2d_projected(&outer, pivot, ux, uy).abs();
                    const REL: f64 = 1e-5;
                    let scale = a_rect.max(a_shoe).max(1.0);
                    let abs_eps = 1e-9 * scale;
                    if a_rect > a_shoe * (1.0 + REL) + abs_eps {
                        let ratio = a_shoe / a_rect.max(1e-12);
                        if ratio >= 0.40 || a_shoe + abs_eps >= a_rect * 0.65 { return Some(a_shoe.max(0.0)); }
                        return Some(a_rect.max(0.0));
                    }
                    if a_shoe > a_rect * (1.0 + REL) + abs_eps { return Some(a_rect); }
                    return Some(a_rect);
                }
            }
            return Some(a_rect);
        }
    }
    let mut outer = sample_wire_polyline_3d(brep, &face.outer_wire);
    trim_almost_closed_polyline(&mut outer, 1e-5);
    if outer.len() < 3 { return None; }
    let (ux, uy) = local_basis_from_normal(face_normal);
    let pivot = outer.first().copied()?;
    let mut a = polygon_area_2d_projected(&outer, pivot, ux, uy).abs();

    // Convex hull cross-check for scrambled wire ordering
    if a * 1e-7 < outer.len() as f64 {
        if let Some(hull_a) = try_boundary_convex_hull_area(brep, &face.outer_wire, pivot, ux, uy) {
            if hull_a > 1e-12 && a < 0.6 * hull_a { return None; }
        }
    }
    if a < 1e-12 {
        let mut bbox_min = outer[0]; let mut bbox_max = outer[0];
        for p in &outer { bbox_min = bbox_min.min(*p); bbox_max = bbox_max.max(*p); }
        let bbox_diag = (bbox_max - bbox_min).length();
        if bbox_diag > 1e-10 && a < 1e-12 * bbox_diag * bbox_diag { return None; }
        else if bbox_diag <= 1e-10 { return None; }
    }
    Some(a.max(0.0))
}

// ══════════════════════════════════════════════════════════════════════════
// Cylinder face area
// ══════════════════════════════════════════════════════════════════════════

fn cylinder_uv_area_gl(uvs: &[DVec2]) -> Option<f64> {
    const NU: usize = 60;
    let n = uvs.len();
    if n < 3 { return None; }
    let two_pi = 2.0 * PI;
    let mut poly = Vec::with_capacity(n);
    poly.push(uvs[0]);
    for i in 1..n { let du = short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x); poly.push(DVec2::new(poly[i - 1].x + du, uvs[i].y)); }
    let v_range_at = |u: f64| -> (f64, f64) {
        let mut v_lo = f64::INFINITY; let mut v_hi = f64::NEG_INFINITY;
        for i in 0..n {
            let j = (i + 1) % n;
            let (u1, v1) = (poly[i].x, poly[i].y); let (u2, v2) = (poly[j].x, poly[j].y);
            if u1 == u2 { continue; }
            let (u_lo_e, u_hi_e, v_lo_e, v_hi_e) = if u1 < u2 { (u1, u2, v1, v2) } else { (u2, u1, v2, v1) };
            if u < u_lo_e - 1e-12 || u > u_hi_e + 1e-12 { continue; }
            let t = if u_hi_e > u_lo_e { (u - u_lo_e) / (u_hi_e - u_lo_e) } else { 0.0 };
            let v = v_lo_e + t * (v_hi_e - v_lo_e);
            v_lo = v_lo.min(v); v_hi = v_hi.max(v);
        }
        (v_lo, v_hi)
    };
    let u_min = poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let u_max = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let u_range = u_max - u_min;
    if !u_min.is_finite() || u_range < 1e-14 { return None; }
    let u_lo = u_min;
    let u_hi = if u_range > two_pi - 0.1 { u_min + two_pi } else { u_max };
    let du = (u_hi - u_lo) / NU as f64;
    const GL_NEG: f64 = -0.5773502691896257; const GL_POS: f64 = 0.5773502691896257;
    let gl_pts = [GL_NEG, GL_POS];
    let mut total = 0.0;
    for i in 0..NU {
        let u_mid = u_lo + (i as f64 + 0.5) * du;
        for &gu in &gl_pts {
            let u = u_mid + gu * du * 0.5;
            let (v_lo, v_hi) = v_range_at(u);
            if v_lo.is_finite() && v_hi > v_lo + 1e-14 { total += (v_hi - v_lo) * du * 0.5; }
        }
    }
    Some(total)
}

fn cylinder_gl_uv_area(uvs: &[DVec2]) -> Option<f64> {
    const N: usize = 60;
    if uvs.len() < 3 { return None; }
    let n = uvs.len();
    const TWO_PI: f64 = 2.0 * PI;
    let unwrapped: Vec<f64> = {
        let mut o = Vec::with_capacity(n); o.push(uvs[0].x);
        for i in 1..n { o.push(o[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x)); }
        o
    };
    let mut umin = f64::INFINITY; let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY; let mut vmax = f64::NEG_INFINITY;
    for (i, uv) in uvs.iter().enumerate() { umin = umin.min(unwrapped[i]); umax = umax.max(unwrapped[i]); vmin = vmin.min(uv.y); vmax = vmax.max(uv.y); }
    if !umin.is_finite() || (umax - umin) < 1e-14 || (vmax - vmin) < 1e-14 { return None; }
    let u_range = umax - umin;
    let u_lo = umin;
    let u_hi = if u_range > TWO_PI - 0.1 { umin + TWO_PI } else { umax };
    let v_lo = vmin; let v_hi = vmax;
    let du = (u_hi - u_lo) / N as f64; let dv = (v_hi - v_lo) / N as f64;
    let cell_area = du * dv / 4.0;
    const GL_NEG: f64 = -0.5773502691896257; const GL_POS: f64 = 0.5773502691896257;
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
                    if winding_number_2d(uvs, DVec2::new(u_mod, v)) != 0 { n_hit += 1; }
                }
            }
            if n_hit > 0 { total += cell_area * n_hit as f64; }
        }
    }
    Some(total)
}

fn cylinder_outer_wire_uv_shoelace_area(brep: &BRep, cyl: &CylindricalSurface, face: &Face) -> Option<f64> {
    let mut pts_3d = sample_wire_polyline_3d_with_n(brep, &face.outer_wire, 512);
    trim_almost_closed_polyline(&mut pts_3d, 1e-5);
    if pts_3d.len() < 3 { return None; }
    let n = pts_3d.len();
    let surf = Surface3::Cylinder(*cyl);
    let uvs: Vec<DVec2> = pts_3d.iter().map(|&p| {
        let proj = closest_point_on_surface(&surf, p, 256);
        DVec2::new(proj.params.0, proj.params.1)
    }).collect();
    let mut area2 = 0.0_f64;
    let unwrapped: Vec<f64> = {
        let mut o = Vec::with_capacity(n); o.push(uvs[0].x);
        for i in 1..n { o.push(o[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x)); }
        o
    };
    for i in 0..n { let j = if i + 1 < n { i + 1 } else { 0 }; area2 += unwrapped[i] * uvs[j].y - unwrapped[j] * uvs[i].y; }
    let uv_area = area2.abs() * 0.5;
    let gl_area = cylinder_uv_area_gl(&uvs).unwrap_or(0.0);
    if gl_area > 0.0 { Some(gl_area * cyl.radius) } else { Some(uv_area * cyl.radius) }
}

fn try_cylinder_trimmed_face_area(cyl: &CylindricalSurface, brep: &BRep, face: &Face, face_flat_idx: usize) -> Option<f64> {
    // Fast path: rectangular UV patch via 2 Lines + 2 Circles
    if face.inner_wires.is_empty() && face.outer_wire.edges.len() == 4 {
        let mut n_lines = 0u32; let mut n_circles = 0u32;
        let mut circle_centers = Vec::new();
        let mut edge_curve_indices = Vec::new();
        let mut valid = true;
        for we in &face.outer_wire.edges {
            if let Some(ci) = brep.tshapes.get(we.idx).and_then(|ts| {
                if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
            }) {
                match ci {
                    Curve3::Line(_) => { n_lines += 1; edge_curve_indices.push(we.idx); }
                    Curve3::Circle(c) => { n_circles += 1; if circle_centers.len() < 2 { circle_centers.push(c.center); } }
                    _ => { valid = false; break; }
                }
            } else { valid = false; break; }
        }
        if valid && n_lines == 2 && n_circles == 2 && circle_centers.len() == 2 {
            let axis = cyl.axis;
            let v0 = (circle_centers[0] - cyl.origin).dot(axis);
            let v1 = (circle_centers[1] - cyl.origin).dot(axis);
            let dv = (v1 - v0).abs();
            let du = 'du: {
                if let Some(range) = brep.tshapes.get(face_flat_idx).and_then(|ts| {
                    if let topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None }
                }) {
                    let du_r = (range[1] - range[0]).abs(); let dv_r = (range[3] - range[2]).abs();
                    if du_r > 1e-14 && dv_r > 1e-14 { break 'du du_r; }
                }
                if edge_curve_indices.len() == 2 && edge_curve_indices[0] == edge_curve_indices[1] { break 'du (2.0 * PI); }
                let x_ax = crate::geom::any_perpendicular(axis);
                let y_ax = axis.cross(x_ax).normalize();
                let mut u_vals = Vec::new();
                for &ci in &edge_curve_indices {
                    if let Some(Curve3::Line(_)) = brep.tshapes.get(ci).and_then(|ts| {
                        if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
                    }) {
                        for we in &face.outer_wire.edges {
                            if we.idx == ci {
                                if let Some(edge) = brep.flat_edges().get(we.idx) {
                                    let vi = if we.forward { edge.0 } else { edge.1 };
                                    if let Some(v) = brep.vertex_point(vi) {
                                        let d = v - cyl.origin;
                                        u_vals.push(d.dot(y_ax).atan2(d.dot(x_ax)));
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                if u_vals.len() == 2 {
                    let du_norm = (u_vals[1] - u_vals[0]).rem_euclid(2.0 * PI);
                    if du_norm > 1e-14 { break 'du du_norm; }
                }
                2.0 * PI
            };
            let r0 = (circle_centers[0] - cyl.origin).cross(axis).length();
            let r1 = (circle_centers[1] - cyl.origin).cross(axis).length();
            if du > 1e-14 && dv > 1e-14 && du.is_finite() && dv.is_finite() && r0 < 1e-8 && r1 < 1e-8 {
                return Some(cyl.radius * du * dv);
            }
        }
    }
    // General case: UV projection + GL integration
    let uv_area = cylinder_outer_wire_uv_shoelace_area(brep, cyl, face)?;
    // Subtract inner wire areas
    let mut total = uv_area;
    for w in &face.inner_wires {
        let area = cylinder_outer_wire_uv_shoelace_area(brep, cyl, &Face {
            outer_wire: w.clone(), inner_wires: vec![], normal: face.normal,
            sample_point: face.sample_point, triangles: Vec::new(),
            mesh_dirty: true, surface_idx: None,
        })?;
        total -= area.min(total);
    }
    Some(total.max(0.0))
}

// ══════════════════════════════════════════════════════════════════════════
// Cone face area
// ══════════════════════════════════════════════════════════════════════════

fn try_cone_trimmed_face_area(cone: &ConicalSurface, brep: &BRep, face: &Face) -> Option<f64> {
    let mut pts_3d = sample_wire_polyline_3d_with_n(brep, &face.outer_wire, 512);
    trim_almost_closed_polyline(&mut pts_3d, 1e-5);
    if pts_3d.len() < 3 { return None; }
    let surf = Surface3::Cone(*cone);
    let uvs: Vec<DVec2> = pts_3d.iter().map(|&p| {
        let proj = closest_point_on_surface(&surf, p, 256);
        DVec2::new(proj.params.0, proj.params.1)
    }).collect();
    let n = uvs.len();
    let mut unwrapped = Vec::with_capacity(n);
    if n > 0 { unwrapped.push(uvs[0].x); }
    for i in 1..n { unwrapped.push(unwrapped[i - 1] + short_delta_on_circle_01(uvs[i - 1].x, uvs[i].x)); }

    // Compute uv area: R(v) = cone.radius + v * tan(α) per the conical surface parameterization
    // Scale factor for area: R(v) / cos(α) — accounts for the slope of the cone.
    let cos_half = cone.half_angle_rad.cos();
    let factor = if cos_half > 1e-12 { 1.0 / cos_half } else { return None; };

    let gl_area = {
        const N: usize = 60;
        let mut umin = f64::INFINITY; let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY; let mut vmax = f64::NEG_INFINITY;
        for (i, uv) in uvs.iter().enumerate() { umin = umin.min(unwrapped[i]); umax = umax.max(unwrapped[i]); vmin = vmin.min(uv.y); vmax = vmax.max(uv.y); }
        if !umin.is_finite() || (umax - umin) < 1e-14 || (vmax - vmin) < 1e-14 { return None; }
        let du = (umax - umin) / N as f64; let dv = (vmax - vmin) / N as f64;
        let cell_area = du * dv / 4.0;
        const GL_NEG: f64 = -0.5773502691896257; const GL_POS: f64 = 0.5773502691896257;
        let gl_pts = [GL_NEG, GL_POS];
        let mut total = 0.0;
        for i in 0..N {
            let u_mid = umin + (i as f64 + 0.5) * du;
            for j in 0..N {
                let v_mid = vmin + (j as f64 + 0.5) * dv;
                let mut n_hit = 0u32;
                for &gu in &gl_pts {
                    let u = u_mid + gu * du * 0.5;
                    for &gv in &gl_pts {
                        let v = v_mid + gv * dv * 0.5;
                        if winding_number_2d(&uvs, DVec2::new(u, v)) != 0 { n_hit += 1; }
                    }
                }
                if n_hit > 0 { total += cell_area * n_hit as f64; }
            }
        }
        total * factor
    };

    Some(gl_area)
}

// ══════════════════════════════════════════════════════════════════════════
// Sphere face area
// ══════════════════════════════════════════════════════════════════════════

fn try_spherical_polygon_great_circle_area(s: &SphericalSurface, brep: &BRep, face: &Face) -> Option<f64> {
    if !face.inner_wires.is_empty() { return None; }
    let n_edges = face.outer_wire.edges.len();
    if n_edges < 3 { return None; }
    let tol = 1e-10;
    let mut verts: Vec<DVec3> = Vec::with_capacity(n_edges + 1);
    for we in &face.outer_wire.edges {
        let ei = we.idx;
        let edge = brep.flat_edges().get(ei).copied()?;
        let curve_idx = brep.tshapes.get(ei).and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { ed.curve.as_ref() } else { None }
        })?;
        match curve_idx {
            Curve3::Circle(c) => { if (c.center - s.center).length() > tol { return None; } }
            _ => return None,
        }
        let vi = if we.forward { edge.0 } else { edge.1 };
        let pt = brep.vertex_point(vi).unwrap_or(DVec3::ZERO);
        if verts.is_empty() || (pt - *verts.last()?).length() > tol { verts.push(pt); }
    }
    if (verts.first()? - verts.last()?).length() > tol { verts.push(*verts.first()?); }
    let n = verts.len() - 1;
    if n < 3 { return None; }
    let mut sum_angles = 0.0;
    for i in 0..n {
        let v_prev = if i > 0 { verts[i - 1] } else { verts[n - 1] };
        let v_curr = verts[i];
        let v_next = verts[i + 1];
        let v_hat = v_curr.normalize();
        let t_in = (v_prev - v_hat * v_prev.dot(v_hat)).normalize();
        let t_out = (v_next - v_hat * v_next.dot(v_hat)).normalize();
        let cos_theta = t_in.dot(t_out).clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        let cross_sign = t_in.cross(t_out).dot(v_hat);
        let interior = if cross_sign < 0.0 { PI - theta } else { PI + theta };
        sum_angles += interior;
    }
    let r2 = s.radius * s.radius;
    let full = 4.0 * PI * r2;
    let mut area = r2 * (sum_angles - (n as f64 - 2.0) * PI);
    area = area.abs();
    if area > full * 0.5 { area = full - area; }
    if let Some(sp) = face.sample_point {
        if !point_in_spherical_polygon_3d(&verts[..n], sp) { area = full - area; }
    }
    if area > 0.0 && area <= full + 1e-12 { Some(area) } else { None }
}

pub fn try_spherical_uv_masked_raster(
    s: &SphericalSurface, brep: &BRep, face: &Face, face_flat_idx: usize, face_normal: DVec3,
) -> Option<Vec<[DVec3; 3]>> {
    let ctx = spherical_holed_uv_mask_setup(s, brep, face)?;
    const N: usize = 30;
    const GL_NEG: f64 = -0.5773502691896257; const GL_POS: f64 = 0.5773502691896257;
    let gl_pts = [GL_NEG, GL_POS];
    let umin = ctx.umin; let umax = ctx.umax; let vmin = ctx.vmin; let vmax = ctx.vmax;
    let du = (umax - umin) / N as f64; let dv = (vmax - vmin) / N as f64;
    let mut tris = Vec::new();
    for i in 0..N {
        let u_mid = umin + (i as f64 + 0.5) * du;
        for j in 0..N {
            let v_mid = vmin + (j as f64 + 0.5) * dv;
            for &gu in &gl_pts {
                let u = u_mid + gu * du * 0.5;
                for &gv in &gl_pts {
                    let v = v_mid + gv * dv * 0.5;
                    if v < 0.0 || v > PI { continue; }
                    let p3d = s.point_at(u, v);
                    if !point_in_spherical_polygon_3d(&ctx.outer_3d, p3d) { continue; }
                    if ctx.inner_3d.iter().any(|h3d| point_in_spherical_polygon_3d(h3d, p3d)) { continue; }
                    tris.push(orient_tri([p3d, s.point_at(u + du, v), s.point_at(u, v + dv)], face_normal));
                    tris.push(orient_tri([s.point_at(u + du, v), s.point_at(u + du, v + dv), s.point_at(u, v + dv)], face_normal));
                }
            }
        }
    }
    if tris.is_empty() { None } else { Some(tris) }
}

// ══════════════════════════════════════════════════════════════════════════
// Generic trimmed face area
// ══════════════════════════════════════════════════════════════════════════

fn try_generic_trimmed_face_area(brep: &BRep, face: &Face, face_flat_idx: usize) -> Option<f64> {
    let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
    })?;
    let surf = &surf_idx;
    let [u0, u1, v0, v1] = surf.default_domain();
    let mut area = 0.0_f64;
    let n_u = 40; let n_v = 40;
    let du = (u1 - u0) / n_u as f64; let dv = (v1 - v0) / n_v as f64;
    for i in 0..n_u {
        for j in 0..n_v {
            let u = u0 + (i as f64 + 0.5) * du;
            let v = v0 + (j as f64 + 0.5) * dv;
            area += surface_normal_jacobian(surf, u, v) * du * dv;
        }
    }
    Some(area)
}
