//! OCCT BRepGProp::VolumeProperties: volume and centroid computation.
//!
//! Split from the original properties.rs.

use glam::DVec3;

use crate::BRep;
use crate::geom::{SphericalSurface, Surface3, SurfaceEval};
use crate::topo::topods;
use crate::topo::topology::Face;
use crate::base::gprop::tri::{
    self, face_flat_iter, face_triangles_pub, sample_wire_polyline_3d,
    trim_almost_closed_polyline, tet_signed_volume,
    point_in_spherical_polygon_3d, spherical_holed_uv_mask_setup, SphereHoledMaskCtx,
};

/// Signed volume of a BRep solid.
pub fn signed_volume(brep: &topods::BRep) -> f64 {
    let mut vol = 0.0;
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        let tris = face_triangles_pub(brep, *fi);
        for [a, b, c] in &tris {
            vol += tet_signed_volume(*a, *b, *c);
        }
    }
    vol
}

/// Absolute volume of a BRep solid.
pub fn volume(brep: &topods::BRep) -> f64 {
    signed_volume(brep).abs()
}

/// Centroid of a BRep solid.
pub fn centroid(brep: &topods::BRep) -> DVec3 {
    let mut total_vol = 0.0;
    let mut center = DVec3::ZERO;
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        let tris = face_triangles_pub(brep, *fi);
        for [a, b, c] in &tris {
            let tv = tet_signed_volume(*a, *b, *c);
            total_vol += tv;
            center += (*a + *b + *c) * tv * 0.25;
        }
    }
    if total_vol.abs() > 1e-15 { center / total_vol } else { DVec3::ZERO }
}

/// Sphere-specific volume using masked parameter integration.
fn sphere_holed_mask_param_volume_sum(s: &SphericalSurface, ctx: &SphereHoledMaskCtx) -> f64 {
    const N: usize = 30;
    const GL_NEG: f64 = -0.5773502691896257;
    const GL_POS: f64 = 0.5773502691896257;
    let gl_pts = [GL_NEG, GL_POS];
    let umin = ctx.umin; let umax = ctx.umax; let vmin = ctx.vmin; let vmax = ctx.vmax;
    let du = (umax - umin) / N as f64; let dv = (vmax - vmin) / N as f64;
    let r3 = s.radius * s.radius * s.radius;
    let cell_area = du * dv / 4.0;
    let mut total = 0.0;
    for i in 0..N {
        let u_mid = umin + (i as f64 + 0.5) * du;
        for j in 0..N {
            let v_mid = vmin + (j as f64 + 0.5) * dv;
            let mut sum = 0.0;
            for &gu in &gl_pts {
                let u = u_mid + gu * du * 0.5;
                for &gv in &gl_pts {
                    let v = v_mid + gv * dv * 0.5;
                    if v < 0.0 || v > std::f64::consts::PI { continue; }
                    let p3d = s.point_at(u, v);
                    if !point_in_spherical_polygon_3d(&ctx.outer_3d, p3d) { continue; }
                    if ctx.inner_3d.iter().any(|h3d| point_in_spherical_polygon_3d(h3d, p3d)) { continue; }
                    sum += v.sin() * v.cos(); // dV = r³·sin(v)·cos(v) dz component
                }
            }
            if sum > 0.0 { total += r3 * cell_area * sum; }
        }
    }
    total
}

fn try_sphere_face_analytic_volume(brep: &BRep, face: &Face, face_flat_idx: usize) -> Option<f64> {
    let surf_idx = brep.tshapes.get(face_flat_idx).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
    })?;
    match &surf_idx {
        Surface3::Sphere(s) => {
            let ctx = spherical_holed_uv_mask_setup(s, brep, face)?;
            let v = sphere_holed_mask_param_volume_sum(s, &ctx);
            if v > 0.0 { Some(v) } else { None }
        }
        _ => None,
    }
}
