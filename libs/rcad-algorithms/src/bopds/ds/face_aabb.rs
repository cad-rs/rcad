//! DS face AABB computation — used to build BVH from DS source shape data.
//!
//! `BOPDS_ShapeInfo::Box()` builds the AABB for each source
//! shape from its sub-shape bounding boxes. Here we compute per-face AABBs
//! from the DS's own vertex and surface data so the BVH uses the same index
//! space as the DS (flat face index).

use glam::DVec3;
use crate::bvh::Aabb;
use crate::bopds::ds::{DS, ShapeOrigin};
use rcad_kernel::geom::{Surface3, SurfaceEval};
use crate::tolerance::TOLERANCE_LINEAR_ULTRA_STRICT;

/// Build a `DsBvh` from all faces of the given origin.
///
/// Returns `None` when the face set is empty (OCCT: no BVH for empty operand).
/// `BOPDS_Iterator::Initialize` builds BVH over source shape
/// bounding boxes using `BOPDS_ShapeInfo::Box()`. Here we use the same index
/// space — the DS face index — so `candidate_pairs` returns DS-compatible indices.
pub fn build_face_bvh(ds: &DS, origin: ShapeOrigin) -> Option<crate::bvh::DsBvh> {
    let faces: Vec<usize> = ds.faces.iter().enumerate()
        .filter(|(_, f)| f.origin == origin)
        .map(|(i, _)| i)
        .collect();
    if faces.is_empty() {
        return None;
    }
    let mut indices = Vec::with_capacity(faces.len());
    let mut aabbs = Vec::with_capacity(faces.len());
    for &fi in &faces {
        indices.push(fi);
        aabbs.push(face_aabb(ds, fi));
    }
    Some(crate::bvh::DsBvh::build(indices, aabbs))
}

/// Compute the AABB of a DS face from its boundary vertices and surface type.
///
/// `BRepBndLib::Add(aF, aBox)` — expands the box by adding
/// face boundary vertices and analytical surface bounds for curved types.
pub fn face_aabb(ds: &DS, fi: usize) -> Aabb {
    let face = &ds.faces[fi];
    let mut aabb = Aabb::empty();

    // Seed AABB from boundary vertices.
    for &vi in &face.boundary_verts {
        if let Some(v) = ds.vertices.get(vi) {
            aabb.expand_point(v.point);
        }
    }

    // Expand for curved surface types (OCCT: BndLib_AddSurface).
    let surf = &face.surface;
    match surf {
        Surface3::Sphere(s) => {
            let r = s.radius.abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
            aabb.expand_point(s.center - DVec3::splat(r));
            aabb.expand_point(s.center + DVec3::splat(r));
        }
        Surface3::Cylinder(c) => {
            let ax = c.axis.normalize_or_zero();
            let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
            let u_dir = ax.cross(perp).normalize_or_zero();
            let v_dir = ax.cross(u_dir).normalize_or_zero();
            let r = c.radius.abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
            let [_, _, v0, v1] = surf.default_domain();
            for &vh in &[v0, v1] {
                for k in 0..8 {
                    let a = std::f64::consts::TAU * k as f64 / 8.0;
                    let p = c.origin + ax * vh + u_dir * r * a.cos() + v_dir * r * a.sin();
                    aabb.expand_point(p);
                }
            }
        }
        Surface3::Cone(c) => {
            let ax = c.axis.normalize_or_zero();
            let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
            let u_dir = ax.cross(perp).normalize_or_zero();
            let v_dir = ax.cross(u_dir).normalize_or_zero();
            let [_, _, v0, v1] = surf.default_domain();
            for &vh in &[v0, v1] {
                let r_at = (c.radius + vh * c.half_angle_rad.tan()).abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
                let center = c.apex + ax * vh;
                for k in 0..8 {
                    let a = std::f64::consts::TAU * k as f64 / 8.0;
                    let p = center + u_dir * r_at * a.cos() + v_dir * r_at * a.sin();
                    aabb.expand_point(p);
                }
            }
        }
        Surface3::Torus(t) => {
            let r_out = t.major_radius.abs() + t.minor_radius.abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
            let ax = t.axis.normalize_or_zero();
            let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
            let u_dir = ax.cross(perp).normalize_or_zero();
            let v_dir = ax.cross(u_dir).normalize_or_zero();
            for k in 0..8 {
                let a = std::f64::consts::TAU * k as f64 / 8.0;
                let c = t.center + u_dir * t.major_radius * a.cos() + v_dir * t.major_radius * a.sin();
                aabb.expand_point(c + ax * t.minor_radius);
                aabb.expand_point(c - ax * t.minor_radius);
            }
        }
        _ => {
            // Plane, BSpline, Bezier: use surface sampling (3x3 grid).
            let domain = surf.default_domain();
            let [u0, u1, v0, v1] = domain;
            for i in 0..=2 {
                for j in 0..=2 {
                    let u = u0 + (u1 - u0) * i as f64 / 2.0;
                    let v = v0 + (v1 - v0) * j as f64 / 2.0;
                    let p = surf.point_at(u, v);
                    if p.is_finite() { aabb.expand_point(p); }
                }
            }
        }
    }

    // Ensure non-zero extent.
    let size = aabb.max - aabb.min;
    if size.x < TOLERANCE_LINEAR_ULTRA_STRICT { aabb.min.x -= TOLERANCE_LINEAR_ULTRA_STRICT; aabb.max.x += TOLERANCE_LINEAR_ULTRA_STRICT; }
    if size.y < TOLERANCE_LINEAR_ULTRA_STRICT { aabb.min.y -= TOLERANCE_LINEAR_ULTRA_STRICT; aabb.max.y += TOLERANCE_LINEAR_ULTRA_STRICT; }
    if size.z < TOLERANCE_LINEAR_ULTRA_STRICT { aabb.min.z -= TOLERANCE_LINEAR_ULTRA_STRICT; aabb.max.z += TOLERANCE_LINEAR_ULTRA_STRICT; }
    aabb
}
