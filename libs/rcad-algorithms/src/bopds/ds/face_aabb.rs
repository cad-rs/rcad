//! DS face AABB computation — used to build BVH from DS source shape data.
//!
//! `BOPDS_ShapeInfo::Box()` builds the AABB for each source
//! shape from its sub-shape bounding boxes. Here we compute per-face AABBs
//! from the DS's own vertex and surface data so the BVH uses the same index
//! space as the DS (flat face index).

use glam::DVec3;
use crate::boptools::bvh::Aabb;
use crate::bopds::ds::{DS, ShapeOrigin};
use crate::tolerance::TOLERANCE_LINEAR_ULTRA_STRICT;

/// Build a `BoxTree` from all faces of the given origin.
///
/// Returns `None` when the face set is empty (OCCT: no BVH for empty operand).
/// `BOPDS_Iterator::Initialize` builds BVH over source shape
/// bounding boxes using `BOPDS_ShapeInfo::Box()`. Here we use the same index
/// space — the DS face index — so `candidate_pairs` returns DS-compatible indices.
pub fn build_face_bvh(ds: &DS, origin: ShapeOrigin) -> Option<crate::boptools::bvh::BoxTree> {
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
    Some(crate::boptools::bvh::BoxTree::build(indices, aabbs))
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

    // ✅ OCCT-aligned: delegate to bnd_lib::surface_bounds (BndLib_AddSurface per-type dispatch).
    let surf_bbox = crate::bnd_lib::surface_bounds(&face.surface, TOLERANCE_LINEAR_ULTRA_STRICT);
    if surf_bbox.is_valid() {
        aabb.min = aabb.min.min(surf_bbox.min);
        aabb.max = aabb.max.max(surf_bbox.max);
    }

    // Ensure non-zero extent.
    let size = aabb.max - aabb.min;
    if size.x < TOLERANCE_LINEAR_ULTRA_STRICT { aabb.min.x -= TOLERANCE_LINEAR_ULTRA_STRICT; aabb.max.x += TOLERANCE_LINEAR_ULTRA_STRICT; }
    if size.y < TOLERANCE_LINEAR_ULTRA_STRICT { aabb.min.y -= TOLERANCE_LINEAR_ULTRA_STRICT; aabb.max.y += TOLERANCE_LINEAR_ULTRA_STRICT; }
    if size.z < TOLERANCE_LINEAR_ULTRA_STRICT { aabb.min.z -= TOLERANCE_LINEAR_ULTRA_STRICT; aabb.max.z += TOLERANCE_LINEAR_ULTRA_STRICT; }
    aabb
}
