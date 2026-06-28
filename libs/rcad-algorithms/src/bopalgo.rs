use glam::DVec3;
use crate::bopds::ds::DS;
use crate::bvh::Aabb;
use crate::classify::{Classification, classify_point};

/// ✅ OCCT-aligned: BOPAlgo_Tools::ClassifyFaces (cxx:1622-1747).
///
/// Classifies result faces relatively draft solids.  For each solid,
/// collects the faces classified as IN into the returned Vec.
///
/// OCCT builds a BVH tree of face boxes (L1670-1680), then runs per-solid
/// classification jobs (L1736-1738, BOPAlgo_FillIn3DParts::Perform).
/// rcad: for each (face, solid) pair with overlapping AABB, calls
/// classify_point against the solid's DS face set.
///
/// Args:
///   the_faces: result face indices, each with sample point at
///     face_samples[i].
///   face_samples: 3D sample points for each face (one per the_faces entry).
///   the_solids: each solid = Vec of shell groups of DS face indices.
///   ds: data structure.
///   aabb_of_face: bounding box for each face (parallel to the_faces).
///   aabb_of_solid: bounding box for each solid.
///
/// Returns: for each solid index, list of result FACE INDICES (values from
///   the_faces, not positions) classified as IN that solid.
pub fn classify_faces(
    the_faces: &[usize],
    face_samples: &[DVec3],
    the_solids: &[Vec<Vec<usize>>],
    ds: &DS,
    aabb_of_face: &[Aabb],
    aabb_of_solid: &[Aabb],
) -> Vec<Vec<usize>> {
    let n_solids = the_solids.len();
    let mut the_in_parts: Vec<Vec<usize>> = vec![Vec::new(); n_solids];

    // Precompute flat DS face sets per solid (for classify_point)
    let solid_faces: Vec<Vec<usize>> = the_solids.iter()
        .map(|shells| shells.iter().flat_map(|sh| sh.iter().copied()).collect())
        .collect();

    for (si, sfaces) in solid_faces.iter().enumerate() {
        if sfaces.is_empty() { continue; }
        let sbox = &aabb_of_solid[si];
        for (pi, &fi) in the_faces.iter().enumerate() {
            if pi >= face_samples.len() { continue; }
            if !aabb_of_face[pi].intersects(sbox) { continue; }
            let class = classify_point(face_samples[pi], sfaces, ds);
            if class == Classification::In {
                the_in_parts[si].push(fi);
            }
        }
    }
    the_in_parts
}
