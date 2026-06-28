use glam::DVec3;
use std::collections::{HashMap, HashSet, VecDeque};
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

/// ✅ OCCT-aligned: BOPAlgo_Tools::TrsfToPoint (cxx:1912-1937).
///
/// Computes a translation from the combined bounding box of two boxes to a point.
/// Returns `Some(translation_vector)` when the point is sufficiently far from the
/// combined box center and the box size is small enough relative to the distance.
/// Returns `None` when the criteria rejects the transformation.
///
/// OCCT parameters:
///   theBox1, theBox2 — bounding boxes to unify.
///   theTrsf        — (output) the transform to fill.
///   thePoint       — target point.
///   theCriteria    — minimal distance criterion.
///
/// rcad: returns Option<DVec3> (the translation vector) instead of bool + gp_Trsf.
pub fn trsf_to_point(
    box1: &crate::bvh::Aabb,
    box2: &crate::bvh::Aabb,
    point: glam::DVec3,
    criteria: f64,
) -> Option<glam::DVec3> {
    // OCCT L1918-1920: Unify two boxes
    let mut a_box = *box1;
    a_box.expand_aabb(box2);

    // OCCT L1922-1923: Compute center of unified box and distance from point
    let a_b_center = (a_box.min + a_box.max) * 0.5;
    let a_pb_dist = (point - a_b_center).length();

    // OCCT L1924-1927: Reject if point is too close to box center
    if a_pb_dist < criteria {
        return None;
    }

    // OCCT L1929-1933: Compute box diagonal length; reject if box is too large
    //   relative to the distance (ratio > 1/criteria)
    let a_b_size = (a_box.max - a_box.min).length();
    if (a_b_size / a_pb_dist) > (1.0 / criteria) {
        return None;
    }

    // OCCT L1935: Set translation from box corner min to the point
    Some(point - a_box.min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fill_map_and_make_blocks_single_edge() {
        let mut m = HashMap::new();
        fill_map(&mut m, 1, 2);
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 1);
        assert!(blocks[0].contains(&1) && blocks[0].contains(&2));
    }

    #[test]
    fn test_make_blocks_two_components() {
        let mut m = HashMap::new();
        fill_map(&mut m, 1, 2);
        fill_map(&mut m, 3, 4);
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 2);
    }

    #[test]
    fn test_make_blocks_chain() {
        let mut m = HashMap::new();
        fill_map(&mut m, 1, 2);
        fill_map(&mut m, 2, 3);
        fill_map(&mut m, 3, 4);
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].len(), 4);
    }

    #[test]
    fn test_make_blocks_isolated() {
        let mut m = HashMap::new();
        m.entry(42).or_default();
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0], vec![42]);
    }

    #[test]
    fn test_make_blocks_empty() {
        let m: HashMap<usize, Vec<usize>> = HashMap::new();
        assert!(make_blocks(&m).is_empty());
    }
}
