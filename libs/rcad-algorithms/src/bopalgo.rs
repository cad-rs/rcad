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

/// ✅ OCCT-aligned: BOPAlgo_Tools::FillMap (hxx:83-102).
///   Adds a bidirectional connection between n1 and n2 in an adjacency map.
///   If either key doesn't exist, it's created.  rcad: uses HashMap<usize, Vec<usize>>.
pub fn fill_map(the_mili: &mut HashMap<usize, Vec<usize>>, n1: usize, n2: usize) {
    the_mili.entry(n1).or_default().push(n2);
    the_mili.entry(n2).or_default().push(n1);
}

/// ✅ OCCT-aligned: BOPAlgo_Tools::MakeBlocks (hxx:45-80).
///   Builds connected components from an adjacency map.
///   Each component is a BFS closure under the adjacency relation.
///   Equivalent to OCCT's template function NCollection_Map + chain building.
///
///   Args:
///     the_mili: adjacency map — for each key, list of adjacent keys.
///
///   Returns: Vec of connected components (each component is a Vec of keys).
pub fn make_blocks(the_mili: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
    let mut fence: HashSet<usize> = HashSet::new();
    let mut blocks: Vec<Vec<usize>> = Vec::new();

    // OCCT L51: for (i = 1; i <= aNb; ++i) — iterate all keys
    let keys: Vec<&usize> = the_mili.keys().collect();
    for &&key in &keys {
        if !fence.insert(key) {
            continue; // OCCT L55: if (!aMFence.Add(n)) continue
        }

        // OCCT L59-61: Start the chain — aChain.Append(n)
        let mut chain: Vec<usize> = Vec::new();
        chain.push(key);

        // OCCT L62-78: BFS-like traversal through adjacency
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(key);
        while let Some(n1) = queue.pop_front() {
            if let Some(adjacent) = the_mili.get(&n1) {
                for &n2 in adjacent {
                    if fence.insert(n2) {
                        chain.push(n2);
                        queue.push_back(n2);
                    }
                }
            }
        }

        blocks.push(chain);
    }

    blocks
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
