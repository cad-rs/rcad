//! BoxTree — simple linear-scan bounding box collection.
use crate::bnd_box::Aabb;

/// Simple linear-scan BVH (no tree structure — O(n) search).
pub struct BoxTree {
    items: Vec<(usize, Aabb)>,
}

impl BoxTree {
    pub fn build(indices: Vec<usize>, aabbs: Vec<Aabb>) -> Self {
        let items = indices.into_iter().zip(aabbs).collect();
        BoxTree { items }
    }

    pub fn find_overlapping(&self, _query: &Aabb) -> Vec<usize> {
        self.items.iter().map(|(i, _)| *i).collect()
    }

    /// Return candidate pairs (all pairs for linear scan).
    pub fn candidate_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for i in 0..self.items.len() {
            for j in (i + 1)..self.items.len() {
                pairs.push((self.items[i].0, self.items[j].0));
            }
        }
        pairs
    }
}
