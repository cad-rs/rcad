//! BoxTree — linear-scan bounding box collection using kernel BndBox.

/// Simple linear-scan BVH using kernel BndBox.
pub struct BoxTree {
    items: Vec<(usize, rcad_kernel::math::bnd::BndBox)>,
}

impl BoxTree {
    pub fn build(indices: Vec<usize>, aabbs: Vec<rcad_kernel::math::bnd::BndBox>) -> Self {
        let items = indices.into_iter().zip(aabbs).collect();
        BoxTree { items }
    }

    pub fn find_overlapping(&self, query: &rcad_kernel::math::bnd::BndBox) -> Vec<usize> {
        self.items.iter().filter_map(|(i, b)| {
            if !b.is_out_point(query.corner_min().unwrap_or_default())
                && !query.is_out_point(b.corner_min().unwrap_or_default())
            {
                Some(*i)
            } else {
                None
            }
        }).collect()
    }

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
