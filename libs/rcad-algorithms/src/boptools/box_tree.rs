//! BVH for DS entity pair culling — OCCT BOPTools_BoxTree / BOPTools_BoxSelector.
//!
//! Uses LBVH (Linear BVH) builder via Morton codes + radix sort,
//! matching OCCT BVH_LinearBuilder + BVH_RadixSorter.
//!
//! OCCT correspondence:
//! - BOPTools_BoxTree           = BoxTree (this file)
//! - BVH_LinearBuilder          = build() method
//! - BVH_RadixSorter            = morton sort inside build()
//! - BOPTools_BoxPairSelector   = self_pairs() method
//! - BOPTools_BoxSelector       = query_aabb()

use crate::bnd_box::Aabb;
use crate::boptools::bvh_tree::BvhNode;
use glam::DVec3;

// ── Morton code LUT (OCCT BVH_RadixSorter.hxx L34-63) ──
const MORTON_LUT: [u32; 256] = [
    0x000000, 0x000001, 0x000008, 0x000009, 0x000040, 0x000041, 0x000048, 0x000049, 0x000200,
    0x000201, 0x000208, 0x000209, 0x000240, 0x000241, 0x000248, 0x000249, 0x001000, 0x001001,
    0x001008, 0x001009, 0x001040, 0x001041, 0x001048, 0x001049, 0x001200, 0x001201, 0x001208,
    0x001209, 0x001240, 0x001241, 0x001248, 0x001249, 0x008000, 0x008001, 0x008008, 0x008009,
    0x008040, 0x008041, 0x008048, 0x008049, 0x008200, 0x008201, 0x008208, 0x008209, 0x008240,
    0x008241, 0x008248, 0x008249, 0x009000, 0x009001, 0x009008, 0x009009, 0x009040, 0x009041,
    0x009048, 0x009049, 0x009200, 0x009201, 0x009208, 0x009209, 0x009240, 0x009241, 0x009248,
    0x009249, 0x040000, 0x040001, 0x040008, 0x040009, 0x040040, 0x040041, 0x040048, 0x040049,
    0x040200, 0x040201, 0x040208, 0x040209, 0x040240, 0x040241, 0x040248, 0x040249, 0x041000,
    0x041001, 0x041008, 0x041009, 0x041040, 0x041041, 0x041048, 0x041049, 0x041200, 0x041201,
    0x041208, 0x041209, 0x041240, 0x041241, 0x041248, 0x041249, 0x048000, 0x048001, 0x048008,
    0x048009, 0x048040, 0x048041, 0x048048, 0x048049, 0x048200, 0x048201, 0x048208, 0x048209,
    0x048240, 0x048241, 0x048248, 0x048249, 0x049000, 0x049001, 0x049008, 0x049009, 0x049040,
    0x049041, 0x049048, 0x049049, 0x049200, 0x049201, 0x049208, 0x049209, 0x049240, 0x049241,
    0x049248, 0x049249, 0x200000, 0x200001, 0x200008, 0x200009, 0x200040, 0x200041, 0x200048,
    0x200049, 0x200200, 0x200201, 0x200208, 0x200209, 0x200240, 0x200241, 0x200248, 0x200249,
    0x201000, 0x201001, 0x201008, 0x201009, 0x201040, 0x201041, 0x201048, 0x201049, 0x201200,
    0x201201, 0x201208, 0x201209, 0x201240, 0x201241, 0x201248, 0x201249, 0x208000, 0x208001,
    0x208008, 0x208009, 0x208040, 0x208041, 0x208048, 0x208049, 0x208200, 0x208201, 0x208208,
    0x208209, 0x208240, 0x208241, 0x208248, 0x208249, 0x209000, 0x209001, 0x209008, 0x209009,
    0x209040, 0x209041, 0x209048, 0x209049, 0x209200, 0x209201, 0x209208, 0x209209, 0x209240,
    0x209241, 0x209248, 0x209249, 0x240000, 0x240001, 0x240008, 0x240009, 0x240040, 0x240041,
    0x240048, 0x240049, 0x240200, 0x240201, 0x240208, 0x240209, 0x240240, 0x240241, 0x240248,
    0x240249, 0x241000, 0x241001, 0x241008, 0x241009, 0x241040, 0x241041, 0x241048, 0x241049,
    0x241200, 0x241201, 0x241208, 0x241209, 0x241240, 0x241241, 0x241248, 0x241249, 0x248000,
    0x248001, 0x248008, 0x248009, 0x248040, 0x248041, 0x248048, 0x248049, 0x248200, 0x248201,
    0x248208, 0x248209, 0x248240, 0x248241, 0x248248, 0x248249, 0x249000, 0x249001, 0x249008,
    0x249009, 0x249040, 0x249041, 0x249048, 0x249049, 0x249200, 0x249201, 0x249208, 0x249209,
    0x249240, 0x249241, 0x249248, 0x249249,
];

/// encode 10-bit xyz → 30-bit morton code (OCCT BVH_RadixSorter.hxx L70-81)
fn encode_morton(vx: u32, vy: u32, vz: u32) -> u32 {
    (MORTON_LUT[(vx & 0xFF) as usize] | (MORTON_LUT[((vx >> 8) & 0x03) as usize] << 24))
        | ((MORTON_LUT[(vy & 0xFF) as usize] | (MORTON_LUT[((vy >> 8) & 0x03) as usize] << 24))
            << 1)
        | ((MORTON_LUT[(vz & 0xFF) as usize] | (MORTON_LUT[((vz >> 8) & 0x03) as usize] << 24))
            << 2)
}

/// MSD radix sort on 30-bit morton codes (OCCT BVH::RadixSorter L218-227)
fn radix_sort_msd(links: &mut [(u32, usize)], start: usize, end: usize, digit: i32) {
    if end - start <= 4 || digit < 0 {
        return;
    }
    let bit = 1u32 << (digit as u32);
    let mut split = start;
    for i in start..end {
        if links[i].0 & bit == 0 {
            links.swap(i, split);
            split += 1;
        }
    }
    if split > start && split < end {
        radix_sort_msd(links, start, split, digit - 1);
        radix_sort_msd(links, split, end, digit - 1);
    } else {
        radix_sort_msd(links, start, end, digit - 1);
    }
}

fn emit_hierarchy(
    links: &[(u32, usize)],
    start: usize,
    end: usize,
    nodes: &mut Vec<BvhNode>,
    next_slot: &mut usize,
) -> usize {
    if end - start <= 4 {
        let idx = *next_slot;
        *next_slot += 1;
        nodes.push(BvhNode::Leaf {
            aabb: Aabb::empty(),
            start,
            end,
        });
        return idx;
    }
    let first = links[start].0;
    let last = links[end - 1].0;
    let diff = first ^ last;
    let hb = (0..=29).rev().find(|&b| diff & (1u32 << b) != 0);
    let split = match hb {
        None => (start + end) / 2,
        Some(bit) => {
            let mask = 1u32 << bit;
            let mut lo = start;
            let mut hi = end;
            while lo < hi {
                let mid = (lo + hi) / 2;
                if links[mid].0 & mask == 0 {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            lo
        }
    };
    if split == start || split == end {
        return emit_hierarchy(links, start, end, nodes, next_slot);
    }
    let idx = *next_slot;
    *next_slot += 1;
    // Reserve a slot: push dummy Internal, fill children later
    nodes.push(BvhNode::Internal {
        aabb: Aabb::empty(),
        left: 0,
        right: 0,
    });
    let left_child = emit_hierarchy(links, start, split, nodes, next_slot);
    let right_child = emit_hierarchy(links, split, end, nodes, next_slot);
    match &mut nodes[idx] {
        BvhNode::Internal { left, right, .. } => {
            *left = left_child;
            *right = right_child;
        }
        _ => unreachable!(),
    }
    idx
}

fn update_bounds(node: usize, nodes: &mut [BvhNode], aabbs: &[Aabb]) {
    // Extract child indices BEFORE recursing (avoids multi-borrow conflict).
    let (child_left, child_right) = match &nodes[node] {
        BvhNode::Internal { left, right, .. } => (*left, *right),
        _ => (usize::MAX, usize::MAX), // leaf marker
    };
    if child_left == usize::MAX {
        // Leaf: compute AABB from sorted aabbs range
        let (start, end) = match &nodes[node] {
            BvhNode::Leaf { start, end, .. } => (*start, *end),
            _ => unreachable!(),
        };
        let mut mn = aabbs[start].min;
        let mut mx = aabbs[start].max;
        let mut gap = aabbs[start].gap;
        for i in (start + 1)..end {
            mn = mn.min(aabbs[i].min);
            mx = mx.max(aabbs[i].max);
            if aabbs[i].gap > gap {
                gap = aabbs[i].gap;
            }
        }
        match &mut nodes[node] {
            BvhNode::Leaf { aabb, .. } => {
                *aabb = Aabb {
                    min: mn,
                    max: mx,
                    gap,
                }
            }
            _ => unreachable!(),
        }
        return;
    }
    update_bounds(child_left, nodes, aabbs);
    update_bounds(child_right, nodes, aabbs);
    let lbb = &nodes[child_left];
    let rbb = &nodes[child_right];
    let mn = lbb.aabb().min.min(rbb.aabb().min);
    let mx = lbb.aabb().max.max(rbb.aabb().max);
    let gap = lbb.aabb().gap.max(rbb.aabb().gap);
    match &mut nodes[node] {
        BvhNode::Internal { aabb, .. } => {
            *aabb = Aabb {
                min: mn,
                max: mx,
                gap,
            }
        }
        _ => unreachable!(),
    }
}

/// BVH for DS entity pair culling — OCCT BOPTools_BoxTree.
pub struct BoxTree {
    nodes: Vec<BvhNode>,
    indices: Vec<usize>,
    aabbs: Vec<Aabb>,
}

impl BoxTree {
    /// LBVH builder — OCCT BVH_LinearBuilder::Build + BVH_RadixSorter::Perform
    pub fn build(indices: Vec<usize>, aabbs: Vec<Aabb>) -> Self {
        let n = indices.len();
        if n == 0 {
            return Self {
                nodes: Vec::new(),
                indices,
                aabbs,
            };
        }

        // Scene AABB
        let mut scene_min = aabbs[0].min;
        let mut scene_max = aabbs[0].max;
        for a in &aabbs {
            scene_min = scene_min.min(a.min);
            scene_max = scene_max.max(a.max);
        }
        const MIN_EXT: f64 = 1e-12;
        let extent = (scene_max - scene_min).max(DVec3::splat(MIN_EXT));
        let inv_extent = DVec3::new(1024.0 / extent.x, 1024.0 / extent.y, 1024.0 / extent.z);

        // Morton codes
        let mut links: Vec<(u32, usize)> = Vec::with_capacity(n);
        for i in 0..n {
            let c = aabbs[i].center();
            let vf = (c - scene_min) * inv_extent;
            let vx = (vf.x as i32).clamp(0, 1023) as u32;
            let vy = (vf.y as i32).clamp(0, 1023) as u32;
            let vz = (vf.z as i32).clamp(0, 1023) as u32;
            links.push((encode_morton(vx, vy, vz), i));
        }

        // Radix sort
        radix_sort_msd(&mut links, 0, n, 29);

        // Reorder
        let mut sorted_idx = Vec::with_capacity(n);
        let mut sorted_abb = Vec::with_capacity(n);
        for &(_, orig_i) in &links {
            sorted_idx.push(indices[orig_i]);
            sorted_abb.push(aabbs[orig_i]);
        }

        // Emit hierarchy
        let mut nodes: Vec<BvhNode> = Vec::new();
        let mut next_slot = 0usize;
        emit_hierarchy(&links, 0, n, &mut nodes, &mut next_slot);

        // Compute AABBs
        for i in 0..nodes.len() {
            update_bounds(i, &mut nodes, &sorted_abb);
        }

        Self {
            nodes,
            indices: sorted_idx,
            aabbs: sorted_abb,
        }
    }

    /// Self-pair query — OCCT BOPTools_PairSelector SetSame(true)
    pub fn self_pairs(&self) -> Vec<(usize, usize)> {
        PairSelector::select(self)
    }

    /// Dual-tree pair query — OCCT BOPTools_PairSelector
    pub fn candidate_pairs(bvh_a: &BoxTree, bvh_b: &BoxTree) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if bvh_a.nodes.is_empty() || bvh_b.nodes.is_empty() {
            return pairs;
        }
        Self::candidate_pairs_node(bvh_a, 0, bvh_b, 0, &mut pairs);
        pairs
    }

    fn candidate_pairs_node(
        bvh_a: &BoxTree,
        na: usize,
        bvh_b: &BoxTree,
        nb: usize,
        pairs: &mut Vec<(usize, usize)>,
    ) {
        if !bvh_a.nodes[na].aabb().intersects(bvh_b.nodes[nb].aabb()) {
            return;
        }
        match (&bvh_a.nodes[na], &bvh_b.nodes[nb]) {
            (
                BvhNode::Leaf {
                    start: sa, end: ea, ..
                },
                BvhNode::Leaf {
                    start: sb, end: eb, ..
                },
            ) => {
                for i in *sa..*ea {
                    for j in *sb..*eb {
                        if bvh_a.aabbs[i].intersects(&bvh_b.aabbs[j]) {
                            pairs.push((bvh_a.indices[i], bvh_b.indices[j]));
                        }
                    }
                }
            }
            (
                BvhNode::Internal {
                    left: la,
                    right: ra,
                    ..
                },
                _,
            ) => {
                Self::candidate_pairs_node(bvh_a, *la, bvh_b, nb, pairs);
                Self::candidate_pairs_node(bvh_a, *ra, bvh_b, nb, pairs);
            }
            (
                _,
                BvhNode::Internal {
                    left: lb,
                    right: rb,
                    ..
                },
            ) => {
                Self::candidate_pairs_node(bvh_a, na, bvh_b, *lb, pairs);
                Self::candidate_pairs_node(bvh_a, na, bvh_b, *rb, pairs);
            }
        }
    }

    /// Query all items whose AABB overlaps the query AABB.
    pub fn query_aabb(&self, query: &Aabb) -> Vec<usize> {
        let mut res = Vec::new();
        if self.nodes.is_empty() {
            return res;
        }
        self.query_aabb_node(0, query, &mut res);
        res
    }

    fn query_aabb_node(&self, node: usize, query: &Aabb, res: &mut Vec<usize>) {
        if !self.nodes[node].aabb().intersects(query) {
            return;
        }
        match &self.nodes[node] {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    if self.aabbs[i].intersects(query) {
                        res.push(self.indices[i]);
                    }
                }
            }
            BvhNode::Internal { left, right, .. } => {
                self.query_aabb_node(*left, query, res);
                self.query_aabb_node(*right, query, res);
            }
        }
    }
}

// ── BOPTools_PairSelector (OCCT L25-108) ──
// OCCT BOPTools_PairSelector<Dimension> with SetSame(true).
// RejectNode  = node-level AABB overlap check (BVH_Box::IsOut)
// RejectElement = element-level AABB check + same-tree dedup (ID1 >= ID2)
// Accept      = store pair as (Element(ID1), Element(ID2))
pub struct PairSelector;

impl PairSelector {
    pub fn select(tree: &BoxTree) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if tree.nodes.is_empty() {
            return pairs;
        }
        Self::self_query(tree, 0, &mut pairs);
        pairs
    }

    // Self-query traversal (OCCT BVH_PairTraverse with same BVH sets)
    fn self_query(tree: &BoxTree, node: usize, pairs: &mut Vec<(usize, usize)>) {
        match &tree.nodes[node] {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    for j in (i + 1)..*end {
                        if !Self::reject_element(tree, tree, i, j, true) {
                            Self::accept(tree, tree, i, j, pairs);
                        }
                    }
                }
            }
            BvhNode::Internal { left, right, .. } => {
                Self::self_query(tree, *left, pairs);
                Self::self_query(tree, *right, pairs);
                Self::cross(tree, *left, *right, pairs);
            }
        }
    }

    // Cross-subtree traversal (OCCT BVH_PairTraverse with two different nodes)
    fn cross(tree: &BoxTree, na: usize, nb: usize, pairs: &mut Vec<(usize, usize)>) {
        // RejectNode: node AABBs don't overlap
        if Self::reject_node(&tree.nodes[na], &tree.nodes[nb]) {
            return;
        }
        match (&tree.nodes[na], &tree.nodes[nb]) {
            (
                BvhNode::Leaf {
                    start: sa, end: ea, ..
                },
                BvhNode::Leaf {
                    start: sb, end: eb, ..
                },
            ) => {
                for i in *sa..*ea {
                    for j in *sb..*eb {
                        if !Self::reject_element(tree, tree, i, j, true) {
                            Self::accept(tree, tree, i, j, pairs);
                        }
                    }
                }
            }
            (
                BvhNode::Internal {
                    left: la,
                    right: ra,
                    ..
                },
                _,
            ) => {
                Self::cross(tree, *la, nb, pairs);
                Self::cross(tree, *ra, nb, pairs);
            }
            (
                _,
                BvhNode::Internal {
                    left: lb,
                    right: rb,
                    ..
                },
            ) => {
                Self::cross(tree, na, *lb, pairs);
                Self::cross(tree, na, *rb, pairs);
            }
        }
    }

    // OCCT RejectNode: BVH_Box::IsOut on node AABBs
    fn reject_node(na: &BvhNode, nb: &BvhNode) -> bool {
        !na.aabb().intersects(nb.aabb())
    }

    // OCCT RejectElement: (mySameBVHs && theID1 >= theID2) || Box(theID1).IsOut(Box(theID2))
    fn reject_element(
        tree_a: &BoxTree,
        tree_b: &BoxTree,
        ia: usize,
        ib: usize,
        same: bool,
    ) -> bool {
        if same && ia >= ib {
            return true;
        }
        !tree_a.aabbs[ia].intersects(&tree_b.aabbs[ib])
    }

    // OCCT Accept: append(Element(theID1), Element(theID2))
    fn accept(
        tree_a: &BoxTree,
        tree_b: &BoxTree,
        ia: usize,
        ib: usize,
        pairs: &mut Vec<(usize, usize)>,
    ) {
        pairs.push((tree_a.indices[ia], tree_b.indices[ib]));
    }
}
