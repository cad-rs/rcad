//! BVH for DS entity pair culling — OCCT BOPTools_BoxTree / BOPTools_BoxSelector.
//!
//! OCCT correspondence:
//! - BOPTools_BoxTree           = BoxTree (this file)
//! - BVH_TreeBase               = BvhTreeBase (bvh_tree.rs)
//! - BVH_LinearBuilder          = build() method
//! - BVH_RadixSorter            = morton sort inside build()
//! - BOPTools_BoxPairSelector   = self_pairs() method
//! - BOPTools_BoxSelector       = query_aabb()

use crate::bop::tools::bvh_tree::BvhTreeBase;
use glam::DVec3;

/// Input bounding box used during BVH build (OCCT Bnd_Box equivalent).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Aabb {
    pub min: DVec3,
    pub max: DVec3,
    pub gap: f64,
}

impl Aabb {
    pub fn empty() -> Self { Self { min: DVec3::splat(f64::INFINITY), max: DVec3::splat(f64::NEG_INFINITY), gap: 0.0 } }
    pub fn center(&self) -> DVec3 { (self.min + self.max) * 0.5 }
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min.x - self.gap <= other.max.x + other.gap && self.max.x + self.gap >= other.min.x - other.gap
            && self.min.y - self.gap <= other.max.y + other.gap && self.max.y + self.gap >= other.min.y - other.gap
            && self.min.z - self.gap <= other.max.z + other.gap && self.max.z + self.gap >= other.min.z - other.gap
    }
}

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

fn encode_morton(vx: u32, vy: u32, vz: u32) -> u32 {
    (MORTON_LUT[(vx & 0xFF) as usize] | (MORTON_LUT[((vx >> 8) & 0x03) as usize] << 24))
        | ((MORTON_LUT[(vy & 0xFF) as usize] | (MORTON_LUT[((vy >> 8) & 0x03) as usize] << 24)) << 1)
        | ((MORTON_LUT[(vz & 0xFF) as usize] | (MORTON_LUT[((vz >> 8) & 0x03) as usize] << 24)) << 2)
}

fn radix_sort_msd(links: &mut [(u32, usize)], start: usize, end: usize, digit: i32) {
    if end - start <= 4 || digit < 0 { return; }
    let bit = 1u32 << (digit as u32);
    let mut split = start;
    for i in start..end {
        if links[i].0 & bit == 0 { links.swap(i, split); split += 1; }
    }
    if split > start && split < end {
        radix_sort_msd(links, start, split, digit - 1);
        radix_sort_msd(links, split, end, digit - 1);
    } else {
        radix_sort_msd(links, start, end, digit - 1);
    }
}

fn emit_hierarchy(
    links: &[(u32, usize)], start: usize, end: usize,
    nodes: &mut Vec<(i32, i32, i32, i32)>, // (isOuter, startPrim, endPrim, level)
    next_slot: &mut usize, depth: i32,
) -> usize {
    if end - start <= 4 {
        let idx = *next_slot; *next_slot += 1;
        nodes.push((0, start as i32, (end - 1) as i32, depth));
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
            let mut lo = start; let mut hi = end;
            while lo < hi { let mid = (lo + hi) / 2; if links[mid].0 & mask == 0 { lo = mid + 1; } else { hi = mid; } }
            lo
        }
    };
    if split == start || split == end {
        return emit_hierarchy(links, start, end, nodes, next_slot, depth);
    }
    let idx = *next_slot; *next_slot += 1;
    nodes.push((0, 0, 0, depth)); // inner node, filled after children
    let left_child = emit_hierarchy(links, start, split, nodes, next_slot, depth + 1);
    let right_child = emit_hierarchy(links, split, end, nodes, next_slot, depth + 1);
    // Store children in inner nodes: left = 2*idx+1, right = 2*idx+2
    // (no explicit child fields needed — heap layout)
    idx
}

/// OCCT BOPTools_BoxSet / BOPTools_BoxTree — BVH with AABB array + LinearBuilder.
pub struct BoxTree {
    pub(crate) tree: BvhTreeBase,
    pub(crate) indices: Vec<usize>,
    pub(crate) aabbs: Vec<Aabb>,
}

impl BoxTree {
    /// LBVH builder — OCCT BVH_LinearBuilder + BVH_RadixSorter.
    pub fn build(indices: Vec<usize>, aabbs: Vec<Aabb>) -> Self {
        let n = indices.len();
        if n == 0 { return Self { tree: BvhTreeBase::new(), indices, aabbs }; }

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

        radix_sort_msd(&mut links, 0, n, 29);

        // Reorder
        let mut sorted_idx = Vec::with_capacity(n);
        let mut sorted_abb = Vec::with_capacity(n);
        for &(_, orig_i) in &links {
            sorted_idx.push(indices[orig_i]);
            sorted_abb.push(aabbs[orig_i]);
        }

        // Emit hierarchy: vec of (isOuter, startPrim, endPrim, level)
        let mut raw_nodes: Vec<(i32, i32, i32, i32)> = Vec::new();
        let mut next_slot = 0usize;
        emit_hierarchy(&links, 0, n, &mut raw_nodes, &mut next_slot, 0);

        // Compute node AABBs (post-order) and push to BvhTreeBase
        let mut tree = BvhTreeBase::new();
        // First pass: create nodes with empty AABBs
        for i in 0..raw_nodes.len() {
            let (is_outer, sp, ep, level) = raw_nodes[i];
            tree.push_node(is_outer != 0, sp, ep, level, DVec3::ZERO, DVec3::ZERO);
        }
        // Second pass: compute AABBs bottom-up
        Self::compute_bounds(0, &raw_nodes, &sorted_abb, &mut tree);

        Self { tree, indices: sorted_idx, aabbs: sorted_abb }
    }

    fn compute_bounds(node: usize, raw: &[(i32, i32, i32, i32)], aabbs: &[Aabb], tree: &mut BvhTreeBase) {
        let (is_outer, sp, ep, _level) = raw[node];
        if is_outer != 0 {
            // Leaf: compute AABB from element boxes
            let mut mn = aabbs[sp as usize].min;
            let mut mx = aabbs[sp as usize].max;
            for i in (sp as usize + 1)..=ep as usize {
                mn = mn.min(aabbs[i].min);
                mx = mx.max(aabbs[i].max);
            }
            tree.min_points[node] = mn;
            tree.max_points[node] = mx;
            return;
        }
        let left = BvhTreeBase::left_child(node);
        let right = BvhTreeBase::right_child(node);
        if left < tree.length() { Self::compute_bounds(left, raw, aabbs, tree); }
        if right < tree.length() { Self::compute_bounds(right, raw, aabbs, tree); }
        let mn = if left < tree.length() { tree.min_points[left].min(tree.min_points.get(right).copied().unwrap_or(DVec3::splat(f64::INFINITY))) } else { DVec3::splat(f64::INFINITY) };
        let mx = if left < tree.length() { tree.max_points[left].max(tree.max_points.get(right).copied().unwrap_or(DVec3::splat(f64::NEG_INFINITY))) } else { DVec3::splat(f64::NEG_INFINITY) };
        tree.min_points[node] = mn;
        tree.max_points[node] = mx;
    }

    /// Self-pair query — OCCT BOPTools_PairSelector with same BVH sets.
    pub fn self_pairs(&self) -> Vec<(usize, usize)> {
        PairSelector::select(self)
    }

    /// Dual-tree pair query — OCCT BOPTools_PairSelector.
    pub fn candidate_pairs(bvh_a: &BoxTree, bvh_b: &BoxTree) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if bvh_a.tree.is_empty() || bvh_b.tree.is_empty() { return pairs; }
        Self::cross_pairs(bvh_a, 0, bvh_b, 0, &mut pairs);
        pairs
    }

    fn cross_pairs(bvh_a: &BoxTree, na: usize, bvh_b: &BoxTree, nb: usize, pairs: &mut Vec<(usize, usize)>) {
        if !Self::aabb_overlap(
            bvh_a.tree.min_point(na), bvh_a.tree.max_point(na),
            bvh_b.tree.min_point(nb), bvh_b.tree.max_point(nb),
        ) { return; }

        let a_is_leaf = bvh_a.tree.is_outer(na);
        let b_is_leaf = bvh_b.tree.is_outer(nb);

        if a_is_leaf && b_is_leaf {
            let sa = bvh_a.tree.beg_primitive(na) as usize;
            let ea = bvh_a.tree.end_primitive(na) as usize;
            let sb = bvh_b.tree.beg_primitive(nb) as usize;
            let eb = bvh_b.tree.end_primitive(nb) as usize;
            for i in sa..=ea {
                for j in sb..=eb {
                    if bvh_a.aabbs[i].intersects(&bvh_b.aabbs[j]) {
                        pairs.push((bvh_a.indices[i], bvh_b.indices[j]));
                    }
                }
            }
            return;
        }

        if !a_is_leaf {
            let la = BvhTreeBase::left_child(na);
            let ra = BvhTreeBase::right_child(na);
            if la < bvh_a.tree.length() { Self::cross_pairs(bvh_a, la, bvh_b, nb, pairs); }
            if ra < bvh_a.tree.length() { Self::cross_pairs(bvh_a, ra, bvh_b, nb, pairs); }
        } else {
            let lb = BvhTreeBase::left_child(nb);
            let rb = BvhTreeBase::right_child(nb);
            if lb < bvh_b.tree.length() { Self::cross_pairs(bvh_a, na, bvh_b, lb, pairs); }
            if rb < bvh_b.tree.length() { Self::cross_pairs(bvh_a, na, bvh_b, rb, pairs); }
        }
    }

    /// Query all items whose AABB overlaps the query AABB — OCCT BOPTools_BoxSelector.
    pub fn query_aabb(&self, query_min: &DVec3, query_max: &DVec3) -> Vec<usize> {
        let mut res = Vec::new();
        if self.tree.is_empty() { return res; }
        self.query_node(0, query_min, query_max, &mut res);
        res
    }

    fn query_node(&self, node: usize, qmin: &DVec3, qmax: &DVec3, res: &mut Vec<usize>) {
        if !Self::aabb_overlap(self.tree.min_point(node), self.tree.max_point(node), qmin, qmax) { return; }
        if self.tree.is_outer(node) {
            let s = self.tree.beg_primitive(node) as usize;
            let e = self.tree.end_primitive(node) as usize;
            for i in s..=e {
                let a = &self.aabbs[i];
                if a.min.x <= qmax.x && a.max.x >= qmin.x
                    && a.min.y <= qmax.y && a.max.y >= qmin.y
                    && a.min.z <= qmax.z && a.max.z >= qmin.z
                { res.push(self.indices[i]); }
            }
            return;
        }
        let l = BvhTreeBase::left_child(node);
        let r = BvhTreeBase::right_child(node);
        if l < self.tree.length() { self.query_node(l, qmin, qmax, res); }
        if r < self.tree.length() { self.query_node(r, qmin, qmax, res); }
    }

    fn aabb_overlap(min1: &DVec3, max1: &DVec3, min2: &DVec3, max2: &DVec3) -> bool {
        min1.x <= max2.x && max1.x >= min2.x
            && min1.y <= max2.y && max1.y >= min2.y
            && min1.z <= max2.z && max1.z >= min2.z
    }
}

// ── OCCT BOPTools_PairSelector — self-pair traversal ──
struct PairSelector;

impl PairSelector {
    fn select(tree: &BoxTree) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if tree.tree.is_empty() { return pairs; }
        Self::self_query(tree, 0, &mut pairs);
        pairs
    }

    fn self_query(tree: &BoxTree, node: usize, pairs: &mut Vec<(usize, usize)>) {
        if tree.tree.is_outer(node) {
            let s = tree.tree.beg_primitive(node) as usize;
            let e = tree.tree.end_primitive(node) as usize;
            for i in s..=e {
                for j in (i + 1)..=e {
                    if tree.aabbs[i].intersects(&tree.aabbs[j]) {
                        pairs.push((tree.indices[i], tree.indices[j]));
                    }
                }
            }
            return;
        }
        let l = BvhTreeBase::left_child(node);
        let r = BvhTreeBase::right_child(node);
        if l < tree.tree.length() { Self::self_query(tree, l, pairs); }
        if r < tree.tree.length() { Self::self_query(tree, r, pairs); }
        Self::cross_nodes(tree, l, r, pairs);
    }

    fn cross_nodes(tree: &BoxTree, na: usize, nb: usize, pairs: &mut Vec<(usize, usize)>) {
        if na >= tree.tree.length() || nb >= tree.tree.length() { return; }
        if !BoxTree::aabb_overlap(
            tree.tree.min_point(na), tree.tree.max_point(na),
            tree.tree.min_point(nb), tree.tree.max_point(nb),
        ) { return; }

        let a_leaf = tree.tree.is_outer(na);
        let b_leaf = tree.tree.is_outer(nb);

        if a_leaf && b_leaf {
            let sa = tree.tree.beg_primitive(na) as usize;
            let ea = tree.tree.end_primitive(na) as usize;
            let sb = tree.tree.beg_primitive(nb) as usize;
            let eb = tree.tree.end_primitive(nb) as usize;
            for i in sa..=ea {
                for j in sb..=eb {
                    if i >= j { continue; }
                    if tree.aabbs[i].intersects(&tree.aabbs[j]) {
                        pairs.push((tree.indices[i], tree.indices[j]));
                    }
                }
            }
            return;
        }

        if !a_leaf {
            Self::cross_nodes(tree, BvhTreeBase::left_child(na), nb, pairs);
            Self::cross_nodes(tree, BvhTreeBase::right_child(na), nb, pairs);
        } else {
            Self::cross_nodes(tree, na, BvhTreeBase::left_child(nb), pairs);
            Self::cross_nodes(tree, na, BvhTreeBase::right_child(nb), pairs);
        }
    }
}
