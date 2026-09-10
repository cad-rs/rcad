//! BVH for DS entity pair culling — OCCT BOPTools_BoxTree / BOPTools_BoxSelector.
//!
//! OCCT correspondence:
//! - BOPTools_BoxSet / BOPTools_BoxTree  = BoxTree
//! - BVH_LinearBuilder                     = build() method
//! - BVH_RadixSorter                       = morton sort inside build()
//! - BOPTools_PairSelector                 = PairSelector
//! - BOPTools_BoxSelector                  = query_aabb()

use crate::bop::tools::bvh_tree::BvhTreeBase;
use glam::DVec3;

/// AABB used during BVH build (OCCT Bnd_Box equivalent via Bnd2BVH).
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

const MORTON_LUT: [u32; 256] = [
    0x000000,0x000001,0x000008,0x000009,0x000040,0x000041,0x000048,0x000049,0x000200,0x000201,
    0x000208,0x000209,0x000240,0x000241,0x000248,0x000249,0x001000,0x001001,0x001008,0x001009,
    0x001040,0x001041,0x001048,0x001049,0x001200,0x001201,0x001208,0x001209,0x001240,0x001241,
    0x001248,0x001249,0x008000,0x008001,0x008008,0x008009,0x008040,0x008041,0x008048,0x008049,
    0x008200,0x008201,0x008208,0x008209,0x008240,0x008241,0x008248,0x008249,0x009000,0x009001,
    0x009008,0x009009,0x009040,0x009041,0x009048,0x009049,0x009200,0x009201,0x009208,0x009209,
    0x009240,0x009241,0x009248,0x009249,0x040000,0x040001,0x040008,0x040009,0x040040,0x040041,
    0x040048,0x040049,0x040200,0x040201,0x040208,0x040209,0x040240,0x040241,0x040248,0x040249,
    0x041000,0x041001,0x041008,0x041009,0x041040,0x041041,0x041048,0x041049,0x041200,0x041201,
    0x041208,0x041209,0x041240,0x041241,0x041248,0x041249,0x048000,0x048001,0x048008,0x048009,
    0x048040,0x048041,0x048048,0x048049,0x048200,0x048201,0x048208,0x048209,0x048240,0x048241,
    0x048248,0x048249,0x049000,0x049001,0x049008,0x049009,0x049040,0x049041,0x049048,0x049049,
    0x049200,0x049201,0x049208,0x049209,0x049240,0x049241,0x049248,0x049249,0x200000,0x200001,
    0x200008,0x200009,0x200040,0x200041,0x200048,0x200049,0x200200,0x200201,0x200208,0x200209,
    0x200240,0x200241,0x200248,0x200249,0x201000,0x201001,0x201008,0x201009,0x201040,0x201041,
    0x201048,0x201049,0x201200,0x201201,0x201208,0x201209,0x201240,0x201241,0x201248,0x201249,
    0x208000,0x208001,0x208008,0x208009,0x208040,0x208041,0x208048,0x208049,0x208200,0x208201,
    0x208208,0x208209,0x208240,0x208241,0x208248,0x208249,0x209000,0x209001,0x209008,0x209009,
    0x209040,0x209041,0x209048,0x209049,0x209200,0x209201,0x209208,0x209209,0x209240,0x209241,
    0x209248,0x209249,0x240000,0x240001,0x240008,0x240009,0x240040,0x240041,0x240048,0x240049,
    0x240200,0x240201,0x240208,0x240209,0x240240,0x240241,0x240248,0x240249,0x241000,0x241001,
    0x241008,0x241009,0x241040,0x241041,0x241048,0x241049,0x241200,0x241201,0x241208,0x241209,
    0x241240,0x241241,0x241248,0x241249,0x248000,0x248001,0x248008,0x248009,0x248040,0x248041,
    0x248048,0x248049,0x248200,0x248201,0x248208,0x248209,0x248240,0x248241,0x248248,0x248249,
    0x249000,0x249001,0x249008,0x249009,0x249040,0x249041,0x249048,0x249049,0x249200,0x249201,
    0x249208,0x249209,0x249240,0x249241,0x249248,0x249249,
];

fn encode_morton(vx: u32, vy: u32, vz: u32) -> u32 {
    (MORTON_LUT[(vx & 0xFF) as usize] | (MORTON_LUT[((vx >> 8) & 0x03) as usize] << 24))
        | ((MORTON_LUT[(vy & 0xFF) as usize] | (MORTON_LUT[((vy >> 8) & 0x03) as usize] << 24)) << 1)
        | ((MORTON_LUT[(vz & 0xFF) as usize] | (MORTON_LUT[((vz >> 8) & 0x03) as usize] << 24)) << 2)
}

fn radix_sort_msd(links: &mut [(u32, usize)], start: usize, end: usize, digit: i32) {
    if end - start <= 4 || digit < 0 { return; }
    let bit = 1u32 << (digit as u32); let mut split = start;
    for i in start..end { if links[i].0 & bit == 0 { links.swap(i, split); split += 1; } }
    if split > start && split < end { radix_sort_msd(links, start, split, digit - 1); radix_sort_msd(links, split, end, digit - 1); }
    else { radix_sort_msd(links, start, end, digit - 1); }
}

fn emit_hierarchy(links: &[(u32, usize)], start: usize, end: usize,
    nodes: &mut Vec<(i32, i32, i32, i32)>, ns: &mut usize, depth: i32) -> usize {
    if end - start <= 4 { let idx = *ns; *ns += 1; nodes.push((1, start as i32, (end - 1) as i32, depth)); return idx; }
    let diff = links[start].0 ^ links[end - 1].0;
    let hb = (0..=29).rev().find(|&b| diff & (1u32 << b) != 0);
    let split = match hb {
        None => (start + end) / 2,
        Some(bit) => { let mask = 1u32 << bit; let mut lo = start; let mut hi = end; while lo < hi { let mid = (lo + hi) / 2; if links[mid].0 & mask == 0 { lo = mid + 1; } else { hi = mid; } } lo }
    };
    if split == start || split == end { let idx = *ns; *ns += 1; nodes.push((1, start as i32, (end - 1) as i32, depth)); return idx; }
    let idx = *ns; *ns += 1; nodes.push((0, 0, 0, depth));
    let lc = emit_hierarchy(links, start, split, nodes, ns, depth + 1);
    let rc = emit_hierarchy(links, split, end, nodes, ns, depth + 1);
    nodes[idx] = (0, lc as i32, rc as i32, depth);
    idx
}

fn compute_bounds(node: usize, raw: &[(i32, i32, i32, i32)], aabbs: &[Aabb], tree: &mut BvhTreeBase) {
    let (is_outer, v1, v2, _) = raw[node];
    if is_outer != 0 {
        let mut mn = aabbs[v1 as usize].min; let mut mx = aabbs[v1 as usize].max;
        for i in (v1 as usize + 1)..=v2 as usize { mn = mn.min(aabbs[i].min); mx = mx.max(aabbs[i].max); }
        tree.min_points[node] = mn; tree.max_points[node] = mx; return;
    }
    compute_bounds(v1 as usize, raw, aabbs, tree); compute_bounds(v2 as usize, raw, aabbs, tree);
    let mn = tree.min_points[v1 as usize].min(tree.min_points[v2 as usize]);
    let mx = tree.max_points[v1 as usize].max(tree.max_points[v2 as usize]);
    tree.min_points[node] = mn; tree.max_points[node] = mx;
}

/// OCCT BOPTools_BoxSet — BVH set with element array + BVH tree.
pub struct BoxTree {
    pub(crate) tree: BvhTreeBase,
    pub(crate) indices: Vec<usize>,
    pub(crate) aabbs: Vec<Aabb>,
}

impl BoxTree {
    pub fn new() -> Self { BoxTree { tree: BvhTreeBase::new(), indices: Vec::new(), aabbs: Vec::new() } }
    pub fn set_size(&mut self, n: usize) { self.indices.reserve(n); self.aabbs.reserve(n); }
    pub fn add(&mut self, idx: usize, aabb: Aabb) { self.indices.push(idx); self.aabbs.push(aabb); }

    pub fn build(&mut self) {
        let n = self.indices.len();
        if n == 0 { self.tree = BvhTreeBase::new(); return; }
        let mut smin = self.aabbs[0].min; let mut smax = self.aabbs[0].max;
        for a in &self.aabbs { smin = smin.min(a.min); smax = smax.max(a.max); }
        let extent = (smax - smin).max(DVec3::splat(1e-12));
        let inv_ext = DVec3::new(1024.0 / extent.x, 1024.0 / extent.y, 1024.0 / extent.z);
        let mut links: Vec<(u32, usize)> = Vec::with_capacity(n);
        for i in 0..n { let c = self.aabbs[i].center(); let vf = (c - smin) * inv_ext; let vx = (vf.x as i32).clamp(0, 1023) as u32; let vy = (vf.y as i32).clamp(0, 1023) as u32; let vz = (vf.z as i32).clamp(0, 1023) as u32; links.push((encode_morton(vx, vy, vz), i)); }
        radix_sort_msd(&mut links, 0, n, 29);
        let mut si = Vec::with_capacity(n); let mut sa = Vec::with_capacity(n);
        for &(_, oi) in &links { si.push(self.indices[oi]); sa.push(self.aabbs[oi]); }
        let mut raw: Vec<(i32, i32, i32, i32)> = Vec::new(); let mut ns = 0;
        emit_hierarchy(&links, 0, n, &mut raw, &mut ns, 0);
        let mut tree = BvhTreeBase::new();
        for i in 0..raw.len() { let (io, v1, v2, lv) = raw[i]; tree.push_node(io != 0, v1, v2, lv, DVec3::ZERO, DVec3::ZERO); }
        compute_bounds(0, &raw, &sa, &mut tree);
        self.tree = tree; self.indices = si; self.aabbs = sa;
    }

    pub fn candidate_pairs(bvh_a: &BoxTree, bvh_b: &BoxTree) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if bvh_a.tree.is_empty() || bvh_b.tree.is_empty() { return pairs; }
        Self::cross_pairs(bvh_a, 0, bvh_b, 0, &mut pairs);
        pairs
    }
    fn cross_pairs(bvh_a: &BoxTree, na: usize, bvh_b: &BoxTree, nb: usize, pairs: &mut Vec<(usize, usize)>) {
        if bvh_a.tree.min_point(na).x > bvh_b.tree.max_point(nb).x || bvh_a.tree.max_point(na).x < bvh_b.tree.min_point(nb).x
            || bvh_a.tree.min_point(na).y > bvh_b.tree.max_point(nb).y || bvh_a.tree.max_point(na).y < bvh_b.tree.min_point(nb).y
            || bvh_a.tree.min_point(na).z > bvh_b.tree.max_point(nb).z || bvh_a.tree.max_point(na).z < bvh_b.tree.min_point(nb).z { return; }
        let a_leaf = bvh_a.tree.is_outer(na); let b_leaf = bvh_b.tree.is_outer(nb);
        if a_leaf && b_leaf {
            for i in bvh_a.tree.beg_primitive(na) as usize..=bvh_a.tree.end_primitive(na) as usize {
                for j in bvh_b.tree.beg_primitive(nb) as usize..=bvh_b.tree.end_primitive(nb) as usize {
                    if bvh_a.aabbs[i].intersects(&bvh_b.aabbs[j]) { pairs.push((bvh_a.indices[i], bvh_b.indices[j])); }
                }
            } return;
        }
        if !a_leaf { Self::cross_pairs(bvh_a, bvh_a.tree.child_left(na), bvh_b, nb, pairs); Self::cross_pairs(bvh_a, bvh_a.tree.child_right(na), bvh_b, nb, pairs); }
        else { Self::cross_pairs(bvh_a, na, bvh_b, bvh_b.tree.child_left(nb), pairs); Self::cross_pairs(bvh_a, na, bvh_b, bvh_b.tree.child_right(nb), pairs); }
    }

    pub fn query_aabb(&self, qmin: &DVec3, qmax: &DVec3) -> Vec<usize> {
        let mut res = Vec::new();
        if self.tree.is_empty() { return res; }
        Self::query_node(self, 0, qmin, qmax, &mut res);
        res
    }
    fn query_node(&self, node: usize, qmin: &DVec3, qmax: &DVec3, res: &mut Vec<usize>) {
        if self.tree.min_point(node).x > qmax.x || self.tree.max_point(node).x < qmin.x
            || self.tree.min_point(node).y > qmax.y || self.tree.max_point(node).y < qmin.y
            || self.tree.min_point(node).z > qmax.z || self.tree.max_point(node).z < qmin.z { return; }
        if self.tree.is_outer(node) {
            for i in self.tree.beg_primitive(node) as usize..=self.tree.end_primitive(node) as usize {
                let a = &self.aabbs[i];
                if a.min.x <= qmax.x && a.max.x >= qmin.x && a.min.y <= qmax.y && a.max.y >= qmin.y && a.min.z <= qmax.z && a.max.z >= qmin.z { res.push(self.indices[i]); }
            } return;
        }
        Self::query_node(self, self.tree.child_left(node), qmin, qmax, res);
        Self::query_node(self, self.tree.child_right(node), qmin, qmax, res);
    }
}

/// OCCT BOPTools_PairSelector — selects overlapping element pairs from BVH trees.
pub struct PairSelector {
    tree: *const BoxTree,
    same: bool,
    pairs: Vec<(usize, usize)>,
}

impl PairSelector {
    pub fn new() -> Self { PairSelector { tree: std::ptr::null(), same: false, pairs: Vec::new() } }
    pub fn set_bvh_sets(&mut self, tree: &BoxTree) { self.tree = tree as *const BoxTree; }
    pub fn set_same(&mut self, is_same: bool) { self.same = is_same; }

    pub fn select(&mut self) {
        self.pairs.clear();
        let tree = unsafe { &*self.tree };
        if tree.tree.is_empty() { return; }
        self.self_query(tree, 0);
    }

    fn self_query(&mut self, tree: &BoxTree, node: usize) {
        if tree.tree.is_outer(node) {
            for i in tree.tree.beg_primitive(node) as usize..=tree.tree.end_primitive(node) as usize {
                for j in (i + 1)..=tree.tree.end_primitive(node) as usize {
                    if tree.aabbs[i].intersects(&tree.aabbs[j]) { self.pairs.push((tree.indices[i], tree.indices[j])); }
                }
            } return;
        }
        self.self_query(tree, tree.tree.child_left(node));
        self.self_query(tree, tree.tree.child_right(node));
        self.cross_nodes(tree, tree.tree.child_left(node), tree.tree.child_right(node));
    }

    fn cross_nodes(&mut self, tree: &BoxTree, na: usize, nb: usize) {
        if na >= tree.tree.length() || nb >= tree.tree.length() || na == nb { return; }
        if tree.tree.min_point(na).x > tree.tree.max_point(nb).x || tree.tree.max_point(na).x < tree.tree.min_point(nb).x
            || tree.tree.min_point(na).y > tree.tree.max_point(nb).y || tree.tree.max_point(na).y < tree.tree.min_point(nb).y
            || tree.tree.min_point(na).z > tree.tree.max_point(nb).z || tree.tree.max_point(na).z < tree.tree.min_point(nb).z { return; }
        let a_leaf = tree.tree.is_outer(na); let b_leaf = tree.tree.is_outer(nb);
        if a_leaf && b_leaf {
            // Every (primitive of na, primitive of nb) pair is distinct: the
            // recursion visits each subtree pair exactly once, so no index
            // dedup applies here (the Morton-sorted leaf ranges need not be
            // ordered, and an `i >= j` skip would drop valid pairs).
            for i in tree.tree.beg_primitive(na) as usize..=tree.tree.end_primitive(na) as usize {
                for j in tree.tree.beg_primitive(nb) as usize..=tree.tree.end_primitive(nb) as usize {
                    if tree.aabbs[i].intersects(&tree.aabbs[j]) { self.pairs.push((tree.indices[i], tree.indices[j])); }
                }
            } return;
        }
        if !a_leaf { let cl = tree.tree.child_left(na); if cl != na { self.cross_nodes(tree, cl, nb); }
                     let cr = tree.tree.child_right(na); if cr != na { self.cross_nodes(tree, cr, nb); } }
        else { let cl = tree.tree.child_left(nb); if cl != nb { self.cross_nodes(tree, na, cl); }
               let cr = tree.tree.child_right(nb); if cr != nb { self.cross_nodes(tree, na, cr); } }
    }

    pub fn sort(&mut self) { self.pairs.sort_unstable_by_key(|&(a, b)| (a, b)); }
    pub fn pairs(&self) -> &[(usize, usize)] { &self.pairs }
}

#[cfg(test)]
mod pair_selector_tests {
    //! The self-pair yield must match the brute-force overlap scan for any
    //! Morton ordering of the primitives (BOPTools_PairSelector contract).

    use super::{Aabb, BoxTree, PairSelector};
    use glam::DVec3;

    #[test]
    fn pair_yield_matches_brute_force() {
        let mut tree = BoxTree::new();
        tree.set_size(35);
        let mut boxes: Vec<(f64, f64, f64)> = Vec::new();
        for k in 0..30 {
            let a = k as f64 * std::f64::consts::TAU / 30.0;
            boxes.push((28.0 * a.cos(), 28.0 * a.sin(), 6.0));
        }
        boxes.push((2.0, 0.0, 5.0));
        boxes.push((0.0, 3.0, 5.0));
        boxes.push((30.0, 5.0, 4.0));
        boxes.push((60.0, 60.0, 2.0));
        boxes.push((15.0, 15.0, 1.0));
        for (i, (cx, cy, r)) in boxes.iter().enumerate() {
            tree.add(
                i,
                Aabb {
                    min: DVec3::new(cx - r, cy - r, 0.0),
                    max: DVec3::new(cx + r, cy + r, 0.0),
                    gap: 0.0,
                },
            );
        }
        tree.build();
        let mut selector = PairSelector::new();
        selector.set_bvh_sets(&tree);
        selector.set_same(true);
        selector.select();
        selector.sort();
        let got = selector.pairs().to_vec();

        let mut want: Vec<(usize, usize)> = Vec::new();
        for i in 0..boxes.len() {
            for j in (i + 1)..boxes.len() {
                let (cx_i, cy_i, r_i) = boxes[i];
                let (cx_j, cy_j, r_j) = boxes[j];
                let dx = cx_i - cx_j;
                let dy = cy_i - cy_j;
                let rr = r_i + r_j;
                if dx * dx + dy * dy <= rr * rr {
                    want.push((i, j));
                }
            }
        }
        // The selector yields each unordered pair once, in arbitrary order.
        let mut got: Vec<(usize, usize)> = got.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
        got.sort_unstable();
        want.sort_unstable();
        assert_eq!(got.len(), want.len(), "pair count mismatch");
        for (g, w) in got.iter().zip(want.iter()) {
            assert_eq!(g, w, "pair mismatch");
        }
    }
}
