// OCCT BVH_TreeBase — bounding volume hierarchy with SOA node storage.
//
// OCCT BVH_Tree.hxx / BVH_TreeBase.hxx
// Nodes stored in parallel arrays (SOA):
//   myNodeInfoBuffer[i].x = isOuter (0 = inner, non-zero = leaf)
//   myNodeInfoBuffer[i].y = childLeft OR startPrim (leaf: first element index)
//   myNodeInfoBuffer[i].z = childRight OR endPrim   (leaf: last element index)
//   myNodeInfoBuffer[i].w = level     (tree depth level)
//   myMinPointBuffer[i]    = node min point
//   myMaxPointBuffer[i]    = node max point

use glam::{DVec3, IVec4};

/// OCCT BVH_TreeBase — SOA node storage for bounding volume hierarchy.
#[derive(Debug, Clone)]
pub(crate) struct BvhTreeBase {
    /// (isOuter, childLeft/startPrim, childRight/endPrim, level) per node
    pub node_info: Vec<IVec4>,
    /// Min corner of each node's bounding box
    pub min_points: Vec<DVec3>,
    /// Max corner of each node's bounding box
    pub max_points: Vec<DVec3>,
}

impl BvhTreeBase {
    pub fn new() -> Self {
        BvhTreeBase { node_info: Vec::new(), min_points: Vec::new(), max_points: Vec::new() }
    }

    pub fn depth(&self) -> i32 {
        if self.node_info.is_empty() { 0 } else { self.node_info.iter().map(|n| n.w).max().unwrap_or(0) }
    }
    pub fn length(&self) -> usize { self.node_info.len() }
    pub fn is_empty(&self) -> bool { self.node_info.is_empty() }
    pub fn min_point(&self, node_idx: usize) -> &DVec3 { &self.min_points[node_idx] }
    pub fn max_point(&self, node_idx: usize) -> &DVec3 { &self.max_points[node_idx] }

    // OCCT: BegPrimitive / EndPrimitive (leaf only)
    pub fn beg_primitive(&self, node_idx: usize) -> i32 { self.node_info[node_idx].y }
    pub fn end_primitive(&self, node_idx: usize) -> i32 { self.node_info[node_idx].z }
    pub fn nb_primitives(&self, node_idx: usize) -> i32 { self.end_primitive(node_idx) - self.beg_primitive(node_idx) + 1 }

    // OCCT: IsOuter
    pub fn is_outer(&self, node_idx: usize) -> bool { self.node_info[node_idx].x != 0 }
    pub fn level(&self, node_idx: usize) -> i32 { self.node_info[node_idx].w }

    // Internal node: explicit child indices (y=left, z=right)
    pub fn child_left(&self, node_idx: usize) -> usize { self.node_info[node_idx].y as usize }
    pub fn child_right(&self, node_idx: usize) -> usize { self.node_info[node_idx].z as usize }

    pub fn push_node(&mut self, is_outer: bool, v1: i32, v2: i32, level: i32,
                     min_p: DVec3, max_p: DVec3) -> usize {
        let idx = self.node_info.len();
        self.node_info.push(IVec4::new(if is_outer { 1 } else { 0 }, v1, v2, level));
        self.min_points.push(min_p);
        self.max_points.push(max_p);
        idx
    }
}
