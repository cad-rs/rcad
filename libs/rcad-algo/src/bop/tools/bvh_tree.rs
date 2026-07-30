//! Shared BvhNode type used by BoxTree (BOPTools_BoxTree).
//!
//! The SAH-based `Bvh` (BVH_Tree / BVH_Builder) moved to
//! `rcad_kernel::math::bvh` — use that for new code.
//!
//! This file retains only the `BvhNode` enum used by `box_tree.rs`.

use crate::bop::tools::bvh::Aabb;

/// BVH node (internal or leaf).
#[derive(Debug, Clone)]
pub(crate) enum BvhNode {
    Leaf { aabb: Aabb, start: usize, end: usize },
    Internal { aabb: Aabb, left: usize, right: usize },
}

impl BvhNode {
    pub(crate) fn aabb(&self) -> &Aabb {
        match self {
            BvhNode::Leaf { aabb, .. } => aabb,
            BvhNode::Internal { aabb, .. } => aabb,
        }
    }
}
