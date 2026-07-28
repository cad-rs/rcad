//! Re-exports: Bnd_Box, BVH_Tree/BVH_Builder, BOPTools_BoxTree
//!
//! | OCCT | rcad |
//! |------|------|
//! | `Bnd_Box` | `bnd_box::Aabb` |
//! | `BVH_Tree` / `BVH_Builder` | `rcad_kernel::math::bvh::Bvh` |
//! | `BOPTools_BoxTree` | `box_tree::BoxTree` |

pub use crate::bnd_box::Aabb;
pub use crate::boptools::box_tree::BoxTree;

/// BVH moved to rcad_kernel::math::bvh. Re-exported for backward compat.
pub use rcad_kernel::math::bvh::{Bvh, BvhStats};
