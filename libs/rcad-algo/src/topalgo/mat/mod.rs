//! OCCT TKTopAlgo/MAT — the MAT graph over the bisecting locus
//! (consumed by MAT2d), 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!
//! Handle mapping (architecture layer): OCCT `occ::handle<MAT_X>` maps to
//! `Arc<RwLock<MatX>>` (Standard_Transient equivalent, see bop/ds/pave.rs
//! SharedPB). OCCT handle `==` compares pointer identity -> `Arc::ptr_eq`.
//! OCCT null handle -> `Option<Arc<...>>`. OCCT stores some back-pointers as
//! raw `void*` to avoid handle cycles; here they are strong handles (the
//! graph is a short-lived structure, the OCCT destructors that break the
//! cycles are translated as `Drop` where present).

use std::sync::{Arc, RwLock};

pub mod mat_arc;
pub mod mat_basic_elt;
pub mod mat_bisector;
pub mod mat_edge;
pub mod mat_graph;
pub mod mat_list_of_bisector_0;
pub mod mat_list_of_edge_0;
pub mod mat_node;
pub mod mat_tlist_node_of_list_of_bisector_0;
pub mod mat_tlist_node_of_list_of_edge_0;
pub mod mat_zone;

pub use mat_arc::{HandleMatArc, MatArc};
pub use mat_basic_elt::{HandleMatBasicElt, MatBasicElt};
pub use mat_bisector::{HandleMatBisector, MatBisector};
pub use mat_edge::{HandleMatEdge, MatEdge};
pub use mat_graph::{HandleMatGraph, MatGraph};
pub use mat_list_of_bisector_0::{HandleMatListOfBisector, MatListOfBisector};
pub use mat_list_of_edge_0::{HandleMatListOfEdge, MatListOfEdge};
pub use mat_node::{HandleMatNode, MatNode};
pub use mat_tlist_node_of_list_of_bisector_0::{
    HandleMatTListNodeOfListOfBisector, MatTListNodeOfListOfBisector,
};
pub use mat_tlist_node_of_list_of_edge_0::{
    HandleMatTListNodeOfListOfEdge, MatTListNodeOfListOfEdge,
};
pub use mat_zone::{HandleMatZone, MatZone};

// OCCT MAT_Side.hxx L17-27
/// Definition on the Left and the Right on the Fig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatSide {
    Left,  // MAT_Left
    Right, // MAT_Right
}

/// OCCT handle `==` compares pointer identity (null == null is true).
/// Arc's PartialEq would delegate to the inner value, so all identity
/// comparisons go through this helper (architecture layer).
pub(crate) fn same_handle<T>(a: &Option<Arc<RwLock<T>>>, b: &Option<Arc<RwLock<T>>>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => Arc::ptr_eq(x, y),
        (None, None) => true,
        _ => false,
    }
}
