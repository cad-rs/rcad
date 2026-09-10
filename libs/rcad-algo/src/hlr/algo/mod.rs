//! OCCT HLRAlgo package (TKHLR) — projector, packed min-max boxes and the
//! core hidden-line data structures.
//!
//! Ported classes (Stages 0-1):
//! - [`edges_block::EdgesBlock`]     (HLRAlgo_EdgesBlock)
//! - [`edges_block::MinMaxIndices`]  (HLRAlgo_EdgesBlock::MinMaxIndices)
//! - [`hlr_algo::HLRAlgo`]           (HLRAlgo static helpers)
//! - [`projector::Projector`]        (HLRAlgo_Projector)
//! - [`edge_status::EdgeStatus`]     (HLRAlgo_EdgeStatus)
//! - [`coincidence::Coincidence`]    (HLRAlgo_Coincidence)
//! - [`intersection::Intersection`]  (HLRAlgo_Intersection)
//! - [`interference::Interference`]  (HLRAlgo_Interference)
//! - [`bi_point::BiPoint`]           (HLRAlgo_BiPoint)
//! - [`wires_block::WiresBlock`]     (HLRAlgo_WiresBlock)
//! - [`poly_data_structs`]           (HLRAlgo_PolyMask / TriangleData /
//!   PolyInternalSegment / PolyInternalNode / PolyHidingData)
//! - [`edge_iterator::EdgeIterator`] (HLRAlgo_EdgeIterator)

pub mod bi_point;
pub mod coincidence;
pub mod edge_iterator;
pub mod edge_status;
pub mod edges_block;
pub mod hlr_algo;
pub mod interference;
pub mod intersection;
pub mod poly_data_structs;
pub mod projector;
pub mod wires_block;

pub use bi_point::{BiPoint, IndicesT, PointsT};
pub use coincidence::Coincidence;
pub use edge_iterator::EdgeIterator;
pub use edge_status::EdgeStatus;
pub use edges_block::{EdgesBlock, MinMaxIndices};
pub use hlr_algo::HLRAlgo;
pub use interference::Interference;
pub use intersection::Intersection;
pub use poly_data_structs::{
    HidingPlaneT, HidingTriangleIndices, NodeData, NodeIndices, PolyHidingData,
    PolyInternalNode, PolyInternalSegment, TriangleData,
};
pub use projector::Projector;
pub use wires_block::WiresBlock;
