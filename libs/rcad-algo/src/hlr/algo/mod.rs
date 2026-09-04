//! OCCT HLRAlgo package (TKHLR) — projector, packed min-max boxes and the
//! core hidden-line data structures.
//!
//! Ported classes (Stage 0):
//! - [`edges_block::EdgesBlock`]     (HLRAlgo_EdgesBlock)
//! - [`edges_block::MinMaxIndices`]  (HLRAlgo_EdgesBlock::MinMaxIndices)
//! - [`hlr_algo::HLRAlgo`]           (HLRAlgo static helpers)
//! - [`projector::Projector`]        (HLRAlgo_Projector)

pub mod edges_block;
pub mod hlr_algo;
pub mod projector;

pub use edges_block::{EdgesBlock, MinMaxIndices};
pub use hlr_algo::HLRAlgo;
pub use projector::Projector;
