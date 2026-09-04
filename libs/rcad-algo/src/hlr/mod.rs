//! OCCT TKHLR port — Hidden Line Removal (HLRAlgo + HLRBRep + HLRTopoBRep +
//! Contap + Intrv + TopBas + TopCnx packages).
//!
//! Runway plan and progress checklist: `rcad/docs/tkhlr-port-plan.md`.
//! Ported so far (Stages 0-1):
//! - [`intrv`]                (Intrv: Interval / Intervals / Position)
//! - [`top_bas`]              (TopBas: TestInterference)
//! - [`top_cnx`]              (TopCnx: EdgeFaceTransition)
//! - [`algo`]                 (HLRAlgo: statics, EdgesBlock, Projector,
//!   EdgeStatus, Coincidence, Intersection, Interference, BiPoint,
//!   WiresBlock, poly data structures, EdgeIterator)

pub mod algo;
pub mod contap;
pub mod intrv;
pub mod top_bas;
pub mod top_cnx;

#[cfg(test)]
pub mod tests;

pub use intrv::{Intervals, Interval as IntrvInterval, Position as IntrvPosition};
pub use top_bas::TestInterference;
pub use top_cnx::EdgeFaceTransition;
pub use algo::bi_point::{BiPoint, IndicesT, PointsT};
pub use algo::coincidence::Coincidence;
pub use algo::edge_iterator::EdgeIterator;
pub use algo::edge_status::EdgeStatus;
pub use algo::edges_block::{EdgesBlock, MinMaxIndices};
pub use algo::hlr_algo::HLRAlgo;
pub use algo::interference::Interference;
pub use algo::intersection::Intersection;
pub use algo::poly_data_structs::{
    HidingPlaneT, HidingTriangleIndices, NodeData, NodeIndices, PolyHidingData,
    PolyInternalNode, PolyInternalSegment, TriangleData,
};
pub use algo::projector::Projector;
pub use algo::wires_block::WiresBlock;
