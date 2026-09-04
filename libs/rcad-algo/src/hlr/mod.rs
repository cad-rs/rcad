//! OCCT TKHLR port — Hidden Line Removal (HLRAlgo + HLRBRep + HLRTopoBRep +
//! Contap + Intrv + TopBas + TopCnx packages).
//!
//! Runway plan and progress checklist: `rcad/docs/tkhlr-port-plan.md`.
//! Ported so far (Stage 0 leaf packages):
//! - [`intrv`]                (Intrv: Interval / Intervals / Position)
//! - [`top_bas`]              (TopBas: TestInterference)
//! - [`top_cnx`]              (TopCnx: EdgeFaceTransition)
//! - [`algo`]                 (HLRAlgo: statics, EdgesBlock, Projector)

pub mod algo;
pub mod intrv;
pub mod top_bas;
pub mod top_cnx;

#[cfg(test)]
pub mod tests;

pub use intrv::{Intervals, Interval as IntrvInterval, Position as IntrvPosition};
pub use top_bas::TestInterference;
pub use top_cnx::EdgeFaceTransition;
pub use algo::edges_block::{EdgesBlock, MinMaxIndices};
pub use algo::hlr_algo::HLRAlgo;
pub use algo::projector::Projector;
