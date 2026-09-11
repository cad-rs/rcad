//! OCCT BRepSweep (TKPrim) — the sweeping primitives package.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKPrim/BRepSweep/
//! Consumers: LocOpe_Prism/Revol (feat/), BRepFeat_MakeDPrism/MakeRevol.
//!
//! Layout (per class; the .lxx inline bodies are folded into their class
//! modules):
//! - `sweep_num_shape` / `sweep_num_shape_iterator` / `sweep_num_shape_tool`
//!   — the TKPrim/Sweep num-simulation layer (the directing-shape services).
//! - `brep_sweep_iterator` / `brep_sweep_tool` / `brep_sweep_builder` — the
//!   generatrix services and the BRep_Builder adapter.
//! - `num_linear_regular_sweep` — the generic sweep engine (the shape grid)
//!   with the `NumLinearRegularSweepSlots` virtual table.
//! - `trsf` — the shared-transform intermediate class (Init / Process /
//!   SetContinuity) and its `BRepSweepTrsfSlots` extension trait.
//! - `translation` / `rotation` — the concrete sweeps (Prism / Revol engines).
//! - `prism` / `revol` — the natural constructors consumed by feat/.
//! - `tool_rehost` — the TKBRep / TKG2d / TKG3d / TKMath leaf re-hosts.

pub mod brep_sweep_builder;
pub mod brep_sweep_iterator;
pub mod brep_sweep_tool;
pub mod make_revol;
pub mod num_linear_regular_sweep;
pub mod prism;
pub mod revol;
pub mod rotation;
pub mod sweep_num_shape;
pub mod sweep_num_shape_iterator;
pub mod sweep_num_shape_tool;
pub mod tool_rehost;
pub mod translation;
pub mod trsf;

pub use brep_sweep_builder::BRepSweepBuilder;
pub use num_linear_regular_sweep::NumLinearRegularSweepSlots;
pub use prism::BRepSweepPrism;
pub use revol::BRepSweepRevol;
pub use rotation::BRepSweepRotation;
pub use sweep_num_shape::SweepNumShape;
pub use translation::BRepSweepTranslation;
