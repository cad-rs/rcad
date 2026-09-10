//! Adaptor3d (TKG3d) — the surface-domain topological tool layer.
//!
//! 1:1 translations:
//! - [`hvertex::HVertex`] — Adaptor3d_HVertex.hxx/.cxx.
//! - [`topol_tool::TopolTool`] — Adaptor3d_TopolTool.hxx/.cxx (the default
//!   domain tool: restriction lines + classification + sampling).

pub mod hvertex;
pub mod topol_tool;

pub use hvertex::HVertex;
pub use topol_tool::{get_cone_apex_param, TopolTool};
