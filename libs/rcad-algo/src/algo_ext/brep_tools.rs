//! `brep_tools` compatibility module.
//!
//! The OCCT test generator emits `rcad_algorithms::brep_tools::extract_solids`
//! / `extract_shells` for `explode ... so` / `explode ... Sh` DRAW commands;
//! keep that module path working by re-exporting the `bool_ops_ext` helpers.

pub use super::bool_ops_ext::{extract_shells, extract_solids};
