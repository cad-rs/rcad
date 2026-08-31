//! OCCT TKShHealing — ShapeFix / ShapeAnalysis / ShapeBuild / ShapeCustom
//! healing packages.
//!
//! Current content migrated from `algo_ext` (legacy `rcad-algorithms`
//! surface): the healing chain, the shape-analysis helpers and the
//! `ShapeCustom::restrict_to_bspline` equivalent — plus the 1:1
//! `ShapeExtend`/`ShapeBuild` foundation classes for the ShapeFix stack
//! (wire_tails heal grid pilot).

pub mod healing;
pub mod shape_analysis;
pub mod shape_build;
pub mod shape_custom;
pub mod shape_extend;

pub use healing::{HealingMode, HealingOptions, HealingReport, analyze_and_heal};
pub use shape_custom::restrict_to_bspline;
