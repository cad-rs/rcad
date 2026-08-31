//! OCCT TKShHealing — ShapeFix / ShapeAnalysis / ShapeBuild / ShapeCustom
//! healing packages.
//!
//! Current content migrated from `algo_ext` (legacy `rcad-algorithms`
//! surface): the healing chain, the shape-analysis helpers and the
//! `ShapeCustom::restrict_to_bspline` equivalent.

pub mod healing;
pub mod shape_analysis;
pub mod shape_custom;

pub use healing::{HealingMode, HealingOptions, HealingReport, analyze_and_heal};
pub use shape_custom::restrict_to_bspline;
