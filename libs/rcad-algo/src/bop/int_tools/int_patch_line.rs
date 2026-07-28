//! IntPatch_Line — intersection line between two surfaces.

/// Type of walking-line algorithm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WLineType { Analytic, Walking, ParamParam, ImpImp, ImpPrm }

/// Walking-line point.
#[derive(Debug, Clone)]
pub struct WLinePnt;

/// Intersection line.
#[derive(Debug, Clone)]
pub struct IntPatchLine;
