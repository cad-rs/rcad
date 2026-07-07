//! OCCT IntPatch_Line.hxx + subtypes (ALine / WLine / GLine / RLine)
//!
//! rcad: single flat struct with line_type discriminator and optional
//! WLine data (point sequence for walking lines).

use super::int_patch_type::IntPatchIType;
use glam::DVec3;

/// A point on a walking line (WLine): 3D point + UV on both surfaces.
#[derive(Debug, Clone, Copy)]
pub struct WLinePnt {
    pub p3d: glam::DVec3,
    pub u1: f64, pub v1: f64,
    pub u2: f64, pub v2: f64,
}

/// OCCT-aligned: IntPatch_Line — intersection line.
///
/// Supports both ALine (analytic, stored as curve) and
/// WLine (walking, stored as point sequence) via wline_pnts.
#[derive(Debug, Clone)]
pub struct IntPatchLine {
    pub line_type: IntPatchIType,
    // ALine/analytic curve data
    pub curve: rcad_kernel::geom::Curve3,
    pub t_range: [f64; 2],
    pub pcurve1: Option<rcad_kernel::geom::Curve2d>,
    pub pcurve2: Option<rcad_kernel::geom::Curve2d>,
    pub tolerance: f64,
    pub tang_tolerance: f64,
    // WLine data (point sequence for walking lines)
    pub wline_pnts: Vec<WLinePnt>,
    // OCCT WLine flags
    pub is_purging_allowed: bool,
    pub wl_type: WLineType,
}

impl IntPatchLine {
    pub fn analytic(line_type: IntPatchIType, curve: rcad_kernel::geom::Curve3,
                    t_range: [f64; 2]) -> Self {
        Self {
            line_type, curve, t_range,
            pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7,
            wline_pnts: Vec::new(),
            is_purging_allowed: false,
            wl_type: WLineType::Unknown,
        }
    }

    pub fn walking(pnts: Vec<WLinePnt>, wl_type: WLineType) -> Self {
        Self {
            line_type: IntPatchIType::Walking,
            curve: rcad_kernel::geom::Curve3::Line(
                rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::X,
                }),
            t_range: [0.0, 1.0],
            pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7,
            wline_pnts: pnts,
            is_purging_allowed: true,
            wl_type,
        }
    }

/// OCCT IntPatch_WLine::IntPatch_WLType
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WLineType { Unknown, ImpImp, ImpPrm, PrmPrm }

impl IntPatchLine {
    pub fn new_analytic(line_type: IntPatchIType, curve: rcad_kernel::geom::Curve3,
                        t_range: [f64; 2]) -> Self {
        Self {
            line_type, curve, t_range,
            pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7,
            wline_pnts: Vec::new(),
            is_purging_allowed: false,
            wl_type: WLineType::Unknown,
        }
    }

    pub fn new_walking(pnts: Vec<WLinePnt>, wl_type: WLineType) -> Self {
        Self {
            line_type: IntPatchIType::Walking,
            curve: rcad_kernel::geom::Curve3::Line(
                rcad_kernel::geom::Line3 { origin: glam::DVec3::ZERO, direction: glam::DVec3::X }),
            t_range: [0.0, 1.0],
            pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7,
            wline_pnts: pnts,
            is_purging_allowed: true,
            wl_type,
        }
    }

    pub fn is_wline(&self) -> bool { !self.wline_pnts.is_empty() }
    pub fn nb_points(&self) -> usize { self.wline_pnts.len() }
    pub fn point(&self, i: usize) -> &WLinePnt { &self.wline_pnts[i] }
}
