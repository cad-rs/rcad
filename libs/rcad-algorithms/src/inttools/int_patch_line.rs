use super::int_patch_type::IntPatchIType;
use glam::DVec3;

#[derive(Debug, Clone, Copy)]
pub struct WLinePnt {
    pub p3d: DVec3,
    pub u1: f64, pub v1: f64,
    pub u2: f64, pub v2: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WLineType { Unknown, ImpImp, ImpPrm, PrmPrm }

#[derive(Debug, Clone)]
pub struct IntPatchLine {
    pub line_type: IntPatchIType,
    pub curve: rcad_kernel::geom::Curve3,
    pub t_range: [f64; 2],
    pub pcurve1: Option<rcad_kernel::geom::Curve2d>,
    pub pcurve2: Option<rcad_kernel::geom::Curve2d>,
    pub tolerance: f64,
    pub tang_tolerance: f64,
    pub wline_pnts: Vec<WLinePnt>,
    pub is_purging_allowed: bool,
    pub wl_type: WLineType,
}

impl IntPatchLine {
    pub fn analytic(lt: IntPatchIType, curve: rcad_kernel::geom::Curve3, tr: [f64; 2]) -> Self {
        Self { line_type: lt, curve, t_range: tr, pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7, wline_pnts: Vec::new(),
            is_purging_allowed: false, wl_type: WLineType::Unknown }
    }
    pub fn walking(pnts: Vec<WLinePnt>, wt: WLineType) -> Self {
        let line = rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X };
        Self { line_type: IntPatchIType::Walking, curve: rcad_kernel::geom::Curve3::Line(line),
            t_range: [0.0, 1.0], pcurve1: None, pcurve2: None,
            tolerance: 1e-7, tang_tolerance: 1e-7, wline_pnts: pnts,
            is_purging_allowed: true, wl_type: wt }
    }
    pub fn is_wline(&self) -> bool { !self.wline_pnts.is_empty() }
    pub fn nb_points(&self) -> usize { self.wline_pnts.len() }
    pub fn point(&self, i: usize) -> &WLinePnt { &self.wline_pnts[i] }
}
