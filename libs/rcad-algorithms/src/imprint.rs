//! Stub: imprint_shape module removed (self-created code).
//! Kept to satisfy existing references until callers are cleaned up.

use rcad_kernel::BRep;

#[derive(Debug, Clone)]
pub struct ImprintResult {
    pub brep: BRep,
    pub seam_edges: Vec<usize>,
}

/// Stub: always returns a no-op result (original BRep unchanged).
pub fn imprint_shape(target: &BRep, _tool: &BRep) -> ImprintResult {
    ImprintResult {
        brep: target.clone(),
        seam_edges: vec![],
    }
}

/// Stub types for gap/overlap detection (deleted).
#[derive(Debug, Clone)]
pub struct Gap;
#[derive(Debug, Clone)]
pub struct Overlap;
#[derive(Debug, Clone)]
pub struct GapOverlapReport;
pub fn detect_gaps_overlaps(
    _a: &BRep, _b: &BRep, _tol: f64,
) -> GapOverlapReport {
    GapOverlapReport
}
pub fn min_distance(_a: &BRep, _b: &BRep) -> f64 { f64::MAX }
