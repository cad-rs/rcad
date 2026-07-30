// OCCT BOPAlgo_CheckResult — result of argument analysis for Boolean Operations.
//
// OCCT BOPAlgo_CheckResult.hxx (standalone class).
// OCCT BOPAlgo_CheckStatus is defined in BOPAlgo_CheckStatus.hxx.

/// OCCT BOPAlgo_CheckStatus — status codes for CheckResult.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckStatus {
    BadType,
    SelfIntersect,
    TooSmallEdge,
    NonRecoverableFace,
    IncompatibilityOfVertex,
    IncompatibilityOfEdge,
    IncompatibilityOfFace,
    GeomAbsC0,
    InvalidCurveOnSurface,
    OperationAborted,
    CheckUnknown,
}

/// OCCT BOPAlgo_CheckResult — one validation result from ArgumentAnalyzer.
#[derive(Debug, Clone)]
pub struct CheckResult {
    pub check_status: CheckStatus,
}
