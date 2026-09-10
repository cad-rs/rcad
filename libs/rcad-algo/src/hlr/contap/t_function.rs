// OCCT Contap_TFunction.hxx L20-26 — the type of the contour function.

/// OCCT Contap_TFunction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TFunction {
    /// OCCT Contap_ContourStd.
    ContourStd,
    /// OCCT Contap_ContourPrs.
    ContourPrs,
    /// OCCT Contap_DraftStd.
    DraftStd,
    /// OCCT Contap_DraftPrs.
    DraftPrs,
}
