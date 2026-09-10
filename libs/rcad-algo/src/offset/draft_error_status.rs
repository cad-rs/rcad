//! OCCT Draft_ErrorStatus (Draft_ErrorStatus.hxx L20-26) — the error status
//! reported by Draft_Modification::Error.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/Draft_ErrorStatus.hxx

/// OCCT Draft_ErrorStatus (Draft_ErrorStatus.hxx L20-26).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DraftErrorStatus {
    NoError,
    FaceRecomputation,
    EdgeRecomputation,
    VertexRecomputation,
}
