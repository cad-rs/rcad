// OCCT BiTgte_ContactType.hxx L24-33 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BiTgte/
//         BiTgte_ContactType.hxx
//
// OCCT inheritance chain (hxx L24): none — BiTgte_ContactType is a plain
// enum.

/// OCCT BiTgte_ContactType (BiTgte_ContactType.hxx L24-33).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BiTgteContactType {
    /// OCCT BiTgte_FaceFace.
    FaceFace,
    /// OCCT BiTgte_FaceEdge.
    FaceEdge,
    /// OCCT BiTgte_FaceVertex.
    FaceVertex,
    /// OCCT BiTgte_EdgeEdge.
    EdgeEdge,
    /// OCCT BiTgte_EdgeVertex.
    EdgeVertex,
    /// OCCT BiTgte_VertexVertex.
    VertexVertex,
}
