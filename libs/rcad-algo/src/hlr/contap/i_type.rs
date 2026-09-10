// OCCT Contap_IType.hxx L20-26 — the type of a Contap contour line.

/// OCCT Contap_IType.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IType {
    /// OCCT Contap_Lin.
    Lin,
    /// OCCT Contap_Circle.
    Circle,
    /// OCCT Contap_Walking.
    Walking,
    /// OCCT Contap_Restriction.
    Restriction,
}
