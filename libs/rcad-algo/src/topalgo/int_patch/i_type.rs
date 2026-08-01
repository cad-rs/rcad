// OCCT IntPatch_IType — the type of geometry of an IntPatch_Line.
//
// OCCT IntPatch_IType.hxx L20-30.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntPatchIType {
    /// Line.
    Lin,
    /// Circle.
    Circle,
    /// Ellipse.
    Ellipse,
    /// Parabola.
    Parabola,
    /// Hyperbola.
    Hyperbola,
    /// Analytic.
    Analytic,
    /// Walking.
    Walking,
    /// Restriction.
    Restriction,
}
