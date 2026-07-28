//! OCCT IntPatch_IType.hxx — type of intersection line

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntPatchIType {
    Unknown,
    Line,
    Circle,
    Ellipse,
    Parabola,
    Hyperbola,
    Walking,
    Restriction,
}
