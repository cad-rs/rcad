//! IntPatch_IType — type of intersection line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntPatchIType {
    Line, Circle, Ellipse, Parabola, Hyperbola, Analytic, Walking, Restricted,
}
