//! OCCT IntPatch_Line.hxx + subtypes (ALine / WLine / GLine / RLine)
//!
//! rcad: single flat struct with line_type discriminator.
//! OCCT has a class hierarchy:
//! - IntPatch_ALine: analytic line (Circle/Ellipse/Hyperbola/Parabola/Line)
//! - IntPatch_WLine: walking line (sequence of 3D points)
//! - IntPatch_GLine: geometric line (with tangent vectors)
//! - IntPatch_RLine: restriction line (boundary intersection)

use super::int_patch_type::IntPatchIType;

/// OCCT IntPatch_Line.hxx — intersection line base type.
#[derive(Debug, Clone)]
pub struct IntPatchLine {
    pub line_type: IntPatchIType,
    pub curve: rcad_kernel::geom::Curve3,
    pub t_range: [f64; 2],
    pub pcurve1: Option<rcad_kernel::geom::Curve2d>,
    pub pcurve2: Option<rcad_kernel::geom::Curve2d>,
    pub tolerance: f64,
    pub tang_tolerance: f64,
}
