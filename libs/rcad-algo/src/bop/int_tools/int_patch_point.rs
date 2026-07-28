//! IntPatch_Point — intersection point between two surfaces.
use glam::DVec3;

#[derive(Debug, Clone)]
pub struct IntPatchPoint {
    pub pnt: DVec3,
    pub param_on_first: (f64, f64),
    pub param_on_second: (f64, f64),
    pub tolerance: f64,
}
impl IntPatchPoint {
    pub fn new() -> Self {
        IntPatchPoint {
            pnt: DVec3::ZERO,
            param_on_first: (0.0, 0.0),
            param_on_second: (0.0, 0.0),
            tolerance: 0.0,
        }
    }
}
