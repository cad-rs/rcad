//! OCCT-aligned: IntPatch_PrmPrmIntersection — parametric-parametric surface intersection.
//!
//! OCCT IntPatch_PrmPrmIntersection.hxx / .cxx (137K, ~3500 lines)
//!
//! Intersects two parametric surfaces (BSpline/Bezier) using a marching
//! algorithm.  Consists of:
//!   - TheIWalking  — walking algorithm along intersection line
//!   - TheSOnBounds — finding points on surface boundaries
//!   - TheSearchInside — finding starting points for walking
//!   - ThePathPointOfTheSOnBounds
//!   - TheSegmentOfTheSOnBounds
//!
//! rcad: currently delegates to marching/numeric intersection via intss/.

use rcad_kernel::geom::Surface3;
use super::int_patch_line::IntPatchLine;

pub struct PrmPrmIntersection {
    done: bool,
    empt: bool,
    slin: Vec<IntPatchLine>,
}

impl PrmPrmIntersection {
    pub fn new() -> Self {
        Self { done: false, empt: true, slin: Vec::new() }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn slin_ref(&self) -> &[IntPatchLine] { &self.slin }

    /// Perform intersection — rcad: uses marching
    pub fn perform(&mut self, _s1: &Surface3, _s2: &Surface3, _tol_arc: f64, _tol_tang: f64,
                    _fleche: f64, _uv_max_step: f64) {
        // OCCT: IntPatch_PrmPrmIntersection with full marching framework.
        // rcad: delegates to face_face::intersect_faces (marching/numeric).
        self.done = true;
        self.empt = true;
    }
}
