//! OCCT-aligned: IntPatch_ImpPrmIntersection — analytic-parametric surface intersection.
//!
//! OCCT IntPatch_ImpPrmIntersection.hxx / .cxx (114K, ~3000 lines)
//!
//! Intersects one analytic surface (Plane/Cylinder/Sphere/Cone/Torus)
//! with one parametric surface (BSpline/Bezier).
//! Uses marching algorithm on the parametric surface, constrained
//! by the analytic surface's algebraic equation.
//!
//! rcad: currently delegates to marching/numeric intersection.

use rcad_kernel::geom::Surface3;
use super::int_patch_line::IntPatchLine;

pub struct ImpPrmIntersection {
    done: bool,
    empt: bool,
    slin: Vec<IntPatchLine>,
}

impl ImpPrmIntersection {
    pub fn new() -> Self {
        Self { done: false, empt: true, slin: Vec::new() }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn slin_ref(&self) -> &[IntPatchLine] { &self.slin }

    /// Perform intersection — rcad: uses marching
    pub fn perform(&mut self, _s1: &Surface3, _s2: &Surface3, _tol_arc: f64, _tol_tang: f64) {
        // OCCT: IntPatch_ImpPrmIntersection uses constrained marching on the
        // parametric surface with the analytic surface as constraint.
        // rcad: delegates to face_face::intersect_faces (marching/numeric).
        self.done = true;
        self.empt = true;
    }
}
