//! OCCT-aligned: IntPatch_PrmPrmIntersection — parametric-parametric surface intersection.
//!
//! OCCT IntPatch_PrmPrmIntersection.hxx / .cxx (137K, ~3500 lines)
//!
//! rcad: delegates to marching infrastructure (face_face::intersect_faces / intss:).

use rcad_kernel::geom::Surface3;
use super::int_patch_line::IntPatchLine;
use super::int_patch_type::IntPatchIType;

pub struct PrmPrmIntersection {
    done: bool,
    empt: bool,
    spnt: Vec<super::int_patch_point::IntPatchPoint>,
    slin: Vec<IntPatchLine>,
}

impl PrmPrmIntersection {
    pub fn new() -> Self {
        Self { done: false, empt: true, spnt: Vec::new(), slin: Vec::new() }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn slin_ref(&self) -> &[IntPatchLine] { &self.slin }

    /// OCCT: Perform(Caro1,Domain1, Caro2,Domain2, TolTang,Epsilon,Deflection,Increment)
    /// rcad: delegates to face_face::intersect_faces (marching/numeric)
    pub fn perform(&mut self, s1: &Surface3, s2: &Surface3,
                   tol_arc: f64, tol_tang: f64, _fleche: f64, _uv_max_step: f64) {
        self.done = false; self.empt = true;
        self.slin.clear(); self.spnt.clear();

        let curves = crate::inttools::face_face::intersect_faces(s1, s2, tol_arc, tol_tang);
        self.slin = curves.into_iter().map(|c| IntPatchLine {
            line_type: IntPatchIType::Walking,
            curve: c.curve, t_range: c.t_range,
            pcurve1: c.pcurve1, pcurve2: c.pcurve2,
            tolerance: c.tolerance, tang_tolerance: c.tang_tolerance,
            wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: WLineType::Unknown,
        }).collect();
        self.empt = self.slin.is_empty();
        self.done = true;
    }
}
