//! OCCT-aligned: IntPatch_ImpPrmIntersection — analytic-parametric surface intersection.
//!
//! OCCT IntPatch_ImpPrmIntersection.hxx / .cxx (114K, ~3000 lines)
//!
//! rcad: delegates to marching. OCCT algorithm structure documented below.

use rcad_kernel::geom::Surface3;
use super::int_patch_line::IntPatchLine;
use super::int_patch_line::WLineType;
use super::int_patch_type::IntPatchIType;
use super::int_surf_quadric::Quadric;

pub struct ImpPrmIntersection {
    done: bool, empt: bool,
    spnt: Vec<super::int_patch_point::IntPatchPoint>,
    slin: Vec<IntPatchLine>,
    my_is_start_pnt: bool, my_u_start: f64, my_v_start: f64,
}

impl ImpPrmIntersection {
    pub fn new() -> Self {
        Self { done: false, empt: true, spnt: Vec::new(), slin: Vec::new(),
            my_is_start_pnt: false, my_u_start: 0.0, my_v_start: 0.0 }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn slin_ref(&self) -> &[IntPatchLine] { &self.slin }

    pub fn set_start_point(&mut self, u: f64, v: f64) {
        self.my_is_start_pnt = true; self.my_u_start = u; self.my_v_start = v;
    }

    /// Perform intersection — rcad: marching
    /// OCCT L617: Set Quadric, Func, SOnBounds, SearchInside, IWalking, Decompose
    pub fn perform(&mut self, s1: &Surface3, s2: &Surface3,
                   tol_arc: f64, tol_tang: f64, _fleche: f64, _pas: f64) {
        self.done = false; self.empt = true; self.slin.clear(); self.spnt.clear();
        let (qsurf, psurf) = match (Quadric::from_surface3(s1), Quadric::from_surface3(s2)) {
            (Some(_), _) => (s1, s2),
            (None, Some(_)) => (s2, s1),
            (None, None) => { self.done = true; return; }
        };
        let curves = crate::inttools::face_face::intersect_faces(qsurf, psurf, tol_arc, tol_tang);
        self.slin = curves.into_iter().map(|c| IntPatchLine {
            line_type: IntPatchIType::Walking, curve: c.curve, t_range: c.t_range,
            pcurve1: c.pcurve1, pcurve2: c.pcurve2,
            tolerance: c.tolerance, tang_tolerance: c.tang_tolerance,
            wline_pnts: Vec::new(), is_purging_allowed: false, wl_type: WLineType::Unknown,
        }).collect();
        self.empt = self.slin.is_empty(); self.done = true;
    }
}
