//! OCCT ChFi3d_Builder_6.cxx (continued) — the three ComputeData and three
//! SimulData overloads of ChFi3d_Builder (Stage 1f).
//!
//! Split from `chfi3d_builder_6` per the 2000-line guideline.  Line anchors
//! refer to OCCT ChFi3d_Builder_6.cxx.
//!
//! Architecture mappings (see chfi3d_builder_6.rs header for the full list):
//!   - `occ::handle<ChFiDS_ElSpine>` maps to `ChFiDSElSpineHandle`
//!     (Arc + RwLock — the walking mutates the shared guide in place).
//!   - `occ::handle<Adaptor3d_TopolTool>` maps to `&BRepTopolTool`
//!     (topalgo/brep_top_adaptor).
//!
//! Cross-file references (Stage 1e second batch, translated in parallel):
//!   BRepBlend_Walking             -> super::brep_blend_walking::BRepBlendWalking
//!   BRepBlend_SurfRstLineBuilder  -> super::brep_blend_surf_rst_line_builder::BRepBlendSurfRstLineBuilder
//!   BRepBlend_RstRstLineBuilder   -> super::brep_blend_rst_rst_line_builder::BRepBlendRstRstLineBuilder
//!   Blend_SurfPointFuncInv /
//!   Blend_SurfCurvFuncInv         -> super::brep_blend_func_inv (family file)
//!   ChFi3d_Builder::SearchFace    -> super::chfi3d_builder_2 (Builder_2.cxx
//!                                    L1523-1682, parallel agent)

use rcad_kernel::geom::Curve2d;
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::Shape;

use super::brep_blend_func_inv::{
    BlendCurvPointFuncInv, BlendFuncInv, BlendSurfCurvFuncInv, BlendSurfPointFuncInv,
};
use super::brep_blend_function::BlendFunction;
use super::brep_blend_line::BRepBlendLine;
use super::brep_blend_point::BlendPoint;
use super::brep_blend_rst_rst_function::BlendRstRstFunction;
use super::brep_blend_rst_rst_line_builder::BRepBlendRstRstLineBuilder;
use super::brep_blend_surf_rst_function::BlendSurfRstFunction;
use super::brep_blend_surf_rst_line_builder::BRepBlendSurfRstLineBuilder;
use super::brep_blend_walking::BRepBlendWalking;
use super::chfi3d::ChFi3dBuilder;
use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_6::{
    chfi3d_fil_common_point, comp_blend_point, is_obst, search_index, update_line,
};
use super::chfi_ds::{ChFiDSElSpine, ChFiDSSpineHandle};
use super::chfi3d_builder_2::BRepTopAdaptorTopolTool;

/// rcad mapping of `occ::handle<ChFiDS_ElSpine>` — a mutable shared guide
/// (SharedStripe precedent).
pub type ChFiDSElSpineHandle = std::sync::Arc<std::sync::RwLock<ChFiDSElSpine>>;

// =========================================================================
// OCCT ChFiDS_ElSpine curve-access methods consumed by ComputeData /
// SimulData (pending the ElSpine curve machinery — ChFiDS_ElSpine.hxx /
// ChFiDS_ElSpine.cxx; the rcad record carries no underlying curve yet).
// Kept next to the consumers per the Stage 1f delegation note.
// =========================================================================
impl ChFiDSElSpine {
    // OCCT ChFiDS_ElSpine::FirstParameter() / LastParameter() are the
    // first_parameter / last_parameter accessors in chfi_ds_spine.rs (1b).

    /// OCCT ChFiDS_ElSpine::FirstParameter(Val) — setter.
    pub fn set_first_parameter(&mut self, val: f64) {
        self.firstparam = val;
    }

    /// OCCT ChFiDS_ElSpine::LastParameter(Val) — setter.
    pub fn set_last_parameter(&mut self, val: f64) {
        self.lastparam = val;
    }

    /// OCCT ChFiDS_ElSpine::IsPeriodic().
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// OCCT ChFiDS_ElSpine::Period() — the parameter period.
    pub fn period(&self) -> f64 {
        self.lastparam - self.firstparam
    }

    /// OCCT ChFiDS_ElSpine::SetOrigin(Ori) — re-origin of the underlying
    /// curve (pending the curve machinery; boundary no-op).
    pub fn set_origin(&mut self, _ori: f64) {}

    /// OCCT ChFiDS_ElSpine::Resolution(Tol3d) — pending the curve
    /// machinery; the input 3d tolerance stands in.
    pub fn resolution(&self, tol3d: f64) -> f64 {
        tol3d
    }

    /// OCCT ChFiDS_ElSpine::D1(U, P, V) — pending the curve machinery
    /// (zero derivatives at the boundary).  gp_Pnt/gp_Vec map to DVec3.
    pub fn d1(&self, _u: f64) -> (glam::DVec3, glam::DVec3) {
        (glam::DVec3::ZERO, glam::DVec3::ZERO)
    }

    /// OCCT ChFiDS_ElSpine::D2(U, P, V1, V2) — pending the curve machinery.
    pub fn d2(&self, _u: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3) {
        (glam::DVec3::ZERO, glam::DVec3::ZERO, glam::DVec3::ZERO)
    }
}

/// OCCT `aHElSpine == HGuide` (handle identity) — rcad boundary substitute
/// on the (firstparam, lastparam) pair; see the offset-guide block in
/// ComputeData.
pub fn elspine_matches_handle(a: &ChFiDSElSpine, hguide: &ChFiDSElSpineHandle) -> bool {
    let hg = hguide.read().expect("elspine lock");
    a.firstparam == hg.firstparam && a.lastparam == hg.lastparam
}

/// The rcad Walking binds the guide as a &Curve3 (BRepBlend_Walking.cxx
/// ctor reads the ElSpine curve).  Pending boundary: the ElSpine curve
/// machinery is deferred (ChFiDSElSpine carries no curve yet), so the
/// ctor binds a straight-line placeholder at the axis origin.
pub fn elspine_guide_curve(hguide: &ChFiDSElSpineHandle) -> rcad_kernel::geom::Curve3 {
    let _ = hguide;
    rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(
        glam::DVec3::ZERO,
        glam::DVec3::X,
    ))
}

/// OCCT SearchFace(Spine, ...) tolerates a null spine (ThreeCorner); the
/// rcad ChFi3d_Builder::search_face (Builder_2) takes a handle — the null
/// case passes an empty Chamf handle (boundary note).
pub fn search_face_handle(spine: Option<&ChFiDSSpineHandle>) -> ChFiDSSpineHandle {
    match spine {
        Some(s) => s.clone(),
        None => ChFiDSSpineHandle::Chamf(super::chfi_ds::ChFiDSChamfSpine::new()),
    }
}

impl ChFi3dBuilder {
    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L1011-1170 — ComputeData (head of the path
    // edge/face for the bypass of obstacle).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn compute_data_surf_rst(
        &mut self,
        data: &mut super::chfi_ds::ChFiDSSurfData,
        hguide: &ChFiDSElSpineHandle,
        lin: &mut Option<BRepBlendLine>,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        pc2: &Curve2d,
        i2: &BRepTopAdaptorTopolTool,
        decroch: &mut bool,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        soldep: &Vector,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
    ) -> bool {
        // OCCT: BRepBlend_SurfRstLineBuilder TheWalk(S1, I1, S2, PC2, I2)
        let mut the_walk = BRepBlendSurfRstLineBuilder::new(s1, i1, s2, pc2, i2);

        data.set_first_extension_value(0.0);
        data.set_last_extension_value(0.0);

        let reverse = !forward || inside;
        let sp_first = hguide.read().expect("elspine lock").first_parameter();
        let sp_last = hguide.read().expect("elspine lock").last_parameter();
        let mut target = sp_last;
        if reverse {
            target = sp_first;
        }
        let targetsov = target;

        let mut ms = max_step;
        let mut again = 0;
        let nbptmin = 3; // jlr
        let mut nbpnt = 1;
        // the initial solution is reframed if necessary.
        let mut par_sol = Vector::new(1, 3);
        let mut new_first = p_first;
        if rec_p || rec_s || rec_rst {
            if !the_walk.perform_first_section(
                func,
                finv,
                finv_p,
                finv_c,
                p_first,
                target,
                &soldep.data.v,
                self.tolapp3d,
                self.tolapp2d,
                tol_guide,
                rec_rst,
                rec_p,
                rec_s,
                &mut new_first,
                &mut par_sol.data.v,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        while again < 2 {
            the_walk.perform(
                func,
                finv,
                finv_p,
                finv_c,
                new_first,
                *last,
                ms,
                self.tolapp3d,
                self.tolapp2d,
                tol_guide,
                &par_sol.data.v,
                fleche,
                appro,
            );

            if !the_walk.is_done() {
                return false;
            }

            if reverse {
                if !the_walk.complete(func, finv, finv_p, finv_c, sp_last) {
                    // OCCT: prints "Not completed" under OCCT_DEBUG only.
                }
            }

            *lin = Some(the_walk.line().clone());
            let line = lin.as_ref().expect("Lin");
            nbpnt = line.nb_points();
            if nbpnt <= 1 && again == 0 {
                again += 1;
                ms /= 50.0;
                target = targetsov;
            } else if nbpnt <= nbptmin && again == 0 {
                again += 1;
                let u1 = line.point(1).parameter();
                let u2 = line.point(nbpnt).parameter();
                ms = (u2 - u1) / (nbptmin as f64 + 1.0);
                target = targetsov;
            } else if nbpnt <= nbptmin {
                return false;
            } else {
                again = 2;
            }
        }
        if forward {
            *decroch = the_walk.decroch_end();
        } else {
            *decroch = the_walk.decroch_start();
        }
        let line = lin.as_ref().expect("Lin");
        *last = line.point(nbpnt).parameter();
        *first = line.point(1).parameter();
        true
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L1177-1341 — ComputeData (heading of the
    // path edge/edge for the bypass of obstacle).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn compute_data_rst_rst(
        &mut self,
        data: &mut super::chfi_ds::ChFiDSSurfData,
        hguide: &ChFiDSElSpineHandle,
        lin: &mut Option<BRepBlendLine>,
        s1: &BRepAdaptorSurface,
        pc1: &Curve2d,
        i1: &BRepTopAdaptorTopolTool,
        decroch1: &mut bool,
        s2: &BRepAdaptorSurface,
        pc2: &Curve2d,
        i2: &BRepTopAdaptorTopolTool,
        decroch2: &mut bool,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        soldep: &Vector,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p1: bool,
        rec_rst1: bool,
        rec_p2: bool,
        rec_rst2: bool,
    ) -> bool {
        // OCCT: BRepBlend_RstRstLineBuilder TheWalk(S1, PC1, I1, S2, PC2, I2)
        let mut the_walk = BRepBlendRstRstLineBuilder::new(s1, pc1, i1, s2, pc2, i2);

        data.set_first_extension_value(0.0);
        data.set_last_extension_value(0.0);

        let reverse = !forward || inside;
        let sp_first = hguide.read().expect("elspine lock").first_parameter();
        let sp_last = hguide.read().expect("elspine lock").last_parameter();
        let mut target = sp_last;
        if reverse {
            target = sp_first;
        }
        let targetsov = target;

        let mut ms = max_step;
        let mut again = 0;
        let nbptmin = 3; // jlr
        let mut nbpnt = 0;
        // the initial solution is reframed if necessary.
        let mut par_sol = Vector::new(1, 2);
        let mut new_first = p_first;
        if rec_p1 || rec_rst1 || rec_p2 || rec_rst2 {
            if !the_walk.perform_first_section(
                func,
                finv1,
                finv_p1,
                finv2,
                finv_p2,
                p_first,
                target,
                &soldep.data.v,
                self.tolapp3d,
                tol_guide,
                rec_rst1,
                rec_p1,
                rec_rst2,
                rec_p2,
                &mut new_first,
                &mut par_sol.data.v,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        while again < 2 {
            the_walk.perform(
                func,
                finv1,
                finv_p1,
                finv2,
                finv_p2,
                new_first,
                *last,
                ms,
                self.tolapp3d,
                tol_guide,
                &par_sol.data.v,
                fleche,
                appro,
            );

            if !the_walk.is_done() {
                return false;
            }

            if reverse {
                if !the_walk.complete(func, finv1, finv_p1, finv2, finv_p2, sp_last) {
                    // OCCT: prints "Not completed" under OCCT_DEBUG only.
                }
            }

            *lin = Some(the_walk.line().clone());
            let line = lin.as_ref().expect("Lin");
            nbpnt = line.nb_points();
            if nbpnt <= 1 && again == 0 {
                again += 1;
                ms /= 50.0;
                target = targetsov;
            } else if nbpnt <= nbptmin && again == 0 {
                again += 1;
                let u1 = line.point(1).parameter();
                let u2 = line.point(nbpnt).parameter();
                ms = (u2 - u1) / (nbptmin as f64 + 1.0);
                target = targetsov;
            } else if nbpnt <= nbptmin {
                return false;
            } else {
                again = 2;
            }
        }
        if forward {
            *decroch1 = the_walk.decroch1_end();
            *decroch2 = the_walk.decroch2_end();
        } else {
            *decroch1 = the_walk.decroch1_start();
            *decroch2 = the_walk.decroch2_start();
        }
        let line = lin.as_ref().expect("Lin");
        *last = line.point(nbpnt).parameter();
        *first = line.point(1).parameter();
        true
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L1348-1500 — SimulData (heading of the path
    // edge/face for the bypass of obstacle in simulation mode).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn simul_data_surf_rst(
        &mut self,
        _data: &mut super::chfi_ds::ChFiDSSurfData,
        hguide: &ChFiDSElSpineHandle,
        lin: &mut Option<BRepBlendLine>,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        pc2: &Curve2d,
        i2: &BRepTopAdaptorTopolTool,
        decroch: &mut bool,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        soldep: &Vector,
        nb_sec_min: i32,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
    ) -> bool {
        // OCCT: BRepBlend_SurfRstLineBuilder TheWalk(S1, I1, S2, PC2, I2)
        let mut the_walk = BRepBlendSurfRstLineBuilder::new(s1, i1, s2, pc2, i2);

        let reverse = !forward || inside;
        let sp_first = hguide.read().expect("elspine lock").first_parameter();
        let sp_last = hguide.read().expect("elspine lock").last_parameter();
        let mut target = sp_last;
        if reverse {
            target = sp_first;
        }
        let targetsov = target;

        let mut ms = max_step;
        let mut again = 0;
        let mut nbpnt = 0;
        // the starting solution is reframed if needed.
        let mut par_sol = Vector::new(1, 3);
        let mut new_first = p_first;
        if rec_p || rec_s || rec_rst {
            if !the_walk.perform_first_section(
                func,
                finv,
                finv_p,
                finv_c,
                p_first,
                target,
                &soldep.data.v,
                self.tolapp3d,
                self.tolapp2d,
                tol_guide,
                rec_rst,
                rec_p,
                rec_s,
                &mut new_first,
                &mut par_sol.data.v,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        while again < 2 {
            the_walk.perform(
                func,
                finv,
                finv_p,
                finv_c,
                new_first,
                *last,
                ms,
                self.tolapp3d,
                self.tolapp2d,
                tol_guide,
                &par_sol.data.v,
                fleche,
                appro,
            );
            if !the_walk.is_done() {
                return false;
            }
            if reverse {
                if !the_walk.complete(func, finv, finv_p, finv_c, sp_last) {
                    // OCCT: prints "Not completed" under OCCT_DEBUG only.
                }
            }
            *lin = Some(the_walk.line().clone());
            let line = lin.as_ref().expect("Lin");
            nbpnt = line.nb_points();
            if nbpnt <= 1 && again == 0 {
                again += 1;
                ms /= 50.0;
                target = targetsov;
            } else if nbpnt <= nb_sec_min && again == 0 {
                again += 1;
                let u1 = line.point(1).parameter();
                let u2 = line.point(nbpnt).parameter();
                ms = (u2 - u1) / (nb_sec_min as f64 + 1.0);
                target = targetsov;
            } else if nbpnt <= nb_sec_min {
                return false;
            } else {
                again = 2;
            }
        }
        if forward {
            *decroch = the_walk.decroch_end();
        } else {
            *decroch = the_walk.decroch_start();
        }
        let line = lin.as_ref().expect("Lin");
        *last = line.point(nbpnt).parameter();
        *first = line.point(1).parameter();
        true
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L1508-1668 — SimulData (heading of path
    // edge/edge for the bypass of obstacle in simulation mode).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn simul_data_rst_rst(
        &mut self,
        _data: &mut super::chfi_ds::ChFiDSSurfData,
        hguide: &ChFiDSElSpineHandle,
        lin: &mut Option<BRepBlendLine>,
        s1: &BRepAdaptorSurface,
        pc1: &Curve2d,
        i1: &BRepTopAdaptorTopolTool,
        decroch1: &mut bool,
        s2: &BRepAdaptorSurface,
        pc2: &Curve2d,
        i2: &BRepTopAdaptorTopolTool,
        decroch2: &mut bool,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        soldep: &Vector,
        nb_sec_min: i32,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p1: bool,
        rec_rst1: bool,
        rec_p2: bool,
        rec_rst2: bool,
    ) -> bool {
        // OCCT: BRepBlend_RstRstLineBuilder TheWalk(S1, PC1, I1, S2, PC2, I2)
        let mut the_walk = BRepBlendRstRstLineBuilder::new(s1, pc1, i1, s2, pc2, i2);

        let reverse = !forward || inside;
        let sp_first = hguide.read().expect("elspine lock").first_parameter();
        let sp_last = hguide.read().expect("elspine lock").last_parameter();
        let mut target = sp_last;
        if reverse {
            target = sp_first;
        }
        let targetsov = target;

        let mut ms = max_step;
        let mut again = 0;
        let mut nbpnt = 0;
        // The initial solution is reframed if necessary.
        let mut par_sol = Vector::new(1, 2);
        let mut new_first = p_first;
        if rec_p1 || rec_rst1 || rec_p2 || rec_rst2 {
            if !the_walk.perform_first_section(
                func,
                finv1,
                finv_p1,
                finv2,
                finv_p2,
                p_first,
                target,
                &soldep.data.v,
                self.tolapp3d,
                tol_guide,
                rec_rst1,
                rec_p1,
                rec_rst2,
                rec_p2,
                &mut new_first,
                &mut par_sol.data.v,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        while again < 2 {
            the_walk.perform(
                func,
                finv1,
                finv_p1,
                finv2,
                finv_p2,
                new_first,
                *last,
                ms,
                self.tolapp3d,
                tol_guide,
                &par_sol.data.v,
                fleche,
                appro,
            );
            if !the_walk.is_done() {
                return false;
            }
            if reverse {
                if !the_walk.complete(func, finv1, finv_p1, finv2, finv_p2, sp_last) {
                    // OCCT: prints "Not completed" under OCCT_DEBUG only.
                }
            }
            *lin = Some(the_walk.line().clone());
            let line = lin.as_ref().expect("Lin");
            nbpnt = line.nb_points();
            if nbpnt <= 1 && again == 0 {
                again += 1;
                ms /= 50.0;
                target = targetsov;
            } else if nbpnt <= nb_sec_min && again == 0 {
                again += 1;
                let u1 = line.point(1).parameter();
                let u2 = line.point(nbpnt).parameter();
                ms = (u2 - u1) / (nb_sec_min as f64 + 1.0);
                target = targetsov;
            } else if nbpnt <= nb_sec_min {
                return false;
            } else {
                again = 2;
            }
        }
        if forward {
            *decroch1 = the_walk.decroch1_end();
            *decroch2 = the_walk.decroch2_end();
        } else {
            *decroch1 = the_walk.decroch1_start();
            *decroch2 = the_walk.decroch2_start();
        }

        let line = lin.as_ref().expect("Lin");
        *last = line.point(nbpnt).parameter();
        *first = line.point(1).parameter();
        true
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L1676-2507 — ComputeData (construction of
    // elementary fillet by path).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn compute_data(
        &mut self,
        data: &mut super::chfi_ds::ChFiDSSurfData,
        hguide: &ChFiDSElSpineHandle,
        spine: Option<&ChFiDSSpineHandle>,
        lin: &mut Option<BRepBlendLine>,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        func: &mut dyn BlendFunction,
        finv: &mut dyn BlendFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tolguide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        soldep: &Vector,
        intf: &mut i32,
        intl: &mut i32,
        gd1: &mut bool,
        gd2: &mut bool,
        gf1: &mut bool,
        gf2: &mut bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
    ) -> bool {
        // Get offset guide if exists
        let mut offset_hguide: Option<ChFiDSElSpine> = None;
        if let Some(spine) = spine {
            if spine.base().my_mode == super::chfi_ds::ChFiDS_ChamfMode::ConstThroatWithPenetrationChamfer {
                // OCCT L1709-1719: find the offset elspine paired with HGuide
                // by handle identity (aHElSpine == HGuide).  rcad boundary:
                // ChFiDS_Spine stores plain ElSpine records while HGuide
                // arrives as a handle — identity is substituted by the
                // (firstparam, lastparam) pair until the spine list carries
                // handles (Stage 1f note).
                for (a, b) in spine.base().elspines.iter().zip(spine.base().offset_elspines.iter()) {
                    if elspine_matches_handle(a, hguide) {
                        offset_hguide = Some(b.clone());
                    }
                }
            }
        }

        // The extrensions are created in case of output of two domains
        // directly and not by path ( too hasardous ).
        data.set_first_extension_value(0.0);
        data.set_last_extension_value(0.0);

        // The eventual faces are restored to test the jump of edge.
        let f1 = s1.face.clone();
        let f2 = s2.face.clone();

        // Path framing variables
        let mut tol_guide = tolguide;
        let nbptmin = 4;

        // OCCT: BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide)
        let guide = elspine_guide_curve(hguide);
        let mut the_walk = BRepBlendWalking::new(&s1.surface, &s2.surface, i1, i2, &guide);

        // Start of removal, 2D path controls
        // that qui s'accomodent mal des surfaces a parametrages non homogenes
        // en u et en v are extinguished.
        the_walk.set_check2d(false);

        let mut ms = max_step;
        let mut nbpnt;
        let mut sp_first = hguide.read().expect("elspine lock").first_parameter();
        let mut sp_last = hguide.read().expect("elspine lock").last_parameter();

        // When the start point is inside, the path goes first to the left
        // to determine the Last for the periodicals.
        let reverse = !forward || inside;
        let mut target;
        if reverse {
            target = sp_first;
            if *intf == 0 {
                target = *last;
            }
        } else {
            target = sp_last + sp_last.abs();
            if *intl == 0 {
                target = *last;
            }
        }

        // In case if the singularity is pre-determined,
        // the path is indicated.
        if let Some(spine) = spine {
            if spine.base().is_tangency_extremity(true) {
                let v = spine.base().first_vertex();
                let e = spine.base().edges(1).clone();
                let param = spine.base().first_parameter();
                let mut bp = BlendPoint::new();
                let brep = self.my_brep.clone();
                if comp_blend_point(&brep, &v, &e, param, &f1, &f2, &mut bp) {
                    let mut vec = Vector::new(1, 4);
                    let (u1, v1) = bp.parameters_on_s1();
                    vec.set(1, u1);
                    vec.set(2, v1);
                    let (u2, v2) = bp.parameters_on_s2();
                    vec.set(3, u2);
                    vec.set(4, v2);
                    func.set_param(param);
                    if func.is_solution(&vec.data.v, self.tolapp3d) {
                        the_walk.add_singular_point(bp.clone());
                    }
                }
            }
            if spine.base().is_tangency_extremity(false) {
                let v = spine.base().last_vertex();
                let e = spine.base().edges(spine.base().nb_edges()).clone();
                let param = spine.base().last_parameter();
                let mut bp = BlendPoint::new();
                let brep = self.my_brep.clone();
                if comp_blend_point(&brep, &v, &e, param, &f1, &f2, &mut bp) {
                    let mut vec = Vector::new(1, 4);
                    let (u1, v1) = bp.parameters_on_s1();
                    vec.set(1, u1);
                    vec.set(2, v1);
                    let (u2, v2) = bp.parameters_on_s2();
                    vec.set(3, u2);
                    vec.set(4, v2);
                    func.set_param(param);
                    if func.is_solution(&vec.data.v, self.tolapp3d) {
                        the_walk.add_singular_point(bp.clone());
                    }
                }
            }
        }

        // The starting solution is reframed if necessary.
        //**********************************************//
        let mut par_sol = Vector::new(1, 4);
        let mut new_first = p_first;
        if rec_on_s1 || rec_on_s2 {
            if !the_walk.perform_first_section_recate(
                func,
                finv,
                p_first,
                target,
                &soldep.data.v,
                self.tolapp3d,
                tol_guide,
                rec_on_s1,
                rec_on_s2,
                &mut new_first,
                &mut par_sol.data.v,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        // First the valid part is calculate, without caring for the extensions.
        //******************************************************************//
        let mut again = 0;
        let mut tchernobyl = false;
        let mut u1sov = 0.0;
        let mut u2sov = 0.0;
        let mut bif = Shape::null();
        // Max step is relevant, but too great, the vector is required to detect
        // the twists.
        if ((*last - *first).abs() <= ms * 5.0)
            && ((*last - *first).abs() >= 0.01 * (new_first - target).abs())
        {
            ms = (*last - *first).abs() * 0.2;
        }

        while again < 3 {
            // Path.
            if again == 0 && ms < 5.0 * tol_guide {
                ms = 5.0 * tol_guide;
            } else {
                if 5.0 * tol_guide > ms {
                    tol_guide = ms / 5.0;
                }
            }
            the_walk.perform(
                func,
                finv,
                new_first,
                target,
                ms,
                self.tolapp3d,
                tol_guide,
                &par_sol.data.v,
                fleche,
                appro,
            );
            if !the_walk.is_done() {
                return false;
            }
            *lin = Some(the_walk.line().clone());
            if hguide.read().expect("elspine lock").is_periodic() && inside {
                {
                    let line = lin.as_ref().expect("Lin");
                    sp_first = line.point(1).parameter();
                }
                sp_last = sp_first + hguide.read().expect("elspine lock").period();
                {
                    // OCCT: HGuide->FirstParameter(SpFirst);
                    //       HGuide->LastParameter(SpLast);
                    //       HGuide->SetOrigin(SpFirst);
                    let mut hg = hguide.write().expect("elspine lock");
                    hg.set_first_parameter(sp_first);
                    hg.set_last_parameter(sp_last);
                    hg.set_origin(sp_first);
                }
                if let Some(off) = offset_hguide.as_mut() {
                    off.set_first_parameter(sp_first);
                    off.set_last_parameter(sp_last);
                    off.set_origin(sp_first);
                }
            }
            let mut complmnt = true;
            if inside {
                complmnt = the_walk.complete(func, finv, sp_last);
            }
            if !complmnt {
                return false;
            }

            // The result is controlled using two criterions :
            // - if there is enough points,
            // - if one has gone far enough.
            let line = lin.as_ref().expect("Lin");
            nbpnt = line.nb_points();
            if nbpnt == 0 {
                return false;
            }
            let fpointpar = line.point(1).parameter();
            let lpointpar = line.point(nbpnt).parameter();
            drop(line);

            let factor = 1.0 / (nbptmin as f64 + 1.0);
            let mut okdeb = forward && !inside;
            let mut okfin = !forward && !inside;
            if !okdeb {
                let narc1 = {
                    let line = lin.as_ref().expect("Lin");
                    line.start_point_on_first().nb_point_on_rst()
                };
                let narc2 = {
                    let line = lin.as_ref().expect("Lin");
                    line.start_point_on_second().nb_point_on_rst()
                };
                okdeb = narc1 > 0 || narc2 > 0 || (fpointpar - *first) < 10.0 * tol_guide;
            }
            if !okfin {
                let narc1 = {
                    let line = lin.as_ref().expect("Lin");
                    line.end_point_on_first().nb_point_on_rst()
                };
                let narc2 = {
                    let line = lin.as_ref().expect("Lin");
                    line.end_point_on_second().nb_point_on_rst()
                };
                okfin = narc1 > 0 || narc2 > 0 || (*last - lpointpar) < 10.0 * tol_guide;
            }
            if !okdeb || !okfin || nbpnt == 1 {
                // It drags, the controls are extended, it is  expected to
                // evaluate a satisfactory maximum step. If it already done,
                // quit.
                if tchernobyl {
                    return false;
                }
                tchernobyl = true;
                the_walk.set_check(false);
                if nbpnt == 1 {
                    ms *= 0.01;
                } else {
                    // EvalStep(Lin)
                    ms = (lpointpar - fpointpar) / nbpnt as f64;
                }
            } else if nbpnt < nbptmin {
                if again == 0 {
                    u1sov = fpointpar;
                    u2sov = lpointpar;
                    ms = (lpointpar - fpointpar) * factor;
                } else if again == 1 {
                    if (fpointpar - u1sov).abs() >= tol_guide
                        || (lpointpar - u2sov).abs() >= tol_guide
                    {
                        ms = (lpointpar - fpointpar) * factor;
                    } else {
                        return false;
                    }
                }
                again += 1;
            } else {
                again = 3;
            }
        }

        if the_walk.twist_on_s1() {
            data.set_twist_on_s1(true);
        }
        if the_walk.twist_on_s2() {
            data.set_twist_on_s2(true);
        }

        // Here there is a more or less presentable result
        // however it covers a the minimum zone.
        // The extensions are targeted.
        //*****************************//

        *gd1 = false;
        *gd2 = false;
        *gf1 = false;
        *gf2 = false;

        let unseulsuffitdeb = *intf >= 2;
        let unseulsuffitfin = *intl >= 2;
        let noproldeb = *intf >= 3;
        let noprolfin = *intl >= 3;

        let rab = 0.03 * (sp_last - sp_first);

        let mut debarc1 = false;
        let mut debarc2 = false;
        let mut debcas1 = false;
        let mut debcas2 = false;
        let mut debobst1 = false;
        let mut debobst2 = false;

        let mut finarc1 = false;
        let mut finarc2 = false;
        let mut fincas1 = false;
        let mut fincas2 = false;
        let mut finobst1 = false;
        let mut finobst2 = false;

        let mut narc1;
        let mut narc2;

        let mut backw_continue_failed = false; // eap
        if reverse && *intf != 0 {
            {
                let line = lin.as_ref().expect("Lin");
                narc1 = line.start_point_on_first().nb_point_on_rst();
                narc2 = line.start_point_on_second().nb_point_on_rst();
            }
            if narc1 != 0 {
                let brep = self.my_brep.clone();
                let line = lin.as_ref().expect("Lin");
                chfi3d_fil_common_point(
                    &brep,
                    line.start_point_on_first(),
                    line.transition_on_s1(),
                    true,
                    data.change_vertex_first_on_s1(),
                    self.tolapp3d,
                );
                drop(line);
                debarc1 = true;
                if !self.search_face(&mut search_face_handle(spine), &data.vertex_first_on_s1().clone(), &f1, &mut bif) {
                    // It is checked if there is not an obstacle.
                    debcas1 = true;
                    if let Some(spine) = spine {
                        if spine.base().is_periodic() {
                            debobst1 = true;
                        } else {
                            debobst1 = is_obst(
                                data.vertex_first_on_s1(),
                                &spine.base().first_vertex(),
                                &self.my_ve_map,
                            );
                        }
                    }
                }
            }
            if narc2 != 0 {
                let brep = self.my_brep.clone();
                let line = lin.as_ref().expect("Lin");
                chfi3d_fil_common_point(
                    &brep,
                    line.start_point_on_second(),
                    line.transition_on_s2(),
                    true,
                    data.change_vertex_first_on_s2(),
                    self.tolapp3d,
                );
                drop(line);
                debarc2 = true;
                if !self.search_face(&mut search_face_handle(spine), &data.vertex_first_on_s2().clone(), &f2, &mut bif) {
                    // It is checked if it is not an obstacle.
                    debcas2 = true;
                    if let Some(spine) = spine {
                        if spine.base().is_periodic() {
                            debobst2 = true;
                        } else {
                            debobst2 = is_obst(
                                data.vertex_first_on_s2(),
                                &spine.base().first_vertex(),
                                &self.my_ve_map,
                            );
                        }
                    }
                }
            }
            let mut oncontinue = !noproldeb && (narc1 != 0 || narc2 != 0);
            if debobst1 || debobst2 {
                oncontinue = false;
            } else if debcas1 && debcas2 {
                oncontinue = false;
            } else if (!debcas1 && debarc1) || (!debcas2 && debarc2) {
                oncontinue = false;
            }

            if oncontinue {
                the_walk.classification_on_s1(!debarc1);
                the_walk.classification_on_s2(!debarc2);
                the_walk.set_check2d(true); // It should be strict (PMN)
                the_walk.continu(func, finv, target);
                the_walk.classification_on_s1(true);
                the_walk.classification_on_s2(true);
                the_walk.set_check2d(false);
                {
                    let line = lin.as_ref().expect("Lin");
                    narc1 = line.start_point_on_first().nb_point_on_rst();
                    narc2 = line.start_point_on_second().nb_point_on_rst();
                }
                //  modified by eap Fri Feb  8 11:43:48 2002 ___BEGIN___
                if !debarc1 {
                    if narc1 == 0 {
                        let line = lin.as_ref().expect("Lin");
                        backw_continue_failed =
                            line.start_point_on_first().parameter_on_guide() > target;
                    } else {
                        let brep = self.my_brep.clone();
                        let line = lin.as_ref().expect("Lin");
                        chfi3d_fil_common_point(
                            &brep,
                            line.start_point_on_first(),
                            line.transition_on_s1(),
                            true,
                            data.change_vertex_first_on_s1(),
                            self.tolapp3d,
                        );
                        drop(line);
                        debarc1 = true;
                        if !self.search_face(&mut search_face_handle(spine), &data.vertex_first_on_s1().clone(), &f1, &mut bif)
                        {
                            // It is checked if it is not an obstacle.
                            debcas1 = true;
                        }
                    }
                }
                if !debarc2 {
                    if narc2 == 0 {
                        let line = lin.as_ref().expect("Lin");
                        backw_continue_failed =
                            line.start_point_on_second().parameter_on_guide() > target;
                    } else {
                        let brep = self.my_brep.clone();
                        let line = lin.as_ref().expect("Lin");
                        chfi3d_fil_common_point(
                            &brep,
                            line.start_point_on_second(),
                            line.transition_on_s2(),
                            true,
                            data.change_vertex_first_on_s2(),
                            self.tolapp3d,
                        );
                        drop(line);
                        debarc2 = true;
                        if !self.search_face(&mut search_face_handle(spine), &data.vertex_first_on_s2().clone(), &f2, &mut bif)
                        {
                            // It is checked if it is not an obstacle.
                            debcas2 = true;
                        }
                    }
                }
                if backw_continue_failed {
                    // if we leave backwContinueFailed as is, we will stop in
                    // this direction but we are to continue if there are no
                    // more faces on the side with arc check this condition
                    let a_cp = if debarc1 {
                        data.vertex_first_on_s1().clone()
                    } else {
                        data.vertex_first_on_s2().clone()
                    };
                    if a_cp.is_on_arc() && bif.is_null() {
                        backw_continue_failed = false;
                    }
                }
            }
        }
        let mut forw_continue_failed = false;
        //  modified by eap Fri Feb  8 11:44:11 2002 ___END___
        if forward && *intl != 0 {
            target = sp_last;
            {
                let line = lin.as_ref().expect("Lin");
                narc1 = line.end_point_on_first().nb_point_on_rst();
                narc2 = line.end_point_on_second().nb_point_on_rst();
            }
            if narc1 != 0 {
                let brep = self.my_brep.clone();
                let line = lin.as_ref().expect("Lin");
                chfi3d_fil_common_point(
                    &brep,
                    line.end_point_on_first(),
                    line.transition_on_s1(),
                    false,
                    data.change_vertex_last_on_s1(),
                    self.tolapp3d,
                );
                drop(line);
                finarc1 = true;
                if !self.search_face(&mut search_face_handle(spine), &data.vertex_last_on_s1().clone(), &f1, &mut bif) {
                    // It is checked if it is not an obstacle.
                    fincas1 = true;
                    if let Some(spine) = spine {
                        finobst1 =
                            is_obst(data.vertex_last_on_s1(), &spine.base().last_vertex(), &self.my_ve_map);
                    }
                }
            }
            if narc2 != 0 {
                let brep = self.my_brep.clone();
                let line = lin.as_ref().expect("Lin");
                chfi3d_fil_common_point(
                    &brep,
                    line.end_point_on_second(),
                    line.transition_on_s2(),
                    false,
                    data.change_vertex_last_on_s2(),
                    self.tolapp3d,
                );
                drop(line);
                finarc2 = true;
                if !self.search_face(&mut search_face_handle(spine), &data.vertex_last_on_s2().clone(), &f2, &mut bif) {
                    // It is checked if it is not an obstacle.
                    fincas2 = true;
                    if let Some(spine) = spine {
                        finobst2 =
                            is_obst(data.vertex_last_on_s2(), &spine.base().last_vertex(), &self.my_ve_map);
                    }
                }
            }
            let mut oncontinue = !noprolfin && (narc1 != 0 || narc2 != 0);
            if finobst1 || finobst2 {
                oncontinue = false;
            } else if fincas1 && fincas2 {
                oncontinue = false;
            } else if (!fincas1 && finarc1) || (!fincas2 && finarc2) {
                oncontinue = false;
            }

            if oncontinue {
                the_walk.classification_on_s1(!finarc1);
                the_walk.classification_on_s2(!finarc2);
                the_walk.set_check2d(true); // It should be strict (PMN)
                the_walk.continu(func, finv, target);
                the_walk.classification_on_s1(true);
                the_walk.classification_on_s2(true);
                the_walk.set_check2d(false);
                {
                    let line = lin.as_ref().expect("Lin");
                    narc1 = line.end_point_on_first().nb_point_on_rst();
                    narc2 = line.end_point_on_second().nb_point_on_rst();
                }
                //  modified by eap Fri Feb  8 11:44:57 2002 ___BEGIN___
                if !finarc1 {
                    if narc1 == 0 {
                        let line = lin.as_ref().expect("Lin");
                        forw_continue_failed = line.end_point_on_first().parameter_on_guide() < target;
                    } else {
                        let brep = self.my_brep.clone();
                        let line = lin.as_ref().expect("Lin");
                        chfi3d_fil_common_point(
                            &brep,
                            line.end_point_on_first(),
                            line.transition_on_s1(),
                            false,
                            data.change_vertex_last_on_s1(),
                            self.tolapp3d,
                        );
                        drop(line);
                        finarc1 = true;
                        if !self.search_face(&mut search_face_handle(spine), &data.vertex_last_on_s1().clone(), &f1, &mut bif) {
                            // It is checked if it is not an obstacle.
                            fincas1 = true;
                        }
                    }
                }
                if !finarc2 {
                    if narc2 == 0 {
                        let line = lin.as_ref().expect("Lin");
                        forw_continue_failed = line.end_point_on_second().parameter_on_guide() < target;
                    } else {
                        let brep = self.my_brep.clone();
                        let line = lin.as_ref().expect("Lin");
                        chfi3d_fil_common_point(
                            &brep,
                            line.end_point_on_second(),
                            line.transition_on_s2(),
                            false,
                            data.change_vertex_last_on_s2(),
                            self.tolapp3d,
                        );
                        drop(line);
                        finarc2 = true;
                        if !self.search_face(&mut search_face_handle(spine), &data.vertex_last_on_s2().clone(), &f2, &mut bif) {
                            // On regarde si ce n'est pas un obstacle.
                            fincas2 = true;
                        }
                    }
                }
                if forw_continue_failed {
                    // if we leave forwContinueFailed as is, we will stop in
                    // this direction but we are to continue if there are no
                    // more faces on the side with arc check this condition
                    let a_cp = if finarc1 {
                        data.vertex_last_on_s1().clone()
                    } else {
                        data.vertex_last_on_s2().clone()
                    };
                    if a_cp.is_on_arc() && bif.is_null() {
                        forw_continue_failed = false;
                    }
                }
                //  modified by eap Fri Feb  8 11:45:10 2002 ___END___
            }
        }
        {
            let line = lin.as_ref().expect("Lin");
            nbpnt = line.nb_points();
            *first = line.point(1).parameter();
            *last = line.point(nbpnt).parameter();
        }

        // ============= INVALIDATION EVENTUELLE =============
        // ------ Preparation des prolongement par plan tangent -----
        if reverse && *intf != 0 {
            *gd1 = debcas1; // skv(occ67)
            *gd2 = debcas2; // skv(occ67)
            if (debarc1 ^ debarc2) && !unseulsuffitdeb && (*first != sp_first) {
                // Case of incomplete path, of course this ends badly :
                // the result is truncated instead of exit.
                let sortie;
                let mut ind;
                if debarc1 {
                    sortie = data.vertex_first_on_s1().parameter();
                } else {
                    sortie = data.vertex_first_on_s2().parameter();
                }
                if sortie - *first > self.tolesp {
                    {
                        let line = lin.as_ref().expect("Lin");
                        ind = search_index(sortie, &line);
                        if line.point(ind).parameter() == sortie {
                            ind -= 1;
                        }
                    }
                    if ind >= 1 {
                        let mut line = lin.as_mut().expect("Lin");
                        line.remove(1, ind);
                        update_line(&mut line, true);
                    }
                    let line = lin.as_ref().expect("Lin");
                    nbpnt = line.nb_points();
                    *first = line.point(1).parameter();
                }
            } else if (*intf >= 5) && !debarc1 && !debarc2 && (*first != sp_first) {
                let sortie = (2.0 * *first + *last) / 3.0;
                let mut ind;
                if sortie - *first > self.tolesp {
                    {
                        let line = lin.as_ref().expect("Lin");
                        ind = search_index(sortie, &line);
                        if line.point(ind).parameter() == sortie {
                            ind -= 1;
                        }
                    }
                    if nbpnt - ind < 3 {
                        ind = nbpnt - 3;
                    }
                    if ind >= 1 {
                        let mut line = lin.as_mut().expect("Lin");
                        line.remove(1, ind);
                        update_line(&mut line, true);
                    }
                    let line = lin.as_ref().expect("Lin");
                    nbpnt = line.nb_points();
                    *first = line.point(1).parameter();
                }
            }
            if *gd1 && *gd2 {
                let line = lin.as_ref().expect("Lin");
                target = (line.point(1).parameter() - rab).min(*first);
                target = target.max(sp_first);
                data.set_first_extension_value((line.point(1).parameter() - target).abs());
            }
            if *intf != 0 && !unseulsuffitdeb {
                // eap
                *intf = ((*gd1 && *gd2) || backw_continue_failed) as i32;
            } else if *intf != 0 && unseulsuffitdeb && (*intf < 5) {
                *intf = (*gd1 || *gd2) as i32;
                // It is checked if there is no new face.
                if *intf != 0 && ((!debcas1 && debarc1) || (!debcas2 && debarc2)) {
                    *intf = 0;
                }
            } else if *intf < 5 {
                *intf = 0;
            }
        }

        if forward && *intl != 0 {
            *gf1 = fincas1; // skv(occ67)
            *gf2 = fincas2; // skv(occ67)
            if (finarc1 ^ finarc2) && !unseulsuffitfin && (*last != sp_last) {
                // Case of incomplete path, of course, this ends badly :
                // the result is truncated instead of exit.
                let sortie;
                let mut ind;
                if finarc1 {
                    sortie = data.vertex_last_on_s1().parameter();
                } else {
                    sortie = data.vertex_last_on_s2().parameter();
                }
                if *last - sortie > self.tolesp {
                    {
                        let line = lin.as_ref().expect("Lin");
                        ind = search_index(sortie, &line);
                        if line.point(ind).parameter() == sortie {
                            ind += 1;
                        }
                    }
                    if ind <= nbpnt {
                        let mut line = lin.as_mut().expect("Lin");
                        line.remove(ind, nbpnt);
                        update_line(&mut line, false);
                    }
                    let line = lin.as_ref().expect("Lin");
                    nbpnt = line.nb_points();
                    *last = line.point(nbpnt).parameter();
                }
            } else if (*intl >= 5) && !finarc1 && !finarc2 && (*last != sp_last) {
                // The same in case when the entire "Lin" is an extension
                let sortie = (*first + 2.0 * *last) / 3.0;
                let mut ind;
                if *last - sortie > self.tolesp {
                    {
                        let line = lin.as_ref().expect("Lin");
                        ind = search_index(sortie, &line);
                        if line.point(ind).parameter() == sortie {
                            ind += 1;
                        }
                    }
                    if ind < 3 {
                        ind = 3;
                    }
                    if ind <= nbpnt {
                        let mut line = lin.as_mut().expect("Lin");
                        line.remove(ind, nbpnt);
                        update_line(&mut line, false);
                    }
                    let line = lin.as_ref().expect("Lin");
                    nbpnt = line.nb_points();
                    *last = line.point(nbpnt).parameter();
                }
            }
            if *gf1 && *gf2 {
                let line = lin.as_ref().expect("Lin");
                target = (line.point(nbpnt).parameter() + rab).max(*last);
                target = target.min(sp_last);
                data.set_last_extension_value((target - line.point(nbpnt).parameter()).abs());
            }

            if *intl != 0 && !unseulsuffitfin {
                // eap
                *intl = ((*gf1 && *gf2) || forw_continue_failed) as i32;
            } else if *intl != 0 && unseulsuffitfin && (*intl < 5) {
                *intl = (*gf1 || *gf2) as i32; // It is checked if there is no new face.
                if *intl != 0 && ((!fincas1 && finarc1) || (!fincas2 && finarc2)) {
                    *intl = 0;
                }
            } else if *intl < 5 {
                *intl = 0;
            }
        }
        true
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L2511-2703 — SimulData (elementary fillet
    // simulation by path with an additional guide).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn simul_data_walking(
        &mut self,
        _data: &mut super::chfi_ds::ChFiDSSurfData,
        hguide: &ChFiDSElSpineHandle,
        additional_hguide: Option<&ChFiDSElSpineHandle>,
        lin: &mut Option<BRepBlendLine>,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        func: &mut dyn BlendFunction,
        finv: &mut dyn BlendFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tolguide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        soldep: &Vector,
        nb_sec_min: i32,
        rec_on_s1: bool,
        rec_on_s2: bool,
    ) -> bool {
        // OCCT: BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide)
        let guide = elspine_guide_curve(hguide);
        let mut the_walk = BRepBlendWalking::new(&s1.surface, &s2.surface, i1, i2, &guide);
        the_walk.set_check2d(false);

        let mut ms = max_step;
        let mut tol_guide = tolguide;
        let mut nbpnt = 0;
        let mut sp_first = hguide.read().expect("elspine lock").first_parameter();
        let mut sp_last = hguide.read().expect("elspine lock").last_parameter();
        let reverse = !forward || inside;
        let mut target;
        if reverse {
            target = sp_first;
        } else {
            target = sp_last;
        }

        let targetsov = target;
        let mut u1sov = 0.0;
        let mut u2sov = 0.0;
        // on recadre la solution de depart a la demande.
        let mut par_sol = Vector::new(1, 4);
        let mut new_first = p_first;
        if rec_on_s1 || rec_on_s2 {
            if !the_walk.perform_first_section_recate(
                func,
                finv,
                p_first,
                target,
                &soldep.data.v,
                self.tolapp3d,
                tol_guide,
                rec_on_s1,
                rec_on_s2,
                &mut new_first,
                &mut par_sol.data.v,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }
        let mut again = 0;
        while again < 3 {
            // When the start point is inside, the path goes first to the left
            // to determine the Last for the periodicals.
            if again == 0 && ms < 5.0 * tol_guide {
                ms = 5.0 * tol_guide;
            } else {
                if 5.0 * tol_guide > ms {
                    tol_guide = ms / 5.0;
                }
            }

            the_walk.perform(
                func,
                finv,
                new_first,
                target,
                ms,
                self.tolapp3d,
                tol_guide,
                &par_sol.data.v,
                fleche,
                appro,
            );

            if !the_walk.is_done() {
                return false;
            }
            *lin = Some(the_walk.line().clone());
            if reverse {
                if hguide.read().expect("elspine lock").is_periodic() {
                    {
                        let line = lin.as_ref().expect("Lin");
                        sp_first = line.point(1).parameter();
                    }
                    sp_last = sp_first + hguide.read().expect("elspine lock").period();
                    {
                        // OCCT: HGuide->FirstParameter(SpFirst);
                        //       HGuide->LastParameter(SpLast);
                        let mut hg = hguide.write().expect("elspine lock");
                        hg.set_first_parameter(sp_first);
                        hg.set_last_parameter(sp_last);
                    }
                    if let Some(add) = additional_hguide {
                        // OCCT: AdditionalHGuide->FirstParameter(SpFirst);
                        //       AdditionalHGuide->LastParameter(SpLast);
                        let mut ah = add.write().expect("elspine lock");
                        ah.set_first_parameter(sp_first);
                        ah.set_last_parameter(sp_last);
                    }
                }
                let mut complmnt = true;
                if inside {
                    complmnt = the_walk.complete(func, finv, sp_last);
                }
                if !complmnt {
                    return false;
                }
            }
            {
                let line = lin.as_ref().expect("Lin");
                nbpnt = line.nb_points();
                let factor = 1.0 / (nb_sec_min as f64 + 1.0);
                if nbpnt == 0 {
                    return false;
                } else if nbpnt == 1 && again == 0 {
                    again += 1;
                    ms *= 0.01;
                    target = targetsov;
                    let u = line.point(1).parameter();
                    u1sov = u;
                    u2sov = u;
                } else if nbpnt < nb_sec_min && again == 0 {
                    again += 1;
                    let u1 = line.point(1).parameter();
                    let u2 = line.point(nbpnt).parameter();
                    u1sov = u1;
                    u2sov = u2;
                    ms = (u2 - u1) * factor;
                    target = targetsov;
                } else if nbpnt < nb_sec_min && again == 1 {
                    let u1 = line.point(1).parameter();
                    let u2 = line.point(nbpnt).parameter();
                    if (u1 - u1sov).abs() >= tol_guide || (u2 - u2sov).abs() >= tol_guide {
                        again += 1;
                        ms /= 100.0;
                        target = targetsov;
                    } else {
                        return false;
                    }
                } else if nbpnt < nb_sec_min {
                    return false;
                } else {
                    again = 3;
                }
            }
        }
        let line = lin.as_ref().expect("Lin");
        *first = line.point(1).parameter();
        *last = line.point(nbpnt).parameter();
        true
    }
}
