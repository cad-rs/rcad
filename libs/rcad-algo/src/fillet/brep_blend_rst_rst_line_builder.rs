//! OCCT BRepBlend_RstRstLineBuilder (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_RstRstLineBuilder.hxx (L65-208) + BRepBlend_RstRstLineBuilder.cxx
//! (L151-1964) + BRepBlend_RstRstLineBuilder.lxx (L24-77).
//!
//! The class builds a BRepBlend_Line between two pcurves (restrictions) from
//! an approached starting solution (hxx L57-64).  The OCCT_DEBUG trace
//! statics (BBPP / tracederiv / Drawsect, cxx L39-147) are debug-only and
//! are not translated.
//!
//! Architecture mappings and pending boundaries (Stage 1e third batch): see
//! brep_blend_surf_rst_line_builder.rs — [`RstArc`] / [`DomainVertex`] /
//! [`DomainTool`] and the [`blend_tool_*`] functions are the shared bridge
//! from that module.  This builder's domains are single-restriction domains
//! (`domain1->Initialize(rst1)` / `domain2->Initialize(rst2)`), so the
//! vertex scans see no vertices until the consumer carries the edge
//! identity (pending boundary).
//!   - PENDING: Blend_RstRstFunction::Decroch(Sol, Tgrst1, Nrrst1, Tgrst2,
//!     Nrrst2) (CheckInside, cxx L1961) is missing from the rcad
//!     BlendRstRstFunction trait (shared trait file outside this batch's
//!     ownership; gap reported) — the stand-in reports Blend_NoDecroch.

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetRoot;

use crate::geomalgo::int_patch::transitions::{
    make_transition, Transition as IntSurfTransition,
};

use super::brep_blend::{BlendDecrochStatus, BlendStatus};
use super::brep_blend_extremity::BRepBlendExtremity;
use super::brep_blend_func_inv::{BlendCurvPointFuncInv, BlendSurfCurvFuncInv};
use super::brep_blend_line::{BRepBlendLine, IntSurfTypeTrans};
use super::brep_blend_point::BlendPoint;
use super::brep_blend_rst_rst_function::BlendRstRstFunction;
use super::brep_blend_surf_rst_line_builder_b::{
    blend_tool_parameter, blend_tool_tolerance, DomainTool, DomainVertex, RstArc,
};
use super::chfi3d::topabs_reverse;
use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_2::{BRepTopAdaptorTopolTool, TopAbsState};

/// OCCT ConvOrToTra(O) (cxx L1764-1771) — the orientation to transition
/// conversion.
fn conv_or_to_tra(
    o: rcad_kernel::topods::Orientation,
) -> IntSurfTypeTrans {
    if o == rcad_kernel::topods::Orientation::Forward {
        IntSurfTypeTrans::In
    } else {
        IntSurfTypeTrans::Out
    }
}

/// OCCT BRepBlend_RstRstLineBuilder — this class processes the data
/// resulting from Blend_CSWalking but it takes in consideration the Surface
/// supporting the curve to detect the breakpoint (hxx L36).  OCCT
/// inheritance: standalone class (no base).
pub struct BRepBlendRstRstLineBuilder<'a> {
    done: bool,
    line: BRepBlendLine,
    sol: Vec<f64>,
    surf1: &'a BRepAdaptorSurface,
    domain1: &'a BRepTopAdaptorTopolTool,
    surf2: &'a BRepAdaptorSurface,
    domain2: &'a BRepTopAdaptorTopolTool,
    rst1: RstArc,
    rst2: RstArc,
    tolpoint3d: f64,
    tolgui: f64,
    pasmax: f64,
    fleche: f64,
    param: f64,
    previous_p: BlendPoint,
    rebrou: bool,
    iscomplete: bool,
    comptra: bool,
    sens: f64,
    decrochdeb: BlendDecrochStatus,
    decrochfin: BlendDecrochStatus,
}

impl<'a> BRepBlendRstRstLineBuilder<'a> {
    /// OCCT BRepBlend_RstRstLineBuilder(Surf1, Rst1, Domain1, Surf2, Rst2,
    /// Domain2) (cxx L151-178).
    pub fn new(
        surf1: &'a BRepAdaptorSurface,
        rst1: &'a Curve2d,
        domain1: &'a BRepTopAdaptorTopolTool,
        surf2: &'a BRepAdaptorSurface,
        rst2: &'a Curve2d,
        domain2: &'a BRepTopAdaptorTopolTool,
    ) -> Self {
        BRepBlendRstRstLineBuilder {
            done: false,
            line: BRepBlendLine::new(),
            sol: vec![0.0; 2], // OCCT: math_Vector sol(1, 2)
            surf1,
            domain1,
            surf2,
            domain2,
            rst1: RstArc::bare(rst1),
            rst2: RstArc::bare(rst2),
            tolpoint3d: 0.0,
            tolgui: 0.0,
            pasmax: 0.0,
            fleche: 0.0,
            param: 0.0,
            previous_p: BlendPoint::new(),
            rebrou: false,
            iscomplete: false,
            comptra: false,
            sens: 0.0,
            decrochdeb: BlendDecrochStatus::NoDecroch,
            decrochfin: BlendDecrochStatus::NoDecroch,
        }
    }

    /// OCCT IsDone() (lxx L27-31).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Line() (lxx L37-44) — throws StdFail_NotDone when not done.
    pub fn line(&self) -> &BRepBlendLine {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_RstRstLineBuilder::Line");
        }
        &self.line
    }

    /// OCCT Decroch1Start() (lxx L50-53).
    pub fn decroch1_start(&self) -> bool {
        self.decrochdeb == BlendDecrochStatus::DecrochRst1
            || self.decrochdeb == BlendDecrochStatus::DecrochBoth
    }

    /// OCCT Decroch1End() (lxx L59-62).
    pub fn decroch1_end(&self) -> bool {
        self.decrochfin == BlendDecrochStatus::DecrochRst1
            || self.decrochfin == BlendDecrochStatus::DecrochBoth
    }

    /// OCCT Decroch2Start() (lxx L68-71).
    pub fn decroch2_start(&self) -> bool {
        self.decrochdeb == BlendDecrochStatus::DecrochRst2
            || self.decrochdeb == BlendDecrochStatus::DecrochBoth
    }

    /// OCCT Decroch2End() (lxx L77-80).
    pub fn decroch2_end(&self) -> bool {
        self.decrochfin == BlendDecrochStatus::DecrochRst2
            || self.decrochfin == BlendDecrochStatus::DecrochBoth
    }

    /// OCCT Perform(Func, Finv1, FinvP1, Finv2, FinvP2, Pdep, Pmax, MaxStep,
    /// Tol3d, TolGuide, ParDep, Fleche, Appro) (cxx L182-281).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        pdep: f64,
        pmax: f64,
        max_step: f64,
        tol3d: f64,
        tolguide: f64,
        par_dep: &[f64],
        fleche: f64,
        appro: bool,
    ) {
        self.done = false;
        self.iscomplete = false;
        self.comptra = false;
        self.line = BRepBlendLine::new();
        self.tolpoint3d = tol3d;
        self.tolgui = tolguide.abs();
        self.fleche = fleche.abs();
        self.rebrou = false;
        self.pasmax = max_step.abs();

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        self.param = pdep;
        func.set_param(self.param);

        if appro {
            let mut siturst1 = TopAbsState::Unknown;
            let mut siturst2 = TopAbsState::Unknown;
            let mut decroch = BlendDecrochStatus::NoDecroch;
            let mut tolerance = vec![0.0; 2];
            let mut infbound = vec![0.0; 2];
            let mut supbound = vec![0.0; 2];
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);
            let mut rsnld = FunctionSetRoot::new(func, &tolerance, 30);

            rsnld.perform(func, par_dep, &infbound, &supbound, false);

            if !rsnld.is_done() {
                return;
            }
            self.sol = rsnld.root();
            if !self.check_inside(func, &mut siturst1, &mut siturst2, &mut decroch) {
                return;
            }
        } else {
            self.sol = par_dep.to_vec();
        }

        let state = self.test_arret(func, false, BlendStatus::Ok);
        if state != BlendStatus::Ok {
            return;
        }
        // Update the line.
        self.line.append(self.previous_p.clone());
        let u = self.previous_p.parameter_on_c1();
        let v = self.previous_p.parameter_on_c2();
        let mut ptf1 = BRepBlendExtremity::new_on_curve(
            self.previous_p.point_on_c1(),
            u,
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        let mut ptf2 = BRepBlendExtremity::new_on_curve(
            self.previous_p.point_on_c2(),
            v,
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        if !self.previous_p.is_tangency_point() {
            ptf1.set_tangent(self.previous_p.tangent_on_c1());
            ptf2.set_tangent(self.previous_p.tangent_on_c2());
        }

        if self.sens > 0.0 {
            self.line.set_start_points(&ptf1, &ptf2);
        } else {
            self.line.set_end_points(&ptf1, &ptf2);
        }

        self.internal_perform(func, finv1, finv_p1, finv2, finv_p2, pmax);
        self.done = true;
    }

    /// OCCT PerformFirstSection(Func, Finv1, FinvP1, Finv2, FinvP2, Pdep,
    /// Pmax, ParDep, Tol3d, TolGuide, RecRst1, RecP1, RecRst2, RecP2, Psol,
    /// ParSol) (cxx L285-534).
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn perform_first_section(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        pdep: f64,
        pmax: f64,
        par_dep: &[f64],
        tol3d: f64,
        tolguide: f64,
        rec_rst1: bool,
        rec_p1: bool,
        rec_rst2: bool,
        rec_p2: bool,
        psol: &mut f64,
        par_sol: &mut [f64],
    ) -> bool {
        self.done = false;
        self.iscomplete = false;
        self.comptra = false;
        self.line = BRepBlendLine::new();
        self.tolpoint3d = tol3d;
        self.tolgui = tolguide.abs();
        self.rebrou = false;

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        let mut recadp1;
        let mut recadp2;
        let mut recadrst1;
        let mut recadrst2;
        let mut wp1;
        let mut wp2;
        let mut wrst1;
        let mut wrst2;
        let mut state = BlendStatus::OnRst12;
        let (mut trst11, mut trst12, mut trst21, mut trst22) = (0.0, 0.0, 0.0, 0.0);
        let mut infbound = vec![0.0; 2];
        let mut supbound = vec![0.0; 2];
        let mut tolerance = vec![0.0; 2];
        let mut solinvp1 = vec![0.0; 2];
        let mut solinvp2 = vec![0.0; 2];
        let mut solinvrst1 = vec![0.0; 3];
        let mut solinvrst2 = vec![0.0; 3];
        let mut vtxp1 = DomainVertex::empty();
        let mut vtxp2 = DomainVertex::empty();
        let mut vtxrst1 = DomainVertex::empty();
        let mut vtxrst2 = DomainVertex::empty();
        let mut is_vtxp1 = false;
        let mut is_vtxp2 = false;
        let mut is_vtxrst1 = false;
        let mut is_vtxrst2 = false;
        wp1 = pmax;
        wp2 = pmax;
        wrst1 = pmax;
        wrst2 = pmax;
        self.param = pdep;
        func.set_param(self.param);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        let mut rsnld = FunctionSetRoot::new(func, &tolerance, 30);
        rsnld.perform(func, par_dep, &infbound, &supbound, false);
        if !rsnld.is_done() {
            return false;
        }
        self.sol = rsnld.root();

        recadrst1 = rec_rst1
            && self.recadre1_func(func, finv1, &mut solinvrst1, &mut is_vtxrst1, &mut vtxrst1);
        if recadrst1 {
            wrst1 = solinvrst1[0];
        }

        recadp1 = rec_p1 && self.recadre1_finv_p(finv_p1, &mut solinvp1, &mut is_vtxp1, &mut vtxp1);
        if recadp1 {
            wp1 = solinvp1[0];
        }

        recadrst2 = rec_rst2
            && self.recadre2_func(func, finv2, &mut solinvrst2, &mut is_vtxrst2, &mut vtxrst2);
        if recadrst2 {
            wrst2 = solinvrst2[0];
        }

        recadp2 = rec_p2 && self.recadre2_finv_p(finv_p2, &mut solinvp2, &mut is_vtxp2, &mut vtxp2);
        if recadp2 {
            wp2 = solinvp2[0];
        }

        if !recadrst1 && !recadp1 && !recadrst2 && !recadp2 {
            return false;
        }

        // it is checked if the contact was lost or domain 1 was left
        if recadp1 && recadrst1 {
            if self.sens * (wrst1 - wp1) > self.tolgui {
                // at first one leaves the domain
                wrst1 = wp1;
                trst12 = solinvp1[1];
                trst11 = blend_tool_parameter(&vtxp1, &self.rst1);
                is_vtxrst2 = is_vtxp1;
                vtxrst2 = vtxp1.clone();
                recadrst1 = false;
            } else {
                // the contact is lost
                trst11 = solinvrst1[2];
                trst12 = solinvrst1[1];
                recadp1 = false;
            }
        } else if recadp1 {
            wrst1 = wp1;
            trst12 = solinvp1[1];
            trst11 = blend_tool_parameter(&vtxp1, &self.rst1);
            is_vtxrst1 = is_vtxp1;
            vtxrst1 = vtxp1.clone();
        } else if recadrst1 {
            trst11 = solinvrst1[2];
            trst12 = solinvrst1[1];
        }

        // it is checked if the contact was lost or domain 2 was left
        if recadp2 && recadrst2 {
            if self.sens * (wrst2 - wp2) > self.tolgui {
                // at first one leaves the domain
                wrst2 = wp2;
                trst21 = solinvp2[1];
                trst22 = blend_tool_parameter(&vtxp2, &self.rst2);
                is_vtxrst2 = is_vtxp2;
                vtxrst2 = vtxp2.clone();
                recadrst2 = false;
            } else {
                trst22 = solinvrst2[2];
                trst21 = solinvrst2[1];
                recadp2 = false;
            }
        } else if recadp2 {
            wrst2 = wp2;
            trst21 = solinvp2[1];
            trst22 = blend_tool_parameter(&vtxp2, &self.rst2);
            is_vtxrst2 = is_vtxp2;
            vtxrst2 = vtxp2.clone();
        } else if recadrst2 {
            trst22 = solinvrst2[2];
            trst21 = solinvrst2[1];
        }

        // it is checked on which curve the contact is lost earlier
        if recadrst1 && recadrst2 {
            if (wrst1 - wrst2).abs() < self.tolgui {
                state = BlendStatus::OnRst12;
                self.param = 0.5 * (wrst1 + wrst2);
                self.sol[0] = trst11;
                self.sol[1] = trst22;
            } else if self.sens * (wrst1 - wrst2) < 0.0 {
                // contact lost on Rst1
                state = BlendStatus::OnRst1;
                self.param = wrst1;
                self.sol[0] = trst11;
                self.sol[1] = trst12;
            } else {
                // contact lost on rst2
                state = BlendStatus::OnRst2;
                self.param = wrst2;
                self.sol[0] = trst21;
                self.sol[1] = trst22;
            }
            func.set_param(self.param);
        } else if recadrst1 {
            // ground on rst1
            state = BlendStatus::OnRst1;
            self.param = wrst1;
            self.sol[0] = trst11;
            self.sol[1] = trst12;
            func.set_param(self.param);
        } else if recadrst2 {
            // ground on rst2
            state = BlendStatus::OnRst2;
            self.param = wrst2;
            self.sol[0] = trst21;
            self.sol[1] = trst22;
            func.set_param(self.param);
        }
        // it is checked on which curves one leaves first
        else if recadp1 && recadp2 {
            if (wrst1 - wrst2).abs() < self.tolgui {
                state = BlendStatus::OnRst12;
                self.param = 0.5 * (wrst1 + wrst2);
                self.sol[0] = trst11;
                self.sol[1] = trst22;
            } else if self.sens * (wrst1 - wrst2) < 0.0 {
                // sol on Rst1
                state = BlendStatus::OnRst1;
                self.param = wrst1;
                self.sol[0] = trst11;
                self.sol[1] = trst12;
            } else {
                // ground on rst2
                state = BlendStatus::OnRst2;
                self.param = wrst2;
                self.sol[0] = trst21;
                self.sol[1] = trst22;
            }
            func.set_param(self.param);
        } else if recadp1 {
            // ground on rst1
            state = BlendStatus::OnRst1;
            self.param = wrst1;
            self.sol[0] = trst11;
            self.sol[1] = trst12;
            func.set_param(self.param);
        } else if recadp2 {
            // ground on rst2
            state = BlendStatus::OnRst2;
            self.param = wrst2;
            self.sol[0] = trst21;
            self.sol[1] = trst22;
            func.set_param(self.param);
        }

        let state = self.test_arret(func, false, state);
        *psol = self.param;
        par_sol.copy_from_slice(&self.sol);
        let _ = state;
        true
    }

    /// OCCT Complete(Func, Finv1, FinvP1, Finv2, FinvP2, Pmin)
    /// (cxx L538-569).
    pub fn complete(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        pmin: f64,
    ) -> bool {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_RstRstLineBuilder::Complete");
        }
        if self.iscomplete {
            return true;
        }
        self.previous_p = if self.sens > 0.0 {
            self.line.point(1).clone()
        } else {
            self.line.point(self.line.nb_points()).clone()
        };
        self.sens = -self.sens;
        self.param = self.previous_p.parameter();
        self.sol[0] = self.previous_p.parameter_on_c1();
        self.sol[1] = self.previous_p.parameter_on_c2();

        self.internal_perform(func, finv1, finv_p1, finv2, finv_p2, pmin);
        self.iscomplete = true;
        true
    }

    /// OCCT InternalPerform(Func, Finv1, FinvP1, Finv2, FinvP2, Bound)
    /// (cxx L576-1135) — algorithm of processing without extremities.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    fn internal_perform(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        bound: f64,
    ) {
        let mut stepw = self.pasmax;
        let nbp = self.line.nb_points();
        if nbp >= 2 {
            // The last step is redone if it is not too small.
            if self.sens < 0.0 {
                stepw = self.line.point(2).parameter() - self.line.point(1).parameter();
            } else {
                stepw = self.line.point(nbp).parameter() - self.line.point(nbp - 1).parameter();
            }
            stepw = stepw.max(100.0 * self.tolgui);
        }
        let mut parprec = self.param;
        if self.sens * (parprec - bound) >= -self.tolgui {
            return;
        }
        let mut state = BlendStatus::OnRst12;
        let (mut trst11, mut trst12, mut trst21, mut trst22) = (0.0, 0.0, 0.0, 0.0);
        let mut situonc1 = TopAbsState::Unknown;
        let mut situonc2 = TopAbsState::Unknown;
        let mut decroch = BlendDecrochStatus::NoDecroch;
        let mut arrive;
        let mut recadp1;
        let mut recadp2;
        let mut recadrst1;
        let mut recadrst2;
        let mut echecrecad;
        let mut wp1;
        let mut wp2;
        let mut wrst1;
        let mut wrst2;
        let mut infbound = vec![0.0; 2];
        let mut supbound = vec![0.0; 2];
        let mut parinit = vec![0.0; 2];
        let mut tolerance = vec![0.0; 2];
        let mut solinvp1 = vec![0.0; 2];
        let mut solinvp2 = vec![0.0; 2];
        let mut solinvrst1 = vec![0.0; 3];
        let mut solinvrst2 = vec![0.0; 3];
        let mut vtxp1 = DomainVertex::empty();
        let mut vtxp2 = DomainVertex::empty();
        let mut vtxrst1 = DomainVertex::empty();
        let mut vtxrst2 = DomainVertex::empty();
        let mut is_vtxp1 = false;
        let mut is_vtxp2 = false;
        let mut is_vtxrst1 = false;
        let mut is_vtxrst2 = false;
        let mut extrst1 = BRepBlendExtremity::new();
        let mut extrst2 = BRepBlendExtremity::new();

        // IntSurf_Transition Tline, Tarc;

        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        let mut rsnld = FunctionSetRoot::new(func, &tolerance, 30);
        parinit = self.sol.clone();

        arrive = false;
        self.param = parprec + self.sens * stepw;
        if self.sens * (self.param - bound) > 0.0 {
            stepw = self.sens * (bound - parprec) * 0.5;
            self.param = parprec + self.sens * stepw;
        }

        while !arrive {
            let mut bonpoint = true;
            func.set_param(self.param);
            rsnld.perform(func, &parinit, &infbound, &supbound, false);

            if rsnld.is_done() {
                self.sol = rsnld.root();
                if !self.check_inside(func, &mut situonc1, &mut situonc2, &mut decroch)
                    && self.line.nb_points() == 1
                {
                    state = BlendStatus::StepTooLarge;
                    bonpoint = false;
                }
            } else {
                state = BlendStatus::StepTooLarge;
                bonpoint = false;
            }
            if bonpoint {
                wp1 = bound;
                wp2 = bound;
                wrst1 = bound;
                wrst2 = bound;
                recadp1 = false;
                recadp2 = false;
                recadrst1 = false;
                recadrst2 = false;
                echecrecad = false;
                if situonc1 != TopAbsState::In {
                    // pb inversion rst/rst
                    recadp1 = self.recadre1_finv_p(finv_p1, &mut solinvp1, &mut is_vtxp1, &mut vtxp1);
                    if recadp1 {
                        wp1 = solinvp1[0];
                    } else {
                        echecrecad = true;
                    }
                }

                if situonc2 != TopAbsState::In {
                    // pb inversion point/surf
                    recadp2 = self.recadre2_finv_p(finv_p2, &mut solinvp2, &mut is_vtxp2, &mut vtxp2);
                    if recadp2 {
                        wp2 = solinvp2[0];
                    } else {
                        echecrecad = true;
                    }
                }

                if decroch == BlendDecrochStatus::DecrochRst1
                    || decroch == BlendDecrochStatus::DecrochBoth
                {
                    // pb inversion rst1/surf1
                    recadrst1 =
                        self.recadre1_func(func, finv1, &mut solinvrst1, &mut is_vtxrst1, &mut vtxrst1);
                    if recadrst1 {
                        wrst1 = solinvrst1[0];
                    } else {
                        echecrecad = true;
                    }
                }

                if decroch == BlendDecrochStatus::DecrochRst2
                    || decroch == BlendDecrochStatus::DecrochBoth
                {
                    // pb inverse rst2/surf2
                    recadrst2 =
                        self.recadre2_func(func, finv2, &mut solinvrst2, &mut is_vtxrst2, &mut vtxrst2);
                    if recadrst2 {
                        wrst2 = solinvrst2[0];
                    } else {
                        echecrecad = true;
                    }
                }

                decroch = BlendDecrochStatus::NoDecroch;
                if recadp1 || recadp2 || recadrst1 || recadrst2 {
                    echecrecad = false;
                }

                if !echecrecad {
                    // it is checked if the contact was lost or domain 1 was left
                    if recadp1 && recadrst1 {
                        if self.sens * (wrst1 - wp1) > self.tolgui {
                            // first one leaves the domain
                            wrst1 = wp1;
                            trst12 = solinvp1[1];
                            trst11 = blend_tool_parameter(&vtxp1, &self.rst1);
                            is_vtxrst2 = is_vtxp1;
                            vtxrst2 = vtxp1.clone();
                            recadrst1 = false;
                        } else {
                            // contact is lost
                            trst11 = solinvrst1[2];
                            trst12 = solinvrst1[1];
                            recadp1 = false;
                        }
                    } else if recadp1 {
                        wrst1 = wp1;
                        trst12 = solinvp1[1];
                        trst11 = blend_tool_parameter(&vtxp1, &self.rst1);
                        is_vtxrst1 = is_vtxp1;
                        vtxrst1 = vtxp1.clone();
                    } else if recadrst1 {
                        trst11 = solinvrst1[2];
                        trst12 = solinvrst1[1];
                    }

                    // it is checked if the contact was lost or domain 2 was left
                    if recadp2 && recadrst2 {
                        if self.sens * (wrst2 - wp2) > self.tolgui {
                            // first one leaves the domain
                            wrst2 = wp2;
                            trst21 = solinvp2[1];
                            trst22 = blend_tool_parameter(&vtxp2, &self.rst2);
                            is_vtxrst2 = is_vtxp2;
                            vtxrst2 = vtxp2.clone();
                            recadrst2 = false;
                        } else {
                            trst22 = solinvrst2[2];
                            trst21 = solinvrst2[1];
                            recadp2 = false;
                        }
                    } else if recadp2 {
                        wrst2 = wp2;
                        trst21 = solinvp2[1];
                        trst22 = blend_tool_parameter(&vtxp2, &self.rst2);
                        is_vtxrst2 = is_vtxp2;
                        vtxrst2 = vtxp2.clone();
                    } else if recadrst2 {
                        trst22 = solinvrst2[2];
                        trst21 = solinvrst2[1];
                    }

                    // it is checked on which curve the contact is lost earlier
                    if recadrst1 && recadrst2 {
                        if (wrst1 - wrst2).abs() < self.tolgui {
                            state = BlendStatus::OnRst12;
                            decroch = BlendDecrochStatus::DecrochBoth;
                            self.param = 0.5 * (wrst1 + wrst2);
                            self.sol[0] = trst11;
                            self.sol[1] = trst22;
                        } else if self.sens * (wrst1 - wrst2) < 0.0 {
                            // contact is lost on Rst1
                            state = BlendStatus::OnRst1;
                            decroch = BlendDecrochStatus::DecrochRst1;
                            self.param = wrst1;
                            self.sol[0] = trst11;
                            self.sol[1] = trst12;
                        } else {
                            // contact is lost on rst2
                            state = BlendStatus::OnRst2;
                            decroch = BlendDecrochStatus::DecrochRst2;
                            self.param = wrst2;
                            self.sol[0] = trst21;
                            self.sol[1] = trst22;
                        }
                        func.set_param(self.param);
                    } else if recadrst1 {
                        // ground on rst1
                        state = BlendStatus::OnRst1;
                        decroch = BlendDecrochStatus::DecrochRst1;
                        self.param = wrst1;
                        self.sol[0] = trst11;
                        self.sol[1] = trst12;
                        func.set_param(self.param);
                    } else if recadrst2 {
                        // ground on rst2
                        state = BlendStatus::OnRst2;
                        decroch = BlendDecrochStatus::DecrochRst2;
                        self.param = wrst2;
                        self.sol[0] = trst21;
                        self.sol[1] = trst22;
                        func.set_param(self.param);
                    }
                    //  it is checked on which curve the contact is lost earlier
                    else if recadp1 && recadp2 {
                        if (wrst1 - wrst2).abs() < self.tolgui {
                            state = BlendStatus::OnRst12;
                            self.param = 0.5 * (wrst1 + wrst2);
                            self.sol[0] = trst11;
                            self.sol[1] = trst22;
                        } else if self.sens * (wrst1 - wrst2) < 0.0 {
                            // ground on Rst1
                            state = BlendStatus::OnRst1;
                            self.param = wrst1;
                            self.sol[0] = trst11;
                            self.sol[1] = trst12;
                        } else {
                            // ground on rst2
                            state = BlendStatus::OnRst2;
                            self.param = wrst2;
                            self.sol[0] = trst21;
                            self.sol[1] = trst22;
                        }
                        func.set_param(self.param);
                    } else if recadp1 {
                        // ground on rst1
                        state = BlendStatus::OnRst1;
                        self.param = wrst1;
                        self.sol[0] = trst11;
                        self.sol[1] = trst12;
                        func.set_param(self.param);
                    } else if recadp2 {
                        // ground on rst2
                        state = BlendStatus::OnRst2;
                        self.param = wrst2;
                        self.sol[0] = trst21;
                        self.sol[1] = trst22;
                        func.set_param(self.param);
                    } else {
                        state = BlendStatus::Ok;
                    }

                    state = self.test_arret(func, true, state);
                } else {
                    // reframing failed. Leave with PointsConfondus
                    state = BlendStatus::SamePoints;
                }
            }

            // The restriction clones below keep the &mut self calls
            // (append / make_extremity) borrow-conflict free.
            let rst1 = self.rst1.clone();
            let rst2 = self.rst2.clone();
            match state {
                BlendStatus::Ok => {
                    // Update the line.
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }
                    parinit = self.sol.clone();
                    parprec = self.param;

                    if self.param == bound {
                        arrive = true;
                        extrst1.set_value_on_curve(
                            self.previous_p.point_on_c1(),
                            self.previous_p.parameter_on_c1(),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        self.make_extremity(
                            &mut extrst2,
                            false,
                            &rst2,
                            self.sol[1],
                            is_vtxrst2,
                            &vtxrst2,
                        );
                        // Show that end is on Bound.
                    } else {
                        self.param = self.param + self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                }

                BlendStatus::StepTooLarge => {
                    stepw = stepw / 2.0;
                    if stepw.abs() < self.tolgui {
                        extrst1.set_value_on_curve(
                            self.previous_p.point_on_c1(),
                            self.previous_p.parameter_on_c1(),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        extrst2.set_value_on_curve(
                            self.previous_p.point_on_c2(),
                            self.previous_p.parameter_on_c2(),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        arrive = true;
                        if self.line.nb_points() >= 2 {
                            // Show that there is a stop during processing
                        }
                    } else {
                        self.param = parprec + self.sens * stepw; // there is no risk to exceed Bound.
                    }
                }

                BlendStatus::StepTooSmall => {
                    // Update the line.
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }
                    parinit = self.sol.clone();
                    parprec = self.param;

                    stepw = (1.5 * stepw).min(self.pasmax);
                    if self.param == bound {
                        arrive = true;
                        extrst1.set_value_on_curve(
                            self.previous_p.point_on_c1(),
                            self.previous_p.parameter_on_c1(),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        self.make_extremity(
                            &mut extrst2,
                            false,
                            &rst2,
                            self.sol[1],
                            is_vtxrst2,
                            &vtxrst2,
                        );
                        // Indicate that end is on Bound.
                    } else {
                        self.param = self.param + self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                }

                BlendStatus::OnRst1 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }
                    self.make_extremity(
                        &mut extrst1,
                        true,
                        &rst1,
                        self.sol[0],
                        is_vtxrst1,
                        &vtxrst1,
                    );
                    self.make_extremity(
                        &mut extrst2,
                        false,
                        &rst2,
                        self.sol[1],
                        is_vtxrst2,
                        &vtxrst2,
                    );
                    arrive = true;
                }

                BlendStatus::OnRst2 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }

                    self.make_extremity(
                        &mut extrst1,
                        true,
                        &rst1,
                        self.sol[0],
                        is_vtxrst1,
                        &vtxrst1,
                    );
                    self.make_extremity(
                        &mut extrst2,
                        false,
                        &rst2,
                        self.sol[1],
                        is_vtxrst2,
                        &vtxrst2,
                    );
                    arrive = true;
                }

                BlendStatus::OnRst12 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }

                    self.make_extremity(
                        &mut extrst1,
                        true,
                        &rst1,
                        self.sol[0],
                        is_vtxrst1,
                        &vtxrst1,
                    );
                    self.make_extremity(
                        &mut extrst2,
                        false,
                        &rst2,
                        self.sol[1],
                        is_vtxrst2,
                        &vtxrst2,
                    );
                    arrive = true;
                }

                BlendStatus::SamePoints => {
                    // Stop
                    extrst1.set_value_on_curve(
                        self.previous_p.point_on_c1(),
                        self.previous_p.parameter_on_c1(),
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    extrst2.set_value_on_curve(
                        self.previous_p.point_on_c2(),
                        self.previous_p.parameter_on_c2(),
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    arrive = true;
                }

                BlendStatus::Backward => {}
            }
            if arrive {
                if self.sens > 0.0 {
                    self.line.set_end_points(&extrst1, &extrst2);
                    self.decrochfin = decroch;
                } else {
                    self.line.set_start_points(&extrst1, &extrst2);
                    self.decrochdeb = decroch;
                }
            }
        }
    }

    /// OCCT Recadre1(Func, Finv, Solinv, IsVtx, Vtx) (cxx L1139-1220).
    fn recadre1_func(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv: &mut dyn BlendSurfCurvFuncInv,
        solinv: &mut [f64],
        is_vtx: &mut bool,
        vtx: &mut DomainVertex,
    ) -> bool {
        let mut toler = vec![0.0; 3];
        let mut infb = vec![0.0; 3];
        let mut supb = vec![0.0; 3];
        finv.get_tolerance(&mut toler, self.tolpoint3d);
        finv.get_bounds(&mut infb, &mut supb);
        solinv[0] = self.param;
        solinv[1] = self.sol[1];
        solinv[2] = self.sol[0];

        // The point where contact is not lost is found
        let mut rsnld = FunctionSetRoot::new(finv, &toler, 30);
        rsnld.perform(finv, solinv, &infb, &supb, false);
        if !rsnld.is_done() {
            // OCCT prints "RSNLD not done" under OCCT_DEBUG.
            return false;
        }

        let root = rsnld.root();
        solinv.copy_from_slice(&root);

        // It is necessary to check if the function value meets the
        // second restriction
        if finv.is_solution(solinv, self.tolpoint3d) {
            let w = solinv[1];
            if w < self.rst2.first_parameter() - toler[1]
                || w > self.rst2.last_parameter() + toler[1]
            {
                return false;
            }

            // it is checked if it is on a Vertex
            let mut domain1 = DomainTool::single(&self.rst1);
            domain1.init_vertex_iterator();
            *is_vtx = !domain1.more_vertex();
            while !*is_vtx {
                *vtx = domain1.vertex();
                if (blend_tool_parameter(vtx, &self.rst1) - solinv[2]).abs()
                    <= blend_tool_tolerance(vtx, &self.rst1)
                {
                    *is_vtx = true;
                } else {
                    domain1.next_vertex();
                    *is_vtx = !domain1.more_vertex();
                }
            }
            if !domain1.more_vertex() {
                *is_vtx = false;
            }
            // The section is recalculated by direct solution, otherwise return
            // incoherences between the parameter and the ground caused by yawn.

            let mut infbound = vec![0.0; 2];
            let mut supbound = vec![0.0; 2];
            let mut parinit = vec![0.0; 2];
            let mut tolerance = vec![0.0; 2];
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);

            let mut rsnld2 = FunctionSetRoot::new(func, &tolerance, 30);
            parinit[0] = solinv[2];
            parinit[1] = solinv[1];
            func.set_param(solinv[0]);
            rsnld2.perform(func, &parinit, &infbound, &supbound, false);
            if !rsnld2.is_done() {
                return false;
            }
            let root = rsnld2.root();
            parinit.copy_from_slice(&root);
            solinv[1] = parinit[1];
            solinv[2] = parinit[0];
            return true;
        }
        false
    }

    /// OCCT Recadre2(Func, Finv, Solinv, IsVtx, Vtx) (cxx L1224-1302).
    fn recadre2_func(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv: &mut dyn BlendSurfCurvFuncInv,
        solinv: &mut [f64],
        is_vtx: &mut bool,
        vtx: &mut DomainVertex,
    ) -> bool {
        let mut toler = vec![0.0; 3];
        let mut infb = vec![0.0; 3];
        let mut supb = vec![0.0; 3];
        finv.get_tolerance(&mut toler, self.tolpoint3d);
        finv.get_bounds(&mut infb, &mut supb);
        solinv[0] = self.param;
        solinv[1] = self.sol[0];
        solinv[2] = self.sol[1];

        let mut rsnld = FunctionSetRoot::new(finv, &toler, 30);
        rsnld.perform(finv, solinv, &infb, &supb, false);
        if !rsnld.is_done() {
            // OCCT prints "RSNLD not done" under OCCT_DEBUG.
            return false;
        }

        let root = rsnld.root();
        solinv.copy_from_slice(&root);

        // It is necessary to check the value of the function
        if finv.is_solution(solinv, self.tolpoint3d) {
            let w = solinv[1];
            if w < self.rst1.first_parameter() - toler[1]
                || w > self.rst1.last_parameter() + toler[1]
            {
                return false;
            }

            let mut domain2 = DomainTool::single(&self.rst2);
            domain2.init_vertex_iterator();
            *is_vtx = !domain2.more_vertex();
            while !*is_vtx {
                *vtx = domain2.vertex();
                if (blend_tool_parameter(vtx, &self.rst2) - solinv[2]).abs()
                    <= blend_tool_tolerance(vtx, &self.rst2)
                {
                    *is_vtx = true;
                } else {
                    domain2.next_vertex();
                    *is_vtx = !domain2.more_vertex();
                }
            }
            if !domain2.more_vertex() {
                *is_vtx = false;
            }
            // The section is recalculated by direct solution, otherwise return
            // incoherences between the parameter and the ground caused by yawn.

            let mut infbound = vec![0.0; 2];
            let mut supbound = vec![0.0; 2];
            let mut parinit = vec![0.0; 2];
            let mut tolerance = vec![0.0; 2];
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);

            let mut rsnld2 = FunctionSetRoot::new(func, &tolerance, 30);
            parinit[0] = solinv[1];
            parinit[1] = solinv[2];
            func.set_param(solinv[0]);
            rsnld2.perform(func, &parinit, &infbound, &supbound, false);
            if !rsnld2.is_done() {
                return false;
            }
            let root = rsnld2.root();
            parinit.copy_from_slice(&root);
            solinv[1] = parinit[0];
            solinv[2] = parinit[1];
            return true;
        }
        false
    }

    /// OCCT Recadre1(FinvP, Solinv, IsVtx, Vtx) (cxx L1306-1375).
    fn recadre1_finv_p(
        &mut self,
        finv_p: &mut dyn BlendCurvPointFuncInv,
        solinv: &mut [f64],
        is_vtx: &mut bool,
        vtx: &mut DomainVertex,
    ) -> bool {
        // One is located on the last or the first point, following the
        // direction of processing.
        let firstrst1 = self.rst1.first_parameter();
        let lastrst1 = self.rst1.last_parameter();
        let mut upoint = firstrst1;

        if (self.sol[0] - firstrst1) > (lastrst1 - self.sol[0]) {
            upoint = lastrst1;
        }
        let p2drst1 = self.rst1.value(upoint);
        let thepoint = self.surf1.surface.point_at(p2drst1.x, p2drst1.y);

        finv_p.set_point(thepoint);
        let mut toler = vec![0.0; 2];
        let mut infb = vec![0.0; 2];
        let mut supb = vec![0.0; 2];
        finv_p.get_tolerance(&mut toler, self.tolpoint3d);
        finv_p.get_bounds(&mut infb, &mut supb);
        solinv[0] = self.param;
        solinv[1] = self.sol[1];

        let mut rsnld = FunctionSetRoot::new(finv_p, &toler, 30);
        rsnld.perform(finv_p, solinv, &infb, &supb, false);
        if !rsnld.is_done() {
            // OCCT prints "RSNLD not done" under OCCT_DEBUG.
            return false;
        }
        let root = rsnld.root();
        solinv.copy_from_slice(&root);

        if finv_p.is_solution(solinv, self.tolpoint3d) {
            let p2drst2 = self.rst2.value(solinv[1]);
            let situ = self.domain2.classify(p2drst2, toler[1], false);
            if situ != TopAbsState::In && situ != TopAbsState::On {
                return false;
            }
            let mut domain1 = DomainTool::single(&self.rst1);
            domain1.init_vertex_iterator();
            *is_vtx = !domain1.more_vertex();
            while !*is_vtx {
                *vtx = domain1.vertex();
                if (blend_tool_parameter(vtx, &self.rst1) - upoint).abs()
                    <= blend_tool_tolerance(vtx, &self.rst1)
                {
                    *is_vtx = true;
                } else {
                    domain1.next_vertex();
                    *is_vtx = !domain1.more_vertex();
                }
            }
            if !domain1.more_vertex() {
                *is_vtx = false;
            }
            return true;
        }
        false
    }

    /// OCCT Recadre2(FinvP, Solinv, IsVtx, Vtx) (cxx L1379-1448).
    fn recadre2_finv_p(
        &mut self,
        finv_p: &mut dyn BlendCurvPointFuncInv,
        solinv: &mut [f64],
        is_vtx: &mut bool,
        vtx: &mut DomainVertex,
    ) -> bool {
        // One is located on the last or the first point, following the
        // direction of processing.
        let firstrst2 = self.rst2.first_parameter();
        let lastrst2 = self.rst2.last_parameter();
        let mut vpoint = firstrst2;

        if (self.sol[1] - firstrst2) > (lastrst2 - self.sol[1]) {
            vpoint = lastrst2;
        }
        let p2drst2 = self.rst2.value(vpoint);
        let thepoint = self.surf2.surface.point_at(p2drst2.x, p2drst2.y);

        finv_p.set_point(thepoint);
        let mut toler = vec![0.0; 2];
        let mut infb = vec![0.0; 2];
        let mut supb = vec![0.0; 2];
        finv_p.get_tolerance(&mut toler, self.tolpoint3d);
        finv_p.get_bounds(&mut infb, &mut supb);
        solinv[0] = self.param;
        solinv[1] = self.sol[0];

        let mut rsnld = FunctionSetRoot::new(finv_p, &toler, 30);
        rsnld.perform(finv_p, solinv, &infb, &supb, false);
        if !rsnld.is_done() {
            // OCCT prints "RSNLD not done" under OCCT_DEBUG.
            return false;
        }
        let root = rsnld.root();
        solinv.copy_from_slice(&root);

        if finv_p.is_solution(solinv, self.tolpoint3d) {
            let p2drst1 = self.rst1.value(solinv[1]);
            let situ = self.domain1.classify(p2drst1, toler[1], false);
            if situ != TopAbsState::In && situ != TopAbsState::On {
                return false;
            }
            let mut domain2 = DomainTool::single(&self.rst2);
            domain2.init_vertex_iterator();
            *is_vtx = !domain2.more_vertex();
            while !*is_vtx {
                *vtx = domain2.vertex();
                if (blend_tool_parameter(vtx, &self.rst2) - vpoint).abs()
                    <= blend_tool_tolerance(vtx, &self.rst2)
                {
                    *is_vtx = true;
                } else {
                    domain2.next_vertex();
                    *is_vtx = !domain2.more_vertex();
                }
            }
            if !domain2.more_vertex() {
                *is_vtx = false;
            }
            return true;
        }
        false
    }

    /// OCCT Transition(OnFirst, Arc, Param, TLine, TArc) (cxx L1452-1514).
    fn transition(
        &self,
        on_first: bool,
        arc: &RstArc,
        param: f64,
        t_line: &mut IntSurfTransition,
        t_arc: &mut IntSurfTransition,
    ) {
        let mut computetranstionaveclacorde = false;
        let tgline;
        let prevprev;

        if self.previous_p.is_tangency_point() {
            if self.line.nb_points() < 2 {
                return;
            }
            computetranstionaveclacorde = true;
            prevprev = if self.sens < 0.0 {
                self.line.point(2).clone()
            } else {
                self.line.point(self.line.nb_points() - 1).clone()
            };
        } else {
            prevprev = BlendPoint::new();
        }
        let (p2d, dp2d) = arc.d1(param);

        let (d1u, d1v);
        if on_first {
            let (_pbid, du, dv) = self.surf1.surface.derivatives(p2d.x, p2d.y);
            d1u = du;
            d1v = dv;
            if !computetranstionaveclacorde {
                tgline = self.previous_p.tangent_on_c1();
            } else {
                tgline = self.previous_p.point_on_c1() - prevprev.point_on_c1();
            }
        } else {
            let (_pbid, du, dv) = self.surf2.surface.derivatives(p2d.x, p2d.y);
            d1u = du;
            d1v = dv;
            if !computetranstionaveclacorde {
                tgline = self.previous_p.tangent_on_c2();
            } else {
                tgline = self.previous_p.point_on_c2() - prevprev.point_on_c2();
            }
        }

        let tgrst = d1u * dp2d.x + d1v * dp2d.y;
        let normale = d1u.cross(d1v);

        make_transition(tgline, tgrst, normale, t_line, t_arc);
    }

    /// OCCT MakeExtremity(Extrem, OnFirst, Arc, Param, IsVtx, Vtx)
    /// (cxx L1521-1585) — produce the extremity of a curve.
    fn make_extremity(
        &mut self,
        extrem: &mut BRepBlendExtremity,
        on_first: bool,
        arc: &RstArc,
        param: f64,
        is_vtx: bool,
        vtx: &DomainVertex,
    ) {
        let mut tline = IntSurfTransition::new();
        let mut tarc = IntSurfTransition::new();
        let mut prm;
        let mut iter = if on_first {
            extrem.set_value_on_curve(
                self.previous_p.point_on_c1(),
                self.sol[0],
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_c1());
            }
            DomainTool::face_domain(self.domain1)
        } else {
            extrem.set_value_on_curve(
                self.previous_p.point_on_c2(),
                self.sol[1],
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            // OCCT literal (cxx L1545): the else branch reads TangentOnC1().
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_c1());
            }
            DomainTool::face_domain(self.domain2)
        };

        iter.init();
        if !is_vtx {
            self.transition(on_first, arc, param, &mut tline, &mut tarc);
            extrem.add_arc(&arc.curve, param, tline, tarc);
        } else {
            extrem.set_vertex(&vtx.hv);
            while iter.more() {
                let arc_cur = iter.value();
                if !arc_cur.is_same(arc) {
                    iter.initialize_arc(&arc_cur);
                    iter.init_vertex_iterator();
                    while iter.more_vertex() {
                        if iter.vertex().is_same(vtx) {
                            prm = blend_tool_parameter(vtx, &arc_cur);
                            self.transition(on_first, &arc_cur, prm, &mut tline, &mut tarc);
                            extrem.add_arc(&arc_cur.curve, prm, tline, tarc);
                        }
                        iter.next_vertex();
                    }
                } else {
                    self.transition(on_first, &arc_cur, param, &mut tline, &mut tarc);
                    extrem.add_arc(&arc_cur.curve, param, tline, tarc);
                }
                iter.next();
            }
        }
    }

    /// OCCT CheckDeflectionOnRst1(CurPoint) (cxx L1589-1673) — controls 3d
    /// of Blend_CSWalking.
    fn check_deflection_on_rst1(&mut self, cur_point: &BlendPoint) -> BlendStatus {
        // rule by tests in U4 corresponds to 11.478
        let cos_ref_3d = 0.98;
        let mut cosi;
        let mut cosi2;
        let curpointistangent = cur_point.is_tangency_point();
        let prevpointistangent = self.previous_p.is_tangency_point();

        let psurf = cur_point.point_on_c1();
        let tgsurf;
        if !curpointistangent {
            tgsurf = cur_point.tangent_on_c1();
        } else {
            tgsurf = DVec3::ZERO;
        }
        let prevp = self.previous_p.point_on_c1();
        let prevtg;
        if !prevpointistangent {
            prevtg = self.previous_p.tangent_on_c1();
        } else {
            prevtg = DVec3::ZERO;
        }
        let norme;
        let mut prevnorme = 0.0;
        let corde = psurf - prevp;
        norme = corde.length_squared();
        if !prevpointistangent {
            prevnorme = prevtg.length_squared();
        }

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // it can be necessary to force the same point
            return BlendStatus::SamePoints;
        }
        if !prevpointistangent {
            if prevnorme <= toler3d * toler3d {
                return BlendStatus::SamePoints;
            }
            cosi = self.sens * corde.dot(prevtg);
            if cosi < 0.0 {
                // angle 3d>pi/2. --> return back
                return BlendStatus::Backward;
            }

            cosi2 = cosi * cosi / prevnorme / norme;
            if cosi2 < cos_ref_3d {
                return BlendStatus::StepTooLarge;
            }
        }

        if !curpointistangent {
            // Check if it is necessary to control the sign of prevtg*Tgsurf
            cosi = self.sens * corde.dot(tgsurf);
            cosi2 = cosi * cosi / tgsurf.length_squared() / norme;
            if cosi2 < cos_ref_3d || cosi < 0.0 {
                return BlendStatus::StepTooLarge;
            }
        }

        if !curpointistangent && !prevpointistangent {
            // Estimation of the current arrow
            let fleche_courante =
                (prevtg.normalize() - tgsurf.normalize()).length_squared() * norme / 64.0;

            if fleche_courante <= 0.25 * self.fleche * self.fleche {
                return BlendStatus::StepTooSmall;
            }
            if fleche_courante > self.fleche * self.fleche {
                // not too great
                return BlendStatus::StepTooLarge;
            }
        }
        BlendStatus::Ok
    }

    /// OCCT CheckDeflectionOnRst2(CurPoint) (cxx L1677-1762) — 3D controls
    /// of Blend_CSWalking.
    fn check_deflection_on_rst2(&mut self, cur_point: &BlendPoint) -> BlendStatus {
        // rule by tests in U4 corresponding to 11.478 d
        let cos_ref_3d = 0.98;
        let mut cosi;
        let mut cosi2;
        let curpointistangent = cur_point.is_tangency_point();
        let prevpointistangent = self.previous_p.is_tangency_point();

        let psurf = cur_point.point_on_c2();
        let tgsurf;
        if !curpointistangent {
            tgsurf = cur_point.tangent_on_c2();
        } else {
            tgsurf = DVec3::ZERO;
        }
        let prevp = self.previous_p.point_on_c2();
        let prevtg;
        if !prevpointistangent {
            prevtg = self.previous_p.tangent_on_c2();
        } else {
            prevtg = DVec3::ZERO;
        }
        let norme;
        let mut prevnorme = 0.0;
        let corde = psurf - prevp;
        norme = corde.length_squared();
        if !prevpointistangent {
            prevnorme = prevtg.length_squared();
        }

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // it can be necessary to force the same point
            return BlendStatus::SamePoints;
        }
        if !prevpointistangent {
            if prevnorme <= toler3d * toler3d {
                return BlendStatus::SamePoints;
            }
            cosi = self.sens * corde.dot(prevtg);
            if cosi < 0.0 {
                // angle 3d>pi/2. --> return back
                return BlendStatus::Backward;
            }

            cosi2 = cosi * cosi / prevnorme / norme;
            if cosi2 < cos_ref_3d {
                return BlendStatus::StepTooLarge;
            }
        }

        if !curpointistangent {
            // Check if it is necessary to control the sign of prevtg*Tgsurf
            cosi = self.sens * corde.dot(tgsurf);
            cosi2 = cosi * cosi / tgsurf.length_squared() / norme;
            if cosi2 < cos_ref_3d || cosi < 0.0 {
                return BlendStatus::StepTooLarge;
            }
        }

        if !curpointistangent && !prevpointistangent {
            // Estimation of the current arrow
            let fleche_courante =
                (prevtg.normalize() - tgsurf.normalize()).length_squared() * norme / 64.0;

            if fleche_courante <= 0.25 * self.fleche * self.fleche {
                return BlendStatus::StepTooSmall;
            }
            if fleche_courante > self.fleche * self.fleche {
                // not too great
                return BlendStatus::StepTooLarge;
            }
        }
        BlendStatus::Ok
    }

    /// OCCT TestArret(Func, TestDeflection, State) (cxx L1775-1916).
    #[allow(unused_assignments)]
    fn test_arret(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        test_deflection: bool,
        state: BlendStatus,
    ) -> BlendStatus {
        let mut tgrst1 = DVec3::ZERO;
        let mut tgrst2 = DVec3::ZERO;
        let mut tg2drst1 = DVec2::ZERO;
        let mut tg2drst2 = DVec2::ZERO;
        let mut trarst1 = IntSurfTypeTrans::Undecided;
        let mut trarst2 = IntSurfTypeTrans::Undecided;

        if func.is_solution(&self.sol, self.tolpoint3d) {
            let curpointistangent = func.is_tangency_point();
            let ptrst1 = func.point_on_rst1();
            let ptrst2 = func.point_on_rst2();
            let pt2drst1 = func.pnt2d_on_rst1();
            let pt2drst2 = func.pnt2d_on_rst2();

            let curpoint;
            if curpointistangent {
                curpoint = BlendPoint::new_on_2_curves_on_surfaces(
                    ptrst1,
                    ptrst2,
                    self.param,
                    pt2drst1.x,
                    pt2drst1.y,
                    pt2drst2.x,
                    pt2drst2.y,
                    self.sol[0],
                    self.sol[1],
                );
            } else {
                tgrst1 = func.tangent_on_rst1();
                tgrst2 = func.tangent_on_rst2();
                tg2drst1 = func.tangent_2d_on_rst1();
                tg2drst2 = func.tangent_2d_on_rst2();
                curpoint = BlendPoint::new_on_2_curves_on_surfaces_with_tangents(
                    ptrst1,
                    ptrst2,
                    self.param,
                    pt2drst1.x,
                    pt2drst1.y,
                    pt2drst2.x,
                    pt2drst2.y,
                    self.sol[0],
                    self.sol[1],
                    tgrst1,
                    tgrst2,
                    tg2drst1,
                    tg2drst2,
                );
            }
            let mut state_rst1;
            let mut state_rst2;
            if test_deflection {
                state_rst1 = self.check_deflection_on_rst1(&curpoint);
                state_rst2 = self.check_deflection_on_rst2(&curpoint);
            } else {
                state_rst1 = BlendStatus::Ok;
                state_rst2 = BlendStatus::Ok;
            }
            if state_rst1 == BlendStatus::Backward {
                state_rst1 = BlendStatus::StepTooLarge;
                self.rebrou = true;
            }
            if state_rst2 == BlendStatus::Backward {
                state_rst2 = BlendStatus::StepTooLarge;
                self.rebrou = true;
            }
            if state_rst1 == BlendStatus::StepTooLarge || state_rst2 == BlendStatus::StepTooLarge {
                return BlendStatus::StepTooLarge;
            }

            if !self.comptra && !curpointistangent {
                let (_p2drstref, tg2drstref) = self.rst1.d1(self.sol[0]);
                let mut testra = tg2drst1.dot(tg2drstref);
                let mut or = self.rst1.orientation();

                if testra.abs() > self.tolpoint3d {
                    if testra < 0.0 {
                        trarst1 = conv_or_to_tra(topabs_reverse(or));
                    } else if testra > 0.0 {
                        trarst1 = conv_or_to_tra(or);
                    }

                    let (_p2drstref, tg2drstref) = self.rst2.d1(self.sol[1]);
                    testra = tg2drst2.dot(tg2drstref);

                    or = self.rst2.orientation();
                    if testra.abs() > self.tolpoint3d {
                        if testra < 0.0 {
                            trarst2 = conv_or_to_tra(topabs_reverse(or));
                        } else if testra > 0.0 {
                            trarst2 = conv_or_to_tra(or);
                        }
                        self.comptra = true;
                        self.line.set_transitions(trarst1, trarst2);
                    }
                }
            }
            if state_rst1 == BlendStatus::Ok || state_rst2 == BlendStatus::Ok {
                self.previous_p = curpoint;
                return state;
            }
            if state_rst1 == BlendStatus::StepTooSmall && state_rst2 == BlendStatus::StepTooSmall {
                self.previous_p = curpoint;
                if state == BlendStatus::Ok {
                    return BlendStatus::StepTooSmall;
                } else {
                    return state;
                }
            }
            if state == BlendStatus::Ok {
                return BlendStatus::SamePoints;
            } else {
                return state;
            }
        }
        BlendStatus::StepTooLarge
    }

    /// OCCT CheckInside(Func, SituOnC1, SituOnC2, Decroch) (cxx L1920-1964).
    fn check_inside(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        situ_on_c1: &mut TopAbsState,
        situ_on_c2: &mut TopAbsState,
        decroch: &mut BlendDecrochStatus,
    ) -> bool {
        //  bool inside = true;
        let mut tolerance = vec![0.0; 2];
        func.get_tolerance(&mut tolerance, self.tolpoint3d);

        // face pcurve 1.
        let mut v = self.sol[0];
        if v < self.rst1.first_parameter() - tolerance[1]
            || v > self.rst1.last_parameter() + tolerance[1]
        {
            *situ_on_c1 = TopAbsState::Out;
        } else if v > self.rst1.first_parameter() && v < self.rst1.last_parameter() {
            *situ_on_c1 = TopAbsState::In;
        } else {
            *situ_on_c1 = TopAbsState::On;
        }

        // face pcurve 2.
        v = self.sol[1];
        if v < self.rst2.first_parameter() - tolerance[1]
            || v > self.rst2.last_parameter() + tolerance[1]
        {
            *situ_on_c2 = TopAbsState::Out;
        } else if v > self.rst2.first_parameter() && v < self.rst2.last_parameter() {
            *situ_on_c2 = TopAbsState::In;
        } else {
            *situ_on_c2 = TopAbsState::On;
        }

        // lost contact
        // OCCT (cxx L1961): Decroch = Func.Decroch(sol, tgrst1, norst1, tgrst2, norst2);
        // PENDING: Blend_RstRstFunction::Decroch is missing from the rcad
        // BlendRstRstFunction trait (shared file outside this batch's
        // ownership — gap reported); the stand-in reports Blend_NoDecroch.
        let tgrst1 = DVec3::ZERO;
        let norst1 = DVec3::ZERO;
        let tgrst2 = DVec3::ZERO;
        let norst2 = DVec3::ZERO;
        let _ = (tgrst1, norst1, tgrst2, norst2);
        *decroch = BlendDecrochStatus::NoDecroch;

        *situ_on_c1 == TopAbsState::In
            && *situ_on_c2 == TopAbsState::In
            && *decroch == BlendDecrochStatus::NoDecroch
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT anchor: the constructor initial state (cxx L151-178).
    #[test]
    fn ctor_anchor() {
        let surf = BRepAdaptorSurface::empty();
        let tool = BRepTopAdaptorTopolTool::default();
        let rst = Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let b = BRepBlendRstRstLineBuilder::new(&surf, &rst, &tool, &surf, &rst, &tool);
        assert!(!b.is_done());
        assert!(!b.decroch1_start());
        assert!(!b.decroch1_end());
        assert!(!b.decroch2_start());
        assert!(!b.decroch2_end());
    }

    /// OCCT anchor: ConvOrToTra (cxx L1764-1771).
    #[test]
    fn conv_or_to_tra_anchor() {
        assert_eq!(
            conv_or_to_tra(rcad_kernel::topods::Orientation::Forward),
            IntSurfTypeTrans::In
        );
        assert_eq!(
            conv_or_to_tra(rcad_kernel::topods::Orientation::Reversed),
            IntSurfTypeTrans::Out
        );
    }
}
