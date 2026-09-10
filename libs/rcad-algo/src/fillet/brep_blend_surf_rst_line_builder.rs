//! OCCT BRepBlend_SurfRstLineBuilder (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_SurfRstLineBuilder.hxx (L65-197) + BRepBlend_SurfRstLineBuilder.cxx
//! (L141-1776) + BRepBlend_SurfRstLineBuilder.lxx (L24-59).
//!
//! The class builds a BRepBlend_Line between a surface and a pcurve on the
//! other surface from an approached starting solution (hxx L56-64).  The
//! OCCT_DEBUG trace statics (BBPP / tracederiv / Drawsect, cxx L40-133) are
//! debug-only and are not translated.
//!
//! Architecture mappings and pending boundaries (Stage 1e third batch):
//!   - `occ::handle<Adaptor3d_Surface>` -> `&BRepAdaptorSurface` (the
//!     consumer-pinned carrier, chfi3d_builder_0/chfi3d_builder_2).
//!   - `occ::handle<Adaptor3d_TopolTool>` -> `&BRepTopAdaptorTopolTool` (the
//!     rcad stub carrier).  OCCT mutates the tool through the
//!     Init/More/Next/Initialize protocol; the rcad carrier is immutable, so
//!     the protocol is an explicit iterator value ([`DomainTool`], shared
//!     with brep_blend_rst_rst_line_builder).
//!   - `occ::handle<Adaptor2d_Curve2d>` -> [`RstArc`] (the pcurve value plus
//!     the owning edge/face identity — the BRepAdaptor_Curve2d view).  The
//!     `rst` member arrives from the consumer as a bare pcurve (no edge
//!     identity): the vertex scans over `domain2->Initialize(rst)`
//!     (cxx L1171-1187, L1265-1281) therefore see no vertices — pending the
//!     edge-identity carrier on the consumer side.
//!   - `occ::handle<Adaptor3d_HVertex>` -> [`DomainVertex`] (the
//!     Adaptor3d_HVertex view + the topological vertex identity used by
//!     BRepTopAdaptor_HVertex::IsSame, BRepTopAdaptor_HVertex.cxx L188-192).
//!   - `BRepBlend_BlendTool` -> the [`blend_tool_*`] free functions
//!     (BRepBlend_BlendTool.cxx L38-142 + .lxx L19-58).
//!   - `math_Vector sol(1, 3)` -> `Vec<f64>` (1-based OCCT access ->
//!     0-based indexing).
//!   - `math_FunctionSetRoot` -> `rcad_kernel::math::FunctionSetRoot` (the
//!     4-argument OCCT Perform maps to perform(..., false)).
//!   - PENDING: BRepTopAdaptor_HVertex::Resolution
//!     (BRepTopAdaptor_HVertex.cxx L48-185 — the UV-refinement chain needs
//!     BRepAdaptor_Surface::U/VResolution + pcurve Resolution) — the bridge
//!     stores the BRep_Tool::Tolerance(V) seed instead.

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetRoot;
use rcad_kernel::topods::Orientation;

use crate::geomalgo::int_patch::transitions::{
    make_transition, Transition as IntSurfTransition,
};

use super::brep_blend::BlendStatus;
use super::brep_blend_extremity::BRepBlendExtremity;
use super::brep_blend_func_inv::{BlendFuncInv, BlendSurfCurvFuncInv, BlendSurfPointFuncInv};
use super::brep_blend_line::{BRepBlendLine, IntSurfTypeTrans};
use super::brep_blend_point::BlendPoint;
use super::brep_blend_surf_rst_function::BlendSurfRstFunction;
use super::brep_blend_surf_rst_line_builder_b::{
    blend_tool_inters, blend_tool_parameter, blend_tool_project, blend_tool_tolerance, DomainTool,
    DomainVertex, RstArc,
};
use super::chfi3d::topabs_reverse;
use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_2::{BRepTopAdaptorTopolTool, TopAbsState};

/// OCCT ConvOrToTra(O) (cxx L1601-1608) — the orientation to transition
/// conversion.
fn conv_or_to_tra(o: Orientation) -> IntSurfTypeTrans {
    if o == Orientation::Forward {
        IntSurfTypeTrans::In
    } else {
        IntSurfTypeTrans::Out
    }
}

/// OCCT BRepBlend_SurfRstLineBuilder — this class processes data resulting
/// from Blend_CSWalking taking in consideration the Surface supporting the
/// curve to detect the breakpoint (hxx L37).  OCCT inheritance: standalone
/// class (no base).
pub struct BRepBlendSurfRstLineBuilder<'a> {
    done: bool,
    line: BRepBlendLine,
    sol: Vec<f64>,
    surf1: &'a BRepAdaptorSurface,
    domain1: &'a BRepTopAdaptorTopolTool,
    surf2: &'a BRepAdaptorSurface,
    rst: RstArc,
    domain2: &'a BRepTopAdaptorTopolTool,
    tolpoint3d: f64,
    tolpoint2d: f64,
    tolgui: f64,
    pasmax: f64,
    fleche: f64,
    param: f64,
    previous_p: BlendPoint,
    rebrou: bool,
    iscomplete: bool,
    comptra: bool,
    sens: f64,
    decrochdeb: bool,
    decrochfin: bool,
}

impl<'a> BRepBlendSurfRstLineBuilder<'a> {
    /// OCCT BRepBlend_SurfRstLineBuilder(Surf1, Domain1, Surf2, Rst, Domain2)
    /// (cxx L195-221).
    pub fn new(
        surf1: &'a BRepAdaptorSurface,
        domain1: &'a BRepTopAdaptorTopolTool,
        surf2: &'a BRepAdaptorSurface,
        rst: &'a Curve2d,
        domain2: &'a BRepTopAdaptorTopolTool,
    ) -> Self {
        BRepBlendSurfRstLineBuilder {
            done: false,
            line: BRepBlendLine::new(),
            sol: vec![0.0; 3], // OCCT: math_Vector sol(1, 3)
            surf1,
            domain1,
            surf2,
            rst: RstArc::bare(rst),
            domain2,
            tolpoint3d: 0.0,
            tolpoint2d: 0.0,
            tolgui: 0.0,
            pasmax: 0.0,
            fleche: 0.0,
            param: 0.0,
            previous_p: BlendPoint::new(),
            rebrou: false,
            iscomplete: false,
            comptra: false,
            sens: 0.0,
            decrochdeb: false,
            decrochfin: false,
        }
    }

    /// OCCT IsDone() (lxx L24-28).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Line() (lxx L34-41) — throws StdFail_NotDone when not done.
    pub fn line(&self) -> &BRepBlendLine {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_SurfRstLineBuilder::Line");
        }
        &self.line
    }

    /// OCCT DecrochStart() (lxx L47-50).
    pub fn decroch_start(&self) -> bool {
        self.decrochdeb
    }

    /// OCCT DecrochEnd() (lxx L55-58).
    pub fn decroch_end(&self) -> bool {
        self.decrochfin
    }

    /// OCCT ArcToRecadre(Sol, PrevIndex, lastpt2d, pt2d, ponarc)
    /// (cxx L141-191) — find a suitable arc; PrevIndex is used to reject an
    /// already tested arc.
    #[allow(clippy::too_many_arguments)]
    pub fn arc_to_recadre(
        &mut self,
        the_sol: &[f64],
        prev_index: i32,
        lastpt2d: &mut DVec2,
        pt2d: &mut DVec2,
        ponarc: &mut f64,
    ) -> i32 {
        let mut index_sol = 0;
        let mut nbarc = 0;
        let byinter = self.line.nb_points() != 0;
        let mut okinter;
        let mut distmin = f64::MAX; // OCCT: RealLast()
        let (mut uprev, mut vprev) = (0.0, 0.0);
        let mut prm;
        let mut dist;

        if byinter {
            (uprev, vprev) = self.previous_p.parameters_on_s();
        }
        *pt2d = DVec2::new(the_sol[0], the_sol[1]);
        *lastpt2d = DVec2::new(uprev, vprev);
        let mut domain1 = DomainTool::face_domain(self.domain1);
        domain1.init();

        while domain1.more() {
            nbarc += 1;
            let arc = domain1.value();
            let mut ok = false;
            okinter = false;
            if byinter {
                let mut prm_loc = 0.0;
                let mut dist_loc = 0.0;
                ok = blend_tool_inters(*pt2d, *lastpt2d, &arc, &mut prm_loc, &mut dist_loc);
                prm = prm_loc;
                dist = dist_loc;
                okinter = ok;
            } else {
                prm = 0.0;
                dist = 0.0;
            }
            if !ok {
                let mut prm_loc = prm;
                let mut dist_loc = dist;
                ok = blend_tool_project(*pt2d, &arc, &mut prm_loc, &mut dist_loc);
                prm = prm_loc;
                dist = dist_loc;
            }

            if ok && nbarc != prev_index {
                if dist < distmin || okinter {
                    distmin = dist;
                    *ponarc = prm;
                    index_sol = nbarc;
                    if okinter && prev_index == 0 {
                        break;
                    }
                }
            }
            domain1.next();
        }
        index_sol
    }

    /// OCCT Perform(Func, Finv, FinvP, FinvC, Pdep, Pmax, MaxStep, Tol3d,
    /// Tol2d, TolGuide, ParDep, Fleche, Appro) (cxx L225-324).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        pdep: f64,
        pmax: f64,
        max_step: f64,
        tol3d: f64,
        tol2d: f64,
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
        self.tolpoint2d = tol2d;
        self.tolgui = tolguide.abs();
        self.fleche = fleche.abs();
        self.rebrou = false;
        self.pasmax = max_step.abs();

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };

        self.param = pdep;
        func.set_param(self.param);

        if appro {
            let mut siturst = TopAbsState::Unknown;
            let mut situs = TopAbsState::Unknown;
            let mut decroch = false;
            let mut tolerance = vec![0.0; 3];
            let mut infbound = vec![0.0; 3];
            let mut supbound = vec![0.0; 3];
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);
            let mut rsnld = FunctionSetRoot::new(func, &tolerance, 30);

            rsnld.perform(func, par_dep, &infbound, &supbound, false);

            if !rsnld.is_done() {
                return;
            }
            self.sol = rsnld.root();
            if !self.check_inside(func, &mut siturst, &mut situs, &mut decroch) {
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
        let (u, v);
        (u, v) = self.previous_p.parameters_on_s();
        //  W = previousP.ParameterOnC();

        let mut ptf1 = BRepBlendExtremity::new_on_surface(
            self.previous_p.point_on_s(),
            u,
            v,
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        let mut ptf2 = BRepBlendExtremity::new_on_surface(
            self.previous_p.point_on_c(),
            u,
            v,
            self.previous_p.parameter(),
            self.tolpoint3d,
        );
        if !self.previous_p.is_tangency_point() {
            ptf1.set_tangent(self.previous_p.tangent_on_s());
            ptf2.set_tangent(self.previous_p.tangent_on_c());
        }
        if self.sens > 0.0 {
            self.line.set_start_points(&ptf1, &ptf2);
        } else {
            self.line.set_end_points(&ptf1, &ptf2);
        }
        self.internal_perform(func, finv, finv_p, finv_c, pmax);
        self.done = true;
    }

    /// OCCT PerformFirstSection(Func, Finv, FinvP, FinvC, Pdep, Pmax, ParDep,
    /// Tol3d, Tol2d, TolGuide, RecRst, RecP, RecS, Psol, ParSol)
    /// (cxx L328-494).
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn perform_first_section(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        pdep: f64,
        pmax: f64,
        par_dep: &[f64],
        tol3d: f64,
        tol2d: f64,
        tolguide: f64,
        rec_rst: bool,
        rec_p: bool,
        rec_s: bool,
        psol: &mut f64,
        par_sol: &mut [f64],
    ) -> bool {
        self.done = false;
        self.iscomplete = false;
        self.comptra = false;
        self.line = BRepBlendLine::new();
        self.tolpoint3d = tol3d;
        self.tolpoint2d = tol2d;
        self.tolgui = tolguide.abs();
        self.rebrou = false;

        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };
        let mut state = BlendStatus::OnRst12;
        let mut trst = 0.0;
        let recads;
        let mut recadrst;
        let recadp;
        let (mut wp, mut wrst, mut ws);
        let (mut u, mut v) = (0.0, 0.0);
        let mut infbound = vec![0.0; 3];
        let mut supbound = vec![0.0; 3];
        let mut tolerance = vec![0.0; 3];
        let mut solinvp = vec![0.0; 3];
        let mut solinvrst = vec![0.0; 4];
        let mut solinvs = vec![0.0; 3];
        let mut vtxp = DomainVertex::empty();
        let mut vtxrst = DomainVertex::empty();
        let mut vtxs = DomainVertex::empty();
        let mut is_vtxp = false;
        let mut is_vtxrst = false;
        let mut is_vtxs = false;
        let mut arc = RstArc::empty();
        wp = pmax;
        wrst = pmax;
        ws = pmax;
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

        recads = rec_s
            && self.recadre_finv_c(finv_c, &mut solinvs, &mut arc, &mut is_vtxs, &mut vtxs);
        if recads {
            ws = solinvs[0];
        }
        recadp = rec_p && self.recadre_finv_p(finv_p, &mut solinvp, &mut is_vtxp, &mut vtxp);
        if recadp {
            wp = solinvp[0];
        }
        recadrst =
            rec_rst && self.recadre_func_finv(func, finv, &mut solinvrst, &mut is_vtxrst, &mut vtxrst);
        if recadrst {
            wrst = solinvrst[1];
        }
        if !recads && !recadp && !recadrst {
            return false;
        }
        if recadp && recadrst {
            if self.sens * (wrst - wp) > self.tolgui {
                // first one leaves the domain
                wrst = wp;
                u = solinvp[1];
                v = solinvp[2];
                trst = blend_tool_parameter(&vtxp, &self.rst);
                is_vtxrst = is_vtxp;
                vtxrst = vtxp.clone();
            } else {
                u = solinvrst[2];
                v = solinvrst[3];
                trst = solinvrst[0];
            }
        } else if recadp {
            wrst = wp;
            u = solinvp[1];
            v = solinvp[2];
            trst = blend_tool_parameter(&vtxp, &self.rst);
            is_vtxrst = is_vtxp;
            vtxrst = vtxp.clone();
            recadrst = true;
        } else if recadrst {
            u = solinvrst[2];
            v = solinvrst[3];
            trst = solinvrst[0];
        }
        if recads && recadrst {
            if (ws - wrst).abs() < self.tolgui {
                state = BlendStatus::OnRst12;
                self.param = 0.5 * (ws + wrst);
                self.sol[0] = u;
                self.sol[1] = v;
                self.sol[2] = solinvs[1];
            } else if self.sens * (ws - wrst) < 0.0 {
                // ground on surf
                state = BlendStatus::OnRst1;
                self.param = ws;
                let p = arc.value(solinvs[2]);
                u = p.x;
                v = p.y;
                self.sol[0] = u;
                self.sol[1] = v;
                self.sol[2] = solinvs[1];
            } else {
                // ground on rst
                state = BlendStatus::OnRst2;
                self.param = wrst;
                self.sol[0] = u;
                self.sol[1] = v;
                self.sol[2] = trst;
            }
            func.set_param(self.param);
        } else if recads {
            // ground on surf
            state = BlendStatus::OnRst1;
            self.param = ws;
            let p = arc.value(solinvs[2]);
            u = p.x;
            v = p.y;
            self.sol[0] = u;
            self.sol[1] = v;
            self.sol[2] = solinvs[1];
            func.set_param(self.param);
        } else if recadrst {
            // ground on rst
            state = BlendStatus::OnRst2;
            self.param = wrst;
            self.sol[0] = u;
            self.sol[1] = v;
            self.sol[2] = trst;
            func.set_param(self.param);
        }
        let state = self.test_arret(func, false, state);
        *psol = self.param;
        par_sol.copy_from_slice(&self.sol);
        let _ = state;
        true
    }

    /// OCCT Complete(Func, Finv, FinvP, FinvC, Pmin) (cxx L498-528).
    pub fn complete(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        pmin: f64,
    ) -> bool {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_SurfRstLineBuilder::Complete");
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
        let (u, v) = self.previous_p.parameters_on_s();
        self.sol[0] = u;
        self.sol[1] = v;
        self.sol[2] = self.previous_p.parameter_on_c();

        self.internal_perform(func, finv, finv_p, finv_c, pmin);
        self.iscomplete = true;
        true
    }

    /// OCCT InternalPerform(Func, Finv, FinvP, FinvC, Bound) (cxx L532-995).
    #[allow(unused_assignments)]
    fn internal_perform(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        bound: f64,
    ) {
        let mut stepw = self.pasmax;
        let nbp = self.line.nb_points();
        if nbp >= 2 {
            // The last step is reproduced if it is not too small.
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
        let mut situonc = TopAbsState::Unknown;
        let mut situons = TopAbsState::Unknown;
        let mut decroch = false;
        let mut arrive;
        let mut recadp;
        let mut recadrst;
        let mut recads;
        let mut echecrecad;
        let mut wp;
        let mut wrst;
        let mut ws;
        let mut u = 0.0;
        let mut v = 0.0;
        let mut trst = 0.0;
        let mut infbound = vec![0.0; 3];
        let mut supbound = vec![0.0; 3];
        let mut parinit = vec![0.0; 3];
        let mut tolerance = vec![0.0; 3];
        let mut solinvp = vec![0.0; 3];
        let mut solinvrst = vec![0.0; 4];
        let mut solinvs = vec![0.0; 3];
        let mut vtxp = DomainVertex::empty();
        let mut vtxrst = DomainVertex::empty();
        let mut vtxs = DomainVertex::empty();
        let mut is_vtxp = false;
        let mut is_vtxrst = false;
        let mut is_vtxs = false;
        let mut extrst = BRepBlendExtremity::new();
        let mut exts = BRepBlendExtremity::new();
        let mut arc = RstArc::empty();

        // IntSurf_Transition Tline,Tarc;

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
                if !self.check_inside(func, &mut situonc, &mut situons, &mut decroch)
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
                wp = bound;
                wrst = bound;
                ws = bound;
                recadp = false;
                recadrst = false;
                recads = false;
                echecrecad = false;
                if situons == TopAbsState::Out || situons == TopAbsState::On {
                    // pb inverse rst/rst
                    recads =
                        self.recadre_finv_c(finv_c, &mut solinvs, &mut arc, &mut is_vtxs, &mut vtxs);
                    if recads {
                        ws = solinvs[0];
                        // It is necessary to reevaluate the deviation (BUC60360)
                        let mut t = DVec3::ZERO;
                        let mut n = DVec3::ZERO;
                        func.set_param(ws);
                        let p = arc.value(solinvs[2]);
                        u = p.x;
                        v = p.y;
                        self.sol[0] = u;
                        self.sol[1] = v;
                        self.sol[2] = solinvs[1];
                        decroch = func.decroch(&self.sol, &mut n, &mut t);
                    } else {
                        echecrecad = true;
                    }
                }
                if situonc == TopAbsState::Out || situonc == TopAbsState::On {
                    // pb inverse point/surf
                    recadp = self.recadre_finv_p(finv_p, &mut solinvp, &mut is_vtxp, &mut vtxp);
                    if recadp {
                        wp = solinvp[0];
                    } else {
                        echecrecad = true;
                    }
                }
                if decroch {
                    // pb inverse rst/surf
                    recadrst =
                        self.recadre_func_finv(func, finv, &mut solinvrst, &mut is_vtxrst, &mut vtxrst);
                    if recadrst {
                        wrst = solinvrst[1];
                    } else {
                        echecrecad = true;
                    }
                }
                decroch = false;
                if recadp || recads || recadrst {
                    echecrecad = false;
                }
                if !echecrecad {
                    if recadp && recadrst {
                        if self.sens * (wrst - wp) > self.tolgui {
                            // first one leaves the domain
                            wrst = wp;
                            u = solinvp[1];
                            v = solinvp[2];
                            trst = blend_tool_parameter(&vtxp, &self.rst);
                            is_vtxrst = is_vtxp;
                            vtxrst = vtxp.clone();
                        } else {
                            decroch = true;
                            u = solinvrst[2];
                            v = solinvrst[3];
                            trst = solinvrst[0];
                        }
                    } else if recadp {
                        wrst = wp;
                        u = solinvp[1];
                        v = solinvp[2];
                        trst = blend_tool_parameter(&vtxp, &self.rst);
                        is_vtxrst = is_vtxp;
                        vtxrst = vtxp.clone();
                        recadrst = true;
                    } else if recadrst {
                        decroch = true;
                        u = solinvrst[2];
                        v = solinvrst[3];
                        trst = solinvrst[0];
                    }
                    if recads && recadrst {
                        if (ws - wrst).abs() < self.tolgui {
                            state = BlendStatus::OnRst12;
                            self.param = 0.5 * (ws + wrst);
                            self.sol[0] = u;
                            self.sol[1] = v;
                            self.sol[2] = solinvs[2];
                        } else if self.sens * (ws - wrst) < 0.0 {
                            // ground on surf
                            decroch = false;
                            state = BlendStatus::OnRst1;
                            self.param = ws;
                            let p = arc.value(solinvs[2]);
                            u = p.x;
                            v = p.y;
                            self.sol[0] = u;
                            self.sol[1] = v;
                            self.sol[2] = solinvs[1];
                        } else {
                            // ground on rst
                            state = BlendStatus::OnRst2;
                            self.param = wrst;
                            self.sol[0] = u;
                            self.sol[1] = v;
                            self.sol[2] = trst;
                        }
                        func.set_param(self.param);
                    } else if recads {
                        // ground on surf
                        state = BlendStatus::OnRst1;
                        self.param = ws;
                        let p = arc.value(solinvs[2]);
                        u = p.x;
                        v = p.y;
                        self.sol[0] = u;
                        self.sol[1] = v;
                        self.sol[2] = solinvs[1];
                        func.set_param(self.param);
                    } else if recadrst {
                        // ground on rst
                        state = BlendStatus::OnRst2;
                        self.param = wrst;
                        self.sol[0] = u;
                        self.sol[1] = v;
                        self.sol[2] = trst;
                        func.set_param(self.param);
                    } else {
                        state = BlendStatus::Ok;
                    }
                    state = self.test_arret(func, true, state);
                } else {
                    // Failed reframing. Leave with PointsConfondus
                    state = BlendStatus::SamePoints;
                }
            }

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
                        exts.set_value(
                            self.previous_p.point_on_s(),
                            self.sol[0],
                            self.sol[1],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        let rst = self.rst.clone();
                        self.make_extremity(
                            &mut extrst,
                            false,
                            &rst,
                            self.sol[2],
                            is_vtxrst,
                            &vtxrst,
                        );
                        // Indicate end on Bound.
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
                        let (pu, pv) = self.previous_p.parameters_on_s();
                        u = pu;
                        v = pv;
                        exts.set_value(
                            self.previous_p.point_on_s(),
                            u,
                            v,
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        extrst.set_value_on_curve(
                            self.previous_p.point_on_c(),
                            self.previous_p.parameter_on_c(),
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        arrive = true;
                        if self.line.nb_points() >= 2 {
                            // Indicate that one stops during the processing
                        }
                    } else {
                        self.param = parprec + self.sens * stepw; // no risk to exceed Bound.
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
                        exts.set_value(
                            self.previous_p.point_on_s(),
                            self.sol[0],
                            self.sol[1],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        let rst = self.rst.clone();
                        self.make_extremity(
                            &mut extrst,
                            false,
                            &rst,
                            self.sol[2],
                            is_vtxrst,
                            &vtxrst,
                        );
                        // Indicate end on Bound.
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
                    self.make_extremity(&mut exts, true, &arc, solinvs[2], is_vtxs, &vtxs);
                    let rst = self.rst.clone();
                    self.make_extremity(
                        &mut extrst,
                        false,
                        &rst,
                        self.sol[2],
                        is_vtxrst,
                        &vtxrst,
                    );
                    arrive = true;
                }

                BlendStatus::OnRst2 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }
                    exts.set_value(
                        self.previous_p.point_on_s(),
                        self.sol[0],
                        self.sol[1],
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    let rst = self.rst.clone();
                    self.make_extremity(
                        &mut extrst,
                        false,
                        &rst,
                        self.sol[2],
                        is_vtxrst,
                        &vtxrst,
                    );
                    arrive = true;
                }

                BlendStatus::OnRst12 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }
                    self.make_extremity(&mut exts, true, &arc, solinvs[0], is_vtxs, &vtxs);
                    let rst = self.rst.clone();
                    self.make_extremity(
                        &mut extrst,
                        false,
                        &rst,
                        self.sol[2],
                        is_vtxrst,
                        &vtxrst,
                    );
                    arrive = true;
                }

                BlendStatus::SamePoints => {
                    // Stop
                    let (pu, pv) = self.previous_p.parameters_on_s();
                    u = pu;
                    v = pv;
                    exts.set_value(
                        self.previous_p.point_on_s(),
                        u,
                        v,
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    extrst.set_value_on_curve(
                        self.previous_p.point_on_c(),
                        self.previous_p.parameter_on_c(),
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    arrive = true;
                }

                BlendStatus::Backward => {}
            }
            if arrive {
                if self.sens > 0.0 {
                    self.line.set_end_points(&exts, &extrst);
                    self.decrochfin = decroch;
                } else {
                    self.line.set_start_points(&exts, &extrst);
                    self.decrochdeb = decroch;
                }
            }
        }
    }

    /// OCCT Recadre(FinvC, Solinv, Arc, IsVtx, Vtx) (cxx L1002-1134) —
    /// reframe section Surface / Restriction on a domain arc.
    fn recadre_finv_c(
        &mut self,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        solinv: &mut [f64],
        arc: &mut RstArc,
        is_vtx: &mut bool,
        vtx: &mut DomainVertex,
    ) -> bool {
        let mut recadre = false;

        let mut pt2d = DVec2::ZERO;
        let mut lastpt2d = DVec2::ZERO;
        let mut pmin = 0.0;
        let mut nbarc;

        let sol = self.sol.clone();
        let mut index_sol = self.arc_to_recadre(&sol, 0, &mut lastpt2d, &mut pt2d, &mut pmin);

        *is_vtx = false;
        if index_sol == 0 {
            return false;
        }

        let mut domain1 = DomainTool::face_domain(self.domain1);
        domain1.init();
        nbarc = 1;
        while nbarc < index_sol {
            nbarc += 1;
            domain1.next();
        }
        *arc = domain1.value();

        finv_c.set_rst(&arc.curve);

        let mut toler = vec![0.0; 3];
        let mut infb = vec![0.0; 3];
        let mut supb = vec![0.0; 3];
        // use reduced Tol argument value to pass testcase
        // blend complex A6 with scale factor of model 0.1 (base scale = 1000)
        // So, here we using 1.0e-5 rather than 1.0e-4 of tolerance of point in 3d
        finv_c.get_tolerance(&mut toler, 0.1 * self.tolpoint3d);
        finv_c.get_bounds(&mut infb, &mut supb);
        solinv[0] = self.param;
        solinv[1] = self.sol[2];
        solinv[2] = pmin;

        let mut rsnld = FunctionSetRoot::new(finv_c, &toler, 30);
        rsnld.perform(finv_c, solinv, &infb, &supb, false);

        if !rsnld.is_done() {
            // OCCT prints "SurfRstLineBuilder : RSNLD not done" under OCCT_DEBUG.
        } else {
            // It is necessary to check the value of the function
            let root = rsnld.root();
            solinv.copy_from_slice(&root);
            recadre = finv_c.is_solution(solinv, self.tolpoint3d);
        }

        // In case of fail, it is checked if another arc
        // can be useful (case of output at the proximity of a vertex)
        if !recadre {
            let sol = self.sol.clone();
            index_sol = self.arc_to_recadre(&sol, index_sol, &mut lastpt2d, &mut pt2d, &mut pmin);
            if index_sol == 0 {
                return false; // No other solution
            }

            let mut domain1 = DomainTool::face_domain(self.domain1);
            domain1.init();
            nbarc = 1;
            while nbarc < index_sol {
                nbarc += 1;
                domain1.next();
            }

            *arc = domain1.value();
            finv_c.set_rst(&arc.curve);

            finv_c.get_tolerance(&mut toler, self.tolpoint3d);
            finv_c.get_bounds(&mut infb, &mut supb);

            solinv[2] = pmin;

            let mut a_rsnld = FunctionSetRoot::new(finv_c, &toler, 30);
            a_rsnld.perform(finv_c, solinv, &infb, &supb, false);

            if !a_rsnld.is_done() {
                // OCCT prints "SurfRstLineBuilder : RSNLD not done" under OCCT_DEBUG.
            } else {
                // It is necessary to check the value of the function
                let root = a_rsnld.root();
                solinv.copy_from_slice(&root);
                recadre = finv_c.is_solution(solinv, self.tolpoint3d);
            }
        }

        if recadre {
            let w = solinv[1];
            if w < self.rst.first_parameter() - toler[1]
                || w > self.rst.last_parameter() + toler[1]
            {
                return false;
            }
            let mut domain1 = DomainTool::face_domain(self.domain1);
            domain1.initialize_arc(arc);
            domain1.init_vertex_iterator();
            *is_vtx = !domain1.more_vertex();
            while !*is_vtx {
                *vtx = domain1.vertex();
                if (blend_tool_parameter(vtx, arc) - solinv[2]).abs()
                    <= blend_tool_tolerance(vtx, arc)
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

    /// OCCT Recadre(Func, Finv, Solinv, IsVtx, Vtx) (cxx L1138-1217) —
    /// reframe section Surface / Restriction by inverse function.
    fn recadre_func_finv(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        solinv: &mut [f64],
        is_vtx: &mut bool,
        vtx: &mut DomainVertex,
    ) -> bool {
        let mut toler = vec![0.0; 4];
        let mut infb = vec![0.0; 4];
        let mut supb = vec![0.0; 4];
        finv.get_tolerance(&mut toler, self.tolpoint3d);
        finv.get_bounds(&mut infb, &mut supb);
        solinv[0] = self.sol[2];
        solinv[1] = self.param;
        solinv[2] = self.sol[0];
        solinv[3] = self.sol[1];

        let mut rsnld = FunctionSetRoot::new(finv, &toler, 30);
        rsnld.perform(finv, solinv, &infb, &supb, false);
        if !rsnld.is_done() {
            // OCCT prints "SurfRstLineBuilder :RSNLD not done" under OCCT_DEBUG.
            return false;
        }
        let root = rsnld.root();
        solinv.copy_from_slice(&root);

        if finv.is_solution(solinv, self.tolpoint3d) {
            let p2d = DVec2::new(solinv[2], solinv[3]);
            let situ = self.domain1.classify(p2d, self.tolpoint2d, false);
            if situ != TopAbsState::In && situ != TopAbsState::On {
                return false;
            }
            let mut domain2 = DomainTool::single(&self.rst);
            domain2.init_vertex_iterator();
            *is_vtx = !domain2.more_vertex();
            while !*is_vtx {
                *vtx = domain2.vertex();
                if (blend_tool_parameter(vtx, &self.rst) - solinv[0]).abs()
                    <= blend_tool_tolerance(vtx, &self.rst)
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
            // The section is recalculated by direct resolution, otherwise
            // incoherences between the parameter and the ground caused by yawn
            // are returned.

            let mut infbound = vec![0.0; 3];
            let mut supbound = vec![0.0; 3];
            let mut parinit = vec![0.0; 3];
            let mut tolerance = vec![0.0; 3];
            func.get_tolerance(&mut tolerance, self.tolpoint3d);
            func.get_bounds(&mut infbound, &mut supbound);

            let mut rsnld2 = FunctionSetRoot::new(func, &tolerance, 30);
            parinit[0] = solinv[2];
            parinit[1] = solinv[3];
            parinit[2] = solinv[0];
            func.set_param(solinv[1]);
            rsnld2.perform(func, &parinit, &infbound, &supbound, false);
            if !rsnld2.is_done() {
                return false;
            }
            let root = rsnld2.root();
            parinit.copy_from_slice(&root);
            solinv[2] = parinit[0];
            solinv[3] = parinit[1];
            solinv[0] = parinit[2];
            return true;
        }
        false
    }

    /// OCCT Recadre(FinvP, Solinv, IsVtx, Vtx) (cxx L1221-1289) — reframe
    /// section Surface / Restriction on an extremity point of the rst.
    fn recadre_finv_p(
        &mut self,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        solinv: &mut [f64],
        is_vtx: &mut bool,
        vtx: &mut DomainVertex,
    ) -> bool {
        // Le point.
        let firstrst = self.rst.first_parameter();
        let lastrst = self.rst.last_parameter();
        let mut wpoint = firstrst;
        if (self.sol[2] - firstrst) > (lastrst - self.sol[2]) {
            wpoint = lastrst;
        }
        let p2drst = self.rst.value(wpoint);
        let thepoint = self.surf2.surface.point_at(p2drst.x, p2drst.y);

        finv_p.set_point(thepoint);
        let mut toler = vec![0.0; 3];
        let mut infb = vec![0.0; 3];
        let mut supb = vec![0.0; 3];
        finv_p.get_tolerance(&mut toler, self.tolpoint3d);
        finv_p.get_bounds(&mut infb, &mut supb);
        solinv[0] = self.param;
        solinv[1] = self.sol[0];
        solinv[2] = self.sol[1];

        let mut rsnld = FunctionSetRoot::new(finv_p, &toler, 30);
        rsnld.perform(finv_p, solinv, &infb, &supb, false);
        if !rsnld.is_done() {
            // OCCT prints "SurfRstLineBuilder :RSNLD not done" under OCCT_DEBUG.
            return false;
        }
        let root = rsnld.root();
        solinv.copy_from_slice(&root);

        if finv_p.is_solution(solinv, self.tolpoint3d) {
            let p2d = DVec2::new(solinv[1], solinv[2]);
            let situ = self.domain1.classify(p2d, self.tolpoint2d, false);
            if situ != TopAbsState::In && situ != TopAbsState::On {
                return false;
            }
            let mut domain2 = DomainTool::single(&self.rst);
            domain2.init_vertex_iterator();
            *is_vtx = !domain2.more_vertex();
            while !*is_vtx {
                *vtx = domain2.vertex();
                if (blend_tool_parameter(vtx, &self.rst) - wpoint).abs()
                    <= blend_tool_tolerance(vtx, &self.rst)
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

    /// OCCT Transition(OnFirst, Arc, Param, TLine, TArc) (cxx L1293-1355).
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
                tgline = self.previous_p.tangent_on_s1();
            } else {
                tgline = self.previous_p.point_on_s() - prevprev.point_on_s();
            }
        } else {
            let (_pbid, du, dv) = self.surf2.surface.derivatives(p2d.x, p2d.y);
            d1u = du;
            d1v = dv;
            if !computetranstionaveclacorde {
                tgline = self.previous_p.tangent_on_s2();
            } else {
                tgline = self.previous_p.point_on_c() - prevprev.point_on_c();
            }
        }

        let tgrst = d1u * dp2d.x + d1v * dp2d.y;
        let normale = d1u.cross(d1v);

        make_transition(tgline, tgrst, normale, t_line, t_arc);
    }

    /// OCCT MakeExtremity(Extrem, OnFirst, Arc, Param, IsVtx, Vtx)
    /// (cxx L1359-1423).
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
            extrem.set_value(
                self.previous_p.point_on_s(),
                self.sol[0],
                self.sol[1],
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_s());
            }
            DomainTool::face_domain(self.domain1)
        } else {
            extrem.set_value_on_curve(
                self.previous_p.point_on_c(),
                self.sol[2],
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_c());
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

    /// OCCT CheckDeflectionOnSurf(CurPoint) (cxx L1427-1511) — controls 3d
    /// of Blend_CSWalking.
    fn check_deflection_on_surf(&mut self, cur_point: &BlendPoint) -> BlendStatus {
        // rule by tests in U4 corresponds to 11.478 d
        let cos_ref_3d = 0.98;
        let mut cosi;
        let mut cosi2;
        let curpointistangent = cur_point.is_tangency_point();
        let prevpointistangent = self.previous_p.is_tangency_point();

        let psurf = cur_point.point_on_s();
        let tgsurf;
        if !curpointistangent {
            tgsurf = cur_point.tangent_on_s();
        } else {
            tgsurf = DVec3::ZERO;
        }
        let prevp = self.previous_p.point_on_s();
        let prevtg;
        if !prevpointistangent {
            prevtg = self.previous_p.tangent_on_s();
        } else {
            prevtg = DVec3::ZERO;
        }
        let norme;
        let mut prevnorme = 0.0;
        let corde = psurf - prevp;
        norme = corde.length_squared();
        //  if(!curpointistangent) curNorme = Tgsurf.SquareMagnitude();
        if !prevpointistangent {
            prevnorme = prevtg.length_squared();
        }

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // it can be necessary to force same point
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
                // not too great :
                return BlendStatus::StepTooLarge;
            }
        }
        BlendStatus::Ok
    }

    /// OCCT CheckDeflectionOnRst(CurPoint) (cxx L1515-1599) — controls 3D
    /// of Blend_CSWalking.
    fn check_deflection_on_rst(&mut self, cur_point: &BlendPoint) -> BlendStatus {
        // rule by tests in U4 corresponds to 11.478 d
        let cos_ref_3d = 0.98;
        let mut cosi;
        let mut cosi2;
        let curpointistangent = cur_point.is_tangency_point();
        let prevpointistangent = self.previous_p.is_tangency_point();

        let psurf = cur_point.point_on_c();
        let tgsurf;
        if !curpointistangent {
            tgsurf = cur_point.tangent_on_c();
        } else {
            tgsurf = DVec3::ZERO;
        }
        let prevp = self.previous_p.point_on_c();
        let prevtg;
        if !prevpointistangent {
            prevtg = self.previous_p.tangent_on_c();
        } else {
            prevtg = DVec3::ZERO;
        }
        let norme;
        let mut prevnorme = 0.0;
        let corde = psurf - prevp;
        norme = corde.length_squared();
        //  if(!curpointistangent) curNorme = Tgsurf.SquareMagnitude();
        if !prevpointistangent {
            prevnorme = prevtg.length_squared();
        }

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // it can be necessary to force same point
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

    /// OCCT TestArret(Func, TestDeflection, State) (cxx L1612-1741).
    #[allow(unused_assignments)]
    fn test_arret(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        test_deflection: bool,
        state: BlendStatus,
    ) -> BlendStatus {
        let mut tgs = DVec3::ZERO;
        let mut tgrst = DVec3::ZERO;
        let mut tg2ds = DVec2::ZERO;
        let mut tg2drst = DVec2::ZERO;
        let mut tras = IntSurfTypeTrans::Undecided;
        let mut trarst = IntSurfTypeTrans::Undecided;

        if func.is_solution(&self.sol, self.tolpoint3d) {
            let curpointistangent = func.is_tangency_point();
            let pts = func.point_on_s();
            let ptrst = func.point_on_rst();
            // OCCT (cxx L1629): pt2drst = Func.Pnt2dOnRst();
            let pt2drst = func.pnt2d_on_rst();
            let curpoint;
            if curpointistangent {
                curpoint = BlendPoint::new_on_surface_curve_on_surface(
                    pts,
                    ptrst,
                    self.param,
                    self.sol[0],
                    self.sol[1],
                    pt2drst.x,
                    pt2drst.y,
                    self.sol[2],
                );
            } else {
                tgs = func.tangent_on_s();
                tgrst = func.tangent_on_rst();
                tg2ds = func.tangent_2d_on_s();
                tg2drst = func.tangent_2d_on_rst();

                curpoint = BlendPoint::new_on_surface_curve_on_surface_with_tangents(
                    pts,
                    ptrst,
                    self.param,
                    self.sol[0],
                    self.sol[1],
                    pt2drst.x,
                    pt2drst.y,
                    self.sol[2],
                    tgs,
                    tgrst,
                    tg2ds,
                    tg2drst,
                );
            }
            let mut state_s;
            let mut state_rst;
            if test_deflection {
                state_s = self.check_deflection_on_surf(&curpoint);
                state_rst = self.check_deflection_on_rst(&curpoint);
            } else {
                state_s = BlendStatus::Ok;
                state_rst = BlendStatus::Ok;
            }
            if state_s == BlendStatus::Backward {
                state_s = BlendStatus::StepTooLarge;
                self.rebrou = true;
            }
            if state_rst == BlendStatus::Backward {
                state_rst = BlendStatus::StepTooLarge;
                self.rebrou = true;
            }
            if state_s == BlendStatus::StepTooLarge || state_rst == BlendStatus::StepTooLarge {
                return BlendStatus::StepTooLarge;
            }

            if !self.comptra && !curpointistangent {
                let mut tgsecs = DVec3::ZERO;
                let mut nors = DVec3::ZERO;
                func.decroch(&self.sol, &mut nors, &mut tgsecs);
                nors = nors.normalize();
                let mut testra = tgsecs.dot(nors.cross(tgs));
                if testra.abs() > self.tolpoint3d {
                    if testra < 0.0 {
                        tras = IntSurfTypeTrans::In;
                    } else if testra > 0.0 {
                        tras = IntSurfTypeTrans::Out;
                    }
                    let (_p2drstref, tg2drstref) = self.rst.d1(self.sol[2]);
                    testra = tg2ds.dot(tg2drstref);
                    let or = self.rst.orientation();
                    if testra.abs() > 1.0e-8 {
                        if testra < 0.0 {
                            trarst = conv_or_to_tra(topabs_reverse(or));
                        } else if testra > 0.0 {
                            trarst = conv_or_to_tra(or);
                        }
                        self.comptra = true;
                        self.line.set_transitions(tras, trarst);
                    }
                }
            }
            if state_s == BlendStatus::Ok || state_rst == BlendStatus::Ok {
                self.previous_p = curpoint;
                return state;
            }
            if state_s == BlendStatus::StepTooSmall && state_rst == BlendStatus::StepTooSmall {
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

    /// OCCT CheckInside(Func, SituOnC, SituOnS, Decroch) (cxx L1745-1776).
    fn check_inside(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        situ_on_c: &mut TopAbsState,
        situ_on_s: &mut TopAbsState,
        decroch: &mut bool,
    ) -> bool {
        let mut tolerance = vec![0.0; 3];
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        // face pcurve.
        let w = self.sol[2];
        if w < self.rst.first_parameter() - tolerance[2]
            || w > self.rst.last_parameter() + tolerance[2]
        {
            *situ_on_c = TopAbsState::Out;
        } else if w > self.rst.first_parameter() && w < self.rst.last_parameter() {
            *situ_on_c = TopAbsState::In;
        } else {
            *situ_on_c = TopAbsState::On;
        }

        // face surface
        let p2d = DVec2::new(self.sol[0], self.sol[1]);
        *situ_on_s = self.domain1.classify(p2d, self.tolpoint2d, false);

        // lost contact
        let mut tgs = DVec3::ZERO;
        let mut nors = DVec3::ZERO;
        *decroch = func.decroch(&self.sol, &mut tgs, &mut nors);

        *situ_on_c == TopAbsState::In && *situ_on_s == TopAbsState::In && !*decroch
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT anchor: the constructor initial state (cxx L195-221).
    #[test]
    fn ctor_anchor() {
        let surf = BRepAdaptorSurface::empty();
        let tool = BRepTopAdaptorTopolTool::default();
        let rst = Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let b = BRepBlendSurfRstLineBuilder::new(&surf, &tool, &surf, &rst, &tool);
        assert!(!b.is_done());
        assert!(!b.decroch_start());
        assert!(!b.decroch_end());
    }

    /// OCCT anchor: ConvOrToTra (cxx L1601-1608).
    #[test]
    fn conv_or_to_tra_anchor() {
        assert_eq!(conv_or_to_tra(Orientation::Forward), IntSurfTypeTrans::In);
        assert_eq!(conv_or_to_tra(Orientation::Reversed), IntSurfTypeTrans::Out);
    }
}
