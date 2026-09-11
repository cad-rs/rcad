//! OCCT BRepBlend_Walking (TKFillet/BRepBlend) — translation of
//! BRepBlend_Walking.hxx (L36-212) and BRepBlend_Walking.cxx
//! (L83-614 here; the solver internals L616-2769 live in
//! [`super::brep_blend_walking_b`]).
//!
//! Architecture mappings: `occ::handle<Adaptor3d_Surface>` -> `&Surface3`;
//! `occ::handle<Adaptor3d_TopolTool>` -> `&BRepTopAdaptorTopolTool` (the
//! concrete rcad TopolTool, chfi3d_builder_2); `occ::handle<ChFiDS_ElSpine>`
//! -> `&Curve3` (the elspine guide curve); `math_Vector sol(1, 4)` ->
//! `Vec<f64>` (0-based storage for the 1-based OCCT indices); the NULL
//! `occ::handle<BRepBlend_Line>` -> `line` + `line_is_null` flag;
//! `NCollection_Sequence<Blend_Point>` -> `Vec<BlendPoint>`;
//! `Blend_Function&` / `Blend_FuncInv&` -> `&mut dyn BlendFunction` /
//! `&mut dyn BlendFuncInv` (trait upcasting feeds the math solver).
//!
//! Pending dependencies (plan 0.6, OCCT-named placeholders + failure path,
//! see brep_blend_walking_b): the restriction-arc iteration of the
//! placeholder `BRepTopAdaptorTopolTool`, BRepBlend_BlendTool::Inters
//! (Geom2dInt_GInter), the ChFiDS_ElSpine vertex list / saved parameters on
//! the rcad guide curve.

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Curve3, Surface3, SurfaceEval as _};
use rcad_kernel::math::cs_lib::{DerivativeStatus, NormalStatus};
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};
use crate::geomalgo::extrema_gen_ext_pc2d::EPCOfExtPC2d;
use crate::geomalgo::int_patch::transitions::Transition;

use super::brep_blend::BlendStatus;
use super::brep_blend_extremity::BRepBlendExtremity;
use super::brep_blend_func_inv::BlendFuncInv;
use super::brep_blend_function::BlendFunction;
use super::brep_blend_line::BRepBlendLine;
use super::brep_blend_point::BlendPoint;
use super::chfi3d_builder_2::{BRepTopAdaptorTopolTool, TopAbsState};
use crate::topalgo::adaptor3d::hvertex::HVertex;

/// OCCT BRepBlend_Walking — this class describes the asynchronous "walking"
/// of a blending surface on a pair of surfaces (BRepBlend_Walking.hxx L31).
pub struct BRepBlendWalking<'a> {
    pub(crate) previous_p: BlendPoint,
    pub(crate) line: BRepBlendLine,
    /// OCCT: `line.IsNull()` — the handle-null state of the line.
    pub(crate) line_is_null: bool,
    pub(crate) sol: Vec<f64>,
    pub(crate) jalons: Vec<BlendPoint>,
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) domain1: &'a BRepTopAdaptorTopolTool,
    pub(crate) domain2: &'a BRepTopAdaptorTopolTool,
    pub(crate) recdomain1: &'a BRepTopAdaptorTopolTool,
    pub(crate) recdomain2: &'a BRepTopAdaptorTopolTool,
    pub(crate) hguide: &'a Curve3,
    pub(crate) to_correct_on_rst1: bool,
    pub(crate) to_correct_on_rst2: bool,
    pub(crate) corrected_param: f64,
    pub(crate) tolpoint3d: f64,
    pub(crate) tolgui: f64,
    pub(crate) pasmax: f64,
    pub(crate) fleche: f64,
    pub(crate) param: f64,
    pub(crate) sens: f64,
    pub(crate) done: bool,
    pub(crate) rebrou: bool,
    pub(crate) iscomplete: bool,
    pub(crate) comptra: bool,
    pub(crate) clason_s1: bool,
    pub(crate) clason_s2: bool,
    pub(crate) check2d: bool,
    pub(crate) check: bool,
    pub(crate) twistflag1: bool,
    pub(crate) twistflag2: bool,
}

impl<'a> BRepBlendWalking<'a> {
    /// OCCT BRepBlend_Walking(Surf1, Surf2, Domain1, Domain2, HGuide)
    /// (BRepBlend_Walking.cxx L83-109).
    pub fn new(
        surf1: &'a Surface3,
        surf2: &'a Surface3,
        domain1: &'a BRepTopAdaptorTopolTool,
        domain2: &'a BRepTopAdaptorTopolTool,
        hguide: &'a Curve3,
    ) -> Self {
        BRepBlendWalking {
            previous_p: BlendPoint::new(),
            line: BRepBlendLine::new(),
            line_is_null: true,
            sol: vec![0.0; 4], // OCCT: math_Vector sol(1, 4)
            jalons: Vec::new(),
            surf1,
            surf2,
            domain1,
            domain2,
            recdomain1: domain1,
            recdomain2: domain2,
            hguide,
            to_correct_on_rst1: false,
            to_correct_on_rst2: false,
            corrected_param: 0.0,
            tolpoint3d: 0.0,
            tolgui: 0.0,
            pasmax: 0.0,
            fleche: 0.0,
            param: 0.0,
            sens: 1.0,
            done: false,
            rebrou: false,
            iscomplete: false,
            comptra: false,
            clason_s1: true,
            clason_s2: true,
            check2d: true,
            check: true,
            twistflag1: false,
            twistflag2: false,
        }
    }

    /// OCCT SetDomainsToRecadre(Domain1, Domain2) (BRepBlend_Walking.cxx
    /// L111-115).
    pub fn set_domains_to_recadre(
        &mut self,
        domain1: &'a BRepTopAdaptorTopolTool,
        domain2: &'a BRepTopAdaptorTopolTool,
    ) {
        self.recdomain1 = domain1;
        self.recdomain2 = domain2;
    }

    /// OCCT AddSingularPoint(P) (BRepBlend_Walking.cxx L117-141).
    pub fn add_singular_point(&mut self, p: BlendPoint) {
        if self.jalons.is_empty() {
            self.jalons.push(p);
        } else {
            let tp = p.parameter();
            let mut ti = self.jalons[0].parameter();
            let mut jj = 1usize;
            let mut ii = 1usize;
            while ii <= self.jalons.len() && tp > ti {
                jj = ii;
                ti = self.jalons[jj].parameter();
                ii += 1;
            }
            if tp > ti {
                // OCCT: jalons.InsertAfter(jj, P).
                self.jalons.insert(jj + 1 - 1, p);
            } else {
                // OCCT: jalons.InsertBefore(jj, P).
                self.jalons.insert(jj - 1, p);
            }
        }
    }

    /// OCCT Perform(Func, FuncInv, Pdep, Pmax, MaxStep, Tol3d, TolGuide,
    /// ParDep, Fleche, Appro) (BRepBlend_Walking.cxx L143-247).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        f: &mut dyn BlendFunction,
        f_inv: &mut dyn BlendFuncInv,
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
        let mut doextremities = true;
        if self.line_is_null {
            // OCCT: line = new BRepBlend_Line();
            self.line = BRepBlendLine::new();
            self.line_is_null = false;
        } else {
            self.line.clear();
            doextremities = false;
        }
        self.tolpoint3d = tol3d;
        self.tolgui = tolguide.abs();
        self.fleche = fleche.abs();
        self.rebrou = false;
        self.pasmax = max_step.abs();
        if pmax - pdep >= 0.0 {
            self.sens = 1.0;
        } else {
            self.sens = -1.0;
        }

        self.param = pdep;
        f.set_param(self.param);

        if appro {
            let mut tolerance = vec![0.0; 4];
            let mut infbound = vec![0.0; 4];
            let mut supbound = vec![0.0; 4];
            f.get_tolerance(&mut tolerance, self.tolpoint3d);
            f.get_bounds(&mut infbound, &mut supbound);
            // OCCT: math_FunctionSetRoot rsnld(Func, tolerance, 30);
            let fswd: &dyn FunctionSetWithDerivatives = &*f;
            let mut rsnld = FunctionSetRoot::new(fswd, &tolerance, 30);

            let fs: &mut dyn FunctionSetWithDerivatives = f;
            rsnld.perform(fs, par_dep, &infbound, &supbound, false);

            if !rsnld.is_done() {
                return;
            }
            self.sol.copy_from_slice(&rsnld.root());

            let situ1;
            let situ2;
            if self.clason_s1 {
                situ1 = self.domain1.classify(
                    glam::DVec2::new(self.sol[0], self.sol[1]),
                    tolerance[0].min(tolerance[1]),
                    false,
                );
            } else {
                situ1 = TopAbsState::In;
            }
            if self.clason_s2 {
                situ2 = self.domain2.classify(
                    glam::DVec2::new(self.sol[2], self.sol[3]),
                    tolerance[2].min(tolerance[3]),
                    false,
                );
            } else {
                situ2 = TopAbsState::In;
            }

            if situ1 != TopAbsState::In || situ2 != TopAbsState::In {
                return;
            }
        } else {
            self.sol.copy_from_slice(par_dep);
        }

        // OCCT: State = TestArret(Func, Blend_OK, false); (TestDeflection =
        // false, TestSolution = true, TestLengthStep = false)
        let state = self.test_arret(f, BlendStatus::Ok, false, true, false);
        if state != BlendStatus::Ok {
            return;
        }
        // Mettre a jour la ligne.
        // Correct first parameter if needed
        if self.to_correct_on_rst1 || self.to_correct_on_rst2 {
            self.previous_p.set_parameter(self.corrected_param);
        }
        self.line.append(self.previous_p.clone());

        if doextremities {
            // OCCT: BRepBlend_Extremity ptf1(P, U, V, Tol) — the 4-argument
            // constructor, Param = 0.
            let mut ptf1 =
                BRepBlendExtremity::new_on_surface(self.previous_p.point_on_s1(), self.sol[0], self.sol[1], 0.0, self.tolpoint3d);
            let mut ptf2 =
                BRepBlendExtremity::new_on_surface(self.previous_p.point_on_s2(), self.sol[2], self.sol[3], 0.0, self.tolpoint3d);
            if !self.previous_p.is_tangency_point() {
                ptf1.set_tangent(self.previous_p.tangent_on_s1());
                ptf2.set_tangent(self.previous_p.tangent_on_s2());
            }

            if self.sens > 0.0 {
                self.line.set_start_points(&ptf1, &ptf2);
            } else {
                self.line.set_end_points(&ptf1, &ptf2);
            }
        }

        self.internal_perform(f, f_inv, pmax);

        self.done = true;
    }

    /// OCCT PerformFirstSection(Func, Pdep, ParDep, Tol3d, TolGuide, Pos1,
    /// Pos2) (BRepBlend_Walking.cxx L249-295) — the first overload without
    /// recadre.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_first_section(
        &mut self,
        f: &mut dyn BlendFunction,
        pdep: f64,
        par_dep: &mut [f64],
        tol3d: f64,
        tolguide: f64,
        pos1: &mut TopAbsState,
        pos2: &mut TopAbsState,
    ) -> bool {
        self.iscomplete = false;
        self.comptra = false;
        self.line = BRepBlendLine::new();
        self.line_is_null = false;
        self.tolpoint3d = tol3d;
        self.tolgui = tolguide.abs();

        *pos1 = TopAbsState::Unknown;
        *pos2 = TopAbsState::Unknown;

        self.param = pdep;
        f.set_param(self.param);

        let mut tolerance = vec![0.0; 4];
        let mut infbound = vec![0.0; 4];
        let mut supbound = vec![0.0; 4];
        f.get_tolerance(&mut tolerance, self.tolpoint3d);
        f.get_bounds(&mut infbound, &mut supbound);
        // OCCT: math_FunctionSetRoot rsnld(Func, tolerance, 30);
        let fswd: &dyn FunctionSetWithDerivatives = &*f;
        let mut rsnld = FunctionSetRoot::new(fswd, &tolerance, 30);

        let fs: &mut dyn FunctionSetWithDerivatives = f;
        rsnld.perform(fs, par_dep, &infbound, &supbound, false);

        if !rsnld.is_done() {
            return false;
        }
        self.sol.copy_from_slice(&rsnld.root());
        par_dep.copy_from_slice(&self.sol);
        *pos1 = self.domain1.classify(
            glam::DVec2::new(self.sol[0], self.sol[1]),
            tolerance[0].min(tolerance[1]),
            false,
        );
        *pos2 = self.domain2.classify(
            glam::DVec2::new(self.sol[2], self.sol[3]),
            tolerance[2].min(tolerance[3]),
            false,
        );
        if *pos1 != TopAbsState::In || *pos2 != TopAbsState::In {
            return false;
        }

        // OCCT: TestArret(Func, Blend_OK, false);
        let _ = self.test_arret(f, BlendStatus::Ok, false, true, false);
        true
    }

    /// OCCT PerformFirstSection(Func, FuncInv, Pdep, Pmax, ParDep, Tol3d,
    /// TolGuide, RecOnS1, RecOnS2, Psol, ParSol) (BRepBlend_Walking.cxx
    /// L297-614) — the second overload with recadre.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_first_section_recate(
        &mut self,
        f: &mut dyn BlendFunction,
        f_inv: &mut dyn BlendFuncInv,
        pdep: f64,
        pmax: f64,
        par_dep: &[f64],
        tol3d: f64,
        tolguide: f64,
        rec_on_s1: bool,
        rec_on_s2: bool,
        psol: &mut f64,
        par_sol: &mut [f64],
    ) -> bool {
        self.iscomplete = false;
        self.comptra = false;
        self.line = BRepBlendLine::new();
        self.line_is_null = false;

        let mut w1;
        let mut w2;
        let extrapol;

        self.tolpoint3d = tol3d;
        self.tolgui = tolguide.abs();
        if pmax - pdep >= 0.0 {
            self.sens = 1.0;
        } else {
            self.sens = -1.0;
        }
        extrapol = (pmax - pdep).abs() / 50.0; // 2%

        let state;

        self.param = pdep;
        f.set_param(self.param);

        let mut tolerance = vec![0.0; 4];
        let mut infbound = vec![0.0; 4];
        let mut supbound = vec![0.0; 4];
        let mut solrst1 = vec![0.0; 4];
        let mut solrst2 = vec![0.0; 4];
        let mut ext1 = BRepBlendExtremity::new();
        let mut ext2 = BRepBlendExtremity::new();
        let mut index1 = 0i32;
        let mut index2 = 0i32;
        let mut isvtx1 = false;
        let mut isvtx2 = false;
        let mut vtx1 = HVertex::new();
        let mut vtx2 = HVertex::new();
        let mut corrected_u = 0.0;
        let mut corrected_v = 0.0;
        let mut corrected_pnt = DVec3::ZERO;

        f.get_tolerance(&mut tolerance, self.tolpoint3d);
        f.get_bounds(&mut infbound, &mut supbound);
        // OCCT: math_FunctionSetRoot rsnld(Func, tolerance, 30);
        let fswd: &dyn FunctionSetWithDerivatives = &*f;
        let mut rsnld = FunctionSetRoot::new(fswd, &tolerance, 30);

        let fs: &mut dyn FunctionSetWithDerivatives = f;
        rsnld.perform(fs, par_dep, &infbound, &supbound, false);

        if !rsnld.is_done() {
            return false;
        }
        self.sol.copy_from_slice(&rsnld.root());

        w1 = pmax;
        w2 = pmax;

        // OCCT passes the member sol by const reference into Recadre; rcad
        // splits the &mut self borrow with a copy.
        let sol_c = self.sol.clone();
        let recad1 = rec_on_s1
            && self.recadre(
                f_inv,
                true,
                &sol_c,
                &mut solrst1,
                &mut index1,
                &mut isvtx1,
                &mut vtx1,
                extrapol,
            );
        if recad1 {
            w1 = solrst1[1];
        }

        let recad2 = rec_on_s2
            && self.recadre(
                f_inv,
                false,
                &sol_c,
                &mut solrst2,
                &mut index2,
                &mut isvtx2,
                &mut vtx2,
                extrapol,
            );
        if recad2 {
            w2 = solrst2[1];
        }

        if !recad1 && !recad2 {
            return false;
        }

        if recad1 && recad2 {
            if (w1 - w2).abs() <= self.tolgui {
                // sol sur 1 et 2 a la fois
                state = BlendStatus::OnRst12;
                self.param = w1;
                par_sol[0] = solrst2[2];
                par_sol[1] = solrst2[3];
                par_sol[2] = solrst1[2];
                par_sol[3] = solrst1[3];
            } else if self.sens * (w2 - w1) < 0.0 {
                // on garde le plus grand
                // sol sur 1
                state = BlendStatus::OnRst1;
                self.param = w1;

                topol_tool_init(self.recdomain1);
                let mut nbarc = 1;
                while nbarc < index1 {
                    nbarc += 1;
                    topol_tool_next(self.recdomain1);
                }
                let p2d = hcurve2d_tool_value(
                    &topol_tool_value(self.recdomain1),
                    solrst1[0],
                );
                par_sol[0] = p2d.x;
                par_sol[1] = p2d.y;
                par_sol[2] = solrst1[2];
                par_sol[3] = solrst1[3];
            } else {
                // sol sur 2
                state = BlendStatus::OnRst2;
                self.param = w2;

                topol_tool_init(self.recdomain2);
                let mut nbarc = 1;
                while nbarc < index2 {
                    nbarc += 1;
                    topol_tool_next(self.recdomain2);
                }
                let p2d = hcurve2d_tool_value(
                    &topol_tool_value(self.recdomain2),
                    solrst2[0],
                );
                par_sol[0] = solrst2[2];
                par_sol[1] = solrst2[3];
                par_sol[2] = p2d.x;
                par_sol[3] = p2d.y;
            }
        } else if recad1 {
            // sol sur 1
            state = BlendStatus::OnRst1;
            self.param = w1;
            topol_tool_init(self.recdomain1);
            let mut nbarc = 1;
            while nbarc < index1 {
                nbarc += 1;
                topol_tool_next(self.recdomain1);
            }
            let p2d = hcurve2d_tool_value(
                &topol_tool_value(self.recdomain1),
                solrst1[0],
            );
            par_sol[0] = p2d.x;
            par_sol[1] = p2d.y;
            par_sol[2] = solrst1[2];
            par_sol[3] = solrst1[3];
            let the_pnt_on_rst =
                hsurface_tool_value(self.surf1, par_sol[0], par_sol[1]);
            let mut corrected_param = self.corrected_param;
            if self.correct_extremity_on_one_rst(
                1,
                par_sol[2],
                par_sol[3],
                self.param,
                the_pnt_on_rst,
                &mut corrected_u,
                &mut corrected_v,
                &mut corrected_pnt,
                &mut corrected_param,
            ) {
                self.to_correct_on_rst1 = true;
                self.corrected_param = corrected_param;
            }
        } else {
            // if (recad2) {
            // sol sur 2
            state = BlendStatus::OnRst2;
            self.param = w2;
            topol_tool_init(self.recdomain2);
            let mut nbarc = 1;
            while nbarc < index2 {
                nbarc += 1;
                topol_tool_next(self.recdomain2);
            }
            let p2d = hcurve2d_tool_value(
                &topol_tool_value(self.recdomain2),
                solrst2[0],
            );
            par_sol[0] = solrst2[2];
            par_sol[1] = solrst2[3];
            par_sol[2] = p2d.x;
            par_sol[3] = p2d.y;
            let the_pnt_on_rst =
                hsurface_tool_value(self.surf2, par_sol[2], par_sol[3]);
            let mut corrected_param = self.corrected_param;
            if self.correct_extremity_on_one_rst(
                2,
                par_sol[0],
                par_sol[1],
                self.param,
                the_pnt_on_rst,
                &mut corrected_u,
                &mut corrected_v,
                &mut corrected_pnt,
                &mut corrected_param,
            ) {
                self.to_correct_on_rst2 = true;
                self.corrected_param = corrected_param;
            }
        }

        *psol = self.param;
        self.sol.copy_from_slice(par_sol);
        f.set_param(self.param);
        // OCCT: State = TestArret(Func, State, false);
        let state = self.test_arret(f, state, false, true, false);
        match state {
            BlendStatus::OnRst1 => {
                self.make_extremity(&mut ext1, true, index1, solrst1[0], isvtx1, &vtx1);
                if self.to_correct_on_rst1 {
                    // OCCT: Ext2.SetValue(CorrectedPnt, CorrectedU,
                    // CorrectedV, tolpoint3d) — 4-argument overload, Param 0.
                    ext2.set_value(corrected_pnt, corrected_u, corrected_v, 0.0, self.tolpoint3d);
                } else {
                    ext2.set_value(
                        self.previous_p.point_on_s2(),
                        self.sol[2],
                        self.sol[3],
                        0.0,
                        self.tolpoint3d,
                    );
                }
            }

            BlendStatus::OnRst2 => {
                if self.to_correct_on_rst2 {
                    // OCCT: Ext1.SetValue(CorrectedPnt, CorrectedU,
                    // CorrectedV, tolpoint3d) — 4-argument overload, Param 0.
                    ext1.set_value(corrected_pnt, corrected_u, corrected_v, 0.0, self.tolpoint3d);
                } else {
                    ext1.set_value(
                        self.previous_p.point_on_s1(),
                        self.sol[0],
                        self.sol[1],
                        0.0,
                        self.tolpoint3d,
                    );
                }
                self.make_extremity(&mut ext2, false, index2, solrst2[0], isvtx2, &vtx2);
            }

            BlendStatus::OnRst12 => {
                self.make_extremity(&mut ext1, true, index1, solrst1[0], isvtx1, &vtx1);
                self.make_extremity(&mut ext2, false, index2, solrst2[0], isvtx2, &vtx2);
            }
            _ => {
                panic!("BRepBlend_Walking::PerformFirstSection : echec");
            }
        }
        if self.sens < 0.0 {
            self.line.set_end_points(&ext1, &ext2);
        } else {
            self.line.set_start_points(&ext1, &ext2);
        }
        true
    }

    /// OCCT Continu(Func, FuncInv, P) (BRepBlend_Walking.cxx L616-641).
    pub fn continu(&mut self, f: &mut dyn BlendFunction, f_inv: &mut dyn BlendFuncInv, p: f64) -> bool {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_Walking::Continu");
        }
        let first_bp = self.line.point(1).clone();
        let last_bp = self.line.point(self.line.nb_points()).clone();

        if p < first_bp.parameter() {
            self.sens = -1.0;
            self.previous_p = first_bp;
        } else if p > last_bp.parameter() {
            self.sens = 1.0;
            self.previous_p = last_bp;
        }

        self.param = self.previous_p.parameter();
        let (u1, v1) = self.previous_p.parameters_on_s1();
        self.sol[0] = u1;
        self.sol[1] = v1;
        let (u2, v2) = self.previous_p.parameters_on_s2();
        self.sol[2] = u2;
        self.sol[3] = v2;

        self.internal_perform(f, f_inv, p);
        true
    }

    /// OCCT Continu(Func, FuncInv, P, OnS1) (BRepBlend_Walking.cxx
    /// L643-741).
    pub fn continu_on_s(
        &mut self,
        f: &mut dyn BlendFunction,
        f_inv: &mut dyn BlendFuncInv,
        p: f64,
        on_s1: bool,
    ) -> bool {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_Walking::Continu");
        }
        let ext1;
        let ext2;
        if self.sens < 0.0 {
            ext1 = self.line.start_point_on_first().clone();
            ext2 = self.line.start_point_on_second().clone();
            if (on_s1 && ext1.nb_point_on_rst() == 0) || (!on_s1 && ext2.nb_point_on_rst() == 0) {
                return false;
            }
            self.previous_p = self.line.point(1).clone();
        } else {
            ext1 = self.line.end_point_on_first().clone();
            ext2 = self.line.end_point_on_second().clone();
            if (on_s1 && ext1.nb_point_on_rst() == 0) || (!on_s1 && ext2.nb_point_on_rst() == 0) {
                return false;
            }
            self.previous_p = self.line.point(self.line.nb_points()).clone();
        }

        let length = self.line.nb_points();
        self.param = self.previous_p.parameter();
        let (u1, v1) = self.previous_p.parameters_on_s1();
        self.sol[0] = u1;
        self.sol[1] = v1;
        let (u2, v2) = self.previous_p.parameters_on_s2();
        self.sol[2] = u2;
        self.sol[3] = v2;

        if on_s1 {
            self.clason_s1 = false;
        } else {
            self.clason_s2 = false;
        }

        self.internal_perform(f, f_inv, p);

        self.clason_s1 = true;
        self.clason_s2 = true;

        let newlength = self.line.nb_points();
        if self.sens < 0.0 {
            if (on_s1 && self.line.start_point_on_second().nb_point_on_rst() == 0)
                || (!on_s1 && self.line.start_point_on_first().nb_point_on_rst() == 0)
            {
                self.line.remove(1, newlength - length);
                self.line.set_start_points(&ext1, &ext2);
                return false;
            }
        } else {
            if (on_s1 && self.line.end_point_on_second().nb_point_on_rst() == 0)
                || (!on_s1 && self.line.end_point_on_first().nb_point_on_rst() == 0)
            {
                self.line.remove(length, newlength);
                self.line.set_end_points(&ext1, &ext2);
                return false;
            }
        }
        true
    }

    /// OCCT Complete(Func, FuncInv, Pmin) (BRepBlend_Walking.cxx L743-772).
    pub fn complete(&mut self, f: &mut dyn BlendFunction, f_inv: &mut dyn BlendFuncInv, pmin: f64) -> bool {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_Walking::Complete");
        }
        if self.iscomplete {
            return true;
        }

        if self.sens > 0.0 {
            self.previous_p = self.line.point(1).clone();
        } else {
            self.previous_p = self.line.point(self.line.nb_points()).clone();
        }

        self.sens = -self.sens;

        self.param = self.previous_p.parameter();
        let (u1, v1) = self.previous_p.parameters_on_s1();
        self.sol[0] = u1;
        self.sol[1] = v1;
        let (u2, v2) = self.previous_p.parameters_on_s2();
        self.sol[2] = u2;
        self.sol[3] = v2;

        self.internal_perform(f, f_inv, pmin);

        self.iscomplete = true;
        true
    }

    /// OCCT ClassificationOnS1(C) — inline header (cxx L774-778).
    pub fn classification_on_s1(&mut self, c: bool) {
        self.clason_s1 = c;
    }

    /// OCCT ClassificationOnS2(C) — inline header (cxx L780-784).
    pub fn classification_on_s2(&mut self, c: bool) {
        self.clason_s2 = c;
    }

    /// OCCT Check2d(C) — inline header (cxx L786-790).
    pub fn set_check2d(&mut self, c: bool) {
        self.check2d = c;
    }

    /// OCCT Check(C) — inline header (cxx L792-796).
    pub fn set_check(&mut self, c: bool) {
        self.check = c;
    }

    /// OCCT TwistOnS1() — inline header.
    pub fn twist_on_s1(&self) -> bool {
        self.twistflag1
    }

    /// OCCT TwistOnS2() — inline header.
    pub fn twist_on_s2(&self) -> bool {
        self.twistflag2
    }

    /// OCCT IsDone() — inline header.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Line() — inline header (throws StdFail_NotDone when not done).
    pub fn line(&self) -> &BRepBlendLine {
        if !self.done {
            panic!("StdFail_NotDone: BRepBlend_Walking::Line");
        }
        &self.line
    }
}
/// OCCT gp_Dir::Angle(Dir) — the angle between two directions in [0, PI].
pub(crate) fn gp_dir_angle(a: DVec3, b: DVec3) -> f64 {
    (a.dot(b) / (a.length() * b.length())).clamp(-1.0, 1.0).acos()
}

/// OCCT gce_MakePln(P1, P2, P3) — the plane through three points; only the
/// Axis().Direction() queries of the OCCT callers are represented.
pub(crate) struct GceMakePln {
    done: bool,
    axis_direction: DVec3,
}

impl GceMakePln {
    /// OCCT gce_MakePln::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT gp_Pln::Axis().Direction() — the plane normal.
    pub(crate) fn axis_direction(&self) -> DVec3 {
        self.axis_direction
    }
}

/// OCCT gp::Resolution() — the resolution below which a vector magnitude is
/// considered null.
pub(crate) const GP_RESOLUTION: f64 = 1.0e-12;

/// OCCT gce_MakePln(P1, P2, P3) (gce_MakePln.cxx — the three-point
/// constructor); the plane normal is the normalized P1P2 x P1P3.
pub(crate) fn gce_make_pln(p1: DVec3, p2: DVec3, p3: DVec3) -> GceMakePln {
    let n = (p2 - p1).cross(p3 - p1);
    if n.length() <= GP_RESOLUTION {
        return GceMakePln {
            done: false,
            axis_direction: DVec3::ZERO,
        };
    }
    GceMakePln {
        done: true,
        axis_direction: n.normalize(),
    }
}

/// OCCT ChFiDS_ElSpine::NbVertices() — PENDING: the vertex list lives on
/// the ChFiDS_ElSpine, not on the rcad guide `Curve3`; reports 0, which
/// selects the OCCT early-out of CorrectExtremityOnOneRst.
pub(crate) fn elspine_nb_vertices(_hguide: &Curve3) -> i32 {
    0
}

/// OCCT ChFiDS_ElSpine::VertexWithTangent(ind) — PENDING (see
/// elspine_nb_vertices); only reached for a non-empty vertex list.
pub(crate) fn elspine_vertex_with_tangent(_hguide: &Curve3, ind: i32) -> (DVec3, DVec3) {
    let _ = ind;
    unreachable!("ChFiDS_ElSpine::VertexWithTangent — pending vertex list on the rcad guide curve")
}

/// OCCT ChFiDS_ElSpine::GetSavedFirstParameter() — the OCCT constructor
/// default is Precision::Infinite() (ChFiDS_ElSpine.cxx L37-45); the saved
/// values are set only by the pending SetOrigin machinery.
pub(crate) fn elspine_get_saved_first_parameter(_hguide: &Curve3) -> f64 {
    rcad_kernel::core::precision::INFINITE_VALUE
}

/// See elspine_get_saved_first_parameter.
pub(crate) fn elspine_get_saved_last_parameter(_hguide: &Curve3) -> f64 {
    rcad_kernel::core::precision::INFINITE_VALUE
}

// =========================================================================
// OCCT BRepBlend_HCurve2dTool (BRepBlend_HCurve2dTool.lxx L8-120) — the
// static dispatchers over the restriction curve, mapped onto the rcad
// Curve2d evaluation trait.
// =========================================================================

/// OCCT BRepBlend_HCurve2dTool::FirstParameter(C).
pub(crate) fn hcurve2d_tool_first_parameter(c: &Curve2d) -> f64 {
    c.default_domain()[0]
}

/// OCCT BRepBlend_HCurve2dTool::LastParameter(C).
pub(crate) fn hcurve2d_tool_last_parameter(c: &Curve2d) -> f64 {
    c.default_domain()[1]
}

/// OCCT BRepBlend_HCurve2dTool::Value(C, U).
pub(crate) fn hcurve2d_tool_value(c: &Curve2d, u: f64) -> DVec2 {
    c.point_at(u)
}

/// OCCT BRepBlend_HCurve2dTool::D0(C, U, P).
pub(crate) fn hcurve2d_tool_d0(c: &Curve2d, u: f64) -> DVec2 {
    c.point_at(u)
}

/// OCCT BRepBlend_HCurve2dTool::D1(C, U, P, V).
pub(crate) fn hcurve2d_tool_d1(c: &Curve2d, u: f64) -> (DVec2, DVec2) {
    (c.point_at(u), c.derivative_at(u))
}

// =========================================================================
// OCCT BRepBlend_HCurveTool (BRepBlend_HCurveTool.lxx) — the static
// dispatchers over the guide curve.
// =========================================================================

/// OCCT BRepBlend_HCurveTool::Period(C) — 2*PI for the periodic conics
/// (the period carrier of the rcad Curve3 is the pending boundary).
pub(crate) fn hcurve_tool_period(c: &Curve3) -> f64 {
    match c {
        Curve3::Circle(_) | Curve3::Ellipse(_) => std::f64::consts::PI * 2.0,
        _ => 0.0,
    }
}

// =========================================================================
// OCCT Adaptor3d_HSurfaceTool (hxx L20-230) — the static dispatchers over
// the support surface, mapped onto the rcad Surface3 evaluation surface.
// =========================================================================

/// OCCT Adaptor3d_HSurfaceTool::Value(S, U, V).
pub(crate) fn hsurface_tool_value(s: &Surface3, u: f64, v: f64) -> DVec3 {
    s.point_at(u, v)
}

/// OCCT Adaptor3d_HSurfaceTool::D1(S, U, V, P, D1U, D1V).
pub(crate) fn hsurface_tool_d1(s: &Surface3, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
    s.derivatives(u, v)
}

/// OCCT Adaptor3d_HSurfaceTool::D2(S, U, V, P, D1U, D1V, D2U, D2V, D2UV).
#[allow(clippy::type_complexity)]
pub(crate) fn hsurface_tool_d2(s: &Surface3, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
    s.derivatives2(u, v)
}

/// OCCT Adaptor3d_HSurfaceTool::FirstUParameter / LastUParameter /
/// FirstVParameter / LastVParameter — the natural UV box.
pub(crate) fn hsurface_tool_parameters(s: &Surface3) -> (f64, f64, f64, f64) {
    let d = s.default_domain();
    (d[0], d[1], d[2], d[3])
}

/// OCCT Adaptor3d_HSurfaceTool::UResolution(S, R).
pub(crate) fn hsurface_tool_u_resolution(s: &Surface3, r: f64) -> f64 {
    s.u_resolution(r)
}

/// OCCT Adaptor3d_HSurfaceTool::VResolution(S, R).
pub(crate) fn hsurface_tool_v_resolution(s: &Surface3, r: f64) -> f64 {
    s.v_resolution(r)
}

/// OCCT Adaptor3d_HSurfaceTool::UPeriod(S) — 2*PI on the revolved quadrics
/// (the analytic UPeriod; periodic BSpline V-periods are the pending
/// boundary of the rcad surface layer).
pub(crate) fn hsurface_tool_u_period(s: &Surface3) -> f64 {
    match s {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
            std::f64::consts::PI * 2.0
        }
        _ => 0.0,
    }
}

/// OCCT Adaptor3d_HSurfaceTool::VPeriod(S) — only the torus is V-periodic
/// among the analytic surfaces (the sphere V-period does not exist in
/// OCCT either).
pub(crate) fn hsurface_tool_v_period(s: &Surface3) -> f64 {
    match s {
        Surface3::Torus(_) => std::f64::consts::PI * 2.0,
        _ => 0.0,
    }
}

// =========================================================================
// OCCT BRepBlend_BlendTool (BRepBlend_BlendTool.cxx L26-142 + .lxx) — the
// static helpers over the restriction arcs.
// =========================================================================

/// OCCT BRepBlend_BlendTool::Project(P, S, C, Paramproj, Dist)
/// (BRepBlend_BlendTool.cxx L36-77) — orthogonal projection of a point on a
/// curve; returns (ok, Paramproj, Dist).
pub(crate) fn blend_tool_project(p: &DVec2, _surf: &Surface3, c: &Curve2d) -> (bool, f64, f64) {
    let mut paramproj = hcurve2d_tool_first_parameter(c);
    let mut p2d = hcurve2d_tool_d0(c, paramproj);
    let mut dist = p2d.distance(*p);

    let t = hcurve2d_tool_last_parameter(c);
    p2d = hcurve2d_tool_d0(c, t);
    if p2d.distance(*p) < dist {
        paramproj = t;
        dist = p2d.distance(*p);
    }

    let eps_x = 1.0e-8;
    let nbu = 20;
    let tol = 1.0e-5;
    // OCCT: Extrema_EPCOfExtPC2d extrema(P, *C, Nbu, epsX, Tol);
    let extrema = EPCOfExtPC2d::with_full_domain(*p, c, nbu, eps_x, tol);
    if !extrema.is_done() {
        return (true, paramproj, dist);
    }

    let nbext = extrema.nb_ext();
    let mut a_dist2 = dist * dist;
    for i in 1..=nbext {
        if extrema.square_distance(i) < a_dist2 {
            a_dist2 = extrema.square_distance(i);
            paramproj = extrema.point(i).parameter();
        }
    }
    dist = a_dist2.sqrt();

    (true, paramproj, dist)
}

/// OCCT BRepBlend_BlendTool::Inters(P1, P2, S, C, Param, Dist)
/// (BRepBlend_BlendTool.cxx L79-116) — intersection of a segment with a
/// curve.  PENDING: Geom2dInt_GInter is not translated yet; the placeholder
/// reports no intersection and the OCCT control flow falls back to
/// Project.
pub(crate) fn blend_tool_inters(
    _p1: &DVec2,
    _p2: &DVec2,
    _surf: &Surface3,
    _c: &Curve2d,
) -> (bool, f64, f64) {
    (false, 0.0, 0.0)
}

/// OCCT BRepBlend_BlendTool::CurveOnSurf(C, S) (.lxx) — identity on the
/// pcurve.
pub(crate) fn blend_tool_curve_on_surf(c: &Curve2d, _surf: &Surface3) -> Curve2d {
    c.clone()
}

/// OCCT BRepBlend_BlendTool::Bounds(A, Ufirst, Ulast) (.cxx L133-141).
pub(crate) fn blend_tool_bounds(a: &Curve2d, ufirst: &mut f64, ulast: &mut f64) {
    *ufirst = hcurve2d_tool_first_parameter(a);
    *ulast = hcurve2d_tool_last_parameter(a);
}

/// OCCT BRepBlend_BlendTool::Tolerance(V, A) (.lxx) — V->Resolution(A).
pub(crate) fn blend_tool_tolerance(v: &HVertex, a: &Curve2d) -> f64 {
    v.resolution(a)
}

/// OCCT BRepBlend_BlendTool::Parameter(V, C) (.lxx) — V->Parameter(C).
pub(crate) fn blend_tool_parameter(v: &HVertex, c: &Curve2d) -> f64 {
    v.parameter(c)
}

// =========================================================================
// OCCT Adaptor3d_TopolTool iteration — PENDING on the placeholder
// BRepTopAdaptorTopolTool (chfi3d_builder_2): the restriction-arc list and
// the vertex iterator are not stored yet.  The helpers below reproduce the
// empty-domain behaviour (More() == false), which drives the OCCT failure
// paths of Recadre / ArcToRecadre and the empty vertex loops.
// =========================================================================

/// OCCT Adaptor3d_TopolTool::Init().
pub(crate) fn topol_tool_init(_t: &BRepTopAdaptorTopolTool) {}

/// OCCT Adaptor3d_TopolTool::More() — the placeholder domain carries no
/// restriction arcs.
pub(crate) fn topol_tool_more(_t: &BRepTopAdaptorTopolTool) -> bool {
    false
}

/// OCCT Adaptor3d_TopolTool::Next().
pub(crate) fn topol_tool_next(_t: &BRepTopAdaptorTopolTool) {}

/// OCCT Adaptor3d_TopolTool::Value() — only reachable when More() is true,
/// which the placeholder domain never reports.
pub(crate) fn topol_tool_value(_t: &BRepTopAdaptorTopolTool) -> Curve2d {
    unreachable!("Adaptor3d_TopolTool::Value: no current restriction (pending BRepTopAdaptor_TopolTool arc iteration)")
}

/// OCCT Adaptor3d_TopolTool::Initialize(C).
pub(crate) fn topol_tool_initialize<'t>(_t: &'t BRepTopAdaptorTopolTool, _c: &Curve2d) {}

/// OCCT Adaptor3d_TopolTool::InitVertexIterator().
pub(crate) fn topol_tool_init_vertex_iterator(_t: &BRepTopAdaptorTopolTool) {}

/// OCCT Adaptor3d_TopolTool::MoreVertex() — the placeholder domain carries
/// no vertices.
pub(crate) fn topol_tool_more_vertex(_t: &BRepTopAdaptorTopolTool) -> bool {
    false
}

/// OCCT Adaptor3d_TopolTool::Vertex() — only reachable when MoreVertex()
/// is true.
pub(crate) fn topol_tool_vertex(_t: &BRepTopAdaptorTopolTool) -> &HVertex {
    unreachable!("Adaptor3d_TopolTool::Vertex: no current vertex (pending vertex iterator)")
}

/// OCCT Adaptor3d_TopolTool::NextVertex().
pub(crate) fn topol_tool_next_vertex(_t: &BRepTopAdaptorTopolTool) {}

/// OCCT Adaptor3d_TopolTool::Identical(V1, V2).
pub(crate) fn topol_tool_identical(_t: &BRepTopAdaptorTopolTool, _v1: &HVertex, _v2: &HVertex) -> bool {
    false
}

// Unused-import guards for the pending-boundary pieces above.
#[allow(dead_code)]
fn _assert_pending_types(_s: NormalStatus, _d: DerivativeStatus, _t: Transition) {}
/// OCCT RecadreIfPeriodic(NewU, NewV, OldU, OldV, UPeriod, VPeriod)
/// (BRepBlend_Walking.cxx L1861-1879).
pub(crate) fn recadre_if_periodic(
    new_u: &mut f64,
    new_v: &mut f64,
    old_u: f64,
    old_v: f64,
    u_period: f64,
    v_period: f64,
) {
    if u_period > 0.0 {
        let sign = if *new_u < old_u { 1.0 } else { -1.0 };
        while (*new_u - old_u).abs() > u_period / 2.0 {
            *new_u += sign * u_period;
        }
    }
    if v_period > 0.0 {
        let sign = if *new_v < old_v { 1.0 } else { -1.0 };
        while (*new_v - old_v).abs() > v_period / 2.0 {
            *new_v += sign * v_period;
        }
    }
}

/// OCCT evalpinit(parinit, previousP, parprec, param, infbound, supbound,
/// classonS1, classonS2) (BRepBlend_Walking.cxx L1881-1930).
pub(crate) fn evalpinit(
    parinit: &mut [f64],
    previous_p: &BlendPoint,
    parprec: f64,
    param: f64,
    infbound: &[f64],
    supbound: &[f64],
    classon_s1: bool,
    classon_s2: bool,
) {
    if previous_p.is_tangency_point() {
        let (u1, v1) = previous_p.parameters_on_s1();
        let (u2, v2) = previous_p.parameters_on_s2();
        parinit[0] = u1;
        parinit[1] = v1;
        parinit[2] = u2;
        parinit[3] = v2;
    } else {
        let mut inside = true;
        let (mut u1, mut v1) = previous_p.parameters_on_s1();
        let (mut u2, mut v2) = previous_p.parameters_on_s2();
        let du1 = previous_p.tangent_2d_on_s1().x;
        let dv1 = previous_p.tangent_2d_on_s1().y;
        let du2 = previous_p.tangent_2d_on_s2().x;
        let dv2 = previous_p.tangent_2d_on_s2().y;
        let step = param - parprec;
        u1 += step * du1;
        v1 += step * dv1;
        if classon_s1 {
            if (u1 < infbound[0]) || (u1 > supbound[0]) {
                inside = false;
            }
            if (v1 < infbound[1]) || (v1 > supbound[1]) {
                inside = false;
            }
        }
        u2 += step * du2;
        v2 += step * dv2;
        if classon_s2 {
            if (u2 < infbound[2]) || (u2 > supbound[2]) {
                inside = false;
            }
            if (v2 < infbound[3]) || (v2 > supbound[3]) {
                inside = false;
            }
        }

        if inside {
            parinit[0] = u1;
            parinit[1] = v1;
            parinit[2] = u2;
            parinit[3] = v2;
        } else {
            // on ne joue pas au plus malin
            let (u1, v1) = previous_p.parameters_on_s1();
            let (u2, v2) = previous_p.parameters_on_s2();
            parinit[0] = u1;
            parinit[1] = v1;
            parinit[2] = u2;
            parinit[3] = v2;
        }
    }
}
