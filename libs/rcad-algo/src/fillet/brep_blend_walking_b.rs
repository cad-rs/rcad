//! OCCT BRepBlend_Walking solver internals (TKFillet/BRepBlend) —
//! BRepBlend_Walking.cxx L798-2769: TestArret, CheckDeflection,
//! ArcToRecadre, Recadre, Transition, MakeExtremity, MakeSingularExtremity,
//! InternalPerform, CorrectExtremityOnOneRst plus the file statics
//! (CosRef3D, RecadreIfPeriodic, evalpinit).
//!
//! Split from `brep_blend_walking` per the 2000-line guideline.
//!
//! Pending dependencies (plan 0.6 — OCCT-named placeholders + failure path):
//! - `BRepBlend_BlendTool::Inters` needs Geom2dInt_GInter (pending);
//!   returns false, and the OCCT control flow falls back to Project.
//! - `Adaptor3d_TopolTool` restriction-arc / vertex iteration is not
//!   carried by the placeholder `BRepTopAdaptorTopolTool` (chfi3d_builder_2);
//!   the iteration helpers report an empty domain, which drives the OCCT
//!   empty-domain paths (ArcToRecadre -> 0, Recadre -> false, empty
//!   vertex loops).
//! - The ChFiDS_ElSpine vertex list and saved parameters are not carried by
//!   the rcad guide `Curve3`; NbVertices reports 0 (the OCCT early-out of
//!   CorrectExtremityOnOneRst) and the saved parameters report the OCCT
//!   constructor default (Precision::Infinite).

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, Curve2dEval as _, CurveEval as _, SurfaceEval as _};
use rcad_kernel::math::cs_lib::{normal_from_derivatives, normal_max_order, DerivativeStatus};
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};
use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};

use crate::geomalgo::int_patch::transitions::{make_transition, Transition};
use rcad_kernel::base::extrema::ExtPS;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;

use super::brep_blend::BlendStatus;
use super::brep_blend_extremity::BRepBlendExtremity;
use super::brep_blend_func_inv::BlendFuncInv;
use super::brep_blend_function::BlendFunction;
use super::brep_blend_line::IntSurfTypeTrans;
use super::brep_blend_point::BlendPoint;
use super::brep_blend_walking::BRepBlendWalking;
use super::chfi3d_builder_2::TopAbsState;
use crate::topalgo::adaptor3d::hvertex::HVertex;

// Support-tool layer shared with the walking front end (declared in
// brep_blend_walking.rs).
use super::brep_blend_walking::{
    blend_tool_bounds, blend_tool_curve_on_surf, blend_tool_inters, blend_tool_parameter,
    blend_tool_project, blend_tool_tolerance, elspine_get_saved_first_parameter,
    elspine_get_saved_last_parameter, elspine_nb_vertices, elspine_vertex_with_tangent,
    gce_make_pln, gp_dir_angle, hcurve2d_tool_d1, hcurve2d_tool_first_parameter,
    hcurve2d_tool_last_parameter, hcurve2d_tool_value, hcurve_tool_period, hsurface_tool_d1,
    hsurface_tool_d2, hsurface_tool_parameters, hsurface_tool_u_period,
    hsurface_tool_u_resolution, hsurface_tool_v_period, hsurface_tool_v_resolution,
    evalpinit, hsurface_tool_value, recadre_if_periodic, topol_tool_identical,
    topol_tool_initialize, topol_tool_init, topol_tool_init_vertex_iterator, topol_tool_more,
    topol_tool_more_vertex, topol_tool_next, topol_tool_next_vertex, topol_tool_value,
    topol_tool_vertex, GP_RESOLUTION,
};

/// OCCT BRepBlend_Walking.cxx L1857: `static const double CosRef3D = 0.88;`
/// (used by InternalPerform for the guide deflection test).
const COS_REF_3D: f64 = 0.88;

/// OCCT `RealLast()` — the sentinel for "no candidate found".
const REAL_LAST: f64 = f64::MAX;

impl BRepBlendWalking<'_> {
    /// OCCT TestArret(Function, State, TestDefl, TestSolu, TestLengthStep)
    /// (BRepBlend_Walking.cxx L798-1000).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn test_arret(
        &mut self,
        f: &mut dyn BlendFunction,
        state: BlendStatus,
        test_defl: bool,
        test_solu: bool,
        test_length_step: bool,
    ) -> BlendStatus {
        let mut v1 = DVec3::ZERO;
        let mut v2 = DVec3::ZERO;
        let mut loctwist1 = false;
        let mut loctwist2 = false;
        let mut tolsolu = self.tolpoint3d;

        if !test_solu {
            tolsolu *= 1000.0; // Ca doit toujours etre bon
        }
        if f.is_solution(&self.sol, tolsolu) {
            let curpointistangent = f.is_tangency_point();
            let pt1 = f.point_on_s1();
            let pt2 = f.point_on_s2();
            let mut curpoint = BlendPoint::new();
            if curpointistangent {
                curpoint.set_value_on_2_surfaces(
                    pt1,
                    pt2,
                    self.param,
                    self.sol[0],
                    self.sol[1],
                    self.sol[2],
                    self.sol[3],
                );
            } else {
                v1 = f.tangent_on_s1();
                v2 = f.tangent_on_s2();
                let v12d = f.tangent_2d_on_s1();
                let v22d = f.tangent_2d_on_s2();
                curpoint.set_value_on_2_surfaces_with_tangents(
                    pt1,
                    pt2,
                    self.param,
                    self.sol[0],
                    self.sol[1],
                    self.sol[2],
                    self.sol[3],
                    v1,
                    v2,
                    v12d,
                    v22d,
                );
                if f.twist_on_s1() {
                    loctwist1 = true;
                }
                if f.twist_on_s2() {
                    loctwist2 = true;
                }
            }

            let mut state1;
            let mut state2;
            if test_defl && self.check {
                // Verification du critere de fleche sur chaque surface
                // et sur la ligne guide

                state1 = self.check_deflection(true, &curpoint);
                state2 = self.check_deflection(false, &curpoint);
            } else {
                state1 = BlendStatus::Ok;
                state2 = BlendStatus::Ok;
                if test_length_step {
                    // On verifie juste que le pas n'est pas trop grand
                    // (Cas des prolongements foireux)
                    let mut inf = vec![0.0; 4];
                    let mut sup = vec![0.0; 4];
                    f.get_bounds(&mut inf, &mut sup);
                    // OCCT: sup -= inf; sup *= 0.05; (Pas max : 5% du domaine)
                    for i in 0..4 {
                        sup[i] = (sup[i] - inf[i]) * 0.05;
                    }

                    let (curparamu1, curparamv1) = curpoint.parameters_on_s1();
                    let (prevparamu1, prevparamv1) = self.previous_p.parameters_on_s1();
                    if (curparamu1 - prevparamu1).abs() > sup[0] {
                        state1 = BlendStatus::StepTooLarge;
                    }
                    if (curparamv1 - prevparamv1).abs() > sup[1] {
                        state1 = BlendStatus::StepTooLarge;
                    }
                    let (curparamu2, curparamv2) = curpoint.parameters_on_s2();
                    let (prevparamu2, prevparamv2) = self.previous_p.parameters_on_s2();
                    if (curparamu2 - prevparamu2).abs() > sup[2] {
                        state2 = BlendStatus::StepTooLarge;
                    }
                    if (curparamv2 - prevparamv2).abs() > sup[3] {
                        state2 = BlendStatus::StepTooLarge;
                    }
                }
            }

            if state1 == BlendStatus::Backward {
                state1 = BlendStatus::StepTooLarge;
                self.rebrou = true;
            }

            if state2 == BlendStatus::Backward {
                state2 = BlendStatus::StepTooLarge;
                self.rebrou = true;
            }

            if state1 == BlendStatus::StepTooLarge || state2 == BlendStatus::StepTooLarge {
                return BlendStatus::StepTooLarge;
            }

            // Ici seulement on peut statuer sur le twist
            // Car les rejet ont ete effectue (BUC60322)
            if loctwist1 {
                self.twistflag1 = true;
            }
            if loctwist2 {
                self.twistflag2 = true;
            }

            if !self.comptra && !curpointistangent {
                let mut tgp1 = DVec3::ZERO;
                let mut tgp2 = DVec3::ZERO;
                let mut nor1 = DVec3::ZERO;
                let mut nor2 = DVec3::ZERO;
                f.tangent(
                    self.sol[0],
                    self.sol[1],
                    self.sol[2],
                    self.sol[3],
                    &mut tgp1,
                    &mut tgp2,
                    &mut nor1,
                    &mut nor2,
                );
                let nor1 = nor1.normalize();
                let nor2 = nor2.normalize();
                let mut testra = tgp1.dot(nor1.cross(v1));
                if testra.abs() > CONFUSION {
                    let mut tras1 = IntSurfTypeTrans::In;
                    if (testra > 0.0 && !loctwist1) || (testra < 0.0 && loctwist1) {
                        tras1 = IntSurfTypeTrans::Out;
                    }

                    testra = tgp2.dot(nor2.cross(v2));
                    if testra.abs() > CONFUSION {
                        let mut tras2 = IntSurfTypeTrans::Out;
                        if (testra > 0.0 && !loctwist2) || (testra < 0.0 && loctwist2) {
                            tras2 = IntSurfTypeTrans::In;
                        }
                        self.comptra = true;
                        self.line.set_transitions(tras1, tras2);
                    }
                }
            }

            if state1 == BlendStatus::Ok || state2 == BlendStatus::Ok {
                self.previous_p = curpoint;
                return state;
            }

            if state1 == BlendStatus::StepTooSmall && state2 == BlendStatus::StepTooSmall {
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
        } else {
            return BlendStatus::StepTooLarge;
        }
    }

    /// OCCT CheckDeflection(OnFirst, CurPoint) (BRepBlend_Walking.cxx
    /// L1002-1116).
    pub(crate) fn check_deflection(&mut self, on_first: bool, cur_point: &BlendPoint) -> BlendStatus {
        // regle par tests dans U4 correspond a 11.478 d
        let cos_ref_3d: f64 = 0.98;

        let cos_ref_2d: f64 = 0.88; // correspond a 25 d

        let curpointistangent = cur_point.is_tangency_point();
        let prevpointistangent = self.previous_p.is_tangency_point();

        let (psurf, tgsurf, prev_p, prev_tg, tolu, tolv);
        if on_first {
            psurf = cur_point.point_on_s1();
            tgsurf = if !curpointistangent {
                cur_point.tangent_on_s1()
            } else {
                DVec3::ZERO
            };
            prev_p = self.previous_p.point_on_s1();
            prev_tg = if !prevpointistangent {
                self.previous_p.tangent_on_s1()
            } else {
                DVec3::ZERO
            };
            tolu = hsurface_tool_u_resolution(self.surf1, self.tolpoint3d);
            tolv = hsurface_tool_v_resolution(self.surf1, self.tolpoint3d);
        } else {
            psurf = cur_point.point_on_s2();
            tgsurf = if !curpointistangent {
                cur_point.tangent_on_s2()
            } else {
                DVec3::ZERO
            };
            prev_p = self.previous_p.point_on_s2();
            prev_tg = if !prevpointistangent {
                self.previous_p.tangent_on_s2()
            } else {
                DVec3::ZERO
            };
            tolu = hsurface_tool_u_resolution(self.surf2, self.tolpoint3d);
            tolv = hsurface_tool_v_resolution(self.surf2, self.tolpoint3d);
        }

        let corde = psurf - prev_p;
        let norme = corde.length_squared();
        //  if(!curpointistangent) curNorme = Tgsurf.SquareMagnitude();
        let prev_norme = if !prevpointistangent {
            prev_tg.length_squared()
        } else {
            0.0
        };

        let toler3d = 0.01 * self.tolpoint3d;
        if norme <= toler3d * toler3d {
            // il faudra peut etre  forcer meme point
            return BlendStatus::SamePoints;
        }
        if !prevpointistangent {
            if prev_norme <= toler3d * toler3d {
                return BlendStatus::SamePoints;
            }
            let cosi = self.sens * corde.dot(prev_tg);
            if cosi < 0.0 {
                // angle 3d>pi/2. --> retour arriere
                return BlendStatus::Backward;
            }

            let cosi2 = cosi * cosi / prev_norme / norme;
            if cosi2 < cos_ref_3d {
                return BlendStatus::StepTooLarge;
            }
        }

        if !curpointistangent {
            // Voir s il faut faire le controle sur le signe de prevtg*Tgsurf
            let cosi = self.sens * corde.dot(tgsurf);
            let cosi2 = cosi * cosi / tgsurf.length_squared() / norme;
            if cosi2 < cos_ref_3d || cosi < 0.0 {
                return BlendStatus::StepTooLarge;
            }
        }

        if self.check2d {
            let (curparamu, curparamv, tgonsurf, prevparamu, prevparamv, previousd2d);
            if on_first {
                let (cu, cv) = cur_point.parameters_on_s1();
                curparamu = cu;
                curparamv = cv;
                tgonsurf = if !curpointistangent {
                    cur_point.tangent_2d_on_s1()
                } else {
                    DVec2::ZERO
                };
                let (pu, pv) = self.previous_p.parameters_on_s1();
                prevparamu = pu;
                prevparamv = pv;
                previousd2d = if !prevpointistangent {
                    self.previous_p.tangent_2d_on_s1()
                } else {
                    DVec2::ZERO
                };
            } else {
                let (cu, cv) = cur_point.parameters_on_s2();
                curparamu = cu;
                curparamv = cv;
                tgonsurf = if !curpointistangent {
                    cur_point.tangent_2d_on_s2()
                } else {
                    DVec2::ZERO
                };
                let (pu, pv) = self.previous_p.parameters_on_s2();
                prevparamu = pu;
                prevparamv = pv;
                previousd2d = if !prevpointistangent {
                    self.previous_p.tangent_2d_on_s2()
                } else {
                    DVec2::ZERO
                };
            }

            let du = curparamu - prevparamu;
            let dv = curparamv - prevparamv;
            let duv = du * du + dv * dv;
            if du.abs() < tolu && dv.abs() < tolv {
                // il faudra peut etre  forcer meme point
                return BlendStatus::SamePoints; // point confondu 2d
            }
            if !prevpointistangent {
                if previousd2d.x.abs() < tolu && previousd2d.y.abs() < tolv {
                    // il faudra peut etre  forcer meme point
                    return BlendStatus::SamePoints; // point confondu 2d
                }
                let cosi = self.sens * (du * previousd2d.x + dv * previousd2d.y);
                if cosi < 0.0 {
                    return BlendStatus::Backward;
                }
            }
            if !curpointistangent {
                // Voir s il faut faire le controle sur le signe de Cosi
                let cosi = self.sens * (du * tgonsurf.x + dv * tgonsurf.y) / tgonsurf.length();
                let cosi2 = cosi * cosi / duv;
                if cosi2 < cos_ref_2d || cosi < 0.0 {
                    return BlendStatus::StepTooLarge;
                }
            }
        }
        if !curpointistangent && !prevpointistangent {
            // Estimation de la fleche courante
            let fleche_courante =
                (prev_tg.normalize() - tgsurf.normalize()).length_squared() * norme / 64.0;

            if fleche_courante <= 0.25 * self.fleche * self.fleche {
                return BlendStatus::StepTooSmall;
            }
            if fleche_courante > self.fleche * self.fleche {
                // pas trop grand : commentaire interessant
                return BlendStatus::StepTooLarge;
            }
        }
        BlendStatus::Ok
    }

    /// OCCT ArcToRecadre(OnFirst, theSol, PrevIndex, lastpt2d, pt2d, ponarc)
    /// (BRepBlend_Walking.cxx L1118-1206).
    pub(crate) fn arc_to_recadre(
        &mut self,
        on_first: bool,
        the_sol: &[f64],
        prev_index: i32,
        lastpt2d: &mut DVec2,
        pt2d: &mut DVec2,
        ponarc: &mut f64,
    ) -> i32 {
        let mut index_sol = 0i32;
        let mut nbarc = 0i32;
        let mut ok;
        let byinter = self.line.nb_points() != 0;
        let mut okinter = false;
        let mut distmin = REAL_LAST;
        let mut uprev = 0.0;
        let mut vprev = 0.0;
        let mut prm = 0.0;
        let mut dist = 0.0;

        if on_first {
            if byinter {
                let (u, v) = self.previous_p.parameters_on_s1();
                uprev = u;
                vprev = v;
            }
            *pt2d = DVec2::new(the_sol[0], the_sol[1]);
        } else {
            if byinter {
                let (u, v) = self.previous_p.parameters_on_s2();
                uprev = u;
                vprev = v;
            }
            *pt2d = DVec2::new(the_sol[2], the_sol[3]);
        }
        *lastpt2d = DVec2::new(uprev, vprev);
        let iter = if on_first { self.recdomain1 } else { self.recdomain2 };
        topol_tool_init(iter);
        while topol_tool_more(iter) {
            nbarc += 1;
            ok = false;
            if on_first {
                if byinter {
                    let r = blend_tool_inters(pt2d, lastpt2d, self.surf1, &topol_tool_value(iter));
                    ok = r.0;
                    okinter = r.0;
                    prm = r.1;
                    dist = r.2;
                }
                if !ok {
                    let r = blend_tool_project(pt2d, self.surf1, &topol_tool_value(iter));
                    ok = r.0;
                    prm = r.1;
                    dist = r.2;
                }
            } else {
                if byinter {
                    let r = blend_tool_inters(pt2d, lastpt2d, self.surf2, &topol_tool_value(iter));
                    ok = r.0;
                    okinter = r.0;
                    prm = r.1;
                    dist = r.2;
                }
                if !ok {
                    let r = blend_tool_project(pt2d, self.surf2, &topol_tool_value(iter));
                    ok = r.0;
                    prm = r.1;
                    dist = r.2;
                }
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
            topol_tool_next(iter);
        }
        index_sol
    }

    /// OCCT Recadre(FuncInv, OnFirst, theSol, solrst, Indexsol, IsVtx, Vtx,
    /// Extrap) (BRepBlend_Walking.cxx L1208-1630).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn recadre(
        &mut self,
        f_inv: &mut dyn BlendFuncInv,
        on_first: bool,
        the_sol: &[f64],
        solrst: &mut [f64],
        indexsol: &mut i32,
        is_vtx: &mut bool,
        vtx: &mut HVertex,
        extrap: f64,
    ) -> bool {
        let mut jalons_trouve = false;
        let mut recadre = true;
        let byinter = self.line.nb_points() != 0;
        let mut le_jalon = 0usize;

        let mut dist;
        let mut prm = 0.0;
        let mut pmin = 0.0;
        let mut vtol;

        let mut toler = vec![0.0; 4];
        let mut infb = vec![0.0; 4];
        let mut supb = vec![0.0; 4];

        let iter = if on_first { self.recdomain1 } else { self.recdomain2 };

        let mut lastpt2d = DVec2::ZERO;
        let mut pt2d = DVec2::ZERO;
        let mut ponarc = 0.0;
        *indexsol = self.arc_to_recadre(on_first, the_sol, 0, &mut lastpt2d, &mut pt2d, &mut ponarc);
        pmin = ponarc;
        *is_vtx = false;
        if *indexsol == 0 {
            return false;
        }

        topol_tool_init(iter);
        let mut nbarc = 1;
        while nbarc < *indexsol {
            nbarc += 1;
            topol_tool_next(iter);
        }

        let mut thearc = topol_tool_value(iter);

        let mut thecur = if on_first {
            blend_tool_curve_on_surf(&thearc, self.surf1)
        } else {
            blend_tool_curve_on_surf(&thearc, self.surf2)
        };

        // Le probleme a resoudre
        f_inv.set_curve_on_surface(on_first, &thecur);
        f_inv.get_bounds(&mut infb, &mut supb);
        infb[1] -= extrap;
        supb[1] += extrap;

        f_inv.get_tolerance(&mut toler, 0.1 * self.tolpoint3d); // Il vaut mieux garder un peu de marge
        // OCCT: math_FunctionSetRoot rsnld(FuncInv, toler, 35);
        let fi: &dyn FunctionSetWithDerivatives = f_inv;
        let mut rsnld = FunctionSetRoot::new(fi, &toler, 35);
        // OCCT: toler *= 10; — Mais on fait les tests correctements
        for t in toler.iter_mut() {
            *t *= 10.0;
        }

        // Calcul d'un point d'init
        let mut ufirst = 0.0;
        let mut ulast = 0.0;
        blend_tool_bounds(&thecur, &mut ufirst, &mut ulast);
        // Pour aider a trouver les coins singuliers on recadre eventuellement
        // le parametre
        if (pmin - ufirst).abs() < (ulast - ufirst).abs() / 1000.0 {
            pmin = ufirst;
        }
        if (pmin - ulast).abs() < (ulast - ufirst).abs() / 1000.0 {
            pmin = ulast;
        }

        if byinter {
            let last_param = self.previous_p.parameter();
            // Verifie que le recadrage n'est pas un jalons
            if !self.jalons.is_empty() {
                let t1;
                let t2;
                let mut cherche = true;
                let mut ii = 0usize;
                if last_param < self.param {
                    t1 = last_param;
                    t2 = self.param;
                } else {
                    t1 = self.param;
                    t2 = last_param;
                }
                while ii < self.jalons.len() && cherche {
                    let t = self.jalons[ii].parameter();
                    if ((t1 < t) && (t2 > t)) || (t == self.param) {
                        jalons_trouve = true;
                        le_jalon = ii;
                        cherche = false; // Ne marche que si l'on sort simultanement
                    } else {
                        cherche = t < t2; // On s'arrete si t>=t2;
                    }
                    ii += 1;
                }
            }
            if !jalons_trouve {
                // Initialisation par Interpolation
                let pnt1;
                let pnt2;
                // OCCT: thecur->D0(pmin, Pnt);
                let pnt = thecur.point_at(pmin);
                if on_first {
                    let (u, v) = self.previous_p.parameters_on_s2();
                    pnt1 = DVec2::new(u, v);
                    pnt2 = DVec2::new(the_sol[2], the_sol[3]);
                } else {
                    let (u, v) = self.previous_p.parameters_on_s1();
                    pnt1 = DVec2::new(u, v);
                    pnt2 = DVec2::new(the_sol[0], the_sol[1]);
                }

                let mut lambda = pnt.distance(lastpt2d);
                if lambda > 1.0e-12 {
                    lambda /= pnt.distance(lastpt2d) + pnt.distance(pt2d);
                } else {
                    lambda = 0.0;
                }
                solrst[0] = pmin;
                solrst[1] = (1.0 - lambda) * last_param + lambda * self.param;
                solrst[2] = (1.0 - lambda) * pnt1.x + lambda * pnt2.x;
                solrst[3] = (1.0 - lambda) * pnt1.y + lambda * pnt2.y;
            }
        } else {
            // sinon on initialise par le dernier point calcule
            solrst[0] = pmin;
            solrst[1] = self.param;
            if on_first {
                solrst[2] = the_sol[2];
                solrst[3] = the_sol[3];
            } else {
                solrst[2] = the_sol[0];
                solrst[3] = the_sol[1];
            }
        }

        if jalons_trouve {
            // On recupere le jalon
            let mon_jalon = self.jalons[le_jalon].clone();
            let periodic;
            let mut uperiod = 0.0;
            let mut vperiod = 0.0;
            solrst[1] = mon_jalon.parameter();
            if on_first {
                let (u, v) = mon_jalon.parameters_on_s2();
                solrst[2] = u;
                solrst[3] = v;
                periodic = self.surf2.is_u_periodic() || self.surf2.is_v_periodic();
            } else {
                let (u, v) = mon_jalon.parameters_on_s1();
                solrst[2] = u;
                solrst[3] = v;
                periodic = self.surf1.is_u_periodic() || self.surf1.is_v_periodic();
            }

            // Recadrage eventuelle pour le cas periodique
            if periodic {
                let surf = if on_first { self.surf2 } else { self.surf1 };

                lastpt2d = thecur.point_at(pmin);

                if surf.is_u_periodic() {
                    uperiod = hsurface_tool_u_period(surf);
                    if solrst[2] - lastpt2d.x > uperiod * 0.6 {
                        solrst[2] -= uperiod;
                    }
                    if solrst[2] - lastpt2d.x < -uperiod * 0.6 {
                        solrst[2] += uperiod;
                    }
                }
                if surf.is_v_periodic() {
                    vperiod = hsurface_tool_v_period(surf);
                    if solrst[3] - lastpt2d.y > vperiod * 0.6 {
                        solrst[3] -= vperiod;
                    }
                    if solrst[3] - lastpt2d.y < -vperiod * 0.6 {
                        solrst[3] += vperiod;
                    }
                }
            }

            // Pour le parametre sur arc il faut projeter...
            pt2d = DVec2::new(solrst[2], solrst[3]);
            let pnt = thecur.point_at(ufirst);
            dist = pt2d.distance(pnt);
            solrst[0] = ufirst;
            let pnt = thecur.point_at(ulast);
            let distaux = pt2d.distance(pnt);
            if distaux < dist {
                solrst[0] = ulast;
                dist = distaux;
            }

            if dist > PCONFUSION {
                prm = pmin;
                let (ok, prm, distaux) = if on_first {
                    blend_tool_project(&pt2d, self.surf1, &thearc)
                } else {
                    blend_tool_project(&pt2d, self.surf2, &thearc)
                };
                let _ = distaux;
                if ok && pt2d.distance(thecur.point_at(prm)) < dist {
                    solrst[0] = prm;
                } else {
                    solrst[0] = pmin;
                }
            }
            // On verifie le jalon
            jalons_trouve = f_inv.is_solution(solrst, self.tolpoint3d);
        }

        if !jalons_trouve {
            // Resolution...
            let fi: &mut dyn FunctionSetWithDerivatives = f_inv;
            rsnld.perform(fi, solrst, &infb, &supb, false);
            if !rsnld.is_done() {
                recadre = false;
            } else {
                let root = rsnld.root();
                solrst.copy_from_slice(&root);
                recadre = f_inv.is_solution(solrst, self.tolpoint3d);
            }
        }

        // En cas d'echecs, on regarde si un autre arc
        // peut faire l'affaire (cas des sorties a proximite d'un vertex)
        dist = (ulast - ufirst) / 100.0;
        if !recadre && ((pmin - ulast).abs() < dist || (pmin - ufirst).abs() < dist) {
            *indexsol = self.arc_to_recadre(on_first, the_sol, *indexsol, &mut lastpt2d, &mut pt2d, &mut pmin);
            if *indexsol == 0 {
                return false;
            }

            topol_tool_init(iter);
            let mut nbarc = 1;
            while nbarc < *indexsol {
                nbarc += 1;
                topol_tool_next(iter);
            }
            thearc = topol_tool_value(iter);

            thecur = if on_first {
                blend_tool_curve_on_surf(&thearc, self.surf1)
            } else {
                blend_tool_curve_on_surf(&thearc, self.surf2)
            };
            solrst[0] = pmin;
            // Le probleme a resoudre
            f_inv.set_curve_on_surface(on_first, &thecur);
            f_inv.get_bounds(&mut infb, &mut supb);
            f_inv.get_tolerance(&mut toler, 0.1 * self.tolpoint3d); // Il vaut mieux garder un peu de marge
            // OCCT: math_FunctionSetRoot aRsnld(FuncInv, toler, 35);
            let fi: &dyn FunctionSetWithDerivatives = f_inv;
            let mut a_rsnld = FunctionSetRoot::new(fi, &toler, 35);
            // OCCT: toler *= 10; — Mais on fait les tests correctements
            for t in toler.iter_mut() {
                *t *= 10.0;
            }
            // Resolution...
            let fi: &mut dyn FunctionSetWithDerivatives = f_inv;
            a_rsnld.perform(fi, solrst, &infb, &supb, false);

            if !a_rsnld.is_done() {
                recadre = false;
            } else {
                let root = a_rsnld.root();
                solrst.copy_from_slice(&root);
                recadre = f_inv.is_solution(solrst, self.tolpoint3d);
            }
        }

        if recadre {
            // Classification topologique
            thecur = if on_first {
                blend_tool_curve_on_surf(&thearc, self.surf1)
            } else {
                blend_tool_curve_on_surf(&thearc, self.surf2)
            };
            blend_tool_bounds(&thecur, &mut ufirst, &mut ulast);

            topol_tool_initialize(iter, &thearc);
            topol_tool_init_vertex_iterator(iter);
            *is_vtx = !topol_tool_more_vertex(iter);
            while !*is_vtx {
                *vtx = *topol_tool_vertex(iter);
                vtol = 0.4 * (ulast - ufirst).abs(); // Un majorant de la tolerance
                if vtol > blend_tool_tolerance(vtx, &thearc).max(toler[0]) {
                    vtol = blend_tool_tolerance(vtx, &thearc).max(toler[0]);
                }
                if (blend_tool_parameter(vtx, &thearc) - solrst[0]).abs() <= vtol {
                    *is_vtx = true; // On est dans la boule du vertex ou
                                     // le vertex est dans la "boule" du recadrage
                } else {
                    topol_tool_next_vertex(iter);
                    *is_vtx = !topol_tool_more_vertex(iter);
                }
            }
            if !topol_tool_more_vertex(iter) {
                *is_vtx = false;
            }
            return true;
        }
        false
    }

    /// OCCT Transition(OnFirst, A, Param, TLine, TArc)
    /// (BRepBlend_Walking.cxx L1632-1770).
    pub(crate) fn transition(
        &self,
        on_first: bool,
        a: &Curve2d,
        param: f64,
        t_line: &mut Transition,
        t_arc: &mut Transition,
    ) {
        let mut computetranstionaveclacorde = false;
        let tgline;
        let prevprev;

        if self.previous_p.is_tangency_point() {
            if self.line.nb_points() < 2 {
                return;
            }
            computetranstionaveclacorde = true;
            if self.sens < 0.0 {
                prevprev = self.line.point(2).clone();
            } else {
                prevprev = self.line.point(self.line.nb_points() - 1).clone();
            }
        } else {
            prevprev = BlendPoint::new();
        }

        let (p2d, dp2d) = hcurve2d_tool_d1(a, param);

        let (d1u, d1v, tgline_dir);
        if on_first {
            let (_pbid, d1u_s, d1v_s) = hsurface_tool_d1(self.surf1, p2d.x, p2d.y);
            if !computetranstionaveclacorde {
                tgline_dir = self.previous_p.tangent_on_s1();
            } else {
                tgline_dir = self.previous_p.point_on_s1() - prevprev.point_on_s1();
            }
            d1u = d1u_s;
            d1v = d1v_s;
        } else {
            let (_pbid, d1u_s, d1v_s) = hsurface_tool_d1(self.surf2, p2d.x, p2d.y);
            if !computetranstionaveclacorde {
                tgline_dir = self.previous_p.tangent_on_s2();
            } else {
                tgline_dir = self.previous_p.point_on_s2() - prevprev.point_on_s2();
            }
            d1u = d1u_s;
            d1v = d1v_s;
        }
        tgline = tgline_dir;

        // OCCT: tgrst.SetLinearForm(dp2d.X(), d1u, dp2d.Y(), d1v);
        let tgrst = d1u * dp2d.x + d1v * dp2d.y;

        // OCCT: CSLib::Normal(d1u, d1v, 1.e-9, stat, thenormal);
        let (thenormal_opt, stat) = normal_from_derivatives(d1u, d1v, 1.0e-9);
        let normale = if stat == DerivativeStatus::Done {
            thenormal_opt.expect("CSLib Normal defined")
        } else {
            let surf = if on_first { self.surf1 } else { self.surf2 };
            // OCCT: NCollection_Array2<gp_Vec> Der(0, 2, 0, 2); filled with
            // the D2 block and the DN(2,1) / DN(1,2) / DN(2,2) derivatives.
            let (pbid, d1u2, d1v2, d2u, d2v, d2uv) = hsurface_tool_d2(surf, p2d.x, p2d.y);
            let _ = pbid;
            let _ = (d1u2, d1v2);
            let mut der = vec![vec![DVec3::ZERO; 3]; 3];
            der[1][0] = d2u;
            der[0][1] = d2v;
            der[1][1] = d2uv;
            der[2][1] = surf.dn(p2d.x, p2d.y, 2, 1);
            der[1][2] = surf.dn(p2d.x, p2d.y, 1, 2);
            der[2][2] = surf.dn(p2d.x, p2d.y, 2, 2);
            let (umin, umax, vmin, vmax) = hsurface_tool_parameters(surf);
            // OCCT: CSLib::Normal(2, Der, 1.e-9, X, Y, UMin, UMax, VMin,
            // VMax, stat, thenormal, iu, iv);
            let (nstat, thenormal_opt, _iu, _iv) =
                normal_max_order(2, &der, 1.0e-9, p2d.x, p2d.y, umin, umax, vmin, vmax);
            let _ = nstat; // OCCT debug-prints InfinityOfSolutions only
            thenormal_opt.unwrap_or(DVec3::ZERO)
        };

        make_transition(tgline, tgrst, normale, t_line, t_arc);
    }

    /// OCCT MakeExtremity(Extrem, OnFirst, Index, Param, IsVtx, Vtx)
    /// (BRepBlend_Walking.cxx L1772-1810).
    pub(crate) fn make_extremity(
        &mut self,
        extrem: &mut BRepBlendExtremity,
        on_first: bool,
        index: i32,
        param: f64,
        is_vtx: bool,
        vtx: &HVertex,
    ) {
        let mut tline = Transition::new();
        let mut tarc = Transition::new();

        if on_first {
            extrem.set_value(
                self.previous_p.point_on_s1(),
                self.sol[0],
                self.sol[1],
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_s1());
            }
        } else {
            extrem.set_value(
                self.previous_p.point_on_s2(),
                self.sol[2],
                self.sol[3],
                self.previous_p.parameter(),
                self.tolpoint3d,
            );
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_s2());
            }
        }

        let iter = if on_first { self.recdomain1 } else { self.recdomain2 };

        topol_tool_init(iter);
        let mut nbarc = 1;

        while nbarc < index {
            nbarc += 1;
            topol_tool_next(iter);
        }

        let arc = topol_tool_value(iter);
        self.transition(on_first, &arc, param, &mut tline, &mut tarc);
        extrem.add_arc(&arc, param, tline, tarc);
        if is_vtx {
            extrem.set_vertex(vtx);
        }
    }

    /// OCCT MakeSingularExtremity(Extrem, OnFirst, Vtx)
    /// (BRepBlend_Walking.cxx L1812-1855).
    pub(crate) fn make_singular_extremity(
        &mut self,
        extrem: &mut BRepBlendExtremity,
        on_first: bool,
        vtx: &HVertex,
    ) {
        let mut tline;
        let mut tarc;
        let mut prm;

        let iter = if on_first { self.recdomain1 } else { self.recdomain2 };
        if on_first {
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_s1());
            }
        } else {
            if !self.previous_p.is_tangency_point() {
                extrem.set_tangent(self.previous_p.tangent_on_s2());
            }
        }

        topol_tool_init(iter);
        extrem.set_vertex(vtx);
        while topol_tool_more(iter) {
            let arc = topol_tool_value(iter);
            topol_tool_initialize(iter, &arc);
            topol_tool_init_vertex_iterator(iter);
            while topol_tool_more_vertex(iter) {
                if topol_tool_identical(iter, vtx, topol_tool_vertex(iter)) {
                    prm = blend_tool_parameter(vtx, &arc);
                    tline = Transition::new();
                    tarc = Transition::new();
                    self.transition(on_first, &arc, prm, &mut tline, &mut tarc);
                    extrem.add_arc(&arc, prm, tline, tarc);
                }
                topol_tool_next_vertex(iter);
            }
            topol_tool_next(iter);
        }
    }
}


impl BRepBlendWalking<'_> {
    /// OCCT InternalPerform(Func, FuncInv, Bound) (BRepBlend_Walking.cxx
    /// L1932-2590).
    pub(crate) fn internal_perform(
        &mut self,
        f: &mut dyn BlendFunction,
        f_inv: &mut dyn BlendFuncInv,
        bound: f64,
    ) {
        let mut stepw = self.pasmax;
        let nbp = self.line.nb_points();
        if nbp >= 2 {
            // On reprend le dernier step s'il n est pas trop petit.
            if self.sens < 0.0 {
                stepw = self.line.point(2).parameter() - self.line.point(1).parameter();
            } else {
                stepw = self.line.point(nbp).parameter() - self.line.point(nbp - 1).parameter();
            }
            stepw = stepw.max(100.0 * self.tolgui);
        }
        let mut parprec = self.param;
        // OCCT: hguide->D1(parprec, PtOnGuide, TgOnGuide);
        let mut pt_on_guide = self.hguide.point_at(parprec);
        let mut tg_on_guide = self.hguide.derivative_at(parprec);
        let mut prev_tg_on_guide = tg_on_guide;

        if self.sens * (parprec - bound) >= -self.tolgui {
            return;
        }
        let mut state = BlendStatus::OnRst12;
        // OCCT: TopAbs_State situ1 = TopAbs_IN, situ2 = TopAbs_IN;
        let mut situ1 = TopAbsState::In;
        let mut situ2 = TopAbsState::In;
        let mut w1 = 0.0f64;
        let mut w2 = 0.0f64;
        let mut index1 = 0i32;
        let mut index2 = 0i32;
        let mut nbarc;
        let mut arrrive;
        let mut recad1;
        let mut recad2;
        let mut control;
        let mut isvtx1 = false;
        let mut isvtx2 = false;
        let mut echecrecad;
        let mut p2d = DVec2::ZERO;
        let mut tolerance = vec![0.0; 4];
        let mut infbound = vec![0.0; 4];
        let mut supbound = vec![0.0; 4];
        let mut parinit = vec![0.0; 4];
        let mut solrst1 = vec![0.0; 4];
        let mut solrst2 = vec![0.0; 4];
        let mut vtx1 = HVertex::new();
        let mut vtx2 = HVertex::new();
        let mut ext1 = BRepBlendExtremity::new();
        let mut ext2 = BRepBlendExtremity::new();

        // IntSurf_Transition Tline,Tarc;

        f.get_tolerance(&mut tolerance, self.tolpoint3d);
        f.get_bounds(&mut infbound, &mut supbound);

        // OCCT: math_FunctionSetRoot rsnld(Func, tolerance, 30);
        let fswd: &dyn FunctionSetWithDerivatives = &*f;
        let mut rsnld = FunctionSetRoot::new(fswd, &tolerance, 30);
        parinit.copy_from_slice(&self.sol);

        arrrive = false;
        self.param = parprec + self.sens * stepw;
        if self.sens * (self.param - bound) > 0.0 {
            stepw = self.sens * (bound - parprec) * 0.5;
            self.param = parprec + self.sens * stepw;
        }

        evalpinit(
            &mut parinit,
            &self.previous_p,
            parprec,
            self.param,
            &infbound,
            &supbound,
            self.clason_s1,
            self.clason_s2,
        );

        while !arrrive {
            // OCCT: hguide->D1(param, PtOnGuide, TgOnGuide);
            pt_on_guide = self.hguide.point_at(self.param);
            tg_on_guide = self.hguide.derivative_at(self.param);
            // Check deflection on guide
            let cosi = prev_tg_on_guide.dot(tg_on_guide);
            let cosi2;
            if cosi < GP_RESOLUTION {
                // angle>=pi/2 or null magnitude
                cosi2 = 0.0;
            } else {
                cosi2 = cosi * cosi / prev_tg_on_guide.length_squared()
                    / tg_on_guide.length_squared();
            }
            if cosi2 < COS_REF_3D {
                // angle 3d too great
                state = BlendStatus::StepTooLarge;
                stepw = stepw / 2.0;
                self.param = parprec + self.sens * stepw; // on ne risque pas de depasser Bound.
                if stepw.abs() < self.tolgui {
                    ext1.set_value(
                        self.previous_p.point_on_s1(),
                        self.sol[0],
                        self.sol[1],
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    ext2.set_value(
                        self.previous_p.point_on_s2(),
                        self.sol[2],
                        self.sol[3],
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    if !self.previous_p.is_tangency_point() {
                        ext1.set_tangent(self.previous_p.tangent_on_s1());
                        ext2.set_tangent(self.previous_p.tangent_on_s2());
                    }
                    arrrive = true;
                }
                continue;
            }
            prev_tg_on_guide = tg_on_guide;
            //////////////////////////

            let mut bonpoint = true;
            f.set_param(self.param);
            let fs: &mut dyn FunctionSetWithDerivatives = f;
            rsnld.perform(fs, &parinit, &infbound, &supbound, false);

            if !rsnld.is_done() {
                state = BlendStatus::StepTooLarge;
                bonpoint = false;
            } else {
                let root = rsnld.root();
                self.sol.copy_from_slice(&root);

                if self.clason_s1 {
                    situ1 = self.domain1.classify(
                        DVec2::new(self.sol[0], self.sol[1]),
                        tolerance[0].min(tolerance[1]),
                        false,
                    );
                } else {
                    situ1 = TopAbsState::In;
                }
                if self.clason_s2 {
                    situ2 = self.domain2.classify(
                        DVec2::new(self.sol[2], self.sol[3]),
                        tolerance[2].min(tolerance[3]),
                        false,
                    );
                } else {
                    situ2 = TopAbsState::In;
                }
            }
            if bonpoint && self.line.nb_points() == 1 && (situ1 != TopAbsState::In || situ2 != TopAbsState::In) {
                state = BlendStatus::StepTooLarge;
                bonpoint = false;
            }
            if bonpoint {
                w1 = bound;
                w2 = bound;
                recad1 = false;
                recad2 = false;
                echecrecad = false;
                control = false;

                // OCCT passes the member sol by const reference into
                // Recadre; rcad splits the &mut self borrow with a copy.
                let sol_c = self.sol.clone();
                if situ1 == TopAbsState::Out || situ1 == TopAbsState::On {
                    // pb inverse sur surf1
                    // Si le recadrage s'effectue dans le sens de la progression
                    // a une tolerance pres, on a pris la mauvaise solution.
                    recad1 = self.recadre(
                        f_inv,
                        true,
                        &sol_c,
                        &mut solrst1,
                        &mut index1,
                        &mut isvtx1,
                        &mut vtx1,
                        0.0,
                    );

                    if recad1 {
                        let wtemp = solrst1[1];
                        if (self.param - wtemp) / self.sens >= -10.0 * self.tolgui {
                            w1 = solrst1[1];
                            control = true;
                        } else {
                            echecrecad = true;
                            recad1 = false;
                            state = BlendStatus::StepTooLarge;
                            bonpoint = false;
                            stepw = stepw / 2.0;
                        }
                    } else {
                        echecrecad = true;
                    }
                }
                if situ2 == TopAbsState::Out || situ2 == TopAbsState::On {
                    // pb inverse sur surf2
                    recad2 = self.recadre(
                        f_inv,
                        false,
                        &sol_c,
                        &mut solrst2,
                        &mut index2,
                        &mut isvtx2,
                        &mut vtx2,
                        0.0,
                    );

                    if recad2 {
                        let wtemp = solrst2[1];
                        if (self.param - wtemp) / self.sens >= -10.0 * self.tolgui {
                            w2 = solrst2[1];
                            control = true;
                        } else {
                            echecrecad = true;
                            recad2 = false;
                            state = BlendStatus::StepTooLarge;
                            bonpoint = false;
                            stepw = stepw / 2.0;
                        }
                    } else {
                        echecrecad = true;
                    }
                }

                // Que faut il controler
                if recad1 && recad2 {
                    if (w1 - w2).abs() <= 10.0 * self.tolgui {
                        // pas besoin de controler les recadrage
                        control = false;
                    } else if self.sens * (w1 - w2) < 0.0 {
                        // sol sur 1 ?
                        recad2 = false;
                    } else {
                        // sol sur 2 ?
                        recad1 = false;
                    }
                }

                // Controle effectif des recadrage
                if control {
                    let situ;
                    if recad1 && self.clason_s2 {
                        situ = self.recdomain2.classify(
                            DVec2::new(solrst1[2], solrst1[3]),
                            tolerance[2].min(tolerance[3]),
                            false,
                        );
                        if situ == TopAbsState::Out {
                            recad1 = false;
                            echecrecad = true;
                        }
                    } else if recad2 && self.clason_s1 {
                        situ = self.recdomain1.classify(
                            DVec2::new(solrst2[2], solrst2[3]),
                            tolerance[0].min(tolerance[1]),
                            false,
                        );
                        if situ == TopAbsState::Out {
                            recad2 = false;
                            echecrecad = true;
                        }
                    }
                }

                if recad1 || recad2 {
                    echecrecad = false;
                }

                if !echecrecad {
                    if recad1 && recad2 {
                        // sol sur 1 et 2 a la fois
                        //  On passe par les arcs , pour ne pas avoir de probleme
                        //  avec les surfaces periodiques.
                        state = BlendStatus::OnRst12;
                        self.param = (w1 + w2) / 2.0;
                        p2d = hcurve2d_tool_value(&topol_tool_value(self.recdomain1), solrst1[0]);
                        self.sol[0] = p2d.x;
                        self.sol[1] = p2d.y;
                        let pnt1 = hsurface_tool_value(self.surf1, self.sol[0], self.sol[1]);
                        p2d = hcurve2d_tool_value(&topol_tool_value(self.recdomain2), solrst2[0]);
                        self.sol[2] = p2d.x;
                        self.sol[3] = p2d.y;
                        let pnt2 = hsurface_tool_value(self.surf2, self.sol[2], self.sol[3]);
                        let tol_prod = 1.0e-5;
                        let mut saved_params = [0.0f64; 2];
                        let mut same_dirs = [false; 2];
                        // OCCT: theElSpine.GetSavedFirstParameter() /
                        // GetSavedLastParameter() — the OCCT constructor
                        // default is Precision::Infinite().
                        saved_params[0] = elspine_get_saved_first_parameter(self.hguide);
                        saved_params[1] = elspine_get_saved_last_parameter(self.hguide);
                        for ind in 0..2 {
                            // OCCT BRepBlend_Walking.cxx L2198:
                            // !Precision::IsInfinite(SavedParams[ind]).
                            if !rcad_kernel::precision::is_infinite_value(saved_params[ind]) {
                                // Check the original first and last parameters
                                // of guide curve for equality to found
                                // parameter <param>:
                                let pnt0 = self.hguide.point_at(saved_params[ind]);
                                let mut dir0 = self.hguide.derivative_at(saved_params[ind]);
                                let length = dir0.length();
                                if length <= GP_RESOLUTION {
                                    continue;
                                }
                                dir0 /= length;
                                let plane = gce_make_pln(pnt0, pnt1, pnt2);
                                if !plane.is_done() {
                                    continue;
                                }
                                let dir_plane = plane.axis_direction();
                                let the_prod = dir0.cross(dir_plane);
                                let prod_mod = the_prod.length();
                                if prod_mod <= tol_prod {
                                    same_dirs[ind] = true;
                                }
                            }
                        }
                        // OCCT Precision::Infinite() (Precision.hxx L350-353).
                        let mut the_param = rcad_kernel::core::precision::INFINITE_VALUE;
                        // Choose the closest parameter
                        if same_dirs[0] && same_dirs[1] {
                            the_param = if (self.param - saved_params[0]).abs()
                                < (self.param - saved_params[1]).abs()
                            {
                                saved_params[0]
                            } else {
                                saved_params[1]
                            };
                        } else if same_dirs[0] {
                            the_param = saved_params[0];
                        } else if same_dirs[1] {
                            the_param = saved_params[1];
                        }

                        let mut new_u = 0.0;
                        let mut new_v = 0.0;
                        let mut new_param = 0.0;
                        let mut new_pnt = DVec3::ZERO;
                        let corrected = self.correct_extremity_on_one_rst(
                            1,
                            self.sol[2],
                            self.sol[3],
                            self.param,
                            pnt1,
                            &mut new_u,
                            &mut new_v,
                            &mut new_pnt,
                            &mut new_param,
                        );
                        if corrected && (self.param - new_param).abs() < (self.param - the_param).abs() {
                            the_param = new_param;
                        }

                        // OCCT BRepBlend_Walking.cxx L2259:
                        // !Precision::IsInfinite(theParam).
                        if !rcad_kernel::precision::is_infinite_value(the_param) {
                            self.param = the_param;
                        }
                    } else if recad1 {
                        // sol sur 1
                        state = BlendStatus::OnRst1;
                        self.param = w1;
                        topol_tool_init(self.recdomain1);
                        nbarc = 1;
                        while nbarc < index1 {
                            nbarc += 1;
                            topol_tool_next(self.recdomain1);
                        }
                        p2d = hcurve2d_tool_value(&topol_tool_value(self.recdomain1), solrst1[0]);
                        self.sol[0] = p2d.x;
                        self.sol[1] = p2d.y;
                        self.sol[2] = solrst1[2];
                        self.sol[3] = solrst1[3];
                        let the_pnt_on_rst = hsurface_tool_value(self.surf1, self.sol[0], self.sol[1]);
                        let mut new_u = 0.0;
                        let mut new_v = 0.0;
                        let mut new_param = 0.0;
                        let mut new_pnt = DVec3::ZERO;
                        let corrected = self.correct_extremity_on_one_rst(
                            1,
                            self.sol[2],
                            self.sol[3],
                            self.param,
                            the_pnt_on_rst,
                            &mut new_u,
                            &mut new_v,
                            &mut new_pnt,
                            &mut new_param,
                        );
                        if corrected {
                            self.param = new_param;
                            self.sol[2] = new_u;
                            self.sol[3] = new_v;
                        }
                    } else if recad2 {
                        // sol sur 2
                        state = BlendStatus::OnRst2;
                        self.param = w2;

                        topol_tool_init(self.recdomain2);
                        nbarc = 1;
                        while nbarc < index2 {
                            nbarc += 1;
                            topol_tool_next(self.recdomain2);
                        }
                        p2d = hcurve2d_tool_value(&topol_tool_value(self.recdomain2), solrst2[0]);
                        self.sol[0] = solrst2[2];
                        self.sol[1] = solrst2[3];
                        self.sol[2] = p2d.x;
                        self.sol[3] = p2d.y;
                        let the_pnt_on_rst = hsurface_tool_value(self.surf2, self.sol[2], self.sol[3]);
                        let mut new_u = 0.0;
                        let mut new_v = 0.0;
                        let mut new_param = 0.0;
                        let mut new_pnt = DVec3::ZERO;
                        let corrected = self.correct_extremity_on_one_rst(
                            2,
                            self.sol[0],
                            self.sol[1],
                            self.param,
                            the_pnt_on_rst,
                            &mut new_u,
                            &mut new_v,
                            &mut new_pnt,
                            &mut new_param,
                        );
                        if corrected {
                            self.param = new_param;
                            self.sol[0] = new_u;
                            self.sol[1] = new_v;
                        }
                    } else {
                        state = BlendStatus::Ok;
                    }

                    let testdefl = true;
                    if recad1 || recad2 {
                        f.set_param(self.param);
                        // Il vaut mieux un pas non orthodoxe que pas de recadrage!! PMN
                        state = self.test_arret(
                            f,
                            state,
                            testdefl && (stepw.abs() > 3.0 * self.tolgui),
                            false,
                            true,
                        );
                    } else {
                        state = self.test_arret(f, state, testdefl, true, false);
                    }
                } else {
                    // Ou bien le pas max est mal regle. On divise.
                    if stepw > 2.0 * self.tolgui {
                        state = BlendStatus::StepTooLarge;
                        // Sinon echec recadrage. On sort avec PointsConfondus
                    } else {
                        state = BlendStatus::SamePoints;
                    }
                }
            }

            match state {
                BlendStatus::Ok => {
                    // Mettre a jour la ligne.
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }

                    parprec = self.param;

                    if self.param == bound {
                        arrrive = true;
                        ext1.set_value(
                            self.previous_p.point_on_s1(),
                            self.sol[0],
                            self.sol[1],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        ext2.set_value(
                            self.previous_p.point_on_s2(),
                            self.sol[2],
                            self.sol[3],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        if !self.previous_p.is_tangency_point() {
                            ext1.set_tangent(self.previous_p.tangent_on_s1());
                            ext2.set_tangent(self.previous_p.tangent_on_s2());
                        }

                        // Indiquer que fin sur Bound.
                    } else {
                        self.param = self.param + self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                    evalpinit(
                        &mut parinit,
                        &self.previous_p,
                        parprec,
                        self.param,
                        &infbound,
                        &supbound,
                        self.clason_s1,
                        self.clason_s2,
                    );
                }

                BlendStatus::StepTooLarge => {
                    stepw = stepw / 2.0;
                    if stepw.abs() < self.tolgui {
                        ext1.set_value(
                            self.previous_p.point_on_s1(),
                            self.sol[0],
                            self.sol[1],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        ext2.set_value(
                            self.previous_p.point_on_s2(),
                            self.sol[2],
                            self.sol[3],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        if !self.previous_p.is_tangency_point() {
                            ext1.set_tangent(self.previous_p.tangent_on_s1());
                            ext2.set_tangent(self.previous_p.tangent_on_s2());
                        }
                        arrrive = true;
                        if self.line.nb_points() >= 2 {
                            // Indiquer qu on s arrete en cours de cheminement
                        }
                    } else {
                        self.param = parprec + self.sens * stepw; // on ne risque pas de depasser Bound.
                        evalpinit(
                            &mut parinit,
                            &self.previous_p,
                            parprec,
                            self.param,
                            &infbound,
                            &supbound,
                            self.clason_s1,
                            self.clason_s2,
                        );
                    }
                }

                BlendStatus::StepTooSmall => {
                    // Mettre a jour la ligne.
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }

                    parprec = self.param;

                    stepw = (1.5 * stepw).min(self.pasmax);
                    if self.param == bound {
                        arrrive = true;
                        ext1.set_value(
                            self.previous_p.point_on_s1(),
                            self.sol[0],
                            self.sol[1],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        ext2.set_value(
                            self.previous_p.point_on_s2(),
                            self.sol[2],
                            self.sol[3],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                        if !self.previous_p.is_tangency_point() {
                            ext1.set_tangent(self.previous_p.tangent_on_s1());
                            ext2.set_tangent(self.previous_p.tangent_on_s2());
                        }
                        // Indiquer que fin sur Bound.
                    } else {
                        self.param = self.param + self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                    evalpinit(
                        &mut parinit,
                        &self.previous_p,
                        parprec,
                        self.param,
                        &infbound,
                        &supbound,
                        self.clason_s1,
                        self.clason_s2,
                    );
                }

                BlendStatus::OnRst1 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }

                    self.make_extremity(&mut ext1, true, index1, solrst1[0], isvtx1, &vtx1);
                    // On blinde le cas singulier ou un des recadrage a planter
                    if self
                        .previous_p
                        .point_on_s1()
                        .distance(self.previous_p.point_on_s2())
                        <= 2.0 * self.tolpoint3d
                    {
                        // OCCT: Ext2.SetValue(P, U, V, Tol) — 4-argument
                        // overload, Param = 0.
                        ext2.set_value(
                            self.previous_p.point_on_s1(),
                            self.sol[2],
                            self.sol[3],
                            0.0,
                            self.tolpoint3d,
                        );
                        if isvtx1 {
                            self.make_singular_extremity(&mut ext2, false, &vtx1);
                        }
                    } else {
                        ext2.set_value(
                            self.previous_p.point_on_s2(),
                            self.sol[2],
                            self.sol[3],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                    }
                    arrrive = true;
                }

                BlendStatus::OnRst2 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }

                    // On blinde le cas singulier ou un des recadrage a plante
                    if self
                        .previous_p
                        .point_on_s1()
                        .distance(self.previous_p.point_on_s2())
                        <= 2.0 * self.tolpoint3d
                    {
                        // OCCT: Ext1.SetValue(P, U, V, Tol) — 4-argument
                        // overload, Param = 0.
                        ext1.set_value(
                            self.previous_p.point_on_s2(),
                            self.sol[0],
                            self.sol[1],
                            0.0,
                            self.tolpoint3d,
                        );
                        if isvtx2 {
                            self.make_singular_extremity(&mut ext1, true, &vtx2);
                        }
                    } else {
                        ext1.set_value(
                            self.previous_p.point_on_s1(),
                            self.sol[0],
                            self.sol[1],
                            self.previous_p.parameter(),
                            self.tolpoint3d,
                        );
                    }
                    self.make_extremity(&mut ext2, false, index2, solrst2[0], isvtx2, &vtx2);
                    arrrive = true;
                }

                BlendStatus::OnRst12 => {
                    if self.sens > 0.0 {
                        self.line.append(self.previous_p.clone());
                    } else {
                        self.line.prepend(self.previous_p.clone());
                    }

                    if (isvtx1 != isvtx2)
                        && self
                            .previous_p
                            .point_on_s1()
                            .distance(self.previous_p.point_on_s2())
                            <= 2.0 * self.tolpoint3d
                    {
                        // On blinde le cas singulier ou un seul recadrage
                        // est reconnu comme vertex.
                        if isvtx1 {
                            isvtx2 = true;
                            vtx2 = vtx1;
                        } else {
                            isvtx1 = true;
                            vtx1 = vtx2;
                        }
                    }

                    self.make_extremity(&mut ext1, true, index1, solrst1[0], isvtx1, &vtx1);
                    self.make_extremity(&mut ext2, false, index2, solrst2[0], isvtx2, &vtx2);
                    arrrive = true;
                }

                BlendStatus::SamePoints => {
                    // On arrete
                    ext1.set_value(
                        self.previous_p.point_on_s1(),
                        self.sol[0],
                        self.sol[1],
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    ext2.set_value(
                        self.previous_p.point_on_s2(),
                        self.sol[2],
                        self.sol[3],
                        self.previous_p.parameter(),
                        self.tolpoint3d,
                    );
                    if !self.previous_p.is_tangency_point() {
                        ext1.set_tangent(self.previous_p.tangent_on_s1());
                        ext2.set_tangent(self.previous_p.tangent_on_s2());
                    }
                    arrrive = true;
                }
                BlendStatus::Backward => {
                    // OCCT `default: break;`
                }
            }
            if arrrive {
                if self.sens > 0.0 {
                    self.line.set_end_points(&ext1, &ext2);
                } else {
                    self.line.set_start_points(&ext1, &ext2);
                }
            }
        }
        // OCCT keeps PtOnGuide for the debug traces; rcad has none.
        let _ = pt_on_guide;
    }

    /// OCCT CorrectExtremityOnOneRst(IndexOfRst, theU, theV, theParam,
    /// thePntOnRst, NewU, NewV, NewPoint, NewParam) const
    /// (BRepBlend_Walking.cxx L2592-2769).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn correct_extremity_on_one_rst(
        &self,
        index_of_rst: i32,
        the_u: f64,
        the_v: f64,
        the_param: f64,
        the_pnt_on_rst: DVec3,
        new_u: &mut f64,
        new_v: &mut f64,
        new_point: &mut DVec3,
        new_param: &mut f64,
    ) -> bool {
        let tol_ang = 0.001; // bug OCC25701

        // OCCT: ChFiDS_ElSpine& theElSpine = *hguide; — the vertex list is
        // not carried by the rcad guide curve (pending); NbVertices
        // reports 0 and the OCCT early-out applies.
        if elspine_nb_vertices(self.hguide) == 0 {
            return false;
        }

        let domain_of_rst = if index_of_rst == 1 { self.recdomain1 } else { self.recdomain2 };
        let surf_of_rst = if index_of_rst == 1 { self.surf1 } else { self.surf2 };
        let another_surf = if index_of_rst == 1 { self.surf2 } else { self.surf1 };

        // Correct point on surface 2
        // First we find right <param>
        let arc = topol_tool_value(domain_of_rst);
        let mut ends = [0.0f64; 2];
        ends[0] = hcurve2d_tool_first_parameter(&arc);
        ends[1] = hcurve2d_tool_last_parameter(&arc);
        let mut global_min_sq_dist = f64::INFINITY; // Precision::Infinite()
        let mut param_on_guide = 0.0;
        let mut point_on_guide = DVec3::ZERO;
        for k in 0..2 {
            let p2d_on_end = hcurve2d_tool_value(&topol_tool_value(domain_of_rst), ends[k]);
            let pnt_on_end = hsurface_tool_value(surf_of_rst, p2d_on_end.x, p2d_on_end.y);
            // OCCT L2643: Extrema_ExtPC projoncurv(PntOnEnd, theElSpine) —
            // the two-arg ctor runs Initialize over the full adaptor range
            // with the default theTolF 1.0e-10, then Perform(PntOnEnd).
            let a_adaptor = GeomCurveAdaptor::new(self.hguide.clone());
            let a_tool = CurveToolHandle::for_curve3(self.hguide, &a_adaptor, &a_adaptor);
            let projoncurv = ExtremaExtPC::new_point_curve(pnt_on_end, &a_tool, 1.0e-10);
            if !projoncurv.is_done() {
                continue;
            }
            let mut min_sq_dist = f64::INFINITY;
            let mut imin = 0usize;
            for ind in 1..=projoncurv.nb_ext() {
                let a_sq_dist = projoncurv.square_distance(ind);
                if a_sq_dist < min_sq_dist {
                    min_sq_dist = a_sq_dist;
                    imin = ind;
                }
            }
            if min_sq_dist < global_min_sq_dist {
                global_min_sq_dist = min_sq_dist;
                param_on_guide = projoncurv.point(imin).param;
                point_on_guide = projoncurv.point(imin).point;
            }
        }
        *new_param = param_on_guide;
        if self.hguide.is_periodic() {
            let period = hcurve_tool_period(self.hguide);
            let sign = if *new_param < the_param { 1.0 } else { -1.0 };
            while (*new_param - the_param).abs() > period / 2.0 {
                *new_param += sign * period;
            }
        }

        // Second we find right point and tangent on guide
        global_min_sq_dist = f64::INFINITY;
        let mut the_ax1_location = DVec3::ZERO;
        let mut the_ax1_direction = DVec3::Z;
        for ind in 1..=elspine_nb_vertices(self.hguide) {
            // OCCT: theElSpine.VertexWithTangent(ind) -> gp_Ax1.
            let (location, direction) = elspine_vertex_with_tangent(self.hguide, ind);
            let a_pnt = location;
            let a_sq_dist = point_on_guide.distance_squared(a_pnt);
            if a_sq_dist < global_min_sq_dist {
                global_min_sq_dist = a_sq_dist;
                the_ax1_location = location;
                the_ax1_direction = direction;
            }
        }
        let pnt0 = the_ax1_location;
        let dir0 = the_ax1_direction;
        // Check new point: is it real solution?
        let old_pon_guide = self.hguide.point_at(the_param);
        let pnt_on_surf2 = hsurface_tool_value(another_surf, the_u, the_v); // old point
        let plane = gce_make_pln(the_pnt_on_rst, old_pon_guide, pnt_on_surf2);
        if !plane.is_done() {
            return false;
        }
        let old_dir = plane.axis_direction();
        let mut angle = gp_dir_angle(old_dir, dir0);
        if angle > std::f64::consts::PI / 2.0 {
            angle = std::f64::consts::PI - angle;
        }
        if angle > tol_ang {
            return false;
        }
        ///////////////////////////////////////
        // Project the point(theU,theV) on the plane(Pnt0,Dir0)
        let a_vec = pnt_on_surf2 - pnt0;
        let a_translation = dir0 * a_vec.dot(dir0);
        let pnt_on_plane = pnt_on_surf2 - a_translation;

        // Check new point again: does point on restriction belong to the plane?
        let plane = gce_make_pln(the_pnt_on_rst, pnt0, pnt_on_plane);
        if !plane.is_done() {
            return false;
        }
        let dir_of_new_plane = plane.axis_direction();
        let mut angle = gp_dir_angle(dir0, dir_of_new_plane);
        if angle > std::f64::consts::PI / 2.0 {
            angle = std::f64::consts::PI - angle;
        }
        if angle > tol_ang {
            return false;
        }
        ////////////////////////////////////////////////////////////////////////

        // Project the point <PntOnPlane> on the surface 2
        // OCCT: Extrema_ExtPS projonsurf(PntOnPlane, *AnotherSurf,
        // PConfusion, PConfusion, Extrema_ExtFlag_MIN);
        let projonsurf = ExtPS::new(pnt_on_plane, another_surf, PCONFUSION, PCONFUSION);
        if projonsurf.is_done() {
            let mut min_sq_dist = f64::INFINITY;
            let mut imin = 0usize;
            for ind in 1..=projonsurf.nb_ext() {
                let a_sq_dist = projonsurf.square_distance(ind);
                if a_sq_dist < min_sq_dist {
                    min_sq_dist = a_sq_dist;
                    imin = ind;
                }
            }
            if imin != 0 {
                let new_pon_surf2 = projonsurf.point(imin);
                *new_point = new_pon_surf2.point;
                *new_u = new_pon_surf2.u;
                *new_v = new_pon_surf2.v;
                let uperiod = if another_surf.is_u_periodic() {
                    hsurface_tool_u_period(another_surf)
                } else {
                    0.0
                };
                let vperiod = if another_surf.is_v_periodic() {
                    hsurface_tool_v_period(another_surf)
                } else {
                    0.0
                };
                recadre_if_periodic(new_u, new_v, the_u, the_v, uperiod, vperiod);
                return true;
            }
        }

        false
    }
}

