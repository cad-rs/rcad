// OCCT Contap_TheSearch = IntStart_SearchOnBoundaries.gxx instantiated for
// Contap (TheVertex = handle<Adaptor3d_HVertex>, TheArc =
// handle<Adaptor2d_Curve2d>, TheArcTool = Contap_HCurve2dTool, TheSOBTool
// = Contap_HContTool, TheTopolTool = Adaptor3d_TopolTool, TheFunction =
// Contap_ArcFunction).
//
// IntStart_SearchOnBoundaries.gxx L68-1233 (the `_1.gxx` body).
// Instantiation notes:
//   - Contap_ArcFunction never sets its IntSurf_Quadric (the default is
//     GeomAbs_OtherSurface, IntSurf_Quadric.cxx L38-46), so the
//     quadric-exact branch guarded by `TypeQuad != GeomAbs_OtherSurface`
//     (gxx L374-473) is dead code for Contap: the engine always walks the
//     math_FunctionAllRoots path.  The guard is translated; the dead body
//     (Adaptor3d_CurveOnSurface + IntCurveSurface_HInter) stays deferred
//     to Stage 3b, like the other HLRBRep-facing deferrals.
//   - IsRegularity / TreatLC call Domain->Edge(): the base Adaptor3d
//     TopolTool has no edge, so both early-out (gxx L1005-1016, L1028-1048).
//     The BRepTopAdaptor branch is Stage 3b/3d (BRepAdaptor_Curve +
//     Extrema_ExtCC over real edges).

use rcad_kernel::math::opt::BrentMinimum;
use rcad_kernel::math::root::{FunctionAllRoots, FunctionSample, FunctionValue, FunctionWithDerivative};
use rcad_kernel::precision::{is_negative_infinite_value, is_positive_infinite_value, CONFUSION, PCONFUSION};

use crate::hlr::contap::arc_function::ArcFunction;
use crate::hlr::contap::domain::ContapDomain;
use crate::hlr::contap::the_path_point_of_the_search::ThePathPointOfTheSearch;
use crate::hlr::contap::the_segment_of_the_search::TheSegmentOfTheSearch;

use crate::hlr::contap::point::Arc;

/// OCCT gxx L105-127 — MinFunction: F(x)^2 as a math_Function.
struct MinFunction<'a> {
    my_func: &'a mut ArcFunction,
}

impl<'a> MinFunction<'a> {
    fn new(the_func: &'a mut ArcFunction) -> Self {
        MinFunction { my_func: the_func }
    }
}

impl FunctionValue for MinFunction<'_> {
    // returns value of the one-dimension-function when parameter
    // is equal to theX (gxx L113-120)
    fn value(&mut self, the_x: f64) -> Option<f64> {
        let the_fval = self.my_func.value(the_x)?;
        Some(the_fval * the_fval)
    }
}

/// OCCT gxx L191-227 — SolInfo, the sortable solution slot.
#[derive(Debug, Clone, Copy)]
struct SolInfo {
    my_math_index: i32,
    my_value: f64,
}

impl SolInfo {
    fn new() -> Self {
        SolInfo {
            my_math_index: -1,
            my_value: f64::MAX,
        }
    }
    /// OCCT Init(const math_FunctionAllRoots&, theIndex) (L200-204).
    fn init_roots(&mut self, the_solution: &FunctionAllRoots, the_index: i32) {
        self.my_math_index = the_index;
        self.my_value = the_solution.get_point(the_index as usize);
    }
    fn value(&self) -> f64 {
        self.my_value
    }
    fn index(&self) -> i32 {
        self.my_math_index
    }
    fn change_value(&mut self) -> &mut f64 {
        &mut self.my_value
    }
}

/// OCCT FindVertex (gxx L131-165).
fn find_vertex(
    a: &Arc,
    domain: &mut dyn ContapDomain,
    func: &mut ArcFunction,
    pnt: &mut Vec<ThePathPointOfTheSearch>,
    toler: f64,
) {
    // Find the vertex of the arc A restriction solutions. It stores
    // Vertex in the list solutions pnt.
    domain.initialize_arc(a);
    domain.init_vertex_iterator();
    while domain.more_vertex() {
        let vtx = domain.vertex();
        let param = crate::hlr::contap::h_cont_tool::parameter(vtx.as_ref(), a.as_ref());

        // Evaluate the function and look compared to tolerance of the
        // Vertex. If distance <= tolerance then add a vertex to the list of
        // solutions. The arc is already assumed in the load function.
        if let Some(valf) = func.value(param) {
            if valf.abs() <= toler {
                let itemp = func.get_state_number();
                pnt.push(ThePathPointOfTheSearch::with_vertex(
                    func.valpoint(itemp),
                    toler,
                    vtx,
                    a.clone(),
                    param,
                ));
                // Solution is added
            }
        }
        domain.next_vertex();
    }
}

/// OCCT IsDegenerated(const IntSurf_Quadric&) (gxx L178-189) — the arc
/// variant (L167-176) needs the Adaptor3d_CurveOnSurface and stays with
/// the dead quadric branch (Stage 3b).
fn is_degenerated(the_quadric: &crate::geomalgo::int_surf::quadric::Quadric) -> bool {
    let type_quad = the_quadric.type_quadric();
    if matches!(type_quad, crate::geomalgo::int_surf::quadric::QuadricType::Cone) {
        let a_cone = the_quadric.cone();
        let a_semi_angle = a_cone.half_angle_rad.abs();
        if a_semi_angle < 0.02 || a_semi_angle > 1.55 {
            return true;
        }
    }
    false
}

/// OCCT BoundedArc (gxx L229-798).
#[allow(clippy::too_many_arguments)]
fn bounded_arc(
    a: &Arc,
    domain: &mut dyn ContapDomain,
    pdeb: f64,
    pfin: f64,
    func: &mut ArcFunction,
    pnt: &mut Vec<ThePathPointOfTheSearch>,
    seg: &mut Vec<TheSegmentOfTheSearch>,
    tol_boundary: f64,
    tol_tangency: f64,
    arcsol: &mut bool,
    recheck_on_regularity: bool,
) {
    // Recherche des points solutions et des bouts d arc solution sur un
    // arc donne. On utilise la fonction math_FunctionAllRoots.
    let mut nbi = 0usize;
    let mut nbp = 0usize;

    let mut pardeb = 0.0;
    let mut parfin = 0.0;
    let mut ideb = 0i32;
    let mut ifin = 0i32;

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@ La Tolerance est asociee a l arc  ( Incoherence avec le cheminement )
    //@@@   ( EpsX ~ 1e-5   et ResolutionU et V ~ 1e-9 )
    //@@@   le vertex trouve ici n'est pas retrouve comme point d arret d une
    //@@@   ligne de cheminement
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    let mut eps_x = 1.0e-10;
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    //  int NbEchant = TheSOBTool::NbSamplesOnArc(A);
    let mut nb_echant = func.nb_samples();
    if nb_echant < 100 {
        nb_echant = 100; //-- lbr le 22 Avril 96
    }
    //-- Toujours des pbs

    //-- Modif 24  Aout 93 -----------------------------
    let mut n_tol_tangency = tol_tangency;
    if (pfin - pdeb) < (tol_tangency * 10.0) {
        n_tol_tangency = (pfin - pdeb) * 0.1;
    }
    if eps_x > (n_tol_tangency + n_tol_tangency) {
        eps_x = n_tol_tangency * 0.1;
    }
    //--------------------------------------------------

    //-------------------------------------------------------------- REJECTIONS le 15 oct 98
    let mut rejection = true;
    let mut maxdr;
    let mut maxr;
    let mut minr;
    let mut dur = (pfin - pdeb) * 0.2;
    {
        minr = f64::MAX;
        maxr = -minr;
        maxdr = -minr;
        let mut ur = pdeb;
        for _ in 1..=6 {
            // double F, D; if (Func.Values(ur, F, D))
            if let Some((f, mut d)) = func.values(ur) {
                if d < 0.0 {
                    d = -d;
                }
                d *= dur + dur;
                if d > maxdr {
                    maxdr = d;
                }
                let lminr = f - d;
                let lmaxr = f + d;
                if lminr < minr {
                    minr = lminr;
                }
                if lmaxr > maxr {
                    maxr = lmaxr;
                }
                if minr < 0.0 && maxr > 0.0 {
                    rejection = false;
                    break;
                }
            }
            ur += dur;
        }
        if rejection {
            dur = 0.001 + maxdr + (maxr - minr) * 0.1;
            minr -= dur;
            maxr += dur;
            if minr < 0.0 && maxr > 0.0 {
                rejection = false;
            }
        }
    }

    *arcsol = false;

    if !rejection {
        let a_quadric = func.quadric().clone();
        let type_quad = a_quadric.type_quadric();

        // IntCurveSurface_HInter IntCS; bool IsIntCSdone = false; — the
        // quadric-exact branch below is dead for Contap (the ArcFunction
        // quadric is never set, so TypeQuad stays OtherSurface); the
        // engine enters the math_FunctionAllRoots path (gxx L475).
        let mut is_int_cs_done = false;
        let _ = &mut is_int_cs_done;
        let mut _params: Vec<f64> = Vec::new(); // NCollection_Sequence<double> Params (dead branch)

        let mut p_sol: Option<FunctionAllRoots> = None;

        let echant = FunctionSample::new(pdeb, pfin, nb_echant);

        let mut maxdist;
        let mut aelargir = true;
        // modified by NIZNHY-PKV Thu Apr 12 09:25:19 2001 f
        maxdist = tol_boundary + tol_tangency;
        // modified by NIZNHY-PKV Thu Apr 12 09:25:23 2001 t
        for i in 1..=nb_echant {
            if !aelargir {
                break;
            }
            let u = echant.get_parameter(i);
            if let Some(dist) = func.value(u) {
                if dist > maxdist || -dist > maxdist {
                    aelargir = false;
                }
            }
        }
        if !(aelargir && maxdist < 0.01) {
            maxdist = tol_boundary;
        }

        // OCCT gxx L374: if (TypeQuad != GeomAbs_OtherSurface) — the exact
        // IntCurveSurface_HInter intersection for canonic CurveOnSurface.
        // Dead for Contap (see the module notes); TypeConS stays
        // GeomAbs_OtherCurve and IsIntCSdone stays false.
        let _type_con_s_is_other = true;
        let _ = is_degenerated(&a_quadric);
        let _ = &_params;

        if !is_int_cs_done {
            let mut sol = FunctionAllRoots::new(func, &echant, eps_x, maxdist, maxdist);
            if !sol.is_done() {
                panic!("Standard_Failure: IntStart_SearchOnBoundaries");
            }
            nbp = sol.nb_points();
            p_sol = Some(sol);
        }

        //
        // jgv: build solution on the whole boundary
        if recheck_on_regularity && nbp > 0 && is_regularity(a, domain) {
            let the_tol = 5.0e-4;
            let mut sol_again = FunctionAllRoots::new(func, &echant, eps_x, the_tol, the_tol);

            if !sol_again.is_done() {
                panic!("Standard_Failure: IntStart_SearchOnBoundaries (SolAgain)");
            }

            let nbi_again = sol_again.nb_intervals();

            if nbi_again > 0 {
                let nb_samples = 10;
                let delta = (pfin - pdeb) / nb_samples as f64;
                let global_tol = the_tol * 10.0;
                let mut sol_on_boundary = true;
                for i in 0..=nb_samples {
                    let a_param = pdeb + i as f64 * delta;
                    let Some(a_value) = func.value(a_param) else {
                        continue;
                    };
                    if a_value.abs() > global_tol {
                        sol_on_boundary = false;
                        break;
                    }
                }

                if sol_on_boundary {
                    for i in 1..=nbi_again {
                        let mut newseg = TheSegmentOfTheSearch::new();
                        newseg.set_value(a.clone());
                        // Recuperer point debut et fin, et leur parametre.
                        let (mut pd, mut pf) = sol_again.get_interval(i);
                        pardeb = pd;
                        parfin = pf;

                        if (pardeb - pdeb).abs() <= PCONFUSION {
                            pardeb = pdeb;
                        }
                        if (parfin - pfin).abs() <= PCONFUSION {
                            parfin = pfin;
                        }
                        pd = pardeb;
                        pf = parfin;
                        let _ = (pd, pf);

                        let (id, ifn) = sol_again.get_interval_state(i);
                        ideb = id;
                        ifin = ifn;

                        let ptdeb = func.valpoint(ideb);
                        let ptfin = func.valpoint(ifin);

                        let mut ranged = 0usize;
                        point_process(&ptdeb, pardeb, a, domain, pnt, the_tol, &mut ranged);
                        newseg.set_limit_point(&pnt[ranged - 1], true);
                        let mut rangef = 0usize;
                        point_process(&ptfin, parfin, a, domain, pnt, the_tol, &mut rangef);
                        newseg.set_limit_point(&pnt[rangef - 1], false);
                        seg.push(newseg);
                    }
                    *arcsol = true;
                    return;
                }
            }
        } // if (RecheckOnRegularity && Nbp > 0 && IsRegularity(A, Domain))
        ////////////////////////////////////////////

        //-- detection du cas ou la fonction est quasi tangente et que les
        //-- zeros sont quasi confondus.
        //-- Dans ce cas on prend le point "milieu"
        //-- On suppose que les solutions sont triees.

        if nbp > 0 {
            let mut a_si: Vec<SolInfo> = vec![SolInfo::new(); nbp];

            for (i, si) in a_si.iter_mut().enumerate() {
                let idx = i as i32 + 1;
                // IsIntCSdone is always false for Contap:
                si.init_roots(p_sol.as_ref().unwrap(), idx);
            }

            a_si.sort_by(|x, y| x.value().partial_cmp(&y.value()).unwrap());

            // modified by NIZNHY-PKV Wed Mar 21 18:34:18 2001 (TreatLC)
            let ip = treat_lc(a, domain, &a_quadric, tol_boundary, pnt);
            if ip != 0 {
                //////////////////////////////////////////////////////////
                // modified by NIZNHY-PKV Wed Mar 21 18:34:23 2001 t
                //
                //  Using of old usual way proposed by Laurent
                //
                for i in 1..nbp {
                    let mut parap1 = a_si[i].value(); // aSI(i + 1)
                    let mut para = a_si[i - 1].value(); // aSI(i)

                    let mut param = (para + parap1) * 0.5;
                    let mut yf = 0.0;
                    let mut ym = 0.0;
                    let mut yl = 0.0;
                    if let (Some(v), Some(vf), Some(vl)) =
                        (func.value(param), func.value(para), func.value(parap1))
                    {
                        ym = v;
                        if ym.abs() < maxdist {
                            let sm = ym.signum();
                            let mut a_tang = true;
                            yf = vf;
                            yl = vl;
                            a_tang = a_tang && yf.abs() < maxdist && yl.abs() < maxdist;
                            if a_tang && is_int_cs_done {
                                // (dead for Contap: TypeConS == GeomAbs_Line
                                //  exact-intersection tangency test)
                                let _ = (sm, &mut a_tang);
                            }
                            if a_tang {
                                //  Modified by skv - Tue Aug 31 12:13:51 2004 OCC569
                                let mut a_tol = tol_boundary * 1000.0;
                                if a_tol > 0.001 {
                                    a_tol = 0.001;
                                }

                                // fix floating point exception 569, chl-922-e9
                                parap1 = if parap1.abs() < 1.0e9 {
                                    parap1
                                } else if parap1 >= 0.0 {
                                    1.0e9
                                } else {
                                    -1.0e9
                                };
                                para = if para.abs() < 1.0e9 {
                                    para
                                } else if para >= 0.0 {
                                    1.0e9
                                } else {
                                    -1.0e9
                                };

                                let a_nb_nodes = ((parap1 - para) / a_tol).ceil() as i32;

                                let mut a_val = f64::MAX;
                                let mut a_val_max = 0.0;
                                let a_delta = (parap1 - para) / (a_nb_nodes as f64 + 1.0);
                                let mut a_cur_par;
                                let mut a_cur_val;

                                for ii in 0..=(a_nb_nodes + 1) {
                                    a_cur_par = if ii < a_nb_nodes + 1 {
                                        para + ii as f64 * a_delta
                                    } else {
                                        parap1
                                    };

                                    if let Some(v) = func.value(a_cur_par) {
                                        a_cur_val = v;
                                        let an_abs_val = a_cur_val.abs();
                                        if an_abs_val < a_val {
                                            a_val = an_abs_val;
                                            param = a_cur_par;
                                        }
                                        if an_abs_val > a_val_max {
                                            a_val_max = an_abs_val;
                                        }
                                    }
                                }
                                // At last, interval got by exact intersection can be considered
                                // as tangent if minimal distance is inside interval and
                                // minimal and maximal values are almost the same
                                if is_int_cs_done && a_nb_nodes > 1 {
                                    a_tang = (param - para).abs() > eps_x
                                        && (parap1 - param).abs() > eps_x
                                        && 0.01 * a_val_max <= a_val;
                                }
                                if a_tang {
                                    *a_si[i - 1].change_value() = pdeb - 1.0;
                                    *a_si[i].change_value() = param;
                                }
                            }
                        }
                    }
                }

                for i in 1..=nbp {
                    let mut para = a_si[i - 1].value();
                    if (para - pdeb) < eps_x || (pfin - para) < eps_x {
                        continue;
                    }

                    let Some(mut dist) = func.value(para) else {
                        continue;
                    };

                    dist = dist.abs();

                    let mut an_indx = -1i32;
                    let a_param = a_si[i - 1].value();
                    if dist < maxdist
                        && !is_int_cs_done
                        && ((a_param - pdeb).abs() <= PCONFUSION
                            || (a_param - pfin).abs() <= PCONFUSION)
                    {
                        an_indx = p_sol.as_ref().unwrap().get_point_state(a_si[i - 1].index() as usize);
                    }

                    let mut a_pnt = if an_indx < 0 {
                        func.last_computed_point()
                    } else {
                        func.valpoint(an_indx)
                    };

                    if dist > 0.1 * CONFUSION {
                        // Precise found points (#27252): refine the parameter
                        // with math_BrentMinimum between the neighbours,
                        // keeping the sort order.
                        let a_f_par = if i == 1 {
                            pdeb
                        } else {
                            (para + a_si[i - 2].value()) / 2.0
                        };
                        let a_l_par = if i == nbp {
                            pfin
                        } else {
                            (para + a_si[i].value()) / 2.0
                        };

                        let mut a_new_func = MinFunction::new(func);
                        let mut a_min = BrentMinimum::new(CONFUSION, 100, 2.0 * CONFUSION);
                        // OCCT: math_BrentMinimum aMin(Precision::Confusion())
                        //   -> Zeps defaults to the machine epsilon; the
                        //   kernel BrentMinimum ctor takes (tol, nbiter,
                        //   zeps) and the A/B/C Perform below.
                        a_min.perform(&mut a_new_func, a_f_par, para, a_l_par);
                        if a_min.is_done() {
                            para = a_min.location();
                            let a_p2d = a.value(para);
                            a_pnt = func.surface().value(a_p2d.x, a_p2d.y);
                        }
                    }

                    let mut range = 0usize;
                    point_process(&a_pnt, para, a, domain, pnt, tol_boundary, &mut range);
                } // end of for (i = 1; i <= Nbp; i++)
            } // end of if(ip)
        } // end of if(Nbp)

        // Pour chaque intervalle trouve faire
        //   Traiter les extremites comme des points
        //   Ajouter intervalle dans la liste des segments

        if !is_int_cs_done {
            nbi = p_sol.as_ref().unwrap().nb_intervals();
        }

        if !recheck_on_regularity && nbp > 0 {
            //--cout<<" Debug : IntStart_SearchOnBoundaries_1.gxx :Nbp>0  0 <- Nbi "<<Nbi<<endl;
            nbi = 0;
        }

        for i in 1..=nbi {
            let mut newseg = TheSegmentOfTheSearch::new();
            newseg.set_value(a.clone());
            // Recuperer point debut et fin, et leur parametre.
            // (IsIntCSdone is false for Contap — the IntCS arm stays in the
            //  dead branch.)
            let (pd, pf) = p_sol.as_ref().unwrap().get_interval(i);
            pardeb = pd;
            parfin = pf;
            let (id, ifn) = p_sol.as_ref().unwrap().get_interval_state(i);
            ideb = id;
            ifin = ifn;

            let ptdeb = func.valpoint(ideb);
            let ptfin = func.valpoint(ifin);

            let mut ranged = 0usize;
            point_process(&ptdeb, pardeb, a, domain, pnt, tol_boundary, &mut ranged);
            newseg.set_limit_point(&pnt[ranged - 1], true);
            let mut rangef = 0usize;
            point_process(&ptfin, parfin, a, domain, pnt, tol_boundary, &mut rangef);
            newseg.set_limit_point(&pnt[rangef - 1], false);
            seg.push(newseg);
        }

        if nbi == 1 {
            if (pardeb - pdeb).abs() < PCONFUSION && (parfin - pfin).abs() < PCONFUSION {
                *arcsol = true;
            }
        }
    }
}

/// OCCT ComputeBoundsfromInfinite (gxx L804-876).
fn compute_bounds_from_infinite(func: &mut ArcFunction, p_deb: &mut f64, p_fin: &mut f64) {
    // - PROVISIONAL - TEMPORARY - NOT GOOD - NYI - TO DO (OCCT comment)
    let nb_echant = 100;
    let _ = nb_echant;

    let u0 = 0.0f64;
    let du = 0.001f64;
    let Some(dist0) = func.value(u0) else {
        *p_deb = 1.0e10;
        *p_fin = -1.0e10;
        return;
    };
    let Some(dist1) = func.value(u0 + du) else {
        *p_deb = 1.0e10;
        *p_fin = -1.0e10;
        return;
    };
    let mut d_dist = dist1 - dist0;
    if d_dist != 0.0 {
        let mut u0 = u0 - du * dist0 / d_dist;
        *p_deb = u0;
        *p_fin = u0;
        let mut umin = u0 - 1.0e5;
        let d0 = func.value(umin).unwrap_or(0.0);
        let d1 = func.value(umin + du).unwrap_or(0.0);
        d_dist = d1 - d0;
        if d_dist != 0.0 {
            umin -= du * d0 / d_dist;
        } else {
            umin -= 10.0;
        }
        let mut umax = u0 + 1.0e8;
        let d0 = func.value(umax).unwrap_or(0.0);
        let d1 = func.value(umax + du).unwrap_or(0.0);
        d_dist = d1 - d0;
        if d_dist != 0.0 {
            umax -= du * d0 / d_dist;
        } else {
            umax += 10.0;
        }
        if umin > u0 {
            umin = u0 - 10.0;
        }
        if umax < u0 {
            umax = u0 + 10.0;
        }

        *p_fin = umax + 10.0 * (umax - umin);
        *p_deb = umin - 10.0 * (umax - umin);
    } else {
        //-- Possibilite de Arc totalement inclu ds Quad
        *p_deb = 1.0e10;
        *p_fin = -1.0e10;
    }
}

/// OCCT PointProcess (gxx L880-1001).
fn point_process(
    pt: &rcad_kernel::geom::Point3,
    para: f64,
    a: &Arc,
    domain: &mut dyn ContapDomain,
    pnt: &mut Vec<ThePathPointOfTheSearch>,
    tol: f64,
    range: &mut usize,
) {
    // Check to see if a solution point is coincident with a vertex.
    let nbsol = pnt.len();

    domain.initialize_arc(a);
    domain.init_vertex_iterator();
    let mut found = false;
    let mut goon = domain.more_vertex();
    while goon {
        let vtx = domain.vertex();
        let dist = (para - crate::hlr::contap::h_cont_tool::parameter(vtx.as_ref(), a.as_ref())).abs();
        let mut toler = crate::hlr::contap::h_cont_tool::tolerance(vtx.as_ref(), a.as_ref());

        if dist <= toler {
            // Locate the vertex in the list of solutions
            let mut k = 1usize;
            found = k > nbsol;
            while !found {
                let ptsol = &pnt[k - 1];
                if !ptsol.is_new() {
                    // jag 940608  if (ptsol.Vertex() == vtx && ptsol.Arc() == A) {
                    let same_arc = std::sync::Arc::ptr_eq(&ptsol.arc(), a);
                    if domain.identical(ptsol.vertex(), &vtx)
                        && same_arc
                        && (ptsol.parameter() - para).abs() <= toler
                    {
                        found = true;
                    } else {
                        k += 1;
                        found = k > nbsol;
                    }
                } else {
                    k += 1;
                    found = k > nbsol;
                }
            }
            if k <= nbsol {
                // We find the vertex
                *range = k;
            } else {
                // Otherwise
                let mut ptsol = ThePathPointOfTheSearch::new();
                ptsol.set_value_vertex(*pt, tol, vtx, a.clone(), para);
                pnt.push(ptsol);
                *range = pnt.len();
            }
            found = true;
            goon = false;
        } else {
            domain.next_vertex();
            goon = domain.more_vertex();
        }
    }

    if !found {
        // No one is falling on a vertex
        // jgv: do not add segment's extremities if they already exist
        let mut found_internal = false;
        for (k, ptsol) in pnt.iter().enumerate() {
            let same_arc = std::sync::Arc::ptr_eq(&ptsol.arc(), a);
            if !same_arc || !ptsol.is_new() {
                // vertex
                continue;
            }
            if (ptsol.parameter() - para).abs() <= PCONFUSION {
                found_internal = true;
                *range = k + 1;
            }
        }
        /////////////////////////////////////////////////////////////

        if !found_internal {
            let mut new_tol = tol;
            new_tol *= 1000.0;
            // if(TOL>0.001) TOL=0.001;
            if new_tol > 0.005 {
                new_tol = 0.005; // #24643
            }

            let mut ptsol = ThePathPointOfTheSearch::new();
            ptsol.set_value(*pt, new_tol, a.clone(), para);
            pnt.push(ptsol);
            *range = pnt.len();
        }
    }
}

/// OCCT IsRegularity (gxx L1005-1016) — the base TopolTool has no edge.
fn is_regularity(_a: &Arc, a_domain: &mut dyn ContapDomain) -> bool {
    let an_e_address = a_domain.edge();
    if an_e_address.is_none() {
        return false;
    }
    // BRep branch: BRep_Tool::HasContinuity(*anE) — Stage 3b/3d.
    false
}

/// OCCT TreatLC (gxx L1020-1123) — the base TopolTool has no edge, so the
//  OCCT body always exits at the first check with code 1.
fn treat_lc(
    _a: &Arc,
    a_domain: &mut dyn ContapDomain,
    _a_quadric: &crate::geomalgo::int_surf::quadric::Quadric,
    _tol_boundary: f64,
    _pnt: &mut Vec<ThePathPointOfTheSearch>,
) -> i32 {
    let an_exit_code = 1;

    let an_e_address = a_domain.edge();
    if an_e_address.is_none() {
        return an_exit_code;
    }

    // BRep branch (Degenerated-TLine-Cylinder tangent-line treatment,
    // BRepAdaptor_Curve + Extrema_ExtCC) — Stage 3b/3d.
    an_exit_code
}

/// OCCT IntStart_SearchOnBoundaries.
#[derive(Clone)]
pub struct TheSearch {
    done: bool,
    all: bool,
    sseg: Vec<TheSegmentOfTheSearch>,
    spnt: Vec<ThePathPointOfTheSearch>,
}

impl TheSearch {
    /// OCCT IntStart_SearchOnBoundaries() (gxx L1127-1131).
    pub fn new() -> Self {
        TheSearch {
            done: false,
            all: false,
            sseg: Vec::new(),
            spnt: Vec::new(),
        }
    }

    /// OCCT Perform(Func, Domain, TolBoundary, TolTangency,
    /// RecheckOnRegularity) (gxx L1135-1232).
    pub fn perform(
        &mut self,
        func: &mut ArcFunction,
        domain: &mut dyn ContapDomain,
        tol_boundary: f64,
        tol_tangency: f64,
        recheck_on_regularity: bool,
    ) {
        self.done = false;
        self.spnt.clear();
        self.sseg.clear();

        let mut arcsol = false;
        let mut p_deb = 0.0;
        let mut p_fin = 0.0;

        domain.init();

        if domain.more() {
            self.all = true;
        } else {
            self.all = false;
        }

        while domain.more() {
            let a = domain.arc();
            if !crate::hlr::contap::h_cont_tool::has_been_seen(a.as_ref()) {
                func.set_arc(a.clone());
                find_vertex(&a, domain, func, &mut self.spnt, tol_boundary);
                let (mut pd, mut pf) = crate::hlr::contap::h_cont_tool::bounds(a.as_ref());
                p_deb = pd;
                p_fin = pf;
                if is_negative_infinite_value(p_deb) || is_positive_infinite_value(p_fin) {
                    compute_bounds_from_infinite(func, &mut pd, &mut pf);
                    p_deb = pd;
                    p_fin = pf;
                }
                bounded_arc(
                    &a,
                    domain,
                    p_deb,
                    p_fin,
                    func,
                    &mut self.spnt,
                    &mut self.sseg,
                    tol_boundary,
                    tol_tangency,
                    &mut arcsol,
                    recheck_on_regularity,
                );
                self.all = self.all && arcsol;
            } else {
                // as it seems we'll never be here, because
                // TheSOBTool::HasBeenSeen(A) always returns FALSE
                let nbfound = self.spnt.len();

                // On recupere les points connus
                let nbknown = crate::hlr::contap::h_cont_tool::nb_points(a.as_ref());
                for _i in 1..=nbknown {
                    // TheSOBTool::Value raises Standard_OutOfRange for the
                    // Contap instantiation (NbPoints() == 0).
                    let (pt, tol, prm) = crate::hlr::contap::h_cont_tool::value(a.as_ref(), 1);
                    if crate::hlr::contap::h_cont_tool::is_vertex(a.as_ref(), 1) {
                        let vtx = std::sync::Arc::new(crate::hlr::contap::h_cont_tool::vertex(
                            a.as_ref(),
                            1,
                        ));
                        self.spnt.push(ThePathPointOfTheSearch::with_vertex(
                            pt, tol, vtx, a.clone(), prm,
                        ));
                    } else {
                        self.spnt
                            .push(ThePathPointOfTheSearch::without_vertex(pt, tol, a.clone(), prm));
                    }
                }
                // On recupere les arcs solutions
                let nbknown = crate::hlr::contap::h_cont_tool::nb_segments(a.as_ref());
                for _i in 1..=nbknown {
                    let mut newseg = TheSegmentOfTheSearch::new();
                    newseg.set_value(a.clone());
                    let (hasf, index) =
                        crate::hlr::contap::h_cont_tool::has_first_point(a.as_ref(), 1);
                    if hasf {
                        newseg.set_limit_point(&self.spnt[nbfound + index as usize - 1], true);
                    }
                    let (hasl, index) =
                        crate::hlr::contap::h_cont_tool::has_last_point(a.as_ref(), 1);
                    if hasl {
                        newseg.set_limit_point(&self.spnt[nbfound + index as usize - 1], false);
                    }
                    self.sseg.push(newseg);
                }
                self.all = self.all
                    && crate::hlr::contap::h_cont_tool::is_all_solution(a.as_ref());
            }
            domain.next();
        }
        self.done = true;
    }

    /// OCCT IsDone (IntStart_SearchOnBoundaries.lxx).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT AllArcSolution.
    pub fn all_arc_solution(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.all
    }

    /// OCCT NbPoints.
    pub fn nb_points(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.spnt.len()
    }

    /// OCCT Point(Index) — 1-based.
    pub fn point(&self, index: usize) -> &ThePathPointOfTheSearch {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        &self.spnt[index - 1]
    }

    /// OCCT NbSegments.
    pub fn nb_segments(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.sseg.len()
    }

    /// OCCT Segment(Index) — 1-based.
    pub fn segment(&self, index: usize) -> &TheSegmentOfTheSearch {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        &self.sseg[index - 1]
    }
}

impl Default for TheSearch {
    fn default() -> Self {
        Self::new()
    }
}

#[allow(unused)]
fn _parity(confusion: f64) -> f64 {
    CONFUSION + confusion
}
