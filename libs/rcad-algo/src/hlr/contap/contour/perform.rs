// OCCT Contap_Contour::Perform(const handle<Adaptor3d_TopolTool>& Domain)
// (Contap_Contour.cxx L1541-1971) — the walking path.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::Point3;
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::precision::{is_negative_infinite_value, is_positive_infinite_value, CONFUSION};
use rcad_kernel::topo::topods::State;

use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::geomalgo::int_patch::imp_prm::i_walking::IWalking;
use crate::geomalgo::int_patch::imp_prm::path_point::{InteriorPoint, PathPoint};
use crate::geomalgo::int_patch::transitions::{make_transition, Transition};
use crate::hlr::contap::contour::functions::{
    compute_internal_points, compute_transition_on_line, compute_tangency, keep_inside_points,
    line_constructor, process_segments,
};
use crate::hlr::contap::domain::ContapDomain;
use crate::hlr::contap::h_cont_tool as hcont;
use crate::hlr::contap::h_curve2d_tool as hcurve2d;
use crate::hlr::contap::i_type::IType;
use crate::hlr::contap::line::Line;
use crate::hlr::contap::point::Point as ContapPoint;
use crate::hlr::contap::surf_props;
use crate::hlr::contap::t_function::TFunction;

use super::Contour;

/// OCCT Perform(Domain) (cxx L1541-1971).
pub(crate) fn perform_domain(c: &mut Contour, domain: &mut dyn ContapDomain) {
    c.set_done(false);
    c.slin_mut().clear();

    let mut k = 0usize;
    let mut u = 0.0;
    let mut v = 0.0;
    let mut pt2d = DVec2::ZERO;
    let mut d2d = DVec2::ZERO;
    let mut ptonsurf = Point3::ZERO;
    let mut d1u = DVec3::ZERO;
    let mut d1v = DVec3::ZERO;
    let mut normale = DVec3::ZERO;
    let mut tgtrst = DVec3::ZERO;
    let mut tgline = DVec3::ZERO;
    let mut currentparam = 0.0;
    let mut tline = Transition::new();
    let mut tarc = Transition::new();

    //  double TolArc = 1.e-5;
    let mut tol_arc = CONFUSION; // cxx L1563

    let surf = c.sfunc_ref().surface().clone();

    let eps_u = surf.u_resolution(CONFUSION);
    let eps_v = surf.v_resolution(CONFUSION);
    let preci = eps_u.min(eps_v);
    //  double Fleche = 5.e-1;
    //  double Pas    = 5.e-2;
    let mut fleche = 0.01;
    let pas = 0.005;
    //-- le 23 janvier 98 0.05 -> 0.01

    //-- ******************************************************************************** Janvier 98
    let mut b1 = BndBox::new();
    let box1_ok;
    {
        let uinf = surf.first_u_parameter();
        let vinf = surf.first_v_parameter();
        let usup = surf.last_u_parameter();
        let vsup = surf.last_v_parameter();

        let uinfinfinite = is_negative_infinite_value(uinf);
        let usupinfinite = is_positive_infinite_value(usup);
        let vinfinfinite = is_negative_infinite_value(vinf);
        let vsupinfinite = is_positive_infinite_value(vsup);

        if uinfinfinite || usupinfinite || vinfinfinite || vsupinfinite {
            box1_ok = false;
        } else {
            bnd_lib_add_surface(&*surf, 1.0e-8, &mut b1);
            box1_ok = true;
        }
    }
    let mut dx;
    let dy;
    let dz;
    if box1_ok {
        let (x0, y0, z0, x1, y1, z1) = b1.get().unwrap_or((0.0, 0.0, 0.0, 1.0, 1.0, 1.0));
        dx = x1 - x0;
        dy = y1 - y0;
        dz = z1 - z0;
    } else {
        dx = 1.0;
        dy = 1.0;
        dz = 1.0;
    }
    if dx < dy {
        dx = dy;
    }
    if dx < dz {
        dx = dz;
    }
    if dx > 10000.0 {
        dx = 10000.0;
    }
    fleche *= dx;
    tol_arc *= dx;
    //-- ********************************************************************************

    // jag 940616  SFunc.Set(1.e-8); // tolerance sur la fonction
    c.my_sfunc().set_tolerance(CONFUSION); // tolerance sur la fonction

    let recheck_on_regularity = true;
    {
        let (f, af, srst, _sins, _slin) = c.parts_mut();
        srst.perform(af, domain, tol_arc, tol_arc, recheck_on_regularity);
    }

    if !c.solrst().is_done() {
        return;
    }

    let nb_point_rst = c.solrst().nb_points();
    let mut seqpdep: Vec<PathPoint> = Vec::new();
    // NCollection_Array1<int> Destination(1, NbPointRst + 1); Init(0) — the
    // array lower bound is 1, upper bound is NbPointRst + 1 (ComputeTangency
    // indexes it 1-based). The slot 0 stays unused.
    let mut destination = vec![0i32; nb_point_rst + 2];
    if nb_point_rst != 0 {
        let (f, _af, srst, _sins, _slin) = c.parts_mut();
        compute_tangency(srst, domain, f, &mut seqpdep, &mut destination[1..=nb_point_rst + 1]);
    }

    // jag 940616  solins.Perform(SFunc,Surf,Domain,1.e-6);
    {
        let (f, _af, _srst, sins, _slin) = c.parts_mut();
        sins.perform(f, &surf, domain, CONFUSION);
    }

    let nb_point_ins = c.solins().nb_points();
    let mut seqpins: Vec<InteriorPoint> = Vec::new();

    if nb_point_ins != 0 {
        let mut b_keep_all_points = false;
        // IFV begin
        if c.solrst().nb_segments() == 0 {
            if c.sfunc_ref().function_type() == TFunction::ContourStd {
                let surf_to_check = c.sfunc_ref().surface().clone();
                if surf_to_check.get_type() == GeomAbsSurfaceType::Torus {
                    let a_tor = surf_to_check.torus();
                    let a_tor_dir = a_tor.axis;
                    let a_proj_dir = c.sfunc_ref().direction();

                    if a_tor_dir.dot(a_proj_dir) < CONFUSION {
                        b_keep_all_points = true;
                    }
                }
            }
        }

        if b_keep_all_points {
            let nbp = c.solins().nb_points();
            for indp in 1..=nbp {
                let pti = c.solins().value(indp);
                seqpins.push(pti.clone());
            }
        }
        // IFV - end
        else {
            let (f, _af, srst, sins, _slin) = c.parts_mut();
            keep_inside_points(sins, srst, f, &mut seqpins);
        }
    }

    if !seqpdep.is_empty() || !seqpins.is_empty() {
        let the_to_fill_holes = true;
        let mut iwalk = IWalking::new(preci, fleche, pas, the_to_fill_holes);
        {
            // OCCT: iwalk.Perform(seqpdep, seqpins, mySFunc, Surf);
            // The rcad IntWalk engine is the Surface3 instantiation; the
            // Contap function carries the adaptor (its Func.Set(Caro) is
            // re-run against the stored adaptor — idempotent, cxx L152).
            let caro = match surf.kernel_surface() {
                Some(s3) => s3,
                None => panic!("Contap_TheIWalking: the adapter does not expose a kernel surface"),
            };
            let domain4 = surf.uv_window();
            let func = c.my_sfunc();
            iwalk.perform(&seqpdep, &seqpins, func, caro, domain4, false);
        }
        if !iwalk.is_done() {
            return;
        }

        let nblines = iwalk.nb_lines();
        // TEMP-DEBUG
        if std::env::var("RCAD_IWALK_DEBUG").is_ok() {
            eprintln!(
                "[CWDBG] nb_point_rst={} nb_point_ins={} seqpdep={} seqpins={} iwalk_done={} nblines={}",
                nb_point_rst,
                nb_point_ins,
                seqpdep.len(),
                seqpins.len(),
                iwalk.is_done(),
                nblines
            );
        }
        for j in 1..=nblines {
            let iwline = iwalk.value(j);
            let nbpts = iwline.nb_points();
            let mut theline = Line::new();
            theline.set_line_on_2s(iwline.line());

            // jag 941018 On calcule une seule fois la transition
            let (tgline_v, kk) = iwline.tangent_vector();
            k = kk as usize;
            // OCCT cxx L1717: iwline->Line()->Value(k) — IntSurf_LineOn2S::Value
            // is 1-based; the rcad LineOn2S::value is 0-based.
            let (u_v, v_v) = iwline.line().value(k - 1).parameters_on_surface(false);
            u = u_v;
            v = v_v;
            tgline = tgline_v;
            let type_trans_on_s = compute_transition_on_line(c.my_sfunc(), u, v, tgline);
            theline.set_transition_on_s(type_trans_on_s);

            //---------------------------------------------------------------------
            //-- On ajoute a la liste des vertex les 1er et dernier points de la  -
            //-- ligne de cheminement si ceux-ci ne sont pas presents             -
            //---------------------------------------------------------------------

            if iwline.has_first_point() {
                let indfirst = iwline.first_point_index() as usize;
                let ppoint = &seqpdep[indfirst - 1];
                let mut themult = ppoint.multiplicity();
                let mut i = nb_point_rst;
                while i >= 1 {
                    if destination[i] == indfirst as i32 {
                        let (mut uu, mut vv) = (0.0, 0.0);
                        ppoint.parameters(themult, &mut uu, &mut vv);
                        u = uu;
                        v = vv;
                        let mut ptdeb = ContapPoint::with_uv(ppoint.value(), u, v);
                        ptdeb.set_parameter(1.0);

                        let pstart = c.solrst().point(i);
                        let currentarc = pstart.arc();
                        currentparam = pstart.parameter();
                        if !iwline.is_tangent_at_begining() {
                            let (p2d_v, d2d_v) = hcurve2d::d1(currentarc.as_ref(), currentparam);
                            pt2d = p2d_v;
                            d2d = d2d_v;
                            surf_props::deriv_and_norm(
                                &*surf,
                                pt2d.x,
                                pt2d.y,
                                &mut ptonsurf,
                                &mut d1u,
                                &mut d1v,
                                &mut normale,
                            );
                            tgtrst = d2d.x * d1u;
                            tgtrst += d2d.y * d1v;

                            make_transition(
                                ppoint.direction_3d(),
                                tgtrst,
                                normale,
                                &mut tline,
                                &mut tarc,
                            );
                        } else {
                            // a voir. En effet, on a cheminer. Si on est sur un
                            // point debut, on sait qu'on rentre dans la matiere
                            tline.set_value_undecided();
                            tarc.set_value_undecided();
                        }

                        ptdeb.set_arc(currentarc.clone(), currentparam, tline, tarc);

                        if !c.solrst().point(i).is_new() {
                            ptdeb.set_vertex(c.solrst().point(i).vertex().clone());
                        }
                        theline.add(ptdeb);
                        themult -= 1;
                    }
                    i -= 1;
                }
            } else {
                let (u_v, v_v) = iwline.value(1).parameters_on_surface(false);
                u = u_v;
                v = v_v;
                let mut ptdeb = ContapPoint::with_uv(theline.point(1).value(), u, v);
                ptdeb.set_parameter(1.0);
                theline.add(ptdeb);
            }

            if iwline.has_last_point() {
                let indlast = iwline.last_point_index() as usize;
                let ppoint = &seqpdep[indlast - 1];
                let mut themult = ppoint.multiplicity();
                let mut i = nb_point_rst;
                while i >= 1 {
                    if destination[i] == indlast as i32 {
                        let (mut uu, mut vv) = (0.0, 0.0);
                        ppoint.parameters(themult, &mut uu, &mut vv);
                        u = uu;
                        v = vv;
                        let mut ptfin = ContapPoint::with_uv(ppoint.value(), u, v);
                        ptfin.set_parameter(nbpts as f64);
                        let pstart = c.solrst().point(i);
                        let currentarc = pstart.arc();
                        currentparam = pstart.parameter();

                        if !iwline.is_tangent_at_end() {
                            let (p2d_v, d2d_v) = hcurve2d::d1(currentarc.as_ref(), currentparam);
                            pt2d = p2d_v;
                            d2d = d2d_v;

                            surf_props::deriv_and_norm(
                                &*surf,
                                pt2d.x,
                                pt2d.y,
                                &mut ptonsurf,
                                &mut d1u,
                                &mut d1v,
                                &mut normale,
                            );
                            tgtrst = d2d.x * d1u;
                            tgtrst += d2d.y * d1v;
                            make_transition(
                                -ppoint.direction_3d(),
                                tgtrst,
                                normale,
                                &mut tline,
                                &mut tarc,
                            );
                        } else {
                            tline.set_value_undecided();
                            tarc.set_value_undecided();
                        }

                        ptfin.set_arc(currentarc.clone(), currentparam, tline, tarc);

                        if !c.solrst().point(i).is_new() {
                            ptfin.set_vertex(c.solrst().point(i).vertex().clone());
                        }
                        theline.add(ptfin);
                        themult -= 1;
                    }
                    i -= 1;
                }
            } else {
                let (u_v, v_v) = iwline.value(nbpts).parameters_on_surface(false);
                u = u_v;
                v = v_v;
                let mut ptfin = ContapPoint::with_uv(theline.point(nbpts).value(), u, v);
                ptfin.set_parameter(nbpts as f64);
                theline.add(ptfin);
            }

            compute_internal_points(&mut theline, c.my_sfunc(), eps_u, eps_v);
            line_constructor(c.slin_mut(), domain, &mut theline, &*surf); //-- lbr
            theline.reset_seq_of_vertex();
            if std::env::var("RCAD_IWALK_DEBUG").is_ok() {
                eprintln!(
                    "[CWDBG] after line_constructor j={j}: contour nb_lines={}",
                    c.slin().len()
                );
            }
        }

        // cxx L1839-1870 — the crossing-vertex multiple marking.
        let nblines = c.slin().len();
        for j in 1..nblines {
            let nbvt1 = c.slin()[j - 1].nb_vertex();
            for ivt1 in 1..=nbvt1 {
                let (on_arc1, pttg1) = {
                    let theli = &c.slin()[j - 1];
                    (theli.vertex(ivt1).is_on_arc(), theli.vertex(ivt1).value())
                };
                if !on_arc1 {
                    for k in (j + 1)..=nblines {
                        let nbvt2 = c.slin()[k - 1].nb_vertex();
                        for ivt2 in 1..=nbvt2 {
                            let (on_arc2, pttg2) = {
                                let theli2 = &c.slin()[k - 1];
                                (theli2.vertex(ivt2).is_on_arc(), theli2.vertex(ivt2).value())
                            };
                            if !on_arc2 && pttg1.distance(pttg2) <= tol_arc {
                                c.slin_mut()[j - 1].vertex_mut(ivt1).set_multiple();
                                c.slin_mut()[k - 1].vertex_mut(ivt2).set_multiple();
                            }
                        }
                    }
                }
            }
        }
    }

    // TEMP-DEBUG
    if std::env::var("RCAD_IWALK_DEBUG").is_ok() {
        for (li, l) in c.slin().iter().enumerate() {
            let vinfo: Vec<String> = (1..=l.nb_vertex())
                .map(|iv| {
                    let vt = l.vertex(iv);
                    format!(
                        "[onarc={} mult={} par={:.3}]",
                        vt.is_on_arc(),
                        vt.is_multiple(),
                        vt.parameter_on_line()
                    )
                })
                .collect();
            eprintln!(
                "[SLIN] {} typ={:?} nbpnts={} trans={:?} vertices={:?}",
                li + 1,
                l.type_contour(),
                if l.type_contour() == crate::hlr::contap::i_type::IType::Walking {
                    l.nb_pnts()
                } else {
                    0
                },
                l.transition_on_s(),
                vinfo
            );
        }
        eprintln!("[SLIN] total={}", c.slin().len());
    }

    // jag 940620 On ajoute le traitement des restrictions solutions.

    if c.solrst().nb_segments() != 0 {
        let (f, _af, srst, _sins, slin) = c.parts_mut();
        process_segments(srst, slin, tol_arc, f, domain);
    }

    // Ajout crad pour depanner CMA en attendant mieux
    if c.solrst().nb_segments() != 0 {
        let nblines = c.slin().len();
        for j in 1..=nblines {
            if c.slin()[j - 1].type_contour() == IType::Walking {
                let nbvt1 = c.slin()[j - 1].nb_vertex();
                for ivt1 in 1..=nbvt1 {
                    let (on_arc, is_multiple, up, vp, ptvt_value, param_on_line) = {
                        let ptvt = c.slin()[j - 1].vertex(ivt1);
                        let (up, vp) = ptvt.parameters();
                        (
                            ptvt.is_on_arc(),
                            ptvt.is_multiple(),
                            up,
                            vp,
                            ptvt.value(),
                            ptvt.parameter_on_line(),
                        )
                    };
                    if !on_arc && !is_multiple {
                        let toproj = DVec2::new(up, vp);
                        for k in 1..=nblines {
                            if c.slin()[k - 1].type_contour() == IType::Restriction {
                                let thearc = c.slin()[k - 1].arc();
                                if let Some((paramproj, ptproj)) =
                                    hcont::project(thearc.as_ref(), toproj)
                                {
                                    let dist = ptproj.distance(toproj);
                                    if dist <= preci {
                                        // Calcul de la transition

                                        let (pt2d_v, d2d_v) =
                                            hcurve2d::d1(thearc.as_ref(), paramproj);
                                        pt2d = pt2d_v;
                                        d2d = d2d_v;

                                        surf_props::deriv_and_norm(
                                            &*surf,
                                            ptproj.x,
                                            ptproj.y,
                                            &mut ptonsurf,
                                            &mut d1u,
                                            &mut d1v,
                                            &mut normale,
                                        );

                                        tgtrst = d2d.x * d1u;
                                        tgtrst += d2d.y * d1v;
                                        let paraml = param_on_line as usize;

                                        if paraml == c.slin()[j - 1].nb_pnts() {
                                            tgline = ptvt_value
                                                - c.slin()[j - 1].point(paraml - 1).value();
                                        } else {
                                            tgline = c.slin()[j - 1].point(paraml + 1).value()
                                                - ptvt_value;
                                        }
                                        make_transition(
                                            tgline,
                                            tgtrst,
                                            normale,
                                            &mut tline,
                                            &mut tarc,
                                        );
                                        {
                                            let ptvt =
                                                c.slin_mut()[j - 1].vertex_mut(ivt1);
                                            ptvt.set_arc(thearc.clone(), paramproj, tline, tarc);
                                            ptvt.set_multiple();
                                        }
                                        let mut ptdeb =
                                            ContapPoint::with_uv(ptonsurf, ptproj.x, ptproj.y);
                                        ptdeb.set_parameter(paramproj);
                                        ptdeb.set_multiple();
                                        c.slin_mut()[k - 1].add(ptdeb);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    c.set_done(true);
}

/// OCCT BndLib_AddSurface::Add(S, Tol, B) — the box feeds only the Fleche/
/// TolArc diagonal scale (cxx L1598).  The kernel BndLib port is a separate
/// runway item; this conservative evaluation samples the surface over its
/// finite UV window (the caller guards the infinite domains).
pub(crate) fn bnd_lib_add_surface(
    s: &dyn crate::hlr::contap::surface_adaptor::SurfaceAdapter,
    tol: f64,
    b: &mut BndBox,
) {
    let nu = 16usize;
    let nv = 16usize;
    let uinf = s.first_u_parameter();
    let usup = s.last_u_parameter();
    let vinf = s.first_v_parameter();
    let vsup = s.last_v_parameter();
    for iu in 0..=nu {
        let u = uinf + (usup - uinf) * (iu as f64) / (nu as f64);
        for iv in 0..=nv {
            let vv = vinf + (vsup - vinf) * (iv as f64) / (nv as f64);
            b.add_point(s.value(u, vv));
        }
    }
    b.enlarge(tol);
}
