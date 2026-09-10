// OCCT Contap_Contour::PerformAna(const handle<Adaptor3d_TopolTool>&
// Domain) (Contap_Contour.cxx L2156-2389) — the analytic path.

use rcad_kernel::precision::{ANGULAR, CONFUSION};

use crate::geomalgo::geom2d_int::Curve2dType;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::hlr::contap::cont_ana::ContAna;
use crate::hlr::contap::contour::functions::{
    compute_transition_on_gp_circle, compute_transition_on_gp_line, line_constructor,
    process_segments, put_points_on_line,
};
use crate::hlr::contap::domain::ContapDomain;
use crate::hlr::contap::i_type::IType;
use crate::hlr::contap::line::Line;
use crate::hlr::contap::t_function::TFunction;

use super::Contour;

/// OCCT PerformAna(Domain) (cxx L2156-2389).
pub(crate) fn perform_ana(c: &mut Contour, domain: &mut dyn ContapDomain) {
    c.set_done(false);
    c.slin_mut().clear();

    let tol_arc = 1.0e-5;

    let surf = c.sfunc_ref().surface().clone();
    let type_func = c.sfunc_ref().function_type();
    let mut perform_sol_rst = true;

    let typ_s = surf.get_type();

    // Contap_ContAna contana; — filled by the dispatch below.
    let mut contana = ContAna::new();

    match typ_s {
        GeomAbsSurfaceType::Plane => {
            let pl = surf.plane();
            match type_func {
                TFunction::ContourStd => {
                    let dirpln = pl.normal;
                    if c.sfunc_ref().direction().dot(dirpln).abs() > ANGULAR {
                        // Precision::Angular() — Aucun point du plan n'est
                        // solution, en particulier aucun point sur restriction.
                        perform_sol_rst = false;
                    }
                }
                TFunction::ContourPrs => {
                    let eye = c.sfunc_ref().eye();
                    // pl.Distance(Eye)
                    let d = (eye - pl.origin).dot(pl.normal).abs();
                    if d > CONFUSION {
                        // Aucun point du plan n'est solution, en particulier
                        // aucun point sur restriction.
                        perform_sol_rst = false;
                    }
                }
                TFunction::DraftStd => {
                    let dirpln = pl.normal;
                    let sina = c.sfunc_ref().angle().sin();
                    // voir SurfFunction
                    if c.sfunc_ref().direction().dot(dirpln) + sina > ANGULAR {
                        perform_sol_rst = false;
                    }
                }
                _ => {}
            }
        }

        GeomAbsSurfaceType::Sphere => {
            let sp = surf.sphere();
            match type_func {
                TFunction::ContourStd => {
                    contana.perform_sphere_dir(&sp, c.sfunc_ref().direction());
                }
                TFunction::ContourPrs => {
                    contana.perform_sphere_eye(&sp, c.sfunc_ref().eye());
                }
                TFunction::DraftStd => {
                    let d = c.sfunc_ref().direction();
                    let a = c.sfunc_ref().angle();
                    contana.perform_sphere_dir_angle(&sp, d, a);
                }
                _ => {}
            }
        }

        GeomAbsSurfaceType::Cylinder => {
            let cy = surf.cylinder();
            match type_func {
                TFunction::ContourStd => {
                    contana.perform_cylinder_dir(&cy, c.sfunc_ref().direction());
                }
                TFunction::ContourPrs => {
                    contana.perform_cylinder_eye(&cy, c.sfunc_ref().eye());
                }
                TFunction::DraftStd => {
                    let d = c.sfunc_ref().direction();
                    let a = c.sfunc_ref().angle();
                    contana.perform_cylinder_dir_angle(&cy, d, a);
                }
                _ => {}
            }
        }

        GeomAbsSurfaceType::Cone => {
            let co = surf.cone();
            match type_func {
                TFunction::ContourStd => {
                    contana.perform_cone_dir(&co, c.sfunc_ref().direction());
                }
                TFunction::ContourPrs => {
                    contana.perform_cone_eye(&co, c.sfunc_ref().eye());
                }
                TFunction::DraftStd => {
                    let d = c.sfunc_ref().direction();
                    let a = c.sfunc_ref().angle();
                    contana.perform_cone_dir_angle(&co, d, a);
                }
                _ => {}
            }
        }

        _ => {}
    }

    if typ_s != GeomAbsSurfaceType::Plane {
        if !contana.is_done() {
            return;
        }

        let nb_cont = contana.nb_contours();

        if contana.nb_contours() == 0 {
            c.set_done(true);
            return;
        }

        let typ_l = contana.type_contour();
        if typ_l == Curve2dType::Circle {
            let mut theline = Line::new();
            theline.set_value_circle(&contana.circle());
            let trans_circle = compute_transition_on_gp_circle(c.my_sfunc(), &contana.circle());
            theline.set_transition_on_s(trans_circle);
            c.slin_mut().push(theline);
        } else if typ_l == Curve2dType::Line {
            for i in 1..=nb_cont {
                let mut theline = Line::new();
                theline.set_value_line(&contana.line(i));
                let trans_line = compute_transition_on_gp_line(c.my_sfunc(), &contana.line(i));
                theline.set_transition_on_s(trans_line);
                c.slin_mut().push(theline);
            }

            /*
            if (typS == GeomAbs_Cone) { ... apex vertex branch (commented out
            in OCCT L2329-2341) ... }
            */
        }
    }

    if perform_sol_rst {
        {
            let (f, af, srst, _sins, _slin) = c.parts_mut();
            srst.perform(af, domain, tol_arc, tol_arc, false);
        }
        if !c.solrst().is_done() {
            return;
        }
        let nb_point_rst = c.solrst().nb_points();

        if nb_point_rst != 0 {
            let solrst = c.solrst().clone();
            let slin = c.slin_mut();
            put_points_on_line(&solrst, &*surf, slin);
        }

        if c.solrst().nb_segments() != 0 {
            let (f, _af, srst, _sins, slin) = c.parts_mut();
            process_segments(srst, slin, tol_arc, f, domain);
        }

        //-- lbr
        let nblinto = c.slin().len();
        let mut seq_to_destroy: Vec<usize> = Vec::new();

        for i in 1..=nblinto {
            if c.slin()[i - 1].type_contour() != IType::Restriction {
                // LineConstructor(slin, Domain, slin.ChangeValue(i), Surf):
                // the OCCT call passes the in-sequence line as L (read-only)
                // and appends the constructed pieces at the end.
                let mut theline = c.slin_mut()[i - 1].clone();
                line_constructor(c.slin_mut(), domain, &mut theline, &*surf);
                seq_to_destroy.push(i);
            }
        }
        for i in (0..seq_to_destroy.len()).rev() {
            let idx = seq_to_destroy[i];
            c.slin_mut().remove(idx - 1);
        }
    }

    c.set_done(true);
}
