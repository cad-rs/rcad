//! OCCT Draft_Modification_1.cxx — the Perform sections "Calculate new
//! edges" (L826-1443), "Calculate new vertices" (L1445-1636) and "the small
//! loop of validation/protection" (L1639-1775).  The split from
//! draft_modification_1.rs follows the single-file <2000-line rule only; the
//! OCCT section comments are the boundaries.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/
//!         Draft_Modification_1.cxx L826-1775

use glam::DVec3;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::geom_api::project_on_surf::ProjectPointOnSurf;
use rcad_kernel::base::proj_lib::geom_adaptor_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{Circle3, Curve3, CurveEval, Plane, Surface3};

use rcad_kernel::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::Orientation;

use crate::bop::int_tools::bean_face_intersector::IntCurveSurfaceHInter;
use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_parameter, brep_tool_tolerance, top_abs_reverse,
    top_exp_vertices_raw,
};
use super::draft_modification::brep_tool_pnt;

use super::draft_error_status::DraftErrorStatus;
use super::draft_modification::DraftModification;
use super::draft_modification_1::dir_is_parallel;
use super::draft_modification_1_b::{
    choose, elclib_circle_parameter, geom_curve2d_reverse, geom_curve_reverse, gp_circ_translate,
    parameter, smart_parameter, ExtremaExtCS, GeomConvertCompCurveToBSplineCurve, GeomIntIntSS,
};
use crate::geomalgo::geom_api_project_point_on_curve::GeomAPIProjectPointOnCurve;

impl DraftModification {
    /// OCCT Perform L826-1443 — "Calculate new edges".  The OCCT null-handle
    /// declarations (NewC/imin) surface as benign unused-assignment paths.
    #[allow(unused_assignments)]
    pub(crate) fn perform_calculate_edges(&mut self) {
        for ii in 1..=self.emap().extent() {
            // OCCT L829-831.
            let the_edge = self.emap().find_key(ii).clone();
            // OCCT L833-838: C = BRep_Tool::Curve(theEdge, L, f, l);
            // C = C->Transformed(L.Transformation()); — identity (#5).
            let c = brep_tool_curve(&the_edge).expect("BRep_Tool::Curve(theEdge)").0;

            // OCCT L840.
            let einf_needs_geometry = {
                let einf = self.emap_mut().change_from_index(ii);
                einf.new_geometry() && einf.geometry().is_none()
            };
            if einf_needs_geometry {
                // OCCT L842-843.
                let mut ptfixe = DVec3::ZERO;
                let is_tangent = self.emap_mut().change_from_index(ii).is_tangent(&mut ptfixe);
                let mut new_c: Option<Curve3> = None;
                if !is_tangent {
                    // OCCT L845-849.
                    let (first_face, second_face) = {
                        let einf = self.emap_mut().change_from_index(ii);
                        (einf.first_face().clone(), einf.second_face().clone())
                    };
                    let mut s1 = self.fmap().find_from_key(&first_face).geometry().cloned();
                    let mut s2 = self.fmap().find_from_key(&second_face).geometry().cloned();

                    // OCCT L851.
                    let mut detrompeur = 0i32;

                    // OCCT L853-862.
                    let fv = top_exp_vertices_raw(&the_edge).0.expect("TopExp::FirstVertex");
                    let lv = top_exp_vertices_raw(&the_edge).1.expect("TopExp::LastVertex");
                    let mut pmin = 0.0f64;
                    let prmfv = brep_tool_parameter(&fv, &the_edge);
                    let prmlv = brep_tool_parameter(&lv, &the_edge);
                    let pfv = c.point_at(prmfv);
                    let d1fv = c.derivative_at(prmfv);
                    let plv = c.point_at(prmlv);

                    // OCCT L864-871.
                    let tol_f1 = brep_tool_tolerance(&first_face);
                    let tol_f2 = brep_tool_tolerance(&second_face);
                    let s1_ref = s1.clone().expect("null S1");
                    let s2_ref = s2.clone().expect("null S2");
                    let proj1 = ProjectPointOnSurf::new_point_tol(pfv, &s1_ref, tol_f1);
                    let proj2 = ProjectPointOnSurf::new_point_tol(plv, &s1_ref, tol_f1);
                    let proj3 = ProjectPointOnSurf::new_point_tol(pfv, &s2_ref, tol_f2);
                    let proj4 = ProjectPointOnSurf::new_point_tol(plv, &s2_ref, tol_f2);

                    // OCCT L873-880.
                    if proj1.is_done() && proj2.is_done()
                        && proj1.lower_distance() <= CONFUSION
                        && proj2.lower_distance() <= CONFUSION
                    {
                        detrompeur = 1;
                    }
                    // OCCT L882-889.
                    if proj3.is_done() && proj4.is_done()
                        && proj3.lower_distance() <= CONFUSION
                        && proj4.lower_distance() <= CONFUSION
                    {
                        detrompeur = 2;
                    }

                    // OCCT L891-894.
                    let mut the_dir_extr = DVec3::ZERO;
                    let mut axis = Plane::new(DVec3::ZERO, DVec3::Z); // OCCT: gp_Ax3 Axis;
                    let mut the_new_curve: Option<Curve3> = None;
                    let mut k_part = false;

                    // OCCT L896-903.
                    if let Some(Surface3::Trimmed(t)) = &s1 {
                        s1 = Some(t.basis.as_ref().clone());
                    }
                    if let Some(Surface3::Trimmed(t)) = &s2 {
                        s2 = Some(t.basis.as_ref().clone());
                    }

                    // OCCT L905-922: the KPart detection.
                    let mut pc1 = true; // KPart on S1
                    if matches!(&s1, Some(Surface3::LinearExtrusion(_)))
                        && matches!(&s2, Some(Surface3::Plane(_)))
                    {
                        k_part = true;
                        axis = match &s2 {
                            Some(Surface3::Plane(p)) => *p,
                            _ => unreachable!(),
                        };
                        the_new_curve = match &s1 {
                            Some(Surface3::LinearExtrusion(le)) => {
                                Some(le.profile.as_ref().clone())
                            }
                            _ => unreachable!(),
                        };
                        the_dir_extr = match &s1 {
                            Some(Surface3::LinearExtrusion(le)) => le.direction,
                            _ => unreachable!(),
                        };
                    } else if matches!(&s2, Some(Surface3::LinearExtrusion(_)))
                        && matches!(&s1, Some(Surface3::Plane(_)))
                    {
                        k_part = true;
                        pc1 = false;
                        axis = match &s1 {
                            Some(Surface3::Plane(p)) => *p,
                            _ => unreachable!(),
                        };
                        the_new_curve = match &s2 {
                            Some(Surface3::LinearExtrusion(le)) => {
                                Some(le.profile.as_ref().clone())
                            }
                            _ => unreachable!(),
                        };
                        the_dir_extr = match &s2 {
                            Some(Surface3::LinearExtrusion(le)) => le.direction,
                            _ => unreachable!(),
                        };
                    }
                    // OCCT L923-936: the circle restriction.
                    let mut a_circ: Option<Circle3> = None;
                    if k_part {
                        // very temporary on circles !!!
                        match &the_new_curve {
                            Some(Curve3::Circle(cir)) => {
                                a_circ = Some(*cir);
                                let ax_of_circ = cir.normal;
                                k_part = dir_is_parallel(ax_of_circ, axis.normal, ANGULAR_TOL);
                            }
                            _ => {
                                k_part = false;
                            }
                        }
                    }

                    // OCCT L938-939.
                    let mut imin = 0usize;
                    let mut i2s = GeomIntIntSS::new();
                    if k_part {
                        // OCCT L941-959: the direct calculation of NewC.
                        let a_circ = a_circ.expect("KPart circle");
                        let a_local_real = (axis.origin - a_circ.center).dot(axis.normal);
                        let cos_ = the_dir_extr.dot(axis.normal);
                        let vv = the_dir_extr * (a_local_real / cos_);
                        // OCCT L947: newC = TheNewCurve->Translated(VV) — the
                        // circle translation.
                        new_c = Some(Curve3::Circle(gp_circ_translate(&a_circ, vv)));
                        // OCCT L949-950: the p-curve.
                        let l2d = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                            origin: glam::DVec2::new(0.0, a_local_real / cos_),
                            direction: glam::DVec2::X, // gp::DX2d()
                        });
                        // OCCT L952-959.
                        let einf = self.emap_mut().change_from_index(ii);
                        if pc1 {
                            *einf.change_first_pc() = Some(l2d);
                        } else {
                            *einf.change_second_pc() = Some(l2d);
                        }
                    } else {
                        // OCCT L963-964.
                        let first_face = self.emap().find_from_index(ii).first_face().clone();
                        let second_face = self.emap().find_from_index(ii).second_face().clone();
                        s1 = self.fmap().find_from_key(&first_face).geometry().cloned();
                        s2 = self.fmap().find_from_key(&second_face).geometry().cloned();

                        // OCCT L970: i2s.Perform(S1, S2, Confusion, true, false, false);
                        let s1r = s1.clone().expect("null S1");
                        let s2r = s2.clone().expect("null S2");
                        i2s.perform(&s1r, &s2r, CONFUSION, true, false, false);

                        // OCCT L972-977.
                        if !i2s.is_done() || i2s.nb_lines() <= 0 {
                            self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                            self.bad_shape_mut().clone_from(&the_edge);
                            return;
                        }

                        // OCCT L979.
                        let mut glob2_min = f64::MAX; // RealLast()

                        // OCCT L984-1158: the non-BSpline branch.
                        let first_is_bspline = matches!(i2s.line(1), Curve3::BSpline(_));
                        if !first_is_bspline {
                            let mut dist2_min = f64::MAX;
                            imin = 0;
                            for i in 1..=i2s.nb_lines() {
                                // OCCT L990-991: TheCurve.Load(i2s.Line(i));
                                // Extrema_ExtPC myExtPC(pfv, TheCurve) — the
                                // two-arg ctor over the full domain, the
                                // default theTolF is 1.0e-10.
                                let the_curve = i2s.line(i);
                                let dom = the_curve.default_domain();
                                let a_adaptor = GeomCurveAdaptor::new(the_curve.clone());
                                let a_tool =
                                    CurveToolHandle::for_curve3(&the_curve, &a_adaptor, &a_adaptor);
                                let my_ext_pc = ExtremaExtPC::new_point_curve(pfv, &a_tool, 1.0e-10);

                                let mut locpmin = 0.0f64;
                                if my_ext_pc.is_done() {
                                    // OCCT L996-1000.
                                    if my_ext_pc.nb_ext() >= 1 {
                                        dist2_min = my_ext_pc.square_distance(1);
                                        locpmin = my_ext_pc.point(1).param;
                                    }
                                    // OCCT L1001-1084.
                                    if my_ext_pc.nb_ext() == 2 && dist2_min > SQUARE_CONFUSION {
                                        // to avoid incorrectly choosing the image
                                        // of the first vertex of the initial edge
                                        let d1_2 = my_ext_pc.square_distance(1);
                                        let d2_2 = my_ext_pc.square_distance(2);
                                        if d1_2 > 1.21 * d2_2 {
                                            dist2_min = my_ext_pc.square_distance(2);
                                            locpmin = my_ext_pc.point(2).param;
                                        } else if d2_2 > 1.21 * d1_2 {
                                            dist2_min = my_ext_pc.square_distance(1);
                                            locpmin = my_ext_pc.point(1).param;
                                        } else {
                                            // OCCT L1019-1083.
                                            let pfvpar = my_ext_pc.point(1).param;
                                            let plvpar = my_ext_pc.point(2).param;
                                            let newc = i2s.line(i);

                                            let pfvprim = newc.point_at(pfvpar);
                                            let plvprim = newc.point_at(plvpar);

                                            // OCCT L1028-1044.
                                            let mut the_surf: Option<Surface3> = None;
                                            if detrompeur == 1 {
                                                if let Some(Surface3::Trimmed(t)) = &s1 {
                                                    s1 = Some(t.basis.as_ref().clone());
                                                }
                                                the_surf = s1.clone();
                                            } else if detrompeur == 2 {
                                                if let Some(Surface3::Trimmed(t)) = &s2 {
                                                    s2 = Some(t.basis.as_ref().clone());
                                                }
                                                the_surf = s2.clone();
                                            }
                                            // OCCT L1045-1082.
                                            if detrompeur != 0 && detrompeur != 4 {
                                                // OCCT L1050-1069: the OCCT text
                                                // tests theSurf but downcasts S2
                                                // — translated literally.
                                                let is_plane =
                                                    matches!(the_surf, Some(Surface3::Plane(_)));
                                                let is_cyl = matches!(
                                                    the_surf,
                                                    Some(Surface3::Cylinder(_))
                                                );
                                                let (mut ul, mut vl) = (0.0f64, 0.0f64);
                                                let (mut uf, mut vf) = (0.0f64, 0.0f64);
                                                let (mut ufprim, mut ulprim) = (0.0f64, 0.0f64);
                                                let (mut vfprim, mut vlprim) = (0.0f64, 0.0f64);
                                                if is_plane {
                                                    let pl = match &s2 {
                                                        Some(Surface3::Plane(p)) => *p,
                                                        _ => unreachable!(),
                                                    };
                                                    (ul, vl) =
                                                        super::draft_modification_1::
                                                            elslib_parameters_plane(&pl, plv);
                                                    (uf, vf) =
                                                        super::draft_modification_1::
                                                            elslib_parameters_plane(&pl, pfv);
                                                    (ulprim, vlprim) = super::draft_modification_1::
                                                        elslib_parameters_plane(&pl, plvprim);
                                                    (ufprim, vfprim) = super::draft_modification_1::
                                                        elslib_parameters_plane(&pl, pfvprim);
                                                } else if is_cyl {
                                                    let cy = match &s2 {
                                                        Some(Surface3::Cylinder(c)) => *c,
                                                        _ => unreachable!(),
                                                    };
                                                    (ul, vl) = super::draft_modification_1::
                                                        elslib_parameters_cyl(&cy, plv);
                                                    (uf, vf) = super::draft_modification_1::
                                                        elslib_parameters_cyl(&cy, pfv);
                                                    (ulprim, vlprim) = super::draft_modification_1::
                                                        elslib_parameters_cyl(&cy, plvprim);
                                                    (ufprim, vfprim) = super::draft_modification_1::
                                                        elslib_parameters_cyl(&cy, pfvprim);
                                                } else {
                                                    detrompeur = 4;
                                                }

                                                // OCCT L1071-1081.
                                                if detrompeur == 1 || detrompeur == 2 {
                                                    let v1 = glam::DVec2::new(
                                                        ul - ufprim,
                                                        vl - vfprim,
                                                    );
                                                    let norm = glam::DVec2::new(
                                                        vf - vfprim,
                                                        ufprim - uf,
                                                    );
                                                    let v2 = glam::DVec2::new(
                                                        ulprim - ufprim,
                                                        vlprim - vfprim,
                                                    );
                                                    if (v1.dot(norm)) * (v2.dot(norm)) < 0.0 {
                                                        dist2_min = my_ext_pc.square_distance(2);
                                                        locpmin = my_ext_pc.point(2).param;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    // OCCT L1085-1098.
                                    if my_ext_pc.nb_ext() == 1
                                        || my_ext_pc.nb_ext() > 2
                                        || detrompeur == 4
                                    {
                                        dist2_min = my_ext_pc.square_distance(1);
                                        locpmin = my_ext_pc.point(1).param;
                                        for j in 2..=my_ext_pc.nb_ext() {
                                            let dist2 = my_ext_pc.square_distance(j);
                                            if dist2 < dist2_min {
                                                dist2_min = dist2;
                                                locpmin = my_ext_pc.point(j).param;
                                            }
                                        }
                                    } else if my_ext_pc.nb_ext() < 1 {
                                        // OCCT L1099-1114:
                                        // myExtPC.TrimmedSquareDistances(dist1_2, dist2_2, p1b, p2b).
                                        let (dist1_2, dist2_2, _p1b, _p2b) =
                                            my_ext_pc.trimmed_square_distances();
                                        if dist1_2 < dist2_2 {
                                            dist2_min = dist1_2;
                                            locpmin = dom[0]; // TheCurve.FirstParameter()
                                        } else {
                                            dist2_min = dist2_2;
                                            locpmin = dom[1]; // TheCurve.LastParameter()
                                        }
                                    }

                                    // OCCT L1116-1121.
                                    if dist2_min < glob2_min {
                                        glob2_min = dist2_min;
                                        imin = i;
                                        pmin = locpmin;
                                    }
                                }
                            }
                            // OCCT L1124-1129.
                            if imin == 0 {
                                self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                                self.bad_shape_mut().clone_from(&the_edge);
                                return;
                            }

                            // OCCT L1131-1139.
                            let mut newc = i2s.line(imin);
                            let newd1 = newc.derivative_at(pmin);
                            let ya_rev = d1fv.dot(newd1) < 0.0;
                            if ya_rev {
                                newc = geom_curve_reverse(&newc);
                            }
                            new_c = Some(newc);

                            // OCCT L1141-1157.
                            let einf = self.emap_mut().change_from_index(ii);
                            if i2s.has_line_on_s1(imin) {
                                let mut pc = i2s.line_on_s1(imin);
                                if ya_rev {
                                    pc = geom_curve2d_reverse(&pc);
                                }
                                *einf.change_first_pc() = Some(pc);
                            }
                            if i2s.has_line_on_s2(imin) {
                                let mut pc = i2s.line_on_s2(imin);
                                if ya_rev {
                                    pc = geom_curve2d_reverse(&pc);
                                }
                                *einf.change_second_pc() = Some(pc);
                            }
                        } else {
                            // OCCT L1159-1327: i2s.Line(1) is a BSpline.
                            // OCCT L1162: Candidates (line index + curve — the
                            // index carries the OCCT handle identity).
                            let mut candidates: Vec<(usize, Curve3)> = Vec::new();
                            let s1_is_cyl_or_con = matches!(
                                &s1,
                                Some(Surface3::Cylinder(_)) | Some(Surface3::Cone(_))
                            );
                            if s1_is_cyl_or_con {
                                for i in 1..=i2s.nb_lines() {
                                    let a_curve = i2s.line(i);
                                    let pnt =
                                        a_curve.point_at(a_curve.default_domain()[0]);
                                    let projector = ProjectPointOnSurf::new_point_tol(
                                        pnt,
                                        s1.as_ref().unwrap(),
                                        CONFUSION,
                                    );
                                    let (u, _v) = projector.lower_distance_parameters();
                                    // OCCT L1173-1189.
                                    let two_pi = 2.0 * std::f64::consts::PI;
                                    if u.abs() <= CONFUSION || (u - two_pi).abs() <= CONFUSION {
                                        candidates.push((i, a_curve));
                                    } else {
                                        let pnt2 =
                                            a_curve.point_at(a_curve.default_domain()[1]);
                                        let projector = ProjectPointOnSurf::new_point_tol(
                                            pnt2,
                                            s1.as_ref().unwrap(),
                                            CONFUSION,
                                        );
                                        let (u2, _v2) = projector.lower_distance_parameters();
                                        if u2.abs() <= CONFUSION
                                            || (u2 - two_pi).abs() <= CONFUSION
                                        {
                                            let rev = geom_curve_reverse(&a_curve);
                                            candidates.push((i, rev));
                                        }
                                    }
                                }
                                // OCCT L1192-1201.
                                if candidates.is_empty() {
                                    for i in 1..=i2s.nb_lines() {
                                        candidates.push((i, i2s.line(i)));
                                    }
                                }
                            } else {
                                // OCCT L1203-1209.
                                for i in 1..=i2s.nb_lines() {
                                    candidates.push((i, i2s.line(i)));
                                }
                            }

                            // OCCT L1211-1230.
                            let mut first_curve: Option<(usize, Curve3)> = None;
                            if candidates.len() > 1 {
                                let mut dist_min = f64::INFINITY; // Precision::Infinite()
                                for (idx, a_curve) in &candidates {
                                    let pnt = a_curve.point_at(a_curve.default_domain()[0]);
                                    let dist = pnt.distance(pfv);
                                    if dist - dist_min < -CONFUSION {
                                        dist_min = dist;
                                        first_curve = Some((*idx, a_curve.clone()));
                                    }
                                }
                            } else if let Some((idx, c1)) = candidates.first() {
                                first_curve = Some((*idx, c1.clone()));
                            }

                            // OCCT L1232-1240: glueing.
                            let mut curves: Vec<(usize, Curve3)> = Vec::new();
                            for i in 1..=i2s.nb_lines() {
                                if first_curve.as_ref().map(|(fi, _)| *fi) != Some(i) {
                                    curves.push((i, i2s.line(i)));
                                }
                            }

                            let mut to_glue: Vec<Curve3> = Vec::new();
                            let fc = first_curve.as_ref().expect("FirstCurve").1.clone();
                            let mut end_point = fc.point_at(fc.default_domain()[1]);
                            let mut added = true;
                            while added {
                                added = false;
                                let mut ci = 0usize;
                                while ci < curves.len() {
                                    let a_curve = curves[ci].1.clone();
                                    let pfirst =
                                        a_curve.point_at(a_curve.default_domain()[0]);
                                    let plast =
                                        a_curve.point_at(a_curve.default_domain()[1]);
                                    if pfirst.distance(end_point) <= CONFUSION {
                                        to_glue.push(a_curve);
                                        end_point = plast;
                                        curves.remove(ci);
                                        added = true;
                                        break;
                                    }
                                    if plast.distance(end_point) <= CONFUSION {
                                        let rev = geom_curve_reverse(&a_curve);
                                        to_glue.push(rev);
                                        end_point = pfirst;
                                        curves.remove(ci);
                                        added = true;
                                        break;
                                    }
                                    ci += 1;
                                }
                            }

                            // OCCT L1274-1279.
                            if first_curve.is_none() {
                                self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                                self.bad_shape_mut().clone_from(&the_edge);
                                return;
                            }

                            // OCCT L1281-1290: the concatenation (the
                            // GeomConvert carrier; the Add tail is deferred
                            // until the GeomConvert batch lands).
                            let fc_curve = first_curve.as_ref().unwrap().1.clone();
                            let first_bspline = match &fc_curve {
                                Curve3::BSpline(b) => b.clone(),
                                _ => panic!("Geom_BSplineCurve cast of FirstCurve"),
                            };
                            let mut concat =
                                GeomConvertCompCurveToBSplineCurve::new(&first_bspline);
                            for g in &to_glue {
                                match g {
                                    Curve3::BSpline(b) => {
                                        concat.add(b, CONFUSION, true);
                                    }
                                    _ => panic!("Geom_BSplineCurve cast of ToGlue"),
                                }
                            }
                            let mut newc = Curve3::BSpline(concat.bspline_curve());

                            // OCCT L1292-1306: TheCurve.Load(newC);
                            // Extrema_ExtPC myExtPC(pfv, TheCurve) — the
                            // two-arg ctor, the default theTolF 1.0e-10.
                            let a_adaptor = GeomCurveAdaptor::new(newc.clone());
                            let a_tool = CurveToolHandle::for_curve3(&newc, &a_adaptor, &a_adaptor);
                            let my_ext_pc = ExtremaExtPC::new_point_curve(pfv, &a_tool, 1.0e-10);
                            let mut dist2_min = f64::MAX;
                            for i in 1..=my_ext_pc.nb_ext() {
                                // OCCT L1297: if (myExtPC.IsMin(i)).
                                if my_ext_pc.is_min(i) {
                                    let dist2 = my_ext_pc.square_distance(i);
                                    if dist2 < dist2_min {
                                        dist2_min = dist2;
                                        pmin = my_ext_pc.point(i).param;
                                    }
                                }
                            }
                            // OCCT L1307-1313.
                            let newd1 = newc.derivative_at(pmin);
                            let ya_rev = d1fv.dot(newd1) < 0.0;
                            if ya_rev {
                                newc = geom_curve_reverse(&newc);
                            }
                            new_c = Some(newc);
                            // OCCT L1314-1326: the commented-out p-curve block.
                        }

                        // OCCT L1329.
                        let einf = self.emap_mut().change_from_index(ii);
                        einf.set_tolerance(einf.tolerance().max(i2s.tol_reached_3d()));
                    }
                    // OCCT L1331: End step KPart.
                } else {
                    // OCCT L1332-1424: case of tangency.
                    let f1 = self.emap().find_from_index(ii).first_face().clone();
                    let f2 = self.emap().find_from_index(ii).second_face().clone();

                    // OCCT L1337-1352.
                    let mut a_local_s1 = self.fmap().find_from_key(&f1).geometry().cloned();
                    let mut a_local_s2 = self.fmap().find_from_key(&f2).geometry().cloned();
                    if a_local_s1.is_none() || a_local_s2.is_none() {
                        self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                        self.bad_shape_mut().clone_from(&the_edge);
                        return;
                    }
                    if let Some(Surface3::Trimmed(t)) = &a_local_s1 {
                        a_local_s1 = Some(t.basis.as_ref().clone());
                    }
                    if let Some(Surface3::Trimmed(t)) = &a_local_s2 {
                        a_local_s2 = Some(t.basis.as_ref().clone());
                    }

                    // OCCT L1354-1404: dirextr.
                    let mut dirextr = DVec3::ZERO;
                    if let Some(Surface3::Cylinder(cyl)) = &a_local_s1 {
                        dirextr = cyl.axis;
                    } else if let Some(Surface3::LinearExtrusion(le)) = &a_local_s1 {
                        dirextr = le.direction;

                        // OCCT L1369-1379: the p-curve on S1.
                        if let Curve3::Circle(cir) = le.profile.as_ref() {
                            let u = elclib_circle_parameter(cir, ptfixe);
                            let pc1 = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                                origin: glam::DVec2::new(u, 0.0),
                                direction: glam::DVec2::Y, // gp::DY2d()
                            });
                            let einf = self.emap_mut().change_from_index(ii);
                            *einf.change_first_pc() = Some(pc1);
                        }
                    } else if let Some(Surface3::Cylinder(cyl)) = &a_local_s2 {
                        dirextr = cyl.axis;
                    } else if let Some(Surface3::LinearExtrusion(le)) = &a_local_s2 {
                        dirextr = le.direction;

                        // OCCT L1394-1403: the p-curve on S2.
                        if let Curve3::Circle(cir) = le.profile.as_ref() {
                            let u = elclib_circle_parameter(cir, ptfixe);
                            let pc2 = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                                origin: glam::DVec2::new(u, 0.0),
                                direction: glam::DVec2::Y,
                            });
                            let einf = self.emap_mut().change_from_index(ii);
                            *einf.change_second_pc() = Some(pc2);
                        }
                    }
                    // OCCT L1405.
                    let mut newc = Curve3::Line(rcad_kernel::geom::Line3::new(ptfixe, dirextr));

                    // OCCT L1407-1423.
                    let d1fv = c.derivative_at(0.0);
                    let newd1 = newc.derivative_at(0.0);
                    let ya_rev = d1fv.dot(newd1) < 0.0;
                    if ya_rev {
                        newc = geom_curve_reverse(&newc);
                        let einf = self.emap_mut().change_from_index(ii);
                        if einf.first_pc().is_some() {
                            let pc = einf.first_pc().cloned().unwrap();
                            *einf.change_first_pc() = Some(geom_curve2d_reverse(&pc));
                        }
                        if einf.second_pc().is_some() {
                            let pc = einf.second_pc().cloned().unwrap();
                            *einf.change_second_pc() = Some(geom_curve2d_reverse(&pc));
                        }
                    }
                    new_c = Some(newc);
                }

                // OCCT L1426-1431.
                if let Some(Curve3::Trimmed(t)) = &new_c {
                    new_c = Some(t.curve.as_ref().clone());
                }
                *self.emap_mut().change_from_index(ii).change_geometry() = new_c;
            } else if !self.emap().find_from_index(ii).new_geometry() {
                // OCCT L1433-1442: set the existing curve 3D.
                let mut cc = c;
                if let Curve3::Trimmed(t) = &cc {
                    cc = t.curve.as_ref().clone();
                }
                *self.emap_mut().change_from_index(ii).change_geometry() = Some(cc);
            }
        }
    }

    /// OCCT Perform L1445-1636 — "Calculate new vertices".
    pub(crate) fn perform_calculate_vertices(&mut self) {
        for ii in 1..=self.vmap().extent() {
            let tvv = self.vmap().find_key(ii).clone();

            // OCCT L1455-1457.
            let vtori = brep_tool_pnt(&tvv);
            let chosen = {
                // The OCCT Choose call passes myFMap/myEMap/Vinf as separate
                // objects; the rcad maps are split-borrowed here and the Vinf
                // element is taken out for the call (architecture difference
                // #1 of draft_modification.rs).
                let (fmap, emap, vmap) = self.draft_maps_split();
                let mut vinf_local = std::mem::take(vmap.change_from_index(ii));
                let chosen = choose(fmap, emap, &tvv, &mut vinf_local);
                *vmap.change_from_index(ii) = vinf_local;
                chosen
            };
            let (hac, has) = match chosen {
                Some(pair) => pair,
                None => {
                    // OCCT L1460-1539: no concerted edge — the alignment of
                    // two consecutive edges.
                    let vinf_edges: Vec<Shape> = {
                        let vinf = self.vmap_mut().change_from_index(ii);
                        vinf.init_edge_iterator();
                        let mut v = Vec::new();
                        while vinf.more_edge() {
                            v.push(vinf.edge());
                            vinf.next_edge();
                        }
                        v
                    };
                    if vinf_edges.is_empty() {
                        // OCCT L1541-1543.
                        self.set_err_stat(DraftErrorStatus::VertexRecomputation);
                        self.bad_shape_mut().clone_from(&tvv);
                        return;
                    }
                    // OCCT L1466-1472.
                    let edg1 = vinf_edges[0].clone();
                    let einf1_geom = self
                        .emap()
                        .find_from_key(&edg1)
                        .geometry()
                        .cloned()
                        .expect("null Geometry (Edg1)");
                    // OCCT L1471-1472: the projection of the original point on
                    // the new curve (patch).
                    let projector = GeomAPIProjectPointOnCurve::new_point_curve(vtori, &einf1_geom);
                    let mut pvt = projector.nearest_point();
                    // OCCT L1484.
                    let dion = pvt.distance_squared(vtori);
                    // OCCT L1485-1521.
            if vinf_edges.len() > 1 {
                        let edg2 = vinf_edges[1].clone();
                        let einf2_geom = self
                            .emap()
                            .find_from_key(&edg2)
                            .geometry()
                            .cloned()
                            .expect("null Geometry (Edg2)");
                        let prm2 = {
                            let vinf = self.vmap_mut().change_from_index(ii);
                            vinf.parameter(&edg2)
                        };
                        let opvt = einf2_geom.point_at(prm2);
                        if opvt.distance_squared(vtori) < dion {
                            pvt = opvt;
                        }
                        // OCCT L1508-1520.
                        let mut done = 0i32;
                        let param = parameter(&einf2_geom, pvt, &mut done);
                        if done != 0 {
                            let einf2 = self.emap().find_from_key(&edg2);
                            let s1 = self
                                .fmap()
                                .find_from_key(einf2.first_face())
                                .geometry()
                                .cloned()
                                .expect("null S1");
                            let s2 = self
                                .fmap()
                                .find_from_key(einf2.second_face())
                                .geometry()
                                .cloned()
                                .expect("null S2");
                            let tol = brep_tool_tolerance(&edg2);
                            let einf2 = self.emap_mut().change_from_key(&edg2);
                            let prm = smart_parameter(einf2, tol, pvt, done, &s1, &s2);
                            let vinf = self.vmap_mut().change_from_index(ii);
                            *vinf.change_parameter(&edg2) = prm;
                        } else {
                            let vinf = self.vmap_mut().change_from_index(ii);
                            *vinf.change_parameter(&edg2) = param;
                        }
                    }
                    // OCCT L1523-1537.
                    {
                        let mut done = 0i32;
                        let param = parameter(&einf1_geom, pvt, &mut done);
                        if done != 0 {
                            let einf1 = self.emap().find_from_key(&edg1);
                            let s1 = self
                                .fmap()
                                .find_from_key(einf1.first_face())
                                .geometry()
                                .cloned()
                                .expect("null S1");
                            let s2 = self
                                .fmap()
                                .find_from_key(einf1.second_face())
                                .geometry()
                                .cloned()
                                .expect("null S2");
                            let tol = brep_tool_tolerance(&edg1);
                            let einf1 = self.emap_mut().change_from_key(&edg1);
                            let prm = smart_parameter(einf1, tol, pvt, done, &s1, &s2);
                            let vinf = self.vmap_mut().change_from_index(ii);
                            *vinf.change_parameter(&edg1) = prm;
                        } else {
                            let vinf = self.vmap_mut().change_from_index(ii);
                            *vinf.change_parameter(&edg1) = param;
                        }
                    }
                    // OCCT L1523: Vinf.ChangeGeometry() = pvt;
                    *self.vmap_mut().change_from_index(ii).change_geometry() = pvt;
                    // OCCT L1538: continue;
                    continue;
                }
            };

            // OCCT L1546-1553.
            let mut myintcs = IntCurveSurfaceHInter::new();
            myintcs.perform(&hac, &has);
            if !myintcs.is_done() {
                self.set_err_stat(DraftErrorStatus::VertexRecomputation);
                self.bad_shape_mut().clone_from(&tvv);
                return;
            }

            // OCCT L1555-1601.
            let pvt;
            let nbsol = myintcs.nb_points();
            if nbsol <= 0 {
                // OCCT L1561-1568: the Extrema_ExtCS fallback (GAP carrier:
                // IsDone=false keeps the OCCT error branch).
                let extr = ExtremaExtCS::new(&hac, &has, PCONFUSION, PCONFUSION);
                if !extr.is_done() || extr.nb_ext() == 0 {
                    self.set_err_stat(DraftErrorStatus::VertexRecomputation);
                    self.bad_shape_mut().clone_from(&tvv);
                    return;
                }
                pvt = DVec3::ZERO;
            } else {
                let mut disref = f64::MAX;
                let mut iref = 0usize;
                for i in 1..=nbsol {
                    let distemp = myintcs.base.point(i).pnt().distance_squared(vtori);
                    if distemp < disref {
                        disref = distemp;
                        iref = i;
                    }
                }
                pvt = myintcs.base.point(iref).pnt();
            }

            // OCCT L1603.
            *self.vmap_mut().change_from_index(ii).change_geometry() = pvt;

            // OCCT L1605-1635.
            let edges: Vec<Shape> = {
                let vinf = self.vmap_mut().change_from_index(ii);
                vinf.init_edge_iterator();
                let mut v = Vec::new();
                while vinf.more_edge() {
                    v.push(vinf.edge());
                    vinf.next_edge();
                }
                v
            };
            for edg in edges {
                let initpar = {
                    let vinf = self.vmap_mut().change_from_index(ii);
                    vinf.parameter(&edg)
                };
                let einf_geom = self
                    .emap()
                    .find_from_key(&edg)
                    .geometry()
                    .cloned()
                    .expect("null Geometry (Edg)");
                // OCCT L1613.
                let mut done = 0i32;
                let param = parameter(&einf_geom, pvt, &mut done);
                if done != 0 {
                    let einf = self.emap().find_from_key(&edg);
                    let s1 = self
                        .fmap()
                        .find_from_key(einf.first_face())
                        .geometry()
                        .cloned()
                        .expect("null S1");
                    let s2 = self
                        .fmap()
                        .find_from_key(einf.second_face())
                        .geometry()
                        .cloned()
                        .expect("null S2");
                    let tol = brep_tool_tolerance(&edg);
                    let einf = self.emap_mut().change_from_key(&edg);
                    let prm = smart_parameter(einf, tol, pvt, done, &s1, &s2);
                    let vinf = self.vmap_mut().change_from_index(ii);
                    *vinf.change_parameter(&edg) = prm;
                } else {
                    // OCCT L1622-1632.
                    if (initpar - param).abs() > PCONFUSION {
                        let a_c = brep_tool_curve(&edg);
                        if let Some((Curve3::Trimmed(_), _, _)) = a_c {
                            self.emap_mut().change_from_key(&edg).set_new_geometry(true);
                        }
                    }
                    let vinf = self.vmap_mut().change_from_index(ii);
                    *vinf.change_parameter(&edg) = param;
                }
            }
        }
    }

    /// OCCT Perform L1639-1775 — the small loop of validation/protection.
    pub(crate) fn perform_validation(&mut self) {
        for i in 1..=self.emap().extent() {
            let edg = self.emap().find_key(i).clone();

            // OCCT L1645-1651.
            let raw = top_exp_vertices_raw(&edg);
            let mut vf = raw.0.unwrap_or_else(Shape::null);
            let mut vl = raw.1.unwrap_or_else(Shape::null);
            if edg.orientation == Orientation::Reversed {
                vf.orientation = top_abs_reverse(vf.orientation);
                vl.orientation = top_abs_reverse(vl.orientation);
            }

            // OCCT L1653-1722.
            if self.vmap().contains(&vf) && self.vmap().contains(&vl) {
                // OCCT L1660-1661: compare the directions of the source edge
                // and of the selected part of the intersection edge.
                let a_par_f = self.vmap_mut().change_from_key(&vf).parameter(&edg);
                let a_par_l = self.vmap_mut().change_from_key(&vl).parameter(&edg);

                if a_par_l < a_par_f {
                    // OCCT L1665-1675.
                    let an_int_curv = self.emap().find_from_key(&edg).geometry().cloned();
                    let (a_scurve, _, _) = brep_tool_curve(&edg).expect("BRep_Tool::Curve(edg)");
                    let a_pf = brep_tool_pnt(&vf);
                    let a_pl = brep_tool_pnt(&vl);
                    let a_dir_of = a_scurve.derivative_at(brep_tool_parameter(&vf, &edg));
                    let a_dir_ol = a_scurve.derivative_at(brep_tool_parameter(&vl, &edg));
                    let an_int_curv = an_int_curv.expect("null anIntCurv");
                    let a_dir_nf = an_int_curv.derivative_at(a_par_f);
                    let a_dir_nl = an_int_curv.derivative_at(a_par_l);

                    // OCCT L1677-1700: the normalizations.
                    let mut a_dir_nf = a_dir_nf;
                    let mut a_sq_magn = a_dir_nf.length_squared();
                    if a_sq_magn > SQUARE_CONFUSION {
                        a_dir_nf /= a_sq_magn.sqrt();
                    }
                    let mut a_dir_nl = a_dir_nl;
                    a_sq_magn = a_dir_nl.length_squared();
                    if a_sq_magn > SQUARE_CONFUSION {
                        a_dir_nl /= a_sq_magn.sqrt();
                    }
                    let mut a_dir_of = a_dir_of;
                    a_sq_magn = a_dir_of.length_squared();
                    if a_sq_magn > SQUARE_CONFUSION {
                        a_dir_of /= a_sq_magn.sqrt();
                    }
                    let mut a_dir_ol = a_dir_ol;
                    a_sq_magn = a_dir_ol.length_squared();
                    if a_sq_magn > SQUARE_CONFUSION {
                        a_dir_ol /= a_sq_magn.sqrt();
                    }

                    // OCCT L1702-1703.
                    let a_cos_f = a_dir_nf.dot(a_dir_of);
                    let a_cos_l = a_dir_nl.dot(a_dir_ol);
                    let a_cos_max = if a_cos_f.abs() > a_cos_l.abs() {
                        a_cos_f
                    } else {
                        a_cos_l
                    };

                    // OCCT L1705-1721.
                    if a_cos_max < 0.0 {
                        let mut an_err = 0i32;
                        let rev = geom_curve_reverse(&an_int_curv);
                        *self.emap_mut().change_from_key(&edg).change_geometry() =
                            Some(rev.clone());
                        let a_par = parameter(&rev, a_pf, &mut an_err);
                        if an_err == 0 {
                            *self.vmap_mut().change_from_key(&vf).change_parameter(&edg) = a_par;
                        }
                        let a_par = parameter(&rev, a_pl, &mut an_err);
                        if an_err == 0 {
                            *self.vmap_mut().change_from_key(&vl).change_parameter(&edg) = a_par;
                        }
                    }
                }
            }

            // OCCT L1724-1732.
            let mut pf = 0.0f64;
            let mut pl = 0.0f64;
            let mut tolerance = 0.0f64;
            if !self.new_parameter(&vf, &edg, &mut pf, &mut tolerance) {
                pf = brep_tool_parameter(&vf, &edg);
            }
            if !self.new_parameter(&vl, &edg, &mut pl, &mut tolerance) {
                pl = brep_tool_parameter(&vl, &edg);
            }
            // OCCT L1733-1766.
            if pl <= pf {
                let the_curve = self
                    .emap()
                    .find_from_key(&edg)
                    .geometry()
                    .cloned()
                    .expect("null theCurve");
                if the_curve.is_closed() {
                    // pf >= pl
                    let dom = the_curve.default_domain();
                    let first_par = dom[0];
                    let last_par = dom[1];
                    let pconf = PCONFUSION;
                    if (pf - last_par).abs() <= pconf {
                        pf = first_par;
                    } else if (pl - first_par).abs() <= pconf {
                        pl = last_par;
                    }
                    if pl <= pf {
                        pl += last_par - first_par;
                    }
                }
                if pl <= pf {
                    self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                    self.bad_shape_mut().clone_from(&edg);
                    return;
                }
            }
            // OCCT L1767-1774.
            if self.vmap().contains(&vf) {
                *self.vmap_mut().change_from_key(&vf).change_parameter(&edg) = pf;
            }
            if self.vmap().contains(&vl) {
                *self.vmap_mut().change_from_key(&vl).change_parameter(&edg) = pl;
            }
        }
    }
}

/// The ANGULAR tolerance re-import (kept local for the section file).
const ANGULAR_TOL: f64 = rcad_kernel::precision::ANGULAR;
const _: Option<&'static str> = Some("BRepAdaptorSurface import anchor");
