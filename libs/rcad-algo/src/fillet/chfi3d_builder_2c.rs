//! OCCT ChFi3d_Builder_2.cxx — 1:1 translation, extremities part (Stage 1f).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_Builder_2.cxx, L1629-2203 (moved here from
//!         chfi3d_builder_2.rs to keep each file under 2000 lines).
//!
//! Coverage: ChFi3d_SingularExtremity, IsFree, ChFi3d_MakeExtremities,
//! ChFi3d_Purge, InsertAfter, RemoveSD, InsertBefore.

use rcad_kernel::geom::Curve2dEval as _;
use rcad_kernel::topo::topods::Shape;
use rcad_kernel::topods;

use super::chfi3d::{chfi3d_index_of_surf_data, chfi3d_index_point_in_ds};
use super::chfi3d_builder_0::{
    chfi3d_bound_surf, chfi3d_compute_arete, chfi3d_compute_pcurv_2pt, chfi3d_enlarge_box_dstr,
    chfi3d_enlarge_box_edge_faces, chfi3d_same_parameter, chfi3d_set_point_tolerance, BndBox,
};
use super::chfi3d_builder_2::chfi3d_coupe_par_plan;
use super::chfi3d_ds::{TopOpeBRepDSCurve, TopOpeBRepDSHDataStructure};
use super::chfi_ds::{
    ChFiDS_State, ChFiDSMap, ChFiDS_CommonPoint, ChFiDSStripe, SharedSurfData,
};

// =========================================================================
// OCCT ChFi3d_Builder_2.cxx L1629-1737 — ChFi3d_SingularExtremity.
// =========================================================================
pub(crate) fn chfi3d_singular_extremity(
    brep: &topods::BRep,
    stripe: &mut ChFiDSStripe,
    dstr: &mut TopOpeBRepDSHDataStructure,
    vtx: &Shape,
    tol3d: f64,
    tol2d: f64,
) {
    let mut tolreached = 0.0f64;
    let mut pardeb = 0.0f64;
    let mut parfin = 0.0f64;
    let mut c3d: Option<rcad_kernel::geom::Curve3> = None;
    let mut pcurv: Option<rcad_kernel::geom::Curve2d> = None;
    // SurfData and its CommonPoints,
    let mut ivtx = 0i32;
    let mut icurv = 0i32;
    let mut isfirst;

    let periodic = stripe
        .spine()
        .map(|s| s.base().is_periodic())
        .unwrap_or(false);
    let fd_arc;
    if periodic {
        isfirst = true;
        fd_arc = stripe
            .set_of_surf_data()
            .first()
            .expect("SetOfSurfData empty")
            .clone();
    } else {
        let mut sens = 0i32;
        let num = chfi3d_index_of_surf_data(vtx, stripe, &mut sens);
        fd_arc = stripe
            .set_of_surf_data()
            .get((num - 1) as usize)
            .expect("SetOfSurfData value")
            .clone();
        isfirst = sens == 1;
    }

    let (cv1, cv2, von_s1, von_s2, surf) = {
        let fd = fd_arc.read().expect("surfdata lock");
        let cv1 = fd.vertex(isfirst, 1).clone();
        let cv2 = fd.vertex(isfirst, 2).clone();
        // Is it always degenerated ?
        let fi1 = fd.interference_on_s1();
        let fi2 = fd.interference_on_s2();
        let pcs1 = fi1.pcurve_on_surf().expect("PCurveOnSurf");
        let pcs2 = fi2.pcurve_on_surf().expect("PCurveOnSurf");
        let (von_s1, von_s2) = if isfirst {
            (
                pcs1.point_at(fi1.parameter_first()),
                pcs2.point_at(fi2.parameter_first()),
            )
        } else {
            (
                pcs1.point_at(fi1.parameter_last()),
                pcs2.point_at(fi2.parameter_last()),
            )
        };
        let surf = dstr.surface(fd.surf()).surface.clone();
        (cv1, cv2, von_s1, von_s2, surf)
    };
    // Is it always degenerated ?
    if cv1.point.distance(cv2.point) <= 0.0 {
        ivtx = chfi3d_index_point_in_ds(&cv1, dstr);

        let (c3d_out, pcurv_out, pardeb_out, parfin_out, tolreached_out) = chfi3d_compute_arete(
            brep,
            &cv1,
            von_s1,
            &cv2,
            von_s2,
            &surf,
            tol3d,
            tol2d,
            0,
        );
        c3d = c3d_out;
        pcurv = Some(pcurv_out);
        pardeb = pardeb_out;
        parfin = parfin_out;
        tolreached = tolreached_out;
        let crv = TopOpeBRepDSCurve::new(c3d.clone(), tolreached);
        icurv = dstr.add_curve(crv);

        stripe.set_curve(icurv, isfirst);
        stripe.set_parameters(isfirst, pardeb, parfin);
        if let Some(pc) = pcurv.clone() {
            stripe.change_pcurve(isfirst, pc);
        }
        stripe.set_index_point(ivtx, isfirst, 1);
        stripe.set_index_point(ivtx, isfirst, 2);

        let periodic = stripe
            .spine()
            .map(|s| s.base().is_periodic())
            .unwrap_or(false);
        if periodic {
            // periodic case : The operation is renewed
            // the curve 3d is not shared.
            // 2 degenerated edges coinciding in 3d
            isfirst = false;
            let fd_arc_last = stripe
                .set_of_surf_data()
                .last()
                .expect("SetOfSurfData empty")
                .clone();
            let (von_s1, von_s2, surf) = {
                let fd = fd_arc_last.read().expect("surfdata lock");
                let fi1 = fd.interference_on_s1();
                let fi2 = fd.interference_on_s2();
                let pcs1 = fi1.pcurve_on_surf().expect("PCurveOnSurf");
                let pcs2 = fi2.pcurve_on_surf().expect("PCurveOnSurf");
                let von_s1 = pcs1.point_at(fi1.parameter_last());
                let von_s2 = pcs2.point_at(fi2.parameter_last());
                let surf = dstr.surface(fd.surf()).surface.clone();
                (von_s1, von_s2, surf)
            };

            let (c3d_out, pcurv_out, pardeb_out, parfin_out, tolreached_out) = chfi3d_compute_arete(
                brep,
                &cv1,
                von_s1,
                &cv2,
                von_s2,
                &surf,
                tol3d,
                tol2d,
                0,
            );
            c3d = c3d_out;
            pcurv = Some(pcurv_out);
            pardeb = pardeb_out;
            parfin = parfin_out;
            tolreached = tolreached_out;
            let crv = TopOpeBRepDSCurve::new(c3d.clone(), tolreached);
            icurv = dstr.add_curve(crv);

            stripe.set_curve(icurv, isfirst);
            stripe.set_parameters(isfirst, pardeb, parfin);
            if let Some(pc) = pcurv.clone() {
                stripe.change_pcurve(isfirst, pc);
            }
            stripe.set_index_point(ivtx, isfirst, 1);
            stripe.set_index_point(ivtx, isfirst, 2);
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_2.cxx L1739-1764 — IsFree.
// =========================================================================
pub(crate) fn is_free(e: &Shape, ef_map: &ChFiDSMap) -> bool {
    if !ef_map.contains(e) {
        return false;
    }
    let mut fref = Shape::null();
    for it in ef_map.find(e).clone() {
        if fref.is_null() {
            fref = it;
        } else if !fref.is_same(&it) {
            return false;
        }
    }
    true
}

// =========================================================================
// OCCT ChFi3d_Builder_2.cxx L1766-2043 — ChFi3d_MakeExtremities.
// =========================================================================
pub(crate) fn chfi3d_make_extremities(
    brep: &topods::BRep,
    stripe: &mut ChFiDSStripe,
    dstr: &mut TopOpeBRepDSHDataStructure,
    ef_map: &ChFiDSMap,
    tol3d: f64,
    tol2d: f64,
) {
    // OCCT: occ::handle<ChFiDS_Spine>& sp = Stripe->ChangeSpine().
    let mut pardeb = 0.0f64;
    let mut parfin = 0.0f64;
    let mut c3d: Option<rcad_kernel::geom::Curve3> = None;
    let mut tolreached = 0.0f64;
    let periodic = stripe
        .spine()
        .map(|s| s.base().is_periodic())
        .unwrap_or(false);
    if periodic {
        let mut b1 = BndBox::default();
        let mut b2 = BndBox::default();
        let sdf_arc = stripe
            .set_of_surf_data()
            .first()
            .expect("SetOfSurfData empty")
            .clone();
        let (cv1, cv2, uvf1, uvf2, surf) = {
            let sdf = sdf_arc.read().expect("surfdata lock");
            let cv1 = sdf.vertex_first_on_s1().clone();
            let cv2 = sdf.vertex_first_on_s2().clone();
            let fi1 = sdf.interference_on_s1();
            let fi2 = sdf.interference_on_s2();
            let uvf1 = fi1
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(fi1.parameter_first());
            let uvf2 = fi2
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(fi2.parameter_first());
            let surf = dstr.surface(sdf.surf()).surface.clone();
            (cv1, cv2, uvf1, uvf2, surf)
        };
        if cv1.point.distance(cv2.point) > 0.0 {
            let (c3d_out, pcurv_out, pardeb_out, parfin_out, tolreached_out) = chfi3d_compute_arete(
                brep, &cv1, uvf1, &cv2, uvf2, &surf, tol3d, tol2d, 0,
            );
            c3d = c3d_out;
            stripe.pcrv1 = Some(pcurv_out);
            pardeb = pardeb_out;
            parfin = parfin_out;
            tolreached = tolreached_out;
            let crv = TopOpeBRepDSCurve::new(c3d.clone(), tolreached);
            stripe.change_first_curve(dstr.add_curve(crv));
            stripe.change_first_parameters(pardeb, parfin);
            let (ip1, ip2) = {
                let sdf = sdf_arc.read().expect("surfdata lock");
                (
                    chfi3d_index_point_in_ds(sdf.vertex_first_on_s1(), dstr),
                    chfi3d_index_point_in_ds(sdf.vertex_first_on_s2(), dstr),
                )
            };
            stripe.change_index_first_point_on_s1(ip1);
            stripe.change_index_first_point_on_s2(ip2);
            let icurv = stripe.first_curve();
            stripe.change_last_parameters(pardeb, parfin);
            stripe.change_last_curve(icurv);
            stripe.change_index_last_point_on_s1(stripe.index_first_point_on_s1());
            stripe.change_index_last_point_on_s2(stripe.index_first_point_on_s2());

            let sdl_arc = stripe
                .set_of_surf_data()
                .last()
                .expect("SetOfSurfData empty")
                .clone();

            // OCCT L1811-1820: ChFi3d_ComputePCurv + the DS curve tolerance.
            {
                let sdl = sdl_arc.read().expect("surfdata lock");
                let fi1 = sdl.interference_on_s1();
                let fi2 = sdl.interference_on_s2();
                let uv1 = fi1
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(fi1.parameter_last());
                let uv2 = fi2
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(fi2.parameter_last());
                let surf = dstr.surface(sdl.surf()).surface.clone();
                let mut pc = chfi3d_compute_pcurv_2pt(uv1, uv2, pardeb, parfin, false);
                let mut tr = tol3d;
                if let Some(c3d) = &c3d {
                    chfi3d_same_parameter(c3d, &mut pc, &surf, tol3d, &mut tr);
                }
                stripe.pcrv2 = Some(pc);
                let oldtol = dstr.change_curve(icurv).tolerance();
                dstr.change_curve(icurv).set_tolerance(oldtol.max(tr));
                drop(sdl);
            }
            let (cv1v, cv2v) = {
                let sdf = sdf_arc.read().expect("surfdata lock");
                (
                    sdf.vertex_first_on_s1().clone(),
                    sdf.vertex_first_on_s2().clone(),
                )
            };
            if cv1v.is_on_arc() {
                chfi3d_enlarge_box_edge_faces(
                    brep,
                    cv1v.arc(),
                    ef_map.find(cv1v.arc()),
                    cv1v.parameter_on_arc(),
                    &mut b1,
                );
            }

            if cv2v.is_on_arc() {
                chfi3d_enlarge_box_edge_faces(
                    brep,
                    cv2v.arc(),
                    ef_map.find(cv2v.arc()),
                    cv2v.parameter_on_arc(),
                    &mut b2,
                );
            }
            chfi3d_enlarge_box_dstr(
                brep,
                dstr,
                Some(stripe),
                &sdf_arc.read().expect("surfdata lock"),
                &mut b1,
                &mut b2,
                true,
            );
            chfi3d_enlarge_box_dstr(
                brep,
                dstr,
                Some(stripe),
                &sdl_arc.read().expect("surfdata lock"),
                &mut b1,
                &mut b2,
                false,
            );
            if !cv1v.is_vertex() {
                chfi3d_set_point_tolerance(dstr, &b1, stripe.index_first_point_on_s1());
            }
            if !cv2v.is_vertex() {
                chfi3d_set_point_tolerance(dstr, &b2, stripe.index_first_point_on_s2());
            }
        } else {
            // Case of the single extremity
            let cv1_is_vertex = {
                let sdf = sdf_arc.read().expect("surfdata lock");
                sdf.vertex_first_on_s1().is_vertex()
            };
            if cv1_is_vertex {
                let vtx = sdf_arc
                    .read()
                    .expect("surfdata lock")
                    .vertex_first_on_s1()
                    .vertex()
                    .clone();
                chfi3d_singular_extremity(brep, stripe, dstr, &vtx, tol3d, tol2d);
            }
            // OCCT CHFI3D_DEB tracing is a compile-time debug aid — dropped.
        }
        return;
    }

    let sddeb_arc = stripe
        .set_of_surf_data()
        .first()
        .expect("SetOfSurfData empty")
        .clone();

    let (cpdeb1, cpdeb2, freedeb0) = {
        let sddeb = sddeb_arc.read().expect("surfdata lock");
        let cpdeb1 = sddeb.vertex_first_on_s1().clone();
        let cpdeb2 = sddeb.vertex_first_on_s2().clone();
        let freedeb = stripe
            .spine()
            .map(|s| s.base().first_status() == ChFiDS_State::FreeBoundary)
            .unwrap_or(false);
        (cpdeb1, cpdeb2, freedeb)
    };
    let mut freedeb = freedeb0;
    if !freedeb && cpdeb1.is_on_arc() && cpdeb2.is_on_arc() {
        freedeb = is_free(cpdeb1.arc(), ef_map) && is_free(cpdeb2.arc(), ef_map);
    }
    if freedeb {
        // OCCT: sp->SetFirstStatus(ChFiDS_FreeBoundary) over ChangeSpine().
        if let Some(sp) = stripe.my_spine.as_mut() {
            sp.base_mut().set_first_status(ChFiDS_State::FreeBoundary);
        }
        let mut b1 = BndBox::default();
        let mut b2 = BndBox::default();
        if cpdeb1.point.distance(cpdeb2.point) > 0.0 {
            let mut plane = false;
            let uv1;
            let uv2;
            {
                let sddeb = sddeb_arc.read().expect("surfdata lock");
                let fi1 = sddeb.interference_on_s1();
                let fi2 = sddeb.interference_on_s2();
                uv1 = fi1
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(fi1.parameter_first());
                uv2 = fi2
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(fi2.parameter_first());
            }
            // The intersection of the fillet by a plane is attempted

            let hconge = chfi3d_bound_surf(dstr, &sddeb_arc.read().expect("surfdata lock"), 1, 2);
            let mut pcrv1 = stripe.pcrv1.clone();
            chfi3d_coupe_par_plan(
                brep,
                &cpdeb1,
                &cpdeb2,
                &hconge,
                uv1,
                uv2,
                tol3d,
                tol2d,
                &mut c3d,
                &mut pcrv1,
                &mut tolreached,
                &mut pardeb,
                &mut parfin,
                &mut plane,
            );
            stripe.pcrv1 = pcrv1;
            if !plane {
                let (uvf1, uvf2, surf) = {
                    let sddeb = sddeb_arc.read().expect("surfdata lock");
                    let fi1 = sddeb.interference_on_s1();
                    let fi2 = sddeb.interference_on_s2();
                    (
                        fi1.pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(fi1.parameter_first()),
                        fi2.pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(fi2.parameter_first()),
                        dstr.surface(sddeb.surf()).surface.clone(),
                    )
                };
                let (c3d_out, pcurv_out, pardeb_out, parfin_out, tolreached_out) =
                    chfi3d_compute_arete(
                        brep, &cpdeb1, uvf1, &cpdeb2, uvf2, &surf, tol3d, tol2d, 0,
                    );
                c3d = c3d_out;
                stripe.pcrv1 = Some(pcurv_out);
                pardeb = pardeb_out;
                parfin = parfin_out;
                tolreached = tolreached_out;
            }
            let crv = TopOpeBRepDSCurve::new(c3d.clone(), tolreached);
            stripe.change_first_curve(dstr.add_curve(crv));
            stripe.change_first_parameters(pardeb, parfin);
            let (ip1, ip2) = {
                let sddeb = sddeb_arc.read().expect("surfdata lock");
                (
                    chfi3d_index_point_in_ds(sddeb.vertex_first_on_s1(), dstr),
                    chfi3d_index_point_in_ds(sddeb.vertex_first_on_s2(), dstr),
                )
            };
            stripe.change_index_first_point_on_s1(ip1);
            stripe.change_index_first_point_on_s2(ip2);
            if cpdeb1.is_on_arc() {
                chfi3d_enlarge_box_edge_faces(
                    brep,
                    cpdeb1.arc(),
                    ef_map.find(cpdeb1.arc()),
                    cpdeb1.parameter_on_arc(),
                    &mut b1,
                );
            }
            if cpdeb2.is_on_arc() {
                chfi3d_enlarge_box_edge_faces(
                    brep,
                    cpdeb2.arc(),
                    ef_map.find(cpdeb2.arc()),
                    cpdeb2.parameter_on_arc(),
                    &mut b2,
                );
            }
            chfi3d_enlarge_box_dstr(
                brep,
                dstr,
                Some(stripe),
                &sddeb_arc.read().expect("surfdata lock"),
                &mut b1,
                &mut b2,
                true,
            );
            if !cpdeb1.is_vertex() {
                chfi3d_set_point_tolerance(dstr, &b1, stripe.index_first_point_on_s1());
            }
            if !cpdeb2.is_vertex() {
                chfi3d_set_point_tolerance(dstr, &b2, stripe.index_first_point_on_s2());
            }
        } else {
            // Case of a singular extremity
            if cpdeb1.is_vertex() {
                let vtx = cpdeb1.vertex().clone();
                chfi3d_singular_extremity(brep, stripe, dstr, &vtx, tol3d, tol2d);
            }
            // OCCT CHFI3D_DEB tracing is a compile-time debug aid — dropped.
        }
    }
    let sdfin_arc = stripe
        .set_of_surf_data()
        .last()
        .expect("SetOfSurfData empty")
        .clone();
    let (cpfin1, cpfin2, freefin0) = {
        let sdfin = sdfin_arc.read().expect("surfdata lock");
        let cpfin1 = sdfin.vertex_last_on_s1().clone();
        let cpfin2 = sdfin.vertex_last_on_s2().clone();
        let freefin = stripe
            .spine()
            .map(|s| s.base().last_status() == ChFiDS_State::FreeBoundary)
            .unwrap_or(false);
        (cpfin1, cpfin2, freefin)
    };
    let mut freefin = freefin0;
    if !freefin && cpfin1.is_on_arc() && cpfin2.is_on_arc() {
        freefin = is_free(cpfin1.arc(), ef_map) && is_free(cpfin2.arc(), ef_map);
    }
    if freefin {
        // OCCT: sp->SetLastStatus(ChFiDS_FreeBoundary) over ChangeSpine().
        if let Some(sp) = stripe.my_spine.as_mut() {
            sp.base_mut().set_last_status(ChFiDS_State::FreeBoundary);
        }
        let mut b1 = BndBox::default();
        let mut b2 = BndBox::default();
        if cpfin1.point.distance(cpfin2.point) > 0.0 {
            let mut plane = false;
            let uv1;
            let uv2;
            {
                let sdfin = sdfin_arc.read().expect("surfdata lock");
                let fi1 = sdfin.interference_on_s1();
                let fi2 = sdfin.interference_on_s2();
                uv1 = fi1
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(fi1.parameter_last());
                uv2 = fi2
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(fi2.parameter_last());
            }
            // Intersection of the fillet by a plane is attempted

            let hconge = chfi3d_bound_surf(dstr, &sdfin_arc.read().expect("surfdata lock"), 1, 2);
            let mut pcrv2 = stripe.pcrv2.clone();
            chfi3d_coupe_par_plan(
                brep,
                &cpfin1,
                &cpfin2,
                &hconge,
                uv1,
                uv2,
                tol3d,
                tol2d,
                &mut c3d,
                &mut pcrv2,
                &mut tolreached,
                &mut pardeb,
                &mut parfin,
                &mut plane,
            );
            stripe.pcrv2 = pcrv2;
            if !plane {
                let (uvl1, uvl2, surf) = {
                    let sdfin = sdfin_arc.read().expect("surfdata lock");
                    let fi1 = sdfin.interference_on_s1();
                    let fi2 = sdfin.interference_on_s2();
                    (
                        fi1.pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(fi1.parameter_last()),
                        fi2.pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(fi2.parameter_last()),
                        dstr.surface(sdfin.surf()).surface.clone(),
                    )
                };
                let (c3d_out, pcurv_out, pardeb_out, parfin_out, tolreached_out) =
                    chfi3d_compute_arete(
                        brep, &cpfin1, uvl1, &cpfin2, uvl2, &surf, tol3d, tol2d, 0,
                    );
                c3d = c3d_out;
                stripe.pcrv2 = Some(pcurv_out);
                pardeb = pardeb_out;
                parfin = parfin_out;
                tolreached = tolreached_out;
            }
            let crv = TopOpeBRepDSCurve::new(c3d.clone(), tolreached);
            stripe.change_last_curve(dstr.add_curve(crv));
            stripe.change_last_parameters(pardeb, parfin);
            let (ip1, ip2) = {
                let sdfin = sdfin_arc.read().expect("surfdata lock");
                (
                    chfi3d_index_point_in_ds(sdfin.vertex_last_on_s1(), dstr),
                    chfi3d_index_point_in_ds(sdfin.vertex_last_on_s2(), dstr),
                )
            };
            stripe.change_index_last_point_on_s1(ip1);
            stripe.change_index_last_point_on_s2(ip2);
            if cpfin1.is_on_arc() {
                chfi3d_enlarge_box_edge_faces(
                    brep,
                    cpfin1.arc(),
                    ef_map.find(cpfin1.arc()),
                    cpfin1.parameter_on_arc(),
                    &mut b1,
                );
            }
            if cpfin2.is_on_arc() {
                chfi3d_enlarge_box_edge_faces(
                    brep,
                    cpfin2.arc(),
                    ef_map.find(cpfin2.arc()),
                    cpfin2.parameter_on_arc(),
                    &mut b2,
                );
            }
            chfi3d_enlarge_box_dstr(
                brep,
                dstr,
                Some(stripe),
                &sdfin_arc.read().expect("surfdata lock"),
                &mut b1,
                &mut b2,
                false,
            );
            if !cpfin1.is_vertex() {
                chfi3d_set_point_tolerance(dstr, &b1, stripe.index_last_point_on_s1());
            }
            if !cpfin2.is_vertex() {
                chfi3d_set_point_tolerance(dstr, &b2, stripe.index_last_point_on_s2());
            }
        } else {
            // Case of the single extremity
            if cpfin1.is_vertex() {
                let vtx = cpfin1.vertex().clone();
                chfi3d_singular_extremity(brep, stripe, dstr, &vtx, tol3d, tol2d);
            }
            // OCCT CHFI3D_DEB tracing is a compile-time debug aid — dropped.
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_2.cxx L2045-2100 — ChFi3d_Purge.
// =========================================================================
pub(crate) fn chfi3d_purge(
    stripe: &mut ChFiDSStripe,
    sd: &SharedSurfData,
    v_ref: &ChFiDS_CommonPoint,
    isfirst: bool,
    ons: i32,
    intf: &mut i32,
    intl: &mut i32,
) {
    if isfirst {
        *intf = 1;
    } else {
        *intl = 1; // End.
    }
    let opp = 3 - ons;
    {
        let sdg = sd.read().expect("surfdata lock");
        if !sdg.vertex(isfirst, opp).is_on_arc() || sdg.twist_on_s1() || sdg.twist_on_s2() {
            let seq = stripe.change_set_of_surf_data();
            if isfirst {
                seq.remove(0);
            } else {
                seq.pop();
            }
            return;
        }
    }
    let mut sdg = sd.write().expect("surfdata lock");
    if ons == 1 {
        sdg.change_index_of_s1(0);
    } else {
        sdg.change_index_of_s2(0);
    }

    *sdg.change_vertex(!isfirst, ons) = v_ref.clone();
    *sdg.change_vertex(isfirst, ons) = v_ref.clone();

    let fi = sdg.change_interference(ons);
    if isfirst {
        fi.set_first_parameter(fi.parameter_last());
    } else {
        fi.set_last_parameter(fi.parameter_first());
    }
    fi.set_line_index(0);
}

// =========================================================================
// OCCT ChFi3d_Builder_2.cxx L2102-2132 — InsertAfter.
// =========================================================================
pub(crate) fn insert_after(
    stripe: &mut ChFiDSStripe,
    ref_: &Option<SharedSurfData>,
    item: &SharedSurfData,
) {
    if let Some(r) = ref_ {
        if std::sync::Arc::ptr_eq(r, item) {
            panic!("Standard_Failure: InsertAfter : twice the same surfdata.");
        }
    }

    let seq = stripe.change_set_of_surf_data();

    if seq.is_empty() || ref_.is_none() {
        seq.insert(0, item.clone());
    }
    for i in 0..seq.len() {
        if ref_
            .as_ref()
            .map_or(false, |r| std::sync::Arc::ptr_eq(&seq[i], r))
        {
            seq.insert(i + 1, item.clone());
            break;
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_2.cxx L2134-2171 — RemoveSD.
// =========================================================================
pub(crate) fn remove_sd(
    stripe: &mut ChFiDSStripe,
    prev: &Option<SharedSurfData>,
    next: &Option<SharedSurfData>,
) {
    let seq = stripe.change_set_of_surf_data();
    if seq.is_empty() {
        return;
    }
    let mut iprev = 0usize;
    let mut inext = 0usize;
    for i in 0..seq.len() {
        if prev
            .as_ref()
            .map_or(false, |p| std::sync::Arc::ptr_eq(&seq[i], p))
        {
            iprev = i + 2; // OCCT 1-based: iprev = i + 1 (then used 1-based).
        }
        if next
            .as_ref()
            .map_or(false, |n| std::sync::Arc::ptr_eq(&seq[i], n))
        {
            inext = i; // OCCT 1-based: inext = i - 1.
            break;
        }
    }
    if prev.is_none() {
        iprev = 1;
    }
    if next.is_none() {
        inext = seq.len(); // OCCT: Seq.Length() (1-based).
    }
    if iprev <= inext {
        // OCCT: Seq.Remove(iprev, inext) — 1-based inclusive; rcad 0-based
        // half-open [iprev-1, inext).
        seq.drain((iprev - 1)..inext);
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_2.cxx L2173-2203 — InsertBefore.
// =========================================================================
pub(crate) fn insert_before(
    stripe: &mut ChFiDSStripe,
    ref_: &Option<SharedSurfData>,
    item: &SharedSurfData,
) {
    if let Some(r) = ref_ {
        if std::sync::Arc::ptr_eq(r, item) {
            panic!("Standard_Failure: InsertBefore : twice the same surfdata.");
        }
    }

    let seq = stripe.change_set_of_surf_data();

    if seq.is_empty() || ref_.is_none() {
        seq.push(item.clone());
    }
    for i in 0..seq.len() {
        if ref_
            .as_ref()
            .map_or(false, |r| std::sync::Arc::ptr_eq(&seq[i], r))
        {
            seq.insert(i, item.clone());
            break;
        }
    }
}
