//! OCCT ChFi3d_Builder_C1.cxx L3771-4529 — PerformMoreSurfdata, 1:1
//! translation (split out of chfi3d_builder_c2.rs for file size).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2dEval as _, Curve3, CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};

use super::chfi3d::{chfi3d_index_of_surf_data, chfi3d_index_point_in_ds, topabs_compose, topabs_reverse, ChFi3dBuilder};
use super::chfi3d_builder_0::{
    chfi3d_compute_arete, chfi3d_fil_curve_in_ds, chfi3d_fil_point_in_ds, chfi3d_project_pcurv,
    GeomAdaptorSurface,
};
use super::chfi3d_builder_c2::{c3d_first, c3d_last, contain_e};
use super::chfi3d_ds::{
    TopOpeBRepDSCurve, TopOpeBRepDSSolidSurfaceInterference, TopOpeBRepDSPoint,
    TopOpeBRepDSSurface,
};

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L3771-4529 — PerformMoreSurfdata.
// Determine intersections at end on several surfdata.
// GAP carriers: the GeomInt_IntSS pcurves (LineOnS1/LineOnS2) and
// TolReached3d are pending (the kernel IntSS returns 3D lines only) — the
// rcad translation stores null pcurves there and keeps the OCCT InterSS
// tolerance; GeomProjLib::Curve2d routes through the chfi3d_project_pcurv
// stand-in.
// =========================================================================
pub fn perform_more_surfdata(this: &mut ChFi3dBuilder, index: usize) {
    let a_tol3d = 1.0e-4;
    // OCCT L3773-3779.
    if this.my_vdata_map.find_from_index(index).is_empty() {
        return;
    }
    let a_stripe = this.my_vdata_map.find_from_index(index)[0].clone();
    let mut st = a_stripe.write().expect("stripe lock");
    let a_spine = st.spine().expect("spine").clone();
    let a_vtx = this.my_vdata_map.find_key(index).clone();
    let brep = this.my_brep.clone();
    let tol2d = this.tol2d;

    let mut a_sens = 0i32;
    let a_ind = chfi3d_index_of_surf_data(&a_vtx, &st, &mut a_sens);
    let is_first = a_sens == 1;
    let a_ind_prev;
    if a_sens == 1 {
        a_ind_prev = a_ind + 1;
    } else {
        a_ind_prev = a_ind - 1;
    }

    // OCCT L3794-3801: aSurfData / aCP1 / aCP2 of anInd.
    let (a_cp1, a_cp2) = {
        let sd = st.my_hdata[(a_ind - 1) as usize].read().expect("surfdata lock");
        (
            sd.vertex(is_first, 1).clone(),
            sd.vertex(is_first, 2).clone(),
        )
    };

    // OCCT L3803-3812: aFace by FindFace; aSurfPrev = BRep_Tool::Surface.
    let mut a_face = Shape::null();
    if !this.find_face(&a_cp1, &a_cp2, &mut a_face) {
        // OCCT: FindFace result is not checked; the surface query would
        // raise on a null face — keep the OCCT raise.
        panic!("Standard_NullObject: PerformMoreSurfdata FindFace");
    }
    let a_surf_prev = brep.face_surface_world(&a_face).expect("aSurfPrev");

    // determination of the two arcs of the face at end (OCCT L3814-3843).
    let ve_list: Vec<Shape> = if this.my_ve_map.contains(&a_vtx) {
        this.my_ve_map.find(&a_vtx).clone()
    } else {
        Vec::new()
    };
    let mut an_arc1 = Shape::null();
    let mut is_found = false;
    for it in &ve_list {
        if is_found {
            break;
        }
        an_arc1 = it.clone();
        if contain_e(&brep, &a_face, &an_arc1) {
            is_found = true;
        }
    }
    let mut an_arc2 = Shape::null();
    is_found = false;
    for it in &ve_list {
        if is_found {
            break;
        }
        an_arc2 = it.clone();
        if contain_e(&brep, &a_face, &an_arc2) && !an_arc2.is_same(&an_arc1) {
            is_found = true;
        }
    }

    // determination of common points aCP1onArc, aCP2onArc, aCP2NotonArc
    // (OCCT L3846-3880).
    let a_surf_data_prev = st.my_hdata[(a_ind_prev - 1) as usize].clone();
    let (b_cp1, b_cp2) = {
        let sd = a_surf_data_prev.read().expect("surfdata lock");
        (
            sd.vertex(is_first, 1).clone(),
            sd.vertex(is_first, 2).clone(),
        )
    };
    let is2nd_cp1_on_arc;
    let a_cp2_on_arc;
    let a_cp2_noton_arc;
    if b_cp1.is_on_arc() && (b_cp1.arc().is_same(&an_arc1) || b_cp1.arc().is_same(&an_arc2)) {
        a_cp2_on_arc = b_cp1.clone();
        a_cp2_noton_arc = b_cp2.clone();
        is2nd_cp1_on_arc = true;
    } else if b_cp2.is_on_arc()
        && (b_cp2.arc().is_same(&an_arc1) || b_cp2.arc().is_same(&an_arc2))
    {
        a_cp2_on_arc = b_cp2.clone();
        a_cp2_noton_arc = b_cp1.clone();
        is2nd_cp1_on_arc = false;
    } else {
        return;
    }

    let is1st_cp1_on_arc;
    let a_cp1_on_arc;
    if a_cp1.point().distance(a_cp2_on_arc.point()) <= a_tol3d {
        a_cp1_on_arc = a_cp2.clone();
        is1st_cp1_on_arc = false;
    } else {
        a_cp1_on_arc = a_cp1.clone();
        is1st_cp1_on_arc = true;
    }
    if !a_cp1_on_arc.is_on_arc() {
        return;
    }

    // determination of neighbor surface (OCCT L3884-3894).
    // OCCT L3886-3893: myListStripe.First()->SetOfSurfData()->Value(anInd).
    // The first list stripe is usually the held a_stripe (same Arc); RwLock
    // is not re-entrant, so the held guard serves that case.
    let list_stripe = &this.my_list_stripe[0];
    let ind_surface = if std::sync::Arc::ptr_eq(list_stripe, &a_stripe) {
        let sd = st.my_hdata[(a_ind - 1) as usize].read().expect("surfdata lock");
        if is1st_cp1_on_arc { sd.index_of_s1 } else { sd.index_of_s2 }
    } else {
        let ls = list_stripe.read().expect("stripe lock");
        let sd = ls.my_hdata[(a_ind - 1) as usize].read().expect("surfdata lock");
        if is1st_cp1_on_arc { sd.index_of_s1 } else { sd.index_of_s2 }
    };
    let a_neighbor_face = this.my_ds.as_ref().expect("DS").shape(ind_surface).clone();

    // calculation of intersections — ChFi3d_ComputeArete (OCCT L3897-3930).
    let a_surf_data = st.my_hdata[(a_ind - 1) as usize].clone();
    let a_surf_lock = a_surf_data.read().expect("surfdata lock");
    let a_surf_of_sd = this
        .my_ds
        .as_ref()
        .expect("DS")
        .surface(a_surf_lock.surf())
        .surface
        .clone();
    let (a_cracc, a_pcurv1, _a_parf, _a_parl, a_tol_reached) = if is_first {
        chfi3d_compute_arete(
            &brep,
            a_surf_lock.vertex_last_on_s1(),
            a_surf_lock
                .interference_on_s1()
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(a_surf_lock.interference_on_s1().parameter_last()),
            a_surf_lock.vertex_last_on_s2(),
            a_surf_lock
                .interference_on_s2()
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(a_surf_lock.interference_on_s2().parameter_last()),
            &a_surf_of_sd,
            a_tol3d,
            tol2d,
            0,
        )
    } else {
        chfi3d_compute_arete(
            &brep,
            a_surf_lock.vertex_first_on_s1(),
            a_surf_lock
                .interference_on_s1()
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(a_surf_lock.interference_on_s1().parameter_first()),
            a_surf_lock.vertex_first_on_s2(),
            a_surf_lock
                .interference_on_s2()
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(a_surf_lock.interference_on_s2().parameter_first()),
            &a_surf_of_sd,
            a_tol3d,
            tol2d,
            0,
        )
    };
    let a_cracc = a_cracc.expect("aCracc");

    // the line of anInd (OCCT L3933-3948).
    let a_fi = if is1st_cp1_on_arc {
        a_surf_lock.interference_on_s1().clone()
    } else {
        a_surf_lock.interference_on_s2().clone()
    };
    let ind_line = a_fi.line_index();
    let a_cline = this
        .my_ds
        .as_ref()
        .expect("DS")
        .curve(a_fi.line_index())
        .curve
        .clone()
        .expect("aCline");
    let a_pcline_on_surf = a_fi.pcurve_on_surf().cloned();
    let a_pcline_on_face = a_fi.pcurve_on_face().cloned();
    drop(a_surf_lock);

    // intersection between the SurfData number anInd and the face aFace
    // (OCCT L3951-3957).
    let a_surf = this
        .my_ds
        .as_ref()
        .expect("DS")
        .surface(st.my_hdata[(a_ind - 1) as usize].read().expect("sd").surf())
        .surface
        .clone();
    // GAP: GeomInt_IntSS pcurves (LineOnS1/LineOnS2) and TolReached3d are
    // pending; the kernel IntSS yields 3D lines and the InterSS tolerance.
    let mut an_inter_ss = rcad_kernel::base::geom_api::int_ss::IntSS::with_surfaces(
        &a_surf_prev,
        &a_surf,
        1.0e-7,
    );
    let mut a_cint1: Option<Curve3> = None;
    let mut a_tolex1 = 0.0f64;
    let mut a_pext1 = DVec3::ZERO;
    let mut a_pext2 = DVec3::ZERO;
    if !an_inter_ss.is_done() {
        return;
    }
    is_found = false;
    for i in 1..=an_inter_ss.nb_lines() {
        if is_found {
            break;
        }
        a_cint1 = Some(an_inter_ss.line(i).clone());
        let c1 = a_cint1.as_ref().expect("aCint1");
        a_tolex1 = 1.0e-7; // OCCT: ChFi3d_EvalTolReached over the InterSS pcurves — GAP.
        a_pext1 = c1.point_at(c3d_first(c1));
        a_pext2 = c1.point_at(c3d_last(c1));
        if a_pext1.distance(a_cp1_on_arc.point()) <= a_tol3d
            || a_pext2.distance(a_cp1_on_arc.point()) <= a_tol3d
        {
            is_found = true;
        }
    }
    if !is_found {
        return;
    }

    // OCCT L4027-4043: aPext selection.
    let mut a_pext = DVec3::ZERO;
    let is_pext_found;
    if a_pext1.distance(a_cp2_on_arc.point()) > a_tol3d
        && a_pext1.distance(a_cp1_on_arc.point()) > a_tol3d
    {
        a_pext = a_pext1;
        is_pext_found = true;
    } else if a_pext2.distance(a_cp2_on_arc.point()) > a_tol3d
        && a_pext2.distance(a_cp1_on_arc.point()) > a_tol3d
    {
        a_pext = a_pext2;
        is_pext_found = true;
    } else {
        is_pext_found = false;
    }

    let mut is_do_second_section = false;
    let mut a_par = 0.0f64;
    if is_pext_found {
        // OCCT L4046-4072: Extrema_ExtPC aPext on aCracc, minimum param.
        let [f0, l0] = a_cracc.default_domain();
        let n = 96usize;
        let mut best = f0;
        let mut bestd = f64::MAX;
        for i in 0..=n {
            let t = f0 + (l0 - f0) * (i as f64) / (n as f64);
            let d = a_cracc.point_at(t).distance(a_pext);
            if d < bestd {
                bestd = d;
                best = t;
            }
        }
        if bestd <= a_tol3d {
            a_par = best;
            is_do_second_section = true;
        }
    }

    // OCCT L4075-4095: orientations and the second surface.
    let an_or_sd1 = {
        st.my_hdata[(a_ind - 1) as usize].read().expect("sd").orientation
    };
    let an_or_sd2 = a_surf_data_prev.read().expect("sd").orientation;
    let a_surf2 = this
        .my_ds
        .as_ref()
        .expect("DS")
        .surface(a_surf_data_prev.read().expect("sd").surf())
        .surface
        .clone();

    // The second section (OCCT L4085-4146).
    let mut a_cint2: Option<Curve3> = None;
    let mut a_tolex2 = 0.0f64;
    let a_tr_cracc: Curve3;
    if is_do_second_section {
        let a_par1 = if a_cracc
            .point_at(c3d_first(&a_cracc))
            .distance(a_cp2_noton_arc.point())
            <= a_tol3d
        {
            c3d_first(&a_cracc)
        } else {
            c3d_last(&a_cracc)
        };
        let (lo, hi) = if a_par1 < a_par {
            (a_par1, a_par)
        } else {
            (a_par, a_par1)
        };
        a_tr_cracc = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
            a_cracc.clone(),
            lo,
            hi,
        ));

        let mut an_inter_ss2 = rcad_kernel::base::geom_api::int_ss::IntSS::with_surfaces(
            &a_surf_prev,
            &a_surf2,
            1.0e-7,
        );
        if !an_inter_ss2.is_done() {
            return;
        }
        is_found = false;
        for i in 1..=an_inter_ss2.nb_lines() {
            if is_found {
                break;
            }
            a_cint2 = Some(an_inter_ss2.line(i).clone());
            let c2 = a_cint2.as_ref().expect("aCint2");
            a_tolex2 = 1.0e-7; // GAP: TolReached2 via the InterSS pcurves.
            a_pext1 = c2.point_at(c3d_first(c2));
            a_pext2 = c2.point_at(c3d_last(c2));
            if a_pext1.distance(a_cp2_on_arc.point()) <= a_tol3d
                || a_pext2.distance(a_cp2_on_arc.point()) <= a_tol3d
            {
                is_found = true;
            }
        }
        if !is_found {
            return;
        }
    } else {
        a_tr_cracc = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
            a_cracc.clone(),
            c3d_first(&a_cracc),
            c3d_last(&a_cracc),
        ));
    }

    // Storage of the data structure (OCCT L4150-4186).
    // GAP: GeomProjLib::Curve2d — the chfi3d_project_pcurv stand-in.
    let a_pcracc_s: Option<rcad_kernel::geom::Curve2d> = chfi3d_project_pcurv(
        &a_tr_cracc,
        &GeomAdaptorSurface::new(a_surf2.clone()),
        a_tol3d,
    )
    .map(|(c, _)| c);

    let (a_fi_prev, ind_shape) = {
        let sd = a_surf_data_prev.read().expect("sd");
        if is2nd_cp1_on_arc {
            (sd.interference_on_s1().clone(), sd.index_of_s1)
        } else {
            (sd.interference_on_s2().clone(), sd.index_of_s2)
        }
    };
    if ind_shape <= 0 {
        return;
    }

    let mut a_cur_orient = this.my_ds.as_ref().expect("DS").shape(ind_shape).orientation;
    a_cur_orient = topabs_compose(a_cur_orient, an_or_sd2);
    a_cur_orient = topabs_compose(topabs_reverse(a_fi_prev.transition()), a_cur_orient);

    // Filling the data structure (OCCT L4190-4222).
    let a_pt_cp1 = TopOpeBRepDSPoint::new(a_cp1_on_arc.point(), a_cp1_on_arc.tolerance());
    let ind_cp1on_arc = this.my_ds.as_mut().expect("DS").add_point(a_pt_cp1);
    let ind_surf1 = {
        st.my_hdata[(a_ind - 1) as usize].read().expect("sd").surf()
    };
    let ind_arc1 = this
        .my_ds
        .as_mut()
        .expect("DS")
        .add_shape(a_cp1_on_arc.arc());
    let ind_sol = st.solid_index();

    let an_interfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
        a_cp1_on_arc.transition_on_arc(),
        ind_arc1,
        ind_cp1on_arc,
        a_cp1_on_arc.parameter_on_arc(),
        false,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_shape_interferences(ind_arc1)
        .push(an_interfp1);

    // OCCT L4199-4205: the SolidSurfaceInterference on the solid.
    let ssi = TopOpeBRepDSSolidSurfaceInterference::new(
        an_or_sd1,
        super::chfi3d_ds::TopOpeBRepDSKind::Solid,
        ind_sol,
        super::chfi3d_ds::TopOpeBRepDSKind::Surface,
        ind_surf1,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_shape_interferences(ind_sol)
        .push(super::chfi3d_ds::TopOpeBRepDSInterference::SolidSurface(ssi));

    // deletion of Surface Data (OCCT L4208-4212).
    st.my_hdata.remove((a_ind - 1) as usize);
    let mut an_ind = a_ind;
    if !is_first {
        an_ind -= 1;
    }
    let a_surf_data = st.my_hdata[(an_ind - 1) as usize].clone();

    // definition of indices of common points in Data Structure
    // (OCCT L4215-4243).
    let ind_cp2on_arc;
    let ind_cp2_noton_arc;
    if is2nd_cp1_on_arc {
        st.set_index_point(
            chfi3d_index_point_in_ds(&a_cp2_on_arc, this.my_ds.as_mut().expect("DS")),
            is_first,
            1,
        );
        st.set_index_point(
            chfi3d_index_point_in_ds(&a_cp2_noton_arc, this.my_ds.as_mut().expect("DS")),
            is_first,
            2,
        );
        if is_first {
            ind_cp2on_arc = st.index_point(is_first, 1);
            ind_cp2_noton_arc = st.index_point(is_first, 2);
        } else {
            ind_cp2on_arc = st.index_point(is_first, 1);
            ind_cp2_noton_arc = st.index_point(is_first, 2);
        }
    } else {
        st.set_index_point(
            chfi3d_index_point_in_ds(&a_cp2_on_arc, this.my_ds.as_mut().expect("DS")),
            is_first,
            2,
        );
        st.set_index_point(
            chfi3d_index_point_in_ds(&a_cp2_noton_arc, this.my_ds.as_mut().expect("DS")),
            is_first,
            1,
        );
        if is_first {
            ind_cp2on_arc = st.index_point(is_first, 2);
            ind_cp2_noton_arc = st.index_point(is_first, 1);
        } else {
            ind_cp2on_arc = st.index_point(is_first, 2);
            ind_cp2_noton_arc = st.index_point(is_first, 1);
        }
    }

    let (a_fi, ind_shape) = {
        let sd = a_surf_data.read().expect("sd");
        if is2nd_cp1_on_arc {
            (sd.interference_on_s1().clone(), sd.index_of_s1)
        } else {
            (sd.interference_on_s2().clone(), sd.index_of_s2)
        }
    };

    // OCCT L4251-4263: aPoint1/aPoint2 evaluations.
    let stemp2 = brep
        .face_surface_world(&this.my_ds.as_ref().expect("DS").shape(ind_shape).clone())
        .expect("Stemp2");
    let a_p2d_last = a_fi
        .pcurve_on_face()
        .expect("aFI PCurveOnFace")
        .point_at(a_fi.parameter_last());
    let a_point2 = stemp2.point_at(a_p2d_last.x, a_p2d_last.y);
    let a_p2d_first = a_fi
        .pcurve_on_face()
        .expect("aFI PCurveOnFace")
        .point_at(a_fi.parameter_first());
    let a_point1 = stemp2.point_at(a_p2d_first.x, a_p2d_first.y);

    let mut an_or_surf = a_cur_orient;
    let mut an_or_face = a_face.orientation;
    let inda_face = this.my_ds.as_mut().expect("DS").add_shape(&a_face);
    let mut ind_point = ind_cp2on_arc;
    let mut ind_point1;
    let mut ind_point2;

    if is_do_second_section {
        // OCCT L4266-4316: storage of aCint2.
        let tpoint = TopOpeBRepDSPoint::new(a_pext, a_tolex2);
        ind_point = this.my_ds.as_mut().expect("DS").add_point(tpoint);
        let ind_curve = this
            .my_ds
            .as_mut()
            .expect("DS")
            .add_curve(TopOpeBRepDSCurve::new(a_cint2.clone(), a_tolex2));
        let c2 = a_cint2.as_ref().expect("aCint2");
        a_pext1 = c2.point_at(c3d_first(c2));
        a_pext2 = c2.point_at(c3d_last(c2));
        if a_pext1.distance(a_pext) <= a_tol3d {
            ind_point1 = ind_point;
            ind_point2 = ind_cp2on_arc;
        } else {
            ind_point1 = ind_cp2on_arc;
            ind_point2 = ind_point;
        }
        // define the orientation of aCint2.
        if a_pext1.distance(a_point2) > a_tol3d && a_pext2.distance(a_point1) > a_tol3d {
            an_or_surf = topabs_reverse(an_or_surf);
        }
        let itfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Forward,
            ind_curve,
            ind_point1,
            c3d_first(c2),
            false,
        );
        let itfp2 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Reversed,
            ind_curve,
            ind_point2,
            c3d_last(c2),
            false,
        );
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_curve_interferences(ind_curve)
            .push(itfp1);
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_curve_interferences(ind_curve)
            .push(itfp2);

        // interference of aCint2 on the SurfData number anIndPrev.
        let prev_surf = this
            .my_ds
            .as_ref()
            .expect("DS")
            .surface(a_surf_data.read().expect("sd").surf())
            .surface
            .clone();
        // GAP: aPCint22 (the InterSS pcurve) is pending — stored null.
        let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
            ind_curve,
            a_surf_data.read().expect("sd").surf(),
            None,
            an_or_surf,
        );
        let _ = prev_surf;
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_surface_interferences(a_surf_data.read().expect("sd").surf())
            .push(itfc);

        // interference of aCint2 on aFace (OCCT L4293-4307).
        if an_or_face == an_or_sd2 {
            an_or_face = topabs_reverse(an_or_surf);
        } else {
            an_or_face = an_or_surf;
        }
        // GAP: aPCint21 pending — stored null.
        let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
            ind_curve,
            inda_face,
            None,
            an_or_face,
        );
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_shape_interferences(inda_face)
            .push(itfc);
    }

    // storage of aTrCracc (OCCT L4319-4390).
    a_pext1 = a_tr_cracc.point_at(c3d_first(&a_tr_cracc));
    a_pext2 = a_tr_cracc.point_at(c3d_last(&a_tr_cracc));
    if a_pext1.distance(a_cp2_noton_arc.point()) <= a_tol3d {
        ind_point1 = ind_cp2_noton_arc;
        ind_point2 = ind_point;
    } else {
        ind_point1 = ind_point;
        ind_point2 = ind_cp2_noton_arc;
    }

    // Define the orientation of aTrCracc (OCCT L4337-4359).
    let is_to_reverse = if is_do_second_section {
        let a_p1 = a_tr_cracc.point_at(c3d_first(&a_tr_cracc));
        let a_p2 = a_tr_cracc.point_at(c3d_last(&a_tr_cracc));
        let a_p3 = a_cint2
            .as_ref()
            .expect("aCint2")
            .point_at(c3d_first(a_cint2.as_ref().expect("aCint2")));
        let a_p4 = a_cint2
            .as_ref()
            .expect("aCint2")
            .point_at(c3d_last(a_cint2.as_ref().expect("aCint2")));
        a_p1.distance(a_p4) > a_tol3d && a_p2.distance(a_p3) > a_tol3d
    } else {
        a_pext1.distance(a_point2) > a_tol3d && a_pext2.distance(a_point1) > a_tol3d
    };
    if is_to_reverse {
        an_or_surf = topabs_reverse(an_or_surf);
    }

    let ind_curve = this
        .my_ds
        .as_mut()
        .expect("DS")
        .add_curve(TopOpeBRepDSCurve::new(Some(a_tr_cracc.clone()), a_tol_reached));
    let itfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
        Orientation::Forward,
        ind_curve,
        ind_point1,
        c3d_first(&a_tr_cracc),
        false,
    );
    let itfp2 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
        Orientation::Reversed,
        ind_curve,
        ind_point2,
        c3d_last(&a_tr_cracc),
        false,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_curve_interferences(ind_curve)
        .push(itfp1);
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_curve_interferences(ind_curve)
        .push(itfp2);

    // interference of aTrCracc on the SurfData number anIndPrev
    // (OCCT L4372-4378).
    let prev_isurf = a_surf_data.read().expect("sd").surf();
    let itfc =
        super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(ind_curve, prev_isurf, a_pcracc_s.clone(), an_or_surf);
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_surface_interferences(prev_isurf)
        .push(itfc);
    st.in_ds(is_first, 1);

    // interference of aTrCracc on the SurfData number anInd
    // (OCCT L4381-4388).
    if an_or_sd1 == an_or_sd2 {
        an_or_surf = topabs_reverse(an_or_surf);
    }
    let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
        ind_curve,
        ind_surf1,
        Some(a_pcurv1.clone()),
        an_or_surf,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_surface_interferences(ind_surf1)
        .push(itfc);

    // storage of aCint1 (OCCT L4391-4440).
    a_pext1 = a_cint1
        .as_ref()
        .expect("aCint1")
        .point_at(c3d_first(a_cint1.as_ref().expect("aCint1")));
    if a_pext1.distance(a_cp1_on_arc.point()) <= a_tol3d {
        ind_point1 = ind_cp1on_arc;
        ind_point2 = ind_point;
    } else {
        ind_point1 = ind_point;
        ind_point2 = ind_cp1on_arc;
    }

    // definition of the orientation of aCint1.
    let a_p1 = a_cint1
        .as_ref()
        .expect("aCint1")
        .point_at(c3d_first(a_cint1.as_ref().expect("aCint1")));
    let a_p2 = a_cint1
        .as_ref()
        .expect("aCint1")
        .point_at(c3d_last(a_cint1.as_ref().expect("aCint1")));
    let a_p3 = a_tr_cracc.point_at(c3d_first(&a_tr_cracc));
    let a_p4 = a_tr_cracc.point_at(c3d_last(&a_tr_cracc));
    if a_p1.distance(a_p4) > a_tol3d && a_p2.distance(a_p3) > a_tol3d {
        an_or_surf = topabs_reverse(an_or_surf);
    }

    let ind_curve = this
        .my_ds
        .as_mut()
        .expect("DS")
        .add_curve(TopOpeBRepDSCurve::new(a_cint1.clone(), a_tolex1));
    let itfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
        Orientation::Forward,
        ind_curve,
        ind_point1,
        c3d_first(a_cint1.as_ref().expect("aCint1")),
        false,
    );
    let itfp2 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
        Orientation::Reversed,
        ind_curve,
        ind_point2,
        c3d_last(a_cint1.as_ref().expect("aCint1")),
        false,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_curve_interferences(ind_curve)
        .push(itfp1);
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_curve_interferences(ind_curve)
        .push(itfp2);

    // interference of aCint1 on the SurfData number anInd (OCCT L4420).
    // GAP: aPCint12 pending — stored null.
    let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
        ind_curve,
        ind_surf1,
        None,
        an_or_surf,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_surface_interferences(ind_surf1)
        .push(itfc);

    // interference of aCint1 on aFace (OCCT L4424-4435).
    an_or_face = a_face.orientation;
    if an_or_face == an_or_sd1 {
        an_or_face = topabs_reverse(an_or_surf);
    } else {
        an_or_face = an_or_surf;
    }
    // GAP: aPCint11 pending — stored null.
    let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
        ind_curve,
        inda_face,
        None,
        an_or_face,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_shape_interferences(inda_face)
        .push(itfc);

    // storage of aCline passing through aCP1onArc and aCP2NotonArc
    // (OCCT L4443-4500).
    let a_tr_cline = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
        a_cline.clone(),
        c3d_first(&a_cline),
        c3d_last(&a_cline),
    ));
    let a_tolerance = this.my_ds.as_ref().expect("DS").curve(ind_line).tolerance();
    let ind_curve = this
        .my_ds
        .as_mut()
        .expect("DS")
        .add_curve(TopOpeBRepDSCurve::new(Some(a_tr_cline.clone()), a_tolerance));

    a_pext1 = a_tr_cline.point_at(c3d_first(&a_tr_cline));
    if a_pext1.distance(a_cp1_on_arc.point()) < a_tol3d {
        ind_point1 = ind_cp1on_arc;
        ind_point2 = ind_cp2_noton_arc;
    } else {
        ind_point1 = ind_cp2_noton_arc;
        ind_point2 = ind_cp1on_arc;
    }
    // definition of the orientation of aTrCline.
    let a_p1 = a_tr_cline.point_at(c3d_first(&a_tr_cline));
    let a_p2 = a_tr_cline.point_at(c3d_last(&a_tr_cline));
    let a_p3 = a_cint1
        .as_ref()
        .expect("aCint1")
        .point_at(c3d_first(a_cint1.as_ref().expect("aCint1")));
    let a_p4 = a_cint1
        .as_ref()
        .expect("aCint1")
        .point_at(c3d_last(a_cint1.as_ref().expect("aCint1")));
    if a_p1.distance(a_p4) > a_tol3d && a_p2.distance(a_p3) > a_tol3d {
        an_or_surf = topabs_reverse(an_or_surf);
    }

    let itfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
        Orientation::Forward,
        ind_curve,
        ind_point1,
        c3d_first(&a_tr_cline),
        false,
    );
    let itfp2 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
        Orientation::Reversed,
        ind_curve,
        ind_point2,
        c3d_last(&a_tr_cline),
        false,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_curve_interferences(ind_curve)
        .push(itfp1);
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_curve_interferences(ind_curve)
        .push(itfp2);

    // interference of aTrCline on the SurfData number anInd
    // (OCCT L4473-4478).
    let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
        ind_curve,
        ind_surf1,
        a_pcline_on_surf.clone(),
        an_or_surf,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_surface_interferences(ind_surf1)
        .push(itfc);

    // interference de ctlin par rapport a Fvoisin (OCCT L4481-4497).
    let ind_shape = this.my_ds.as_mut().expect("DS").add_shape(&a_neighbor_face);
    an_or_face = a_neighbor_face.orientation;
    if an_or_face == an_or_sd1 {
        an_or_face = topabs_reverse(an_or_surf);
    } else {
        an_or_face = an_or_surf;
    }
    let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
        ind_curve,
        ind_shape,
        a_pcline_on_face.clone(),
        an_or_face,
    );
    this.my_ds
        .as_mut()
        .expect("DS")
        .change_shape_interferences(ind_shape)
        .push(itfc);
}
