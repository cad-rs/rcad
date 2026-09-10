//! OCCT ChFi3d_Builder_CnCrn.cxx L1237-3927 — `PerformMoreThreeCorner`
//! (Stage 1f).  Process case of a top with n edges.
//!
//! The helper statics and the OCCT-form carriers live in
//! `chfi3d_builder_cncrn.rs`.  Borrow structure: OCCT holds one
//! `DStr = myDS->ChangeDS()` for the whole body; rcad splits the body into
//! two `dstr` scopes so the mid-function `PerformOneCorner` call (OCCT
//! L1676) can take `&mut self`.

use std::sync::{Arc, RwLock};

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::BRepTool as _;
use rcad_kernel::topods::{self, Orientation, Shape};

use super::chfi3d::ChFi3dBuilder;
use super::chfi3d::{chfi3d_index_of_surf_data, is_tangent_faces, topabs_reverse};
use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_cherche_bords_libres, chfi3d_compute_arete,
    chfi3d_couture_on_vertex, chfi3d_fil_curve_in_ds, chfi3d_fil_point_in_ds,
    chfi3d_fil_vertex_in_ds, chfi3d_nb_not_degenerated_edges, topexp_face_edges, GeomAdaptorSurface,
};
use super::chfi3d_builder_cncrn::{
    calcul_batten, calcul_c2d_on_face, calcul_droite, calcul_param, calcul_p2d_on_surf,
    cherche_edge1, chfi3d_angle_edge, chfi3d_cherche_edge, chfi3d_cherche_face1, chfi3d_cherche_vertex,
    chfi3d_edge_common_faces, chfi3d_search_fd, cp_change_vertex, curve_hermite, geom_lib_build_curve3d,
    indices, orientation_arete_vive_consecutive, orientation_ic_non_vive, orientation_icplus_non_vive,
    parametre_plate, perform_two_corner_same_ext, plate_orientation, remove_sd, remove_surf_data,
    surf_index, Adaptor3dCurveOnSurface, ChFi3dSurfType, GeomPlateBuildPlateSurface,
    GeomPlateCurveConstraint, GeomPlateMakeApprox, GeomPlatePlateG0Criterion, Geom2dAdaptorCurve,
    BRepAdaptorCurve,
};
use super::chfi3d_ds::{
    TopOpeBRepDSCurve, TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference, TopOpeBRepDSKind,
    TopOpeBRepDSPoint, TopOpeBRepDSSolidSurfaceInterference, TopOpeBRepDSSurface,
};
use super::chfi_ds::{ChFiDS_State, ChFiDSRegul, ChFiDSStripe, SharedStripe};

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L1237-3927 — PerformMoreThreeCorner.
// =========================================================================

impl ChFi3dBuilder {
    // OCCT ChFi3d_Builder_CnCrn.cxx L1237-3927
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::too_many_lines)]
    pub fn perform_more_three_corner(&mut self, jndex: usize, nconges: i32) {
        let pi = std::f64::consts::PI;
        // rcad architecture: the BRep context is carried beside the DS
        // facade (OCCT reads geometry through TopoDS handles).
        let brep = self.my_brep.clone();
        // OCCT L1246: const TopoDS_Vertex& V1 = myVDataMap.FindKey(Jndex);
        let v1 = self.my_vdata_map.find_key(jndex).clone();
        // OCCT L1350/L1445: myVDataMap(Jndex) — the stripe list of the top.
        let vdata: Vec<SharedStripe> = self.my_vdata_map.find_from_index(jndex).clone();
        // OCCT L1332: myVFMap(V1).
        let vf_list: Vec<Shape> = if self.my_vf_map.contains(&v1) {
            self.my_vf_map.find(&v1).clone()
        } else {
            Vec::new()
        };
        let tolapp3d = self.tolapp3d;
        let tol2d = self.tol2d;
        let angular = self.angular;

        //    ========================================
        //             Initialisations
        //     ========================================
        // OCCT L1251-1263
        let mut nedge = chfi3d_nb_not_degenerated_edges(&v1, &self.my_ve_map) as i32;
        let (bordlibre, edgelibre1, edgelibre2) = chfi3d_cherche_bords_libres(&self.my_ve_map, &v1);
        let mut droit = false;
        if bordlibre {
            nedge = (nedge - 2) / 2 + 2;
            let angedg = chfi3d_angle_edge(&brep, &v1, &edgelibre1, &edgelibre2).abs();
            droit = (angedg - pi).abs() < 0.01;
        } else {
            nedge /= 2;
        }
        // OCCT L1264-1287 — the (0, size) arrays.
        let size = (nedge * 2) as usize;
        let mut cd: Vec<SharedStripe> = (0..=size)
            .map(|_| Arc::new(RwLock::new(ChFiDSStripe::default())))
            .collect();
        let mut jf = vec![0i32; size + 1];
        let mut index_arr = vec![0i32; size + 1];
        let mut indice_arr = vec![0i32; size + 1];
        let mut sens = vec![0i32; size + 1];
        let mut indcurve3d = vec![0i32; size + 1];
        let mut numfa = vec![vec![0i32; size + 1]; size + 1];
        let mut order = vec![0i32; size + 1];
        let mut i_arr = vec![vec![0i32; size + 1]; size + 1];
        let mut indpoint = vec![vec![0i32; 2]; size + 1];
        let mut oksea = vec![false; size + 1];
        let mut sharp = vec![false; size + 1];
        let mut regul = vec![false; size + 1];
        let mut p_arr = vec![vec![0.0f64; size + 1]; size + 1];
        let mut errapp = vec![0.0f64; size + 1];
        let mut evive: Vec<Shape> = (0..=size).map(|_| Shape::null()).collect();
        let mut fvive: Vec<Vec<Shape>> = (0..=size)
            .map(|_| (0..=size).map(|_| Shape::null()).collect())
            .collect();
        let mut ponctuel = vec![false; size + 1];
        let mut samedge = vec![false; size + 1];
        let mut moresurf = vec![false; size + 1];
        let mut libre = vec![false; size + 1];
        let mut tangentregul = vec![false; size + 1];
        let mut isg1 = vec![false; size + 1];
        // OCCT L1288-1304 — the init loop.
        for ind in 0..=size {
            indpoint[ind][0] = 0;
            indpoint[ind][1] = 0;
            indice_arr[ind] = 0;
            oksea[ind] = false;
            sharp[ind] = false;
            regul[ind] = false;
            ponctuel[ind] = false;
            samedge[ind] = false;
            moresurf[ind] = false;
            libre[ind] = false;
            tangentregul[ind] = false;
            isg1[ind] = false;
        }
        // OCCT L1305-1326 — the scalar locals.
        let mut face = Shape::null();
        let mut icplus: i32;
        #[allow(unused_assignments)]
        let mut icmoins: i32;
        let mut icplus2: i32;
        let mut index = 0i32;
        let mut indice;
        let mut n3d = 0i32;
        let mut ivtx = 0i32;
        let mut nb;
        let mut sameside = false;
        let mut trouve;
        let mut isfirst;
        let mut pardeb;
        let mut parfin;
        let mut xdir;
        let mut ydir;
        let tolapp = 1.0e-4f64;
        let mut maxapp = 0.0f64;
        let mut maxapp1 = 0.0f64;
        let mut avedev = 0.0f64;
        let mut curv3d: Option<rcad_kernel::geom::Curve3> = None;
        let mut regular = ChFiDSRegul::new();
        let mut fproj: Vec<Shape> = Vec::new();
        let mut num;
        let mut ecur = Shape::null();
        let mut f1 = Shape::null();
        let mut f2 = Shape::null();
        let mut sum_face_normal_at_v1 = DVec3::ZERO; // is used to define Plate orientation

        // it is determined if there is a sewing edge
        // OCCT L1329-1340
        let mut couture = false;
        let mut facecouture = Shape::null();
        let mut edgecouture = Shape::null();
        for itf in &vf_list {
            if couture {
                break;
            }
            let fcur = itf.clone();
            let (cout, ecout) = chfi3d_couture_on_vertex(&brep, &fcur, &v1);
            couture = cout;
            if couture {
                edgecouture = ecout;
                facecouture = fcur;
            }
        }

        // unused surfdata are removed
        // OCCT L1343
        remove_surf_data(
            &brep,
            &self.my_vdata_map,
            &self.my_ef_map,
            &edgecouture,
            &facecouture,
            &v1,
        );

        // parse edges and faces
        // OCCT L1345-1540 — the first DStr scope (parse loop + walk).
        trouve = false;
        let mut enext = Shape::null();
        let mut vv = Shape::null();
        let mut fcur = Shape::null();
        let mut fnext = Shape::null();
        {
            let dstr = self.my_ds.as_mut().expect("DS");
            // OCCT L1350-1373
            for it in &vdata {
                if trouve {
                    break;
                }
                cd[0] = it.clone();
                let mut sense = 0i32;
                index_arr[0] =
                    chfi3d_index_of_surf_data(&v1, &it.read().expect("stripe lock"), &mut sense);
                sens[0] = sense;
                numfa[0][1] = surf_index(&cd, 0, index_arr[0], ChFi3dSurfType::FACE2);
                numfa[1][0] = numfa[0][1];
                fcur = dstr.shape(numfa[0][1]).clone();
                fvive[0][1] = fcur.clone();
                fvive[1][0] = fcur.clone();
                jf[0] = 2;
                let st = it.read().expect("stripe lock");
                let spine = st.spine().expect("spine");
                ecur = if sens[0] == 1 {
                    spine.base().edges(1).clone()
                } else {
                    spine.base().edges(spine.base().nb_edges()).clone()
                };
                drop(st);
                evive[0] = ecur.clone();
                chfi3d_cherche_edge(&brep, &v1, &evive, &fcur, &mut enext, &mut vv);
                trouve = !enext.is_null();
            }
            // find sum of all face normals at V1
            // OCCT L1375
            super::chfi3d_builder_cncrn::summarize_normal(
                &brep,
                &v1,
                &fcur,
                &ecur,
                &mut sum_face_normal_at_v1,
            );

            // OCCT L1377-1540 — the edge/face walk.
            let mut nbcouture = 0;
            let mut ii = 1i32;
            while ii < nedge {
                if fcur.is_same(&facecouture) && nbcouture == 0 {
                    enext = edgecouture.clone();
                    nbcouture += 1;
                } else {
                    chfi3d_cherche_edge(&brep, &v1, &evive, &fcur, &mut enext, &mut vv);
                }
                if enext.is_null() {
                    panic!("PerformMoreThreeCorner: pb in the parsing of edges and faces");
                }
                if enext.is_same(&edgelibre1) || enext.is_same(&edgelibre2) {
                    // OCCT L1396-1413 — the first free-border edge.
                    cd[ii as usize] = Arc::new(RwLock::new(ChFiDSStripe::default()));
                    index_arr[ii as usize] = 0;
                    sens[ii as usize] = -1;
                    let mut vref = brep.first_vertex(&enext);
                    if vref.is_same(&v1) {
                        sens[ii as usize] = 1;
                    }
                    sharp[ii as usize] = true;
                    evive[ii as usize] = enext.clone();
                    jf[ii as usize] = 0;
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ii, &mut icplus_l, &mut icmoins_l);
                    fvive[ii as usize][icplus_l as usize] = fcur.clone();
                    fvive[icplus_l as usize][ii as usize] = fcur.clone();
                    numfa[ii as usize][icplus_l as usize] = dstr.add_shape(&fcur);
                    numfa[icplus_l as usize][ii as usize] = numfa[ii as usize][icplus_l as usize];
                    ii += 1;
                    // OCCT L1414-1438 — the second free-border edge.
                    if enext.is_same(&edgelibre1) {
                        ecur = edgelibre2.clone();
                    } else {
                        ecur = edgelibre1.clone();
                    }
                    // OCCT L1422: ChFi3d_edge_common_faces(myEFMap(Ecur),
                    // Fcur, Fcur) — both outputs alias Fcur; the final value
                    // is the second output (F2 = F1 when not found).
                    let mut out1 = Shape::null();
                    let mut out2 = Shape::null();
                    chfi3d_edge_common_faces(self.my_ef_map.find(&ecur), &mut out1, &mut out2);
                    fcur = out2;
                    indices(nedge, ii, &mut icplus_l, &mut icmoins_l);
                    fvive[ii as usize][icplus_l as usize] = fcur.clone();
                    fvive[icplus_l as usize][ii as usize] = fcur.clone();
                    numfa[ii as usize][icplus_l as usize] = dstr.add_shape(&fcur);
                    numfa[icplus_l as usize][ii as usize] = numfa[ii as usize][icplus_l as usize];
                    cd[ii as usize] = Arc::new(RwLock::new(ChFiDSStripe::default()));
                    index_arr[ii as usize] = 0;
                    sens[ii as usize] = -1;
                    vref = brep.first_vertex(&ecur);
                    if vref.is_same(&v1) {
                        sens[ii as usize] = 1;
                    }
                    sharp[ii as usize] = true;
                    evive[ii as usize] = ecur.clone();
                    jf[ii as usize] = 0;
                } else {
                    // it is found if Enext is in the map of stripes
                    // OCCT L1441-1461
                    let mut ee = Shape::null();
                    trouve = false;
                    for it in &vdata {
                        if trouve {
                            break;
                        }
                        let mut sense = 0i32;
                        index =
                            chfi3d_index_of_surf_data(&v1, &it.read().expect("stripe lock"), &mut sense);
                        let st = it.read().expect("stripe lock");
                        let spine = st.spine().expect("spine");
                        ee = if sense == 1 {
                            spine.base().edges(1).clone()
                        } else {
                            spine.base().edges(spine.base().nb_edges()).clone()
                        };
                        if enext.is_same(&ee) {
                            cd[ii as usize] = it.clone();
                            sens[ii as usize] = sense;
                            index_arr[ii as usize] = index;
                            trouve = true;
                        }
                    }
                    if trouve {
                        sharp[ii as usize] = false;
                        evive[ii as usize] = enext.clone();
                    } else {
                        //      edge ii is alive
                        // OCCT L1471-1485
                        cd[ii as usize] = Arc::new(RwLock::new(ChFiDSStripe::default()));
                        index_arr[ii as usize] = 0;
                        sens[ii as usize] = -1;
                        let vref = brep.first_vertex(&enext);
                        if vref.is_same(&v1) {
                            sens[ii as usize] = 1;
                        }
                        sharp[ii as usize] = true;
                        evive[ii as usize] = enext.clone();
                        jf[ii as usize] = 0;
                    }
                    // Face Fnext!=Fcur containing Enext
                    // OCCT L1486-1493
                    fnext = fcur.clone();
                    chfi3d_cherche_face1(self.my_ef_map.find(&enext), &fcur, &mut fnext);
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ii, &mut icplus_l, &mut icmoins_l);
                    fvive[ii as usize][icplus_l as usize] = fnext.clone();
                    fvive[icplus_l as usize][ii as usize] = fnext.clone();
                    numfa[ii as usize][icplus_l as usize] = dstr.add_shape(&fnext);
                    numfa[icplus_l as usize][ii as usize] = numfa[ii as usize][icplus_l as usize];
                    // OCCT L1494-1534 — the cts16288 jf/numfa re-update.
                    let mut numface1 = 0i32;
                    let mut numface2 = 0i32;
                    if trouve {
                        // it is checked if numfa corresponds to IndexOfS1 or
                        // IndexOfS2 (cts16288)
                        numface2 =
                            surf_index(&cd, ii as usize, index_arr[ii as usize], ChFi3dSurfType::FACE2);
                        if numface2 == numfa[ii as usize][icplus_l as usize] {
                            jf[ii as usize] = 2;
                        } else {
                            numface1 =
                                surf_index(&cd, ii as usize, index_arr[ii as usize], ChFi3dSurfType::FACE1);
                            if numface1 == numfa[ii as usize][icplus_l as usize] {
                                jf[ii as usize] = 1;
                            } else {
                                if numface1 == numfa[icmoins_l as usize][ii as usize] {
                                    jf[ii as usize] = 2;
                                    fvive[ii as usize][icplus_l as usize] = dstr.shape(numface2).clone();
                                    fvive[icplus_l as usize][ii as usize] = dstr.shape(numface2).clone();
                                    let shape2 = dstr.shape(numface2).clone();
                                    numfa[ii as usize][icplus_l as usize] = dstr.add_shape(&shape2);
                                    numfa[icplus_l as usize][ii as usize] =
                                        numfa[ii as usize][icplus_l as usize];
                                }
                                if numface2 == numfa[icmoins_l as usize][ii as usize] {
                                    jf[ii as usize] = 1;
                                    fvive[ii as usize][icplus_l as usize] = dstr.shape(numface1).clone();
                                    fvive[icplus_l as usize][ii as usize] = dstr.shape(numface1).clone();
                                    let shape1 = dstr.shape(numface1).clone();
                                    numfa[ii as usize][icplus_l as usize] = dstr.add_shape(&shape1);
                                    numfa[icplus_l as usize][ii as usize] =
                                        numfa[ii as usize][icplus_l as usize];
                                }
                            }
                        }
                    }
                    ecur = enext.clone();
                    fcur = fnext.clone();
                    // find sum of all face normales at V1
                    // OCCT L1538
                    super::chfi3d_builder_cncrn::summarize_normal(
                        &brep,
                        &v1,
                        &fcur,
                        &ecur,
                        &mut sum_face_normal_at_v1,
                    );
                }
                ii += 1;
            }
        } // end of the first DStr scope (OCCT L1345-1540)

        // mise a jour du tableau regul
        // OCCT L1541-1557
        for ic in 0..nedge {
            if sharp[ic as usize] {
                ecur = evive[ic as usize].clone();
                if !ecur.is_same(&edgecouture) {
                    chfi3d_edge_common_faces(self.my_ef_map.find(&ecur), &mut f1, &mut f2);
                    //  Modified by Sergey KHROMOV - Fri Dec 21 18:11:02 2001
                    regul[ic as usize] = is_tangent_faces(
                        &brep,
                        &ecur,
                        &f1,
                        &f2,
                        crate::geomalgo::gtests_stubs::GeomAbsShape::G1,
                    );
                    //  Modified by Sergey KHROMOV - Fri Dec 21 18:11:07 2001
                }
            }
        }
        // it is checked if a regular edge is not tangent to another edge
        // in case if it is not considered regular (cts60072)
        // OCCT L1558-1581
        for ic in 0..nedge {
            if regul[ic as usize] {
                trouve = false;
                let ereg = evive[ic as usize].clone();
                let mut ind = 0i32;
                while ind < nedge && !trouve {
                    if ind != ic {
                        let ecur2 = evive[ind as usize].clone();
                        let ang = chfi3d_angle_edge(&brep, &v1, &ecur2, &ereg).abs();
                        if ang < 0.01 || (ang - pi).abs() < 0.01 {
                            regul[ic as usize] = false;
                            tangentregul[ic as usize] = true;
                            trouve = true;
                        }
                    }
                    ind += 1;
                }
            }
        }

        // variable deuxconges allows detecting cases when there is a top
        // with n edges and two fillets on two tangent edges that are not
        // free borders the connecting curves start from the fillet and end
        // on top
        // OCCT L1583-1604
        let mut deuxconges = false;
        let deuxcgnontg;
        trouve = false;
        if nconges == 2 {
            let mut ic = 0i32;
            while ic < nedge && !trouve {
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                if !sharp[ic as usize] && !sharp[icplus_l as usize] {
                    let e1 = evive[ic as usize].clone();
                    let e2 = evive[icplus_l as usize].clone();
                    deuxconges = chfi3d_angle_edge(&brep, &v1, &e1, &e2).abs() < 0.01;
                    trouve = deuxconges;
                }
                ic += 1;
            }
        }

        // variable deuxconges is used in the special case when there are
        // two fillets and if two other living edges are tangent (cts60072)
        // OCCT L1606-1626
        if nconges == 2 && nedge == 4 {
            let mut ic = 0i32;
            while ic < nedge && !deuxconges {
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                if sharp[ic as usize] && sharp[icplus_l as usize] {
                    let e1 = evive[ic as usize].clone();
                    let e2 = evive[icplus_l as usize].clone();
                    if !e1.is_same(&edgelibre1)
                        && !e1.is_same(&edgelibre2)
                        && !e2.is_same(&edgelibre1)
                        && !e2.is_same(&edgelibre2)
                    {
                        let ang = chfi3d_angle_edge(&brep, &v1, &e1, &e2).abs();
                        deuxconges = ang < 0.01 || (ang - pi).abs() < 0.01;
                    }
                }
                ic += 1;
            }
        }

        // OCCT L1628
        deuxcgnontg = nconges == 2 && nedge == 3 && !deuxconges; // pro12305

        if deuxconges {
            for ic in 0..nedge {
                regul[ic as usize] = false;
            }
        }

        // Detect case of 3 edges & 2 conges: OnSame + OnDiff
        // (eap, Arp 9 2002, occ266)
        // OCCT L1638-1670
        let mut is_on_same_diff = false;
        if deuxcgnontg {
            let mut is_on_same = false;
            let mut is_on_diff = false;
            for ic in 0..nedge {
                if sharp[ic as usize] {
                    continue;
                }
                let stat = if sens[ic as usize] == 1 {
                    cd[ic as usize]
                        .read()
                        .expect("stripe lock")
                        .spine()
                        .expect("spine")
                        .base()
                        .first_status()
                } else {
                    cd[ic as usize]
                        .read()
                        .expect("stripe lock")
                        .spine()
                        .expect("spine")
                        .base()
                        .last_status()
                };

                if stat == ChFiDS_State::OnSame {
                    is_on_same = true;
                } else if stat == ChFiDS_State::OnDiff {
                    is_on_diff = true;
                }
            }
            is_on_same_diff = is_on_same && is_on_diff;
        }
        if is_on_same_diff {
            // OCCT L1671-1677
            self.perform_one_corner(jndex, true);
        }

        // ============ second DStr scope (OCCT L1684-3927) ============
        let dstr = self.my_ds.as_mut().expect("DS");

        // if the commonpoint is on an edge that does not have a
        // vertex at the extremity, Evive is found anew
        // Fvive is found anew if it does not correspond
        // to two faces adjacent to Evive (cts16288)
        // OCCT L1684-1783
        if !deuxconges && !is_on_same_diff {
            for ic in 0..nedge {
                if sharp[ic as usize] {
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                    let icp = ic as usize;
                    let icpl = icplus_l as usize;
                    let icm = icmoins_l as usize;
                    let arc = evive[icp].clone();
                    let mut angedg = pi;
                    let mut vcom = Shape::null();
                    if !sharp[icpl] {
                        isfirst = sens[icpl] == 1;
                        let jfp = 3 - jf[icpl];
                        let cp1 = {
                            let st = cd[icpl].read().expect("stripe lock");
                            let fd =
                                st.set_of_surf_data()[(index_arr[icpl] - 1) as usize].clone();
                            cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, jfp)
                        };
                        if cp1.is_on_arc() {
                            chfi3d_cherche_vertex(&arc, cp1.arc(), &mut vcom, &mut trouve);
                            if trouve {
                                angedg = chfi3d_angle_edge(&brep, &vcom, &arc, cp1.arc()).abs();
                            }
                            if !cp1.arc().is_same(&arc) && (angedg - pi).abs() < 0.01 {
                                evive[icp] = cp1.arc().clone();
                                chfi3d_edge_common_faces(self.my_ef_map.find(cp1.arc()), &mut f1, &mut f2);
                                if !fvive[icp][icpl].is_same(&f1) && !fvive[icp][icpl].is_same(&f2) {
                                    if fvive[icp][icm].is_same(&f2) {
                                        fvive[icp][icpl] = f1.clone();
                                        fvive[icpl][icp] = f1.clone();
                                        numfa[icp][icpl] = dstr.add_shape(&f1);
                                        numfa[icpl][icp] = dstr.add_shape(&f1);
                                    } else {
                                        fvive[icp][icpl] = f2.clone();
                                        fvive[icpl][icp] = f2.clone();
                                        numfa[icp][icpl] = dstr.add_shape(&f2);
                                        numfa[icpl][icp] = dstr.add_shape(&f2);
                                    }
                                }
                                samedge[icp] = true;
                                p_arr[icp][icpl] = cp1.parameter_on_arc();
                                p_arr[icp][icm] = cp1.parameter_on_arc();
                                i_arr[icp][icpl] = 1;
                            }
                        }
                    }
                    if !sharp[icm] {
                        isfirst = sens[icm] == 1;
                        let cp2 = {
                            let st = cd[icm].read().expect("stripe lock");
                            let fd = st.set_of_surf_data()[(index_arr[icm] - 1) as usize].clone();
                            cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, jf[icm])
                        };
                        if cp2.is_on_arc() {
                            angedg = pi;
                            chfi3d_cherche_vertex(&arc, cp2.arc(), &mut vcom, &mut trouve);
                            if trouve {
                                angedg = chfi3d_angle_edge(&brep, &vcom, &arc, cp2.arc()).abs();
                            }
                            if !cp2.arc().is_same(&arc) && (angedg - pi).abs() < 0.01 {
                                evive[icp] = cp2.arc().clone();
                                chfi3d_edge_common_faces(self.my_ef_map.find(cp2.arc()), &mut f1, &mut f2);
                                if !fvive[icp][icm].is_same(&f1) && !fvive[icp][icm].is_same(&f2) {
                                    if fvive[icp][icpl].is_same(&f2) {
                                        fvive[icp][icm] = f1.clone();
                                        numfa[icp][icm] = dstr.add_shape(&f1);
                                        fvive[icm][icp] = f1.clone();
                                        numfa[icm][icp] = dstr.add_shape(&f1);
                                    } else {
                                        fvive[icp][icm] = f2.clone();
                                        numfa[icp][icm] = dstr.add_shape(&f2);
                                        fvive[icm][icp] = f2.clone();
                                        numfa[icm][icp] = dstr.add_shape(&f2);
                                    }
                                }
                                samedge[icp] = true;
                                p_arr[icp][icm] = cp2.parameter_on_arc();
                                p_arr[icp][icpl] = cp2.parameter_on_arc();
                                i_arr[icp][icm] = 1;
                            }
                        }
                    }
                }
            }
        }

        // the first free edge is restored if it exists
        // OCCT L1785-1796
        trouve = false;
        {
            let mut ic = 0i32;
            while ic < nedge && !trouve {
                let ecom = evive[ic as usize].clone();
                if ecom.is_same(&edgelibre1) || ecom.is_same(&edgelibre2) {
                    libre[ic as usize] = true;
                    trouve = true;
                }
                ic += 1;
            }
        }

        // determine the minimum recoil distance that can't be exceeded
        // OCCT L1798-1830
        let mut distmini = false;
        let som = brep.vertex_position(&v1);
        let mut pic;
        let mut p2 = DVec2::ZERO;
        let mut edgemin;
        let mut distmin = 1.0e30f64;
        for ic in 0..nedge {
            if sharp[ic as usize] {
                edgemin = evive[ic as usize].clone();
            } else {
                let st = cd[ic as usize].read().expect("stripe lock");
                let spine = st.spine().expect("spine");
                edgemin = if sens[ic as usize] == 1 {
                    spine.base().edges(1).clone()
                } else {
                    spine.base().edges(spine.base().nb_edges()).clone()
                };
            }
            let v = brep.first_vertex(&edgemin);
            let v2 = brep.last_vertex(&edgemin);
            let dst = brep.vertex_position(&v).distance(brep.vertex_position(&v2)) / 1.5;
            if dst < distmin {
                distmin = dst;
            }
        }

        //  calculate intersections between stripes and determine the
        //  parameters on each pcurve
        // OCCT L1832-2000
        let mut inters = true;
        for ic in 0..nedge {
            let mut icplus_l = 0i32;
            let mut icmoins_l = 0i32;
            indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
            let icp = ic as usize;
            let icpl = icplus_l as usize;
            let _ = icmoins_l;
            if sharp[icp] || sharp[icpl] {
                oksea[icp] = false;
            } else {
                let mut jf1 = 0i32;
                let mut jfp = 0i32;
                let mut i1 = 0i32;
                let mut i2 = 0i32;
                let mut pa1 = 0.0f64;
                let mut pa2 = 0.0f64;
                // if two edges are tangent the intersection is not
                // attempted (cts60046)
                let angedg = chfi3d_angle_edge(&brep, &v1, &evive[icp], &evive[icpl]).abs();
                let ok = if (angedg - pi).abs() > 0.01 {
                    chfi3d_search_fd(
                        &brep,
                        dstr,
                        &cd[icp],
                        &cd[icpl],
                        sens[icp],
                        sens[icpl],
                        &mut i1,
                        &mut i2,
                        &mut pa1,
                        &mut pa2,
                        index_arr[icp],
                        index_arr[icpl],
                        &mut face,
                        &mut sameside,
                        &mut jf1,
                        &mut jfp,
                    )
                } else {
                    false
                };
                // if there is an intersection it is checked if surfdata with
                // the intersection corresponds to the first or the last
                // if this is not the case, the surfdata are removed from SD
                if ok {
                    if i1 != index_arr[icp] {
                        let ideb;
                        let ifin;
                        if sens[icp] == 1 {
                            ideb = index_arr[icp];
                            ifin = i1 - 1;
                        } else {
                            ifin = index_arr[icp];
                            ideb = i1 + 1;
                        }
                        if i1 < index_arr[icp] {
                            nb = index_arr[icp];
                            while nb >= i1 {
                                let iface = if (3 - jf1) == 1 {
                                    surf_index(&cd, icp, nb, ChFi3dSurfType::FACE1)
                                } else {
                                    surf_index(&cd, icp, nb, ChFi3dSurfType::FACE2)
                                };
                                fproj.push(dstr.shape(iface).clone());
                                nb -= 1;
                            }
                        }
                        if i1 > index_arr[icp] {
                            nb = index_arr[icp];
                            while nb <= i1 {
                                let iface = if (3 - jf1) == 1 {
                                    surf_index(&cd, icp, nb, ChFi3dSurfType::FACE1)
                                } else {
                                    surf_index(&cd, icp, nb, ChFi3dSurfType::FACE2)
                                };
                                fproj.push(dstr.shape(iface).clone());
                                nb += 1;
                            }
                        }
                        let strip = cd[icp].clone();
                        let mut w = strip.write().expect("stripe lock");
                        remove_sd(&mut w, ideb, ifin);
                        drop(w);
                        let mut sense = 0i32;
                        num = chfi3d_index_of_surf_data(
                            &v1,
                            &cd[icp].read().expect("stripe lock"),
                            &mut sense,
                        );
                        index_arr[icp] = num;
                        i1 = num;
                    }
                    if i2 != index_arr[icpl] {
                        let ideb;
                        let ifin;
                        if sens[icpl] == 1 {
                            ideb = index_arr[icpl];
                            ifin = i2 - 1;
                        } else {
                            ifin = index_arr[icpl];
                            ideb = i2 + 1;
                        }

                        if i2 < index_arr[icpl] {
                            nb = i2;
                            while nb <= index_arr[icpl] {
                                let iface = if (3 - jfp) == 1 {
                                    surf_index(&cd, icpl, nb, ChFi3dSurfType::FACE1)
                                } else {
                                    surf_index(&cd, icpl, nb, ChFi3dSurfType::FACE2)
                                };
                                fproj.push(dstr.shape(iface).clone());
                                nb += 1;
                            }
                        }
                        if i2 > index_arr[icpl] {
                            nb = i2;
                            while nb >= index_arr[icpl] {
                                let iface = if (3 - jfp) == 1 {
                                    surf_index(&cd, icpl, nb, ChFi3dSurfType::FACE1)
                                } else {
                                    surf_index(&cd, icpl, nb, ChFi3dSurfType::FACE2)
                                };
                                fproj.push(dstr.shape(iface).clone());
                                nb -= 1;
                            }
                        }
                        let strip = cd[icpl].clone();
                        let mut w = strip.write().expect("stripe lock");
                        remove_sd(&mut w, ideb, ifin);
                        drop(w);
                        let mut sense = 0i32;
                        num = chfi3d_index_of_surf_data(
                            &v1,
                            &cd[icpl].read().expect("stripe lock"),
                            &mut sense,
                        );
                        index_arr[icpl] = num;
                        i2 = num;
                    }
                    calcul_p2d_on_surf(
                        &cd[icp].read().expect("stripe lock"),
                        jf1,
                        i1,
                        pa1,
                        &mut p2,
                    );
                    indice = surf_index(&cd, icp, i1, ChFi3dSurfType::ChFiSURFACE);
                    pic = dstr.surface(indice).surface.point_at(p2.x, p2.y);
                    if pic.distance(som) > distmin {
                        distmini = true;
                    }
                    jf[icp] = jf1;
                    i_arr[icp][icpl] = i1;
                    i_arr[icpl][icp] = i2;
                    p_arr[icp][icpl] = pa1;
                    p_arr[icpl][icp] = pa2;
                }
                oksea[icp] = ok;
            }
            if !oksea[icp] {
                inters = false;
            }
        }

        // case if there are only intersections
        // the parametres on Pcurves are the extremities of the stripe
        // OCCT L2002-2106
        let mut para = 0.0f64;
        if !inters {
            for ic in 0..nedge {
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                let mut icplus2_l = 0i32;
                indices(nedge, icplus_l, &mut icplus2_l, &mut icmoins_l);
                let icp = ic as usize;
                let icpl = icplus_l as usize;
                if !oksea[icp] {
                    if sharp[icp] {
                        if !samedge[icp] {
                            para = brep_tool_parameter(&brep, &v1, &evive[icp]);
                            p_arr[icp][icpl] = para;
                            i_arr[icp][icpl] = 1;
                        }
                    } else {
                        isfirst = sens[icp] == 1;
                        let mut sense = 0i32;
                        i_arr[icp][icpl] = chfi3d_index_of_surf_data(
                            &v1,
                            &cd[icp].read().expect("stripe lock"),
                            &mut sense,
                        );
                        if oksea[icmoins_l as usize] {
                            para = p_arr[icp][icmoins_l as usize];
                            p_arr[icp][icpl] = para;
                        } else {
                            calcul_param(
                                &cd[icp].read().expect("stripe lock"),
                                jf[icp],
                                i_arr[icp][icpl],
                                isfirst,
                                &mut para,
                            );
                            p_arr[icp][icpl] = para;
                        }
                    }
                    if sharp[icpl] {
                        if !samedge[icpl] {
                            para = brep_tool_parameter(&brep, &v1, &evive[icpl]);
                            p_arr[icpl][icp] = para;
                            i_arr[icpl][icp] = 1;
                        }
                    } else {
                        isfirst = sens[icpl] == 1;
                        let mut sense = 0i32;
                        i_arr[icpl][icp] = chfi3d_index_of_surf_data(
                            &v1,
                            &cd[icpl].read().expect("stripe lock"),
                            &mut sense,
                        );
                        if oksea[icpl] {
                            para = p_arr[icpl][icplus2_l as usize];
                            p_arr[icpl][icp] = para;
                        } else {
                            let jfp = 3 - jf[icpl];
                            calcul_param(
                                &cd[icpl].read().expect("stripe lock"),
                                jfp,
                                i_arr[icpl][icp],
                                isfirst,
                                &mut para,
                            );
                            p_arr[icpl][icp] = para;
                        }
                    }
                }
            }

            //  calculate max distance to the top at each point
            // OCCT L2065-2106
            let mut dist1 = vec![0.0f64; size + 1];
            let mut dist2 = vec![0.0f64; size + 1];
            let mut distance = 0.0f64;
            let sommet = brep.vertex_position(&v1);
            if !deuxconges {
                for ic in 0..nedge {
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                    let icp = ic as usize;
                    let icm = icmoins_l as usize;
                    let icpl = icplus_l as usize;
                    if sharp[icp] {
                        dist1[icp] = 0.0;
                        dist2[icp] = 0.0;
                    } else {
                        let jfp = 3 - jf[icp];
                        calcul_p2d_on_surf(
                            &cd[icp].read().expect("stripe lock"),
                            jfp,
                            i_arr[icp][icm],
                            p_arr[icp][icm],
                            &mut p2,
                        );
                        indice = surf_index(&cd, icp, i_arr[icp][icm], ChFi3dSurfType::ChFiSURFACE);
                        pic = dstr.surface(indice).surface.point_at(p2.x, p2.y);
                        dist1[icp] = sommet.distance(pic);
                        if dist1[icp] > distance {
                            distance = dist1[icp];
                        }

                        calcul_p2d_on_surf(
                            &cd[icp].read().expect("stripe lock"),
                            jf[icp],
                            i_arr[icp][icpl],
                            p_arr[icp][icpl],
                            &mut p2,
                        );
                        indice = surf_index(&cd, icp, i_arr[icp][icpl], ChFi3dSurfType::ChFiSURFACE);
                        pic = dstr.surface(indice).surface.point_at(p2.x, p2.y);
                        dist2[icp] = sommet.distance(pic);
                        if dist2[icp] > distance {
                            distance = dist2[icp];
                        }
                    }
                }
            }

            //  offset of parameters and removal of intersection points
            //  too close to the top
            // OCCT L2108-2196
            let mut ec;
            let mut dist;
            if !deuxconges && !deuxcgnontg {
                for ic in 0..nedge {
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                    let icp = ic as usize;
                    let icm = icmoins_l as usize;
                    let icpl = icplus_l as usize;
                    if sharp[icp] {
                        let c = BRepAdaptorCurve::initialize(&brep, &evive[icp]);
                        // to pass from 3D distance to a parametric distance
                        if !tangentregul[icp] {
                            ec = distance * 100.0 * c.resolution(0.01);
                        } else {
                            ec = 0.0;
                        }
                        if brep.first_vertex(&evive[icp]).is_same(&v1) {
                            para = p_arr[icp][icm] + ec;
                            p_arr[icp][icm] = para;
                        } else {
                            para = p_arr[icp][icm] - ec;
                            p_arr[icp][icm] = para;
                        }
                        // it is necessary to be on to remain on the edge
                        p_arr[icp][icpl] = p_arr[icp][icm];
                    } else if !distmini {
                        dist = dist1[icp];
                        if !oksea[icm] || (oksea[icm] && distance > 1.3 * dist) {
                            ec = distance - dist;
                            if oksea[icm] {
                                oksea[icm] = false;
                                inters = false;
                            }
                            if sens[icp] == 1 {
                                para = p_arr[icp][icm] + ec;
                                p_arr[icp][icm] = para;
                            } else {
                                para = p_arr[icp][icm] - ec;
                                p_arr[icp][icm] = para;
                            }
                        }
                        dist = dist2[icp];
                        if !oksea[icp] || (oksea[icp] && distance > 1.3 * dist) {
                            if oksea[icp] {
                                oksea[icp] = false;
                                inters = false;
                            }
                            if nconges != 1 {
                                let parold = p_arr[icp][icpl];
                                let parnew = p_arr[icp][icm];
                                if sens[icp] == 1 {
                                    if parnew > parold {
                                        p_arr[icp][icpl] = p_arr[icp][icm];
                                    }
                                } else if parnew < parold {
                                    p_arr[icp][icpl] = p_arr[icp][icm];
                                }
                            }
                        }
                    }
                }
            }
        }

        // it is attempted to limit the edge by a commonpoint
        // OCCT L2198-2317
        let mut tolcp = 0.0f64;
        let sommet = brep.vertex_position(&v1);
        if !deuxconges {
            for ic in 0..nedge {
                if sharp[ic as usize] {
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                    let icp = ic as usize;
                    let icpl = icplus_l as usize;
                    let icm = icmoins_l as usize;
                    let c = BRepAdaptorCurve::initialize(&brep, &evive[icp]);
                    let pe = c.value(p_arr[icp][icpl]);
                    let ds = pe.distance(sommet);
                    let mut cp1: Option<super::chfi_ds::ChFiDS_CommonPoint> = None;
                    let mut cp2: Option<super::chfi_ds::ChFiDS_CommonPoint> = None;
                    if !sharp[icpl] {
                        isfirst = sens[icpl] == 1;
                        let jfp = 3 - jf[icpl];
                        let st = cd[icpl].read().expect("stripe lock");
                        let fd = st.set_of_surf_data()[(i_arr[icpl][icp] - 1) as usize].clone();
                        cp1 = Some(cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, jfp));
                    }
                    if !sharp[icm] {
                        isfirst = sens[icm] == 1;
                        let st = cd[icm].read().expect("stripe lock");
                        let fd = st.set_of_surf_data()[(i_arr[icm][icp] - 1) as usize].clone();
                        cp2 = Some(cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, jf[icm]));
                    }
                    let d1 = cp1
                        .as_ref()
                        .map(|c| c.point().distance(sommet))
                        .unwrap_or(0.0);
                    let d2 = cp2
                        .as_ref()
                        .map(|c| c.point().distance(sommet))
                        .unwrap_or(0.0);
                    let mut samecompoint = false;
                    if !sharp[icm] && !sharp[icpl] {
                        samecompoint = cp1.as_ref().expect("cp1").point().distance(
                            cp2.as_ref().expect("cp2").point(),
                        ) < tolapp;
                    }
                    if (ds < d1 || ds < d2) && !samecompoint {
                        // step back till Common Points
                        // without leaving the Edge ??
                        if d2 < d1 && cp1.as_ref().expect("cp1").is_on_arc() {
                            // cp1 is chosen
                            p_arr[icp][icm] = cp1.as_ref().expect("cp1").parameter_on_arc();
                            p_arr[icp][icpl] = p_arr[icp][icm];
                            isfirst = sens[icpl] == 1;
                            let jfp = 3 - jf[icpl];
                            calcul_param(
                                &cd[icpl].read().expect("stripe lock"),
                                jfp,
                                i_arr[icpl][icp],
                                isfirst,
                                &mut para,
                            );
                            p_arr[icpl][icp] = para;
                            let cp1t = cp1.as_ref().expect("cp1").tolerance();
                            if cp1t > tolcp && cp1t < 1.0 {
                                tolcp = cp1t;
                            }
                        } else if cp2.as_ref().expect("cp2").is_on_arc() {
                            // cp2 is chosen
                            p_arr[icp][icm] = cp2.as_ref().expect("cp2").parameter_on_arc();
                            p_arr[icp][icpl] = p_arr[icp][icm];
                            isfirst = sens[icm] == 1;
                            calcul_param(
                                &cd[icm].read().expect("stripe lock"),
                                jf[icm],
                                i_arr[icm][icp],
                                isfirst,
                                &mut para,
                            );
                            p_arr[icm][icp] = para;
                            let cp2t = cp2.as_ref().expect("cp2").tolerance();
                            if cp2t > tolcp && cp2t < 1.0 {
                                tolcp = cp2t;
                            }
                        }
                    } else {
                        // step back till Common Point only if it is very close
                        if !sharp[icpl] {
                            let cp1v = cp1.as_ref().expect("cp1");
                            if (cp1v.point().distance(pe) < cp1v.tolerance()
                                || samecompoint
                                || nconges == 1)
                                && cp1v.is_on_arc()
                            {
                                // it is very close to cp1
                                p_arr[icp][icm] = cp1v.parameter_on_arc();
                                ponctuel[icp] = true;
                                p_arr[icp][icpl] = p_arr[icp][icm];
                                isfirst = sens[icpl] == 1;
                                let jfp = 3 - jf[icpl];
                                calcul_param(
                                    &cd[icpl].read().expect("stripe lock"),
                                    jfp,
                                    i_arr[icpl][icp],
                                    isfirst,
                                    &mut para,
                                );
                                p_arr[icpl][icp] = para;
                                if cp1v.tolerance() > tolcp && cp1v.tolerance() < 1.0 {
                                    tolcp = cp1v.tolerance();
                                }
                            }
                        }
                        if !sharp[icm] {
                            let cp2v = cp2.as_ref().expect("cp2");
                            if (cp2v.point().distance(pe) < cp2v.tolerance()
                                || samecompoint
                                || nconges == 1)
                                && cp2v.is_on_arc()
                            {
                                // it is very close to cp2
                                ponctuel[icm] = true;
                                p_arr[icp][icm] = cp2v.parameter_on_arc();
                                p_arr[icp][icpl] = p_arr[icp][icm];
                                isfirst = sens[icm] == 1;
                                calcul_param(
                                    &cd[icm].read().expect("stripe lock"),
                                    jf[icm],
                                    i_arr[icm][icp],
                                    isfirst,
                                    &mut para,
                                );
                                p_arr[icm][icp] = para;
                                if cp2v.tolerance() > tolcp && cp2v.tolerance() < 1.0 {
                                    tolcp = cp2v.tolerance();
                                }
                            }
                        }
                    }
                }
            }
        }

        // in case of a free border the parameter corresponding
        // to the common point on the free edge is chosen.
        // OCCT L2319-2367
        for ic in 0..nedge {
            if evive[ic as usize].is_same(&edgelibre1) || evive[ic as usize].is_same(&edgelibre2) {
                let indic;
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                let icp = ic as usize;
                if libre[icp] {
                    indic = icmoins_l;
                } else {
                    indic = icplus_l;
                }
                if !sharp[indic as usize] {
                    isfirst = sens[indic as usize] == 1;
                    let mut cp1 = {
                        let st = cd[indic as usize].read().expect("stripe lock");
                        let fd =
                            st.set_of_surf_data()[(index_arr[indic as usize] - 1) as usize].clone();
                        cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, 1)
                    };
                    trouve = false;
                    if cp1.is_on_arc() {
                        if cp1.arc().is_same(&evive[icp]) {
                            p_arr[icp][icmoins_l as usize] = cp1.parameter_on_arc();
                            p_arr[icp][icplus_l as usize] = cp1.parameter_on_arc();
                            trouve = true;
                        }
                    }
                    if !trouve {
                        let st = cd[indic as usize].read().expect("stripe lock");
                        let fd =
                            st.set_of_surf_data()[(index_arr[indic as usize] - 1) as usize].clone();
                        cp1 = cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, 2);
                        if cp1.is_on_arc() {
                            if cp1.arc().is_same(&evive[icp]) {
                                p_arr[icp][icmoins_l as usize] = cp1.parameter_on_arc();
                                p_arr[icp][icplus_l as usize] = cp1.parameter_on_arc();
                            }
                        }
                    }
                }
            }
        }

        // if ic is a regular edge, one finds edge indfin which is not
        // a regular edge, and construtc a curve 3d
        // between edges (or stripes ) icmoins and indfin.
        // Then this courbe3d is projected on all faces (nbface) that
        // separate icmoins and indfin
        // OCCT L2369-2540
        let mut nbface = 0i32;
        let mut error = 0.0f64;
        let mut proj2d1: Vec<Option<rcad_kernel::geom::Curve2d>> = (0..=size).map(|_| None).collect();
        let mut proj2d2: Vec<Option<rcad_kernel::geom::Curve2d>> = (0..=size).map(|_| None).collect();
        let mut cproj1: Vec<Option<rcad_kernel::geom::Curve3>> = (0..=size).map(|_| None).collect();
        let mut cproj2: Vec<Option<rcad_kernel::geom::Curve3>> = (0..=size).map(|_| None).collect();
        if !deuxconges {
            let mut ic = 0i32;
            while ic < nedge {
                if regul[ic as usize] {
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                    let mut icplus2_l = 0i32;
                    indices(nedge, icplus_l, &mut icplus2_l, &mut icmoins_l);
                    let icp = ic as usize;
                    let icpl = icplus_l as usize;
                    let icm = icmoins_l as usize;
                    let mut indfin = icpl as i32;
                    trouve = false;
                    let mut ii = icplus_l;
                    while !trouve {
                        if !regul[ii as usize] {
                            indfin = ii;
                            trouve = true;
                        }
                        if ii == nedge - 1 {
                            ii = 0;
                        } else {
                            ii += 1;
                        }
                    }
                    let mut indfinplus = 0i32;
                    let mut indfinmoins = 0i32;
                    indices(nedge, indfin, &mut indfinplus, &mut indfinmoins);
                    let _ = indfinplus;
                    let mut lface: Vec<Shape> = Vec::new();
                    if !sharp[icm] {
                        let ilin = if jf[icm] == 1 {
                            surf_index(&cd, icm, i_arr[icm][icp], ChFi3dSurfType::FACE1)
                        } else {
                            surf_index(&cd, icm, i_arr[icm][icp], ChFi3dSurfType::FACE2)
                        };
                        lface.push(dstr.shape(ilin).clone());
                    } else {
                        lface.push(fvive[icp][icm].clone());
                    }
                    if indfin > icmoins_l {
                        nbface = indfin - icmoins_l;
                    } else {
                        nbface = nedge - (icmoins_l - indfin);
                    }
                    let mut ledge: Vec<Shape> = Vec::new();
                    let mut seqpr: Vec<f64> = Vec::new();
                    ii = ic;
                    for nf in 1..=(nbface - 1) {
                        let mut iimoins = 0i32;
                        let mut iiplus = 0i32;
                        indices(nedge, ii, &mut iiplus, &mut iimoins);
                        ledge.push(evive[ii as usize].clone());
                        seqpr.push(p_arr[ii as usize][iiplus as usize]);
                        if nf != nbface - 1 {
                            lface.push(fvive[ii as usize][iiplus as usize].clone());
                        }
                        if ii == nedge - 1 {
                            ii = 0;
                        } else {
                            ii += 1;
                        }
                    }
                    if !sharp[indfin as usize] {
                        let jfp = 3 - jf[indfin as usize];
                        let ilin = if jfp == 1 {
                            surf_index(
                                &cd,
                                indfin as usize,
                                i_arr[indfin as usize][indfinmoins as usize],
                                ChFi3dSurfType::FACE1,
                            )
                        } else {
                            surf_index(
                                &cd,
                                indfin as usize,
                                i_arr[indfin as usize][indfinmoins as usize],
                                ChFi3dSurfType::FACE2,
                            )
                        };
                        lface.push(dstr.shape(ilin).clone());
                    } else {
                        lface.push(fvive[indfin as usize][indfinmoins as usize].clone());
                    }
                    // OCCT L2439: Epj — the projected-edge sequence.
                    let mut epj: Vec<Shape> = Vec::new();
                    let mut pr: Vec<Option<rcad_kernel::geom::Curve2d>> = Vec::new();
                    let mut cr: Vec<Option<rcad_kernel::geom::Curve3>> = Vec::new();
                    curve_hermite(
                        &brep,
                        dstr,
                        &cd[icm].read().expect("stripe lock"),
                        jf[icm],
                        i_arr[icm][icp],
                        p_arr[icm][icp],
                        sens[icm],
                        sharp[icm],
                        &evive[icm],
                        &cd[indfin as usize].read().expect("stripe lock"),
                        jf[indfin as usize],
                        i_arr[indfin as usize][indfinmoins as usize],
                        p_arr[indfin as usize][indfinmoins as usize],
                        sens[indfin as usize],
                        sharp[indfin as usize],
                        &evive[indfin as usize],
                        nbface,
                        &mut ledge,
                        &lface,
                        &mut pr,
                        &mut cr,
                        &mut epj,
                        &mut seqpr,
                        &mut error,
                    );
                    ii = ic;
                    for ind in 1..=(nbface - 1) {
                        let mut iimoins = 0i32;
                        let mut iiplus = 0i32;
                        indices(nedge, ii, &mut iiplus, &mut iimoins);
                        p_arr[ii as usize][iiplus as usize] = seqpr[(ind - 1) as usize];
                        p_arr[ii as usize][iimoins as usize] = seqpr[(ind - 1) as usize];
                        proj2d1[ii as usize] = pr[(ind - 1) as usize].clone();
                        proj2d2[ii as usize] = pr[ind as usize].clone();
                        cproj1[ii as usize] = cr[(ind - 1) as usize].clone();
                        cproj2[ii as usize] = cr[ind as usize].clone();
                        if ii == nedge - 1 {
                            ii = 0;
                        } else {
                            ii += 1;
                        }
                    }
                    if !sharp[icm] && !sharp[indfin as usize] {
                        ii = icmoins_l;
                        while ii != indfin {
                            isg1[ii as usize] = true;
                            if ii == nedge - 1 {
                                ii = 0;
                            } else {
                                ii += 1;
                            }
                        }
                    }
                    ic += nbface - 1;
                }
                ic += 1;
            }
        }

        // case when the connecting curve between ic and icplus crosses many
        // faces
        // OCCT L2542-2640
        let mut ecom: Vec<Shape> = Vec::new();
        let mut eproj: Vec<Shape> = Vec::new();
        let mut parcom: Vec<f64> = Vec::new();
        if !deuxconges {
            for ic in 0..nedge {
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                let icp = ic as usize;
                let icpl = icplus_l as usize;
                let _ = icmoins_l;
                if !oksea[icp] {
                    let mut iface1 = numfa[icp][icpl];
                    let mut iface2 = numfa[icpl][icp];
                    if !sharp[icp] {
                        iface1 = if jf[icp] == 1 {
                            surf_index(&cd, icp, i_arr[icp][icpl], ChFi3dSurfType::FACE1)
                        } else {
                            surf_index(&cd, icp, i_arr[icp][icpl], ChFi3dSurfType::FACE2)
                        };
                    }
                    let face1 = dstr.shape(iface1).clone();

                    if !sharp[icpl] {
                        iface2 = if jf[icpl] == 1 {
                            surf_index(&cd, icpl, i_arr[icpl][icp], ChFi3dSurfType::FACE2)
                        } else {
                            surf_index(&cd, icpl, i_arr[icpl][icp], ChFi3dSurfType::FACE1)
                        };
                    }
                    let face2 = dstr.shape(iface2).clone();
                    if !face1.is_same(&face2) {
                        if fproj.is_empty() {
                            fproj.push(face1.clone());
                            fproj.push(face2.clone());
                        }
                        moresurf[icp] = true;
                        nbface = fproj.len() as i32;
                        if !fproj[(nbface - 1) as usize].is_same(&face2) {
                            fproj.remove((nbface - 1) as usize);
                            fproj.push(face2.clone());
                        }
                        if !fproj[0].is_same(&face1) {
                            fproj.remove(0);
                            fproj.insert(0, face1.clone());
                        }
                        let mut edge = Shape::null();
                        for nb in 1..=(nbface - 1) {
                            cherche_edge1(
                                &brep,
                                &fproj[(nb - 1) as usize],
                                &fproj[nb as usize],
                                &mut edge,
                            );
                            ecom.push(edge.clone());
                            para = brep_tool_parameter(&brep, &brep.first_vertex(&edge), &edge);
                            parcom.push(para);
                        }
                        let mut pr: Vec<Option<rcad_kernel::geom::Curve2d>> = Vec::new();
                        let mut cr: Vec<Option<rcad_kernel::geom::Curve3>> = Vec::new();
                        curve_hermite(
                            &brep,
                            dstr,
                            &cd[icp].read().expect("stripe lock"),
                            jf[icp],
                            i_arr[icp][icpl],
                            p_arr[icp][icpl],
                            sens[icp],
                            sharp[icp],
                            &evive[icp],
                            &cd[icpl].read().expect("stripe lock"),
                            jf[icpl],
                            i_arr[icpl][icp],
                            p_arr[icpl][icp],
                            sens[icpl],
                            sharp[icpl],
                            &evive[icpl],
                            nbface,
                            &mut ecom,
                            &fproj,
                            &mut pr,
                            &mut cr,
                            &mut eproj,
                            &mut parcom,
                            &mut error,
                        );
                        ecom.push(ecom[(nbface - 2) as usize].clone());
                        parcom.push(parcom[(nbface - 2) as usize]);
                    }
                }
            }
        }

        // case when two fillets have the same commonpoints
        // one continues then by intersection
        // it is checked if the extremities of the intersection coincide
        // with commonpoints
        // OCCT L2642-2703
        let mut intersection = false;
        if nconges == 2 && !deuxconges {
            let mut p1;
            let mut p2v;
            let mut p3;
            let mut p4;
            let mut ic1 = 0i32;
            let mut ic2 = 0i32;
            trouve = false;
            {
                let mut ic = 0i32;
                while ic < nedge && !trouve {
                    if !sharp[ic as usize] {
                        ic1 = ic;
                        trouve = true;
                    }
                    ic += 1;
                }
            }
            for ic in 0..nedge {
                if !sharp[ic as usize] && ic != ic1 {
                    ic2 = ic;
                }
            }
            let mut jfp = 3 - jf[ic1 as usize];
            let mut icplus_l = 0i32;
            let mut icmoins_l = 0i32;
            indices(nedge, ic1, &mut icplus_l, &mut icmoins_l);
            calcul_p2d_on_surf(
                &cd[ic1 as usize].read().expect("stripe lock"),
                jfp,
                i_arr[ic1 as usize][icmoins_l as usize],
                p_arr[ic1 as usize][icmoins_l as usize],
                &mut p2,
            );
            indice = surf_index(
                &cd,
                ic1 as usize,
                i_arr[ic1 as usize][icmoins_l as usize],
                ChFi3dSurfType::ChFiSURFACE,
            );
            p1 = dstr.surface(indice).surface.point_at(p2.x, p2.y);

            calcul_p2d_on_surf(
                &cd[ic1 as usize].read().expect("stripe lock"),
                jf[ic1 as usize],
                i_arr[ic1 as usize][icplus_l as usize],
                p_arr[ic1 as usize][icplus_l as usize],
                &mut p2,
            );
            indice = surf_index(
                &cd,
                ic1 as usize,
                i_arr[ic1 as usize][icplus_l as usize],
                ChFi3dSurfType::ChFiSURFACE,
            );
            p2v = dstr.surface(indice).surface.point_at(p2.x, p2.y);

            jfp = 3 - jf[ic2 as usize];
            indices(nedge, ic2, &mut icplus_l, &mut icmoins_l);
            calcul_p2d_on_surf(
                &cd[ic2 as usize].read().expect("stripe lock"),
                jfp,
                i_arr[ic2 as usize][icmoins_l as usize],
                p_arr[ic2 as usize][icmoins_l as usize],
                &mut p2,
            );
            indice = surf_index(
                &cd,
                ic2 as usize,
                i_arr[ic2 as usize][icmoins_l as usize],
                ChFi3dSurfType::ChFiSURFACE,
            );
            p3 = dstr.surface(indice).surface.point_at(p2.x, p2.y);

            calcul_p2d_on_surf(
                &cd[ic2 as usize].read().expect("stripe lock"),
                jf[ic2 as usize],
                i_arr[ic2 as usize][icplus_l as usize],
                p_arr[ic2 as usize][icplus_l as usize],
                &mut p2,
            );
            indice = surf_index(
                &cd,
                ic2 as usize,
                i_arr[ic2 as usize][icplus_l as usize],
                ChFi3dSurfType::ChFiSURFACE,
            );
            p4 = dstr.surface(indice).surface.point_at(p2.x, p2.y);
            intersection = (p1.distance(p4) <= 1.0e-7 || p1.distance(p3) <= 1.0e-7)
                && (p2v.distance(p4) <= 1.0e-7 || p2v.distance(p3) <= 1.0e-7);
            if intersection {
                let mut introuve = false;
                perform_two_corner_same_ext(
                    &brep,
                    dstr,
                    &cd[ic1 as usize],
                    index_arr[ic1 as usize],
                    sens[ic1 as usize],
                    &cd[ic2 as usize],
                    index_arr[ic2 as usize],
                    sens[ic2 as usize],
                    &mut introuve,
                );
                if introuve {
                    return;
                }
            }
        }

        // declaration for plate
        // GeomPlate_BuildPlateSurface PSurf(3,10,3,tol2d,tolesp,angular);
        //
        // Sence of Plate parameters and their preferable values :
        // degree is total order of ordinary or mixed derivatives:
        // dS/dU, dS/dV have degree 1, d2S/dU2, d2S/dV2, d2S/(dUdV) have
        // degree 2 nbiter - number of iterations, when surface from
        // previous iteration uses as initial surface for next one
        // practically this process does not converge, using "bad" initial
        // surface leads to much more "bad" solution. constr is order of
        // constraint: 0 - G0, 1 - G1 ... Using constraint order > 0 very
        // often causes unpredictable undulations of solution
        // OCCT L2705-2717
        let degree = 3i32;
        let nbcurvpnt = 10i32;
        let nbiter = 1i32;
        let constr = 1i32; // G1
        let mut psurf =
            GeomPlateBuildPlateSurface::new(degree, nbcurvpnt, nbiter, tol2d, tolapp3d, angular);
        // calculation of curves on surface for each stripe
        // OCCT L2719-2851
        for ic in 0..nedge {
            let mut p2d1 = DVec2::ZERO;
            let mut p2d2 = DVec2::ZERO;
            if !sharp[ic as usize] {
                n3d += 1;
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                let icp = ic as usize;
                let icm = icmoins_l as usize;
                let icpl = icplus_l as usize;
                let jfp = 3 - jf[icp];
                calcul_p2d_on_surf(
                    &cd[icp].read().expect("stripe lock"),
                    jfp,
                    i_arr[icp][icm],
                    p_arr[icp][icm],
                    &mut p2d1,
                );
                calcul_p2d_on_surf(
                    &cd[icp].read().expect("stripe lock"),
                    jf[icp],
                    i_arr[icp][icpl],
                    p_arr[icp][icpl],
                    &mut p2d2,
                );
                //      if (i[ic][icplus]!=  i[ic][icmoins]) std::cout<<"probleme surface"<<std::endl;
                indice = surf_index(&cd, icp, i_arr[icp][icpl], ChFi3dSurfType::ChFiSURFACE);
                let surf_payload = dstr.surface(indice).surface.clone();
                let asurf = GeomAdaptorSurface::new(surf_payload.clone());
                // calculation of curve 2d
                xdir = p2d2.x - p2d1.x;
                ydir = p2d2.y - p2d1.y;
                let l0 = (xdir * xdir + ydir * ydir).sqrt();
                let dir_len = if l0 > 0.0 { l0 } else { 1.0 };
                let l = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                    origin: p2d1,
                    direction: DVec2::new(xdir / dir_len, ydir / dir_len),
                });
                let pcurve = rcad_kernel::geom::Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                    curve: Box::new(l),
                    t_min: 0.0,
                    t_max: l0,
                });
                // OCCT: Acurv = new Geom2dAdaptor_Curve(pcurve);
                //       CurvOnS(Acurv, Asurf); HCons = new
                //       Adaptor3d_CurveOnSurface(CurvOnS);
                let acurv = Geom2dAdaptorCurve::load(pcurve.clone(), 0.0, l0);
                let hcons = Adaptor3dCurveOnSurface::new(acurv, asurf);
                // Order.SetValue(ic,1);
                order[icp] = constr;
                let cont = GeomPlateCurveConstraint::new(hcons, order[icp], nbcurvpnt, tolapp3d, angular, 0.1);
                psurf.add(cont);

                // calculate indexes of points and of the curve for the DS
                isfirst = sens[icp] == 1;
                // OCCT: GeomLib::BuildCurve3d(tolapp, CurvOnS, fp, lp,
                //       Curv3d, maxapp, avedev); — CurvOnS rebuilt from the
                // same pcurve/Asurf pair (OCCT shares the handle).
                let curv_on_s = {
                    let acurv = Geom2dAdaptorCurve::load(pcurve.clone(), 0.0, l0);
                    Adaptor3dCurveOnSurface::new(acurv, GeomAdaptorSurface::new(surf_payload.clone()))
                };
                geom_lib_build_curve3d(
                    tolapp,
                    &curv_on_s,
                    curv_on_s.first_parameter(),
                    curv_on_s.last_parameter(),
                    &mut curv3d,
                    &mut maxapp,
                    &mut avedev,
                );
                let tcurv3d = TopOpeBRepDSCurve::new(curv3d.clone(), maxapp);
                indcurve3d[n3d as usize] = dstr.add_curve(tcurv3d);
                // OCCT: point1 = CurvOnS.Value(fp); point2 =
                //       CurvOnS.Value(lp).
                let point1 = curv_on_s.value(curv_on_s.first_parameter());
                let point2 = curv_on_s.value(curv_on_s.last_parameter());

                let tpoint1 = TopOpeBRepDSPoint::new(point1, maxapp);
                let tpoint2 = TopOpeBRepDSPoint::new(point2, maxapp);
                errapp[icp] = maxapp;
                if ic == 0 {
                    // it is necessary to create two points
                    indpoint[icp][0] = dstr.add_point(tpoint1);
                    indpoint[icp][1] = dstr.add_point(tpoint2);
                } else {
                    // probably the points are already on the fillet
                    // (previous intersection...)
                    trouve = false;
                    let mut ii_found = 0usize;
                    for ii in 0..ic as usize {
                        if trouve {
                            break;
                        }
                        if !sharp[ii] {
                            let tpt = dstr.point(indpoint[ii][1]);
                            if point1.distance(tpt.point()) < 1.0e-4 {
                                trouve = true;
                                ii_found = ii;
                            }
                        }
                    }
                    if trouve {
                        indpoint[icp][0] = indpoint[ii_found][1];
                    } else {
                        indpoint[icp][0] = dstr.add_point(tpoint1);
                    }

                    trouve = false;
                    for ii in 0..ic as usize {
                        if trouve {
                            break;
                        }
                        if !sharp[ii] {
                            let tpt = dstr.point(indpoint[ii][0]);
                            if point2.distance(tpt.point()) < 1.0e-4 {
                                trouve = true;
                                ii_found = ii;
                            }
                        }
                    }
                    if trouve {
                        indpoint[icp][1] = indpoint[ii_found][0];
                    } else {
                        indpoint[icp][1] = dstr.add_point(tpoint2);
                    }
                }

                //   update of the stripe
                // OCCT L2820-2849
                let isurf1 = 3 - jf[icp];
                let isurf2 = jf[icp];
                let (fp, lp) = match &pcurve {
                    rcad_kernel::geom::Curve2d::Trimmed(t) => (t.t_min, t.t_max),
                    _ => (0.0, 0.0),
                };
                {
                    let mut st = cd[icp].write().expect("stripe lock");
                    if isurf1 == 2 {
                        st.set_orientation(Orientation::Reversed, isfirst);
                    }
                    st.set_curve(indcurve3d[n3d as usize], isfirst);
                    st.set_index_point(indpoint[icp][0], isfirst, isurf1);
                    st.set_index_point(indpoint[icp][1], isfirst, isurf2);
                    st.set_parameters(isfirst, fp, lp);
                    let sd = st.set_of_surf_data()[(i_arr[icp][icm] - 1) as usize].clone();
                    let mut fd = sd.write().expect("surfdata lock");
                    let mut cp1 = super::chfi_ds::ChFiDS_CommonPoint::default();
                    let mut cp2 = super::chfi_ds::ChFiDS_CommonPoint::default();
                    cp1.set_point(point1);
                    cp2.set_point(point2);
                    *fd.change_vertex(isfirst, isurf1) = cp1;
                    *fd.change_vertex(isfirst, isurf2) = cp2;
                    fd.change_interference(isurf1)
                        .set_parameter(isfirst, p_arr[icp][icm]);
                    fd.change_interference(isurf2)
                        .set_parameter(isfirst, p_arr[icp][icpl]);
                    st.change_pcurve(isfirst, pcurve.clone());
                }
            }
        }

        // calculate the indices of points for living edges
        // OCCT L2853-2905
        for ic in 0..nedge {
            if sharp[ic as usize] {
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                let icp = ic as usize;
                let icpl = icplus_l as usize;
                let icm = icmoins_l as usize;
                let c = BRepAdaptorCurve::initialize(&brep, &evive[icp]);
                let pe = c.value(p_arr[icp][icpl]);
                let edge_tol = evive[icp].as_edge().map(|e| e.tolerance).unwrap_or(0.0);
                let tpe = TopOpeBRepDSPoint::new(pe, edge_tol);
                if deuxconges {
                    ivtx = dstr.add_shape(&v1);
                    indpoint[icp][0] = ivtx;
                    indpoint[icp][1] = ivtx;
                }
                if !sharp[icpl] {
                    isfirst = sens[icpl] == 1;
                    let jfp = 3 - jf[icpl];
                    let st = cd[icpl].read().expect("stripe lock");
                    let fd = st.set_of_surf_data()[(i_arr[icpl][icp] - 1) as usize].clone();
                    let cp = cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, jfp);
                    if cp.point().distance(pe) <= 1.0e-4f64.max(tolcp) {
                        // edge was limited by the 1st CommonPoint of CD[icplus]
                        indpoint[icp][0] = indpoint[icpl][0];
                        indpoint[icp][1] = indpoint[icpl][0];
                    }
                }
                if !sharp[icm] {
                    isfirst = sens[icm] == 1;
                    let st = cd[icm].read().expect("stripe lock");
                    let fd = st.set_of_surf_data()[(i_arr[icm][icp] - 1) as usize].clone();
                    let cp = cp_change_vertex(&fd.read().expect("surfdata lock"), isfirst, jf[icm]);
                    if cp.point().distance(pe) <= 1.0e-4f64.max(tolcp) {
                        // edge was limited by the 2nd CommonPoint of CD[icmoins]
                        if indpoint[icp][0] == 0 {
                            indpoint[icp][0] = indpoint[icm][1];
                            indpoint[icp][1] = indpoint[icm][1];
                        }
                    }
                }
                if indpoint[icp][0] == 0 {
                    indpoint[icp][0] = dstr.add_point(tpe);
                    indpoint[icp][1] = indpoint[icp][0];
                }
            }
        }

        // calculation of intermediary curves connecting two stripes in case
        // if there is no intersection. The curve is a straight line,
        // projection or batten
        // OCCT L2907-3260
        let mut raccordbatten;
        if !inters {
            for ic in 0..nedge {
                if !oksea[ic as usize] && !moresurf[ic as usize] && !libre[ic as usize] {
                    let mut icplus_l = 0i32;
                    let mut icmoins_l = 0i32;
                    indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                    let icp = ic as usize;
                    let icpl = icplus_l as usize;
                    let icm = icmoins_l as usize;
                    raccordbatten = false;
                    if !regul[icp] {
                        raccordbatten = true;
                        if regul[icpl] {
                            raccordbatten = false;
                        }
                    }
                    n3d += 1;
                    let mut curv2d1: Option<rcad_kernel::geom::Curve2d> = None;
                    let mut curv2d2: Option<rcad_kernel::geom::Curve2d> = None;
                    let mut pcurve: Option<rcad_kernel::geom::Curve2d> = None;
                    let mut curveint: Option<rcad_kernel::geom::Curve3> = None;

                    // return the 1st curve 2d
                    // and the 1st connection point
                    // OCCT L2937-2955
                    if sharp[icp] {
                        if let Some((c, _, _)) = brep.curve_on_surface(&evive[icp], &fvive[icp][icpl]) {
                            curv2d1 = Some(c);
                        }
                    } else {
                        calcul_c2d_on_face(
                            &cd[icp].read().expect("stripe lock"),
                            jf[icp],
                            i_arr[icp][icpl],
                            &mut curv2d1,
                        );
                    }

                    if curv2d1.is_none() {
                        panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                    }
                    let mut p2d1 = curv2d1
                        .as_ref()
                        .expect("curv2d1")
                        .point_at(p_arr[icp][icpl]);

                    // recuperation de la deuxieme courbe 2d
                    // et du deuxieme point de raccordement
                    // OCCT L2957-2975
                    if sharp[icpl] {
                        if let Some((c, _, _)) = brep.curve_on_surface(&evive[icpl], &fvive[icp][icpl]) {
                            curv2d2 = Some(c);
                        }
                    } else {
                        let jfp = 3 - jf[icpl];
                        calcul_c2d_on_face(
                            &cd[icpl].read().expect("stripe lock"),
                            jfp,
                            i_arr[icpl][icp],
                            &mut curv2d2,
                        );
                    }
                    if curv2d2.is_none() {
                        panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                    }
                    let mut p2d2 = curv2d2
                        .as_ref()
                        .expect("curv2d2")
                        .point_at(p_arr[icpl][icp]);

                    let face_surf = brep
                        .face_surface(&fvive[icp][icpl])
                        .cloned()
                        .expect("face surface");
                    let asurf = GeomAdaptorSurface::new(face_surf.clone());
                    let tolu = asurf.u_resolution(1.0e-3);
                    let tolv = asurf.v_resolution(1.0e-3);
                    let ratio = if tolu > tolv { tolu / tolv } else { tolv / tolu };

                    // in case of a sewing edge the parameters are reframed
                    // OCCT L2990-3024
                    if couture {
                        let mut xx;
                        let pi1 = 0.0 <= p2d1.x && p2d1.x <= pi;
                        let pi2 = 0.0 <= p2d2.x && p2d2.x <= pi;

                        if evive[icp].is_same(&edgecouture) {
                            xx = p2d1.x;
                            if pi2 && !pi1 {
                                xx -= 2.0 * pi;
                            }
                            if !pi2 && pi1 {
                                xx += 2.0 * pi;
                            }
                            p2d1.x = xx;
                        }
                        if evive[icpl].is_same(&edgecouture) {
                            xx = p2d2.x;
                            if pi2 && !pi1 {
                                xx += 2.0 * pi;
                            }
                            if !pi2 && pi1 {
                                xx -= 2.0 * pi;
                            }
                            p2d2.x = xx;
                        }
                    }
                    xdir = p2d2.x - p2d1.x;
                    ydir = p2d2.y - p2d1.y;

                    let l0 = (xdir * xdir + ydir * ydir).sqrt();
                    if l0 < 1.0e-7 || ponctuel[icp] {
                        // unused connection
                        // OCCT L3029-3047
                        n3d -= 1;
                        ponctuel[icp] = true;
                        if !deuxconges {
                            if sharp[icpl] && indpoint[icpl][0] == 0 {
                                indpoint[icpl][0] = indpoint[icp][1];
                                indpoint[icpl][1] = indpoint[icp][1];
                            }
                            if sharp[icp] && indpoint[icp][0] == 0 {
                                indpoint[icp][0] = indpoint[icm][1];
                                indpoint[icp][1] = indpoint[icm][1];
                            }
                        }
                    } else {
                        // the connection is a straight line, projection or
                        // batten
                        // OCCT L3048-3139
                        if ratio > 10.0 && nconges == 1 {
                            raccordbatten = true;
                        }
                        if ratio > 10.0 && raccordbatten {
                            calcul_droite(p2d1, xdir, ydir, &mut pcurve);
                            raccordbatten = false;
                        } else if !raccordbatten {
                            // the projected curves are returned
                            if regul[icp] {
                                if cproj2[icp].is_none() {
                                    raccordbatten = true;
                                } else {
                                    pcurve = proj2d2[icp].clone();
                                    curveint = cproj2[icp].clone();
                                    maxapp1 = 1.0e-6;
                                }
                            } else if cproj1[ic as usize + 1].is_none() {
                                raccordbatten = true;
                            } else {
                                pcurve = proj2d1[ic as usize + 1].clone();
                                curveint = cproj1[ic as usize + 1].clone();
                                maxapp1 = 1.0e-6;
                            }
                        }
                        let mut contraint1 = true;
                        let mut contraint2 = true;
                        if raccordbatten {
                            let inverseic;
                            let inverseicplus;
                            if sharp[icp] {
                                inverseic = brep.first_vertex(&evive[icp]).is_same(&v1);
                            } else {
                                inverseic = sens[icp] == 1;
                            }
                            if sharp[icpl] {
                                inverseicplus = brep.first_vertex(&evive[icpl]).is_same(&v1);
                            } else {
                                inverseicplus = sens[icpl] == 1;
                            }
                            if evive[icp].is_same(&edgelibre1) || evive[icp].is_same(&edgelibre2) {
                                contraint1 = false;
                            }
                            if evive[icpl].is_same(&edgelibre1) || evive[icpl].is_same(&edgelibre2) {
                                contraint2 = false;
                            }
                            calcul_batten(
                                &asurf,
                                &fvive[icp][icpl],
                                &brep,
                                xdir,
                                ydir,
                                p2d1,
                                p2d2,
                                contraint1,
                                contraint2,
                                &curv2d1,
                                &curv2d2,
                                p_arr[icp][icpl],
                                p_arr[icpl][icp],
                                inverseic,
                                inverseicplus,
                                &mut pcurve,
                            );
                        }

                        // construction of borders for Plate
                        // OCCT L3141-3144
                        let pcurve_v = pcurve.clone().expect("pcurve");
                        let (pcf, pcl) = match &pcurve_v {
                            rcad_kernel::geom::Curve2d::Trimmed(t) => (t.t_min, t.t_max),
                            _ => (0.0, 0.0),
                        };
                        let acurv = Geom2dAdaptorCurve::load(pcurve_v.clone(), pcf, pcl);
                        let curv_on_s = Adaptor3dCurveOnSurface::new(acurv, asurf);

                        // constraints G1 are set if edges ic and icplus are
                        // not both alive
                        // OCCT L3146-3164
                        order[n3d as usize] = 0;
                        if !sharp[icp] && !sharp[icpl] {
                            order[n3d as usize] = 1;
                        }
                        if !contraint1 && !sharp[icpl] {
                            order[n3d as usize] = 1;
                        }
                        if !contraint2 && !sharp[icp] {
                            order[n3d as usize] = 1;
                        }
                        if tangentregul[icp] || tangentregul[icpl] {
                            order[n3d as usize] = 1;
                        }
                        if isg1[icp] {
                            order[n3d as usize] = 1;
                        }
                        let cont = GeomPlateCurveConstraint::new(
                            curv_on_s,
                            order[n3d as usize],
                            10,
                            tolapp3d,
                            angular,
                            0.1,
                        );
                        psurf.add(cont);

                        // calculation of curve 3d if it is not a projection
                        // OCCT L3173-3186
                        if curveint.is_none() {
                            let curv_on_s2 = {
                                let acurv = Geom2dAdaptorCurve::load(pcurve_v.clone(), pcf, pcl);
                                Adaptor3dCurveOnSurface::new(acurv, GeomAdaptorSurface::new(face_surf.clone()))
                            };
                            geom_lib_build_curve3d(
                                tolapp,
                                &curv_on_s2,
                                curv_on_s2.first_parameter(),
                                curv_on_s2.last_parameter(),
                                &mut curv3d,
                                &mut maxapp1,
                                &mut avedev,
                            );
                            pardeb = pcf;
                            parfin = pcl;
                            if let Some(c3) = curv3d.clone() {
                                curveint = Some(rcad_kernel::geom::Curve3::Trimmed(
                                    rcad_kernel::geom::TrimmedCurve3 {
                                        curve: Box::new(c3),
                                        first: pardeb,
                                        last: parfin,
                                    },
                                ));
                            }
                        }

                        // storage in the DS
                        // OCCT L3188-3230
                        let tcurv3d = TopOpeBRepDSCurve::new(curveint.clone(), maxapp1);
                        indcurve3d[n3d as usize] = dstr.add_curve(tcurv3d);
                        pardeb = curveint
                            .as_ref()
                            .map(|c| c.default_domain()[0])
                            .unwrap_or(0.0);
                        parfin = curveint
                            .as_ref()
                            .map(|c| c.default_domain()[1])
                            .unwrap_or(0.0);
                        if sharp[icpl] && indpoint[icpl][0] == 0 {
                            // it is necessary to initialize
                            // indpoint[icplus][0] and indpoint[icplus][1]
                            let point2 = curveint
                                .as_ref()
                                .map(|c| c.point_at(parfin))
                                .unwrap_or(DVec3::ZERO);
                            let tpoint2 = TopOpeBRepDSPoint::new(point2, maxapp);
                            indpoint[icpl][0] = dstr.add_point(tpoint2);
                            indpoint[icpl][1] = indpoint[icpl][0];
                        }
                        let mut isvt1 = false;
                        let mut isvt2 = false;
                        if deuxconges {
                            isvt1 = sharp[icp];
                            isvt2 = sharp[icpl];
                        }
                        let interfp1 = chfi3d_fil_point_in_ds(
                            Orientation::Forward,
                            indcurve3d[n3d as usize],
                            indpoint[icp][1],
                            pardeb,
                            isvt1,
                        );
                        let interfp2 = chfi3d_fil_point_in_ds(
                            Orientation::Reversed,
                            indcurve3d[n3d as usize],
                            indpoint[icpl][0],
                            parfin,
                            isvt2,
                        );
                        dstr
                            .change_curve_interferences(indcurve3d[n3d as usize])
                            .push(interfp1);
                        dstr
                            .change_curve_interferences(indcurve3d[n3d as usize])
                            .push(interfp2);
                        if !isvt1 {
                            let tpt1_tol = dstr.point(indpoint[icp][1]).tolerance() + maxapp1;
                            dstr.change_point(indpoint[icp][1]).set_tolerance(tpt1_tol);
                        }
                        if !isvt2 {
                            let tpt2_tol = dstr.point(indpoint[icpl][0]).tolerance() + maxapp1;
                            dstr.change_point(indpoint[icpl][0]).set_tolerance(tpt2_tol);
                        }

                        // calculate orientation of the curve
                        // OCCT L3232-3253
                        let mut orinterf = Orientation::Forward;
                        if !sharp[icp] {
                            orientation_ic_non_vive(
                                &cd[icp].read().expect("stripe lock"),
                                jf[icp],
                                i_arr[icp][icpl],
                                sens[icp],
                                &mut orinterf,
                            );
                        } else if !sharp[icpl] {
                            orientation_icplus_non_vive(
                                &cd[icpl].read().expect("stripe lock"),
                                jf[icpl],
                                i_arr[icpl][icp],
                                sens[icpl],
                                &mut orinterf,
                            );
                        } else {
                            orientation_arete_vive_consecutive(
                                &brep,
                                &fvive[icp][icpl],
                                &evive[icp],
                                &v1,
                                &mut orinterf,
                            );
                        }
                        let interfc = chfi3d_fil_curve_in_ds(
                            indcurve3d[n3d as usize],
                            numfa[icp][icpl],
                            pcurve.clone(),
                            orinterf,
                        );
                        dstr.change_shape_interferences(numfa[icp][icpl]).push(interfc);
                    }
                } // end of processing by edge
            } // end of the loop on edges
        } // end of processing for intermediary curves

        //  storage in the DS of curves projected on several faces
        // OCCT L3262-3405
        for ic in 0..nedge {
            if moresurf[ic as usize] {
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                let icp = ic as usize;
                let icpl = icplus_l as usize;
                let _ = icmoins_l;
                let mut ind = 0i32; // must be initialized because of possible use, see L2249
                let mut orvt;
                let mut oredge = Orientation::Forward;
                let mut indpoint1;
                let mut indpoint2;
                // OCCT L3277-3280 — tpt1/tpt2 tolerance updates.  Note the
                // OCCT L3280 reads tpt1 (not tpt2) for the second update.
                let tpt1_tol = dstr.point(indpoint[icp][1]).tolerance() + error;
                dstr.change_point(indpoint[icp][1]).set_tolerance(tpt1_tol);
                let tpt2_tol = dstr.point(indpoint[icp][1]).tolerance() + error;
                dstr.change_point(indpoint[icpl][0]).set_tolerance(tpt2_tol);
                for nb in 1..=nbface {
                    orvt = Orientation::Reversed;
                    let ecom_nb = ecom[(nb - 1) as usize].clone();
                    let vf = brep.first_vertex(&ecom_nb);
                    let vl = brep.last_vertex(&ecom_nb);
                    let pf = brep.vertex_position(&vf);
                    let pl = brep.vertex_position(&vl);
                    para = parcom[(nb - 1) as usize];
                    let pcom = brep
                        .edge_curve_world(&ecom_nb)
                        .map(|(c, _)| c.point_at(para))
                        .unwrap_or(DVec3::ZERO);
                    if pf.distance(brep.vertex_position(&v1)) < pl.distance(brep.vertex_position(&v1)) {
                        orvt = Orientation::Forward;
                    }
                    if !eproj[(nb - 1) as usize].is_null() {
                        n3d += 1;
                        let (proj, up1, up2) = match brep.curve_on_surface(
                            &eproj[(nb - 1) as usize],
                            &fproj[(nb - 1) as usize],
                        ) {
                            Some(v) => v,
                            None => {
                                panic!("Standard_ConstructionError: Failed to get p-curve of edge")
                            }
                        };
                        let proj2d =
                            rcad_kernel::geom::Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                                curve: Box::new(proj),
                                t_min: up1,
                                t_max: up2,
                            });
                        let (projc, up1c, up2c) = match brep.edge_curve_world(&eproj[(nb - 1) as usize]) {
                            Some((c, r)) => (c, r[0], r[1]),
                            None => panic!(
                                "Standard_ConstructionError: Failed to get 3D curve of edge"
                            ),
                        };
                        let cproj =
                            rcad_kernel::geom::Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
                                curve: Box::new(projc),
                                first: up1c,
                                last: up2c,
                            });
                        pardeb = cproj.default_domain()[0];
                        parfin = cproj.default_domain()[1];
                        let p1 = cproj.point_at(pardeb);
                        let p2v = cproj.point_at(parfin);
                        if p1.distance(dstr.point(indpoint[icp][1]).point()) < 1.0e-3 {
                            indpoint1 = indpoint[icp][1];
                        } else {
                            indpoint1 = ind;
                        }
                        if p2v.distance(dstr.point(indpoint[icpl][0]).point()) < 1.0e-3 {
                            indpoint2 = indpoint[icpl][0];
                        } else {
                            let tpoint2 = TopOpeBRepDSPoint::new(p2v, error);
                            indpoint2 = dstr.add_point(tpoint2);
                            ind = indpoint2;
                        }
                        let asurf = GeomAdaptorSurface::new(
                            brep.face_surface(&fproj[(nb - 1) as usize])
                                .cloned()
                                .expect("face surface"),
                        );
                        let acurv = Geom2dAdaptorCurve::load(proj2d.clone(), up1, up2);
                        let curv_on_s = Adaptor3dCurveOnSurface::new(acurv, asurf);
                        order[n3d as usize] = 1;
                        let cont = GeomPlateCurveConstraint::new(
                            curv_on_s,
                            order[n3d as usize],
                            10,
                            tolapp3d,
                            angular,
                            0.1,
                        );
                        psurf.add(cont);
                        let tcurv3d = TopOpeBRepDSCurve::new(Some(cproj.clone()), error);
                        indcurve3d[n3d as usize] = dstr.add_curve(tcurv3d);
                        let interfp1 = chfi3d_fil_point_in_ds(
                            Orientation::Forward,
                            indcurve3d[n3d as usize],
                            indpoint1,
                            pardeb,
                            false,
                        );
                        let interfp2 = chfi3d_fil_point_in_ds(
                            Orientation::Reversed,
                            indcurve3d[n3d as usize],
                            indpoint2,
                            parfin,
                            false,
                        );
                        dstr
                            .change_curve_interferences(indcurve3d[n3d as usize])
                            .push(interfp1);
                        dstr
                            .change_curve_interferences(indcurve3d[n3d as usize])
                            .push(interfp2);
                        num = dstr.add_shape(&fproj[(nb - 1) as usize]);
                        // OCCT L3350-3358 — explore the face's FORWARD
                        // oriented edges for Ecom.
                        let mut f_fwd = fproj[(nb - 1) as usize].clone();
                        f_fwd.orientation = Orientation::Forward;
                        for we in topexp_face_edges(&brep, &f_fwd) {
                            if ecom_nb.is_same(&we) {
                                oredge = we.orientation;
                                break;
                            }
                        }

                        // calculation of the orientation
                        // OCCT L3360-3383
                        let orinterf;
                        if p1.distance(pcom) > 1.0e-4 {
                            if orvt == Orientation::Forward {
                                orinterf = oredge;
                            } else {
                                orinterf = topabs_reverse(oredge);
                            }
                        } else if orvt == Orientation::Forward {
                            orinterf = topabs_reverse(oredge);
                        } else {
                            orinterf = oredge;
                        }
                        let interfc = chfi3d_fil_curve_in_ds(
                            indcurve3d[n3d as usize],
                            num,
                            Some(proj2d.clone()),
                            orinterf,
                        );
                        dstr.change_shape_interferences(num).push(interfc);
                    }
                    indice = ind;
                    if nb != nbface {
                        if eproj[(nb - 1) as usize].is_null() {
                            indice = indpoint[icp][1];
                        }
                        if eproj[nb as usize].is_null() {
                            indice = indpoint[icpl][0];
                        }
                        indice_arr[n3d as usize] = indice;
                        let iarc1 = dstr.add_shape(&ecom_nb);
                        let interfp1 =
                            chfi3d_fil_point_in_ds(orvt, iarc1, indice, parcom[(nb - 1) as usize], false);
                        dstr.change_shape_interferences(iarc1).push(interfp1);
                    }
                }
            }
        }

        // case when two free borders are tangent
        // OCCT L3407-3488
        if droit {
            for ic in 0..nedge {
                let mut icplus_l = 0i32;
                let mut icmoins_l = 0i32;
                indices(nedge, ic, &mut icplus_l, &mut icmoins_l);
                let icp = ic as usize;
                let icpl = icplus_l as usize;
                let icm = icmoins_l as usize;
                let mut indpoint1;
                let mut indpoint2;
                let mut isvt1 = false;
                let mut isvt2 = false;
                let ecur2 = evive[icp].clone();
                if ecur2.is_same(&edgelibre1) || ecur2.is_same(&edgelibre2) {
                    n3d += 1;
                    let (curve2d, ufirst0, ulast0) =
                        match brep.curve_on_surface(&evive[icp], &fvive[icp][icpl]) {
                            Some(v) => v,
                            None => {
                                panic!("Standard_ConstructionError: Failed to get p-curve of edge")
                            }
                        };
                    let (curve, ufirst1, ulast1) = match brep.edge_curve_world(&evive[icp]) {
                        Some((c, r)) => (c, r[0], r[1]),
                        None => {
                            panic!("Standard_ConstructionError: Failed to get 3D curve of edge")
                        }
                    };
                    let mut ufirst = ufirst0;
                    let mut ulast = ulast0;
                    let _ = (ufirst1, ulast1);
                    let mut ctrim: Option<rcad_kernel::geom::Curve3>;
                    let mut ctrim2d: Option<rcad_kernel::geom::Curve2d>;
                    if brep.first_vertex(&evive[icp]).is_same(&v1) {
                        ctrim = Some(rcad_kernel::geom::Curve3::Trimmed(
                            rcad_kernel::geom::TrimmedCurve3 {
                                curve: Box::new(curve.clone()),
                                first: ufirst,
                                last: p_arr[icp][icm],
                            },
                        ));
                        ctrim2d = Some(rcad_kernel::geom::Curve2d::Trimmed(
                            rcad_kernel::geom::TrimmedCurve2 {
                                curve: Box::new(curve2d.clone()),
                                t_min: ufirst,
                                t_max: p_arr[icp][icm],
                            },
                        ));
                        indpoint1 = dstr.add_shape(&v1);
                        isvt1 = true;
                        indpoint2 = indpoint[icp][1];
                    } else {
                        ctrim = Some(rcad_kernel::geom::Curve3::Trimmed(
                            rcad_kernel::geom::TrimmedCurve3 {
                                curve: Box::new(curve.clone()),
                                first: p_arr[icp][icm],
                                last: ulast,
                            },
                        ));
                        ctrim2d = Some(rcad_kernel::geom::Curve2d::Trimmed(
                            rcad_kernel::geom::TrimmedCurve2 {
                                curve: Box::new(curve2d.clone()),
                                t_min: p_arr[icp][icm],
                                t_max: ulast,
                            },
                        ));
                        indpoint2 = dstr.add_shape(&v1);
                        isvt2 = true;
                        indpoint1 = indpoint[icp][1];
                    }
                    if libre[icp] {
                        if brep.first_vertex(&evive[icp]).is_same(&v1) {
                            // OCCT: ctrim->Reverse(); ctrim2d->Reverse();
                            ctrim = ctrim.take().map(|c| reverse_trimmed3(&c));
                            ctrim2d = ctrim2d.take().map(|c| reverse_trimmed2(&c));
                            indpoint2 = dstr.add_shape(&v1);
                            isvt2 = true;
                            isvt1 = false;
                            indpoint1 = indpoint[icp][1];
                        }
                    } else if brep.last_vertex(&evive[icp]).is_same(&v1) {
                        ctrim = ctrim.take().map(|c| reverse_trimmed3(&c));
                        ctrim2d = ctrim2d.take().map(|c| reverse_trimmed2(&c));
                        indpoint1 = dstr.add_shape(&v1);
                        isvt1 = true;
                        isvt2 = false;
                        indpoint2 = indpoint[icp][1];
                    }
                    let ctrim_v = ctrim.take().expect("ctrim");
                    let ctrim2d_v = ctrim2d.take().expect("ctrim2d");
                    ufirst = ctrim_v.default_domain()[0];
                    ulast = ctrim_v.default_domain()[1];
                    let asurf = GeomAdaptorSurface::new(
                        brep.face_surface(&fvive[icp][icpl])
                            .cloned()
                            .expect("face surface"),
                    );
                    let acurv = Geom2dAdaptorCurve::load(ctrim2d_v.clone(), ufirst, ulast);
                    let curv_on_s = Adaptor3dCurveOnSurface::new(acurv, asurf);
                    order[n3d as usize] = 0;
                    let cont = GeomPlateCurveConstraint::new(
                        curv_on_s,
                        order[n3d as usize],
                        10,
                        tolapp3d,
                        angular,
                        0.1,
                    );
                    psurf.add(cont);
                    let tcurv3d = TopOpeBRepDSCurve::new(Some(ctrim_v.clone()), 1.0e-4);
                    indcurve3d[n3d as usize] = dstr.add_curve(tcurv3d);
                    let interfp1 = chfi3d_fil_point_in_ds(
                        Orientation::Forward,
                        indcurve3d[n3d as usize],
                        indpoint1,
                        ufirst,
                        isvt1,
                    );
                    let interfp2 = chfi3d_fil_point_in_ds(
                        Orientation::Reversed,
                        indcurve3d[n3d as usize],
                        indpoint2,
                        ulast,
                        isvt2,
                    );
                    dstr
                        .change_curve_interferences(indcurve3d[n3d as usize])
                        .push(interfp1);
                    dstr
                        .change_curve_interferences(indcurve3d[n3d as usize])
                        .push(interfp2);
                }
            }
        }

        // OCCT L3494 — PSurf.Perform();
        psurf.perform();

        // OCCT L3505 — if (PSurf.IsDone())
        if psurf.is_done() {
            // OCCT L3506-3892 — the plate branch.
            let nbcarreau = 9i32;
            let degmax = 8i32;
            let gp_plate = psurf.surface().expect("pending GeomPlate Surface");

            let mut s2d: Vec<DVec2> = Vec::new();
            let mut s3d: Vec<DVec3> = Vec::new();
            psurf.disc2d_contour(4, &mut s2d);
            psurf.disc3d_contour(4, 0, &mut s3d);
            let seuil = tolapp.max(10.0 * psurf.g0_error());
            let critere = GeomPlatePlateG0Criterion;
            let mapp = GeomPlateMakeApprox::new(&gp_plate, &critere, tolapp, nbcarreau, degmax);
            let surf = mapp.surface();
            let coef = 1.1f64;
            let mut apperror = mapp.criterion_error() * coef;

            //  Storage of the surface plate and corresponding curves in the
            //  DS
            // OCCT L3541-3542
            let tsurf = TopOpeBRepDSSurface::new(
                surf.clone().expect("pending GeomPlate MakeApprox surface"),
                mapp.approx_error(),
            );
            let isurf = dstr.add_surface(tsurf);
            // lbo : historique QDF.
            // OCCT L3544-3549
            let evi_key = v1.ptr_id();
            if !self.my_evi_map.contains_key(&evi_key) {
                self.my_evi_map.insert(evi_key, Vec::new());
            }
            self.my_evi_map.get_mut(&evi_key).expect("evi map").push(isurf);

            // OCCT L3551 — SolInd = CD.Value(0)->SolidIndex();
            let sol_ind = cd[0].read().expect("stripe lock").solid_index();

            // in case when one rereads at top, it is necessary that
            // alive edges that arrive at the top should be removed from the
            // DS. For this they are stored in the DS with their inverted
            // orientation
            // OCCT L3559-3593
            if deuxconges {
                for ic in 0..nedge {
                    if !sharp[ic as usize] {
                        let st = cd[ic as usize].read().expect("stripe lock");
                        let spine = st.spine().expect("spine");
                        let nbedge = spine.base().nb_edges();
                        let arcspine = if sens[ic as usize] == 1 {
                            spine.base().edges(1).clone()
                        } else {
                            spine.base().edges(nbedge).clone()
                        };
                        drop(st);
                        let iarcspine = dstr.add_shape(&arcspine);
                        // OCCT L3578-3585 — explore the FORWARD-oriented
                        // arcspine vertices for V1.
                        let mut ovtx = Orientation::Forward;
                        let (vfirst, vlast) = (
                            brep.first_vertex(&arcspine),
                            brep.last_vertex(&arcspine),
                        );
                        if v1.is_same(&vfirst) {
                            ovtx = Orientation::Forward;
                        } else if v1.is_same(&vlast) {
                            ovtx = Orientation::Reversed;
                        }
                        ovtx = topabs_reverse(ovtx);
                        let parvtx = brep_tool_parameter(&brep, &v1, &arcspine);
                        let interfv = chfi3d_fil_vertex_in_ds(ovtx, iarcspine, ivtx, parvtx);
                        dstr.change_shape_interferences(iarcspine).push(interfv);
                    }
                }
            }

            // calculate orientation of Plate orplate corresponding to
            // surfdata calculation corresponding to the first stripe
            // OCCT L3595-3617
            let mut icplus_l = 0i32;
            let mut icmoins_l = 0i32;
            indices(nedge, 0, &mut icplus_l, &mut icmoins_l);
            isfirst = sens[0] == 1;
            let fd = cd[0]
                .read()
                .expect("stripe lock")
                .set_of_surf_data()[(i_arr[0][icplus_l as usize] - 1) as usize]
                .clone();
            indice = fd.read().expect("surfdata lock").surf();
            let orsurfdata = fd.read().expect("surfdata lock").orientation();
            let surf_ref = surf.as_ref().expect("pending GeomPlate surface");
            let curves2d_vec: Vec<Option<rcad_kernel::geom::Curve2d>> =
                (1..=n3d).map(|k| psurf.curves2d_value(k)).collect();
            let orplate = plate_orientation(surf_ref, &curves2d_vec, sum_face_normal_at_v1);

            //  creation of solidinterderence for Plate
            // OCCT L3619-3626
            let ssi = TopOpeBRepDSSolidSurfaceInterference::new(
                orplate,
                TopOpeBRepDSKind::Solid,
                sol_ind,
                TopOpeBRepDSKind::Surface,
                isurf,
            );
            dstr
                .change_shape_interferences(sol_ind)
                .push(TopOpeBRepDSInterference::SolidSurface(ssi));

            // calculate orientation orien of pcurves of Plate
            // the curves from ic to icplus the pcurves of Plate
            // all have the same orientation
            // OCCT L3628-3674
            let ishape1;
            let ishape2;
            let mut trafil1 = Orientation::Forward;
            let mut trafil2 = Orientation::Forward;
            let fi1_transition;
            let fi2_transition;
            {
                let fdg = fd.read().expect("surfdata lock");
                ishape1 = fdg.index_of_s1();
                ishape2 = fdg.index_of_s2();
                fi1_transition = fdg.interference_on_s1().transition();
                fi2_transition = fdg.interference_on_s2().transition();
            }
            if ishape1 != 0 {
                if ishape1 > 0 {
                    trafil1 = dstr.shape(ishape1).orientation;
                }
                trafil1 = super::chfi3d_builder_cncrn::topabs_compose(trafil1, orsurfdata);
                trafil1 = super::chfi3d_builder_cncrn::topabs_compose(
                    topabs_reverse(fi1_transition),
                    trafil1,
                );
                trafil2 = topabs_reverse(trafil1);
            } else {
                if ishape2 > 0 {
                    trafil2 = dstr.shape(ishape2).orientation;
                }
                trafil2 = super::chfi3d_builder_cncrn::topabs_compose(trafil2, orsurfdata);
                trafil2 = super::chfi3d_builder_cncrn::topabs_compose(
                    topabs_reverse(fi2_transition),
                    trafil2,
                );
                trafil1 = topabs_reverse(trafil2);
            }
            let orpcurve;
            let orien;
            {
                let st0 = cd[0].read().expect("stripe lock");
                if isfirst {
                    orpcurve = super::chfi3d_builder_cncrn::topabs_compose(
                        topabs_reverse(trafil1),
                        st0.first_pcurve_orientation(),
                    );
                } else {
                    orpcurve = super::chfi3d_builder_cncrn::topabs_compose(
                        trafil1,
                        st0.last_pcurve_orientation(),
                    );
                }
            }
            if orsurfdata == orplate {
                orien = topabs_reverse(orpcurve);
            } else {
                orien = orpcurve;
            }

            // OCCT L3676-3730 — the free-border curves on the Plate.
            if !droit {
                for ic in 0..=nedge {
                    if libre[ic as usize] {
                        let mut icplus21 = 0i32;
                        let mut icplus_l2 = 0i32;
                        let mut icmoins_l2 = 0i32;
                        indices(nedge, ic, &mut icplus_l2, &mut icmoins_l2);
                        indices(nedge, icplus_l2, &mut icplus21, &mut icmoins_l2);
                        let icp = ic as usize;
                        let icpl = icplus_l2 as usize;
                        let mut uv1 = DVec2::ZERO;
                        let mut uv2 = DVec2::ZERO;
                        let bcurv1 = BRepAdaptorCurve::initialize(&brep, &evive[icp]);
                        let bcurv2 = BRepAdaptorCurve::initialize(&brep, &evive[icpl]);
                        let par1 = p_arr[icp][icpl];
                        let par2 = p_arr[icpl][icp];
                        let (ptic, _) = bcurv1.d1(par1);
                        let (pticplus, _) = bcurv2.d1(par2);
                        parametre_plate(n3d, &psurf, surf_ref, ptic, apperror, &mut uv1);
                        parametre_plate(n3d, &psurf, surf_ref, pticplus, apperror, &mut uv2);
                        let to3d = 1.0e-3f64;
                        let to2d = 1.0e-6f64;
                        let mut cp1 = super::chfi_ds::ChFiDS_CommonPoint::default();
                        let mut cp2 = super::chfi_ds::ChFiDS_CommonPoint::default();
                        cp1.set_arc(1.0e-3, evive[icp].clone(), par1, Orientation::Forward);
                        cp1.set_point(ptic);
                        cp2.set_arc(1.0e-3, evive[icpl].clone(), par2, Orientation::Forward);
                        cp2.set_point(pticplus);
                        let (c3d, c2d, param1, param2, tolreached) = chfi3d_compute_arete(
                            &brep,
                            &cp1,
                            uv1,
                            &cp2,
                            uv2,
                            surf_ref,
                            to3d,
                            to2d,
                            0,
                        );
                        let tcurv3d = TopOpeBRepDSCurve::new(c3d, tolreached);
                        let ind1 = indpoint[icp][0];
                        let ind2 = indpoint[icpl][0];
                        let indcurv = dstr.add_curve(tcurv3d);
                        let interfp1 =
                            chfi3d_fil_point_in_ds(Orientation::Forward, indcurv, ind1, param1, false);
                        let interfp2 =
                            chfi3d_fil_point_in_ds(Orientation::Reversed, indcurv, ind2, param2, false);
                        dstr.change_curve_interferences(indcurv).push(interfp1);
                        dstr.change_curve_interferences(indcurv).push(interfp2);
                        let interfc = chfi3d_fil_curve_in_ds(indcurv, isurf, Some(c2d), orien);
                        dstr.change_surface_interferences(isurf).push(interfc);
                    }
                }
            }

            //  stockage des courbes relatives aux stripes
            // OCCT L3732-3774
            n3d = 0;
            for ic in 0..nedge {
                if !sharp[ic as usize] {
                    n3d += 1;
                    let mut icplus_l2 = 0i32;
                    let mut icmoins_l2 = 0i32;
                    indices(nedge, ic, &mut icplus_l2, &mut icmoins_l2);
                    let icp = ic as usize;
                    isfirst = sens[icp] == 1;
                    //   calculate curves interference relative to stripes

                    apperror = mapp.criterion_error() * coef;
                    let (pardeb_s, parfin_s) = {
                        let st = cd[icp].read().expect("stripe lock");
                        let pc = if isfirst {
                            st.first_pcurve()
                        } else {
                            st.last_pcurve()
                        };
                        match pc {
                            Some(c) => (c.default_domain()[0], c.default_domain()[1]),
                            None => (0.0, 0.0),
                        }
                    };
                    pardeb = pardeb_s;
                    parfin = parfin_s;

                    let interfp1 = chfi3d_fil_point_in_ds(
                        Orientation::Forward,
                        indcurve3d[n3d as usize],
                        indpoint[icp][0],
                        pardeb,
                        false,
                    );
                    let interfp2 = chfi3d_fil_point_in_ds(
                        Orientation::Reversed,
                        indcurve3d[n3d as usize],
                        indpoint[icp][1],
                        parfin,
                        false,
                    );
                    dstr
                        .change_curve_interferences(indcurve3d[n3d as usize])
                        .push(interfp1);
                    dstr
                        .change_curve_interferences(indcurve3d[n3d as usize])
                        .push(interfp2);
                    let tcourb_tol = dstr.curve(indcurve3d[n3d as usize]).tolerance();
                    dstr
                        .change_curve(indcurve3d[n3d as usize])
                        .set_tolerance(tcourb_tol + errapp[icp] + apperror);
                    let tpt1_tol = dstr.point(indpoint[icp][0]).tolerance() + apperror;
                    dstr.change_point(indpoint[icp][0]).set_tolerance(tpt1_tol);
                    let tpt2_tol = dstr.point(indpoint[icp][1]).tolerance() + apperror;
                    dstr.change_point(indpoint[icp][1]).set_tolerance(tpt2_tol);

                    // calculate surfaceinterference
                    // OCCT L3764-3767
                    let interfc = chfi3d_fil_curve_in_ds(
                        indcurve3d[n3d as usize],
                        isurf,
                        psurf.curves2d_value(n3d),
                        orien,
                    );
                    dstr.change_surface_interferences(isurf).push(interfc);
                    regular.set_curve(indcurve3d[n3d as usize]);
                    regular.set_s1(isurf, false);
                    indice = cd[icp]
                        .read()
                        .expect("stripe lock")
                        .set_of_surf_data()[(i_arr[icp][icplus_l2 as usize] - 1) as usize]
                        .read()
                        .expect("surfdata lock")
                        .surf();
                    regular.set_s2(indice, false);
                    self.my_regul.push(regular);
                }
            }

            // storage of connection curves
            // OCCT L3776-3836
            for ic in 0..nedge {
                let mut icplus_l2 = 0i32;
                let mut icmoins_l2 = 0i32;
                indices(nedge, ic, &mut icplus_l2, &mut icmoins_l2);
                let icp = ic as usize;
                if !oksea[icp] {
                    if sharp[icp] && !deuxconges {
                        // limitation of the alive edge
                        let vd = brep.first_vertex(&evive[icp]);
                        let vf = brep.last_vertex(&evive[icp]);
                        let pf = brep.vertex_position(&vd);
                        let pl = brep.vertex_position(&vf);
                        let sommet1 = brep.vertex_position(&v1);
                        let ori = if pf.distance(sommet1) < pl.distance(sommet1) {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                        let iarc1 = dstr.add_shape(&evive[icp]);
                        let interfp1 = chfi3d_fil_point_in_ds(
                            ori,
                            iarc1,
                            indpoint[icp][1],
                            p_arr[icp][icplus_l2 as usize],
                            false,
                        );
                        dstr
                            .change_shape_interferences_of(&evive[icp])
                            .push(interfp1);
                    }

                    if !ponctuel[icp] && !libre[icp] {
                        // actual connection
                        if !moresurf[icp] {
                            n3d += 1;
                            let tcourb_tol = dstr.curve(indcurve3d[n3d as usize]).tolerance();
                            dstr
                                .change_curve(indcurve3d[n3d as usize])
                                .set_tolerance(tcourb_tol + apperror);
                            if !deuxconges {
                                let tpt11_tol = dstr.point(indpoint[icp][1]).tolerance() + apperror;
                                dstr.change_point(indpoint[icp][1]).set_tolerance(tpt11_tol);
                                let tpt21_tol = dstr.point(indpoint[icplus_l2 as usize][0]).tolerance()
                                    + apperror;
                                dstr
                                    .change_point(indpoint[icplus_l2 as usize][0])
                                    .set_tolerance(tpt21_tol);
                            }
                            let interfc = chfi3d_fil_curve_in_ds(
                                indcurve3d[n3d as usize],
                                isurf,
                                psurf.curves2d_value(n3d),
                                orien,
                            );
                            dstr.change_surface_interferences(isurf).push(interfc);
                            if order[n3d as usize] == 1 {
                                regular.set_curve(indcurve3d[n3d as usize]);
                                regular.set_s1(isurf, false);
                                regular.set_s2(numfa[icp][icplus_l2 as usize], true);
                                self.my_regul.push(regular);
                            }
                        }
                    }
                }
            }

            // storage of curves projected on several faces
            // OCCT L3838-3871
            for ic in 0..nedge {
                let mut icplus_l2 = 0i32;
                let mut icmoins_l2 = 0i32;
                indices(nedge, ic, &mut icplus_l2, &mut icmoins_l2);
                let icp = ic as usize;
                let _ = icplus_l2;
                let _ = icmoins_l2;
                if moresurf[icp] {
                    for nb in 1..=nbface {
                        if !eproj[(nb - 1) as usize].is_null() {
                            n3d += 1;
                            let tcourb_tol = dstr.curve(indcurve3d[n3d as usize]).tolerance();
                            dstr
                                .change_curve(indcurve3d[n3d as usize])
                                .set_tolerance(tcourb_tol + apperror);
                            if indice_arr[n3d as usize] != 0 {
                                let tpt_tol =
                                    dstr.point(indice_arr[n3d as usize]).tolerance() + apperror;
                                dstr
                                    .change_point(indice_arr[n3d as usize])
                                    .set_tolerance(tpt_tol);
                            }
                            let interfc = chfi3d_fil_curve_in_ds(
                                indcurve3d[n3d as usize],
                                isurf,
                                psurf.curves2d_value(n3d),
                                orien,
                            );
                            dstr.change_surface_interferences(isurf).push(interfc);
                            if order[n3d as usize] == 1 {
                                regular.set_curve(indcurve3d[n3d as usize]);
                                regular.set_s1(isurf, false);
                                let fshape = dstr.add_shape(&fproj[(nb - 1) as usize]);
                                regular.set_s2(fshape, true);
                                self.my_regul.push(regular);
                            }
                        }
                    }
                }
            }

            // storage of curves in case of tangent free borders
            // OCCT L3873-3891
            if droit {
                for ic in 0..nedge {
                    let mut icplus_l2 = 0i32;
                    let mut icmoins_l2 = 0i32;
                    indices(nedge, ic, &mut icplus_l2, &mut icmoins_l2);
                    let icp = ic as usize;
                    let _ = (icplus_l2, icmoins_l2);
                    let ecom = evive[icp].clone();
                    if ecom.is_same(&edgelibre1) || ecom.is_same(&edgelibre2) {
                        n3d += 1;
                        let tcourb_tol = dstr.curve(indcurve3d[n3d as usize]).tolerance();
                        dstr
                            .change_curve(indcurve3d[n3d as usize])
                            .set_tolerance(tcourb_tol + apperror);
                        let interfc = chfi3d_fil_curve_in_ds(
                            indcurve3d[n3d as usize],
                            isurf,
                            psurf.curves2d_value(n3d),
                            orien,
                        );
                        dstr.change_surface_interferences(isurf).push(interfc);
                    }
                }
            }
        } else {
            // there is only one partial result
            // OCCT L3893-3926
            self.done = false;
            self.hasresult = true;
            for ic in 0..nedge {
                let mut icplus_l2 = 0i32;
                let mut icmoins_l2 = 0i32;
                indices(nedge, ic, &mut icplus_l2, &mut icmoins_l2);
                let icp = ic as usize;
                if !oksea[icp] {
                    if sharp[icp] && !deuxconges {
                        // limitation of the alive edge
                        let vd = brep.first_vertex(&evive[icp]);
                        let vf = brep.last_vertex(&evive[icp]);
                        let pf = brep.vertex_position(&vd);
                        let pl = brep.vertex_position(&vf);
                        let sommet1 = brep.vertex_position(&v1);
                        let ori = if pf.distance(sommet1) < pl.distance(sommet1) {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                        let iarc1 = dstr.add_shape(&evive[icp]);
                        let interfp1 = chfi3d_fil_point_in_ds(
                            ori,
                            iarc1,
                            indpoint[icp][1],
                            p_arr[icp][icplus_l2 as usize],
                            false,
                        );
                        dstr
                            .change_shape_interferences_of(&evive[icp])
                            .push(interfp1);
                    }
                }
            }
        }
    }
}

/// OCCT Geom_TrimmedCurve::Reverse — pending; the trimmed range is kept and
/// the basis curve reversed (builder_0 reverse_curve).
fn reverse_trimmed3(c: &rcad_kernel::geom::Curve3) -> rcad_kernel::geom::Curve3 {
    if let rcad_kernel::geom::Curve3::Trimmed(t) = c {
        let basis = super::chfi3d_builder_0::reverse_curve(t.curve.as_ref());
        return rcad_kernel::geom::Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(basis),
            first: t.first,
            last: t.last,
        });
    }
    c.clone()
}

/// OCCT Geom2d_TrimmedCurve::Reverse — pending; the trimmed range is kept.
fn reverse_trimmed2(c: &rcad_kernel::geom::Curve2d) -> rcad_kernel::geom::Curve2d {
    if let rcad_kernel::geom::Curve2d::Trimmed(t) = c {
        return rcad_kernel::geom::Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
            curve: t.curve.clone(),
            t_min: t.t_min,
            t_max: t.t_max,
        });
    }
    c.clone()
}
