//! OCCT ChFi3d_FilBuilder_C3.cxx (TKFillet/ChFi3d) — 1:1 translation.
//!
//! Source: C:/Users/lilu/works/OCCT/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_FilBuilder_C3.cxx (L1-1331).
//!
//! Coverage:
//!   - SearchPivot (L78-105, file static)
//!   - SearchFD (L109-211, file static)
//!   - ToricCorner (L219-233, file static)
//!   - ChFi3d_FilBuilder::PerformThreeCorner (L241-1331, the override of the
//!     pure virtual ChFi3d_Builder::PerformThreeCorner)
//!
//! Architecture mapping: OCCT `ChFi3d_FilBuilder::PerformThreeCorner` reads
//! the ChFi3d_Builder inherited members plus the FilBuilder-only
//! `myShape` (BlendFunc_SectionShape) consumed by func.Set(myShape).  The
//! rcad override is a free function over `&mut ChFi3dBuilder` taking the
//! FilBuilder member as an explicit argument (the virtual dispatch from
//! ChFi3d_Builder::PerformFilletOnVertex, ChFi3d_Builder.cxx L891, lands in
//! chfi3d.rs `perform_three_corner_pending`).
//!
//! Small ChFi3d_Builder_0.cxx helpers missing from the existing translation
//! (ChFi3d_ExtrSpineCarac L773-836, ChFi3d_CircularSpine L846-881,
//! ChFi3d_Spine L888-905) are translated here with anchors.  GAP carriers
//! (outside TKFillet or pending architecture): BlendFunc_EvolRad / Law_S
//! (the variable-radius branch), the ChFiDS_ElSpine curve storage accessors.

use std::sync::{Arc, RwLock};

use glam::{DVec2, DVec3};
use rcad_kernel::base::int_ana::intersect_plane_plane_intana;
use rcad_kernel::geom::{BezierCurve3, Circle3, Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Plane, Surface3, SurfaceEval as _};
use rcad_kernel::math::el::elclib_line_value;
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::{Orientation, Shape, TShape};

use super::chfi3d::{chfi3d_index_of_surf_data, chfi3d_index_point_in_ds, next_side, ChFi3dBuilder};
use super::chfi3d_builder_0::{
    chfi3d_compute_arete, chfi3d_enlarge_box_dstr, chfi3d_same_parameter, chfi3d_set_point_tolerance,
    vec_angle, BRepAdaptorSurface, BndBox, P_CONFUSION,
};
use super::chfi3d_builder_2::BRepTopAdaptorTopolTool;
use super::chfi3d_builder_2b::{GeomFillBoundary, GeomFillConstrainedFilling};
use super::chfi3d_builder_cncrn::chfi3d_is_in_front;
use super::chfi3d_builder_6b::{elspine_guide_curve, ChFiDSElSpineHandle};
use super::chfi3d_ds::{TopOpeBRepDSCurve, TopOpeBRepDSHDataStructure};
use super::chfi_ds::{ChFiDSElSpine, ChFiDSStripe, SharedStripe};
use super::chfi_kpart_gp::{elslib_torus_d1, GpAx3};
use super::brep_blend_func_consrad::{BlendFuncConstRad, BlendFuncConstRadInv};
use super::brep_blend_line::BRepBlendLine;

use super::chfi3d_filbuilder_c2::{
    chfi3d_mkbound_adaptor_pcurve_two_points, chfi3d_mkbound_c2d_tangents,
    chfi3d_mkbound_surface_two_points,
};

/// OCCT ChFi3d_FilBuilder_C3.cxx L78-105 — SearchPivot (file static).
fn search_pivot(s: &[i32; 3], u: &[[f64; 3]; 3], t: f64) -> i32 {
    for i in 0..=2i32 {
        // OCCT L83-90: bondeb.
        let bondeb = if s[((i + 1) % 3) as usize] == 1 {
            u[((i + 1) % 3) as usize][i as usize] - u[((i + 1) % 3) as usize][((i + 2) % 3) as usize] >= -t
        } else {
            u[((i + 1) % 3) as usize][i as usize] - u[((i + 1) % 3) as usize][((i + 2) % 3) as usize] <= t
        };
        // OCCT L91-98: bonfin.
        let bonfin = if s[((i + 2) % 3) as usize] == 1 {
            u[((i + 2) % 3) as usize][i as usize] - u[((i + 2) % 3) as usize][((i + 1) % 3) as usize] >= -t
        } else {
            u[((i + 2) % 3) as usize][i as usize] - u[((i + 2) % 3) as usize][((i + 1) % 3) as usize] <= t
        };
        if bondeb && bonfin {
            return i;
        }
    }
    -1
}

/// OCCT ChFi3d_FilBuilder_C3.cxx L109-211 — SearchFD (file static).
#[allow(clippy::too_many_arguments)]
fn search_fd(
    brep: &rcad_kernel::topods::BRep,
    dstr: &TopOpeBRepDSHDataStructure,
    cd1: &SharedStripe,
    cd2: &SharedStripe,
    sens1: i32,
    sens2: i32,
    i1: &mut i32,
    i2: &mut i32,
    p1: &mut f64,
    p2: &mut f64,
    ind1: i32,
    ind2: i32,
    face: &mut Shape,
    sameside: &mut bool,
    jf1: &mut i32,
    jf2: &mut i32,
) -> bool {
    let mut found = false;
    let mut id1 = ind1;
    let mut id2 = ind2;
    let mut if1 = ind1;
    let mut if2 = ind2;
    let l1 = {
        let g = cd1.read().expect("stripe lock");
        g.set_of_surf_data().len() as i32
    };
    let l2 = {
        let g = cd2.read().expect("stripe lock");
        g.set_of_surf_data().len() as i32
    };
    let mut fini1 = false;
    let mut fini2 = false;
    let mut visavis = false;
    let vtx = Shape::null();
    while !found {
        // OCCT L136-160.
        let mut i = id1;
        while (i * sens1) <= (if1 * sens1) && !found && !fini2 {
            if chfi3d_is_in_front(
                brep, dstr, cd1, cd2, i, if2, sens1, sens2, p1, p2, face, sameside, jf1, jf2,
                &mut visavis, &vtx, false, false,
            ) {
                *i1 = i;
                *i2 = if2;
                found = true;
            }
            i += sens1;
        }
        // OCCT L161-169.
        if !fini1 {
            if1 += sens1;
            if if1 < 1 || if1 > l1 {
                if1 -= sens1;
                fini1 = true;
            }
        }
        // OCCT L171-195.
        let mut i = id2;
        while (i * sens2) <= (if2 * sens2) && !found && !fini1 {
            if chfi3d_is_in_front(
                brep, dstr, cd1, cd2, if1, i, sens1, sens2, p1, p2, face, sameside, jf1, jf2,
                &mut visavis, &vtx, false, false,
            ) {
                *i1 = if1;
                *i2 = i;
                found = true;
            }
            i += sens2;
        }
        // OCCT L196-204.
        if !fini2 {
            if2 += sens2;
            if if2 < 1 || if2 > l2 {
                if2 -= sens2;
                fini2 = true;
            }
        }
        // OCCT L205-208.
        if fini1 && fini2 {
            break;
        }
    }
    found
}

/// OCCT ChFi3d_FilBuilder_C3.cxx L219-233 — ToricCorner (file static).
/// Test if this is a particular case of a torus corner (or spherical
/// limited by isos).
fn toric_corner(brep: &rcad_kernel::topods::BRep, f: &Shape, rd: f64, rf: f64, v: DVec3) -> bool {
    // OCCT L221-224.
    if (rd - rf).abs() > P_CONFUSION {
        return false;
    }
    // OCCT L225-229.
    let bs = BRepAdaptorSurface::initialize(brep, f);
    if bs.get_type() != super::chfi3d_builder_0::GeomAbsSurfaceType::Plane {
        return false;
    }
    // OCCT L230-232: the plane XDirection/YDirection dot v.
    let pl = match &bs.surface {
        Surface3::Plane(p) => p.clone(),
        _ => panic!("Standard_NoSuchObject: ToricCorner plane expected"),
    };
    let scal1 = pl.u_dir.dot(v).abs();
    let scal2 = pl.v_dir.dot(v).abs();
    scal1 <= P_CONFUSION && scal2 <= P_CONFUSION
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L773-836 — ChFi3d_ExtrSpineCarac.
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub(crate) fn chfi3d_extr_spine_carac(
    dstr: &TopOpeBRepDSHDataStructure,
    cd: &SharedStripe,
    i: i32,
    p: f64,
    jf: i32,
    sens: i32,
    out_p: &mut DVec3,
    out_v: &mut DVec3,
    out_r: &mut f64,
) {
    // Attention for approximated surfaces it is assumed that the parameters
    // of the pcurve are the same as of the elspine used for its
    // construction. (OCCT L779-782)
    let fffil: Surface3 = {
        let cdg = cd.read().expect("stripe lock");
        dstr
            .surface(cdg.set_of_surf_data()[(i - 1) as usize].read().expect("surfdata lock").surf())
            .surface
            .clone()
    };
    let pp: DVec2 = {
        let cdg = cd.read().expect("stripe lock");
        cdg.set_of_surf_data()[(i - 1) as usize]
            .read()
            .expect("surfdata lock")
            .interference(jf)
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(p)
    };
    // OCCT L786: P = fffil->Value(pp.X(), pp.Y()).
    *out_p = fffil.point_at(pp.x, pp.y);
    match &fffil {
        // OCCT L788-793: the cylinder branch; V = D1V.
        Surface3::Cylinder(cyl) => {
            *out_r = cyl.radius;
            let (_pbid, _vbid, v) = super::chfi_kpart::elslib_cylinder_d1(
                pp.x,
                pp.y,
                cyl.origin,
                cyl.ref_dir,
                cyl.axis,
                cyl.radius,
            );
            *out_v = v;
        }
        // OCCT L794-799: the torus branch; V = D1U.
        Surface3::Torus(tor) => {
            *out_r = tor.minor_radius;
            let ax3 = GpAx3::new_pn_vx(tor.center, tor.axis, tor.ref_dir);
            let (pbid, v, _vbid) = elslib_torus_d1(pp.x, pp.y, &ax3, tor.major_radius, tor.minor_radius);
            let _ = pbid;
            *out_v = v;
        }
        // OCCT L800-823: the default branch through the elspine.
        _ => {
            let (nbelspine, fsp, hels) = {
                let cdg = cd.read().expect("stripe lock");
                let sp = cdg.spine().expect("null spine");
                let fsp = sp.down_cast_fil().cloned();
                let nbelspine = sp.base().nb_edges() as i32;
                let hels = if nbelspine == 1 {
                    sp.base().el_spine_of_index(1).cloned()
                } else {
                    sp.base().el_spine_of_param(p).cloned()
                };
                (nbelspine, fsp, hels)
            };
            let _ = nbelspine;
            if let Some(fsp) = &fsp {
                if fsp.is_constant() {
                    // OCCT L812-815: R = fsp->Radius().
                    *out_r = fsp.radius();
                } else {
                    // OCCT L817: R = fsp->Law(hels)->Value(p).
                    // GAP carrier: ChFiDS_FilSpine::Law / Law_Function
                    // (TKMath/Law) is pending; the panic preserves the OCCT
                    // path boundary for the variable-radius corner.
                    panic!("GAP: ChFiDS_FilSpine::Law pending (OCCT ChFi3d_Builder_0.cxx L817)");
                }
            }
            // OCCT L820: hels->D1(p, Pbid, V).
            // Architecture: rcad ChFiDSElSpine carries no curve yet; the
            // guide-curve carrier of chfi3d_builder_6b.rs supplies the
            // walking guide (same convention as ComputeData).
            let hels_handle: ChFiDSElSpineHandle = Arc::new(RwLock::new(hels.expect("ElSpine")));
            let guide = elspine_guide_curve(&hels_handle);
            *out_v = guide.derivative_at(p);
        }
    }
    // OCCT L825-832.
    *out_v = out_v.normalize();
    if sens == 1 {
        *out_v = -*out_v;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L846-881 — ChFi3d_CircularSpine.
// =========================================================================
pub(crate) fn chfi3d_circular_spine(
    w_first: &mut f64,
    w_last: &mut f64,
    pdeb: DVec3,
    vdeb: DVec3,
    pfin: DVec3,
    vfin: DVec3,
    rad: f64,
) -> Option<Curve3> {
    // OCCT L854: gp_Pln Pl1(Pdeb, gp_Dir(Vdeb)), Pl2(Pfin, gp_Dir(Vfin)).
    let pl1 = Plane::new(pdeb, vdeb);
    let pl2 = Plane::new(pfin, vfin);
    // OCCT L855: IntAna_QuadQuadGeo LInt(Pl1, Pl2, Angular, Confusion).
    let lint = intersect_plane_plane_intana(&pl1, &pl2);
    match lint {
        rcad_kernel::base::int_ana::PlnPlnResult::Line(li) => {
            // OCCT L858-861: li = LInt.Line(1); cendeb/cenfin = ElCLib
            // projections of Pdeb/Pfin on li.
            let tdeb = (pdeb - li.origin).dot(li.direction);
            let cendeb = elclib_line_value(tdeb, li.origin, li.direction);
            let tfin = (pfin - li.origin).dot(li.direction);
            let cenfin = elclib_line_value(tfin, li.origin, li.direction);
            let vvdeb = pdeb - cendeb;
            let vvfin = pfin - cenfin;
            let dddeb = vvdeb.normalize();
            let ddfin = vvfin.normalize();
            // OCCT L866-870.
            if vdeb.cross(vvdeb).dot(vfin.cross(vvfin)) > 0.0 {
                return None;
            }
            // OCCT L871-874: gp_Ax2 circax2(cendeb, dddeb ^ ddfin, dddeb).
            let normal = dddeb.cross(ddfin).normalize();
            let x_dir = dddeb;
            let ccc = Circle3 {
                center: cendeb,
                normal,
                x_dir,
                y_dir: normal.cross(x_dir),
                radius: rad,
            };
            // OCCT L875-876.
            *w_first = 0.0;
            *w_last = vec_angle(dddeb, ddfin);
            Some(Curve3::Circle(ccc))
        }
        // OCCT L880: return null handle.
        _ => None,
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L888-905 — ChFi3d_Spine (the 4-pole Bezier
// guideline).
// =========================================================================
pub(crate) fn chfi3d_spine(pd: DVec3, vd: &mut DVec3, pf: DVec3, vf: &mut DVec3, r: f64) -> Curve3 {
    // OCCT L893: fac = 0.5 * tan((PI - vd.Angle(vf)) * 0.5).
    let fac = 0.5 * ((std::f64::consts::PI - vec_angle(*vd, *vf)) * 0.5).tan();
    let mut pol = [DVec3::ZERO; 4];
    pol[0] = pd;
    // OCCT L894-895: vd.Multiply(fac * R).
    *vd = *vd * (fac * r);
    pol[1] = pd + *vd;
    pol[3] = pf;
    // OCCT L897-898: vf.Multiply(fac * R).
    *vf = *vf * (fac * r);
    pol[2] = pf + *vf;
    Curve3::Bezier(BezierCurve3 {
        control_points: pol.to_vec(),
        weights: vec![1.0, 1.0, 1.0, 1.0],
    })
}

/// OCCT ChFiDS_ElSpine::Resolution (ChFiDS_ElSpine.cxx L134-137) = the
/// guide-curve resolution.  Architecture carrier: rcad ChFiDSElSpine
/// carries no curve yet (see elspine_guide_curve, chfi3d_builder_6b.rs);
/// the unit-speed identity keeps the OCCT call form until the curve
/// storage lands.
pub(crate) fn elspine_resolution(hguide: &ChFiDSElSpineHandle, r3d: f64) -> f64 {
    let _ = hguide;
    r3d
}

// =========================================================================
// OCCT ChFi3d_FilBuilder_C3.cxx L241-1331 — PerformThreeCorner.
// =========================================================================
#[allow(clippy::too_many_lines)]
pub fn perform_three_corner(
    fb: &mut ChFi3dBuilder,
    my_shape: super::brep_blend_func::BlendFuncSectionShape,
    jndex: usize,
) {
    let brep = fb.my_brep.clone();
    // OCCT L249: TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
    let mut dstr = fb.my_ds.take().expect("DS");
    // OCCT L250: const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Jndex);
    let vtx = fb.my_vdata_map.find_key(jndex).clone();
    // OCCT L252-270: Index[3], pivot, deb, fin, ii, jj, kk, pivdif,
    // c1pointu, c1toric, c1spheric, CD[3], face[3], jf[3][3], sameside[3],
    // oksea[3], i[3][3], sens[3], p[3][3], filling.
    let mut index = [0i32; 3];
    let mut pivot;
    let mut deb;
    let mut fin;
    let mut pivdif = true;
    let mut c1pointu = false;
    let mut c1toric = false;
    let mut c1spheric = false;
    let mut face = [Shape::null(), Shape::null(), Shape::null()];
    let mut jf = [[0i32; 3]; 3];
    let mut sameside = [false; 3];
    let mut oksea = [false; 3];
    let mut i_arr = [[0i32; 3]; 3];
    let mut sens = [0i32; 3];
    let mut p_arr = [[0.0f64; 3]; 3];
    let mut filling = false;

    // OCCT L272-276: the three stripes of the vertex.
    let vdata: Vec<SharedStripe> = fb.my_vdata_map.find_from_index(jndex).clone();
    let cd: Vec<SharedStripe> = vdata.iter().take(3).cloned().collect();
    for ii in 0..3usize {
        if ii >= cd.len() {
            break;
        }
        let stg = cd[ii].read().expect("stripe lock");
        index[ii] = chfi3d_index_of_surf_data(&vtx, &stg, &mut sens[ii]);
    }
    // It is checked if one of CD is not present twice in which case it is
    // necessary to modify the return of IndexOfSurfData that takes the
    // first solution. (OCCT L277-294)
    if cd.len() >= 2 && Arc::ptr_eq(&cd[0], &cd[1]) {
        index[1] = {
            let g = cd[1].read().expect("stripe lock");
            g.set_of_surf_data().len() as i32
        };
        sens[1] = -1;
    } else if cd.len() >= 3 && Arc::ptr_eq(&cd[1], &cd[2]) {
        index[2] = {
            let g = cd[2].read().expect("stripe lock");
            g.set_of_surf_data().len() as i32
        };
        sens[2] = -1;
    } else if cd.len() >= 3 && Arc::ptr_eq(&cd[0], &cd[2]) {
        index[2] = {
            let g = cd[2].read().expect("stripe lock");
            g.set_of_surf_data().len() as i32
        };
        sens[2] = -1;
    }
    // OCCT L295-339: the three SearchFD calls.  The rcad rows are split
    // into bindings so the disjoint out-parameters borrow separately.
    let [i_row0, i_row1, i_row2] = &mut i_arr;
    let [p_row0, p_row1, p_row2] = &mut p_arr;
    let [j_row0, j_row1, j_row2] = &mut jf;
    oksea[2] = if cd.len() >= 2 {
        search_fd(
            &brep, &dstr, &cd[0], &cd[1], sens[0], sens[1], &mut i_row0[1], &mut i_row1[0],
            &mut p_row0[1], &mut p_row1[0], index[0], index[1], &mut face[2], &mut sameside[2],
            &mut j_row0[1], &mut j_row1[0],
        )
    } else {
        false
    };
    oksea[1] = if cd.len() >= 3 {
        search_fd(
            &brep, &dstr, &cd[0], &cd[2], sens[0], sens[2], &mut i_row0[2], &mut i_row2[0],
            &mut p_row0[2], &mut p_row2[0], index[0], index[2], &mut face[1], &mut sameside[1],
            &mut j_row0[2], &mut j_row2[0],
        )
    } else {
        false
    };
    oksea[0] = if cd.len() >= 3 {
        search_fd(
            &brep, &dstr, &cd[1], &cd[2], sens[1], sens[2], &mut i_row1[2], &mut i_row2[1],
            &mut p_row1[2], &mut p_row2[1], index[1], index[2], &mut face[0], &mut sameside[0],
            &mut j_row1[2], &mut j_row2[1],
        )
    } else {
        false
    };
    //
    // Analyze concavities of 3 fillets :
    //        - 2 concavities identic and 1 inverted.
    //        - 3 concavities identic
    // (OCCT L340-391)
    if oksea[2] && oksea[1] && !sameside[2] && !sameside[1] {
        pivot = 0;
        deb = 1;
        fin = 2;
    } else if oksea[2] && oksea[0] && !sameside[2] && !sameside[0] {
        pivot = 1;
        deb = 2;
        fin = 0;
    } else if oksea[1] && oksea[0] && !sameside[1] && !sameside[0] {
        pivot = 2;
        deb = 0;
        fin = 1;
    } else if oksea[0] && oksea[1] && oksea[2] {
        // 3 concavities identic. (OCCT L363-386)
        pivot = search_pivot(&sens, &p_arr, fb.tol2d);
        if pivot < 0 {
            fb.perform_more_three_corner(jndex, 3);
            fb.my_ds = Some(dstr);
            return;
        } else {
            deb = (pivot + 1) % 3;
            fin = (pivot + 2) % 3;
        }
        pivdif = false;
        if (p_arr[0][1] - p_arr[0][2]).abs() <= fb.tol2d
            && (p_arr[1][0] - p_arr[1][2]).abs() <= fb.tol2d
            && (p_arr[2][0] - p_arr[2][1]).abs() <= fb.tol2d
        {
            c1pointu = true;
        }
    } else {
        // OCCT L387-391.
        fb.perform_more_three_corner(jndex, 3);
        fb.my_ds = Some(dstr);
        return;
    }
    // OCCT L392-410.
    let (ifacdeb, ifacfin) = {
        let cdd = cd[deb as usize].read().expect("stripe lock");
        let cdf = cd[fin as usize].read().expect("stripe lock");
        let ifacdeb = cdd.set_of_surf_data()[(i_arr[deb as usize][pivot as usize] - 1) as usize]
            .read()
            .expect("surfdata lock")
            .index_of(3 - jf[deb as usize][pivot as usize]);
        let ifacfin = cdf.set_of_surf_data()[(i_arr[fin as usize][pivot as usize] - 1) as usize]
            .read()
            .expect("surfdata lock")
            .index_of(3 - jf[fin as usize][pivot as usize]);
        (ifacdeb, ifacfin)
    };
    if ifacfin != ifacdeb {
        fb.perform_more_three_corner(jndex, 3);
        fb.my_ds = Some(dstr);
        return;
    }
    if i_arr[pivot as usize][deb as usize] != i_arr[pivot as usize][fin as usize] {
        fb.perform_more_three_corner(jndex, 3);
        fb.my_ds = Some(dstr);
        return;
    }

    // OCCT L412-442: the c1pointu toric/spheric analysis.
    let (mut rdeb, mut rfin, mut rdp, mut rfp);
    let mut pdeb = DVec3::ZERO;
    let mut pfin = DVec3::ZERO;
    let mut pdp = DVec3::ZERO;
    let mut pfp = DVec3::ZERO;
    let mut vdeb = DVec3::ZERO;
    let mut vfin = DVec3::ZERO;
    let mut vdp = DVec3::ZERO;
    let mut vfp = DVec3::ZERO;
    rdeb = 0.0;
    rfin = 0.0;
    rdp = 0.0;
    rfp = 0.0;
    if c1pointu {
        let mut qv = [DVec3::ZERO; 3];
        let mut qr = [0.0f64; 3];
        let mut pbid = DVec3::ZERO;
        for ii in 0..=2i32 {
            let jj = (ii + 1) % 3;
            let kk = (ii + 2) % 3;
            chfi3d_extr_spine_carac(
                &dstr,
                &cd[jj as usize],
                i_arr[jj as usize][ii as usize],
                p_arr[jj as usize][ii as usize],
                1,
                sens[jj as usize],
                &mut pbid,
                &mut qv[ii as usize],
                &mut qr[ii as usize],
            );
        }
        let mut ii = 0i32;
        while ii <= 2 && !c1toric {
            let jj = (ii + 1) % 3;
            let kk = (ii + 2) % 3;
            if toric_corner(&brep, &face[ii as usize], qr[jj as usize], qr[kk as usize], qv[ii as usize]) {
                c1toric = true;
                pivot = ii;
                deb = jj;
                fin = kk;
            }
            ii += 1;
        }
        if !c1toric {
            // OCCT L438-441.
            c1spheric =
                (qr[0] - qr[1]).abs() < fb.tolapp3d && (qr[0] - qr[2]).abs() < fb.tolapp3d;
        }
    }

    //  Previously to avoid loops the points were always located inside ...
    // (OCCT L444-485)
    let jjjd = jf[deb as usize][pivot as usize];
    let jjjf = jf[fin as usize][pivot as usize];
    chfi3d_extr_spine_carac(
        &dstr,
        &cd[deb as usize],
        i_arr[deb as usize][pivot as usize],
        p_arr[deb as usize][pivot as usize],
        jjjd,
        sens[deb as usize],
        &mut pdeb,
        &mut vdeb,
        &mut rdeb,
    );
    chfi3d_extr_spine_carac(
        &dstr,
        &cd[fin as usize],
        i_arr[fin as usize][pivot as usize],
        p_arr[fin as usize][pivot as usize],
        jjjf,
        sens[fin as usize],
        &mut pfin,
        &mut vfin,
        &mut rfin,
    );
    chfi3d_extr_spine_carac(
        &dstr,
        &cd[pivot as usize],
        i_arr[pivot as usize][deb as usize],
        p_arr[pivot as usize][deb as usize],
        0,
        sens[pivot as usize],
        &mut pdp,
        &mut vdp,
        &mut rdp,
    );
    chfi3d_extr_spine_carac(
        &dstr,
        &cd[pivot as usize],
        i_arr[pivot as usize][fin as usize],
        p_arr[pivot as usize][fin as usize],
        0,
        sens[pivot as usize],
        &mut pfp,
        &mut vfp,
        &mut rfp,
    );
    // in cas of allsame it is checked that points on the face are not
    // too close, which can stop the work. (OCCT L486-520)
    if !pivdif {
        let mut ptestdeb = DVec3::ZERO;
        let mut ptestfin = DVec3::ZERO;
        let mut bidvec = DVec3::ZERO;
        let mut bidr = 0.0f64;
        chfi3d_extr_spine_carac(
            &dstr,
            &cd[deb as usize],
            i_arr[deb as usize][pivot as usize],
            p_arr[deb as usize][pivot as usize],
            jf[deb as usize][fin as usize],
            sens[deb as usize],
            &mut ptestdeb,
            &mut bidvec,
            &mut bidr,
        );
        chfi3d_extr_spine_carac(
            &dstr,
            &cd[fin as usize],
            i_arr[fin as usize][pivot as usize],
            p_arr[fin as usize][pivot as usize],
            jf[fin as usize][deb as usize],
            sens[fin as usize],
            &mut ptestfin,
            &mut bidvec,
            &mut bidr,
        );
        let distest = ptestdeb.distance(ptestfin);
        if distest < (rdp + rfp) * 0.05 {
            filling = true;
        }
        if distest < (rdp + rfp) * 0.005 {
            c1pointu = true;
        }
    }

    // OCCT L522-533.
    if !c1pointu {
        if !pivdif {
            c1pointu = (p_arr[deb as usize][pivot as usize] - p_arr[deb as usize][fin as usize]).abs()
                <= fb.tol2d
                && (p_arr[fin as usize][pivot as usize] - p_arr[fin as usize][deb as usize]).abs()
                    <= fb.tol2d;
        }
        if (p_arr[pivot as usize][deb as usize] - p_arr[pivot as usize][fin as usize]).abs()
            <= fb.tol2d
        {
            c1toric = toric_corner(&brep, &face[pivot as usize], rdeb, rfin, vdp);
        }
    }
    // there is a pivot, the start and the end CD (finally !?!) :
    // ------------------------------------------------------------- (OCCT
    // L534-546)
    let fddeb = cd[deb as usize].read().expect("stripe lock").set_of_surf_data()
        [(i_arr[deb as usize][pivot as usize] - 1) as usize]
        .clone();
    let fdfin = cd[fin as usize].read().expect("stripe lock").set_of_surf_data()
        [(i_arr[fin as usize][pivot as usize] - 1) as usize]
        .clone();
    let fdpiv = cd[pivot as usize].read().expect("stripe lock").set_of_surf_data()
        [(i_arr[pivot as usize][deb as usize] - 1) as usize]
        .clone();

    // HSurfaces and other suitable tools are constructed.
    // ---------------------------------------------------------- (OCCT
    // L548-595)
    let ofac = face[pivot as usize].orientation;
    let fac = BRepAdaptorSurface::initialize(&brep, &face[pivot as usize]);
    let (ppp1, ppp2) = {
        let fdpiv_g = fdpiv.read().expect("surfdata lock");
        // OCCT L554-556: bid1 = InterferenceOnS1(); ppp1 = PCurveOnSurf
        // Value(FirstParameter).
        let bid1 = fdpiv_g.interference_on_s1();
        let ppp1 = bid1
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(bid1.parameter_first());
        // OCCT L557-559: bid2 = InterferenceOnS2(); ppp2 = Value(Last).
        let bid2 = fdpiv_g.interference_on_s2();
        let ppp2 = bid2
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(bid2.parameter_last());
        (ppp1, ppp2)
    };
    let (uu1, uu2, vv1, vv2) = (ppp1.x, ppp2.x, ppp1.y, ppp2.y);
    // OCCT L561-566: GeomAdaptor_Surface gasurf(DStr.Surface(...).Surface(),
    // uu1, uu2, vv1, vv2) — the rcad adaptor carries the bounds on its
    // fields.
    let pivot_surf: Surface3 = dstr
        .surface(fdpiv.read().expect("surfdata lock").surf())
        .surface
        .clone();
    let mut gasurf = BRepAdaptorSurface::initialize_surface(pivot_surf.clone());
    gasurf.ufirst = uu1;
    gasurf.ulast = uu2;
    gasurf.vfirst = vv1;
    gasurf.vlast = vv2;
    // OCCT L567-593: the parametric-range enlargements per surface kind.
    let styp = gasurf.get_type();
    if styp == super::chfi3d_builder_0::GeomAbsSurfaceType::Cylinder {
        let h = vv2 - vv1;
        gasurf.vfirst = vv1 - 0.5 * h;
        gasurf.vlast = vv2 + 0.5 * h;
        // OCCT L573-577: gasurf.Load(surf, uu1, uu2, vv1, vv2) — the rcad
        // adaptor re-Load is the field update above.
    } else if styp == super::chfi3d_builder_0::GeomAbsSurfaceType::Torus {
        let h = uu2 - uu1;
        gasurf.ufirst = uu1 - 0.1 * h;
        gasurf.ulast = uu2 + 0.1 * h;
        // OCCT L584-588: gasurf.Load(...).
    } else if styp == super::chfi3d_builder_0::GeomAbsSurfaceType::BezierSurface
        || styp == super::chfi3d_builder_0::GeomAbsSurfaceType::BSplineSurface
    {
        // OCCT L592: gasurf.Load(surf) — the full-range load is the rcad
        // adaptor's default.
    }
    // OCCT L595: Surf = new GeomAdaptor_Surface(gasurf).
    let surf = gasurf.clone();
    // OCCT L599-602: bidsurf over Fac->GeomSurfaceOriginal(); IFac/ISurf
    // TopolTools.
    let bidsurf = BRepAdaptorSurface::initialize_surface(fac.surface.clone());
    let mut ifac = BRepTopAdaptorTopolTool::default();
    ifac.initialize_brep(&bidsurf);
    let mut isurf = BRepTopAdaptorTopolTool::default();
    isurf.initialize_brep(&surf);
    // OCCT L603-608: the corner stripe + the coin surfdata.
    let corner = Arc::new(RwLock::new(ChFiDSStripe::default()));
    let coin = Arc::new(RwLock::new(super::chfi_ds::ChFiDSSurfData::default()));
    {
        let mut cg = corner.write().expect("stripe lock");
        let cornerset = cg.change_set_of_surf_data();
        cornerset.clear();
        cornerset.push(coin.clone());
    }
    // OCCT L609-646: choix + orientations.
    let mut o1 = face[pivot as usize].orientation;
    let mut o2 = fdpiv.read().expect("surfdata lock").orientation();
    let oo1 = o1;
    let oo2 = o2;
    let os1 = cd[deb as usize].read().expect("stripe lock").orientation_on_face1();
    let os2 = cd[deb as usize].read().expect("stripe lock").orientation_on_face2();
    let mut choix = cd[deb as usize].read().expect("stripe lock").choix();
    if jf[deb as usize][fin as usize] == 1 {
        choix = next_side(&mut o1, &mut o2, os1, os2, choix);
        if sens[deb as usize] == 1 {
            if choix % 2 == 1 {
                choix += 1;
            } else {
                choix -= 1;
            }
        }
    } else {
        choix = next_side(&mut o2, &mut o1, os1, os2, -choix);
        if sens[deb as usize] == -1 {
            if choix % 2 == 1 {
                choix += 1;
            } else {
                choix -= 1;
            }
        }
    }

    // OCCT L648-673: pfac1/vfac1, pfac2/vfac2, psurf1, psurf2.
    let (pfac1, vfac1) = {
        let cdd = cd[deb as usize].read().expect("stripe lock");
        let pc = cdd.set_of_surf_data()[(i_arr[deb as usize][pivot as usize] - 1) as usize]
            .read()
            .expect("surfdata lock")
            .interference(jf[deb as usize][fin as usize])
            .pcurve_on_face()
            .expect("PCurveOnFace")
            .clone();
        let par = p_arr[deb as usize][pivot as usize];
        (pc.point_at(par), pc.derivative_at(par))
    };
    let (pfac2, vfac2) = {
        let cdf = cd[fin as usize].read().expect("stripe lock");
        let pc = cdf.set_of_surf_data()[(i_arr[fin as usize][pivot as usize] - 1) as usize]
            .read()
            .expect("surfdata lock")
            .interference(jf[fin as usize][deb as usize])
            .pcurve_on_face()
            .expect("PCurveOnFace")
            .clone();
        let par = p_arr[fin as usize][pivot as usize];
        (pc.point_at(par), pc.derivative_at(par))
    };
    let psurf1 = {
        let fdpiv_g = fdpiv.read().expect("surfdata lock");
        fdpiv_g
            .interference(jf[pivot as usize][deb as usize])
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(p_arr[pivot as usize][deb as usize])
    };
    let psurf2 = {
        let fdpiv_g = fdpiv.read().expect("surfdata lock");
        fdpiv_g
            .interference(jf[pivot as usize][fin as usize])
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(p_arr[pivot as usize][fin as usize])
    };

    // OCCT L675: done = false.
    let mut done = false;

    // OCCT L685-712: the toric corner — ChFiKPart_ComputeData::ComputeCorner
    // (the toric or spheric overload, ChFiKPart_ComputeData.cxx L641-712).
    if c1toric {
        // Marshaling: S1 = Fac (face[pivot]); S2 = Surf (the GeomAdaptor
        // over the pivot fillet surface) — the rcad backend reads the
        // surface off a face-shaped wrapper (the rcad stand-in for the
        // GeomAdaptor handle argument).
        let s2_shape = Shape::new(
            Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
                my_shapes: vec![],
                flags: rcad_kernel::topo::topods::tshape_flags::DEFAULT,
                surface: Some(surf.surface.clone()),
                surface_location: 0,
                outer_wire: Shape::null(),
                inner_wires: vec![],
                sample_point: None,
                uv_domain: None,
                internal_vertices: vec![],
                tolerance: 0.0,
                natural_restriction: false,
            })),
            0,
            Orientation::Forward,
        );
        done = super::chfi_kpart::compute_data_compute_corner_cyl(
            &mut dstr,
            &mut coin.write().expect("surfdata lock"),
            &face[pivot as usize],
            &s2_shape,
            oo1,
            oo2,
            o1,
            o2,
            rdeb,
            rdp,
            pfac1,
            pfac2,
            psurf1,
            psurf2,
        );
    } else if c1spheric {
        // OCCT L713-736: the spheric corner — ChFiKPart_ComputeData::
        // ComputeCorner (the non-iso-pcurve spheric overload,
        // ChFiKPart_ComputeData.cxx L716-730).
        let s2_shape = Shape::new(
            Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
                my_shapes: vec![],
                flags: rcad_kernel::topo::topods::tshape_flags::DEFAULT,
                surface: Some(surf.surface.clone()),
                surface_location: 0,
                outer_wire: Shape::null(),
                inner_wires: vec![],
                sample_point: None,
                uv_domain: None,
                internal_vertices: vec![],
                tolerance: 0.0,
                natural_restriction: false,
            })),
            0,
            Orientation::Forward,
        );
        done = super::chfi_kpart::compute_data_compute_corner_sphere(
            &mut dstr,
            &mut coin.write().expect("surfdata lock"),
            &face[pivot as usize],
            &s2_shape,
            oo1,
            oo2,
            o1,
            o2,
            rdp,
            pfac1,
            psurf1,
            psurf2,
        );
    } else if c1pointu {
        // OCCT L737-740.
        filling = true;
    }
    // OCCT L741-1015.
    if !done {
        if !filling {
            // Calculate a guideline (OCCT L750-787).
            let radpondere = (rdp + rfp) / 2.0;
            let mut locfleche = fb.fleche;
            let mut w_first = 0.0f64;
            let mut w_last = 0.0f64;
            let spinecoin = chfi3d_circular_spine(
                &mut w_first, &mut w_last, pdeb, vdeb, pfin, vfin, radpondere,
            );
            let spinecoin = match spinecoin {
                Some(sc) => {
                    // OCCT L774-776.
                    locfleche = radpondere * (w_last - w_first) * fb.fleche;
                    sc
                }
                None => {
                    // This is a bad case when the intersection of section
                    // planes is done out of the sector. (OCCT L767-771)
                    w_first = 0.0;
                    w_last = 1.0;
                    chfi3d_spine(pdeb, &mut vdeb, pfin, &mut vfin, radpondere)
                }
            };
            let pasmax = (w_last - w_first) * 0.05;
            // OCCT L778-784: cornerspine = new ChFiDS_ElSpine;
            // SetCurve(spinecoin); FirstParameter(WFirst - pasmax);
            // LastParameter(WLast + pasmax).  The ChFiDSElSpine curve field
            // (ChFiDS_ElSpine.hxx L148, ChFiDS_ElSpine.cxx L299-302) carries
            // the spinecoin curve; the walking guide reads it through
            // elspine_guide_curve.
            let cornerspine = Arc::new(RwLock::new(ChFiDSElSpine {
                curve: Some(spinecoin),
                firstparam: w_first - pasmax,
                lastparam: w_last + pasmax,
                firstpnt: DVec3::ZERO,
                firsttgt: DVec3::ZERO,
                lastpnt: DVec3::ZERO,
                lasttgt: DVec3::ZERO,
                vertices_with_tangents: Vec::new(),
                period: 0.0,
                periodic: false,
                // OCCT ChFiDS_ElSpine.cxx L42-43: pfirstsav/plastsav ctor
                // defaults are Precision::Infinite().
                pfirstsav: rcad_kernel::core::precision::INFINITE_VALUE,
                plastsav: rcad_kernel::core::precision::INFINITE_VALUE,
                next: None,
                previous: None,
            }));
            // OCCT L784: NullSpine.
            let null_spine: Option<&super::chfi_ds::ChFiDSSpineHandle> = None;
            // OCCT L788-792: math_Vector Soldep(1, 4).
            let mut soldep = Vector::new(1, 4);
            soldep.set(1, pfac1.x);
            soldep.set(2, pfac1.y);
            soldep.set(3, psurf1.x);
            soldep.set(4, psurf1.y);

            let mut gd1 = false;
            let mut gd2 = false;
            let mut gf1 = false;
            let mut gf2 = false;
            let mut lin: Option<BRepBlendLine> = None;
            let mut ffi = w_first;
            let mut lla = w_last + pasmax;

            if (rdeb - rfin).abs() <= fb.tolapp3d {
                // OCCT L798-836: the constant-radius ComputeData.
                // Architecture: BRepBlend_ConstRad is a typedef of
                // BlendFunc_ConstRad (BRepBlend_ConstRad.hxx); the rcad
                // BlendFuncConstRad implements the BlendFunction trait.
                let guide = elspine_guide_curve(&cornerspine);
                let mut func = BlendFuncConstRad::new(&fac.surface, &surf.surface, &guide);
                func.set(rdeb, choix);
                func.set_section_shape(my_shape);
                let mut finv = BlendFuncConstRadInv::new(&fac.surface, &surf.surface, &guide);
                finv.set(rdeb, choix);
                let tol_guide = elspine_resolution(&cornerspine, fb.tolapp3d);
                let mut intf = 3i32;
                let mut intl = 3i32;
                done = fb.compute_data(
                    &mut coin.write().expect("surfdata lock"),
                    &cornerspine,
                    null_spine,
                    &mut lin,
                    &fac,
                    &ifac,
                    &surf,
                    &isurf,
                    &mut func,
                    &mut finv,
                    ffi,
                    pasmax,
                    locfleche,
                    tol_guide,
                    &mut ffi,
                    &mut lla,
                    false,
                    false,
                    true,
                    &soldep,
                    &mut intf,
                    &mut intl,
                    &mut gd1,
                    &mut gd2,
                    &mut gf1,
                    &mut gf2,
                    false,
                    true,
                );
                if done && gf2 {
                    let line_guard = lin.clone().expect("Lin");
                    done = fb.complete_data_function(
                        &mut coin.write().expect("surfdata lock"),
                        &mut func,
                        &line_guard,
                        &fac,
                        Some(&surf),
                        ofac,
                        gd1,
                        false,
                        gf1,
                        false,
                        false,
                    );
                    filling = !done;
                } else {
                    filling = true;
                }
            } else {
                // OCCT L851-903: the variable-radius ComputeData through
                // Law_S + BRepBlend_EvolRad (= typedef BlendFunc_EvolRad).
                // GAP carrier: BlendFunc_EvolRad (OCCT BlendFunc/
                // BlendFunc_EvolRad.cxx) and Law_S (TKMath/Law/Law_S.cxx)
                // are not translated; the panic preserves the OCCT path
                // boundary until the variable-radius batch lands.
                panic!("GAP: BlendFunc_EvolRad / Law_S pending (OCCT ChFi3d_FilBuilder_C3.cxx L853-903)");
            }
        }

        if filling {
            // the contour to be fillet consists of straight lines uv in
            // beginning and end of two pcurves ... (OCCT L917-1015)
            let mut pcurve_on_face: Option<Curve2d> = None;
            let bfac = if !c1pointu {
                // OCCT L924-934.
                let b = chfi3d_mkbound_c2d_tangents(
                    &fac,
                    &mut pcurve_on_face,
                    sens[deb as usize],
                    pfac1,
                    vfac1,
                    sens[fin as usize],
                    pfac2,
                    vfac2,
                    fb.tolapp3d,
                    2.0e-4,
                );
                b
            } else {
                GeomFillBoundary
            };
            // OCCT L935-943: the 3d tangents at the pivot interference
            // parameters.
            let vp1 = {
                let kkk = {
                    let cdd = cd[deb as usize].read().expect("stripe lock");
                    cdd.set_of_surf_data()[(i_arr[deb as usize][pivot as usize] - 1) as usize]
                        .read()
                        .expect("surfdata lock")
                        .interference(jf[deb as usize][pivot as usize])
                        .line_index()
                };
                let curv = dstr
                    .curve(kkk)
                    .curve
                    .clone()
                    .expect("DS curve null");
                curv.derivative_at(p_arr[deb as usize][pivot as usize])
            };
            let vp2 = {
                let kkk = {
                    let cdf = cd[fin as usize].read().expect("stripe lock");
                    cdf.set_of_surf_data()[(i_arr[fin as usize][pivot as usize] - 1) as usize]
                        .read()
                        .expect("surfdata lock")
                        .interference(jf[fin as usize][pivot as usize])
                        .line_index()
                };
                let curv = dstr
                    .curve(kkk)
                    .curve
                    .clone()
                    .expect("DS curve null");
                curv.derivative_at(p_arr[fin as usize][pivot as usize])
            };
            let _ = (vp1, vp2);
            // OCCT L944-947: Bpiv = ChFi3d_mkbound(Surf, PCurveOnPiv,
            // psurf1, psurf2, tolapp3d, 2.e-4, false).
            let mut pcurve_on_piv: Option<Curve2d> = None;
            let bpiv = chfi3d_mkbound_adaptor_pcurve_two_points(
                &surf,
                &mut pcurve_on_piv,
                psurf1,
                psurf2,
                fb.tolapp3d,
                2.0e-4,
                false,
            );
            // OCCT L948-954: pardeb2/parfin2.
            let pardeb2 = if c1pointu {
                p_arr[deb as usize][fin as usize]
            } else {
                p_arr[deb as usize][pivot as usize]
            };
            let parfin2 = if c1pointu {
                p_arr[fin as usize][deb as usize]
            } else {
                p_arr[fin as usize][pivot as usize]
            };
            // OCCT L955-978: pdeb1/pdeb2/pfin1/pfin2.
            let pdeb1 = {
                let cdd = cd[deb as usize].read().expect("stripe lock");
                cdd.set_of_surf_data()[(i_arr[deb as usize][pivot as usize] - 1) as usize]
                    .read()
                    .expect("surfdata lock")
                    .interference(jf[deb as usize][pivot as usize])
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(p_arr[deb as usize][pivot as usize])
            };
            let pdeb2 = {
                let cdd = cd[deb as usize].read().expect("stripe lock");
                cdd.set_of_surf_data()[(i_arr[deb as usize][pivot as usize] - 1) as usize]
                    .read()
                    .expect("surfdata lock")
                    .interference(jf[deb as usize][fin as usize])
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(pardeb2)
            };
            let pfin1 = {
                let cdf = cd[fin as usize].read().expect("stripe lock");
                cdf.set_of_surf_data()[(i_arr[fin as usize][pivot as usize] - 1) as usize]
                    .read()
                    .expect("surfdata lock")
                    .interference(jf[fin as usize][pivot as usize])
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(p_arr[fin as usize][pivot as usize])
            };
            let pfin2 = {
                let cdf = cd[fin as usize].read().expect("stripe lock");
                cdf.set_of_surf_data()[(i_arr[fin as usize][pivot as usize] - 1) as usize]
                    .read()
                    .expect("surfdata lock")
                    .interference(jf[fin as usize][deb as usize])
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(parfin2)
            };
            // OCCT L979-985.
            let sdeb: Surface3 = {
                let cdd = cd[deb as usize].read().expect("stripe lock");
                dstr
                    .surface(
                        cdd.set_of_surf_data()[(i_arr[deb as usize][pivot as usize] - 1) as usize]
                            .read()
                            .expect("surfdata lock")
                            .surf(),
                    )
                    .surface
                    .clone()
            };
            let sfin: Surface3 = {
                let cdf = cd[fin as usize].read().expect("stripe lock");
                dstr
                    .surface(
                        cdf.set_of_surf_data()[(i_arr[fin as usize][pivot as usize] - 1) as usize]
                            .read()
                            .expect("surfdata lock")
                            .surf(),
                    )
                    .surface
                    .clone()
            };
            let bdeb = chfi3d_mkbound_surface_two_points(&sdeb, pdeb1, pdeb2, fb.tolapp3d, 2.0e-4);
            let bfin = chfi3d_mkbound_surface_two_points(&sfin, pfin1, pfin2, fb.tolapp3d, 2.0e-4);

            // OCCT L987-995: GeomFill_ConstrainedFilling fil(11, 20).
            let mut fil = GeomFillConstrainedFilling::new(11, 20);
            if c1pointu {
                fil.init3(&bpiv, &bfin, &bdeb, true);
            } else {
                fil.init4(&bpiv, &bfin, &bfac, &bdeb, true);
            }

            // OCCT L997-1010: Surfcoin = fil.Surface(); Surfcoin->VReverse();
            // done = CompleteData(coin, Surfcoin, Fac, PCurveOnFace, Surf,
            // PCurveOnPiv, fdpiv->Orientation(), false, false, false, false,
            // false).
            done = match fil.surface() {
                Some(surfcoin) => {
                    // OCCT L998: Surfcoin->VReverse() — the V parametrization
                    // flip is a pending Surface3 operation (the branch only
                    // becomes reachable once GeomFill lands, see
                    // chfi3d_builder_2b.rs L1464 for the same convention).
                    let _ = &surfcoin;
                    let fdpiv_orientation = fdpiv.read().expect("surfdata lock").orientation();
                    fb.complete_data_surfcoin(
                        &mut coin.write().expect("surfdata lock"),
                        &surfcoin,
                        &fac,
                        pcurve_on_face.as_ref(),
                        &surf,
                        pcurve_on_piv.as_ref(),
                        fdpiv_orientation,
                        false,
                        false,
                        false,
                        false,
                        false,
                    )
                }
                None => {
                    // GAP: GeomFill_ConstrainedFilling::Surface pending
                    // (chfi3d_builder_2b.rs) — the OCCT flow always produces
                    // a surface; report the not-done path.
                    false
                }
            };
        }
    }
    // OCCT L1017-1022.
    let mut p1deb;
    let mut p2deb;
    let mut p1fin;
    let mut p2fin;
    if !c1pointu {
        p_arr[deb as usize][fin as usize] = p_arr[deb as usize][pivot as usize];
        p_arr[fin as usize][deb as usize] = p_arr[fin as usize][pivot as usize];
    }

    // OCCT L1024-1283.
    if done {
        // Update of 4 Stripes and the DS
        // ------------------------------------- (OCCT L1026-1054)
        let (pf1, pf2): (super::chfi_ds::ChFiDS_CommonPoint, super::chfi_ds::ChFiDS_CommonPoint);
        {
            let mut coing = coin.write().expect("surfdata lock");
            // OCCT L1036-1037: Pf1/Pf2.
            pf1 = coing.vertex_first_on_s1().clone();
            pf2 = coing.vertex_first_on_s2().clone();
            // OCCT L1038-1042: Pl1 = ChangeVertexLastOnS1(); when c1pointu
            // the slot is overwritten with VertexFirstOnS1.
            if c1pointu {
                let first = coing.vertex_first_on_s1().clone();
                *coing.change_vertex_last_on_s1() = first;
            }
        }
        let pl1 = {
            let coing = coin.read().expect("surfdata lock");
            coing.vertex_last_on_s1().clone()
        };
        let pl2 = {
            let coing = coin.read().expect("surfdata lock");
            coing.vertex_last_on_s2().clone()
        };

        // OCCT L1045-1054: Bnd_Box bf1/bl1/bf2/bl2 + the pointers
        // pbf1/pbl1/pbf2/pbl2 (pbl1 aliases pbf1 when c1pointu).  The rcad
        // encoding resolves the alias per use site (c1pointu branch).
        let mut bf1 = BndBox::default();
        let mut bl1 = BndBox::default();
        let mut bf2 = BndBox::default();
        let mut bl2 = BndBox::default();
        bf1.add(pf1.point());
        bf2.add(pf2.point());
        if c1pointu {
            bf1.add(pl1.point());
        } else {
            bl1.add(pl1.point());
        }
        bl2.add(pl2.point());

        // the start corner,
        // ----------------------- (OCCT L1056-1107)
        let if1 = chfi3d_index_point_in_ds(&pf1, &mut dstr);
        let if2 = chfi3d_index_point_in_ds(&pf2, &mut dstr);
        let il1 = if c1pointu {
            if1
        } else {
            chfi3d_index_point_in_ds(&pl1, &mut dstr)
        };
        let il2 = chfi3d_index_point_in_ds(&pl2, &mut dstr);
        let (pp1, pp2) = {
            let coing = coin.read().expect("surfdata lock");
            (
                coing
                    .interference_on_s1()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(coing.interference_on_s1().parameter_first()),
                coing
                    .interference_on_s2()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(coing.interference_on_s2().parameter_first()),
            )
        };
        {
            let mut coing = coin.write().expect("surfdata lock");
            if c1pointu {
                coing.change_index_of_s1(0);
            } else {
                coing.change_index_of_s1(dstr.add_shape(&face[pivot as usize]));
            }
            coing.change_index_of_s2(-(fdpiv.read().expect("surfdata lock").surf()));
        }
        let coin_surf: Surface3 = dstr
            .surface(coin.read().expect("surfdata lock").surf())
            .surface
            .clone();
        let (c3d, _pcurv, cp1deb, cp2deb, tolreached) =
            chfi3d_compute_arete(&brep, &pf1, pp1, &pf2, pp2, &coin_surf, fb.tolapp3d, fb.tol2d, 0);
        let tcurv1 = TopOpeBRepDSCurve::new(c3d, tolreached);
        let icf = dstr.add_curve(tcurv1);
        let mut regdeb = super::chfi_ds::ChFiDSRegul::default();
        regdeb.set_curve(icf);
        regdeb.set_s1(coin.read().expect("surfdata lock").surf(), false);
        regdeb.set_s2(fddeb.read().expect("surfdata lock").surf(), false);
        fb.my_regul.push(regdeb);
        {
            let mut cornerg = corner.write().expect("stripe lock");
            cornerg.change_first_curve(icf);
            cornerg.change_first_parameters(cp1deb, cp2deb);
            cornerg.change_index_first_point_on_s1(if1);
            cornerg.change_index_first_point_on_s2(if2);
        }
        {
            let coing = coin.read().expect("surfdata lock");
            // OCCT L1107: EnlargeBox(DStr, corner, coin, *pbf1, *pbf2, true).
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                Some(&corner.read().expect("stripe lock")),
                &coing,
                &mut bf1,
                &mut bf2,
                true,
            );
        }

        // OCCT L1109-1134: the last corner.
        let (pp1, pp2) = {
            let coing = coin.read().expect("surfdata lock");
            (
                coing
                    .interference_on_s1()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(coing.interference_on_s1().parameter_last()),
                coing
                    .interference_on_s2()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(coing.interference_on_s2().parameter_last()),
            )
        };
        let (c3d, _pcurv, cp1fin, cp2fin, tolreached) =
            chfi3d_compute_arete(&brep, &pl1, pp1, &pl2, pp2, &coin_surf, fb.tolapp3d, fb.tol2d, 0);
        let tcurv2 = TopOpeBRepDSCurve::new(c3d, tolreached);
        let icl = dstr.add_curve(tcurv2);
        let mut regfin = super::chfi_ds::ChFiDSRegul::default();
        regfin.set_curve(icl);
        regfin.set_s1(coin.read().expect("surfdata lock").surf(), false);
        regfin.set_s2(fdfin.read().expect("surfdata lock").surf(), false);
        fb.my_regul.push(regfin);
        {
            let mut cornerg = corner.write().expect("stripe lock");
            cornerg.change_last_curve(icl);
            cornerg.change_last_parameters(cp1fin, cp2fin);
            cornerg.change_index_last_point_on_s1(il1);
            cornerg.change_index_last_point_on_s2(il2);
        }
        {
            let coing = coin.read().expect("surfdata lock");
            // OCCT L1134: EnlargeBox(DStr, corner, coin, *pbl1, *pbl2,
            // false) — pbl1 is pbf1 when c1pointu.
            if c1pointu {
                chfi3d_enlarge_box_dstr(
                    &brep,
                    &dstr,
                    Some(&corner.read().expect("stripe lock")),
                    &coing,
                    &mut bf1,
                    &mut bl2,
                    false,
                );
            } else {
                chfi3d_enlarge_box_dstr(
                    &brep,
                    &dstr,
                    Some(&corner.read().expect("stripe lock")),
                    &coing,
                    &mut bl1,
                    &mut bl2,
                    false,
                );
            }
        }

        // then CornerData of the beginning,
        // -------------------------------- (OCCT L1136-1181)
        let mut isfirst = sens[deb as usize] == 1;
        let mut rev = jf[deb as usize][fin as usize] == 2;
        let mut isurf1 = 1i32;
        let mut isurf2 = 2i32;
        let mut parpp1 = p_arr[deb as usize][fin as usize];
        let mut parpp2 = p_arr[deb as usize][pivot as usize];
        if rev {
            isurf1 = 2;
            isurf2 = 1;
            parpp1 = p_arr[deb as usize][pivot as usize];
            parpp2 = p_arr[deb as usize][fin as usize];
            cd[deb as usize]
                .write()
                .expect("stripe lock")
                .set_orientation(Orientation::Reversed, isfirst);
        }
        let (pp1, pp2) = {
            let fddebg = fddeb.read().expect("surfdata lock");
            (
                fddebg
                    .interference_on_s1()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(parpp1),
                fddebg
                    .interference_on_s2()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(parpp2),
            )
        };
        {
            let mut cddeb_g = cd[deb as usize].write().expect("stripe lock");
            cddeb_g.set_curve(icf, isfirst);
            cddeb_g.set_index_point(if1, isfirst, isurf1);
            cddeb_g.set_index_point(if2, isfirst, isurf2);
            cddeb_g.set_parameters(isfirst, cp1deb, cp2deb);
        }
        {
            let mut fddebg = fddeb.write().expect("surfdata lock");
            *fddebg.change_vertex(isfirst, isurf1) = pf1.clone();
            *fddebg.change_vertex(isfirst, isurf2) = pf2.clone();
            fddebg
                .change_interference_on_s1()
                .set_parameter(isfirst, parpp1);
            fddebg
                .change_interference_on_s2()
                .set_parameter(isfirst, parpp2);
        }
        // OCCT L1160-1173: tcdeb = DStr.ChangeCurve(Icf); crefdeb =
        // tcdeb.Curve(); pp1/pp2 are the interference values of L1150-1151;
        // ChFi3d_ComputePCurv(crefdeb, pp1, pp2, CD[deb]->ChangePCurve(
        // isfirst), DStr.Surface(fddeb->Surf()).Surface(), P1deb, P2deb,
        // tolapp3d, tolrdeb, rev); tcdeb.Tolerance(max(tolrdeb, old)).
        let tolrdeb = {
            let sd_surf: Surface3 = dstr
                .surface(fddeb.read().expect("surfdata lock").surf())
                .surface
                .clone();
            let mut tolrdeb = 0.0f64;
            let (crefdeb, tc_tol) = {
                let tcdeb = dstr.change_curve(icf);
                (tcdeb.curve().cloned(), tcdeb.tolerance())
            };
            if let Some(crefdeb) = &crefdeb {
                // OCCT: CD[deb]->ChangePCurve(isfirst) is the out slot; the
                // ComputePCurv assigns a fresh curve into it.
                let mut pcurv = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                    origin: DVec2::ZERO,
                    direction: DVec2::Y,
                });
                super::chfi3d_filbuilder_c2::chfi3d_compute_pcurv_c3d(
                    crefdeb,
                    pp1,
                    pp2,
                    &mut pcurv,
                    &sd_surf,
                    cp1deb,
                    cp2deb,
                    fb.tolapp3d,
                    &mut tolrdeb,
                    rev,
                );
                cd[deb as usize]
                    .write()
                    .expect("stripe lock")
                    .change_pcurve(isfirst, pcurv);
            }
            tolrdeb.max(tc_tol)
        };
        dstr.change_curve(icf).set_tolerance(tolrdeb);
        if rev {
            let fddebg = fddeb.read().expect("surfdata lock");
            chfi3d_enlarge_box_dstr(
                &brep, &dstr, Some(&cd[deb as usize].read().expect("stripe lock")), &fddebg,
                &mut bf2, &mut bf1, isfirst,
            );
        } else {
            let fddebg = fddeb.read().expect("surfdata lock");
            chfi3d_enlarge_box_dstr(
                &brep, &dstr, Some(&cd[deb as usize].read().expect("stripe lock")), &fddebg,
                &mut bf1, &mut bf2, isfirst,
            );
        }

        // then the end CornerData,
        // ------------------------ (OCCT L1183-1230)
        isfirst = sens[fin as usize] == 1;
        rev = jf[fin as usize][deb as usize] == 2;
        isurf1 = 1;
        isurf2 = 2;
        let mut parpp1 = p_arr[fin as usize][deb as usize];
        let mut parpp2 = p_arr[fin as usize][pivot as usize];
        if rev {
            isurf1 = 2;
            isurf2 = 1;
            parpp1 = p_arr[fin as usize][pivot as usize];
            parpp2 = p_arr[fin as usize][deb as usize];
            cd[fin as usize]
                .write()
                .expect("stripe lock")
                .set_orientation(Orientation::Reversed, isfirst);
        }
        let (pp1, pp2) = {
            let fdfing = fdfin.read().expect("surfdata lock");
            (
                fdfing
                    .interference_on_s1()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(parpp1),
                fdfing
                    .interference_on_s2()
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(parpp2),
            )
        };
        {
            let mut cdfin_g = cd[fin as usize].write().expect("stripe lock");
            cdfin_g.set_curve(icl, isfirst);
            cdfin_g.set_index_point(il1, isfirst, isurf1);
            cdfin_g.set_index_point(il2, isfirst, isurf2);
            cdfin_g.set_parameters(isfirst, cp1fin, cp2fin);
        }
        {
            let mut fdfing = fdfin.write().expect("surfdata lock");
            *fdfing.change_vertex(isfirst, isurf1) = pl1.clone();
            *fdfing.change_vertex(isfirst, isurf2) = pl2.clone();
            fdfing
                .change_interference_on_s1()
                .set_parameter(isfirst, parpp1);
            fdfing
                .change_interference_on_s2()
                .set_parameter(isfirst, parpp2);
        }
        // OCCT L1209-1222: tcfin + ChFi3d_ComputePCurv(creffin, pp1, pp2,
        // CD[fin]->ChangePCurve(isfirst), ..., P1fin, P2fin, tolapp3d,
        // tolrfin, rev) + tcfin.Tolerance(max(tolrfin, old)).
        let tolrfin = {
            let sd_surf: Surface3 = dstr
                .surface(fdfin.read().expect("surfdata lock").surf())
                .surface
                .clone();
            let mut tolrfin = 0.0f64;
            let (creffin, tc_tol) = {
                let tcfin = dstr.change_curve(icl);
                (tcfin.curve().cloned(), tcfin.tolerance())
            };
            if let Some(creffin) = &creffin {
                // OCCT: the out slot of ChangePCurve(isfirst).
                let mut pcurv = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                    origin: DVec2::ZERO,
                    direction: DVec2::Y,
                });
                super::chfi3d_filbuilder_c2::chfi3d_compute_pcurv_c3d(
                    creffin,
                    pp1,
                    pp2,
                    &mut pcurv,
                    &sd_surf,
                    cp1fin,
                    cp2fin,
                    fb.tolapp3d,
                    &mut tolrfin,
                    rev,
                );
                cd[fin as usize]
                    .write()
                    .expect("stripe lock")
                    .change_pcurve(isfirst, pcurv);
            }
            tolrfin.max(tc_tol)
        };
        dstr.change_curve(icl).set_tolerance(tolrfin);
        if rev {
            let fdfing = fdfin.read().expect("surfdata lock");
            // OCCT L1225: EnlargeBox(DStr, CD[fin], fdfin, *pbl2, *pbl1,
            // isfirst) — pbl1 is pbf1 when c1pointu.
            if c1pointu {
                chfi3d_enlarge_box_dstr(
                    &brep, &dstr, Some(&cd[fin as usize].read().expect("stripe lock")), &fdfing,
                    &mut bl2, &mut bf1, isfirst,
                );
            } else {
                chfi3d_enlarge_box_dstr(
                    &brep, &dstr, Some(&cd[fin as usize].read().expect("stripe lock")), &fdfing,
                    &mut bl2, &mut bl1, isfirst,
                );
            }
        } else {
            let fdfing = fdfin.read().expect("surfdata lock");
            // OCCT L1229: EnlargeBox(DStr, CD[fin], fdfin, *pbl1, *pbl2,
            // isfirst).
            if c1pointu {
                chfi3d_enlarge_box_dstr(
                    &brep, &dstr, Some(&cd[fin as usize].read().expect("stripe lock")), &fdfing,
                    &mut bf1, &mut bl2, isfirst,
                );
            } else {
                chfi3d_enlarge_box_dstr(
                    &brep, &dstr, Some(&cd[fin as usize].read().expect("stripe lock")), &fdfing,
                    &mut bl1, &mut bl2, isfirst,
                );
            }
        }

        // anf finally the pivot.
        // ------------------ (OCCT L1232-1276)
        let isfirst = sens[pivot as usize] == 1;
        let rev = jf[pivot as usize][deb as usize] == 2;
        let mut isurf1 = 1i32;
        let mut isurf2 = 2i32;
        if rev {
            isurf1 = 2;
            isurf2 = 1;
            cd[pivot as usize]
                .write()
                .expect("stripe lock")
                .set_orientation(Orientation::Reversed, isfirst);
        }
        let (iccoinpiv, fi_first, fi_last, fi_pcurve) = {
            let mut coing = coin.write().expect("surfdata lock");
            let fi = coing.change_interference_on_s2();
            let iccoinpiv = fi.line_index();
            let fi_first = fi.parameter_first();
            let fi_last = fi.parameter_last();
            let fi_pcurve = fi.change_pcurve_on_face().clone();
            (iccoinpiv, fi_first, fi_last, fi_pcurve)
        };
        {
            cd[pivot as usize]
                .write()
                .expect("stripe lock")
                .set_curve(iccoinpiv, isfirst);
        }
        let ccoinpiv = dstr.curve(iccoinpiv).curve.clone();
        // OCCT L1252-1259: ChFi3d_SameParameter(Ccoinpiv, C2dOnPiv, Spiv,
        // fi.FirstParameter(), fi.LastParameter(), tolapp3d, tolr) +
        // TCcoinpiv.Tolerance(max(TCcoinpiv.Tolerance(), tolr)).
        let tolr = {
            let mut c2d_on_piv = fi_pcurve.clone();
            let spiv: Surface3 = dstr
                .surface(fdpiv.read().expect("surfdata lock").surf())
                .surface
                .clone();
            let mut tolr = 0.0f64;
            if let (Some(ccoinpiv), Some(c2d)) = (&ccoinpiv, c2d_on_piv.as_mut()) {
                chfi3d_same_parameter(ccoinpiv, c2d, &spiv, fb.tolapp3d, &mut tolr);
            }
            tolr.max(dstr.curve(iccoinpiv).tolerance())
        };
        dstr.change_curve(iccoinpiv).set_tolerance(tolr);
        if let Some(c2d) = fi_pcurve {
            cd[pivot as usize]
                .write()
                .expect("stripe lock")
                .change_pcurve(isfirst, c2d);
        }
        {
            let mut cdpiv = cd[pivot as usize].write().expect("stripe lock");
            cdpiv.set_index_point(if2, isfirst, isurf1);
            cdpiv.set_index_point(il2, isfirst, isurf2);
            cdpiv.set_parameters(isfirst, fi_first, fi_last);
        }
        {
            let mut fdpivg = fdpiv.write().expect("surfdata lock");
            *fdpivg.change_vertex(isfirst, isurf1) = pf2.clone();
            *fdpivg.change_vertex(isfirst, isurf2) = pl2.clone();
            fdpivg
                .change_interference(isurf1)
                .set_parameter(isfirst, p_arr[pivot as usize][deb as usize]);
            fdpivg
                .change_interference(isurf2)
                .set_parameter(isfirst, p_arr[pivot as usize][fin as usize]);
        }
        {
            // OCCT L1268: CD[pivot]->InDS(isfirst) — filDS already does it
            // from the corner.
            cd[pivot as usize]
                .write()
                .expect("stripe lock")
                .in_ds(isfirst, 0);
        }
        if rev {
            let fdpivg = fdpiv.read().expect("surfdata lock");
            // OCCT L1271: EnlargeBox(DStr, CD[pivot], fdpiv, *pbl2, *pbf2,
            // isfirst) — pbl2 carries no alias.
            chfi3d_enlarge_box_dstr(
                &brep, &dstr, Some(&cd[pivot as usize].read().expect("stripe lock")), &fdpivg,
                &mut bl2, &mut bf2, isfirst,
            );
        } else {
            let fdpivg = fdpiv.read().expect("surfdata lock");
            // OCCT L1275: EnlargeBox(DStr, CD[pivot], fdpiv, *pbf2, *pbl2,
            // isfirst).
            chfi3d_enlarge_box_dstr(
                &brep, &dstr, Some(&cd[pivot as usize].read().expect("stripe lock")), &fdpivg,
                &mut bf2, &mut bl2, isfirst,
            );
        }

        // To end the tolerances of points are rescaled. (OCCT L1278-1282)
        chfi3d_set_point_tolerance(&mut dstr, &bf1, if1);
        chfi3d_set_point_tolerance(&mut dstr, &bf2, if2);
        if c1pointu {
            chfi3d_set_point_tolerance(&mut dstr, &bf1, il1);
        } else {
            chfi3d_set_point_tolerance(&mut dstr, &bl1, il1);
        }
        chfi3d_set_point_tolerance(&mut dstr, &bl2, il2);
        p1deb = cp1deb;
        p2deb = cp2deb;
        p1fin = cp1fin;
        p2fin = cp2fin;
    } else {
        p1deb = 0.0;
        p2deb = 0.0;
        p1fin = 0.0;
        p2fin = 0.0;
    }
    let _ = (p1deb, p2deb, p1fin, p2fin);

    // The data corners are truncated and index is updated.
    //---------------------------------------------------- (OCCT L1285-1318)
    {
        let mut cdd = cd[deb as usize].write().expect("stripe lock");
        if i_arr[deb as usize][pivot as usize] < index[deb as usize] {
            let seq = cdd.change_set_of_surf_data();
            seq.drain((i_arr[deb as usize][pivot as usize] as usize)..(index[deb as usize] as usize));
            index[deb as usize] = i_arr[deb as usize][pivot as usize];
        } else if i_arr[deb as usize][pivot as usize] > index[deb as usize] {
            let seq = cdd.change_set_of_surf_data();
            seq.drain((index[deb as usize] as usize)..((i_arr[deb as usize][pivot as usize] - 1) as usize));
            i_arr[deb as usize][pivot as usize] = index[deb as usize];
        }
    }
    {
        let mut cdf = cd[fin as usize].write().expect("stripe lock");
        if i_arr[fin as usize][pivot as usize] < index[fin as usize] {
            let seq = cdf.change_set_of_surf_data();
            seq.drain((i_arr[fin as usize][pivot as usize] as usize)..(index[fin as usize] as usize));
            index[fin as usize] = i_arr[fin as usize][pivot as usize];
        } else if i_arr[fin as usize][pivot as usize] > index[fin as usize] {
            let seq = cdf.change_set_of_surf_data();
            seq.drain((index[fin as usize] as usize)..((i_arr[fin as usize][pivot as usize] - 1) as usize));
            i_arr[fin as usize][pivot as usize] = index[fin as usize];
        }
    }
    // it is necessary to take into account mutant corners. (OCCT L1308-1318)
    {
        let mut cdp = cd[pivot as usize].write().expect("stripe lock");
        if i_arr[pivot as usize][deb as usize] < index[pivot as usize] {
            let seq = cdp.change_set_of_surf_data();
            seq.drain((i_arr[pivot as usize][deb as usize] as usize)..(index[pivot as usize] as usize));
            index[pivot as usize] = i_arr[pivot as usize][deb as usize];
        } else if i_arr[pivot as usize][deb as usize] > index[pivot as usize] {
            let seq = cdp.change_set_of_surf_data();
            seq.drain((index[pivot as usize] as usize)..((i_arr[pivot as usize][deb as usize] - 1) as usize));
            i_arr[pivot as usize][deb as usize] = index[pivot as usize];
        }
    }
    // OCCT L1319-1326.
    if !fb.my_evi_map.contains_key(&vtx.ptr_id()) {
        fb.my_evi_map.insert(vtx.ptr_id(), Vec::new());
    }
    let coin_surf = coin.read().expect("surfdata lock").surf();
    fb.my_evi_map
        .get_mut(&vtx.ptr_id())
        .expect("EVIMap")
        .push(coin_surf);
    {
        let solid_index = cd[pivot as usize].read().expect("stripe lock").solid_index();
        corner.write().expect("stripe lock").set_solid_index(solid_index);
    }
    fb.my_list_stripe.push(corner.clone());

    fb.my_ds = Some(dstr);
}
