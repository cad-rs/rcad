//! OCCT ChFi3d_Builder corner-tail machinery — 1:1 translation.
//!
//! Despite the module name (kept from the port-plan batch naming), the two
//! member functions translated here live in OCCT
//! `ChFi3d_Builder_C1.cxx`:
//!   - ChFi3d_Builder::PerformIntersectionAtEnd (C1.cxx L1828-3763)
//!   - ChFi3d_Builder::PerformMoreSurfdata      (C1.cxx L3771-4529)
//! plus the file statics they use:
//!   - cherche_face (C1.cxx L1660-1696), containV / containE
//!     (C1.cxx L1738-1786), IsShrink (C1.cxx L1792-1822)
//!   - ChFi3d_nbface (Builder_0.cxx L5468-5489)
//! and the GeomLib (TKGeomBase) extension primitives invoked from both the
//! corner tails and `ChFi3d_ExtendSurface` (chfi3d_builder_c1.rs):
//!   - GeomLib::ExtendSurfByLength  (GeomLib.cxx L1485-1962)
//!   - GeomLib::ExtendCurveToPoint  (GeomLib.cxx L1269-1410)
//!   - BSplCLib::TangExtendToConstraint (BSplCLib.cxx L3794-4160)
//!   - PLib::CoefficientsPoles, dim version (PLib.cxx L1522-1608)
//!
//! GAP carriers (outside-package dependencies, OCCT failure path kept):
//!   - GeomInt_IntSS pcurves (LineOnS1/LineOnS2, TolReached3d): the rcad
//!     kernel IntSS (GeomAPI_IntSS stand-in) returns 3D lines only.  Where
//!     OCCT consumes the pcurves, the rcad translation carries None and the
//!     dependent extend-branches keep the OCCT not-done behavior.
//!   - GeomLib KPart extension (ExtendKPart, GeomLib.cxx L1420+) and
//!     GeomConvert_ApproxSurface: a non-BSpline/Bezier surface is returned
//!     unextended.
//!   - GeomConvert_CompCurveToBSplineCurve + ComputeLambda
//!     (GeomLib.cxx L126-200): ExtendCurveToPoint keeps the curve
//!     unextended (the branch is only reachable through the IntSS-pcurve
//!     GAP above).

use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{
    BSplineSurface, BezierSurface, Curve2dEval as _, Curve3, CurveEval as _, Surface3,
    SurfaceEval as _,
};
use rcad_kernel::math::bspl_lib::{
    eval_flat, increase_degree as bspl_increase_degree, knot_sequence as bspl_knot_sequence,
    remove_knot as bspl_remove_knot,
};
use rcad_kernel::math::math_matrix::{Matrix, Vector};
use rcad_kernel::math::plib::hermite_coefficients;
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};

use super::chfi3d::{
    chfi3d_index_of_surf_data, chfi3d_index_point_in_ds, is_tangent_faces, next_side,
    topabs_compose, topabs_reverse, ChFi3dBuilder,
};
use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_bound_fac, chfi3d_bound_surf, chfi3d_boite,
    chfi3d_compute_arete, chfi3d_compute_curves, chfi3d_compute_pcurv_2pt, chfi3d_couture,
    chfi3d_couture_on_vertex, chfi3d_eval_tol_reached, chfi3d_nb_not_degenerated_edges,
    chfi3d_project_pcurv, topexp_face_edges, topexp_face_vertices,
    topexp_vertices, BRepAdaptorSurface, GeomAdaptorSurface,
};
use super::chfi3d_builder_c1::{
    chfi3d_recale, compute_curve2d, inters_update_on_same, project_point_on_curve2d,
};
use super::chfi3d_builder_c2_geomlib::geom_lib_extend_surf_by_length;
use super::chfi3d_builder_cncrn::{
    cherche_edge1, chfi3d_angle_edge, chfi3d_cherche_edge, chfi3d_cherche_element,
    chfi3d_cherche_face1, chfi3d_cherche_vertex, chfi3d_edge_common_faces,
};
use super::chfi3d_builder_0::chfi3d_edge_state;
use super::chfi3d_ds::{
    TopOpeBRepDSCurve, TopOpeBRepDSCurvePointInterference, TopOpeBRepDSSolidSurfaceInterference,
    TopOpeBRepDSPoint, TopOpeBRepDSSurface,
};
use super::chfi_ds::{ChFiDS_CommonPoint, ChFiDS_State, ChFiDSSurfData, SharedSurfData};

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5468-5489 — ChFi3d_nbface.
// =========================================================================
pub(crate) fn chfi3d_nbface(map_vf: &[Shape]) -> i32 {
    let mut nface = 0i32;
    for (fj, cur) in map_vf.iter().enumerate() {
        let fj = (fj + 1) as i32;
        let mut kf = 1i32;
        let mut counted = true;
        for jt in map_vf.iter().take((fj - 1) as usize) {
            if cur.is_same(jt) {
                counted = false;
                break;
            }
            kf += 1;
        }
        if counted && kf == fj {
            nface += 1;
        }
    }
    nface
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L1660-1696 — cherche_face: find face F
// belonging to the map, different from faces F1 F2 F3 and containing edge E.
// =========================================================================
pub(crate) fn cherche_face(
    brep: &rcad_kernel::topods::BRep,
    map: &[Shape],
    e: &Shape,
    f1: &Shape,
    f2: &Shape,
    f3: &Shape,
    f: &mut Shape,
) {
    let mut trouve = false;
    for fcur in map {
        if trouve {
            break;
        }
        if !fcur.is_same(f1) && !fcur.is_same(f2) && !fcur.is_same(f3) {
            for ecur in topexp_face_edges(brep, fcur) {
                if e.is_same(&ecur) {
                    *f = fcur.clone();
                    trouve = true;
                    break;
                }
            }
        }
    }
    if f.is_null() {
        panic!("Standard_ConstructionError: Failed to find face.");
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L1738-1756 — containV: true if vertex V
// belongs to face F1.
// =========================================================================
pub(crate) fn contain_v(brep: &rcad_kernel::topods::BRep, f1: &Shape, v: &Shape) -> bool {
    topexp_face_vertices(brep, f1).iter().any(|vcur| vcur.is_same(v))
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L1760-1786 — containE: true if edge E belongs
// to face F1.
// =========================================================================
pub(crate) fn contain_e(brep: &rcad_kernel::topods::BRep, f1: &Shape, e: &Shape) -> bool {
    topexp_face_edges(brep, f1).iter().any(|ecur| ecur.is_same(e))
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L1792-1822 — IsShrink (chfi3d_is_shrink):
// check if U (if <is_u>) or V of points of <pc> is within <tol> from
// <param>, points between <pf> and <pl>.
// =========================================================================
pub(crate) fn chfi3d_is_shrink(
    pc: &rcad_kernel::geom::Curve2d,
    pf: f64,
    pl: f64,
    param: f64,
    is_u: bool,
    tol: f64,
) -> bool {
    let coord = |p: DVec2| if is_u { p.x } else { p.y };
    match pc {
        // GeomAbs_Line
        rcad_kernel::geom::Curve2d::Line(_) => {
            let p1 = pc.point_at(pf);
            let p2 = pc.point_at(pl);
            (coord(p1) - param).abs() <= tol && (coord(p2) - param).abs() <= tol
        }
        // GeomAbs_BezierCurve / GeomAbs_BSplineCurve: math_FunctionSample
        // (Pf, Pl, 10) — 11 sample points.
        rcad_kernel::geom::Curve2d::BSpline(_) | rcad_kernel::geom::Curve2d::Bezier(_) => {
            for i in 0..=10i32 {
                let t = pf + (pl - pf) * (i as f64) / 10.0;
                let p = pc.point_at(t);
                if (coord(p) - param).abs() > tol {
                    return false;
                }
            }
            true
        }
        _ => false,
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L1828-3763 — PerformIntersectionAtEnd.
// Intersection at end of fillet with at least two faces.
// =========================================================================
/// The myEFMap(E) list (an empty stand-in keeps the OCCT null-map case).
pub(crate) fn ef_list(efmap: &super::chfi_ds::ChFiDSMap, e: &Shape) -> Vec<Shape> {
    if efmap.contains(e) {
        efmap.find(e).clone()
    } else {
        Vec::new()
    }
}

/// OCCT Geom_Curve::FirstParameter (the trim range for a trimmed carrier).
pub(crate) fn c3d_first(c: &Curve3) -> f64 {
    c.default_domain()[0]
}

/// OCCT Geom_Curve::LastParameter.
pub(crate) fn c3d_last(c: &Curve3) -> f64 {
    c.default_domain()[1]
}

/// OCCT ChFi3d_Builder_C1.cxx L3373/L3391: Extrema_ExtPC ext(pext, cad,
/// tolpt); par1 = ext.Point(1).Parameter() — the first extremum parameter
/// over the loaded (cad) curve; Point(1) raises StdFail_NotDone when the
/// extrema were not found, mirrored by the None return the call sites
/// expect() on.
fn extrema_ext_pc_first_parameter(c: &Curve3, p: DVec3, tol: f64) -> Option<f64> {
    // OCCT L3372: cad.Load(csau) — the GeomAdaptor_Curve full-range load.
    let a_adaptor = GeomCurveAdaptor::new(c.clone());
    let a_tool = CurveToolHandle::for_curve3(c, &a_adaptor, &a_adaptor);
    let ext = ExtremaExtPC::new_point_curve(p, &a_tool, tol);
    if ext.is_done() && ext.nb_ext() >= 1 {
        Some(ext.point(1).param)
    } else {
        None
    }
}

pub fn perform_intersection_at_end(this: &mut ChFi3dBuilder, index: usize) {
    let pi = std::f64::consts::PI;
    // OCCT L1837: TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
    // OCCT L1839: const int nn = 15.
    let nn = 15usize;
    // OCCT L1841: stripe = It.Value() (first of myVDataMap(Index)).
    let stripe = match this.my_vdata_map.find_from_index(index).first() {
        Some(s) => s.clone(),
        None => return,
    };
    let mut st = stripe.write().expect("stripe lock");
    let spine = st.spine().expect("spine").clone();
    // OCCT L1838-1843: SeqFil = stripe->ChangeSetOfSurfData()->ChangeSequence().
    // OCCT L1843: const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index).
    let vtx = this.my_vdata_map.find_key(index).clone();
    let brep = this.my_brep.clone();
    let tolapp3d = this.tolapp3d;
    let tol2d = this.tol2d;

    // determine the number of faces and edges (OCCT L1845-1860).
    let vf_list: Vec<Shape> = if this.my_vf_map.contains(&vtx) {
        this.my_vf_map.find(&vtx).clone()
    } else {
        Vec::new()
    };
    let mut nface = chfi3d_nbface(&vf_list);
    let (bordlibre, edgelibre1, edgelibre2) =
        super::chfi3d_builder_0::chfi3d_cherche_bords_libres(&this.my_ve_map, &vtx);
    let mut nbarete = chfi3d_nb_not_degenerated_edges(&vtx, &this.my_ve_map) as i32;
    if bordlibre {
        nbarete = (nbarete - 2) / 2 + 2;
    } else {
        nbarete /= 2;
    }

    // it is determined if there is an edge of sewing and it face
    // (OCCT L1863-1875).
    let mut facecouture = Shape::null();
    let mut edgecouture = Shape::null();
    let mut couture = false;
    for itf in &vf_list {
        if couture {
            break;
        }
        let fcur = itf.clone();
        let (cout, ecout) = chfi3d_couture_on_vertex(&brep, &fcur, &vtx);
        if cout {
            couture = true;
            edgecouture = ecout;
            facecouture = fcur;
        }
    }

    // it is determined if one of edges adjacent to the fillet is regular
    // (OCCT L1877-1896).
    let mut reg1 = false;
    let mut reg2 = false;
    let mut eadj1 = Shape::null();
    let mut eadj2 = Shape::null();
    let nbsurf = st.my_hdata.len();
    let nbedge = spine.base().nb_edges();
    let mut sens = 0i32;
    let mut num = chfi3d_index_of_surf_data(&vtx, &st, &mut sens);
    let isfirst = sens == 1;
    let mut num1;
    let edge_spine = if isfirst {
        num1 = num + 1;
        spine.base().edges(1).clone()
    } else {
        num1 = num - 1;
        spine.base().edges(nbedge).clone()
    };
    let state = if isfirst {
        spine.base().first_status()
    } else {
        spine.base().last_status()
    };

    if nbsurf != nbedge && nbsurf != 1 {
        // two faces common to the edge are found (OCCT L1898-1940).
        let mut f1 = Shape::null();
        let mut f2 = Shape::null();
        chfi3d_edge_common_faces(
            &ef_list(&this.my_ef_map, &edge_spine),
            &mut f1,
            &mut f2,
        );
        if f1.is_same(&facecouture) {
            eadj1 = edgecouture.clone();
        } else {
            let mut vbid1 = Shape::null();
            chfi3d_cherche_element(&brep, &vtx, &edge_spine, &f1, &mut eadj1, &mut vbid1);
        }
        let mut fga = Shape::null();
        let mut fdr = Shape::null();
        chfi3d_edge_common_faces(&ef_list(&this.my_ef_map, &eadj1), &mut fga, &mut fdr);
        // Modified by Sergey KHROMOV: reg1 = ChFi3d::IsTangentFaces(Eadj1, Fga, Fdr).
        reg1 = is_tangent_faces(
            &brep,
            &eadj1,
            &fga,
            &fdr,
            crate::geomalgo::gtests_stubs::GeomAbsShape::G1,
        );
        if f2.is_same(&facecouture) {
            eadj2 = edgecouture.clone();
        } else {
            let mut vbid1 = Shape::null();
            chfi3d_cherche_element(&brep, &vtx, &edge_spine, &f2, &mut eadj2, &mut vbid1);
        }
        chfi3d_edge_common_faces(&ef_list(&this.my_ef_map, &eadj2), &mut fga, &mut fdr);
        reg2 = is_tangent_faces(
            &brep,
            &eadj2,
            &fga,
            &fdr,
            crate::geomalgo::gtests_stubs::GeomAbsShape::G1,
        );

        if reg1 || reg2 {
            let mut compoint1 = false;
            let mut compoint2 = false;
            let (cp1, cp2) = {
                let sd = st.my_hdata[(num1 - 1) as usize].read().expect("surfdata lock");
                (
                    sd.vertex(isfirst, 1).clone(),
                    sd.vertex(isfirst, 2).clone(),
                )
            };
            if cp1.is_on_arc() && (cp1.arc().is_same(&eadj1) || cp1.arc().is_same(&eadj2)) {
                compoint1 = true;
            }
            if cp2.is_on_arc() && (cp2.arc().is_same(&eadj1) || cp2.arc().is_same(&eadj2)) {
                compoint2 = true;
            }
            if compoint1 && compoint2 {
                st.my_hdata.remove((num - 1) as usize);
                let mut sens2 = 0i32;
                let num2 = chfi3d_index_of_surf_data(&vtx, &st, &mut sens2);
                num = num2;
                if isfirst {
                    num1 = num + 1;
                } else {
                    num1 = num - 1;
                }
                reg1 = false;
                reg2 = false;
            }
        }
    }

    // there is only one face at end if FindFace is true and if the face is
    // not the face with sewing edge (OCCT L1942-1956).
    let mut face = Shape::null();
    let fd_lock = st.my_hdata[(num - 1) as usize].clone();
    let mut fd = fd_lock.write().expect("surfdata lock");
    let mut cv1 = fd.vertex(isfirst, 1).clone();
    let mut cv2 = fd.vertex(isfirst, 2).clone();
    let mut onecorner = false;
    if this.find_face(&cv1, &cv2, &mut face) {
        if !couture {
            onecorner = true;
        } else if !face.is_same(&facecouture) {
            onecorner = true;
        }
    }
    // OCCT L1958-1963: if (onecorner) { if (MoreSurfdata(Index)) { ... } }.
    // rcad architecture: the stripe write lock is held here and the
    // MoreSurfdata query reads the same stripe (RwLock is not re-entrant),
    // so the check runs on the held guard (the translated MoreSurfdata
    // tail reads the stripe surfdata count).
    if onecorner && st.my_hdata.len() > 1 {
        drop(fd);
        drop(st);
        super::chfi3d_builder_c2b::perform_more_surfdata(this, index);
        return;
    }
    if !onecorner && (reg1 || reg2) && !couture && state != ChFiDS_State::OnSame {
        // OCCT L1963-1967: PerformMoreThreeCorner(Index, 1); return.
        drop(fd);
        drop(st);
        this.perform_more_three_corner(index, 1);
        return;
    }

    // calculate the orientation of curves at end (OCCT L1967-2058).
    let hgs = chfi3d_bound_surf(this.my_ds.as_ref().expect("DS"), &fd, 1, 2);
    let mut fi1 = fd.interference_on_s1().clone();
    let mut fi2 = fd.interference_on_s2().clone();
    let orsurfdata = fd.orientation();
    let mut isurf_prev = 0i32;
    let mut isurf = fd.surf();
    let sdprev_opt: Option<SharedSurfData> = if num1 > 0 && num1 <= st.my_hdata.len() as i32 {
        Some(st.my_hdata[(num1 - 1) as usize].clone())
    } else {
        None
    };
    if let Some(sdprev) = &sdprev_opt {
        isurf_prev = sdprev.read().expect("surfdata lock").surf();
    }

    let tolpt = 1.0e-4;
    let mut tolreached = 0.0f64;
    let mut orface = Orientation::Forward;
    let mut orien = Orientation::Forward;

    st.set_index_point(
        chfi3d_index_point_in_ds(&cv1, this.my_ds.as_mut().expect("DS")),
        isfirst,
        1,
    );
    st.set_index_point(
        chfi3d_index_point_in_ds(&cv2, this.my_ds.as_mut().expect("DS")),
        isfirst,
        2,
    );

    let ishape1 = fd.index_of_s1;
    let mut trafil1 = Orientation::Forward;
    if ishape1 != 0 {
        if ishape1 > 0 {
            trafil1 = this
                .my_ds
                .as_ref()
                .expect("DS")
                .shape(ishape1)
                .orientation;
        }
        // OCCT L2043: trafil1 = TopAbs::Compose(trafil1, Fd->Orientation()).
        trafil1 = topabs_compose(trafil1, fd.orientation());
        // OCCT L2045: trafil1 = TopAbs::Compose(TopAbs::Reverse(Fi1.Transition()), trafil1).
        trafil1 = topabs_compose(topabs_reverse(fi1.transition()), trafil1);
    }
    // OCCT L2060-2064.
    let orcourbe = if !isfirst {
        trafil1
    } else {
        topabs_reverse(trafil1)
    };

    // eap, Apr 22 2002, occ 293 — variables to show OnSame situation
    // (OCCT L2085-2098).
    let mut is_on_same1 = false;
    let mut is_on_same2 = false;
    let mut is_shrink = false;
    let mut is_u_shrink = false;
    let mut check_shr_param = 0.0f64;
    let mut prev_sd_param = 0.0f64;
    let mut mid_p2d = DVec2::ZERO;
    let mut mid_ipoint = 0i32;

    // find Fi1,Fi2 lengths used to extend ChFi surface and by the way
    // define necessity to check shrink (OCCT L2101-2124).
    let p2d1 = fi1
        .pcurve_on_surf()
        .expect("Fi1 PCurveOnSurf")
        .point_at(fi1.parameter(isfirst));
    let p2d2 = fi1
        .pcurve_on_surf()
        .expect("Fi1 PCurveOnSurf")
        .point_at(fi1.parameter(!isfirst));
    let a_p1 = hgs.surface().point_at(p2d1.x, p2d1.y);
    let a_p2 = hgs.surface().point_at(p2d2.x, p2d2.y);
    let fi1_length = a_p1.distance(a_p2);
    let mut check_shrink = fi1_length <= super::chfi3d_builder_0::P_CONFUSION;

    let p2d3 = fi2
        .pcurve_on_surf()
        .expect("Fi2 PCurveOnSurf")
        .point_at(fi2.parameter(isfirst));
    let p2d4 = fi2
        .pcurve_on_surf()
        .expect("Fi2 PCurveOnSurf")
        .point_at(fi2.parameter(!isfirst));
    let a_p1 = hgs.surface().point_at(p2d3.x, p2d3.y);
    let a_p2 = hgs.surface().point_at(p2d4.x, p2d4.y);
    let fi2_length = a_p1.distance(a_p2);
    check_shrink = check_shrink || (fi2_length <= super::chfi3d_builder_0::P_CONFUSION);

    if check_shrink {
        // OCCT L2126-2141.
        if (p2d2.y - p2d4.y).abs() <= super::chfi3d_builder_0::P_CONFUSION {
            is_u_shrink = false;
            check_shr_param = p2d2.y;
        } else if (p2d2.x - p2d4.x).abs() <= super::chfi3d_builder_0::P_CONFUSION {
            is_u_shrink = true;
            check_shr_param = p2d2.x;
        } else {
            check_shrink = false;
        }
    }

    /**********************************************************************/
    // find faces intersecting with the fillet and edges limiting
    // intersections — nbface is the nb of faces intersected, Face[i]
    // contains the faces to intersect, Edge[i] the limiting edges
    // (OCCT L2128-2169).
    /**********************************************************************/
    let mut nb = 1usize;
    let mut nbface;
    let mut edge: Vec<Shape> = (0..=nn).map(|_| Shape::null()).collect();
    let mut e = Shape::null();
    let mut ei;
    let mut edgesau = Shape::null();
    let mut facesau = Shape::null();
    let mut oneintersection1 = false;
    let mut oneintersection2 = false;
    let mut face_arr: Vec<Shape> = (0..=nn).map(|_| Shape::null()).collect();
    let mut f3;
    let mut v = Shape::null();
    let mut findonf1 = false;
    let mut findonf2 = false;
    let f1 = this
        .my_ds
        .as_ref()
        .expect("DS")
        .shape(fd.index_of_s1)
        .clone();
    let f2 = this
        .my_ds
        .as_ref()
        .expect("DS")
        .shape(fd.index_of_s2)
        .clone();
    f3 = f1.clone();
    if couture || bordlibre {
        nface += 1;
    }
    if nface == 3 {
        nbface = 2;
    } else {
        nbface = (nface - 2) as usize;
    }
    if !cv1.is_on_arc() || !cv2.is_on_arc() {
        // OCCT L2165-2169: PerformMoreThreeCorner(Index, 1); return.
        drop(fd);
        drop(st);
        this.perform_more_three_corner(index, 1);
        return;
    }

    edge[0] = cv1.arc().clone();
    edge[nbface] = cv2.arc().clone();
    // processing of a fillet arriving on a vertex: edge contained in CV.Arc
    // is not inevitably good — the edge concerned by the intersection is
    // found (OCCT L2171-2232).
    if cv1.is_vertex() {
        let mut trouve = false;
        v = cv1.vertex().clone();
        let ve_list: Vec<Shape> = if this.my_ve_map.contains(&v) {
            this.my_ve_map.find(&v).clone()
        } else {
            Vec::new()
        };
        for it3 in &ve_list {
            if trouve {
                break;
            }
            e = it3.clone();
            if !e.is_same(&edge[0]) && contain_e(&brep, &f1, &e) {
                trouve = true;
            }
        }
        let (ev1, ev2) = topexp_vertices(&edge[0]);
        let vt;
        if v.is_same(&ev1) {
            vt = ev2.clone();
        } else {
            vt = ev1.clone();
        }
        let dist1 = brep.vertex_position(&vt).distance(brep.vertex_position(&vtx));
        let (ev3, ev4) = topexp_vertices(&e);
        let vt;
        if v.is_same(&ev3) {
            vt = ev4.clone();
        } else {
            vt = ev3.clone();
        }
        let dist2 = brep.vertex_position(&vt).distance(brep.vertex_position(&vtx));
        if dist2 < dist1 && trouve {
            edge[0] = e.clone();
            let ori = if ev2.is_same(&ev3) || ev1.is_same(&ev4) {
                cv1.transition_on_arc()
            } else {
                topabs_reverse(cv1.transition_on_arc())
            };
            let par = brep_tool_parameter(&brep, &v, &edge[0]);
            let tol = cv1.tolerance();
            cv1.set_arc(tol, edge[0].clone(), par, ori);
            *fd.change_vertex(isfirst, 1) = cv1.clone();
        }
    }
    if cv2.is_vertex() {
        // OCCT L2233-2290 — the same treatment for CV2 (Edge[2] slot).
        let mut trouve = false;
        v = cv2.vertex().clone();
        let ve_list: Vec<Shape> = if this.my_ve_map.contains(&v) {
            this.my_ve_map.find(&v).clone()
        } else {
            Vec::new()
        };
        for it3 in &ve_list {
            if trouve {
                break;
            }
            e = it3.clone();
            if !e.is_same(&edge[2]) && contain_e(&brep, &f2, &e) {
                trouve = true;
            }
        }
        let (ev1, ev2) = topexp_vertices(&edge[2]);
        let vt;
        if v.is_same(&ev1) {
            vt = ev2.clone();
        } else {
            vt = ev1.clone();
        }
        let dist1 = brep.vertex_position(&vt).distance(brep.vertex_position(&vtx));
        let (ev3, ev4) = topexp_vertices(&e);
        let vt;
        if v.is_same(&ev3) {
            vt = ev4.clone();
        } else {
            vt = ev3.clone();
        }
        let dist2 = brep.vertex_position(&vt).distance(brep.vertex_position(&vtx));
        if dist2 < dist1 && trouve {
            edge[2] = e.clone();
            let ori = if ev2.is_same(&ev3) || ev1.is_same(&ev4) {
                cv2.transition_on_arc()
            } else {
                topabs_reverse(cv2.transition_on_arc())
            };
            let par = brep_tool_parameter(&brep, &v, &edge[2]);
            let tol = cv2.tolerance();
            cv2.set_arc(tol, edge[2].clone(), par, ori);
            *fd.change_vertex(isfirst, 2) = cv2.clone();
        }
    }
    if !onecorner {
        // If there is a regular edge, the faces adjacent to it are not in
        // Fd->IndexOfS1 or Fd->IndexOfS2 (OCCT L2292-2318).
        if nface == 3 {
            // nface = 3: a top with 3 edges and a fillet whose common
            // points are on different faces (OCCT L2320-2350).
            if cv1.is_vertex() {
                findonf1 = true;
            }
            if cv2.is_vertex() {
                findonf2 = true;
            }
            if !findonf1 {
                let ev = topexp_vertices(&edge[0]);
                if !ev.0.is_same(&vtx) && !ev.1.is_same(&vtx) {
                    findonf1 = true;
                }
            }
            if !findonf2 {
                let ev = topexp_vertices(&edge[2]);
                if !ev.0.is_same(&vtx) && !ev.1.is_same(&vtx) {
                    findonf2 = true;
                }
            }

            // detect and process OnSame situation (OCCT L2352-2374).
            if state == ChFiDS_State::OnSame {
                let mut three_e: [Shape; 3] = [Shape::null(), Shape::null(), Shape::null()];
                let mut vb2 = Shape::null();
                chfi3d_cherche_element(&brep, &vtx, &edge_spine, &f1, &mut three_e[0], &mut vb2);
                chfi3d_cherche_element(&brep, &vtx, &edge_spine, &f2, &mut three_e[1], &mut vb2);
                three_e[2] = edge_spine.clone();
                if chfi3d_edge_state(&three_e, &this.my_ef_map, &brep) == ChFiDS_State::OnSame {
                    is_on_same1 = true;
                    nb = 1;
                    edge[0] = three_e[0].clone();
                    chfi3d_cherche_face1(
                        &ef_list(&this.my_ef_map, &edge[0]),
                        &f1,
                        &mut face_arr[0],
                    );
                    if findonf2 {
                        findonf1 = true; // not to look for Face[0] again
                    } else {
                        edge[1] = cv2.arc().clone();
                    }
                } else {
                    is_on_same2 = true;
                }
            }

            // findonf1 findonf2 show if F1 and/or F2 are adjacent to many
            // faces at end — the faces at end and intersected edges are
            // found (OCCT L2377-2400).
            if findonf1 && !is_on_same1 {
                v = if cv1.transition_on_arc() == Orientation::Forward {
                    brep.first_vertex(cv1.arc())
                } else {
                    brep.last_vertex(cv1.arc())
                };
                chfi3d_cherche_face1(&ef_list(&this.my_ef_map, cv1.arc()), &f1, &mut face_arr[0]);
                nb = 1;
                ei = edge[0].clone();
                while !v.is_same(&vtx) {
                    let mut v2l = Shape::null();
                    chfi3d_cherche_element(&brep, &v, &ei, &f1, &mut e, &mut v2l);
                    v = v2l;
                    ei = e.clone();
                    chfi3d_cherche_face1(&ef_list(&this.my_ef_map, &e), &f1, &mut face_arr[nb]);
                    cherche_edge1(&brep, &face_arr[nb - 1], &face_arr[nb], &mut edge[nb]);
                    nb += 1;
                    if nb >= nn {
                        panic!("Standard_Failure: IntersectionAtEnd : the max number of faces reached");
                    }
                }
                if !findonf2 {
                    edge[nb] = cv2.arc().clone();
                }
            }
            if findonf2 && !is_on_same2 {
                if !findonf1 {
                    nb = 1;
                }
                v = vtx.clone();
                let vfin = if cv2.transition_on_arc() == Orientation::Forward {
                    brep.last_vertex(cv2.arc())
                } else {
                    brep.first_vertex(cv2.arc())
                };
                if !findonf1 {
                    chfi3d_cherche_face1(
                        &ef_list(&this.my_ef_map, cv1.arc()),
                        &f1,
                        &mut face_arr[nb - 1],
                    );
                }
                let mut v2l = Shape::null();
                chfi3d_cherche_element(&brep, &v, &edge_spine, &f2, &mut e, &mut v2l);
                ei = e.clone();
                v = v2l;
                while !v.is_same(&vfin) {
                    let mut v2l = Shape::null();
                    chfi3d_cherche_element(&brep, &v, &ei, &f2, &mut e, &mut v2l);
                    ei = e.clone();
                    v = v2l;
                    chfi3d_cherche_face1(&ef_list(&this.my_ef_map, &e), &f2, &mut face_arr[nb]);
                    cherche_edge1(&brep, &face_arr[nb - 1], &face_arr[nb], &mut edge[nb]);
                    nb += 1;
                    if nb >= nn {
                        panic!("Standard_Failure: IntersectionAtEnd : the max number of faces reached");
                    }
                }
                edge[nb] = cv2.arc().clone();
            }
            if is_on_same2 {
                cherche_edge1(&brep, &face_arr[nb - 1], &f2, &mut edge[nb]);
                face_arr[nb] = f2.clone();
            }

            nbface = nb;
        } else {
            // this is the case when a top has more than three edges — the
            // faces and edges concerned are found (OCCT L2393-2467).
            let mut trouve;
            let mut possible1 = false;
            let mut possible2 = false;
            trouve = false;
            nb = 0;
            let ev0 = topexp_vertices(cv1.arc());
            if ev0.0.is_same(&vtx) || ev0.1.is_same(&vtx) {
                possible1 = true;
            }
            let ev2 = topexp_vertices(cv2.arc());
            if ev2.0.is_same(&vtx) || ev2.1.is_same(&vtx) {
                possible2 = true;
            }
            if (possible1 && possible2) || (!possible1 && !possible2) || nbarete > 4 {
                while !trouve {
                    nb += 1;
                    if nb >= nn {
                        panic!("Standard_Failure: IntersectionAtEnd : the max number of faces reached");
                    }
                    if nb != 1 {
                        f3 = face_arr[nb - 2].clone();
                    }
                    face_arr[nb - 1] = f3.clone();
                    if cv1.arc().is_same(&edgelibre1) {
                        let mut fl = face_arr[nb - 1].clone();
                        cherche_face(&brep, &vf_list, &edgelibre2, &f1, &f2, &f3, &mut fl);
                        face_arr[nb - 1] = fl;
                    } else if cv1.arc().is_same(&edgelibre2) {
                        let mut fl = face_arr[nb - 1].clone();
                        cherche_face(&brep, &vf_list, &edgelibre1, &f1, &f2, &f3, &mut fl);
                        face_arr[nb - 1] = fl;
                    } else {
                        let mut fl = face_arr[nb - 1].clone();
                        cherche_face(&brep, &vf_list, &edge[nb - 1], &f1, &f2, &f3, &mut fl);
                        face_arr[nb - 1] = fl;
                    }
                    // OCCT passes tabedg(0, nn); its live entries mirror
                    // the Edge array bookkeeping.
                    let mut tabedg: Vec<Shape> = (0..=nn).map(|_| Shape::null()).collect();
                    tabedg[0] = edge[0].clone();
                    tabedg[nbface] = edge[nbface].clone();
                    let mut vl = Shape::null();
                    let mut el = Shape::null();
                    chfi3d_cherche_edge(&brep, &vtx, &tabedg, &face_arr[nb - 1], &mut el, &mut vl);
                    edge[nb] = el;
                    if edge[nb].is_same(cv2.arc()) {
                        trouve = true;
                    }
                }
                nbface = nb;
            } else {
                // OCCT L2462-2467: IntersectMoreCorner(Index); return.
                drop(fd);
                drop(st);
                this.perform_more_three_corner(index, 1);
                return;
            }
            if nbarete == 4 {
                // if two consecutive edges are G1 there is only one face of
                // intersection (OCCT L2469-2503).
                let mut ang1 = 0.0f64;
                let mut vcom = Shape::null();
                let mut trouve_v = false;
                chfi3d_cherche_vertex(&edge[0], &edge[1], &mut vcom, &mut trouve_v);
                if vcom.is_same(&vtx) {
                    ang1 = chfi3d_angle_edge(&brep, &vtx, &edge[0], &edge[1]);
                }
                if (ang1 - pi).abs() < 0.01 {
                    oneintersection1 = true;
                    facesau = face_arr[0].clone();
                    edgesau = edge[1].clone();
                    face_arr[0] = face_arr[1].clone();
                    edge[1] = edge[2].clone();
                    nbface = 1;
                }

                if !oneintersection1 {
                    let mut vcom = Shape::null();
                    let mut trouve_v = false;
                    chfi3d_cherche_vertex(&edge[1], &edge[2], &mut vcom, &mut trouve_v);
                    if vcom.is_same(&vtx) {
                        ang1 = chfi3d_angle_edge(&brep, &vtx, &edge[1], &edge[2]);
                    }
                    if (ang1 - pi).abs() < 0.01 {
                        oneintersection2 = true;
                        facesau = face_arr[1].clone();
                        edgesau = edge[1].clone();
                        edge[1] = edge[2].clone();
                        nbface = 1;
                    }
                }
            } else if nbarete == 5 {
                // pro15368 (OCCT L2506-2530): Modified by Sergey KHROMOV.
                let is_tangent0 = is_tangent_faces(
                    &brep,
                    &edge[0],
                    &f1,
                    &face_arr[0],
                    crate::geomalgo::gtests_stubs::GeomAbsShape::G1,
                );
                let is_tangent1 = is_tangent_faces(
                    &brep,
                    &edge[1],
                    &face_arr[0],
                    &face_arr[1],
                    crate::geomalgo::gtests_stubs::GeomAbsShape::G1,
                );
                let is_tangent2 = is_tangent_faces(
                    &brep,
                    &edge[2],
                    &face_arr[1],
                    &face_arr[2],
                    crate::geomalgo::gtests_stubs::GeomAbsShape::G1,
                );
                if (is_tangent0 || is_tangent2) && is_tangent1 {
                    facesau = face_arr[0].clone();
                    edgesau = edge[0].clone();
                    nbface = 1;
                    edge[1] = edge[3].clone();
                    face_arr[0] = face_arr[2].clone();
                    oneintersection1 = true;
                }
            }
        }
    } else {
        // onecorner (OCCT L2381-2386).
        nbface = 1;
        face_arr[0] = face.clone();
        edge[1] = edge[2].clone();
    }

    // OCCT L2390-2404: loop declarations.
    let mut pardeb = [0.0f64; 4];
    let mut parfin = [0.0f64; 4];
    let mut pfil1 = DVec2::ZERO;
    let mut pfac1 = DVec2::ZERO;
    let mut pfil2 = DVec2::ZERO;
    let mut pfac2 = DVec2::ZERO;
    let mut pint = DVec2::ZERO;
    let mut pfildeb = DVec2::ZERO;
    let mut hc1: Option<rcad_kernel::geom::Curve2d> = None;
    let mut hc2: Option<rcad_kernel::geom::Curve2d> = None;
    let mut proledge = vec![0i32; nn + 1];
    let mut prolface = vec![0i32; nn + 2]; // last prolface[nn] is for Fd
    let mut shrink = vec![0i32; nn + 1];
    let mut faceprol_surf: Vec<Option<Surface3>> = (0..=nn).map(|_| None).collect();
    let mut indcurve = vec![0i32; nn + 1];
    let mut indpoint1 = 0i32;
    let mut indpoint2 = 0i32;
    let mut interfedge: Vec<Option<super::chfi3d_ds::TopOpeBRepDSInterference>> =
        (0..=nn).map(|_| None).collect();
    let mut interf_pc: Vec<Option<super::chfi3d_ds::TopOpeBRepDSInterference>> =
        (0..nn).map(|_| None).collect();
    let mut interf_ps: Vec<Option<super::chfi3d_ds::TopOpeBRepDSInterference>> =
        (0..nn).map(|_| None).collect();
    let mut paredge1;
    let mut paredge2 = 0.0f64;
    let mut tolex = 1.0e-4;
    let mut extend = false;
    let mut sfacemoins1: Option<Surface3> = None;
    let mut sface: Option<Surface3> = None;
    let mut cint: Option<Curve3> = None;
    let mut c2dint1: Option<rcad_kernel::geom::Curve2d> = None;
    let mut c2dint2: Option<rcad_kernel::geom::Curve2d> = None;
    let mut hgs = hgs;
    let mut bs: Option<BRepAdaptorSurface> = None;
    let mut ps: Option<rcad_kernel::geom::Curve2d> = None;
    let mut pc: Option<rcad_kernel::geom::Curve2d> = None;
    let mut interfp1: Option<super::chfi3d_ds::TopOpeBRepDSInterference> = None;
    let mut interfp2: Option<super::chfi3d_ds::TopOpeBRepDSInterference> = None;
    let mut interfc: Option<super::chfi3d_ds::TopOpeBRepDSInterference> = None;
    let mut dist = 0.0f64;

    /***************************************************************************/
    // calculate intersection of the fillet and each face and storage in
    // the DS (OCCT L2406-2432).
    /***************************************************************************/
    for nb0 in 1..=nbface {
        prolface[nb0 - 1] = 0;
        proledge[nb0 - 1] = 0;
        shrink[nb0 - 1] = 0;
    }
    proledge[nbface] = 0;
    prolface[nn] = 0;
    if oneintersection1 || oneintersection2 {
        // OCCT: faceprol[1] = facesau (the saved face of the
        // one-intersection shortcuts).
        faceprol_surf[1] = if facesau.is_null() {
            None
        } else {
            brep.face_surface_world(&facesau)
        };
    }
    if !is_on_same1 && !is_on_same2 {
        check_shrink = false;
    }
    // in OnSame situation we need intersect Fd with Edge[0] or Edge[nbface]
    // as well.
    nb = if is_on_same1 { 0 } else { 1 };
    let mut inters_on_same_failed = false;

    // OCCT L2433: for (; nb <= nbface; nb++).
    loop {
        if nb > nbface {
            break;
        }
        extend = false;
        let e2 = edge[nb].clone();
        let f = if nb == 0 {
            f1.clone()
        } else {
            let fcur = face_arr[nb - 1].clone();
            if prolface[nb - 1] == 0 {
                // OCCT: faceprol[nb-1] = F (the face keeps its identity; the
                // extended surface is tracked in faceprol_surf — GAP note:
                // the OCCT MakeFace copy of the extended surface into a new
                // TopoDS_Face cannot be materialized in the rcad BRep
                // without a face-replacement API).
                faceprol_surf[nb - 1] = sfacemoins1.clone();
            }
            fcur
        };

        if f.is_null() {
            panic!("Standard_NullObject: IntersectionAtEnd : Trying to intersect with NULL face");
        }

        sfacemoins1 = brep.face_surface_world(&f);

        // determine intersections of edges and the fillet to find
        // limitations of intersections face - fillet (OCCT L2455-2475).
        if nb == 1 {
            hc1 = brep
                .curve_on_surface(&edge[0], &face_arr[0])
                .map(|(c, _, _)| c);
            if is_on_same1 {
                // update interference param on Fi1 and point of CV1
                // (OCCT L2477-2533).
                if prolface[0] != 0 {
                    bs = Some(match &faceprol_surf[0] {
                        Some(s) => BRepAdaptorSurface::initialize_surface(s.clone()),
                        None => BRepAdaptorSurface::initialize(&brep, &face_arr[0]),
                    });
                } else {
                    bs = Some(BRepAdaptorSurface::initialize(&brep, &face_arr[0]));
                }
                let c3df = this
                    .my_ds
                    .as_ref()
                    .expect("DS")
                    .curve(fi1.line_index())
                    .curve
                    .clone();
                let mut ufi = fi2.parameter(isfirst);
                let mut fi = fi1.clone();
                let mut pfac1_out = pfac1;
                let bs_r = bs.as_mut().expect("Bs");
                if !inters_update_on_same(
                    &hgs,
                    bs_r,
                    c3df.as_ref().expect("c3df"),
                    fi1.parameter_first(),
                    fi1.parameter_last(),
                    &f1,
                    &face_arr[0],
                    &edge[0],
                    &vtx,
                    isfirst,
                    10.0 * tolapp3d, // in
                    &mut fi,
                    &mut cv1,
                    &mut pfac1_out,
                    &mut ufi,
                    &brep,
                ) {
                    panic!("Standard_Failure: IntersectionAtEnd: pb intersection Face - Fi");
                }
                fi1 = fi;
                pfac1 = pfac1_out;
                if inters_on_same_failed {
                    // probable at fillet building — look for paredge2
                    // (OCCT L2523-2533).
                    let proj_target = pfac1;
                    paredge2 = match (&c2dint2, &hc1) {
                        (Some(c2), _) => project_point_on_curve2d(proj_target, c2),
                        (None, Some(c1)) => project_point_on_curve2d(proj_target, c1),
                        _ => 0.0,
                    };
                }
                // update stripe point (OCCT L2535-2538).
                let tpoint = TopOpeBRepDSPoint::new(cv1.point(), tolapp3d);
                indpoint1 = this.my_ds.as_mut().expect("DS").add_point(tpoint);
                st.set_index_point(indpoint1, isfirst, 1);
                // reset arc of CV1 (OCCT L2540-2545).
                let (vert1, vert2) = topexp_vertices(&edge[0]);
                let arc_ori = if vtx.is_same(&vert1) {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                cv1.set_arc(tolapp3d, edge[0].clone(), paredge2, arc_ori);
                *fd.change_vertex(isfirst, 1) = cv1.clone();
            } else if hc1.is_none() {
                // curve 2d not found: Sfacemoins1 is extended and projection
                // is done there (OCCT L2546-2572).
                let mut s_ext = sfacemoins1.clone().expect("face surface");
                let mut prol0 = 0i32;
                super::chfi3d_builder_c1::chfi3d_extend_surface(&mut s_ext, &mut prol0);
                prolface[0] = prol0;
                sfacemoins1 = Some(s_ext);
                if prolface[0] != 0 {
                    extend = true;
                    // OCCT: BRE.MakeFace(faceprol[0], Sfacemoins1, F.Location(), tol).
                    faceprol_surf[0] = sfacemoins1.clone();
                    if !is_on_same1 {
                        // Extrema_ExtPS ext(CV1.Point(), Asurf, tol, tol, MIN).
                        let s_ext_ref = sfacemoins1.clone().expect("s_ext");
                        let [uf0, ul0, vf0, vl0] = s_ext_ref.default_domain();
                        let tol = super::chfi3d_builder_0::P_CONFUSION;
                        let ext = rcad_kernel::base::extrema::ExtPS::with_domain(
                            cv1.point(),
                            &s_ext_ref,
                            uf0,
                            ul0,
                            vf0,
                            vl0,
                            tol,
                            tol,
                        );
                        if ext.nb_ext() > 0 {
                            let pon = ext.point(1);
                            pfac1 = DVec2::new(pon.u, pon.v);
                        }
                    }
                }
            } else {
                pfac1 = hc1
                    .as_ref()
                    .expect("Hc1")
                    .point_at(cv1.parameter_on_arc());
            }
            paredge1 = cv1.parameter_on_arc();
            if fi1.line_index() != 0 {
                pfil1 = fi1
                    .pcurve_on_surf()
                    .expect("pfil1")
                    .point_at(fi1.parameter(isfirst));
            } else {
                pfil1 = fi1
                    .pcurve_on_surf()
                    .expect("pfil1")
                    .point_at(fi1.parameter(!isfirst));
            }
            pfildeb = pfil1;
        } else {
            // OCCT L2578-2582.
            pfil1 = pfil2;
            paredge1 = paredge2;
            pfac1 = pint;
        }

        if nb != nbface || is_on_same2 {
            // intersect the limiting edge with the fillet surface
            // (OCCT L2585-2640).
            let mut nbp;
            let (c_edge, range) = brep
                .edge_curve_world(&e2)
                .expect("Standard_NullObject: edge curve");
            // OCCT: Ctrim = new Geom_TrimmedCurve(C, Ubid, Vbid); Utrim/Vtrim
            // from the basis curve; the restricted domain is adjusted for
            // the periodic case.
            let (utrim, vtrim) = {
                let [f0, l0] = c_edge.default_domain();
                (f0, l0)
            };
            let mut ubid = range[0];
            let mut vbid = range[1];
            if c_edge.is_periodic() {
                let period = vtrim - utrim;
                if ubid > period {
                    ubid = (utrim + vtrim) / 2.0;
                    vbid = vtrim;
                } else {
                    ubid = utrim;
                    vbid = (utrim + vtrim) / 2.0;
                }
            } else {
                ubid = utrim;
                vbid = vtrim;
            }
            // OCCT: inters.Perform(HC, HGs) over GeomAdaptor_Curve(C, Ubid, Vbid).
            let mut inters = crate::geomalgo::int_patch::int_cs::IntCurveSurface::new();
            let perform_curve = cint.clone().unwrap_or_else(|| c_edge.clone());
            let mut perform_domain = [ubid, vbid];
            let mut inters_done;
            {
                inters.perform(
                    &perform_curve,
                    hgs.surface(),
                    perform_domain,
                );
                inters_done = inters.is_done() && inters.nb_points() != 0;
            }
            // OCCT L2592-2637: extend surface of conge when nothing found.
            if prolface[nn] == 0 && !inters_done {
                let s1 = this.my_ds.as_ref().expect("DS").surface(fd.surf()).surface.clone();
                // OCCT down_cast<Geom_BoundedSurface> — BSpline/Bezier only
                // (GeomLib extension); other bounded kinds are GAP.
                let mut s1v = s1;
                let extended = geom_lib_extend_surf_by_length(
                    &mut s1v,
                    0.5 * (fi1_length.max(fi2_length)),
                    1,
                    false,
                    !isfirst,
                );
                if extended {
                    prolface[nn] = 1;
                    if st.is_in_ds(!isfirst) == 0 {
                        hgs = GeomAdaptorSurface::new(s1v.clone());
                        let mut inters2 = crate::geomalgo::int_patch::int_cs::IntCurveSurface::new();
                        inters2.perform(&perform_curve, hgs.surface(), perform_domain);
                        if inters2.is_done() && inters2.nb_points() != 0 {
                            let stoler = this
                                .my_ds
                                .as_ref()
                                .expect("DS")
                                .surface(isurf)
                                .tolerance();
                            let new_isurf = this
                                .my_ds
                                .as_mut()
                                .expect("DS")
                                .add_surface(TopOpeBRepDSSurface::new(s1v.clone(), stoler));
                            fd.change_surf(new_isurf);
                            // update history (OCCT L2617-2634).
                            let key = edge_spine.ptr_id();
                            if let Some(li) = this.my_evi_map.get_mut(&key) {
                                if let Some(pos) = li.iter().position(|&x| x == isurf) {
                                    li.remove(pos);
                                }
                                li.push(fd.surf());
                            } else {
                                this.my_evi_map.insert(key, vec![fd.surf()]);
                            }
                            isurf = fd.surf();
                            inters = inters2;
                            inters_done = true;
                        }
                    }
                }
            }
            // OCCT L2639-2672: the BSpline/Bezier edge-curve fallback —
            // extend the support faces and intersect them (GeomInt_IntSS).
            if !inters_done {
                let is_bezier_or_bspline = matches!(
                    &c_edge,
                    rcad_kernel::geom::Curve3::BSpline(_) | rcad_kernel::geom::Curve3::Bezier(_)
                );
                if is_bezier_or_bspline {
                    // Sface = BRep_Tool::Surface(Face[nb]); extend it.
                    let mut sface_v = brep.face_surface_world(&face_arr[nb]);
                    if let Some(sv) = sface_v.as_mut() {
                        let mut prol_nb = prolface[nb];
                        super::chfi3d_builder_c1::chfi3d_extend_surface(sv, &mut prol_nb);
                        prolface[nb] = prol_nb;
                    }
                    sface = sface_v.clone();
                    if nb != 0 && prolface[nb - 1] == 0 {
                        if let Some(sv) = sfacemoins1.as_mut() {
                            let mut prol_prev = 0i32;
                            super::chfi3d_builder_c1::chfi3d_extend_surface(sv, &mut prol_prev);
                            prolface[nb - 1] = prol_prev;
                            if prol_prev != 0 {
                                faceprol_surf[nb - 1] = sfacemoins1.clone();
                            }
                        }
                    } else if let Some(sv) = sfacemoins1.as_mut() {
                        let mut prol_tmp = 0i32;
                        super::chfi3d_builder_c1::chfi3d_extend_surface(sv, &mut prol_tmp);
                    }
                    // OCCT: GeomInt_IntSS InterSS(Sfacemoins1, Sface, 1.e-7,
                    // true, true, true).  GAP: the rcad kernel IntSS
                    // (GeomAPI_IntSS stand-in) returns 3D lines only —
                    // LineOnS1/LineOnS2 and TolReached3d are pending; the
                    // 3D line is consumed and the pcurves stay null (the
                    // dependent extend-branches keep the OCCT null-pcurve
                    // outcome).
                    if let (Some(sm1), Some(sf1)) = (sfacemoins1.clone(), sface.clone()) {
                        let mut inter_ss =
                            rcad_kernel::base::geom_api::int_ss::IntSS::with_surfaces(
                                &sm1, &sf1, 1.0e-7,
                            );
                        if inter_ss.is_done() {
                            let mut trouve = false;
                            for i in 1..=inter_ss.nb_lines() {
                                if trouve {
                                    break;
                                }
                                extend = true;
                                cint = Some(inter_ss.line(i).clone());
                                // GAP: C2dint1/C2dint2 stay null (see above).
                                perform_domain = c_edge.default_domain()[..2]
                                    .try_into()
                                    .unwrap_or([ubid, vbid]);
                                let mut inters3 =
                                    crate::geomalgo::int_patch::int_cs::IntCurveSurface::new();
                                inters3.perform(&perform_curve, hgs.surface(), perform_domain);
                                inters = inters3;
                                trouve = inters.is_done() && inters.nb_points() != 0;
                                // OCCT: tolex = InterSS.TolReached3d() — GAP,
                                // the InterSS tolerance stands in.
                                tolex = 1.0e-7;
                            }
                            if trouve {
                                inters_done = true;
                            }
                        }
                    }
                }
            }
            // OCCT L2674-2725.
            if inters.is_done() || inters_done {
                nbp = inters.nb_points();
                if nbp == 0 {
                    if nb == 0 || nb == nbface {
                        inters_on_same_failed = true;
                    } else {
                        // OCCT: PerformMoreThreeCorner(Index, 1); return.
                        drop(fd);
                        drop(st);
                        this.perform_more_three_corner(index, 1);
                        return;
                    }
                } else {
                    let p_vtx = brep.vertex_position(&vtx);
                    let mut distmin = p_vtx.distance(inters.point(1).pnt());
                    nbp = 1;
                    for i in 2..=inters.nb_points() {
                        dist = p_vtx.distance(inters.point(i).pnt());
                        if dist < distmin {
                            distmin = dist;
                            nbp = i;
                        }
                    }
                    let pt = inters.point(nbp);
                    pfil2 = DVec2::new(pt.u(), pt.v());
                    paredge2 = pt.w();
                    if !extend {
                        // pcurves of the edge on both faces (OCCT L2700-2712).
                        let Some((cfm1, _, _)) = brep.curve_on_surface(&e2, &f) else {
                            panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                        };
                        let Some((cf1, _, _)) = brep.curve_on_surface(&e2, &face_arr[nb]) else {
                            panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                        };
                        pfac2 = cfm1.point_at(paredge2);
                        pint = cf1.point_at(paredge2);
                    } else if c2dint1.is_none() || c2dint2.is_none() {
                        panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                    } else {
                        pfac2 = c2dint1
                            .as_ref()
                            .expect("C2dint1")
                            .point_at(paredge2);
                        pint = c2dint2.as_ref().expect("C2dint2").point_at(paredge2);
                    }
                }
            } else {
                panic!("Standard_Failure: IntersectionAtEnd: pb intersection Face cb");
            }
        } else {
            // OCCT L2727-2790 — the last face: Hc2 branch.
            hc2 = brep
                .curve_on_surface(&e2, &face_arr[nbface - 1])
                .map(|(c, _, _)| c);
            if hc2.is_none() {
                // curve 2d is not found: Sfacemoins1 is extended and
                // CV2.Point() is projected there.
                let mut s_ext = sfacemoins1.clone().expect("face surface");
                let mut prol0 = prolface[0];
                super::chfi3d_builder_c1::chfi3d_extend_surface(&mut s_ext, &mut prol0);
                prolface[0] = prol0;
                sfacemoins1 = Some(s_ext);
                if prolface[0] != 0 {
                    extend = true;
                    faceprol_surf[nb - 1] = sfacemoins1.clone();
                    let s_ext_ref = sfacemoins1.clone().expect("s_ext");
                    let [uf0, ul0, vf0, vl0] = s_ext_ref.default_domain();
                    let tol = super::chfi3d_builder_0::P_CONFUSION;
                    let ext = rcad_kernel::base::extrema::ExtPS::with_domain(
                        cv2.point(),
                        &s_ext_ref,
                        uf0,
                        ul0,
                        vf0,
                        vl0,
                        tol,
                        tol,
                    );
                    if ext.nb_ext() > 0 {
                        let pon = ext.point(1);
                        pfac2 = DVec2::new(pon.u, pon.v);
                    }
                }
            } else {
                pfac2 = hc2
                    .as_ref()
                    .expect("Hc2")
                    .point_at(cv2.parameter_on_arc());
            }
            paredge2 = cv2.parameter_on_arc();
            if fi2.line_index() != 0 {
                pfil2 = fi2
                    .pcurve_on_surf()
                    .expect("pfil2")
                    .point_at(fi2.parameter(isfirst));
            } else {
                pfil2 = fi2
                    .pcurve_on_surf()
                    .expect("pfil2")
                    .point_at(fi2.parameter(!isfirst));
            }
        }
        if nb == 0 {
            // found paredge1 on Edge[0] in OnSame situation on F1.
            nb += 1;
            continue;
        }

        if nb == nbface && is_on_same2 {
            // update interference param on Fi2 and point of CV2
            // (OCCT L2800-2852).
            if prolface[nb - 1] != 0 {
                bs = Some(match &faceprol_surf[nb - 1] {
                    Some(s) => BRepAdaptorSurface::initialize_surface(s.clone()),
                    None => BRepAdaptorSurface::initialize(&brep, &face_arr[nb - 1]),
                });
            } else {
                bs = Some(BRepAdaptorSurface::initialize(&brep, &face_arr[nb - 1]));
            }
            let c3df = this
                .my_ds
                .as_ref()
                .expect("DS")
                .curve(fi2.line_index())
                .curve
                .clone();
            let mut ufi = fi1.parameter(isfirst);
            let mut fi = fi2.clone();
            let mut pfac2_out = pfac2;
            let bs_r = bs.as_mut().expect("Bs");
            if !inters_update_on_same(
                &hgs,
                bs_r,
                c3df.as_ref().expect("c3df"),
                fi2.parameter_first(),
                fi2.parameter_last(),
                &f2,
                &f,
                &edge[nb],
                &vtx,
                isfirst,
                10.0 * tolapp3d, // in
                &mut fi,
                &mut cv2,
                &mut pfac2_out,
                &mut ufi,
                &brep,
            ) {
                panic!("Standard_Failure: IntersectionAtEnd: pb intersection Face - Fi");
            }
            fi2 = fi;
            pfac2 = pfac2_out;
            if inters_on_same_failed {
                // look for paredge2 (OCCT L2824-2833).
                let proj_target = pfac2;
                paredge2 = if extend {
                    match &c2dint2 {
                        Some(c2) => project_point_on_curve2d(proj_target, c2),
                        None => 0.0,
                    }
                } else {
                    match brep.curve_on_surface(&e2, &face_arr[nbface - 1]) {
                        Some((c, _, _)) => project_point_on_curve2d(proj_target, &c),
                        None => 0.0,
                    }
                };
            }
            // update stripe point (OCCT L2835-2839).
            let tpoint = TopOpeBRepDSPoint::new(cv2.point(), tolapp3d);
            indpoint2 = this.my_ds.as_mut().expect("DS").add_point(tpoint);
            st.set_index_point(indpoint2, isfirst, 2);
            // reset arc of CV2 (OCCT L2841-2846).
            let (vert1, vert2) = topexp_vertices(&edge[nbface]);
            let arc_ori = if vtx.is_same(&vert1) {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            cv2.set_arc(tolapp3d, edge[nbface].clone(), paredge2, arc_ori);
            *fd.change_vertex(isfirst, 2) = cv2.clone();
        }

        // OCCT L2848-2858: Bs.Initialize(faceprol[nb-1] / Face[nb-1]).
        if prolface[nb - 1] != 0 {
            bs = Some(match &faceprol_surf[nb - 1] {
                Some(s) => BRepAdaptorSurface::initialize_surface(s.clone()),
                None => BRepAdaptorSurface::initialize(&brep, &face_arr[nb - 1]),
            });
        } else {
            bs = Some(BRepAdaptorSurface::initialize(&brep, &face_arr[nb - 1]));
        }

        // offset of parameters if they are not in the same period — the
        // OCC354 comment block is elided in OCCT itself (OCCT L2861-2908).
        pardeb[0] = pfil1.x;
        pardeb[1] = pfil1.y;
        pardeb[2] = pfac1.x;
        pardeb[3] = pfac1.y;
        parfin[0] = pfil2.x;
        parfin[1] = pfil2.y;
        parfin[2] = pfac2.x;
        parfin[3] = pfac2.y;

        let (uu1, uu2, vv1, vv2) = chfi3d_boite(pfac1, pfac2);
        let bs_r = bs.as_mut().expect("Bs");
        chfi3d_bound_fac(bs_r, uu1, uu2, vv1, vv2, true);

        ////////////////////////////////////////////////////////////////////////
        // calculate intersections face - fillet (OCCT L2915-2930).
        ////////////////////////////////////////////////////////////////////////
        let bs_view = GeomAdaptorSurface::new(bs_r.surface.clone());
        let computed = chfi3d_compute_curves(
            &hgs,
            &bs_view,
            pardeb,
            parfin,
            tolapp3d,
            tol2d,
            &mut tolreached,
        );
        let (cc, pc_new, ps_new, tolreached_v) = match computed {
            Some(res) => (res.c3d, res.pc2, res.pc1, res.tolreached),
            None => {
                // OCCT L2923-2927: PerformMoreThreeCorner(Index, 1); return.
                drop(fd);
                drop(st);
                this.perform_more_three_corner(index, 1);
                return;
            }
        };
        tolreached = tolreached_v;
        pc = Some(pc_new);
        ps = Some(ps_new);

        // storage of information in the data structure — evaluate
        // tolerances (OCCT L2935-2960).
        let p1 = c3d_first(&cc);
        let p2 = c3d_last(&cc);
        let p2d_a;
        let p2d_b;
        {
            // Pc->D0(p1, p2d1); Pc->D0(p2, p2d2).
            let pc_r = pc.as_ref().expect("Pc");
            p2d_a = pc_r.point_at(p1);
            p2d_b = pc_r.point_at(p2);
        }
        let bs_r = bs.as_ref().expect("Bs");
        let pp1 = hgs.surface().point_at(pardeb[0], pardeb[1]);
        let pp2 = hgs.surface().point_at(parfin[0], parfin[1]);
        let pp3 = bs_r.surface.point_at(pardeb[2], pardeb[3]);
        let pp4 = bs_r.surface.point_at(parfin[2], parfin[3]);
        let pp7 = bs_r.surface.point_at(p2d_a.x, p2d_a.y);
        let pp8 = bs_r.surface.point_at(p2d_b.x, p2d_b.y);
        let ps_r = ps.as_ref().expect("Ps");
        let ps_a = ps_r.point_at(p1);
        let ps_b = ps_r.point_at(p2);
        let pp5 = hgs.surface().point_at(ps_a.x, ps_a.y);
        let pp6 = hgs.surface().point_at(ps_b.x, ps_b.y);
        let to1 = (pp1.distance(pp5) + pp3.distance(pp7)).max(tolreached);
        let to2 = (pp2.distance(pp6) + pp4.distance(pp8)).max(tolreached);

        //////////////////////////////////////////////////////////////////////
        // storage in the DS of the intersection curve (OCCT L2963-3016).
        //////////////////////////////////////////////////////////////////////
        let mut isvtx1 = false;
        let mut isvtx2 = false;

        if nb == 1 {
            indpoint1 = st.index_point(isfirst, 1);
            if !cv1.is_vertex() {
                let tpt = this.my_ds.as_mut().expect("DS").change_point(indpoint1);
                tpt.set_tolerance(tpt.tolerance().max(to1));
            } else {
                isvtx1 = true;
            }
        }
        if nb == nbface {
            indpoint2 = st.index_point(isfirst, 2);
            if !cv2.is_vertex() {
                let tpt = this.my_ds.as_mut().expect("DS").change_point(indpoint2);
                tpt.set_tolerance(tpt.tolerance().max(to2));
            } else {
                isvtx2 = true;
            }
        } else {
            let point = cc.point_at(c3d_last(&cc));
            let tpoint = TopOpeBRepDSPoint::new(point, to2);
            indpoint2 = this.my_ds.as_mut().expect("DS").add_point(tpoint);
        }

        if nb != 1 {
            let tpt = this.my_ds.as_mut().expect("DS").change_point(indpoint1);
            tpt.set_tolerance(tpt.tolerance().max(to1));
        }
        indcurve[nb - 1] = this
            .my_ds
            .as_mut()
            .expect("DS")
            .add_curve(TopOpeBRepDSCurve::new(Some(cc.clone()), tolreached));

        interfp1 = Some(super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Forward,
            indcurve[nb - 1],
            indpoint1,
            c3d_first(&cc),
            isvtx1,
        ));
        interfp2 = Some(super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Reversed,
            indcurve[nb - 1],
            indpoint2,
            c3d_last(&cc),
            isvtx2,
        ));

        this.my_ds
            .as_mut()
            .expect("DS")
            .change_curve_interferences(indcurve[nb - 1])
            .push(interfp1.take().expect("Interfp1"));
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_curve_interferences(indcurve[nb - 1])
            .push(interfp2.take().expect("Interfp2"));

        //////////////////////////////////////////////////////////////////////
        // storage for the face (OCCT L3019-3050).
        //////////////////////////////////////////////////////////////////////
        let mut ori = Orientation::Forward;
        orface = face_arr[nb - 1].orientation;
        if orface == orsurfdata {
            orien = topabs_reverse(orcourbe);
        } else {
            orien = orcourbe;
        }
        // limitation of edges of faces.
        if nb == 1 {
            let iarc1 = this.my_ds.as_mut().expect("DS").add_shape(&edge[0]);
            let itf = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                cv1.transition_on_arc(),
                iarc1,
                indpoint1,
                paredge1,
                isvtx1,
            );
            interfedge[0] = Some(itf);
        }
        if nb == nbface {
            let iarc2 = this.my_ds.as_mut().expect("DS").add_shape(&edge[nb]);
            let itf = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                cv2.transition_on_arc(),
                iarc2,
                indpoint2,
                paredge2,
                isvtx2,
            );
            interfedge[nb] = Some(itf);
        }

        // OCCT L3052-3290: storage of the limiting curve on the next edge
        // (possibly extended).
        if nb != nbface || oneintersection1 || oneintersection2 {
            // orientation of the limitation (OCCT L3057-3103).
            if nface == 3 {
                let (ev1, ev2) = topexp_vertices(&edge[nb]);
                let o1 = contain_v(&brep, &f1, &ev1) || contain_v(&brep, &f2, &ev1);
                let o2 = contain_v(&brep, &f1, &ev2) || contain_v(&brep, &f2, &ev2);
                if o1 {
                    ori = Orientation::Forward;
                } else if o2 {
                    ori = Orientation::Reversed;
                } else {
                    panic!("Standard_Failure: IntersectionAtEnd : pb orientation");
                }
                if contain_v(&brep, &f1, &ev1) && contain_v(&brep, &f1, &ev2) {
                    let d1 = brep.vertex_position(&ev1).distance(brep.vertex_position(&vtx));
                    let d2 = brep.vertex_position(&ev2).distance(brep.vertex_position(&vtx));
                    ori = if d1 < d2 {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                }
                if contain_v(&brep, &f2, &ev1) && contain_v(&brep, &f2, &ev2) {
                    let d1 = brep.vertex_position(&ev1).distance(brep.vertex_position(&vtx));
                    let d2 = brep.vertex_position(&ev2).distance(brep.vertex_position(&vtx));
                    ori = if d1 < d2 {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                }
            } else {
                let (evf, _) = topexp_vertices(&edge[nb]);
                ori = if evf.is_same(&vtx) {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
            }
            if !extend && !(oneintersection1 || oneintersection2) {
                let iarc2 = this.my_ds.as_mut().expect("DS").add_shape(&edge[nb]);
                let itf = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                    ori,
                    iarc2,
                    indpoint2,
                    paredge2,
                    false,
                );
                interfedge[nb] = Some(itf);
            } else {
                if !(oneintersection1 || oneintersection2) {
                    proledge[nb] = 1;
                }
                // OCCT L3120-3231: the extended-curve limitation.
                let mut pext = brep.vertex_position(&vtx);
                // OCCT: cad.Load(cint / edgesau curve); Extrema_ExtPC.
                let csau: Option<Curve3> = if !(oneintersection1 || oneintersection2) {
                    cint.clone()
                } else {
                    brep.edge_curve_world(&edgesau).map(|(c, _)| c)
                };
                let csau = match csau {
                    Some(c) => c,
                    None => panic!("Standard_NullObject: IntersectionAtEnd extension curve"),
                };
                // Extrema_ExtPC ext(pext, cad, tolpt).
                let par1 = extrema_ext_pc_first_parameter(&csau, pext, tolpt)
                    .expect("StdFail_NotDone: IntersectionAtEnd ExtPC");
                let mut par_vtx = par1;
                let mut par2;
                let mut vtx1 = false;
                let mut vtx2 = false;
                let mut ind;
                if oneintersection1 || oneintersection2 {
                    if oneintersection2 {
                        pext = cv2.point();
                        ind = indpoint2;
                    } else {
                        pext = cv1.point();
                        ind = indpoint1;
                    }
                    par2 = extrema_ext_pc_first_parameter(&csau, pext, tolpt)
                        .expect("StdFail_NotDone: IntersectionAtEnd ExtPC2");
                } else {
                    par2 = paredge2;
                    ind = indpoint2;
                }
                let (mut par_lo, mut par_hi, indp1, indp2) = if par1 > par2 {
                    vtx2 = true;
                    (par2, par1, ind, this.my_ds.as_mut().expect("DS").add_shape(&vtx))
                } else {
                    vtx1 = true;
                    (
                        par1,
                        par2,
                        this.my_ds.as_mut().expect("DS").add_shape(&vtx),
                        ind,
                    )
                };
                let mut ct = rcad_kernel::geom::Curve3::Trimmed(
                    rcad_kernel::geom::TrimmedCurve3::new(csau.clone(), par_lo, par_hi),
                );
                // orientation of the trimmed limitation (OCCT L3164-3176).
                let c_end1 = csau.point_at(par_lo);
                let c_end2 = csau.point_at(par_hi);
                let cc_b1 = cc.point_at(cc.default_domain()[0]);
                let cc_b2 = cc.point_at(cc.default_domain()[1]);
                let mut orient = if c_end2.distance(cc_b1) < tolpt || c_end1.distance(cc_b2) < tolpt {
                    orien
                } else {
                    topabs_reverse(orien)
                };
                let mut c2d_use1: Option<rcad_kernel::geom::Curve2d> = None;
                let mut indice;
                if oneintersection1 || oneintersection2 {
                    indice = this.my_ds.as_mut().expect("DS").add_shape(&face_arr[0]);
                    if extend {
                        // OCCT: DStr.SetNewSurface(Face[0], Sfacemoins1) —
                        // GAP: the DS new-surface map is pending (D6); the
                        // face keeps its original surface payload.
                        // ComputeCurve2d(ct, faceprol[0], C2dint1).
                        let faceprol_face = if face_arr[0].is_null() {
                            Shape::null()
                        } else {
                            face_arr[0].clone()
                        };
                        let mut c2d_out = None;
                        compute_curve2d(&ct, &faceprol_face, &mut c2d_out, &brep);
                        c2d_use1 = c2d_out;
                    } else {
                        let mut a_local_edge = edgesau.clone();
                        if a_local_edge.orientation != orient {
                            a_local_edge.orientation = topabs_reverse(a_local_edge.orientation);
                        }
                        c2d_use1 = brep
                            .curve_on_surface(&a_local_edge, &face_arr[0])
                            .map(|(c, _, _)| c);
                    }
                } else {
                    indice = this.my_ds.as_mut().expect("DS").add_shape(&face_arr[nb - 1]);
                    // OCCT: DStr.SetNewSurface(Face[nb-1], Sfacemoins1) —
                    // GAP: DS new-surface map pending (D6).
                }
                //// for periodic 3d curves //// (OCCT L3188-3200): guarded
                // by !C2dint1.IsNull() — naturally skipped under the IntSS
                // pcurve GAP.
                if csau.is_periodic() && c2d_use1.is_some() {
                    // OCCT: P2d = BRep_Tool::Parameters(Vtx, Face[0]) —
                    // GAP: vertex-UV parameters pending; the periodic
                    // parameter shift keeps the curve parameters.
                    let _ = &c2d_use1;
                }

                ct = rcad_kernel::geom::Curve3::Trimmed(
                    rcad_kernel::geom::TrimmedCurve3::new(csau.clone(), par_lo, par_hi),
                );
                if oneintersection1 || oneintersection2 {
                    tolex = 10.0 * brep.tolerance(&edgesau);
                }
                // OCCT L3210-3221: extend → ChFi3d_EvalTolReached (needs the
                // IntSS pcurves — GAP: tolex is kept).
                let _ = &sface;
                let indcurv = this
                    .my_ds
                    .as_mut()
                    .expect("DS")
                    .add_curve(TopOpeBRepDSCurve::new(Some(ct.clone()), tolex));
                let itfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                    Orientation::Forward,
                    indcurv,
                    indp1,
                    par_lo,
                    vtx1,
                );
                let itfp2 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                    Orientation::Reversed,
                    indcurv,
                    indp2,
                    par_hi,
                    vtx2,
                );
                this.my_ds
                    .as_mut()
                    .expect("DS")
                    .change_curve_interferences(indcurv)
                    .push(itfp1);
                this.my_ds
                    .as_mut()
                    .expect("DS")
                    .change_curve_interferences(indcurv)
                    .push(itfp2);

                let interfc_v = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
                    indcurv,
                    indice,
                    c2d_use1.clone(),
                    orient,
                );
                interfc = Some(interfc_v);
                this.my_ds
                    .as_mut()
                    .expect("DS")
                    .change_shape_interferences(indice)
                    .push(interfc.clone().expect("Interfc"));
                if oneintersection1 || oneintersection2 {
                    indice = this.my_ds.as_mut().expect("DS").add_shape(&facesau);
                    if facesau.orientation == face_arr[0].orientation {
                        orient = topabs_reverse(orient);
                    }
                    let mut c2d_use2: Option<rcad_kernel::geom::Curve2d> = None;
                    if extend {
                        // OCCT: ComputeCurve2d(ct, faceprol[1], C2dint2).
                        let mut c2d_out = None;
                        compute_curve2d(&ct, &facesau, &mut c2d_out, &brep);
                        c2d_use2 = c2d_out;
                    } else {
                        let mut a_local_edge = edgesau.clone();
                        if a_local_edge.orientation != orient {
                            a_local_edge.orientation = topabs_reverse(a_local_edge.orientation);
                        }
                        c2d_use2 = brep
                            .curve_on_surface(&a_local_edge, &facesau)
                            .map(|(c, _, _)| c);
                    }
                    let interfc_v = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
                        indcurv,
                        indice,
                        c2d_use2,
                        orient,
                    );
                    interfc = Some(interfc_v);
                    if !bordlibre {
                        this.my_ds
                            .as_mut()
                            .expect("DS")
                            .change_shape_interferences(indice)
                            .push(interfc.expect("Interfc"));
                    }
                } else {
                    indice = this.my_ds.as_mut().expect("DS").add_shape(&face_arr[nb]);
                    // OCCT: DStr.SetNewSurface(Face[nb], Sface) — GAP (D6).
                    if face_arr[nb].orientation == face_arr[nb - 1].orientation {
                        orient = topabs_reverse(orient);
                    }
                    let interfc_v = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
                        indcurv,
                        indice,
                        c2dint2.clone(),
                        orient,
                    );
                    interfc = Some(interfc_v);
                    if !bordlibre {
                        this.my_ds
                            .as_mut()
                            .expect("DS")
                            .change_shape_interferences(indice)
                            .push(interfc.expect("Interfc"));
                    }
                }
                let _ = par_vtx;
                let _ = &mut par_lo;
                let _ = &mut par_hi;
            }
        }

        // OCCT L3293-3352: the shrink detection and storage.
        let p2d_a;
        let p2d_b;
        {
            let pc_r = pc.as_ref().expect("Pc");
            p2d_a = pc_r.point_at(p1);
            p2d_b = pc_r.point_at(p2);
        }
        if check_shrink
            && chfi3d_is_shrink(
                ps.as_ref().expect("Ps"),
                p1,
                p2,
                check_shr_param,
                is_u_shrink,
                super::chfi3d_builder_0::P_CONFUSION,
            )
        {
            shrink[nb - 1] = 1;
            // store section face-chamf curve for previous SurfData
            // (OCCT L3300-3320).
            if !is_shrink {
                // first time
                if let Some(sdprev) = &sdprev_opt {
                    let fi = sdprev.read().expect("surfdata lock").interference_on_s1().clone();
                    let uv = fi
                        .pcurve_on_surf()
                        .expect("prevSD Fi PCurveOnSurf")
                        .point_at(fi.parameter(isfirst));
                    prev_sd_param = if is_u_shrink { uv.x } else { uv.y };
                }
            }
            let mut uv1 = p2d_a;
            let mut uv2 = p2d_b;
            if is_u_shrink {
                uv1.x = prev_sd_param;
                uv2.x = prev_sd_param;
            } else {
                uv1.y = prev_sd_param;
                uv2.y = prev_sd_param;
            }
            // OCCT: ChFi3d_ComputePCurv(Cc, UV1, UV2, Ps, prev surf, p1, p2,
            // tolapp3d, aTolreached).
            if let Some(sdprev) = &sdprev_opt {
                let prev_surf = this
                    .my_ds
                    .as_ref()
                    .expect("DS")
                    .surface(sdprev.read().expect("surfdata lock").surf())
                    .surface
                    .clone();
                let ps_new = chfi3d_compute_pcurv_2pt(uv1, uv2, p1, p2, false);
                ps = Some(ps_new);
                let atolreached = tolapp3d;
                let tcurv = this.my_ds.as_mut().expect("DS").change_curve(indcurve[nb - 1]);
                tcurv.set_tolerance(tcurv.tolerance().max(atolreached));

                let itfps = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
                    indcurve[nb - 1],
                    isurf_prev,
                    ps.clone(),
                    orcourbe,
                );
                interf_ps[nb - 1] = Some(itfps);
                this.my_ds
                    .as_mut()
                    .expect("DS")
                    .change_surface_interferences(isurf_prev)
                    .push(interf_ps[nb - 1].clone().expect("InterfPS"));

                if is_on_same2 {
                    mid_p2d = p2d_b;
                    mid_ipoint = indpoint2;
                } else if !is_shrink {
                    mid_p2d = p2d_a;
                    mid_ipoint = indpoint1;
                }
                is_shrink = true;
            }
        }

        // OCCT L3355-3360: the face pcurve interference + loop tail.
        let indice = this.my_ds.as_mut().expect("DS").add_shape(&face_arr[nb - 1]);
        let itfpc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
            indcurve[nb - 1],
            indice,
            pc.clone(),
            orien,
        );
        interf_pc[nb - 1] = Some(itfpc);
        if shrink[nb - 1] == 0 {
            let itfps = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
                indcurve[nb - 1],
                isurf,
                ps.clone(),
                orcourbe,
            );
            interf_ps[nb - 1] = Some(itfps);
        }
        indpoint1 = indpoint2;

        nb += 1;
    } // end loop on faces being intersected with ChFi

    // OCCT L3601-3604.
    if is_on_same1 {
        cv1.reset();
        *fd.change_vertex(isfirst, 1) = cv1.clone();
    }
    if is_on_same2 {
        cv2.reset();
        *fd.change_vertex(isfirst, 2) = cv2.clone();
    }

    // OCCT L3606-3624: final storage loop.
    for nb0 in 1..=nbface {
        let indice = this.my_ds.as_mut().expect("DS").add_shape(&face_arr[nb0 - 1]);
        if let Some(itfpc) = interf_pc[nb0 - 1].take() {
            this.my_ds
                .as_mut()
                .expect("DS")
                .change_shape_interferences(indice)
                .push(itfpc);
        }
        if shrink[nb0 - 1] == 0 {
            if let Some(itfps) = interf_ps[nb0 - 1].take() {
                this.my_ds
                    .as_mut()
                    .expect("DS")
                    .change_surface_interferences(isurf)
                    .push(itfps);
            }
        }
        if proledge[nb0 - 1] == 0 {
            if let Some(itfe) = interfedge[nb0 - 1].take() {
                // OCCT: DStr.ChangeShapeInterferences(Edge[nb-1]).Append(...).
                let ish = this.my_ds.as_mut().expect("DS").add_shape(&edge[nb0 - 1]);
                this.my_ds
                    .as_mut()
                    .expect("DS")
                    .change_shape_interferences(ish)
                    .push(itfe);
            }
        }
    }
    if let Some(itfe) = interfedge[nbface].take() {
        // OCCT: DStr.ChangeShapeInterferences(Edge[nbface]).Append(...).
        let ish = this.my_ds.as_mut().expect("DS").add_shape(&edge[nbface]);
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_shape_interferences(ish)
            .push(itfe);
    }

    // OCCT L3627-3760: shrink tail / InDS.
    if !is_shrink {
        st.in_ds(isfirst, 1);
    } else {
        // compute curves for the !<isfirst> end of <Fd> and the <isfirst>
        // end of previous <SurfData> (OCCT L3637-3660).
        let mut uv1 = mid_p2d;
        let mut uv2 = mid_p2d;
        let uv;
        if is_on_same1 {
            let v = fi1
                .pcurve_on_surf()
                .expect("Fi1 PCurveOnSurf")
                .point_at(fi1.parameter(!isfirst));
            uv = v;
            uv2 = v;
        } else {
            let v = fi2
                .pcurve_on_surf()
                .expect("Fi2 PCurveOnSurf")
                .point_at(fi2.parameter(!isfirst));
            uv = uv1;
            uv1 = v;
        }
        let a_surf = this
            .my_ds
            .as_ref()
            .expect("DS")
            .surface(fd.surf())
            .surface
            .clone();
        // ChFi3d_ComputeArete for Fd.
        let (c3d, ps_new, parf, parl, atolreached) = chfi3d_compute_arete(
            &brep,
            &cv1,
            uv1,
            &cv2,
            uv2,
            &a_surf,
            tolapp3d,
            tol2d,
            0,
        );
        let c3d = c3d.expect("ChFi3d_ComputeArete C3d");
        ps = Some(ps_new);
        let p1 = parf;
        let p2 = parl;
        let _ = uv;

        indpoint1 = mid_ipoint;
        indpoint2 = mid_ipoint;
        if is_on_same1 {
            let point = c3d.point_at(p2);
            let tpoint = TopOpeBRepDSPoint::new(point, atolreached);
            indpoint2 = this.my_ds.as_mut().expect("DS").add_point(tpoint);
        } else {
            let point = c3d.point_at(p1);
            let tpoint = TopOpeBRepDSPoint::new(point, atolreached);
            indpoint1 = this.my_ds.as_mut().expect("DS").add_point(tpoint);
        }

        // for SDprev (OCCT L3683-3710).
        let icurv = this
            .my_ds
            .as_mut()
            .expect("DS")
            .add_curve(TopOpeBRepDSCurve::new(Some(c3d.clone()), atolreached));
        let itfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Forward,
            icurv,
            indpoint1,
            p1,
            false,
        );
        let itfp2 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Reversed,
            icurv,
            indpoint2,
            p2,
            false,
        );
        let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
            icurv,
            isurf,
            ps.clone(),
            orcourbe,
        );
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_curve_interferences(icurv)
            .push(itfp1);
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_curve_interferences(icurv)
            .push(itfp2);
        this.my_ds
            .as_mut()
            .expect("DS")
            .change_surface_interferences(isurf)
            .push(itfc);

        if let Some(sdprev) = &sdprev_opt {
            let a_surf_prev = this
                .my_ds
                .as_ref()
                .expect("DS")
                .surface(sdprev.read().expect("surfdata lock").surf())
                .surface
                .clone();
            let mut uv1p = uv1;
            let mut uv2p = uv2;
            if is_u_shrink {
                uv1p.x = prev_sd_param;
                uv2p.x = prev_sd_param;
            } else {
                uv1p.y = prev_sd_param;
                uv2p.y = prev_sd_param;
            }
            let pc_new = chfi3d_compute_pcurv_2pt(uv1p, uv2p, p1, p2, false);
            let itfc = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
                icurv,
                isurf_prev,
                Some(pc_new),
                topabs_reverse(orcourbe),
            );
            this.my_ds
                .as_mut()
                .expect("DS")
                .change_surface_interferences(isurf_prev)
                .push(itfc);
        }

        // to process properly this case in ChFi3d_FilDS() (OCCT L3712-3730).
        st.in_ds(isfirst, 2);
        let fi_idx = if is_on_same1 { 2 } else { 1 };
        fd.change_interference(fi_idx).lineindex = 0;
        if let Some(sdprev) = &sdprev_opt {
            let mut sdprev_w = sdprev.write().expect("surfdata lock");
            let cp_prev1 = sdprev_w.change_vertex(isfirst, fi_idx).clone();
            let cp_last1 = fd.change_vertex(isfirst, fi_idx).clone();
            let mut cp_last2 = fd.change_vertex(!isfirst, fi_idx).clone();
            if cp_prev1.is_on_arc() {
                *fd.change_vertex(isfirst, fi_idx) = cp_prev1.clone();
                sdprev_w.change_vertex(isfirst, fi_idx).reset();
                sdprev_w
                    .change_vertex(isfirst, fi_idx)
                    .set_point(cp_last1.point());
                cp_last2.reset();
                cp_last2.set_point(cp_last1.point());
                *fd.change_vertex(!isfirst, fi_idx) = cp_last2;
            }
        }

        // in shrink case, self intersection is possible at <midIpoint>
        // (OCCT L3733-3758).  GAP: the generic Geom2dInt_GInter chain is
        // pending — geom2d_int_g_inter covers the conic pairs only, so the
        // tolerance widening below keeps the OCCT not-done outcome.
        let mut nb_sh = 0usize;
        while nb_sh < nbface {
            let probe = if is_on_same1 {
                shrink[nb_sh + 1] == 1
            } else {
                shrink[nb_sh] == 0
            };
            if probe {
                break;
            }
            nb_sh += 1;
        }
        let nb_sh = nb_sh.min(nn - 1);
        let _ = nb_sh;
        {
            let pds = this.my_ds.as_mut().expect("DS").change_point(mid_ipoint);
            let _ = (pds.point(), pds.tolerance());
        }
        let _ = c3d;
    }

    // write the mutated CommonPoints back (the OCCT references CV1/CV2 are
    // in-place).
    *fd.change_vertex(isfirst, 1) = cv1.clone();
    *fd.change_vertex(isfirst, 2) = cv2.clone();
}
