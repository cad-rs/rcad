//! OCCT BRepCheck_Wire::SelfIntersect (Wire.cxx L1074-1742).
//!
//! Split out of `brep_check_wire.rs` (the 2000-line file limit) and
//! registered like the other BRepCheck files in the module registrar
//! (`brep_check_analyzer.rs`).  Runs the general 2D curve/curve intersector
//! `Geom2dInt_GInter` — the rcad real body is
//! `crate::geomalgo::geom2d_int::GInter` (the `IntCurve_IntCurveCurveGen`
//! instantiation, `Geom2dInt_GInter_0.cxx` + `IntCurve_IntCurveCurveGen.gxx`).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, CurveEval, Line3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, ShapeType};

use crate::geomalgo::geom2d_int::{elclib2d, Curve2dAdaptor, GInter};
use crate::geomalgo::int_res2d::{Domain as Res2dDomain, Position};
use crate::topalgo::brep_class::bnd_lib_add2d_curve::add_2d_curve;

use super::brep_check_result::{
    brep_check_add, brep_tool_tolerance_vertex, iterator_subshapes, BRepCheckStatus,
};
use super::brep_check_wire::{
    explorer_of_edge_vertices, face_surface_adaptor, uv_points, BRepCheckWire, ShapeSet,
};

/// OCCT `Standard::IsEqual(v1, v2)` (Standard_Real.hxx L148-151) —
/// `|v1 - v2| < RealSmall()` (= DBL_MIN, Standard_Real.hxx L143-146).
/// Used by BRepCheck_Wire::SelfIntersect L1444.
fn standard_is_equal(the_value1: f64, the_value2: f64) -> bool {
    (the_value1 - the_value2).abs() < f64::MIN_POSITIVE
}

/// OCCT `BRep_Tool::Parameter(V, E)` — the vertex parameter on the edge
/// (BRep_Tool.cxx L1503-1512; the 3-argument overload it calls lives at
/// L1340-1497). The rcad form reads the `vertex_params` slot of the edge; the
/// OCCT `Standard_NoSuchObject` ("BRep_Tool:: no parameter on edge") maps to
/// the panic.
fn vertex_parameter_on_edge(brep: &BRep, v: &Shape, e: &Shape, f: &Shape) -> f64 {
    match brep.parameter_on_edge(v, e, f) {
        Some(p) => p,
        None => panic!("Standard_NoSuchObject: BRep_Tool::Parameter(V, E) not found"),
    }
}

/// OCCT gp_Dir2d::Angle(theOther) (gp_Dir2d.cxx L26-65) — the signed angle
/// between two unit 2D directions, in [-PI, PI].
fn gp_dir2d_angle(me: DVec2, the_other: DVec2) -> f64 {
    // OCCT L34-35.
    let cosinus = me.dot(the_other);
    let sinus = me.x * the_other.y - me.y * the_other.x;
    // OCCT L36-64.
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        if sinus > 0.0 {
            cosinus.acos()
        } else {
            -cosinus.acos()
        }
    } else if cosinus > 0.0 {
        sinus.asin()
    } else if sinus > 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        -std::f64::consts::PI - sinus.asin()
    }
}

/// OCCT gp_Dir2d::IsParallel(theOther, theAngularTolerance)
/// (gp_Dir2d.hxx L422-431).
fn gp_dir2d_is_parallel(me: DVec2, the_other: DVec2, the_angular_tolerance: f64) -> bool {
    let mut an_ang = gp_dir2d_angle(me, the_other);
    if an_ang < 0.0 {
        an_ang = -an_ang;
    }
    an_ang <= the_angular_tolerance || std::f64::consts::PI - an_ang <= the_angular_tolerance
}

impl BRepCheckWire {
    /// OCCT BRepCheck_Wire::SelfIntersect (Wire.cxx L1074-1742).
    ///
    /// `c1` / `c2` / `my_domain1` / `ip_param_on_first` mirror the OCCT locals
    /// declared at L1082-1091 and assigned inside the loop heads; Rust needs
    /// an initializer the first assignment always overwrites.
    #[allow(unused_assignments)]
    pub fn self_intersect(
        &mut self,
        brep: &BRep,
        f: &Shape,
        ret_e1: &mut Shape,
        ret_e2: &mut Shape,
        update: bool,
    ) -> BRepCheckStatus {
        // OCCT L1079-1080: aHList = myMap(myShape) / aStatusList.
        let my_shape = self.base.my_shape.clone();
        // OCCT L1092: tolint = 1.e-10.
        let tolint = 1e-10;
        // OCCT L1093-1094: HS = new BRepAdaptor_Surface();
        // HS->Initialize(F, false).
        let hs = face_surface_adaptor(brep, f);
        // OCCT L1096-1102: EMap <- the direct edge children of myShape.
        let mut e_map = ShapeSet::new();
        for it1 in iterator_subshapes(brep, &my_shape) {
            if it1.shape_type() == ShapeType::Edge {
                e_map.add(&it1);
            }
        }
        // OCCT L1104: Nbedges = EMap.Extent().
        let nbedges = e_map.extent();
        // OCCT L1105-1112: the empty-wire guard.
        if nbedges == 0 {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("SelfIntersect: myShape must be bound");
                brep_check_add(lst, BRepCheckStatus::EmptyWire);
            }
            return BRepCheckStatus::EmptyWire;
        }

        // OCCT L1114-1116: tabDom / tabCur / boxes.
        // `IntRes2d_Domain* tabDom = new IntRes2d_Domain[Nbedges]` — the
        // default constructed domain (no first/last point), never written for
        // index 0 (see the i == 1 pass below).
        let mut tab_dom: Vec<Res2dDomain> =
            (0..nbedges).map(|_| Res2dDomain::infinite()).collect();
        let mut tab_cur: Vec<Option<Curve2d>> = vec![None; nbedges];
        let mut boxes: Vec<BndBox2d> = (0..nbedges).map(|_| BndBox2d::new()).collect();

        // OCCT L1090-1091: C1 / C2 (Geom2dAdaptor_Curve), myDomain1, and the
        // single Geom2dInt_GInter Inter reused by every Perform.
        let mut c1: Option<Curve2d> = None;
        let mut c2: Option<Curve2d> = None;
        let mut my_domain1 = Res2dDomain::infinite();
        let mut inter = GInter::new();

        // OCCT L1118: for (i = 1; i <= Nbedges; i++)
        for i in 0..nbedges {
            // OCCT L1120: const TopoDS_Edge& E1 = EMap.FindKey(i).
            let e1 = e_map.items()[i].clone();
            if i == 0 {
                // OCCT L1123: pcu = BRep_Tool::CurveOnSurface(E1, F, first1, last1).
                let Some((pcu, first1_raw, last1_raw)) = brep.curve_on_surface(&e1, f) else {
                    // OCCT L1124-1133: the null-pcurve guard.
                    *ret_e1 = e1;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("SelfIntersect: myShape must be bound");
                        brep_check_add(lst, BRepCheckStatus::SelfIntersectingWire);
                    }
                    return BRepCheckStatus::SelfIntersectingWire;
                };
                // OCCT L1135: C1.Load(pcu).
                c1 = Some(pcu.clone());
                // OCCT L1136-1144: clamp onto the adaptor range when C1 is
                // not periodic (to avoid exception in Segment if C1 is
                // BSpline - IFV).
                let mut first1 = first1_raw;
                let mut last1 = last1_raw;
                if !Curve2dEval::is_periodic(&pcu) {
                    let dom = Curve2dEval::default_domain(&pcu);
                    if dom[0] > first1 {
                        first1 = dom[0];
                    }
                    if dom[1] < last1 {
                        last1 = dom[1];
                    }
                }
                // OCCT L1146-1147: BRep_Tool::UVPoints(E1, F, pfirst1, plast1);
                // myDomain1.SetValues(pfirst1, first1, tolint, plast1, last1, tolint).
                let (pfirst1, plast1) = uv_points(&brep, &e1, f, &pcu, first1_raw, last1_raw);
                my_domain1 = Res2dDomain::bounded(pfirst1, first1, tolint, plast1, last1, tolint);
                // OCCT L1149: BndLib_Add2dCurve::Add(C1, first1, last1,
                // Precision::PConfusion(), boxes(1)).
                add_2d_curve(
                    &pcu,
                    first1,
                    last1,
                    rcad_kernel::precision::PCONFUSION,
                    &mut boxes[i],
                );
            } else {
                // OCCT L1151-1156: C1.Load(tabCur(i)); myDomain1 = tabDom[i-1].
                c1 = tab_cur[i].clone();
                my_domain1 = tab_dom[i].clone();
            }
            let c1_ref = c1
                .as_ref()
                .expect("SelfIntersect: C1 must be loaded (BRep_Tool::CurveOnSurface)");

            // OCCT L1159: Inter.Perform(C1, myDomain1, tolint, tolint) — the
            // self-intersection of C1.
            inter.perform_cd(c1_ref, &my_domain1, tolint, tolint);

            // OCCT L1161-1225: the self-intersection point loop.
            if inter.is_done() {
                // OCCT L1163.
                let nbp = inter.nb_points();
                for p in 1..=nbp {
                    // OCCT L1167-1169.
                    let ip = inter.point(p);
                    let tr1 = ip.transition_of_first().clone();
                    let tr2 = ip.transition_of_second().clone();
                    if tr1.position_on_curve() == Position::Middle
                        || tr2.position_on_curve() == Position::Middle
                    {
                        // OCCT L1171-1188: check the point against the true
                        // (3d) tolerances — a point inside a vertex tolerance
                        // is a correct intersection.
                        let mut localok = false;
                        // OCCT L1175: ConS = BRep_Tool::Curve(E1, L, f, l).
                        let con_s = brep.edge_curve_world(&e1);
                        let p3d;
                        if let Some((a_curve, _range)) = con_s.as_ref() {
                            // OCCT L1178-1179: P3d = ConS->Value(...);
                            // P3d.Transform(L.Transformation()).
                            p3d = a_curve.point_at(ip.param_on_first());
                        } else {
                            // OCCT L1184-1185: gp_Pnt2d aP2d = C1.Value(...);
                            // P3d = HS->Value(aP2d.X(), aP2d.Y()).
                            let a_p2d = c1_ref.point_at(ip.param_on_first());
                            p3d = hs
                                .as_ref()
                                .expect("SelfIntersect: BRepAdaptor_Surface must be initialized")
                                .point_at(a_p2d.x, a_p2d.y);
                        }
                        // OCCT L1187-1196: TopExp_Explorer ExplVtx(E1, VERTEX).
                        for vtt in explorer_of_edge_vertices(&e1) {
                            let p3dvtt = brep.vertex_position(&vtt);
                            let tolvtt = brep_tool_tolerance_vertex(brep, &vtt);
                            let tolvtt = tolvtt * tolvtt;
                            let p3dvtt_distance_p3d = p3dvtt.distance_squared(p3d);
                            if p3dvtt_distance_p3d <= tolvtt {
                                localok = true;
                                break;
                            }
                        }
                        // OCCT L1199-1211.
                        if !localok {
                            *ret_e1 = e1;
                            if update {
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(&my_shape)
                                    .expect("SelfIntersect: myShape must be bound");
                                brep_check_add(lst, BRepCheckStatus::SelfIntersectingWire);
                            }
                            return BRepCheckStatus::SelfIntersectingWire;
                        }
                    }
                }
            }

            // OCCT L1227: for (j = i + 1; j <= Nbedges; j++)
            for j in (i + 1)..nbedges {
                // OCCT L1229.
                let e2 = e_map.items()[j].clone();
                if i == 0 {
                    // OCCT L1247: tabCur(j) = BRep_Tool::CurveOnSurface(E2, F, first2, last2).
                    let pc2 = brep.curve_on_surface(&e2, f);
                    tab_cur[j] = pc2.as_ref().map(|(c, _f, _l)| c.clone());
                    // OCCT L1248-1281.
                    let has_range = match pc2.as_ref() {
                        Some((_, f2, l2)) => *l2 > *f2,
                        None => false,
                    };
                    if has_range {
                        let (c2c, first2_raw, last2_raw) = match pc2 {
                            Some(v) => v,
                            None => unreachable!(),
                        };
                        // OCCT L1251: C2.Load(tabCur(j)).
                        c2 = Some(c2c.clone());
                        // OCCT L1252-1262: the periodic guard.
                        let mut first2 = first2_raw;
                        let mut last2 = last2_raw;
                        if !Curve2dEval::is_periodic(&c2c) {
                            let dom = Curve2dEval::default_domain(&c2c);
                            if dom[0] > first2 {
                                first2 = dom[0];
                            }
                            if dom[1] < last2 {
                                last2 = dom[1];
                            }
                        }
                        // OCCT L1264-1266: BRep_Tool::UVPoints(E2, F, pfirst2,
                        // plast2); tabDom[j-1].SetValues(...).
                        let (pfirst2, plast2) =
                            uv_points(&brep, &e2, f, &c2c, first2_raw, last2_raw);
                        tab_dom[j] =
                            Res2dDomain::bounded(pfirst2, first2, tolint, plast2, last2, tolint);
                        // OCCT L1268.
                        add_2d_curve(
                            &c2c,
                            first2,
                            last2,
                            rcad_kernel::precision::PCONFUSION,
                            &mut boxes[j],
                        );
                    } else if pc2.is_none() {
                        // OCCT L1271-1279: tabCur(j) is null.
                        return BRepCheckStatus::NoCurveOnSurface;
                    } else {
                        // OCCT L1271-1280: last2 <= first2.
                        return BRepCheckStatus::InvalidRange;
                    }
                } else {
                    // OCCT L1281-1286: C2.Load(tabCur(j)).
                    c2 = tab_cur[j].clone();
                }
                let c2_ref = c2
                    .as_ref()
                    .expect("SelfIntersect: C2 must be loaded (BRep_Tool::CurveOnSurface)");

                // OCCT L1288-1291: the box rejection.
                if boxes[i].is_out_box(&boxes[j]) {
                    continue;
                }
                // OCCT L1292-1296: the same-edge rejection.
                if e1.is_same(&e2) {
                    continue;
                }

                // OCCT L1302: Inter.Perform(C1, myDomain1, C2, tabDom[j-1],
                // tolint, tolint).
                inter.perform_cd_cd(c1_ref, &my_domain1, c2_ref, &tab_dom[j], tolint, tolint);

                // OCCT L1304-1732.
                if inter.is_done() {
                    // OCCT L1306-1310.
                    let mut common_vertices: Vec<Shape> = Vec::new();
                    let mut v_map = ShapeSet::new();
                    // OCCT L1312-1315: it(E1) — Vmap.
                    for v in explorer_of_edge_vertices(&e1) {
                        v_map.add(&v);
                    }
                    // OCCT L1317-1323: it(E2) — the common vertices.
                    for v in explorer_of_edge_vertices(&e2) {
                        if v_map.contains(&v) {
                            common_vertices.push(v);
                        }
                    }
                    // OCCT L1326-1328.
                    let nbp = inter.nb_points();
                    let nbs = inter.nb_segments();
                    let mut ip_param_on_first = 0.;
                    let mut ip_param_on_second = 0.;

                    // OCCT L1332-1521: the intersection points.
                    for p in 1..=nbp {
                        // OCCT L1334-1338.
                        let ip = inter.point(p);
                        ip_param_on_first = ip.param_on_first();
                        ip_param_on_second = ip.param_on_second();
                        let tr1 = ip.transition_of_first().clone();
                        let tr2 = ip.transition_of_second().clone();
                        if tr1.position_on_curve() != Position::Middle
                            && tr2.position_on_curve() != Position::Middle
                        {
                            continue;
                        }
                        // OCCT L1339-1343.
                        let mut localok = false;
                        // OCCT L1345-1350.
                        let con_s = brep.edge_curve_world(&e1);
                        let con_s2 = brep.edge_curve_world(&e2);
                        let (f1, l1) = con_s
                            .as_ref()
                            .map(|(_, r)| (r[0], r[1]))
                            .unwrap_or((0.0, 0.0));
                        let (f2, l2) = con_s2
                            .as_ref()
                            .map(|(_, r)| (r[0], r[1]))
                            .unwrap_or((0.0, 0.0));
                        // OCCT L1351-1358 (gka: protect against working out of
                        // the edge range).
                        if f1 - ip_param_on_first > rcad_kernel::precision::p_confusion()
                            || ip_param_on_first - l1 > rcad_kernel::precision::p_confusion()
                            || f2 - ip_param_on_second > rcad_kernel::precision::p_confusion()
                            || ip_param_on_second - l2 > rcad_kernel::precision::p_confusion()
                        {
                            continue;
                        }
                        // OCCT L1359.
                        let mut tolvtt = 0.;
                        // OCCT L1361-1386: P3d / P3d2.
                        let p3d;
                        if let Some((a_curve, _range)) = con_s.as_ref() {
                            p3d = a_curve.point_at(ip_param_on_first);
                        } else {
                            let a_p2d = c1_ref.point_at(ip_param_on_first);
                            p3d = hs
                                .as_ref()
                                .expect("SelfIntersect: BRepAdaptor_Surface must be initialized")
                                .point_at(a_p2d.x, a_p2d.y);
                        }
                        let p3d2;
                        if let Some((a_curve, _range)) = con_s2.as_ref() {
                            p3d2 = a_curve.point_at(ip_param_on_second);
                        } else {
                            let a_p2d = c2_ref.point_at(ip_param_on_second);
                            p3d2 = hs
                                .as_ref()
                                .expect("SelfIntersect: BRepAdaptor_Surface must be initialized")
                                .point_at(a_p2d.x, a_p2d.y);
                        }
                        // OCCT L1387-1402: the common-vertex tolerance test.
                        for vtt in common_vertices.iter() {
                            let p3dvtt = brep.vertex_position(vtt);
                            tolvtt = brep_tool_tolerance_vertex(brep, vtt);
                            tolvtt = 1.1 * tolvtt;
                            tolvtt = tolvtt * tolvtt;
                            let p3dvtt_distance_p3d = p3dvtt.distance_squared(p3d);
                            let p3dvtt_distance_p3d2 = p3dvtt.distance_squared(p3d2);
                            if p3dvtt_distance_p3d <= tolvtt && p3dvtt_distance_p3d2 <= tolvtt {
                                localok = true;
                                break;
                            }
                        }

                        // OCCT L1404-1412: check the maximum yawn between the
                        // two edges.
                        if !localok && !common_vertices.is_empty() {
                            // OCCT L1422-1424.
                            let mut distauvtxleplusproche;
                            let mut v_para_on_edge1;
                            let mut v_para_on_edge2;
                            let mut vertex_le_plus_proche;
                            let _ = tolvtt;
                            v_para_on_edge1 = 0.;
                            v_para_on_edge2 = 0.;
                            distauvtxleplusproche = rcad_kernel::REAL_LAST;
                            vertex_le_plus_proche = DVec3::ZERO;
                            // OCCT L1425-1439: find the nearest common vertex.
                            for vtt in common_vertices.iter() {
                                let p3dvtt = brep.vertex_position(vtt);
                                let disptvtx = p3d.distance(p3dvtt);
                                if disptvtx < distauvtxleplusproche {
                                    vertex_le_plus_proche = p3dvtt;
                                    distauvtxleplusproche = disptvtx;
                                    v_para_on_edge1 = vertex_parameter_on_edge(brep, vtt, &e1, f);
                                    v_para_on_edge2 = vertex_parameter_on_edge(brep, vtt, &e2, f);
                                } else if standard_is_equal(distauvtxleplusproche, disptvtx) {
                                    // OCCT L1444-1456: eap — case of a closed edge.
                                    let new_v_para_on_edge1 =
                                        vertex_parameter_on_edge(brep, vtt, &e1, f);
                                    let new_v_para_on_edge2 =
                                        vertex_parameter_on_edge(brep, vtt, &e2, f);
                                    if (ip_param_on_first - v_para_on_edge1).abs()
                                        + (ip_param_on_second - v_para_on_edge2).abs()
                                        > (ip_param_on_first - new_v_para_on_edge1).abs()
                                            + (ip_param_on_second - new_v_para_on_edge2).abs()
                                    {
                                        vertex_le_plus_proche = p3dvtt;
                                        v_para_on_edge1 = new_v_para_on_edge1;
                                        v_para_on_edge2 = new_v_para_on_edge2;
                                    }
                                }
                            }
                            // OCCT L1459-1467: patch for the extra-ordinary
                            // situation (e.g. tolerance(v) == 0.).
                            if vertex_le_plus_proche.distance(p3d) <= GP_RESOLUTION
                                || vertex_le_plus_proche.distance(p3d2) <= GP_RESOLUTION
                            {
                                localok = true;
                            } else {
                                // OCCT L1468-1470.
                                let mut lig = Line3::new(
                                    vertex_le_plus_proche,
                                    p3d - vertex_le_plus_proche,
                                );
                                // OCCT L1471-1474.
                                let du1 = 0.1 * (ip_param_on_first - v_para_on_edge1);
                                let du2 = 0.1 * (ip_param_on_second - v_para_on_edge2);
                                let mut maxd1 = 0.;
                                let mut maxd2 = 0.;
                                // OCCT L1476-1497: edge 1.
                                localok = true;
                                let tole1 = brep.tolerance(&e1);
                                for k in 2..9 {
                                    if !localok {
                                        break;
                                    }
                                    let u = v_para_on_edge1 + k as f64 * du1;
                                    let p1;
                                    if let Some((a_curve, _range)) = con_s.as_ref() {
                                        p1 = a_curve.point_at(u);
                                    } else {
                                        let a_p2d = c1_ref.point_at(u);
                                        p1 = hs
                                            .as_ref()
                                            .expect(
                                                "SelfIntersect: BRepAdaptor_Surface must be initialized",
                                            )
                                            .point_at(a_p2d.x, a_p2d.y);
                                    }
                                    let d1 = lig.distance(p1);
                                    if d1 > maxd1 {
                                        maxd1 = d1;
                                    }
                                    if d1 > tole1 * 2.0 {
                                        localok = false;
                                    }
                                }
                                // OCCT L1498-1501: Lig.SetDirection(aTmpDir)
                                // with aTmpDir = P3d2 - VertexLePlusProche.
                                lig = Line3::new(
                                    vertex_le_plus_proche,
                                    p3d2 - vertex_le_plus_proche,
                                );
                                // OCCT L1501-1519: edge 2.
                                let tole2 = brep.tolerance(&e2);
                                for k in 2..9 {
                                    if !localok {
                                        break;
                                    }
                                    let u = v_para_on_edge2 + k as f64 * du2;
                                    let p2;
                                    if let Some((a_curve, _range)) = con_s2.as_ref() {
                                        p2 = a_curve.point_at(u);
                                    } else {
                                        let a_p2d = c2_ref.point_at(u);
                                        p2 = hs
                                            .as_ref()
                                            .expect(
                                                "SelfIntersect: BRepAdaptor_Surface must be initialized",
                                            )
                                            .point_at(a_p2d.x, a_p2d.y);
                                    }
                                    let d2 = lig.distance(p2);
                                    if d2 > maxd2 {
                                        maxd2 = d2;
                                    }
                                    if d2 > tole2 * 2.0 {
                                        localok = false;
                                    }
                                }
                            }
                        }

                        // OCCT L1565-1580.
                        if !localok {
                            *ret_e1 = e1.clone();
                            *ret_e2 = e2.clone();
                            if update {
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(&my_shape)
                                    .expect("SelfIntersect: myShape must be bound");
                                brep_check_add(lst, BRepCheckStatus::SelfIntersectingWire);
                            }
                            return BRepCheckStatus::SelfIntersectingWire;
                        }
                    }

                    // OCCT L1584-1729: the intersection segments.
                    for s in 1..=nbs {
                        // OCCT L1586-1588.
                        let seg = inter.segment(s);
                        if !(seg.has_first_point() && seg.has_last_point()) {
                            continue;
                        }
                        // OCCT L1590-1597.
                        let mut localok = false;
                        let pseg = [seg.first_point().clone(), seg.last_point().clone()];
                        // OCCT L1600-1602: at least one of the extremities of
                        // the segment must be inside the tolerance of a common
                        // vertex.
                        for k in 0..2 {
                            // OCCT L1604-1610.
                            ip_param_on_first = pseg[k].param_on_first();
                            ip_param_on_second = pseg[k].param_on_second();
                            let tr1 = pseg[k].transition_of_first().clone();
                            let tr2 = pseg[k].transition_of_second().clone();
                            let a_pcr1 = tr1.position_on_curve();
                            let a_pcr2 = tr2.position_on_curve();
                            if a_pcr1 != Position::Middle && a_pcr2 != Position::Middle {
                                // OCCT L1613-1617.
                                let a_ct1 = c1_ref.get_type();
                                let a_ct2 = c2_ref.get_type();
                                if a_ct1 == crate::geomalgo::geom2d_int::Curve2dType::Line
                                    && a_ct2 == crate::geomalgo::geom2d_int::Curve2dType::Line
                                {
                                    // OCCT L1618-1662: check for the two
                                    // lines coincidence.
                                    let a_par_t = 0.43213918;
                                    //
                                    let a_tol_e1 = brep.tolerance(&e1);
                                    let a_tol_e2 = brep.tolerance(&e2);
                                    let a_tol2 = a_tol_e1 + a_tol_e2;
                                    let a_tol2 = a_tol2 * a_tol2;
                                    //
                                    let a_l1 = c1_ref.line();
                                    let a_l2 = c2_ref.line();
                                    //
                                    let a_t11 = pseg[0].param_on_first();
                                    let a_t12 = pseg[1].param_on_first();
                                    let a_t21 = pseg[0].param_on_second();
                                    let a_t22 = pseg[1].param_on_second();
                                    //
                                    let a_t1m = (1.0 - a_par_t) * a_t11 + a_par_t * a_t12;
                                    let a_p1m = c1_ref.point_at(a_t1m);
                                    //
                                    // gp_Lin2d::SquareDistance(gp_Pnt2d)
                                    // (gp_Lin2d.hxx L240-247):
                                    // ((P - Location) ^ Direction)^2.
                                    let a_coord = a_p1m - a_l2.origin;
                                    let a_d = a_coord.x * a_l2.direction.y
                                        - a_coord.y * a_l2.direction.x;
                                    let a_d2 = a_d * a_d;
                                    if a_d2 < a_tol2 {
                                        // OCCT L1637: ElCLib::Parameter(aL2, aP1m).
                                        let a_t2m = elclib2d::line_parameter(
                                            a_l2.origin,
                                            a_l2.direction,
                                            a_p1m,
                                        );
                                        if a_t2m > a_t21 && a_t2m < a_t22 {
                                            // OCCT L1640-1643.
                                            if gp_dir2d_is_parallel(
                                                a_l1.direction,
                                                a_l2.direction,
                                                rcad_kernel::ANGULAR,
                                            ) {
                                                localok = false;
                                                break;
                                            }
                                        }
                                    }
                                }
                                // OCCT L1665-1667.
                                localok = true;
                                break;
                            }
                            // OCCT L1669-1700.
                            let con_s = brep.edge_curve_world(&e1);
                            let con_s2 = brep.edge_curve_world(&e2);
                            let p3d;
                            if let Some((a_curve, _range)) = con_s.as_ref() {
                                p3d = a_curve.point_at(ip_param_on_first);
                            } else {
                                let a_p2d = c1_ref.point_at(ip_param_on_first);
                                p3d = hs
                                    .as_ref()
                                    .expect(
                                        "SelfIntersect: BRepAdaptor_Surface must be initialized",
                                    )
                                    .point_at(a_p2d.x, a_p2d.y);
                            }
                            let p3d2;
                            if let Some((a_curve, _range)) = con_s2.as_ref() {
                                p3d2 = a_curve.point_at(ip_param_on_second);
                            } else {
                                let a_p2d = c2_ref.point_at(ip_param_on_second);
                                p3d2 = hs
                                    .as_ref()
                                    .expect(
                                        "SelfIntersect: BRepAdaptor_Surface must be initialized",
                                    )
                                    .point_at(a_p2d.x, a_p2d.y);
                            }
                            // OCCT L1701-1717.
                            for vtt in common_vertices.iter() {
                                let p3dvtt = brep.vertex_position(vtt);
                                let mut tolvtt = brep_tool_tolerance_vertex(brep, vtt);
                                tolvtt = 1.1 * tolvtt;
                                tolvtt = tolvtt * tolvtt;
                                let p3dvtt_distance_p3d = p3dvtt.distance_squared(p3d);
                                let p3dvtt_distance_p3d2 = p3dvtt.distance_squared(p3d2);
                                if p3dvtt_distance_p3d <= tolvtt
                                    && p3dvtt_distance_p3d2 <= tolvtt
                                {
                                    localok = true;
                                    break;
                                }
                            }
                            if localok {
                                break;
                            }
                        }
                        // OCCT L1719-1730.
                        if !localok {
                            *ret_e1 = e1.clone();
                            *ret_e2 = e2.clone();
                            if update {
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(&my_shape)
                                    .expect("SelfIntersect: myShape must be bound");
                                brep_check_add(lst, BRepCheckStatus::SelfIntersectingWire);
                            }
                            return BRepCheckStatus::SelfIntersectingWire;
                        }
                    }
                }
            }
        }

        // OCCT L1736-1741.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("SelfIntersect: myShape must be bound");
            brep_check_add(lst, BRepCheckStatus::NoError);
        }
        BRepCheckStatus::NoError
    }
}
