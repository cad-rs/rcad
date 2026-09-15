//! OCCT BRepLib::SameParameter (ModelingAlgorithms/TKTopAlgo/BRepLib/
//! BRepLib.cxx L1032-1740) — the 1:1 translation:
//! - static `EvalTol` (cxx L1034-1066),
//! - static `ComputeTol` (cxx L1070-1188),
//! - static `GetCurve3d` (cxx L1192-1218),
//! - `UpdateVTol` (cxx L1222-1233),
//! - `BRepLib::SameParameter(E, Tol)` — the void overload (cxx L1237-1247),
//! - `BRepLib::SameParameter(E, Tol, NewTol, IsUseOldEdge)` — the 4-arg
//!   engine (cxx L1251-1740).
//!
//! Architecture difference (the rcad pool model, topods::BRep): the OCCT
//! BRep_Builder edits the shared edge TShape in place through the handle;
//! rcad mirrors this through the BRep pool (`BRep::edge_mut_inplace`) —
//! the engine takes `&mut BRep` where OCCT mutates, `&BRep` where OCCT only
//! reads.  These are free functions (the OCCT statics of BRepLib.cxx); the
//! pool-free callers adopt the edge graph into a standalone pool with the
//! Arc SHARED (topo_builder::brep_from_shape), so the in-place writes stay
//! observable through the caller's own Shape handle.
//!
//! Encodings recorded along the way:
//! - `occ::handle<GeomAdaptor_Curve> HC` -> [`GeomCurveAdaptor`] (the
//!   `(Curve3, first, last)` window carrier),
//! - `occ::handle<Geom2dAdaptor_Curve> HC2d` -> [`Geom2dCurveAdaptor`],
//! - `occ::handle<GeomAdaptor_Surface> HS` -> [`GeomSurfaceAdaptor`],
//! - `handle(Geom2d_BSplineCurve) bs2d` ->
//!   [`Geom2dBSplineCurve`] (the OCCT-faithful knot/mult carrier); the
//!   adaptor reads go through the legacy kernel `Curve2d` encoding
//!   (`geom2d_bspline_to_curve2d`), the pcurve write-backs carry the same
//!   value,
//! - the NIZHNY-OCC486 `m_TrimmedPeriodical` quirk (cxx L1304-1311) reads
//!   `Curve3::Trimmed`.

use std::sync::Arc;

use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_locate_ext_pc::LocateExtPC;
use rcad_kernel::base::geom2d_convert::Geom2dBSplineCurve;
use rcad_kernel::base::proj_lib::adaptor::{Adaptor2dCurve2d, Adaptor3dCurve, Adaptor3dSurface};
use rcad_kernel::base::proj_lib::{
    CurveType, Geom2dCurveAdaptor, GeomCurveAdaptor, GeomSurfaceAdaptor,
};
use rcad_kernel::core::precision::{is_infinite_value, CONFUSION, PCONFUSION};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::{bspl_lib, GeomAbsShape};
use rcad_kernel::topo::topods::{
    edge_data_pool_free, shape_is_in_pool, BRep, CurveRepresentation, TEdgeData, TShape,
    tshape_flags,
};
use rcad_kernel::topo_shape::Shape;

use rcad_kernel::base::geom2d_convert::bspline_curve::SmoothShape as SmoothShapeQ;

use crate::geomalgo::approx_curvilinear_parameter::ApproxCurvilinearParameter;
use crate::geomalgo::approx_same_parameter::ApproxSameParameter;
use crate::geomalgo::geom_lib_same_range::{c0_bspline_to_c1_bspline_curve, same_range};
use crate::topalgo::brep_lib_validate_edge::{
    brep_check_prec_curve, brep_check_prec_surface, GeomAdaptorCurve, GeomAdaptorSurface,
};

/// OCCT BRepLib.cxx L1294: `const int NCONTROL = 22`.
const NCONTROL: i32 = 22;
/// OCCT BRepLib.cxx L1349: `const double BigError = 1.e10`.
const BIG_ERROR: f64 = 1.0e10;

// =========================================================================
// static EvalTol (BRepLib.cxx L1034-1066).
// =========================================================================

/// OCCT static EvalTol(pc, s, gac, tol, tolbail) — samples five pcurve
/// points, projects them on the 3d curve and reports through `the_tolbail`
/// the maximal deviation; returns `ok > 2`.
fn eval_tol(
    the_pc: &Curve2d,
    the_s: &Surface3,
    the_gac: &GeomCurveAdaptor,
    the_tol: f64,
    the_tolbail: &mut f64,
) -> bool {
    // OCCT L1040: int ok = 0.
    let mut a_ok = 0;
    // OCCT L1041-1042.
    let a_f = the_gac.first;
    let a_l = the_gac.last;
    // OCCT L1043-1044: Extrema_LocateExtPC Projector;
    //                   Projector.Initialize(gac, f, l, tol).
    let a_cth = CurveToolHandle::for_curve3(&the_gac.curve, the_gac, the_gac);
    let mut a_projector = LocateExtPC::new();
    a_projector.initialize(&a_cth, a_f, a_l, the_tol);
    // OCCT L1047: tolbail = tol.
    *the_tolbail = the_tol;
    for a_i in 1..=5 {
        let mut a_t = a_i as f64 / 6.0;
        a_t = (1.0 - a_t) * a_f + a_t * a_l;
        // OCCT L1052: pc->Value(t).Coord(u, v).
        let a_uv = Curve2dEval::point_at(the_pc, a_t);
        // OCCT L1053: p = s->Value(u, v).
        let a_p = SurfaceEval::point_at(the_s, a_uv.x, a_uv.y);
        // OCCT L1054-1055.
        a_projector.perform(a_p, a_t);
        if a_projector.is_done() {
            let a_dist2 = a_projector.square_distance();
            if a_dist2 > *the_tolbail * *the_tolbail {
                *the_tolbail = a_dist2.sqrt();
            }
            a_ok += 1;
        }
    }
    // OCCT L1065.
    a_ok > 2
}

// =========================================================================
// static ComputeTol (BRepLib.cxx L1070-1188).
// =========================================================================

/// OCCT static ComputeTol(c3d, c2d, surf, nbp) — the maximal deviation of
/// the pcurve from the 3d curve over NCONTROL samples, with the
/// out-of-surface `dapp` estimate and the small-distance analytic branch.
fn compute_tol(
    the_c3d: &dyn Adaptor3dCurve,
    the_c2d: &dyn Adaptor2dCurve2d,
    the_surf: &dyn Adaptor3dSurface,
    the_nbp: i32,
) -> f64 {
    // OCCT L1077-1078: NCollection_Array1<double> dist(1, nbp + 10);
    //                  dist.Init(-1.).
    let mut a_dist = vec![-1.0f64; the_nbp as usize + 10];

    // OCCT L1081-1085.
    let a_uf = the_surf.first_u_parameter();
    let a_ul = the_surf.last_u_parameter();
    let a_vf = the_surf.first_v_parameter();
    let a_vl = the_surf.last_v_parameter();
    let a_du = 0.01 * (a_ul - a_uf);
    let a_dv = 0.01 * (a_vl - a_vf);
    let a_is_u_periodic = the_surf.is_u_periodic();
    let a_is_v_periodic = the_surf.is_v_periodic();
    let a_dsdu = 1.0 / the_surf.u_resolution(1.0);
    let a_dsdv = 1.0 / the_surf.v_resolution(1.0);
    let mut a_d2 = 0.0f64;
    let a_first = the_c3d.first_parameter();
    let a_last = the_c3d.last_parameter();
    let mut a_dapp = -1.0f64;
    for a_i in 0..=the_nbp {
        let the_t = a_i as f64 / the_nbp as f64;
        let the_u = a_first * (1.0 - the_t) + a_last * the_t;
        let a_pc3d = the_c3d.value(the_u);
        let a_puv = the_c2d.value(the_u);
        if !a_is_u_periodic {
            if a_puv.x < a_uf - a_du {
                a_dapp = a_dapp.max(a_dsdu * (a_uf - a_puv.x));
                continue;
            } else if a_puv.x > a_ul + a_du {
                a_dapp = a_dapp.max(a_dsdu * (a_puv.x - a_ul));
                continue;
            }
        }
        if !a_is_v_periodic {
            if a_puv.y < a_vf - a_dv {
                a_dapp = a_dapp.max(a_dsdv * (a_vf - a_puv.y));
                continue;
            } else if a_puv.y > a_vl + a_dv {
                a_dapp = a_dapp.max(a_dsdv * (a_puv.y - a_vl));
                continue;
            }
        }
        let a_pcons = the_surf.value(a_puv.x, a_puv.y);
        if is_infinite_value(a_pcons.x)
            || is_infinite_value(a_pcons.y)
            || is_infinite_value(a_pcons.z)
        {
            // OCCT L1126: d2 = Precision::Infinite(); break.
            a_d2 = f64::INFINITY;
            break;
        }
        let a_temp = a_pc3d.distance_squared(a_pcons);

        // OCCT L1131: dist(i + 1) = temp (the 1-based array slot).
        a_dist[a_i as usize] = a_temp;

        a_d2 = a_d2.max(a_temp);
    }

    // OCCT L1136-1139.
    if is_infinite_value(a_d2) {
        return a_d2;
    }

    // OCCT L1141-1145.
    a_d2 = a_d2.sqrt();
    if a_dapp > a_d2 {
        return a_dapp;
    }

    // OCCT L1147-1151.
    let mut a_ana = false;
    let mut a_big_d2 = 0.0f64;
    let mut a_n1 = 0;
    let mut a_n2 = 0;
    let mut a_n3 = 0;

    for a_i in 1..=the_nbp + 10 {
        if a_dist[(a_i - 1) as usize] > 0.0 {
            if a_dist[(a_i - 1) as usize] < 1.0 {
                a_n1 += 1;
            } else {
                a_n2 += 1;
            }
        }
    }

    // OCCT L1168-1182.
    if a_n1 > a_n2 && a_n2 != 0 {
        a_n3 = 100 * a_n2 / (a_n1 + a_n2);
    }
    if a_n3 < 10 && a_n3 != 0 {
        a_ana = true;
        for a_i in 1..=the_nbp + 10 {
            let a_d = a_dist[(a_i - 1) as usize];
            if a_d > 0.0 && a_d < 1.0 {
                a_big_d2 = a_big_d2.max(a_d);
            }
        }
    }

    // OCCT L1184-1187.
    a_d2 = if !a_ana {
        1.5 * a_d2
    } else {
        1.5 * a_big_d2.sqrt()
    };
    a_d2 = a_d2.max(1.0e-7);
    a_d2
}

// =========================================================================
// static GetCurve3d (BRepLib.cxx L1192-1218).
// =========================================================================

/// OCCT static GetCurve3d(theEdge, theC3d, theF3d, theL3d, theLoc3d,
/// theCList) — walks the edge curve representations for the first
/// `IsCurve3D()` GCurve and reads Curve3D/First/Last/Location.
///
/// Architecture difference: the rcad Curve3D representation is the edge
/// curve slot (`TEdgeData::curve` + `TEdgeData::range`); the location
/// travels world-space (the edge location is applied by the kernel read,
/// identity in the offset pipelines), so `theLoc3d` is not returned and the
/// `theCList` out-parameter is the edge representations list the engine
/// iterates through `edge_data`.
fn get_curve3d(the_brep: &BRep, the_edge: &Shape) -> (Option<Curve3>, f64, f64) {
    if shape_is_in_pool(the_brep, the_edge) {
        let a_te = the_brep.edge(the_edge.clone());
        (a_te.curve.clone(), a_te.range[0], a_te.range[1])
    } else {
        let a_te = edge_data_pool_free(the_edge).expect("BRep_Tool: edge without a TEdgeData");
        (a_te.curve.clone(), a_te.range[0], a_te.range[1])
    }
}

// =========================================================================
// UpdateVTol (BRepLib.cxx L1222-1233).
// =========================================================================

/// OCCT UpdateVTol(theV1, theV2, theTol) — raises both vertex tolerances
/// through the builder (BRep_TVertex::UpdateTolerance keeps the max).
fn update_vtol(the_v1: Option<Shape>, the_v2: Option<Shape>, the_tol: f64) {
    // OCCT L1224: BRep_Builder aB.
    if let Some(mut a_v1) = the_v1 {
        if !a_v1.is_null() {
            crate::brep_algo::tool::builder_update_vertex_tol(&mut a_v1, the_tol);
        }
    }
    if let Some(mut a_v2) = the_v2 {
        if !a_v2.is_null() {
            crate::brep_algo::tool::builder_update_vertex_tol(&mut a_v2, the_tol);
        }
    }
}

// =========================================================================
// SameParameter — the void overload (BRepLib.cxx L1237-1247).
// =========================================================================

/// OCCT BRepLib::SameParameter(theEdge, theTolerance) — runs the 4-arg
/// engine with IsUseOldEdge = true and raises the edge vertices to the
/// resulting tolerance.
pub fn same_parameter(the_brep: &mut BRep, the_edge: &Shape, the_tolerance: f64) {
    // OCCT L1239-1240.
    let mut a_new_tol = -1.0f64;
    same_parameter_with_result(the_brep, the_edge, the_tolerance, &mut a_new_tol, true);
    // OCCT L1241-1246.
    if a_new_tol > 0.0 {
        let (a_v1, a_v2) = crate::brep_algo::tool::top_exp_vertices_raw(the_edge);
        update_vtol(a_v1, a_v2, a_new_tol);
    }
}

// =========================================================================
// SameParameter — the 4-arg engine (BRepLib.cxx L1251-1740).
// =========================================================================

/// OCCT BRepLib::SameParameter(theEdge, theTolerance, theNewTol,
/// IsUseOldEdge) — enforces the same parametrization of the edge pcurves
/// with the 3d curve; returns the processed edge (the copy when
/// IsUseOldEdge is false, null when nothing was done).
pub fn same_parameter_with_result(
    the_brep: &mut BRep,
    the_edge: &Shape,
    the_tolerance: f64,
    the_new_tol: &mut f64,
    is_use_old_edge: bool,
) -> Shape {
    // OCCT L1256-1259: if (BRep_Tool::SameParameter(theEdge)) return null.
    if edge_data(the_brep, the_edge).same_parameter {
        return Shape::null();
    }
    // OCCT L1260-1264: f3d = l3d = 0.; L3d; C3d; CList; GetCurve3d(...).
    let (mut a_c3d, mut f3d, mut l3d) = get_curve3d(the_brep, the_edge);
    if a_c3d.is_none() {
        return Shape::null();
    }

    // OCCT L1270-1272: BRep_Builder B; TopoDS_Edge aNE; handle(BRep_TEdge) aNTE.
    let a_ne: Shape;
    if is_use_old_edge {
        // OCCT L1273-1277: aNE = theEdge; aNTE = theEdge.TShape().
        a_ne = the_edge.clone();
    } else {
        // OCCT L1281: aNE = TopoDS::Edge(theEdge.EmptyCopied()) — the copy is
        // modified a little bit later, so copy anyway.
        a_ne = the_brep.empty_copied(the_edge);
        // OCCT L1282: GetCurve3d(aNE, ...) — C3d pointer and CList differ
        // after copying.
        let (a_c3d_copy, a_f3d, a_l3d) = get_curve3d(the_brep, &a_ne);
        a_c3d = a_c3d_copy;
        f3d = a_f3d;
        l3d = a_l3d;
        // OCCT L1285-1289: TopoDS_Iterator sit(theEdge); add the vertices
        // from the old edge to the new ones (B.Add(aNE, sit.Value())).
        let (a_v_old_first, a_v_old_last) = {
            let a_ed_old = edge_data(the_brep, the_edge);
            (a_ed_old.first.clone(), a_ed_old.last.clone())
        };
        for a_v in [&a_v_old_first, &a_v_old_last] {
            if a_v.is_null() {
                continue;
            }
            // BRep_Builder::Add(E, V) — the extremity assignment by
            // orientation (the tool::builder_add_edge_vertex encoding).
            let a_ed = the_brep.edge_mut_inplace(a_ne.clone());
            match a_v.orientation {
                rcad_kernel::topo::topods::Orientation::Reversed => a_ed.last = a_v.clone(),
                _ => a_ed.first = a_v.clone(),
            }
        }
    }

    // OCCT L1292: NCollection_List<...>::Iterator It(CList) — the rcad walk
    // is the index loop over the aNE representations below.
    let a_nb_repr = edge_data(the_brep, &a_ne).representations.len();

    // OCCT L1296-1301: HC / HC2d / HS handles + GAC / GAC2d / GAS refs — the
    // rcad adaptor values start at the OCCT default-constructed state
    // (GeomAdaptor_Curve::Reset) and are rebuilt on each Load below.
    let mut a_hc = {
        let mut a_gac = GeomCurveAdaptor::new(Curve3::Line(rcad_kernel::geom::Line3::new(
            glam::DVec3::ZERO,
            glam::DVec3::X,
        )));
        a_gac.reset();
        a_gac
    };
    let mut a_hc2d = Geom2dCurveAdaptor::new(Curve2d::Line(rcad_kernel::geom::Line2d::new(
        glam::DVec2::ZERO,
        glam::DVec2::X,
    )));
    let mut a_hs = GeomSurfaceAdaptor::new(Surface3::Plane(rcad_kernel::geom::Plane::new(
        glam::DVec3::ZERO,
        glam::DVec3::Z,
    )));

    // OCCT L1304-1310 (NIZHNY-OCC486): m_TrimmedPeriodical.
    let mut m_trimmed_periodical = false;
    let a_c3d_val = a_c3d.clone().unwrap();
    if let Curve3::Trimmed(a_tc) = &a_c3d_val {
        let a_gt_c = a_tc.curve.as_ref().clone();
        m_trimmed_periodical = CurveEval::is_periodic(&a_gt_c);
    }

    // OCCT L1313-1332.
    if !CurveEval::is_periodic(&a_c3d_val) {
        let [a_udeb, a_ufin] = CurveEval::default_domain(&a_c3d_val);
        if !m_trimmed_periodical {
            if a_udeb > f3d {
                f3d = a_udeb;
            }
            if l3d > a_ufin {
                l3d = a_ufin;
            }
        }
    }
    // OCCT L1333-1336: if (!L3d.IsIdentity()) C3d = Transformed(...) — the
    // rcad GetCurve3d reads world-space (the location is already applied),
    // so the branch is the recorded no-op.
    // OCCT L1337: GAC.Load(C3d, f3d, l3d).
    a_hc.load_with_range(a_c3d_val.clone(), f3d, l3d);

    // OCCT L1339: double Prec_C3d = BRepCheck::PrecCurve(GAC).
    let a_prec_c3d = brep_check_prec_curve(&GeomAdaptorCurve::new(a_c3d_val.clone(), f3d, l3d));

    // OCCT L1341-1349.
    let mut is_same_p = true;
    let mut maxdist = 0.0f64;
    let an_edge_tol = edge_data(the_brep, &a_ne).tolerance;
    let same_range_flag = edge_data(the_brep, &a_ne).same_range;
    let mut ya_pcu = false;

    // OCCT L1350-1352: It.Initialize(CList); while (It.More()).
    for a_it in 0..a_nb_repr {
        // OCCT L1355-1360: isANA / isBSP / GCurve / PC[2] / S.
        let mut a_is_ana = false;
        let mut a_is_bsp = false;
        let a_gcurve = edge_data(the_brep, &a_ne).representations[a_it].clone();
        let mut a_pc: [Option<Curve2d>; 2] = [None, None];
        let mut a_s: Option<Surface3> = None;
        if gcurve_is_curve_on_surface(&a_gcurve) {
            // OCCT L1362-1369: YaPCu = true; PC[0] = PCurve(); PCLoc; S =
            // Surface().  The rcad representation carries no pcurve
            // location (locations travel world-space) — the transform
            // branch is the recorded no-op.
            ya_pcu = true;
            a_pc[0] = gcurve_pcurve(&a_gcurve);
            a_s = gcurve_surface(the_brep, &a_gcurve);

            // OCCT L1371: GAS.Load(S).
            if let Some(a_s_val) = a_s.clone() {
                a_hs = GeomSurfaceAdaptor::new(a_s_val);
            }
            // OCCT L1372-1375: closed-surface PC[1].
            if gcurve_is_curve_on_closed_surface(&a_gcurve) {
                a_pc[1] = gcurve_pcurve2(&a_gcurve);
            }

            // OCCT L1377: the GCurve range (the SameRange source window).
            let (a_gf, a_gl) = gcurve_range(the_brep, &a_ne, &a_gcurve).unwrap_or((0.0, 0.0));
            // OCCT L1378: TolSameRange.
            let a_tol_same_range = a_hc.resolution(the_tolerance).max(PCONFUSION);

            // OCCT L1379-1381: for (i = 0; i < 2; i++).
            for a_i in 0..2usize {
                let mut a_update_pc = false;
                // OCCT L1381: curPC = PC[i]; L1383-1386: null check breaks.
                let mut a_cur_pc = match a_pc[a_i].clone() {
                    Some(a_c) => a_c,
                    None => break,
                };
                if !same_range_flag {
                    // OCCT L1387-1392: GeomLib::SameRange(TolSameRange, PC[i],
                    // GCurve->First(), GCurve->Last(), f3d, l3d, curPC);
                    // updatepc = (curPC != PC[i]) — the rcad value semantics
                    // always carry a fresh curve; writing the (possibly
                    // identical) value back is the OCCT same-handle no-op.
                    a_cur_pc =
                        same_range(a_tol_same_range, a_pc[a_i].as_ref().unwrap(), a_gf, a_gl, f3d, l3d);
                    a_update_pc = true;
                }
                let mut a_good_pc = true;
                // OCCT L1394: GAC2d.Load(curPC, f3d, l3d).
                a_hc2d = Geom2dCurveAdaptor::with_range(a_cur_pc.clone(), f3d, l3d);

                // OCCT L1396: double error = ComputeTol(HC, HC2d, HS, NCONTROL).
                let mut a_error = compute_tol(&a_hc, &a_hc2d, &a_hs, NCONTROL);

                // OCCT L1398-1402.
                if a_error > BIG_ERROR {
                    maxdist = a_error;
                    break;
                }

                // OCCT L1404-1623: the C0 BSpline block (the closing brace
                // is after the IsBad branch, as in OCCT).
                if a_hc2d.get_type() == CurveType::BSpline
                    && a_hc2d.continuity() == GeomAbsShape::C0
                {
                    // OCCT L1406-1409.
                    let a_u_resol = a_hs.u_resolution(the_tolerance);
                    let a_v_resol = a_hs.v_resolution(the_tolerance);
                    let mut a_tol_conf2d = a_u_resol.min(a_v_resol);
                    a_tol_conf2d = a_tol_conf2d.max(PCONFUSION);
                    // OCCT L1410-1415: bs2d = GAC2d.BSpline(); bs2dsov =
                    // bs2d; fC0/lC0; repar = true; OriginPoint = D0(fC0).
                    let a_legacy = a_hc2d.bspline().expect("Geom2dAdaptor_Curve::BSpline");
                    // Architecture note: the legacy kernel Curve2d encoding
                    // carries no periodic flag — the OCCT-faithful carrier
                    // is rebuilt non-periodic (the offset-pipeline pcurves
                    // are non-periodic; the SetOrigin arms below stay
                    // structural).
                    let mut a_bs2d = Geom2dBSplineCurve::from_bspline2(&a_legacy, false);
                    let a_bs2dsov = a_bs2d.clone();
                    let a_fc0 = a_bs2d.first_parameter();
                    let a_lc0 = a_bs2d.last_parameter();
                    let mut a_repar = true;
                    let mut a_origin_point = a_bs2d.eval_d0(a_fc0);
                    // OCCT L1416.
                    a_bs2d = c0_bspline_to_c1_bspline_curve(&a_bs2d, a_tol_conf2d);
                    a_is_bsp = true;

                    // OCCT L1419-1442: the periodic SetOrigin pass.  The
                    // rebuilt rcad carrier is non-periodic (arch. note
                    // above) — the arm is preserved with its OCCT anchor.
                    if a_bs2d.is_periodic() {
                        let a_new_origin_point = a_bs2d.eval_d0(a_bs2d.first_parameter());
                        if (a_origin_point.x - a_new_origin_point.x).abs() > PCONFUSION
                            || (a_origin_point.y - a_new_origin_point.y).abs() > PCONFUSION
                        {
                            for a_index in 1..=a_bs2d.nb_knots() {
                                let a_knot_point = a_bs2d.eval_d0(a_bs2d.knot(a_index));
                                if (a_origin_point.x - a_knot_point.x).abs() > PCONFUSION
                                    || (a_origin_point.y - a_knot_point.y).abs() > PCONFUSION
                                {
                                    continue;
                                }
                                a_bs2d.set_origin(a_index);
                                break;
                            }
                        }
                    }

                    // OCCT L1444-1507: the tolbail branch.
                    if matches!(a_bs2d.continuity(), SmoothShapeQ::C0) {
                        let mut a_tolbail = 0.0f64;
                        if eval_tol(
                            // OCCT L1446: EvalTol(curPC, ...) — the current
                            // (SameRange-reparametrized) pcurve, not the
                            // original PC[i].
                            &a_cur_pc,
                            a_s.as_ref().unwrap(),
                            &a_hc,
                            the_tolerance,
                            &mut a_tolbail,
                        ) {
                            // OCCT L1449: bs2d = bs2dsov.
                            a_bs2d = a_bs2dsov.clone();
                            let a_u_resbail = a_hs.u_resolution(a_tolbail);
                            let a_v_resbail = a_hs.v_resolution(a_tolbail);
                            let mut a_tol2dbail = a_u_resbail.min(a_v_resbail);
                            // OCCT L1453.
                            a_origin_point = a_bs2d.eval_d0(a_bs2d.first_parameter());

                            // OCCT L1455-1465: the minimal pole distance.
                            let a_nbp = a_bs2d.nb_poles_curve();
                            let mut a_p = a_bs2d.pole(1);
                            let mut a_d = f64::INFINITY;
                            for a_ip in 2..=a_nbp {
                                let a_p1 = a_bs2d.pole(a_ip);
                                a_d = a_d.min(a_p.distance_squared(a_p1));
                                a_p = a_p1;
                            }
                            let a_d = a_d.sqrt() * 0.1;

                            // OCCT L1467.
                            a_tol2dbail = a_tol2dbail.min(a_d).max(a_tol_conf2d);

                            // OCCT L1469.
                            a_bs2d = c0_bspline_to_c1_bspline_curve(&a_bs2d, a_tol2dbail);

                            // OCCT L1471-1494: the periodic SetOrigin pass
                            // (structural, see the arch. note above).
                            if a_bs2d.is_periodic() {
                                let a_new_origin_point =
                                    a_bs2d.eval_d0(a_bs2d.first_parameter());
                                if (a_origin_point.x - a_new_origin_point.x).abs() > PCONFUSION
                                    || (a_origin_point.y - a_new_origin_point.y).abs() > PCONFUSION
                                {
                                    for a_index in 1..=a_bs2d.nb_knots() {
                                        let a_knot_point =
                                            a_bs2d.eval_d0(a_bs2d.knot(a_index));
                                        if (a_origin_point.x - a_knot_point.x).abs() > PCONFUSION
                                            || (a_origin_point.y - a_knot_point.y).abs()
                                                > PCONFUSION
                                        {
                                            continue;
                                        }
                                        a_bs2d.set_origin(a_index);
                                        break;
                                    }
                                }
                            }

                            // OCCT L1496-1501.
                            if matches!(a_bs2d.continuity(), SmoothShapeQ::C0) {
                                a_good_pc = true;
                                a_bs2d = a_bs2dsov.clone();
                                a_repar = false;
                            }
                        } else {
                            // OCCT L1503-1506.
                            a_good_pc = false;
                        }
                    }

                    // OCCT L1509-1535: the repar block.
                    if a_good_pc && a_repar {
                        // OCCT L1513-1516: Knots = bs2d->Knots();
                        // Reparametrize(fC0, lC0, Knots); bs2d->SetKnots(Knots).
                        let mut a_knots: Vec<f64> =
                            (1..=a_bs2d.nb_knots()).map(|a_k| a_bs2d.knot(a_k)).collect();
                        bspl_lib::reparametrize(a_fc0, a_lc0, &mut a_knots);
                        a_bs2d = geom2d_bspline_set_knots(&a_bs2d, &a_knots);
                        // OCCT L1517-1518: GAC2d.Load(bs2d, f3d, l3d); curPC = bs2d.
                        a_hc2d = Geom2dCurveAdaptor::with_range(
                            geom2d_bspline_to_curve2d(&a_bs2d),
                            f3d,
                            l3d,
                        );
                        a_cur_pc = geom2d_bspline_to_curve2d(&a_bs2d);
                        let a_update_pcsov = a_update_pc;
                        a_update_pc = true;

                        // OCCT L1522-1534.
                        let a_error1 = compute_tol(&a_hc, &a_hc2d, &a_hs, NCONTROL);
                        if a_error1 > a_error {
                            a_bs2d = a_bs2dsov.clone();
                            a_hc2d = Geom2dCurveAdaptor::with_range(
                                geom2d_bspline_to_curve2d(&a_bs2d),
                                f3d,
                                l3d,
                            );
                            a_cur_pc = geom2d_bspline_to_curve2d(&a_bs2d);
                            a_update_pc = a_update_pcsov;
                            a_is_ana = true;
                        } else {
                            a_error = a_error1;
                        }
                    }

                    // OCCT L1537-1582: the IsBad detection.
                    let mut a_cont = smooth_to_geomabs(a_bs2d.continuity());
                    let mut a_is_bad = false;

                    if (a_cont as u8) > (GeomAbsShape::C0 as u8)
                        && a_error > (1.0e-3f64).max(the_tolerance)
                    {
                        let a_nb_knots = a_bs2d.nb_knots();
                        let a_critratio = 10.0f64;
                        let mut a_dtprev = a_bs2d.knot(2) - a_bs2d.knot(1);
                        let mut a_dtratio = 1.0f64;
                        let mut a_dtmin = a_dtprev;
                        let mut a_dtcur = 0.0f64;
                        for a_j in 2..a_nb_knots {
                            a_dtcur = a_bs2d.knot(a_j + 1) - a_bs2d.knot(a_j);
                            a_dtmin = a_dtmin.min(a_dtcur);

                            if a_is_bad {
                                continue;
                            }

                            if a_dtcur > a_dtprev {
                                a_dtratio = a_dtcur / a_dtprev;
                            } else {
                                a_dtratio = a_dtprev / a_dtcur;
                            }
                            if a_dtratio > a_critratio {
                                a_is_bad = true;
                            }
                            a_dtprev = a_dtcur;
                        }
                        if a_is_bad {
                            // To avoid failures in Approx_CurvilinearParameter
                            // (OCCT L1576): bs2d->Resolution(..., dtcur).
                            a_dtcur =
                                geom2d_bspline_resolution(&a_bs2d, (1.0e-3f64).max(the_tolerance));
                            if a_dtmin < a_dtcur {
                                a_is_bad = false;
                            }
                        }
                    }

                    // OCCT L1584-1621: the IsBad branch — reparametrize the
                    // pcurve by its curve length.
                    if a_is_bad {
                        if (a_cont as u8) > (GeomAbsShape::C2 as u8) {
                            a_cont = GeomAbsShape::C2;
                        }
                        let mut a_maxdeg = a_bs2d.degree() as i32;
                        if a_maxdeg == 1 {
                            a_maxdeg = 14;
                        }
                        let mut a_app_cur_par = ApproxCurvilinearParameter::new_curve_on_surface(
                            Arc::new(Geom2dCurveAdaptor::with_range(a_cur_pc.clone(), f3d, l3d)),
                            Arc::new(a_hs.clone()),
                            (1.0e-3f64).max(the_tolerance),
                            a_cont,
                            a_maxdeg,
                            10,
                        );
                        if a_app_cur_par.is_done() || a_app_cur_par.has_result() {
                            a_bs2d = a_app_cur_par.curve2d1().unwrap().clone();
                            a_hc2d = Geom2dCurveAdaptor::with_range(
                                geom2d_bspline_to_curve2d(&a_bs2d),
                                f3d,
                                l3d,
                            );
                            a_cur_pc = geom2d_bspline_to_curve2d(&a_bs2d);

                            // OCCT L1610-1619.
                            if (a_bs2d.first_parameter() - a_fc0).abs() > a_tol_same_range
                                || (a_bs2d.last_parameter() - a_lc0).abs() > a_tol_same_range
                            {
                                let mut a_knots: Vec<f64> = (1..=a_bs2d.nb_knots())
                                    .map(|a_k| a_bs2d.knot(a_k))
                                    .collect();
                                bspl_lib::reparametrize(a_fc0, a_lc0, &mut a_knots);
                                a_bs2d = geom2d_bspline_set_knots(&a_bs2d, &a_knots);
                                a_hc2d = Geom2dCurveAdaptor::with_range(
                                    geom2d_bspline_to_curve2d(&a_bs2d),
                                    f3d,
                                    l3d,
                                );
                                a_cur_pc = geom2d_bspline_to_curve2d(&a_bs2d);
                            }
                        }
                    }
                }

                // OCCT L1625-1701: the Approx_SameParameter tail.
                if a_good_pc {
                    // OCCT L1628.
                    let a_tol = if a_is_ana && a_is_bsp { 1.0e-7 } else { the_tolerance };
                    // OCCT L1627-1631: Approx_SameParameter SameP(HC, HC2d,
                    // HS, aTol).
                    let mut a_same_p = ApproxSameParameter::new(
                        &a_hc.curve,
                        a_hc.first,
                        a_hc.last,
                        &a_cur_pc,
                        a_s.as_ref().unwrap(),
                        a_tol,
                    );

                    if a_same_p.is_same_parameter() {
                        // OCCT L1633-1647.
                        maxdist = maxdist.max(a_same_p.tol_reached());
                        if a_update_pc {
                            if a_i == 0 {
                                gcurve_set_pcurve(the_brep, &a_ne, a_it, Some(a_cur_pc.clone()));
                            } else {
                                gcurve_set_pcurve2(the_brep, &a_ne, a_it, Some(a_cur_pc.clone()));
                            }
                        }
                    } else if a_same_p.is_done() {
                        // OCCT L1648-1672.
                        let a_tol_reached = a_same_p.tol_reached();
                        if a_tol_reached <= a_error {
                            a_cur_pc = a_same_p.curve2d();
                            a_update_pc = true;
                            maxdist = maxdist.max(a_tol_reached);
                        } else {
                            maxdist = maxdist.max(a_error);
                        }
                        if a_update_pc {
                            if a_i == 0 {
                                gcurve_set_pcurve(the_brep, &a_ne, a_it, Some(a_cur_pc.clone()));
                            } else {
                                gcurve_set_pcurve2(the_brep, &a_ne, a_it, Some(a_cur_pc.clone()));
                            }
                        }
                    } else {
                        // OCCT L1673-1696: Approx_SameParameter has failed —
                        // re-run SameRange over the ORIGINAL pcurve and take
                        // it as-is.
                        a_cur_pc = same_range(
                            a_tol_same_range,
                            a_pc[a_i].as_ref().unwrap(),
                            a_gf,
                            a_gl,
                            f3d,
                            l3d,
                        );

                        if a_i == 0 {
                            gcurve_set_pcurve(the_brep, &a_ne, a_it, Some(a_cur_pc.clone()));
                        } else {
                            gcurve_set_pcurve2(the_brep, &a_ne, a_it, Some(a_cur_pc.clone()));
                        }

                        is_same_p = false;
                    }
                } else {
                    // OCCT L1698-1701.
                    is_same_p = false;
                }

                // OCCT L1703-1714 (OCC5898).
                if !is_same_p {
                    let a_prec_surf =
                        brep_check_prec_surface(&GeomAdaptorSurface::new(a_s.clone().unwrap()));
                    let a_cur_tol = an_edge_tol + a_prec_c3d.max(a_prec_surf);
                    if a_cur_tol >= a_error {
                        maxdist = maxdist.max(an_edge_tol);
                        is_same_p = true;
                    }
                }
            }
        }
        // OCCT L1717: It.Next().
    }

    // OCCT L1719-1720.
    builder_range(the_brep, &a_ne, f3d, l3d);
    builder_same_range(the_brep, &a_ne, true);
    if is_same_p {
        // OCCT L1721-1737: reduce the edge tolerance when every
        // representation was processed (except the plane-associated ones
        // not stored in the edge).
        if ya_pcu {
            // Avoid setting too small tolerances.
            maxdist = maxdist.max(CONFUSION);
            *the_new_tol = maxdist;
            // OCCT L1733-1734: aNTE->Modified(true); aNTE->Tolerance(maxdist).
            let a_nte = the_brep.edge_mut_inplace(a_ne.clone());
            a_nte.flags |= tshape_flags::MODIFIED;
            a_nte.tolerance = maxdist;
        }
        // OCCT L1736.
        builder_same_parameter(the_brep, &a_ne, true);
    }

    // OCCT L1739.
    a_ne
}

// =========================================================================
// File-local re-hosts (mirroring build_curves3d.rs; private there, so the
// engine carries its own copies — the same OCCT bodies).
// =========================================================================

/// The guarded edge-data read: the pool walk for an in-pool edge,
/// `edge_data_pool_free` for the offset-engine edges living outside the
/// pool (no OCCT counterpart — OCCT has a single representation).
fn edge_data<'a>(the_brep: &'a BRep, the_e: &'a Shape) -> &'a TEdgeData {
    if shape_is_in_pool(the_brep, the_e) {
        the_brep.edge(the_e.clone())
    } else {
        edge_data_pool_free(the_e).expect("BRep_Tool: edge without a TEdgeData")
    }
}

/// OCCT BRep_CurveRepresentation::IsCurveOnSurface().
fn gcurve_is_curve_on_surface(a_cr: &CurveRepresentation) -> bool {
    matches!(
        a_cr,
        CurveRepresentation::CurveOnSurface { .. } | CurveRepresentation::CurveOnClosedSurface { .. }
    )
}

/// OCCT BRep_CurveOnSurface::IsCurveOnClosedSurface().
fn gcurve_is_curve_on_closed_surface(a_cr: &CurveRepresentation) -> bool {
    matches!(a_cr, CurveRepresentation::CurveOnClosedSurface { .. })
}

/// OCCT BRep_GCurve::First()/Last() — the representation range.
fn gcurve_range(the_brep: &BRep, the_e: &Shape, a_cr: &CurveRepresentation) -> Option<(f64, f64)> {
    match a_cr {
        CurveRepresentation::Curve3D { .. } => {
            let a_ed = edge_data(the_brep, the_e);
            Some((a_ed.range[0], a_ed.range[1]))
        }
        CurveRepresentation::CurveOnSurface { range, .. } => Some((range[0], range[1])),
        CurveRepresentation::CurveOnClosedSurface { range, .. } => Some((range[0], range[1])),
        CurveRepresentation::CurveOn2Surfaces { .. } => None,
    }
}

/// OCCT BRep_CurveOnSurface::PCurve().
fn gcurve_pcurve(a_cr: &CurveRepresentation) -> Option<Curve2d> {
    match a_cr {
        CurveRepresentation::CurveOnSurface { pcurve, .. } => Some(pcurve.clone()),
        CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => Some(pcurve1.clone()),
        _ => None,
    }
}

/// OCCT BRep_CurveOnClosedSurface::PCurve2().
fn gcurve_pcurve2(a_cr: &CurveRepresentation) -> Option<Curve2d> {
    match a_cr {
        CurveRepresentation::CurveOnClosedSurface { pcurve2, .. } => Some(pcurve2.clone()),
        _ => None,
    }
}

/// OCCT BRep_GCurve::Surface() — the support surface resolved from the
/// owning face key (the same encoding as build_curves3d.rs L819-833).
fn gcurve_surface(the_brep: &BRep, a_cr: &CurveRepresentation) -> Option<Surface3> {
    let a_face_key = match a_cr {
        CurveRepresentation::CurveOnSurface { face, .. } => *face,
        CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
        _ => return None,
    };
    let a_ts = the_brep
        .tshapes
        .iter()
        .find(|a_ts| std::sync::Arc::as_ptr(a_ts) as u64 == a_face_key.0)?;
    match a_ts.as_ref() {
        TShape::Face(a_fd) => a_fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_CurveOnSurface::PCurve(C) setter — replaces the pcurve at the
/// iterator slot `the_index` (and the pcurves-map entry under its face key).
fn gcurve_set_pcurve(the_brep: &mut BRep, the_e: &Shape, the_index: usize, the_pc: Option<Curve2d>) {
    let the_pc = match the_pc {
        Some(the_pc) => the_pc,
        None => return,
    };
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    let a_face_key = match &a_ed.representations[the_index] {
        CurveRepresentation::CurveOnSurface { face, .. } => *face,
        CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
        _ => return,
    };
    match &mut a_ed.representations[the_index] {
        CurveRepresentation::CurveOnSurface { pcurve, .. } => *pcurve = the_pc.clone(),
        CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => *pcurve1 = the_pc.clone(),
        _ => {}
    }
    if let Some(a_entry) = a_ed.pcurves.get_mut(&a_face_key) {
        a_entry.0 = the_pc;
    }
}

/// OCCT BRep_CurveOnClosedSurface::PCurve2(C) setter.
fn gcurve_set_pcurve2(
    the_brep: &mut BRep,
    the_e: &Shape,
    the_index: usize,
    the_pc: Option<Curve2d>,
) {
    let the_pc = match the_pc {
        Some(the_pc) => the_pc,
        None => return,
    };
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    match &mut a_ed.representations[the_index] {
        CurveRepresentation::CurveOnClosedSurface { pcurve2, .. } => *pcurve2 = the_pc,
        _ => {}
    }
}

/// OCCT BRep_Builder::Range(E, First, Last) — the 3d range of the edge.
fn builder_range(the_brep: &mut BRep, the_e: &Shape, the_first: f64, the_last: f64) {
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    a_ed.range = [the_first, the_last];
}

/// OCCT BRep_Builder::SameParameter(E, B).
fn builder_same_parameter(the_brep: &mut BRep, the_e: &Shape, the_b: bool) {
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    a_ed.same_parameter = the_b;
}

/// OCCT BRep_Builder::SameRange(E, B).
fn builder_same_range(the_brep: &mut BRep, the_e: &Shape, the_b: bool) {
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    a_ed.same_range = the_b;
}

// =========================================================================
// Geom2d_BSplineCurve value-semantics bridges (architecture notes inline).
// =========================================================================

/// The OCCT GeomAbs_Shape ordering on the rcad Geom2dBSplineCurve::Continuity
/// payload (SmoothShape carries no G1/G2 and its declaration order is not
/// the OCCT ordinal — the math GeomAbsShape is the ordered carrier).
fn smooth_to_geomabs(
    the_s: rcad_kernel::base::geom2d_convert::bspline_curve::SmoothShape,
) -> GeomAbsShape {
    match the_s {
        rcad_kernel::base::geom2d_convert::bspline_curve::SmoothShape::C0 => GeomAbsShape::C0,
        rcad_kernel::base::geom2d_convert::bspline_curve::SmoothShape::C1 => GeomAbsShape::C1,
        rcad_kernel::base::geom2d_convert::bspline_curve::SmoothShape::C2 => GeomAbsShape::C2,
        rcad_kernel::base::geom2d_convert::bspline_curve::SmoothShape::C3 => GeomAbsShape::C3,
        rcad_kernel::base::geom2d_convert::bspline_curve::SmoothShape::CN => GeomAbsShape::CN,
    }
}

/// OCCT Geom2d_BSplineCurve::SetKnots(K) (Geom2d_BSplineCurve.cxx L930-936)
/// — CheckCurveData + myKnots copy + updateKnots; the rcad value encoding
/// rebuilds the curve over the new knot table (poles/weights/mults/degree/
/// periodicity preserved).
fn geom2d_bspline_set_knots(the_bs: &Geom2dBSplineCurve, the_knots: &[f64]) -> Geom2dBSplineCurve {
    let a_poles: Vec<glam::DVec2> = (1..=the_bs.nb_poles_curve())
        .map(|a_i| the_bs.pole(a_i))
        .collect();
    let a_mults: Vec<i32> = (1..=the_bs.nb_knots())
        .map(|a_i| the_bs.multiplicity(a_i))
        .collect();
    if the_bs.is_rational() {
        let a_weights: Vec<f64> = (1..=the_bs.nb_poles_curve())
            .map(|a_i| the_bs.weight(a_i))
            .collect();
        Geom2dBSplineCurve::new_rational(
            a_poles,
            a_weights,
            the_knots.to_vec(),
            a_mults,
            the_bs.degree(),
            the_bs.is_periodic(),
        )
    } else {
        Geom2dBSplineCurve::new(
            a_poles,
            the_knots.to_vec(),
            a_mults,
            the_bs.degree(),
            the_bs.is_periodic(),
        )
    }
}

/// OCCT Geom2d_BSplineCurve -> the legacy kernel Curve2d encoding (the
/// adaptor read bridge; the same packing as the kernel c1_concat
/// `to_bspline2` helper).
fn geom2d_bspline_to_curve2d(the_bs: &Geom2dBSplineCurve) -> Curve2d {
    let mut a_flat_knots: Vec<f64> = Vec::new();
    for a_i in 1..=the_bs.nb_knots() {
        let a_k = the_bs.knot(a_i);
        let a_m = the_bs.multiplicity(a_i);
        for _ in 0..a_m {
            a_flat_knots.push(a_k);
        }
    }
    Curve2d::BSpline(rcad_kernel::geom::BSplineCurve2 {
        degree: the_bs.degree(),
        knots: a_flat_knots,
        control_points: (1..=the_bs.nb_poles_curve())
            .map(|a_i| the_bs.pole(a_i))
            .collect(),
        weights: (1..=the_bs.nb_poles_curve())
            .map(|a_i| the_bs.weight(a_i))
            .collect(),
        // Geom2d_BSplineCurve::IsPeriodic() — carried from the kernel curve.
        is_periodic: the_bs.is_periodic(),
    })
}

/// OCCT Geom2d_BSplineCurve::Resolution(Tolerance3D, U) — the BSplCLib
/// case-2 `Resolution` arm (BSplCLib.cxx L4329-4445): the per-pole value
/// is the MANHATTAN sum of the two coordinate differences maximized over
/// the wrapped pole pairs (rational: normalized by the minimal weight),
/// then UTolerance = Tolerance3D / (Degree * max).  This is the JOINT 2d
/// result — NOT the min of two dim-1 resolutions (the joint
/// max_derivative is at least each per-axis one, so the joint parametric
/// tolerance is at most each per-axis one).
fn geom2d_bspline_resolution(the_bs: &Geom2dBSplineCurve, the_tolerance3d: f64) -> f64 {
    let mut a_flat_knots: Vec<f64> = Vec::new();
    for a_i in 1..=the_bs.nb_knots() {
        let a_k = the_bs.knot(a_i);
        for _ in 0..the_bs.multiplicity(a_i) {
            a_flat_knots.push(a_k);
        }
    }
    // The OCCT flat PA layout: x0, y0, x1, y1, ...
    let mut a_poles_xy: Vec<f64> = Vec::with_capacity(2 * the_bs.nb_poles_curve() as usize);
    for a_i in 1..=the_bs.nb_poles_curve() {
        let a_p = the_bs.pole(a_i);
        a_poles_xy.push(a_p.x);
        a_poles_xy.push(a_p.y);
    }
    let a_weights: Option<Vec<f64>> = if the_bs.is_rational() {
        Some(
            (1..=the_bs.nb_poles_curve())
                .map(|a_i| the_bs.weight(a_i))
                .collect(),
        )
    } else {
        None
    };
    bspl_lib::resolution(
        2,
        &a_poles_xy,
        a_weights.as_deref(),
        &a_flat_knots,
        the_bs.degree(),
        the_tolerance3d,
    )
}

// =========================================================================
// Tests (hand-derived; see the module docs).
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use glam::{DVec2, DVec3};
    use rcad_kernel::geom::{BSplineCurve2, BSplineCurve3, Line2d, Line3, Plane};

    /// Builds a pool edge: a straight X segment (BSpline degree 1 on
    /// [0, 2]) on the plane z = 0 whose pcurve is the SAME segment in UV
    /// over the SAME range — the exact-sameparameter construction, with the
    /// sameparameter/samerange flags reset (the OCCT tolerance pass of
    /// BRepOffsetAPI_ThruSections L1053-1055 does exactly that before
    /// calling SameParameter).
    fn make_exact_edge(brep: &mut BRep) -> Shape {
        let a_face = brep.add_tface(
            Some(Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))),
            Shape::null(),
            Vec::new(),
            None,
            None,
            Vec::new(),
            false,
        );
        let a_face_key = (a_face.ptr_id(), a_face.location);
        let v1 = brep.add_tvertex(DVec3::ZERO);
        let v2 = brep.add_tvertex(DVec3::new(2.0, 0.0, 0.0));
        let bs = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 2.0, 2.0],
            control_points: vec![DVec3::ZERO, DVec3::new(2.0, 0.0, 0.0)],
            weights: vec![1.0; 2],
            is_periodic: false,
        };
        let a_edge = brep.add_tedge(Some(Curve3::BSpline(bs)), v1, v2, [0.0, 2.0]);
        let a_pcurve = Curve2d::BSpline(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 2.0, 2.0],
            control_points: vec![DVec2::ZERO, DVec2::new(2.0, 0.0)],
            weights: vec![1.0; 2],
            is_periodic: false,
        });
        {
            let a_ed = brep.edge_mut_inplace(a_edge.clone());
            a_ed.representations.push(CurveRepresentation::CurveOnSurface {
                face: a_face_key,
                pcurve: a_pcurve,
                range: [0.0, 2.0],
            });
            a_ed.same_parameter = false;
            a_ed.same_range = false;
        }
        a_edge
    }

    /// The engine on the exact straight edge: the edge comes out flagged
    /// sameparameter/samerange and the tolerance lands at the OCCT floor
    /// (max(maxdist, Precision::Confusion()), cxx L1731).
    #[test]
    fn engine_flags_exact_straight_edge() {
        let mut brep = BRep::new();
        let a_edge = make_exact_edge(&mut brep);
        let mut a_new_tol = -1.0f64;
        let a_ne = same_parameter_with_result(&mut brep, &a_edge, 1.0e-5, &mut a_new_tol, true);
        assert!(!a_ne.is_null(), "the engine must return the processed edge");
        let a_ed = brep.edge(a_ne.clone());
        assert!(a_ed.same_range, "SameRange flag");
        assert!(a_ed.same_parameter, "SameParameter flag");
        assert!(
            (a_ed.tolerance - CONFUSION).abs() < 1.0e-12,
            "edge tolerance {} == Precision::Confusion",
            a_ed.tolerance
        );
        // OCCT L1719: B.Range(aNE, f3d, l3d) — the 3d window [0, 2].
        assert_eq!(a_ed.range, [0.0, 2.0]);
        assert!((a_new_tol - CONFUSION).abs() < 1.0e-12);
    }

    /// A pcurve with a DIFFERENT parameter window ([10, 14] vs the 3d
    /// [0, 2]) on the same plane: the engine must rebuild it through
    /// GeomLib::SameRange onto [f3d, l3d] and flag the edge.
    #[test]
    fn engine_reparametrizes_offset_window_pcurve() {
        let mut brep = BRep::new();
        let a_face = brep.add_tface(
            Some(Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))),
            Shape::null(),
            Vec::new(),
            None,
            None,
            Vec::new(),
            false,
        );
        let a_face_key = (a_face.ptr_id(), a_face.location);
        let v1 = brep.add_tvertex(DVec3::ZERO);
        let v2 = brep.add_tvertex(DVec3::new(2.0, 0.0, 0.0));
        let bs = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 2.0, 2.0],
            control_points: vec![DVec3::ZERO, DVec3::new(2.0, 0.0, 0.0)],
            weights: vec![1.0; 2],
            is_periodic: false,
        };
        let a_edge = brep.add_tedge(Some(Curve3::BSpline(bs)), v1, v2, [0.0, 2.0]);
        // pcurve on [10, 14]: the abscissa s carries the same arc position
        // (s - 10) as the 3d curve at u = s - 10.
        let a_pcurve = Curve2d::BSpline(BSplineCurve2 {
            degree: 1,
            knots: vec![10.0, 10.0, 14.0, 14.0],
            control_points: vec![DVec2::ZERO, DVec2::new(2.0, 0.0)],
            weights: vec![1.0; 2],
            is_periodic: false,
        });
        {
            let a_ed = brep.edge_mut_inplace(a_edge.clone());
            a_ed.representations.push(CurveRepresentation::CurveOnSurface {
                face: a_face_key,
                pcurve: a_pcurve,
                range: [10.0, 14.0],
            });
            a_ed.same_parameter = false;
            a_ed.same_range = false;
        }

        let mut a_new_tol = -1.0f64;
        let a_ne = same_parameter_with_result(&mut brep, &a_edge, 1.0e-5, &mut a_new_tol, true);
        assert!(!a_ne.is_null());
        let a_ed = brep.edge(a_ne.clone());
        assert!(a_ed.same_range && a_ed.same_parameter);
        assert_eq!(a_ed.range, [0.0, 2.0]);
        assert!(a_new_tol >= CONFUSION, "theNewTol must be set");
        // The stored pcurve must now be defined on [0, 2] (SameRange
        // rebuild) and pass through the segment ends.  The representation
        // First/Last range is NOT rewritten by the OCCT engine (the
        // BRep_CurveOnSurface::PCurve(C) setter BRep_CurveOnSurface.cxx
        // L81-84 only swaps the handle; only the edge 3d range is set at
        // cxx L1719) — the rcad assertion mirrors that.
        let a_cr = &a_ed.representations[0];
        if let CurveRepresentation::CurveOnSurface { pcurve, range, .. } = a_cr {
            assert_eq!(*range, [10.0, 14.0], "the representation range stays");
            let a_p0 = Curve2dEval::point_at(pcurve, 0.0);
            let a_p1 = Curve2dEval::point_at(pcurve, 2.0);
            assert!((a_p0 - DVec2::ZERO).length() < 1.0e-6);
            assert!((a_p1 - DVec2::new(2.0, 0.0)).length() < 1.0e-6);
        } else {
            panic!("expected a curve-on-surface representation");
        }
        let _ = Line2d::new(DVec2::ZERO, DVec2::X);
        let _ = Line3::new(DVec3::ZERO, DVec3::X);
    }

    /// The OCCT empty-copied arm: with IsUseOldEdge = false the engine
    /// processes a fresh copy (cxx L1281) and leaves the original flags
    /// untouched.
    #[test]
    fn engine_copies_edge_when_not_using_old_edge() {
        let mut brep = BRep::new();
        let a_edge = make_exact_edge(&mut brep);
        let mut a_new_tol = -1.0f64;
        let a_ne = same_parameter_with_result(&mut brep, &a_edge, 1.0e-5, &mut a_new_tol, false);
        assert!(!a_ne.is_null());
        assert!(a_ne.ptr_id() != a_edge.ptr_id(), "a fresh copy must be returned");
        let a_ed_new = brep.edge(a_ne.clone());
        assert!(a_ed_new.same_parameter && a_ed_new.same_range);
        let a_ed_old = brep.edge(a_edge.clone());
        assert!(!a_ed_old.same_parameter, "the original edge stays unflagged");
        // The copy carries the vertices (cxx L1285-1289).
        assert!(!a_ed_new.first.is_null() && !a_ed_new.last.is_null());
    }

    /// The void overload raises the vertex tolerances to theNewTol
    /// (cxx L1241-1246 through UpdateVTol).
    #[test]
    fn void_overload_updates_vertices() {
        let mut brep = BRep::new();
        let a_edge = make_exact_edge(&mut brep);
        let a_v_first = brep.edge(a_edge.clone()).first.clone();
        let a_tol_before = rcad_kernel::topo::topods::BRepTool::vertex_tolerance(&brep, &a_v_first);
        same_parameter(&mut brep, &a_edge, 1.0e-5);
        let a_ed = brep.edge(a_edge.clone());
        assert!(a_ed.same_parameter);
        let a_tol_after =
            rcad_kernel::topo::topods::BRepTool::vertex_tolerance(&brep, &a_v_first);
        assert!(
            a_tol_after >= a_tol_before && a_tol_after >= CONFUSION - 1.0e-15,
            "vertex tolerance raised: {} -> {}",
            a_tol_before,
            a_tol_after
        );
    }

    /// The 2d Resolution joint-norm form (BSplCLib.cxx L4413-4445): a
    /// degree-1 curve with poles (0,0), (1,100) answers
    /// Tolerance3D/(Degree*(|dx|+|dy|)) = Tolerance3D/101 — the Manhattan
    /// sum over both coordinates, NOT the per-axis min (Tolerance3D/100
    /// from the y extent) the retired min-of-two-dim-1 body produced.
    #[test]
    fn geom2d_bspline_resolution_is_the_joint_manhattan_norm() {
        let a_bs = Geom2dBSplineCurve::new(
            vec![DVec2::ZERO, DVec2::new(1.0, 100.0)],
            vec![0.0, 1.0],
            vec![2, 2],
            1,
            false,
        );
        let a_u = geom2d_bspline_resolution(&a_bs, 1.0e-3);
        let a_expected = 1.0e-3 / 101.0;
        assert!(
            (a_u - a_expected).abs() < 1.0e-18,
            "joint 2d Resolution: {} vs {}",
            a_u,
            a_expected
        );
        assert!(
            a_u < 1.0e-5,
            "the joint form is tighter than the per-axis min {}",
            1.0e-5
        );
    }
}
