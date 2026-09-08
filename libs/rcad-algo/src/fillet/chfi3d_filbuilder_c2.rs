//! OCCT ChFi3d_FilBuilder_C2.cxx (TKFillet/ChFi3d) — 1:1 translation.
//!
//! Source: C:/Users/lilu/works/OCCT/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_FilBuilder_C2.cxx (L1-1153).
//!
//! Coverage:
//!   - ToricRotule (L87-121, file static)
//!   - RemoveSD (L123-139, file static)
//!   - ChFi3d_FilBuilder::PerformTwoCorner (L143-1153, the override of the
//!     pure virtual ChFi3d_Builder::PerformTwoCorner)
//!
//! Architecture mapping: OCCT `ChFi3d_FilBuilder::PerformTwoCorner` reads
//! only members inherited from ChFi3d_Builder (myVDataMap, myDS, myEFMap,
//! myEVIMap, myListStripe, myRegul, done, tolapp3d, tol2d), so the rcad
//! override is a free function over `&mut ChFi3dBuilder` (the rcad base
//! struct).  The virtual dispatch from ChFi3d_Builder::PerformFilletOnVertex
//! (ChFi3d_Builder.cxx L868) lands in chfi3d.rs `perform_two_corner_pending`.
//!
//! Small ChFi3d_Builder_0.cxx helpers missing from the existing translation
//! (ChFi3d_Coefficient, ChFi3d_BuildPCurve, the ChFi3d_mkbound family and
//! the ChFi3d_ComputePCurv Geom_Curve overload) are translated here with
//! anchors.  The GeomFill filling leaf classes (GeomFill_Boundary,
//! GeomFill_BoundWithSurf, GeomFill_SimpleBound, GeomFill_ConstrainedFilling
//! ::Surface) and Geom2dConvert are outside TKFillet and pending — GAP
//! carriers below preserve the OCCT failure paths.

use std::sync::{Arc, RwLock};

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{BezierCurve2, Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};

use super::chfi3d::{chfi3d_index_of_surf_data, chfi3d_index_point_in_ds, next_side, topabs_reverse, ChFi3dBuilder};
use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_compute_arete, chfi3d_compute_pcurv_2pt, chfi3d_same_parameter,
    vec_angle, BRepAdaptorSurface,
};
use super::chfi3d_builder_2::BRepAdaptorCurve2d;
use super::chfi3d_builder_2b::{GeomFillBoundary, GeomFillConstrainedFilling};
use super::chfi3d_builder_cncrn::chfi3d_is_in_front;
use super::chfi3d_ds::TopOpeBRepDSCurve;
use super::chfi_ds::{ChFiDS_CommonPoint, ChFiDS_State, ChFiDSStripe, SharedStripe, SharedSurfData};

/// OCCT ChFi3d_Builder.hxx L196 — PerformTwoCornerbyInter is a non-virtual
/// ChFi3d_Builder method defined in ChFi3d_Builder_C2.cxx L102-1043 (the
/// base Builder file, outside this F1 batch).  GAP carrier: the OCCT failure
/// path (Standard_False = the caller falls through to the cornerdata
/// construction) is reported until the owning batch lands.
impl ChFi3dBuilder {
    pub fn perform_two_cornerby_inter(&mut self, _index: usize) -> bool {
        // GAP carrier: OCCT ChFi3d_Builder_C2.cxx L108 (ChFi3d_Builder::
        // PerformTwoCornerbyInter) not translated; returns false = the OCCT
        // not-done path consumed at ChFi3d_FilBuilder_C2.cxx L475/L501.
        false
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1451-1465 — ChFi3d_Coefficient.
// =========================================================================
pub(crate) fn chfi3d_coefficient(v3d: DVec3, d1u: DVec3, d1v: DVec3) -> (f64, f64) {
    let aa = d1u.dot(d1u);
    let bb = d1u.dot(d1v);
    let cc = d1v.dot(d1v);
    let dd = d1u.dot(v3d);
    let ee = d1v.dot(v3d);
    let delta = aa * cc - bb * bb;
    let du = (dd * cc - ee * bb) / delta;
    let dv = (aa * ee - bb * dd) / delta;
    (du, dv)
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1866-1898 — ChFi3d_BuildPCurve (the 2d
// point/tangent core; a 4-pole Bezier).
// =========================================================================
pub(crate) fn chfi3d_build_pc_pts(
    p1: DVec2,
    d1: DVec2,
    p2: DVec2,
    d2: DVec2,
    redresse: bool,
) -> Curve2d {
    // OCCT L1869-1871: gp_Vec2d vref(p1, p2); gp_Dir2d dref(vref).
    let vref = p2 - p1;
    let dref = vref.normalize_or_zero();
    let mref = vref.length();
    // OCCT L1872-1881: redresse.
    let mut d1 = d1;
    let mut d2 = d2;
    if redresse {
        if d1.dot(dref) < 0.0 {
            d1 = -d1;
        }
        if d2.dot(dref) > 0.0 {
            d2 = -d2;
        }
    }
    // OCCT L1883-1893: the cubic through the two points.
    let mut pol = [DVec2::ZERO; 4];
    pol[0] = p1;
    pol[3] = p2;
    let lambda1 = d2.dot(d1).abs().max(dref.dot(d1).abs());
    let lambda1 = (0.5 * mref * lambda1).max(1.0e-5);
    pol[1] = p1 + lambda1 * d1;
    let lambda2 = d1.dot(d2).abs().max(dref.dot(d2).abs());
    let lambda2 = (0.5 * mref * lambda2).max(1.0e-5);
    pol[2] = p2 + lambda2 * d2;
    Curve2d::Bezier(BezierCurve2 {
        control_points: pol.to_vec(),
        weights: vec![1.0, 1.0, 1.0, 1.0],
    })
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1900-1937 — ChFi3d_BuildPCurve (the 2d-tangent
// overload; resolution-scaled in, poles scaled back).
// =========================================================================
pub(crate) fn chfi3d_build_pc_uv(
    surf: &BRepAdaptorSurface,
    p1: DVec2,
    v1: DVec2,
    p2: DVec2,
    v2: DVec2,
    redresse: bool,
) -> Curve2d {
    // OCCT L1906-1917: resolution scaling.
    let ures = surf.u_resolution(1.0);
    let vres = surf.v_resolution(1.0);
    let invures = 1.0 / ures;
    let invvres = 1.0 / vres;
    let pp1 = DVec2::new(invures * p1.x, invvres * p1.y);
    let pp2 = DVec2::new(invures * p2.x, invvres * p2.y);
    let vv1 = DVec2::new(invures * v1.x, invvres * v1.y);
    let vv2 = DVec2::new(invures * v2.x, invvres * v2.y);
    // OCCT L1918-1919: gp_Dir2d d1(vv1), d2(vv2).
    let d1 = vv1.normalize_or_zero();
    let d2 = vv2.normalize_or_zero();
    // OCCT L1920-1930: the Bezier + pole rescale.
    let g2dc = chfi3d_build_pc_pts(pp1, d1, pp2, d2, redresse);
    if let Curve2d::Bezier(pc) = &g2dc {
        let nbp = pc.control_points.len();
        let mut poles = pc.control_points.clone();
        for ip in 0..nbp {
            poles[ip] = DVec2::new(ures * poles[ip].x, vres * poles[ip].y);
        }
        return Curve2d::Bezier(BezierCurve2 {
            control_points: poles,
            weights: vec![1.0; nbp],
        });
    }
    g2dc
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1939-1970 — ChFi3d_BuildPCurve (the 3d-tangent
// overload; projects the 3d tangents through ChFi3d_Coefficient).
// =========================================================================
pub(crate) fn chfi3d_build_pc_3d(
    surf: &BRepAdaptorSurface,
    p1: DVec2,
    v1: DVec3,
    p2: DVec2,
    v2: DVec3,
    redresse: bool,
) -> Curve2d {
    // OCCT L1942-1949.
    let (_, d1u1, d1v1) = surf.surface.derivatives(p1.x, p1.y);
    let (du1, dv1) = chfi3d_coefficient(v1, d1u1, d1v1);
    let vv1 = DVec2::new(du1, dv1);
    let (_, d1u2, d1v2) = surf.surface.derivatives(p2.x, p2.y);
    let (du2, dv2) = chfi3d_coefficient(v2, d1u2, d1v2);
    let vv2 = DVec2::new(du2, dv2);
    // OCCT L1950-1961: Vref + redresse.
    let pp1 = surf.surface.point_at(p1.x, p1.y);
    let pp2 = surf.surface.point_at(p2.x, p2.y);
    let vref = pp2 - pp1;
    let mut vv1 = vv1;
    let mut vv2 = vv2;
    if redresse {
        if vref.dot(v1) < 0.0 {
            vv1 = -vv1;
        }
        if vref.dot(v2) > 0.0 {
            vv2 = -vv2;
        }
    }
    // OCCT L1963.
    chfi3d_build_pc_uv(surf, p1, vv1, p2, vv2, false)
}

// =========================================================================
// OCCT ChFi3d_Builder_0.hxx L163-172 / Builder_0.cxx L1771-1788 —
// ChFi3d_mkbound (Adaptor3d_Surface + 2d-tangent pair overload).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub(crate) fn chfi3d_mkbound_c2d_tangents(
    fac: &BRepAdaptorSurface,
    curv: &mut Option<Curve2d>,
    sens1: i32,
    pfac1: DVec2,
    vfac1: DVec2,
    sens2: i32,
    pfac2: DVec2,
    vfac2: DVec2,
    t3d: f64,
    ta: f64,
) -> GeomFillBoundary {
    // OCCT L1773-1782: v1/v2 reversed per sens.
    let mut v1 = vfac1;
    if sens1 == 1 {
        v1 = -v1;
    }
    let mut v2 = vfac2;
    if sens2 == 1 {
        v2 = -v2;
    }
    // OCCT L1783: curv = ChFi3d_BuildPCurve(Fac, pfac1, v1, pfac2, v2, false).
    *curv = Some(chfi3d_build_pc_uv(fac, pfac1, v1, pfac2, v2, false));
    // OCCT L1784: return ChFi3d_mkbound(Fac, curv, t3d, ta).
    chfi3d_mkbound_pcurve(fac, curv.as_ref(), t3d, ta)
}

/// OCCT ChFi3d_Builder_0.hxx L185-191 / Builder_0.cxx L1809-1816 —
/// ChFi3d_mkbound (Geom_Surface + two 2d points overload; marshals to the
/// adaptor overload).
pub(crate) fn chfi3d_mkbound_surface_two_points(
    s: &Surface3,
    p1: DVec2,
    p2: DVec2,
    t3d: f64,
    ta: f64,
) -> GeomFillBoundary {
    // OCCT L1812-1815: HS = GeomAdaptor_Surface(s); mkbound(HS, p1, p2, ...).
    let hs = BRepAdaptorSurface::initialize_surface(s.clone());
    chfi3d_mkbound_adaptor_two_points(&hs, p1, p2, t3d, ta)
}

/// OCCT ChFi3d_Builder_0.hxx L192-198 / Builder_0.cxx L1818-1828 —
/// ChFi3d_mkbound (Adaptor3d_Surface + two 2d points overload; a 2-pole
/// Bezier pcurve is built and assigned).
pub(crate) fn chfi3d_mkbound_adaptor_two_points(
    hs: &BRepAdaptorSurface,
    p1: DVec2,
    p2: DVec2,
    t3d: f64,
    ta: f64,
) -> GeomFillBoundary {
    // OCCT L1821-1824: curv = Geom2d_BezierCurve(pol).
    let curv = Curve2d::Bezier(BezierCurve2 {
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    });
    // OCCT L1825: return ChFi3d_mkbound(HS, curv, t3d, ta, isfreeboundary).
    chfi3d_mkbound_pcurve(hs, Some(&curv), t3d, ta)
}

/// OCCT ChFi3d_Builder_0.hxx L199-204 / Builder_0.cxx L1830-1843 —
/// ChFi3d_mkbound (Adaptor3d_Surface + pcurve overload).  GAP leaf: the
/// GeomFill_BoundWithSurf / GeomFill_SimpleBound construction is outside
/// TKFillet and pending; the empty boundary carrier keeps the call sites
/// typed and the filling machinery reports not-done downstream.
pub(crate) fn chfi3d_mkbound_pcurve(
    _hs: &BRepAdaptorSurface,
    _curv: Option<&Curve2d>,
    _t3d: f64,
    _ta: f64,
) -> GeomFillBoundary {
    GeomFillBoundary
}

/// OCCT ChFi3d_Builder_0.hxx L205-212 / Builder_0.cxx L1845-1857 —
/// ChFi3d_mkbound (Adaptor3d_Surface + pcurve + two 2d points overload).
#[allow(clippy::too_many_arguments)]
pub(crate) fn chfi3d_mkbound_adaptor_pcurve_two_points(
    fac: &BRepAdaptorSurface,
    curv: &mut Option<Curve2d>,
    p1: DVec2,
    p2: DVec2,
    t3d: f64,
    ta: f64,
    isfreeboundary: bool,
) -> GeomFillBoundary {
    // OCCT L1848-1852: curv = Geom2d_BezierCurve(pol).
    *curv = Some(Curve2d::Bezier(BezierCurve2 {
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    }));
    // OCCT L1853.
    let _ = isfreeboundary;
    chfi3d_mkbound_pcurve(fac, curv.as_ref(), t3d, ta)
}

/// OCCT GeomFill_SimpleBound (TKGeomAlgo/GeomFill/GeomFill_SimpleBound.hxx
/// L40-90) — a GeomFill_Boundary over a trimmed 3d curve.  GAP leaf (the
/// GeomFill boundary hierarchy is pending); the carrier keeps the
/// construction site typed.
pub(crate) fn geom_fill_simple_bound(_hcurve: &Curve3, _tol3d: f64, _tol2d: f64) -> GeomFillBoundary {
    GeomFillBoundary
}

/// OCCT GeomFill_Boundary::Bounds(f, l) / D1(u, p, v) (GeomFill_Boundary.hxx)
/// — GAP leaf accessors over the pending boundary object.
impl GeomFillBoundary {
    pub fn bounds(&self, _fb: &mut f64, _lb: &mut f64) {}

    pub fn d1(&self, _u: f64, _p: &mut DVec3, _v: &mut DVec3) {}
}

/// OCCT Geom2dConvert::CurveToBSplineCurve (TKGeomBase/Geom2dConvert.cxx
/// L181-380) — GAP leaf (outside TKFillet, pending).  The panic marks the
/// carrier boundary on the pcpivot-is-not-a-bspline branch (C2 L641-648).
pub(crate) fn geom2d_convert_curve_to_bspline_curve(_trc: &Curve2d) -> Curve2d {
    panic!("GAP: Geom2dConvert::CurveToBSplineCurve pending (TKGeomBase/Geom2dConvert.cxx L181-380)");
}

/// OCCT Geom2dConvert::SplitBSplineCurve (TKGeomBase/Geom2dConvert.cxx
/// L102/L146) — GAP leaf (outside TKFillet, pending).
pub(crate) fn geom2d_convert_split_bspline_curve(_bspl: &Curve2d, _u1: f64, _u2: f64, _tol: f64) -> Curve2d {
    panic!("GAP: Geom2dConvert::SplitBSplineCurve pending (TKGeomBase/Geom2dConvert.cxx L102-176)");
}

/// OCCT BRepAdaptor_Surface::IsUClosed() — GAP leaf (the rcad adaptor in
/// chfi3d_builder_0.rs carries no closed flag yet); the false result keeps
/// the OCCT no-adjustment path of the OCC26173 block (C2 L956-971).
pub(crate) fn brep_adaptor_surface_is_u_closed(_s: &BRepAdaptorSurface) -> bool {
    false
}

/// OCCT ChFi3d_Builder_0.cxx L1669-1688 — ChFi3d_ComputePCurv (Geom_Curve
/// overload; marshals to the adaptor overload L1648-1660, i.e. the 2-point
/// pcurve plus ChFi3d_SameParameter).  Consumed by PerformTwoCorner L754 and
/// PerformThreeCorner L1163.
pub(crate) fn chfi3d_compute_pcurv_c3d(
    c3d: &Curve3,
    uv1: DVec2,
    uv2: DVec2,
    pcurv: &mut Curve2d,
    s: &Surface3,
    pardeb: f64,
    parfin: f64,
    tol3d: f64,
    tolreached: &mut f64,
    reverse: bool,
) {
    // OCCT L1677-1679: hs = GeomAdaptor_Surface(S); hc = GeomAdaptor_Curve(
    // C3d, Pardeb, Parfin) — the rcad SameParameter carrier takes the curve
    // and surface directly.
    *pcurv = chfi3d_compute_pcurv_2pt(uv1, uv2, pardeb, parfin, reverse);
    chfi3d_same_parameter(c3d, pcurv, s, tol3d, tolreached);
}

/// The oriented edges of a face (OCCT TopExp_Explorer(F, TopAbs_EDGE)
/// marshaling: the face outer wire + inner wires carry the oriented edges).
pub(crate) fn face_edges(brep: &rcad_kernel::topods::BRep, face: &Shape) -> Vec<Shape> {
    use rcad_kernel::topo::topods::TShape;
    let mut out = Vec::new();
    let Some(fd) = face.as_face() else {
        return out;
    };
    for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
        if let Some(ts) = brep.tshapes.get(w.index) {
            if let TShape::Wire(wd) = ts.as_ref() {
                out.extend(wd.edges.iter().cloned());
            }
        }
    }
    out
}

// =========================================================================
// OCCT ChFi3d_FilBuilder_C2.cxx L87-121 — ToricRotule (file static).
// Test if it is a particular case of torus routine.  Three planes with two
// constant incident fillets of the same radius and the third face
// perpendicular to two others are required.
// =========================================================================
fn toric_rotule(
    fac: &BRepAdaptorSurface,
    s1: &BRepAdaptorSurface,
    s2: &BRepAdaptorSurface,
    c1: &SharedStripe,
    c2: &SharedStripe,
) -> bool {
    use super::chfi3d_builder_0::GeomAbsSurfaceType;
    // OCCT L94: double tolesp = 1.e-7;
    let tolesp = 1.0e-7;

    // OCCT L96-97: down_cast<ChFiDS_FilSpine>(cN->Spine()).
    let sp1 = {
        let g = c1.read().expect("stripe lock");
        g.spine().and_then(|sp| sp.down_cast_fil().cloned())
    };
    let sp2 = {
        let g = c2.read().expect("stripe lock");
        g.spine().and_then(|sp| sp.down_cast_fil().cloned())
    };
    // OCCT L98-101.
    let (Some(sp1), Some(sp2)) = (sp1, sp2) else {
        return false;
    };
    // OCCT L102-105.
    if !sp1.is_constant() || !sp2.is_constant() {
        return false;
    }
    // OCCT L106-110.
    if fac.get_type() != GeomAbsSurfaceType::Plane
        || s1.get_type() != GeomAbsSurfaceType::Plane
        || s2.get_type() != GeomAbsSurfaceType::Plane
    {
        return false;
    }
    // OCCT L111-113: fac.Plane().Position().Direction() normals.  The
    // plane types were checked above (OCCT would throw in Plane()).
    let plane_normal = |s: &BRepAdaptorSurface| -> DVec3 {
        match &s.surface {
            Surface3::Plane(p) => p.normal,
            _ => panic!("Standard_NoSuchObject: ToricRotule plane expected"),
        }
    };
    let df = plane_normal(fac);
    let ds1 = plane_normal(s1);
    let ds2 = plane_normal(s2);
    // OCCT L114-117.
    if df.dot(ds1).abs() >= tolesp || df.dot(ds2).abs() >= tolesp {
        return false;
    }
    // OCCT L118-120.
    let r1 = sp1.radius();
    let r2 = sp2.radius();
    (r1 - r2).abs() < tolesp
}

// =========================================================================
// OCCT ChFi3d_FilBuilder_C2.cxx L123-139 — RemoveSD (file static).
// =========================================================================
fn remove_sd(stripe: &SharedStripe, num1: i32, num2: i32) {
    // OCCT L125-126: the sequence of surfdata of the stripe (1-based).
    let mut g = stripe.write().expect("stripe lock");
    let seq = g.change_set_of_surf_data();
    if seq.is_empty() {
        return;
    }
    if num1 == num2 {
        // OCCT L133: Seq.Remove(num1).
        seq.remove((num1 - 1) as usize);
    } else {
        // OCCT L137: Seq.Remove(num1, num2) — inclusive range.
        seq.drain(((num1 - 1) as usize)..(num2 as usize));
    }
}

// =========================================================================
// OCCT ChFi3d_FilBuilder_C2.cxx L143-1153 — PerformTwoCorner.
// =========================================================================
#[allow(clippy::too_many_lines)]
pub fn perform_two_corner(fb: &mut ChFi3dBuilder, index: usize) {
    let brep = fb.my_brep.clone();
    // OCCT L150: done = false;
    fb.done = false;
    // OCCT L151: const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);
    let vtx = fb.my_vdata_map.find_key(index).clone();
    // OCCT L152: TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
    let mut dstr = fb.my_ds.take().expect("DS");
    // OCCT L153-154: It.Initialize(myVDataMap(Index)).
    let vdata: Vec<SharedStripe> = fb.my_vdata_map.find_from_index(index).clone();
    // OCCT L155-168: st1, st2, Sens1, Sens2, Isd1, Isd2, sd1, sd2, SeqFil1,
    // SeqFil2, surf1, surf2, OkinterCC, Okvisavis, SameSide, IFaCo1, IFaCo2,
    // UIntPC1, UIntPC2, FaCo, E1, E2, nbsurf1, nbsurf2, deb/fin, parE1/2.
    let mut sens1 = 0i32;
    let mut sens2 = 0i32;
    let mut isd1;
    let mut isd2;
    let mut u_intpc1 = 0.0f64;
    let mut u_intpc2 = 0.0f64;
    let mut faco = Shape::null();
    let mut same_side = false;
    let mut ifa_co1 = 0i32;
    let mut ifa_co2 = 0i32;
    let mut okvisavis = false;
    // OCCT L416-417: TopoDS_Face FF1, FF2, F, FaPiv; pctrans = FORWARD.
    let mut fapiv = Shape::null();
    let mut pctrans = Orientation::Forward;
    // OCCT L375-379: resetcp1/resetcp2, pivot, yapiv.
    let mut resetcp1 = false;
    let mut resetcp2 = false;
    // OCCT L605/L587: PCurveOnFace / PCurveOnPiv out slots.
    let mut pcurve_on_face: Option<Curve2d> = None;
    let mut pcurve_on_piv: Option<Curve2d> = None;
    let mut bid = 1i32;

    // the first
    //---------- (OCCT L173-177)
    let st1 = vdata[0].clone();
    isd1 = {
        let st1g = st1.read().expect("stripe lock");
        chfi3d_index_of_surf_data(&vtx, &st1g, &mut sens1)
    };

    // the second
    //---------- (OCCT L179-191)
    let st2 = vdata[1].clone();
    if Arc::ptr_eq(&st2, &st1) {
        sens2 = -1;
        isd2 = {
            let st2g = st2.read().expect("stripe lock");
            st2g.set_of_surf_data().len() as i32
        };
    } else {
        isd2 = {
            let st2g = st2.read().expect("stripe lock");
            chfi3d_index_of_surf_data(&vtx, &st2g, &mut sens2)
        };
    }

    // If two edges to rounded are tangent GeomPlate is called
    // (OCCT L193-237)
    let e1 = {
        let st1g = st1.read().expect("stripe lock");
        let sp = st1g.spine().expect("null spine").base().clone();
        if sens1 == 1 {
            sp.edges(1).clone()
        } else {
            sp.edges(sp.nb_edges()).clone()
        }
    };
    let e2 = {
        let st2g = st2.read().expect("stripe lock");
        let sp = st2g.spine().expect("null spine").base().clone();
        if sens2 == 1 {
            sp.edges(1).clone()
        } else {
            sp.edges(sp.nb_edges()).clone()
        }
    };

    // OCCT L213-218: BRepAdaptor_Curve BCurv1/2 + BRep_Tool::Parameter +
    // BRepLProp_CLProps(BCurv, par, 1, 1.e-4).
    let par_e1 = brep_tool_parameter(&brep, &vtx, &e1);
    let par_e2 = brep_tool_parameter(&brep, &vtx, &e2);
    let (bc1, _, _) = crate::brep_algo::tool::brep_tool_curve(&e1).expect("edge without 3d curve");
    let (bc2, _, _) = crate::brep_algo::tool::brep_tool_curve(&e2).expect("edge without 3d curve");
    // OCCT L219-221: CL1.Tangent(dir1); CL2.Tangent(dir2).
    let mut dir1 = bc1.derivative_at(par_e1).normalize();
    let mut dir2 = bc2.derivative_at(par_e2).normalize();
    // OCCT L222-229.
    if sens1 == -1 {
        dir1 = -dir1;
    }
    if sens2 == -1 {
        dir2 = -dir2;
    }
    // OCCT L230-237: ang1 = Abs(dir1.Angle(dir2)).
    let ang1 = vec_angle(dir1, dir2).abs();
    if ang1 < std::f64::consts::PI / 180.0 {
        fb.perform_more_three_corner(index, 2);
        fb.done = true;
        fb.my_ds = Some(dstr);
        return;
    }

    // OCCT L239-254.
    let okintercc = chfi3d_is_in_front(
        &brep, &dstr, &st1, &st2, isd1, isd2, sens1, sens2, &mut u_intpc1, &mut u_intpc2,
        &mut faco, &mut same_side, &mut ifa_co1, &mut ifa_co2, &mut okvisavis, &vtx, true, false,
    );

    let mut trouve = false;
    if !okvisavis {
        // one is not limited to the first or the last surfdata
        // to find the opposing data (OCCT L260-293)
        let nbsurf1 = {
            let st1g = st1.read().expect("stripe lock");
            st1g.set_of_surf_data().len() as i32
        };
        let nbsurf2 = {
            let st2g = st2.read().expect("stripe lock");
            st2g.set_of_surf_data().len() as i32
        };
        let mut deb1 = 1i32;
        let mut deb2 = 1i32;
        let mut fin1 = 1i32;
        let mut fin2 = 1i32;
        if nbsurf1 != 1 {
            if sens1 == 1 {
                deb1 = 1;
                fin1 = 2;
            } else {
                deb1 = nbsurf1 - 1;
                fin1 = nbsurf1;
            }
        }
        if nbsurf2 != 1 {
            if sens2 == 1 {
                deb2 = 1;
                fin2 = 2;
            } else {
                deb2 = nbsurf2 - 1;
                fin2 = nbsurf2;
            }
        }

        // OCCT L295-320.
        let mut i1 = deb1;
        while i1 <= fin1 && !trouve {
            isd1 = i1;
            let mut i2 = deb2;
            while i2 <= fin2 && !trouve {
                isd2 = i2;
                let _unused = chfi3d_is_in_front(
                    &brep, &dstr, &st1, &st2, isd1, isd2, sens1, sens2, &mut u_intpc1,
                    &mut u_intpc2, &mut faco, &mut same_side, &mut ifa_co1, &mut ifa_co2,
                    &mut okvisavis, &vtx, true, false,
                );
                trouve = okvisavis;
                i2 += 1;
            }
            i1 += 1;
        }
        // OCCT L321-326.
        if !trouve {
            fb.perform_more_three_corner(index, 2);
            fb.done = true;
            fb.my_ds = Some(dstr);
            return;
        } else {
            // OCCT L329-344.
            if sens1 == 1 && isd1 != 1 {
                remove_sd(&st1, 1, 1);
            }
            if sens1 != 1 && isd1 != nbsurf1 {
                remove_sd(&st1, fin1, fin1);
            }
            if sens2 == 1 && isd2 != 1 {
                remove_sd(&st2, 1, 1);
            }
            if sens2 != 1 && isd2 != nbsurf2 {
                remove_sd(&st2, fin2, fin2);
            }
        }
        // OCCT L346-347.
        isd1 = {
            let st1g = st1.read().expect("stripe lock");
            chfi3d_index_of_surf_data(&vtx, &st1g, &mut sens1)
        };
        isd2 = {
            let st2g = st2.read().expect("stripe lock");
            chfi3d_index_of_surf_data(&vtx, &st2g, &mut sens2)
        };
    }
    // OCCT L350-356.
    let ifa_arc1 = 3 - ifa_co1;
    let ifa_arc2 = 3 - ifa_co2;
    let seqfil1: Vec<SharedSurfData> = {
        let st1g = st1.read().expect("stripe lock");
        st1g.set_of_surf_data().clone()
    };
    let seqfil2: Vec<SharedSurfData> = {
        let st2g = st2.read().expect("stripe lock");
        st2g.set_of_surf_data().clone()
    };
    let sd1 = seqfil1[(isd1 - 1) as usize].clone();
    let surf1: Surface3 = dstr
        .surface(sd1.read().expect("surfdata lock").surf())
        .surface
        .clone();
    let sd2 = seqfil2[(isd2 - 1) as usize].clone();
    let surf2: Surface3 = dstr
        .surface(sd2.read().expect("surfdata lock").surf())
        .surface
        .clone();
    let ofa_co = faco.orientation;

    // The concavities are analyzed and the opposite face and the
    // eventual intersection of 2 pcurves on this face are found.
    // (OCCT L358-367)
    let isfirst1 = sens1 == 1;
    let isfirst2 = sens2 == 1;
    let (stat1, stat2) = {
        let st1g = st1.read().expect("stripe lock");
        let st2g = st2.read().expect("stripe lock");
        (
            st1g.spine().expect("null spine").base().status(isfirst1),
            st2g.spine().expect("null spine").base().status(isfirst2),
        )
    };
    let c1biseau = stat1 == ChFiDS_State::AllSame;
    let c1rotule = stat1 == ChFiDS_State::OnSame && stat2 == ChFiDS_State::OnSame;

    // It is checked if the fillets have a commonpoint on a common arc.
    // This edge is the pivot of the bevel or the knee. (OCCT L369-393)
    let (cp1, cp2): (ChFiDS_CommonPoint, ChFiDS_CommonPoint) = {
        let mut sd1g = sd1.write().expect("surfdata lock");
        let mut sd2g = sd2.write().expect("surfdata lock");
        (
            sd1g.change_vertex(isfirst1, ifa_arc1).clone(),
            sd2g.change_vertex(isfirst2, ifa_arc2).clone(),
        )
    };

    let mut pivot = Shape::null();
    let mut yapiv = false;
    if cp1.is_on_arc() {
        // OCCT L380-383.
        pivot = cp1.arc().clone();
    } else {
        // OCCT L384-389.
        fb.perform_more_three_corner(index, 2);
        fb.done = true;
        fb.my_ds = Some(dstr);
        return;
    }
    if cp1.is_on_arc() && cp2.is_on_arc() {
        // OCCT L390-393.
        yapiv = pivot.is_same(cp2.arc());
    }
    // OCCT L394-405: Hpivot + the sameparam test on the pivot curve.
    let hpivot: Option<Curve3> = if yapiv {
        crate::brep_algo::tool::brep_tool_curve(&pivot).map(|(c, _, _)| c)
    } else {
        None
    };
    let parcp1 = cp1.parameter_on_arc();
    let parcp2 = cp2.parameter_on_arc();
    let mut sameparam = false;
    if yapiv {
        if let Some(pc) = &hpivot {
            // OCCT L402-404: tst1 = Hpivot->Value(parCP1) etc.
            let tst1 = pc.point_at(parcp1);
            let tst2 = pc.point_at(parcp2);
            sameparam = tst1.distance(tst2) <= fb.tolapp3d;
        }
    }
    // OCCT L406-414: HFaCo/HBRS1/HBRS2 handles; BRFaCo.Initialize(FaCo).
    let hfacov = BRepAdaptorSurface::initialize(&brep, &faco);

    // OCCT L416-428: FF1/FF2 from the DS; nullity -> PerformMoreThreeCorner;
    // BRS1/BRS2 initialization.
    let ff1 = dstr
        .shape(sd1.read().expect("surfdata lock").index_of(ifa_arc1))
        .clone();
    let ff2 = dstr
        .shape(sd2.read().expect("surfdata lock").index_of(ifa_arc2))
        .clone();
    if ff1.is_null() || ff2.is_null() {
        fb.perform_more_three_corner(index, 2);
        fb.done = true;
        fb.my_ds = Some(dstr);
        return;
    }
    let hbrs1 = BRepAdaptorSurface::initialize(&brep, &ff1);
    let hbrs2 = BRepAdaptorSurface::initialize(&brep, &ff2);

    if yapiv {
        // OCCT L430-452: the pivot must be shared by FF1 and FF2 through
        // myEFMap(pivot).
        let mut ok1 = false;
        let mut ok2 = false;
        let faces = if fb.my_ef_map.contains(&pivot) {
            fb.my_ef_map.find(&pivot).clone()
        } else {
            Vec::new()
        };
        for kt in &faces {
            if !ok1 && ff1.is_same(kt) {
                ok1 = true;
            }
            if !ok2 && ff2.is_same(kt) {
                ok2 = true;
            }
        }
        if !ok1 || !ok2 {
            fb.perform_more_three_corner(index, 2);
            fb.done = true;
            fb.my_ds = Some(dstr);
            return;
        }
    }

    // bevel
    //------ (OCCT L458-514)
    let mut done = false;
    if c1biseau {
        // OCCT L469.
        done = fb.perform_two_cornerby_inter(index);
        if !done {
            // OCCT L477-480.
            fb.perform_more_three_corner(index, 2);
            fb.done = true;
            fb.my_ds = Some(dstr);
            return;
        }
    } else if c1rotule {
        // save. (OCCT L483-491)
        let (cp11, cp12, cp21, cp22) = {
            let sd1g = sd1.read().expect("surfdata lock");
            let sd2g = sd2.read().expect("surfdata lock");
            (
                sd1g.vertex(isfirst1, 1).clone(),
                sd1g.vertex(isfirst1, 2).clone(),
                sd2g.vertex(isfirst2, 1).clone(),
                sd2g.vertex(isfirst2, 2).clone(),
            )
        };
        let (intf11, intf12, intf21, intf22) = {
            let sd1g = sd1.read().expect("surfdata lock");
            let sd2g = sd2.read().expect("surfdata lock");
            (
                sd1g.interference_on_s1().clone(),
                sd1g.interference_on_s2().clone(),
                sd2g.interference_on_s1().clone(),
                sd2g.interference_on_s2().clone(),
            )
        };
        // OCCT L496.
        done = fb.perform_two_cornerby_inter(index);
        if !done {
            // restore (OCCT L501-513)
            {
                let mut sd1g = sd1.write().expect("surfdata lock");
                let mut sd2g = sd2.write().expect("surfdata lock");
                *sd1g.change_vertex(isfirst1, 1) = cp11;
                *sd1g.change_vertex(isfirst1, 2) = cp12;
                *sd2g.change_vertex(isfirst2, 1) = cp21;
                *sd2g.change_vertex(isfirst2, 2) = cp22;
                *sd1g.change_interference_on_s1() = intf11;
                *sd1g.change_interference_on_s2() = intf12;
                *sd2g.change_interference_on_s1() = intf21;
                *sd2g.change_interference_on_s2() = intf22;
            }
            done = false;
        }
    }

    if !c1biseau && !done {
        // new cornerdata is created
        //------------------------------- (OCCT L518-525)
        let corner = Arc::new(RwLock::new(ChFiDSStripe::default()));
        let coin = Arc::new(RwLock::new(super::chfi_ds::ChFiDSSurfData::default()));
        {
            // OCCT L521-525: cornerset = new HSequence; Append(coin).
            let mut cg = corner.write().expect("stripe lock");
            let cornerset = cg.change_set_of_surf_data();
            cornerset.clear();
            cornerset.push(coin.clone());
        }

        if same_side {
            // OCCT L529: ToricRotule(BRFaCo, BRS1, BRS2, st1, st2).
            if toric_rotule(&hfacov, &hbrs1, &hbrs2, &st1, &st2) {
                // Direct construction.
                // --------------------- (OCCT L531-565)
                let mut ori = ofa_co;
                let ori_s = {
                    let st1g = st1.read().expect("stripe lock");
                    st1g.orientation_on_s(ifa_co1)
                };
                let mut off1 = ff1.orientation;
                let ori_sff1 = {
                    let st1g = st1.read().expect("stripe lock");
                    st1g.orientation_on_s(ifa_arc1)
                };
                // OCCT L539-540: bid = 1; bid = ChFi3d::NextSide(ori, OFF1,
                // oriS, oriSFF1, bid).
                bid = next_side(&mut ori, &mut off1, ori_s, ori_sff1, bid);
                let mut op1 = Orientation::Forward;
                let mut op2 = Orientation::Forward;
                if yapiv {
                    // OCCT L544.
                    bid = super::chfi3d::concave_side(&brep, &ff1, &ff2, &pivot, &mut op1, &mut op2);
                }
                // OCCT L546-547.
                op1 = topabs_reverse(op1);
                op2 = topabs_reverse(op2);
                // OCCT L551.
                let radius = {
                    let st1g = st1.read().expect("stripe lock");
                    st1g
                        .spine()
                        .and_then(|sp| sp.down_cast_fil().cloned())
                        .expect("fil spine")
                        .radius()
                };
                // OCCT L552-561: ChFiKPart_ComputeData::ComputeCorner(DStr,
                // coin, HFaCo, HBRS1, HBRS2, OFaCo, ori, op1, op2, radius) —
                // the toric rotule overload (ChFiKPart_ComputeData.cxx
                // L734-762); the rcad backend takes the face shapes behind
                // the adaptors (HFaCo->FaCo, HBRS1->FF1, HBRS2->FF2).
                done = super::chfi_kpart::compute_data_compute_corner_rotule(
                    &mut dstr,
                    &mut coin.write().expect("surfdata lock"),
                    &faco,
                    &ff1,
                    &ff2,
                    ofa_co,
                    ori,
                    op1,
                    op2,
                    radius,
                );
            } else {
                // Construction by filling remplissage
                // ---------------------------- (OCCT L566-682)
                let (u_pcarc1, p2da1, p2df1, p2dfac1, v2dfac1);
                let (u_pcarc2, p2da2, p2df2, p2dfac2, v2dfac2);
                {
                    let sd1g = sd1.read().expect("surfdata lock");
                    let sd2g = sd2.read().expect("surfdata lock");
                    // OCCT L574-577.
                    u_pcarc1 = sd1g.interference(ifa_arc1).parameter(isfirst1);
                    p2da1 = sd1g
                        .interference(ifa_arc1)
                        .pcurve_on_surf()
                        .expect("PCurveOnSurf")
                        .point_at(u_pcarc1);
                    p2df1 = sd1g
                        .interference(ifa_co1)
                        .pcurve_on_surf()
                        .expect("PCurveOnSurf")
                        .point_at(u_pcarc1);
                    let (p, v) = {
                        let pc = sd1g
                            .interference(ifa_co1)
                            .pcurve_on_face()
                            .expect("PCurveOnFace");
                        (pc.point_at(u_pcarc1), pc.derivative_at(u_pcarc1))
                    };
                    p2dfac1 = p;
                    v2dfac1 = v;
                    // OCCT L578-581.
                    u_pcarc2 = sd2g.interference(ifa_arc2).parameter(isfirst2);
                    p2da2 = sd2g
                        .interference(ifa_arc2)
                        .pcurve_on_surf()
                        .expect("PCurveOnSurf")
                        .point_at(u_pcarc2);
                    p2df2 = sd2g
                        .interference(ifa_co2)
                        .pcurve_on_surf()
                        .expect("PCurveOnSurf")
                        .point_at(u_pcarc2);
                    let (p, v) = {
                        let pc = sd2g
                            .interference(ifa_co2)
                            .pcurve_on_face()
                            .expect("PCurveOnFace");
                        (pc.point_at(u_pcarc2), pc.derivative_at(u_pcarc2))
                    };
                    p2dfac2 = p;
                    v2dfac2 = v;
                }
                // OCCT L585-586.
                let b1 = chfi3d_mkbound_surface_two_points(&surf1, p2df1, p2da1, fb.tolapp3d, 2.0e-4);
                let b2 = chfi3d_mkbound_surface_two_points(&surf2, p2df2, p2da2, fb.tolapp3d, 2.0e-4);
                // OCCT L587-597: Bfac = ChFi3d_mkbound(HFaCo, PCurveOnFace,
                // Sens1, p2dfac1, v2dfac1, Sens2, p2dfac2, v2dfac2, ...).
                let bfac = chfi3d_mkbound_c2d_tangents(
                    &hfacov,
                    &mut pcurve_on_face,
                    sens1,
                    p2dfac1,
                    v2dfac1,
                    sens2,
                    p2dfac2,
                    v2dfac2,
                    fb.tolapp3d,
                    2.0e-4,
                );
                // OCCT L598: GeomFill_ConstrainedFilling fil(8, 20).
                let mut fil = GeomFillConstrainedFilling::new(8, 20);
                if sameparam {
                    // OCCT L599-601: fil.Init(Bfac, B2, B1, true).
                    fil.init3(&bfac, &b2, &b1, true);
                } else {
                    // OCCT L603-608: HPivTrim + Bpiv + fil.Init(Bfac, B2,
                    // Bpiv, B1, true).
                    let hpivtrim = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
                        curve: Box::new(hpivot.clone().expect("pivot curve")),
                        first: parcp1.min(parcp2),
                        last: parcp1.max(parcp2),
                    });
                    let bpiv = geom_fill_simple_bound(&hpivtrim, fb.tolapp3d, 2.0e-4);
                    fil.init4(&bfac, &b2, &bpiv, &b1, true);
                    // OCCT L609-629: pivotverslebas decides FaPiv.
                    let pivotverslebas = {
                        // OCCT L610-615: Hpivot->D1(parCP1, bidon, dArc);
                        // B1->Bounds(fb1, lb1); B1->D1(lb1, bidon, dcf).
                        let d_arc = hpivot
                            .as_ref()
                            .expect("pivot curve")
                            .derivative_at(parcp1);
                        let mut fb1 = 0.0;
                        let mut lb1 = 0.0;
                        let mut bidon = DVec3::ZERO;
                        let mut dcf = DVec3::ZERO;
                        b1.bounds(&mut fb1, &mut lb1);
                        b1.d1(lb1, &mut bidon, &mut dcf);
                        // OCCT L616.
                        d_arc.dot(dcf) <= 0.0
                    };
                    // OCCT L617: pcfalenvers.
                    let pcfalenvers = parcp1 > parcp2;
                    if (pivotverslebas && !pcfalenvers) || (!pivotverslebas && pcfalenvers) {
                        // OCCT L620-623.
                        fapiv = ff2.clone();
                        resetcp2 = true;
                    } else {
                        // OCCT L625-629.
                        fapiv = ff1.clone();
                        resetcp1 = true;
                    }
                    // OCCT L630-631.
                    fapiv.orientation = Orientation::Forward;
                    let mut pcpivot = BRepAdaptorCurve2d::new();
                    pcpivot.initialize(&brep, &pivot, &fapiv);
                    // OCCT L632-640: the orientation of the pivot edge inside
                    // FaPiv.
                    let fapiv_edges = face_edges(&brep, &fapiv);
                    let n_edges = fapiv_edges.len();
                    for (i_ex, cur) in fapiv_edges.iter().enumerate() {
                        if cur.is_same(&pivot) {
                            pctrans = cur.orientation;
                            break;
                        } else if i_ex + 1 >= n_edges {
                            // OCCT: the explorer runs out (no break).
                            let _ = i_ex;
                        }
                    }
                    // OCCT L641-656: build PCurveOnPiv from the pcpivot.
                    {
                        let (base_pc, _, _) = brep
                            .curve_on_surface(&pivot, &fapiv)
                            .expect("pcpivot");
                        if !matches!(base_pc, Curve2d::BSpline(_)) {
                            // OCCT L643-647: Geom2d_TrimmedCurve + 
                            // Geom2dConvert::CurveToBSplineCurve.
                            let trc = Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                                curve: Box::new(base_pc.clone()),
                                t_min: parcp1.min(parcp2),
                                t_max: parcp1.max(parcp2),
                            });
                            pcurve_on_piv = Some(geom2d_convert_curve_to_bspline_curve(&trc));
                        } else {
                            // OCCT L649-656: Geom2dConvert::SplitBSplineCurve.
                            pcurve_on_piv = Some(geom2d_convert_split_bspline_curve(
                                &base_pc,
                                parcp1.min(parcp2),
                                parcp1.max(parcp2),
                                fb.tol2d,
                            ));
                        }
                    }
                    // OCCT L657-659: BSplCLib::Reparametrize(0, 1, kk).
                    if let Some(Curve2d::BSpline(bs)) = pcurve_on_piv.as_mut() {
                        let mut kk = bs.knots.clone();
                        rcad_kernel::math::bspl_lib::reparametrize(0.0, 1.0, &mut kk);
                        bs.knots = kk;
                    }
                    // OCCT L660-664.
                    if pcfalenvers {
                        if let Some(pc) = pcurve_on_piv.as_mut() {
                            *pc = rcad_kernel::geom::reverse_curve2d(pc);
                        }
                        pctrans = topabs_reverse(pctrans);
                    }
                }
                // OCCT L666-678: Surfcoin = fil.Surface(); done =
                // CompleteData(coin, Surfcoin, HFaCo, PCurveOnFace, HFaPiv,
                // PCurveOnPiv, OFaCo, true, false, false, false, false).
                done = match fil.surface() {
                    Some(surfcoin) => {
                        // The HFaPiv adaptor marshals the FaPiv face; when
                        // sameparam the OCCT passes the still-null HFaPiv and
                        // PC2 is None so the backend never reads it.
                        let hfapiv = if fapiv.is_null() {
                            None
                        } else {
                            Some(BRepAdaptorSurface::initialize(&brep, &fapiv))
                        };
                        let empty = BRepAdaptorSurface::initialize_surface(surfcoin.clone());
                        fb.complete_data_surfcoin(
                            &mut coin.write().expect("surfdata lock"),
                            &surfcoin,
                            &hfacov,
                            pcurve_on_face.as_ref(),
                            hfapiv.as_ref().unwrap_or(&empty),
                            pcurve_on_piv.as_ref(),
                            ofa_co,
                            true,
                            false,
                            false,
                            false,
                            false,
                        )
                    }
                    None => {
                        // GAP: GeomFill_ConstrainedFilling::Surface pending
                        // (chfi3d_builder_2b.rs) — the OCCT flow always
                        // produces a surface; report the not-done path.
                        false
                    }
                };
            }
            // OCCT L686-856: if (done) — update 3 CornerData and the DS.
            if done {
                // OCCT L690-705: the CP1/CP2 reset to a point-only
                // commonpoint.
                if resetcp1 {
                    let pjyl = cp1.point();
                    let tolsav = cp1.tolerance();
                    let mut sd1g = sd1.write().expect("surfdata lock");
                    let cp = sd1g.change_vertex(isfirst1, ifa_arc1);
                    cp.reset();
                    cp.set_point(pjyl);
                    cp.set_tolerance(tolsav);
                } else if resetcp2 {
                    let pjyl = cp2.point();
                    let tolsav = cp2.tolerance();
                    let mut sd2g = sd2.write().expect("surfdata lock");
                    let cp = sd2g.change_vertex(isfirst2, ifa_arc2);
                    cp.reset();
                    cp.set_point(pjyl);
                    cp.set_tolerance(tolsav);
                }
                // OCCT L706-713: Pf1/Pf2/Pl1/Pl2; Pf2 = CP1; Pl2 = CP2.
                let cp1_now = sd1
                    .read()
                    .expect("surfdata lock")
                    .vertex(isfirst1, ifa_arc1)
                    .clone();
                let cp2_now = sd2
                    .read()
                    .expect("surfdata lock")
                    .vertex(isfirst2, ifa_arc2)
                    .clone();
                let (pf1, pl1): (ChFiDS_CommonPoint, ChFiDS_CommonPoint);
                {
                    let mut coing = coin.write().expect("surfdata lock");
                    pf1 = coing.vertex_first_on_s1().clone();
                    pl1 = coing.vertex_last_on_s1().clone();
                    *coing.change_vertex_first_on_s2() = cp1_now.clone();
                    *coing.change_vertex_last_on_s2() = cp2_now.clone();
                }
                let pf2 = cp1_now;
                let pl2 = cp2_now;

                // the corner to start,
                // ----------------------- (OCCT L715-728)
                let if1 = chfi3d_index_point_in_ds(&pf1, &mut dstr);
                let if2 = chfi3d_index_point_in_ds(&pf2, &mut dstr);
                let il1 = chfi3d_index_point_in_ds(&pl1, &mut dstr);
                let il2 = if sameparam {
                    if2
                } else {
                    chfi3d_index_point_in_ds(&pl2, &mut dstr)
                };

                let coin_surf: Surface3 = dstr
                    .surface(coin.read().expect("surfdata lock").surf())
                    .surface
                    .clone();
                // OCCT L730-749.
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
                let (c3d, mut pcurv, p1deb, p2deb, mut tolreached) = chfi3d_compute_arete(
                    &brep, &pf1, pp1, &pf2, pp2, &coin_surf, fb.tolapp3d, fb.tol2d, 0,
                );
                // OCCT L750-763: par1 + pp1/pp2 + ChFi3d_ComputePCurv on the
                // st1 side.
                let (par1, pp1_s1, pp2_s1) = {
                    let sd1g = sd1.read().expect("surfdata lock");
                    let par1 = sd1g.interference(ifa_arc1).parameter(isfirst1);
                    (
                        par1,
                        sd1g
                            .interference(ifa_co1)
                            .pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(par1),
                        sd1g
                            .interference(ifa_arc1)
                            .pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(par1),
                    )
                };
                let sd1_surf: Surface3 = dstr
                    .surface(sd1.read().expect("surfdata lock").surf())
                    .surface
                    .clone();
                let mut tolr1 = 0.0f64;
                if let Some(c3d) = &c3d {
                    chfi3d_compute_pcurv_c3d(
                        c3d,
                        pp1_s1,
                        pp2_s1,
                        &mut pcurv,
                        &sd1_surf,
                        p1deb,
                        p2deb,
                        fb.tolapp3d,
                        &mut tolr1,
                        false,
                    );
                }
                tolreached = tolreached.max(tolr1);
                // OCCT L764-773: Tcurv1 + Icf + regdeb + the corner entries.
                let tcurv1 = TopOpeBRepDSCurve::new(c3d, tolreached);
                let icf = dstr.add_curve(tcurv1);
                let mut regdeb = super::chfi_ds::ChFiDSRegul::default();
                regdeb.set_curve(icf);
                regdeb.set_s1(coin.read().expect("surfdata lock").surf(), false);
                regdeb.set_s2(sd1.read().expect("surfdata lock").surf(), false);
                fb.my_regul.push(regdeb);
                {
                    let mut cornerg = corner.write().expect("stripe lock");
                    cornerg.change_first_curve(icf);
                    cornerg.change_first_parameters(p1deb, p2deb);
                    cornerg.change_index_first_point_on_s1(if1);
                    cornerg.change_index_first_point_on_s2(if2);
                }

                // OCCT L775-791: last-parameter ComputeArete.
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
                let (c3d, mut pcurv_l, p1fin, p2fin, mut tolreached) = chfi3d_compute_arete(
                    &brep, &pl1, pp1, &pl2, pp2, &coin_surf, fb.tolapp3d, fb.tol2d, 0,
                );
                // OCCT L792-805: par2 + ChFi3d_ComputePCurv on the st2 side.
                let (par2, pp1_s2, pp2_s2) = {
                    let sd2g = sd2.read().expect("surfdata lock");
                    let par2 = sd2g.interference(ifa_arc2).parameter(isfirst2);
                    (
                        par2,
                        sd2g
                            .interference(ifa_co2)
                            .pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(par2),
                        sd2g
                            .interference(ifa_arc2)
                            .pcurve_on_surf()
                            .expect("PCurveOnSurf")
                            .point_at(par2),
                    )
                };
                let sd2_surf: Surface3 = dstr
                    .surface(sd2.read().expect("surfdata lock").surf())
                    .surface
                    .clone();
                let mut tolr2 = 0.0f64;
                if let Some(c3d) = &c3d {
                    chfi3d_compute_pcurv_c3d(
                        c3d,
                        pp1_s2,
                        pp2_s2,
                        &mut pcurv_l,
                        &sd2_surf,
                        p1fin,
                        p2fin,
                        fb.tolapp3d,
                        &mut tolr2,
                        false,
                    );
                }
                tolreached = tolreached.max(tolr2);
                // OCCT L806-815.
                let tcurv2 = TopOpeBRepDSCurve::new(c3d, tolreached);
                let icl = dstr.add_curve(tcurv2);
                let mut regfin = super::chfi_ds::ChFiDSRegul::default();
                regfin.set_curve(icl);
                regfin.set_s1(coin.read().expect("surfdata lock").surf(), false);
                regfin.set_s2(sd2.read().expect("surfdata lock").surf(), false);
                fb.my_regul.push(regfin);
                {
                    let mut cornerg = corner.write().expect("stripe lock");
                    cornerg.change_last_curve(icl);
                    cornerg.change_last_parameters(p1fin, p2fin);
                    cornerg.change_index_last_point_on_s1(il1);
                    cornerg.change_index_last_point_on_s2(il2);
                }

                // OCCT L817-827: the S1/S2 shape indices of the coin + the
                // pivot transition + the solid index.
                {
                    let mut coing = coin.write().expect("surfdata lock");
                    coing.change_index_of_s1(dstr.add_shape(&faco));
                    if sameparam {
                        coing.change_index_of_s2(0);
                    } else {
                        coing.change_index_of_s2(dstr.add_shape(&fapiv));
                        coing.change_interference_on_s2().set_transition(pctrans);
                    }
                    let solid_index = {
                        let st1g = st1.read().expect("stripe lock");
                        st1g.solid_index()
                    };
                    corner.write().expect("stripe lock").set_solid_index(solid_index);
                }
                let _ = (&pcurv, &pcurv_l);

                // then the starting Stripe,
                // ------------------------ (OCCT L829-841)
                {
                    let mut st1g = st1.write().expect("stripe lock");
                    st1g.set_curve(icf, isfirst1);
                    st1g.set_index_point(if1, isfirst1, ifa_co1);
                    st1g.set_index_point(if2, isfirst1, ifa_arc1);
                    st1g.set_parameters(isfirst1, p1deb, p2deb);
                    let mut sd1g = sd1.write().expect("surfdata lock");
                    *sd1g.change_vertex(isfirst1, ifa_co1) = pf1.clone();
                    *sd1g.change_vertex(isfirst1, ifa_arc1) = pf2.clone();
                    sd1g.change_interference(ifa_co1).set_parameter(isfirst1, par1);
                    if ifa_co1 == 2 {
                        st1g.set_orientation(Orientation::Reversed, isfirst1);
                    }
                }

                // then the end Stripe,
                // ------------------------- (OCCT L843-855)
                {
                    let mut st2g = st2.write().expect("stripe lock");
                    st2g.set_curve(icl, isfirst2);
                    st2g.set_index_point(il1, isfirst2, ifa_co2);
                    st2g.set_index_point(il2, isfirst2, ifa_arc2);
                    st2g.set_parameters(isfirst2, p1fin, p2fin);
                    let mut sd2g = sd2.write().expect("surfdata lock");
                    *sd2g.change_vertex(isfirst2, ifa_co2) = pl1.clone();
                    *sd2g.change_vertex(isfirst2, ifa_arc2) = pl2.clone();
                    sd2g.change_interference(ifa_co2).set_parameter(isfirst2, par2);
                    if ifa_co2 == 2 {
                        st2g.set_orientation(Orientation::Reversed, isfirst2);
                    }
                }
            }
        } else {
            // it is necessary to make difference with (OCCT L861-1144)
            if !okintercc {
                // OCCT L866.
                panic!("Standard_Failure: TwoCorner : No intersection pc pc");
            }
            // OCCT L868-906: the OnSame/OnDiff assignment bundle.
            let (stsam, stdif, sdsam, sddif, uintpcsam, uintpcdif, ifacosam, ifacodif, ifaopsam, ifaopdif, isfirstsam, isfirstdif);
            if stat1 == ChFiDS_State::OnSame && stat2 == ChFiDS_State::OnDiff {
                stsam = st1.clone();
                sdsam = sd1.clone();
                uintpcsam = u_intpc1;
                ifacosam = ifa_co1;
                ifaopsam = ifa_arc1;
                isfirstsam = isfirst1;
                stdif = st2.clone();
                sddif = sd2.clone();
                uintpcdif = u_intpc2;
                ifacodif = ifa_co2;
                ifaopdif = ifa_arc2;
                isfirstdif = isfirst2;
            } else if stat1 == ChFiDS_State::OnDiff && stat2 == ChFiDS_State::OnSame {
                stsam = st2.clone();
                sdsam = sd2.clone();
                uintpcsam = u_intpc2;
                ifacosam = ifa_co2;
                ifaopsam = ifa_arc2;
                isfirstsam = isfirst2;
                stdif = st1.clone();
                sddif = sd1.clone();
                uintpcdif = u_intpc1;
                ifacodif = ifa_co1;
                ifaopdif = ifa_arc1;
                isfirstdif = isfirst1;
            } else {
                // OCCT L905.
                panic!("Standard_Failure: TwoCorner : Config unknown");
            }
            // It is checked if surface ondiff has a point on arc from the
            // side opposed to the common face and if this arc is connected
            // to the base face opposed to common face of the surface
            // onsame. (OCCT L907-928)
            let cpopdif = {
                let mut sddifg = sddif.write().expect("surfdata lock");
                sddifg.change_vertex(isfirstdif, ifaopdif).clone()
            };
            if !cpopdif.is_on_arc() {
                // OCCT L913.
                panic!("Standard_Failure: TwoCorner : No point on restriction on surface OnDiff");
            }
            let arcopdif = cpopdif.arc().clone();
            let fopsam = dstr
                .shape(sdsam.read().expect("surfdata lock").index_of(ifaopsam))
                .clone();
            {
                let fopsam_edges = face_edges(&brep, &fopsam);
                let n = fopsam_edges.len();
                for (i_ex, cur) in fopsam_edges.iter().enumerate() {
                    if cur.is_same(&arcopdif) {
                        break;
                    } else if i_ex + 1 >= n {
                        // OCCT L924-927 (the !ex.More() branch).
                        panic!("Standard_Failure: TwoCorner : No common face to loop the contour");
                    }
                }
            }
            // OCCT L932-938.
            let (ppopsam, ppcosam) = {
                let sdsamg = sdsam.read().expect("surfdata lock");
                (
                    sdsamg
                        .interference(ifaopsam)
                        .pcurve_on_surf()
                        .expect("PCurveOnSurf")
                        .point_at(uintpcsam),
                    sdsamg
                        .interference(ifacosam)
                        .pcurve_on_surf()
                        .expect("PCurveOnSurf")
                        .point_at(uintpcsam),
                )
            };
            let surfsam: Surface3 = dstr
                .surface(sdsam.read().expect("surfdata lock").surf())
                .surface
                .clone();
            let hsurfsam = BRepAdaptorSurface::initialize_surface(surfsam.clone());
            let mut pcsurfsam: Option<Curve2d> = None;
            let bsam = chfi3d_mkbound_adaptor_two_points(&hsurfsam, ppopsam, ppcosam, fb.tolapp3d, 2.0e-4);
            // OCCT L939-945.
            let (upcopdif, ppopdif, ppcodif) = {
                let sddifg = sddif.read().expect("surfdata lock");
                let upcopdif = sddifg.interference(ifaopdif).parameter(isfirstdif);
                let ppopdif = sddifg
                    .interference(ifaopdif)
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(upcopdif);
                let ppcodif = sddifg
                    .interference(ifacodif)
                    .pcurve_on_surf()
                    .expect("PCurveOnSurf")
                    .point_at(uintpcdif);
                (upcopdif, ppopdif, ppcodif)
            };
            let surfdif: Surface3 = dstr
                .surface(sddif.read().expect("surfdata lock").surf())
                .surface
                .clone();
            let hsurfdif = BRepAdaptorSurface::initialize_surface(surfdif.clone());
            let mut pcsurfdif: Option<Curve2d> = None;
            let bdif = chfi3d_mkbound_adaptor_two_points(&hsurfdif, ppcodif, ppopdif, fb.tolapp3d, 2.0e-4);
            // OCCT L946-952.
            let ppfacsam = {
                let sdsamg = sdsam.read().expect("surfdata lock");
                sdsamg
                    .interference(ifaopsam)
                    .pcurve_on_face()
                    .expect("PCurveOnFace")
                    .point_at(uintpcsam)
            };
            let vvfacsam = {
                let line_index = {
                    let sdsamg = sdsam.read().expect("surfdata lock");
                    sdsamg.interference(ifaopsam).line_index()
                };
                let curv = dstr
                    .curve(line_index)
                    .curve
                    .clone()
                    .expect("DS curve null");
                curv.derivative_at(uintpcsam)
            };
            // OCCT L953-954: PCArcFac.D0(cpopdif.ParameterOnArc(), ppfacdif).
            let mut pcarcfac = BRepAdaptorCurve2d::new();
            pcarcfac.initialize(&brep, &arcopdif, &fopsam);
            let mut ppfacdif = pcarcfac.value(cpopdif.parameter_on_arc());
            // jgv for OCC26173 (OCCT L955-972)
            let surfopsam = BRepAdaptorSurface::initialize(&brep, &fopsam);
            if brep_adaptor_surface_is_u_closed(&surfopsam) {
                let u_period = surfopsam.ulast - surfopsam.ufirst;
                if (ppfacsam.x - ppfacdif.x).abs() > u_period / 2.0 {
                    if ppfacdif.x < ppfacsam.x {
                        ppfacdif.x += u_period;
                    } else {
                        ppfacdif.x -= u_period;
                    }
                }
            }
            // OCCT L973-974: CArcFac.D1(cpopdif.ParameterOnArc(), PPfacdif,
            // VVfacdif).
            let vvfacdif = {
                let (c3d, _, _) =
                    crate::brep_algo::tool::brep_tool_curve(&arcopdif).expect("edge 3d curve");
                c3d.derivative_at(cpopdif.parameter_on_arc())
            };
            // OCCT L975-980.
            let hbrfopsam = BRepAdaptorSurface::initialize(&brep, &fopsam);
            let pcfopsam =
                chfi3d_build_pc_3d(&hbrfopsam, ppfacsam, vvfacsam, ppfacdif, vvfacdif, true);
            let bfac = chfi3d_mkbound_pcurve(&hbrfopsam, Some(&pcfopsam), fb.tolapp3d, 2.0e-4);
            // OCCT L981-983.
            let mut fil = GeomFillConstrainedFilling::new(8, 20);
            fil.init3(&bsam, &bdif, &bfac, true);
            let surfcoin = fil.surface();
            // OCCT L984-997.
            let osurfsam = {
                let sdsamg = sdsam.read().expect("surfdata lock");
                sdsamg.orientation()
            };
            done = match surfcoin {
                Some(surfcoin) => fb.complete_data_surfcoin(
                    &mut coin.write().expect("surfdata lock"),
                    &surfcoin,
                    &hsurfsam,
                    pcsurfsam.as_ref(),
                    &hbrfopsam,
                    None,
                    osurfsam,
                    true,
                    false,
                    false,
                    false,
                    false,
                ),
                None => {
                    // GAP: GeomFill_ConstrainedFilling::Surface pending.
                    false
                }
            };
            if !done {
                // OCCT L1003.
                panic!("Standard_Failure: concavites inverted : fail");
            }
            // Update 3 CornerData and the DS
            // ---------------------------------------- (OCCT L1008-1143)
            let (pf1, pl1): (ChFiDS_CommonPoint, ChFiDS_CommonPoint);
            let pf2: ChFiDS_CommonPoint;
            let pl2: ChFiDS_CommonPoint;
            {
                let mut coing = coin.write().expect("surfdata lock");
                // OCCT L1014-1018: Pf2 = Pl2 = cpopdif.
                pf1 = coing.vertex_first_on_s1().clone();
                pl1 = coing.vertex_last_on_s1().clone();
                *coing.change_vertex_first_on_s2() = cpopdif.clone();
                *coing.change_vertex_last_on_s2() = cpopdif.clone();
                pf2 = cpopdif.clone();
                pl2 = cpopdif.clone();
            }

            // the corner to start,
            // ----------------------- (OCCT L1010-1024)
            let if1 = chfi3d_index_point_in_ds(&pf1, &mut dstr);
            let if2 = chfi3d_index_point_in_ds(&pf2, &mut dstr);
            let il1 = chfi3d_index_point_in_ds(&pl1, &mut dstr);
            let il2 = if2;

            let coin_surf: Surface3 = dstr
                .surface(coin.read().expect("surfdata lock").surf())
                .surface
                .clone();
            // OCCT L1026-1045.
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
            let (c3d, _pcurv, p1deb, p2deb, mut tolreached) = chfi3d_compute_arete(
                &brep, &pf1, pp1, &pf2, pp2, &coin_surf, fb.tolapp3d, fb.tol2d, 0,
            );
            // OCCT L1046-1049: ChFi3d_SameParameter(HC3d, pcFopsam,
            // HBRFopsam, tolapp3d, tolr1).
            let mut tolr1 = 0.0f64;
            let mut pcfopsam = pcfopsam.clone();
            if let Some(c3d) = &c3d {
                chfi3d_same_parameter(c3d, &mut pcfopsam, &hbrfopsam.surface, fb.tolapp3d, &mut tolr1);
            }
            tolreached = tolreached.max(tolr1);
            // OCCT L1050-1051.
            let tcurv1 = TopOpeBRepDSCurve::new(c3d, tolreached);
            let icf = dstr.add_curve(tcurv1);
            // place the pcurve on face in the DS (OCCT L1052-1061)
            let (opc_fopsam, ifopsam) = {
                let sdsamg = sdsam.read().expect("surfdata lock");
                let mut o = sdsamg.interference(ifaopsam).transition();
                let i = sdsamg.index_of(ifaopsam);
                if isfirstsam {
                    o = topabs_reverse(o);
                }
                (o, i)
            };
            let interf = super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(
                icf,
                ifopsam,
                Some(pcfopsam.clone()),
                opc_fopsam,
            );
            dstr.change_shape_interferences(ifopsam).push(interf);

            // OCCT L1063-1070.
            let mut regdeb = super::chfi_ds::ChFiDSRegul::default();
            regdeb.set_curve(icf);
            regdeb.set_s1(coin.read().expect("surfdata lock").surf(), false);
            regdeb.set_s2(ifopsam, true);
            fb.my_regul.push(regdeb);
            {
                let mut cornerg = corner.write().expect("stripe lock");
                cornerg.change_first_curve(icf);
                cornerg.change_first_parameters(p1deb, p2deb);
                cornerg.change_index_first_point_on_s1(if1);
                cornerg.change_index_first_point_on_s2(if2);
            }

            // OCCT L1072-1088.
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
            let (c3d, _pcurv, p1fin, p2fin, mut tolreached) = chfi3d_compute_arete(
                &brep, &pl1, pp1, &pl2, pp2, &coin_surf, fb.tolapp3d, fb.tol2d, 0,
            );
            // OCCT L1089-1092.
            let mut tolr2 = 0.0f64;
            if let (Some(c3d), Some(pcs)) = (&c3d, pcsurfdif.as_mut()) {
                chfi3d_same_parameter(c3d, pcs, &hsurfdif.surface, fb.tolapp3d, &mut tolr2);
            }
            tolreached = tolreached.max(tolr2);
            // OCCT L1093-1098.
            let tcurv2 = TopOpeBRepDSCurve::new(c3d, tolreached);
            let icl = dstr.add_curve(tcurv2);
            let mut regfin = super::chfi_ds::ChFiDSRegul::default();
            regfin.set_curve(icl);
            regfin.set_s1(coin.read().expect("surfdata lock").surf(), false);
            regfin.set_s2(sddif.read().expect("surfdata lock").surf(), false);
            fb.my_regul.push(regfin);
            {
                let mut cornerg = corner.write().expect("stripe lock");
                cornerg.change_last_curve(icl);
                cornerg.change_last_parameters(p1fin, p2fin);
                cornerg.change_index_last_point_on_s1(il1);
                cornerg.change_index_last_point_on_s2(il2);
            }

            // OCCT L1104-1105.
            {
                let mut coing = coin.write().expect("surfdata lock");
                coing.change_index_of_s1(-(sdsam.read().expect("surfdata lock").surf()));
                coing.change_index_of_s2(0);
            }

            // OCCT L1107: corner->SetSolidIndex(stsam->SolidIndex()).
            {
                let solid_index = stsam.read().expect("stripe lock").solid_index();
                corner.write().expect("stripe lock").set_solid_index(solid_index);
            }

            // then Stripe OnSame
            // --------------------- (OCCT L1109-1125)
            let (intcoin1_first, intcoin1_last, intcoin1_line, intcoin1_pc) = {
                let coing = coin.read().expect("surfdata lock");
                let intcoin1 = coing.interference_on_s1();
                (
                    intcoin1.parameter_first(),
                    intcoin1.parameter_last(),
                    intcoin1.line_index(),
                    intcoin1.pcurve_on_face().cloned(),
                )
            };
            {
                let mut stsamg = stsam.write().expect("stripe lock");
                stsamg.set_curve(intcoin1_line, isfirstsam);
                stsamg.in_ds(isfirstsam, 0); // filDS already works from the corner.
                if let Some(pc) = intcoin1_pc {
                    stsamg.change_pcurve(isfirstsam, pc);
                }
                stsamg.set_index_point(if1, isfirstsam, ifaopsam);
                stsamg.set_index_point(il1, isfirstsam, ifacosam);
                stsamg.set_parameters(isfirstsam, intcoin1_first, intcoin1_last);
                let mut sdsamg = sdsam.write().expect("surfdata lock");
                *sdsamg.change_vertex(isfirstsam, ifaopsam) = pf1.clone();
                *sdsamg.change_vertex(isfirstsam, ifacosam) = pl1.clone();
                sdsamg
                    .change_interference_on_s1()
                    .set_parameter(isfirstsam, uintpcsam);
                sdsamg
                    .change_interference_on_s2()
                    .set_parameter(isfirstsam, uintpcsam);
                if ifaopsam == 2 {
                    stsamg.set_orientation(Orientation::Reversed, isfirstsam);
                }
            }

            // then Stripe OnDiff
            // --------------------- (OCCT L1127-1140)
            {
                let mut stdifg = stdif.write().expect("stripe lock");
                stdifg.set_curve(icl, isfirstdif);
                if let Some(pcs) = pcsurfdif.clone() {
                    stdifg.change_pcurve(isfirstdif, pcs);
                }
                stdifg.set_index_point(il2, isfirstdif, ifaopdif);
                stdifg.set_index_point(il1, isfirstdif, ifacodif);
                stdifg.set_parameters(isfirstdif, p1fin, p2fin);
                let mut sddifg = sddif.write().expect("surfdata lock");
                *sddifg.change_vertex(isfirstdif, ifaopdif) = pl2.clone();
                *sddifg.change_vertex(isfirstdif, ifacodif) = pl1.clone();
                sddifg
                    .change_interference(ifacodif)
                    .set_parameter(isfirstdif, uintpcdif);
                if ifaopdif == 1 {
                    stdifg.set_orientation(Orientation::Reversed, isfirstdif);
                }
            }
        }

        // OCCT L1145-1152: myEVIMap + myListStripe.Append(corner).
        if !fb.my_evi_map.contains_key(&vtx.ptr_id()) {
            fb.my_evi_map.insert(vtx.ptr_id(), Vec::new());
        }
        let coin_surf = coin.read().expect("surfdata lock").surf();
        fb.my_evi_map
            .get_mut(&vtx.ptr_id())
            .expect("EVIMap")
            .push(coin_surf);
        fb.my_list_stripe.push(corner.clone());
    }

    fb.my_ds = Some(dstr);
}
