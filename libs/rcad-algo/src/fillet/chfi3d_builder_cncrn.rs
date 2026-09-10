//! OCCT ChFi3d_Builder_CnCrn.cxx — 1:1 translation (Stage 1f).
//!
//! Source: src/ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx
//! (3,927 lines).  This file holds:
//!   - the file-scope static helpers L128-1230,
//!   - the Builder_0.cxx / ChFi3d.cxx helpers that CnCrn depends on and that
//!     have no rcad home yet (cherche_face1 / cherche_element / cherche_edge
//!     / edge_common_faces / AngleEdge / cherche_vertex / SameSide /
//!     SearchFD / IsInFront / IntTraces) — pending relocation into
//!     `chfi3d_builder_0.rs` when that stage covers Builder_0.cxx
//!     L900-1232 / L4549-4663 / L5334-5760 and ChFi3d.cxx L622-645,
//!   - OCCT-form carriers for the TKGeomAlgo / TKTopAlgo classes the corner
//!     pipeline calls whose rcad translations have not landed yet
//!     (FairCurve_Batten, BRepAlgo_NormalProjection, Extrema_ExtCC/ExtPC,
//!     Geom2dInt_GInter, GeomLib::BuildCurve3d, GeomPlate_MakeApprox,
//!     GeomPlate_PlateG0Criterion, BndLib_Add2dCurve, BRepAdaptor_Curve,
//!     Adaptor3d_CurveOnSurface, the GeomPlate_BuildPlateSurface curve
//!     path).
//!
//! `PerformMoreThreeCorner` (L1237-3927) lives in
//! `chfi3d_builder_cncrn_b.rs` (file-size rule).

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema::POnCurve;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_cc::ExtremaExtCC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::Curve2d;
use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::BRepTool as _;
use rcad_kernel::topods::{self, GeomAbsShape, Orientation, Shape, ShapeType};

use crate::brep_algo::normal_projection::BRepAlgoNormalProjection;
use super::chfi3d::topabs_reverse;
use super::chfi3d::chfi3d_index_of_surf_data;
use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_bound_surf, chfi3d_compute_curves, topexp_face_edges,
    topexp_vertices, vec_angle, GeomAdaptorSurface, GeomAbsSurfaceType, P_CONFUSION,
};
use super::chfi3d_ds::TopOpeBRepDSHDataStructure;
use super::chfi_ds::{ChFiDS_CommonPoint, ChFiDSMap, ChFiDSStripeMap, ChFiDSSurfData, SharedStripe};

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L126-146 — Indices.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L128-146
pub fn indices(n: i32, ic: i32, icplus: &mut i32, icmoins: &mut i32) {
    if ic == (n - 1) {
        *icplus = 0;
    } else {
        *icplus = ic + 1;
    }
    if ic == 0 {
        *icmoins = n - 1;
    } else {
        *icmoins = ic - 1;
    }
}

/// OCCT ChFiDS_StripeMap::operator()(V1) — the stripe list of a vertex key
/// (rcad ChFiDSStripeMap is index-addressed; the key scan is the D6
/// architecture shim).
fn stripe_map_find<'a>(map: &'a ChFiDSStripeMap, v1: &Shape) -> &'a [SharedStripe] {
    for (i, k) in map.keys().iter().enumerate() {
        if k.is_same(v1) {
            return map.find_from_index(i + 1);
        }
    }
    &[]
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L150-164 — Calcul_Param.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L150-164
pub fn calcul_param(
    stripe: &super::chfi_ds::ChFiDSStripe,
    jfposit: i32,
    indice: i32,
    isfirst: bool,
    param: &mut f64,
) {
    if jfposit == 2 {
        *param = stripe.set_of_surf_data()[(indice - 1) as usize]
            .read()
            .expect("surfdata lock")
            .interference_on_s2()
            .parameter(isfirst);
    } else {
        *param = stripe.set_of_surf_data()[(indice - 1) as usize]
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .parameter(isfirst);
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L168-183 — Calcul_P2dOnSurf.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L168-183
pub fn calcul_p2d_on_surf(
    stripe: &super::chfi_ds::ChFiDSStripe,
    jfposit: i32,
    indice: i32,
    param: f64,
    p2: &mut DVec2,
) {
    if jfposit == 1 {
        if let Some(pc) = stripe.set_of_surf_data()[(indice - 1) as usize]
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .pcurve_on_surf()
            .cloned()
        {
            *p2 = pc.point_at(param);
        }
    } else if let Some(pc) = stripe.set_of_surf_data()[(indice - 1) as usize]
        .read()
        .expect("surfdata lock")
        .interference_on_s2()
        .pcurve_on_surf()
        .cloned()
    {
        *p2 = pc.point_at(param);
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L187-201 — Calcul_C2dOnFace.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L187-201
pub fn calcul_c2d_on_face(
    stripe: &super::chfi_ds::ChFiDSStripe,
    jfposit: i32,
    indice: i32,
    c2d: &mut Option<rcad_kernel::geom::Curve2d>,
) {
    let fd = stripe.set_of_surf_data()[(indice - 1) as usize]
        .read()
        .expect("surfdata lock");
    if jfposit == 1 {
        *c2d = fd.interference_on_s1().pcurve_on_face().cloned();
    } else {
        *c2d = fd.interference_on_s2().pcurve_on_face().cloned();
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L205-218 — Calcul_Orientation.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L205-218
pub fn calcul_orientation(
    stripe: &super::chfi_ds::ChFiDSStripe,
    jfposit: i32,
    indice: i32,
    orient: &mut Orientation,
) {
    let fd = stripe.set_of_surf_data()[(indice - 1) as usize]
        .read()
        .expect("surfdata lock");
    if jfposit == 1 {
        *orient = fd.interference_on_s1().transition();
    } else {
        *orient = fd.interference_on_s2().transition();
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L222-238 — RemoveSD.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L222-238
pub fn remove_sd(stripe: &mut super::chfi_ds::ChFiDSStripe, num1: i32, num2: i32) {
    let seq = stripe.change_set_of_surf_data();
    if seq.is_empty() {
        return;
    }
    if num1 == num2 {
        // OCCT NCollection_Sequence::Remove(index) — 1-based.
        seq.remove((num1 - 1) as usize);
    } else {
        // OCCT NCollection_Sequence::Remove(from, to) — inclusive, 1-based.
        seq.drain(((num1 - 1) as usize)..(num2 as usize));
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L245-274 — cherche_edge1: find common edge
// of faces F1 and F2.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L245-274
pub fn cherche_edge1(brep: &topods::BRep, f1: &Shape, f2: &Shape, edge: &mut Shape) {
    let mut trouve = false;
    let map_e1 = topexp_face_edges(brep, f1);
    let map_e2 = topexp_face_edges(brep, f2);
    for ecur1 in &map_e1 {
        if trouve {
            break;
        }
        for ecur2 in &map_e2 {
            if trouve {
                break;
            }
            if ecur2.is_same(ecur1) {
                *edge = ecur1.clone();
                trouve = true;
            }
        }
    }
    if edge.is_null() {
        panic!("Standard_ConstructionError: Failed to find edge");
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5334-5354 — ChFi3d_cherche_face1 (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L5334-5354
pub fn chfi3d_cherche_face1(map: &[Shape], f1: &Shape, f: &mut Shape) {
    let mut trouve = false;
    for fcur in map {
        if trouve {
            break;
        }
        let fcur = fcur.clone();
        if !fcur.is_same(f1) {
            *f = fcur;
            trouve = true;
        }
    }
    if f.is_null() {
        panic!("Standard_ConstructionError: Failed to find face");
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5361-5403 — ChFi3d_cherche_element (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L5361-5403
pub fn chfi3d_cherche_element(
    brep: &topods::BRep,
    v: &Shape,
    e1: &Shape,
    f1: &Shape,
    e: &mut Shape,
    vtx: &mut Shape,
) {
    let mut trouve = false;
    let map_e = topexp_face_edges(brep, f1);
    for ecur in &map_e {
        if trouve {
            break;
        }
        if !ecur.is_same(e1) {
            let (v1, v2) = topexp_vertices(ecur);
            if !v1.is_null() && !v2.is_null() {
                if v1.is_same(v) {
                    *vtx = v2;
                    *e = ecur.clone();
                    trouve = true;
                } else if v2.is_same(v) {
                    *vtx = v1;
                    *e = ecur.clone();
                    trouve = true;
                }
            }
        }
    }
    if e.is_null() {
        panic!("Standard_ConstructionError: Failed to find element");
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5410-5461 — ChFi3d_cherche_edge (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L5410-5461
pub fn chfi3d_cherche_edge(
    brep: &topods::BRep,
    v: &Shape,
    e1: &[Shape],
    f1: &Shape,
    e: &mut Shape,
    vtx: &mut Shape,
) {
    let mut trouve = false;
    let map_e = topexp_face_edges(brep, f1);
    for ecur in &map_e {
        if trouve {
            break;
        }
        let mut same = false;
        for ei in e1 {
            if ecur.is_same(ei) {
                same = true;
            }
        }
        if !same {
            let (v1, v2) = topexp_vertices(ecur);
            if !v1.is_null() && !v2.is_null() {
                if v1.is_same(v) {
                    *vtx = v2;
                    *e = ecur.clone();
                    trouve = true;
                } else if v2.is_same(v) {
                    *vtx = v1;
                    *e = ecur.clone();
                    trouve = true;
                }
            }
        }
    }
    if e.is_null() {
        panic!("Standard_ConstructionError: Failed to find edge");
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5498-5521 — ChFi3d_edge_common_faces (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L5498-5521
pub fn chfi3d_edge_common_faces(map_ef: &[Shape], f1: &mut Shape, f2: &mut Shape) {
    if map_ef.is_empty() {
        return;
    }
    *f1 = map_ef[0].clone();
    let mut trouve = false;
    for f in map_ef {
        if trouve {
            break;
        }
        if !f.is_same(f1) {
            *f2 = f.clone();
            trouve = true;
        }
    }
    if !trouve {
        *f2 = f1.clone();
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5527-5549 — ChFi3d_AngleEdge (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L5527-5549
pub fn chfi3d_angle_edge(brep: &topods::BRep, vtx: &Shape, e1: &Shape, e2: &Shape) -> f64 {
    let bc1 = BRepAdaptorCurve::initialize(brep, e1);
    let bc2 = BRepAdaptorCurve::initialize(brep, e2);
    let par_e1 = brep_tool_parameter(brep, vtx, e1);
    let par_e2 = brep_tool_parameter(brep, vtx, e2);
    let (_p1, mut dir1) = bc1.d1(par_e1);
    let (_p2, mut dir2) = bc2.d1(par_e2);
    let first1 = brep.first_vertex(e1);
    let first2 = brep.first_vertex(e2);
    if !vtx.is_same(&first1) {
        dir1 = -dir1;
    }
    if !vtx.is_same(&first2) {
        dir2 = -dir2;
    }
    vec_angle(dir1, dir2).abs()
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5716-5748 — ChFi3d_cherche_vertex (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L5716-5748
pub fn chfi3d_cherche_vertex(e1: &Shape, e2: &Shape, vertex: &mut Shape, trouve: &mut bool) {
    *trouve = false;
    let (vf1, vl1) = topexp_vertices(e1);
    let (vf2, vl2) = topexp_vertices(e2);
    for vcur1 in [&vf1, &vl1] {
        if *trouve {
            break;
        }
        for vcur2 in [&vf2, &vl2] {
            if *trouve {
                break;
            }
            if vcur2.is_same(vcur1) {
                *vertex = vcur1.clone();
                *trouve = true;
            }
        }
    }
}

// =========================================================================
// OCCT TopAbs::Compose(A, B) (TopAbs_Orientation.hxx L69-95).
// =========================================================================

// OCCT TopAbs_Orientation.hxx L69-95
pub fn topabs_compose(a: Orientation, b: Orientation) -> Orientation {
    if a == Orientation::Forward {
        b
    } else if a == Orientation::Reversed {
        topabs_reverse(b)
    } else {
        a
    }
}

// =========================================================================
// OCCT ChFi3d.cxx L622-645 — ChFi3d::SameSide (pending relocation).
// =========================================================================

// OCCT ChFi3d.cxx L622-645
pub fn chfi3d_same_side(
    or: Orientation,
    or_save1: Orientation,
    or_save2: Orientation,
    or_face1: Orientation,
    or_face2: Orientation,
) -> bool {
    let o1 = if or == or_face1 {
        or_save1
    } else {
        topabs_reverse(or_save1)
    };
    let o2 = if or == or_face2 {
        or_save2
    } else {
        topabs_reverse(or_save2)
    };
    o1 == o2
}

/// OCCT `TopExp::MapShapes(S, TopAbs_EDGE, M)` (TopExp.cxx L75-88) — fills an
/// `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>` with the
/// unique edges of `S` in exploration order.  rcad encoding: a `Vec<Shape>`
/// deduplicated by the ShapeMapHasher key (TShape pointer, Location).
pub(crate) fn top_exp_map_shapes_edges(the_s: &Shape) -> Vec<Shape> {
    let mut a_map: Vec<Shape> = Vec::new();
    for an_edge in crate::brep_algo::tool::explorer(the_s, ShapeType::Edge, ShapeType::Shape) {
        let a_key = (an_edge.ptr_id(), an_edge.location);
        if !a_map
            .iter()
            .any(|k| (k.ptr_id(), k.location) == a_key)
        {
            a_map.push(an_edge);
        }
    }
    a_map
}

// =========================================================================
// OCCT-form carriers for the not-yet-translated TKGeomAlgo / TKTopAlgo
// dependencies of ChFi3d_Builder_CnCrn.cxx.  Each carrier keeps the OCCT
// call surface; bodies marked "pending" hold deterministic neutral values
// until the dedicated module translations land.
// =========================================================================

/// OCCT BRepAdaptor_Curve — the edge 3D curve adaptor.
/// Pending: full BRepAdaptor_Curve translation; served over the BRep edge
/// curve + CurveEval (precedent: fillet/chfi2d_ana_fillet_algo.rs).
pub struct BRepAdaptorCurve {
    curve: Option<rcad_kernel::geom::Curve3>,
    first: f64,
    last: f64,
}

impl BRepAdaptorCurve {
    pub fn initialize(brep: &topods::BRep, e: &Shape) -> Self {
        match brep.edge_curve_world(e) {
            Some((c, r)) => BRepAdaptorCurve {
                curve: Some(c),
                first: r[0],
                last: r[1],
            },
            None => BRepAdaptorCurve {
                curve: None,
                first: 0.0,
                last: 0.0,
            },
        }
    }

    /// OCCT BRepAdaptor_Curve::Value(U).
    pub fn value(&self, u: f64) -> DVec3 {
        match &self.curve {
            Some(c) => c.point_at(u),
            None => DVec3::ZERO,
        }
    }

    /// OCCT BRepAdaptor_Curve::D1(U, P, V).
    pub fn d1(&self, u: f64) -> (DVec3, DVec3) {
        match &self.curve {
            Some(c) => (c.point_at(u), c.derivative_at(u)),
            None => (DVec3::ZERO, DVec3::ZERO),
        }
    }

    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// OCCT BRepAdaptor_Curve::Resolution(R3d) — pending; neutral 0 keeps
    /// the CnCrn recoil offset at 0 until the curve-resolution machinery
    /// lands.
    pub fn resolution(&self, _r3d: f64) -> f64 {
        0.0
    }
}

/// OCCT Geom2dAdaptor_Curve — the pcurve adaptor (Load with optional trim).
pub struct Geom2dAdaptorCurve {
    curve: Option<rcad_kernel::geom::Curve2d>,
    first: f64,
    last: f64,
}

impl Geom2dAdaptorCurve {
    /// OCCT Geom2dAdaptor_Curve::Load(C, U1, U2).
    pub fn load(c: rcad_kernel::geom::Curve2d, u1: f64, u2: f64) -> Self {
        Geom2dAdaptorCurve {
            curve: Some(c),
            first: u1,
            last: u2,
        }
    }

    /// OCCT Geom2dAdaptor_Curve::Value(U) / D0.
    pub fn value(&self, u: f64) -> DVec2 {
        match &self.curve {
            Some(c) => c.point_at(u),
            None => DVec2::ZERO,
        }
    }

    /// OCCT Geom2dAdaptor_Curve::D1(U, P, V).
    pub fn d1(&self, u: f64) -> (DVec2, DVec2) {
        match &self.curve {
            Some(c) => (c.point_at(u), c.derivative_at(u)),
            None => (DVec2::ZERO, DVec2::ZERO),
        }
    }

    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    pub fn last_parameter(&self) -> f64 {
        self.last
    }
}

/// OCCT Adaptor3d_CurveOnSurface — the pcurve-surfaced 3D curve
/// (pending dedicated translation; served over the pcurve adaptor + the
/// GeomAdaptor surface evaluation).
pub struct Adaptor3dCurveOnSurface {
    pcurve: Geom2dAdaptorCurve,
    surf: GeomAdaptorSurface,
}

impl Adaptor3dCurveOnSurface {
    pub fn new(pcurve: Geom2dAdaptorCurve, surf: GeomAdaptorSurface) -> Self {
        Adaptor3dCurveOnSurface { pcurve, surf }
    }

    /// OCCT Adaptor3d_CurveOnSurface::FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.pcurve.first_parameter()
    }

    /// OCCT Adaptor3d_CurveOnSurface::LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.pcurve.last_parameter()
    }

    /// OCCT Adaptor3d_CurveOnSurface::Value(U).
    pub fn value(&self, u: f64) -> DVec3 {
        let uv = self.pcurve.value(u);
        self.surf.value(uv.x, uv.y)
    }
}

/// OCCT FairCurve_AnalysisCode — pending translation (FairCurve module).
#[derive(Debug, Clone, Copy)]
pub struct FairCurveAnalysisCode;

/// OCCT FairCurve_Batten — pending TKFairCurve translation.  The carrier
/// keeps the OCCT call surface; Compute() reports failure so the OCCT
/// fallback in CalculBatten (CalculDroite) executes, exactly as the OCCT
/// code does when the batten solver fails.
pub struct FairCurveBatten {
    #[allow(dead_code)]
    point1: DVec2,
    #[allow(dead_code)]
    point2: DVec2,
    #[allow(dead_code)]
    height: f64,
    #[allow(dead_code)]
    free_sliding: bool,
    #[allow(dead_code)]
    angle1: Option<f64>,
    #[allow(dead_code)]
    angle2: Option<f64>,
    #[allow(dead_code)]
    constraint_order1: Option<i32>,
    #[allow(dead_code)]
    constraint_order2: Option<i32>,
}

impl FairCurveBatten {
    /// OCCT FairCurve_Batten::FairCurve_Batten(P1, P2, Height).
    pub fn new(p1: DVec2, p2: DVec2, height: f64) -> Self {
        FairCurveBatten {
            point1: p1,
            point2: p2,
            height,
            free_sliding: false,
            angle1: None,
            angle2: None,
            constraint_order1: None,
            constraint_order2: None,
        }
    }

    /// OCCT FairCurve_Batten::SetFreeSliding.
    pub fn set_free_sliding(&mut self, free_sliding: bool) {
        self.free_sliding = free_sliding;
    }

    /// OCCT FairCurve_Batten::SetAngle1.
    pub fn set_angle1(&mut self, angle: f64) {
        self.angle1 = Some(angle);
    }

    /// OCCT FairCurve_Batten::SetAngle2.
    pub fn set_angle2(&mut self, angle: f64) {
        self.angle2 = Some(angle);
    }

    /// OCCT FairCurve_Batten::SetConstraintOrder1.
    pub fn set_constraint_order1(&mut self, order: i32) {
        self.constraint_order1 = Some(order);
    }

    /// OCCT FairCurve_Batten::SetConstraintOrder2.
    pub fn set_constraint_order2(&mut self, order: i32) {
        self.constraint_order2 = Some(order);
    }

    /// OCCT FairCurve_Batten::Compute(AnalysisCode, NbIterations,
    /// Tolerance) — pending TKFairCurve translation; reports not-done so
    /// the OCCT CalculBatten fallback (CalculDroite) runs.
    pub fn compute(&mut self, _iana: &mut FairCurveAnalysisCode, _iters: i32, _tol: f64) -> bool {
        false
    }

    /// OCCT FairCurve_Batten::Curve().
    pub fn curve(&self) -> Option<rcad_kernel::geom::Curve2d> {
        None
    }

    /// OCCT FairCurve_Batten::Dump — debug only, pending.
    pub fn dump(&self) {}
}

/// OCCT Bnd_Box2d — the 2D bounding box (minimal OCCT call surface).
#[derive(Debug, Clone, Copy)]
pub struct BndBox2d {
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
}

impl BndBox2d {
    pub fn new() -> Self {
        BndBox2d {
            xmin: f64::MAX,
            ymin: f64::MAX,
            xmax: f64::MIN,
            ymax: f64::MIN,
        }
    }

    /// OCCT Bnd_Box2d::Update(x, y, X, Y).
    pub fn update(&mut self, x: f64, y: f64, big_x: f64, big_y: f64) {
        self.xmin = self.xmin.min(x);
        self.ymin = self.ymin.min(y);
        self.xmax = self.xmax.max(big_x);
        self.ymax = self.ymax.max(big_y);
    }

    /// OCCT Bnd_Box2d::Get(x, y, X, Y).
    pub fn get(&self) -> (f64, f64, f64, f64) {
        (self.xmin, self.ymin, self.xmax, self.ymax)
    }
}

/// OCCT BndLib_Add2dCurve::Add(C, Tol, B) — pending TKTopAlgo/BndLib 2D
/// translation; served by dense sampling of the pcurve adaptor.
pub fn bnd_lib_add2d_curve(acur: &Geom2dAdaptorCurve, _tol: f64, b: &mut BndBox2d) {
    let (f, l) = (acur.first_parameter(), acur.last_parameter());
    let n = 32;
    for k in 0..=n {
        let t = f + (l - f) * (k as f64) / (n as f64);
        let p = acur.value(t);
        b.update(p.x, p.y, p.x, p.y);
    }
}

/// OCCT GeomLProp_CLProps2d — the point/tangent evaluator on a 2D curve.
pub struct GeomLPropCLProps2d {
    #[allow(dead_code)]
    p: DVec2,
    tangent: DVec2,
}

impl GeomLPropCLProps2d {
    /// OCCT GeomLProp_CLProps2d(C, U, N, Tolerance).
    pub fn new(c: &rcad_kernel::geom::Curve2d, u: f64, _n: i32, _tolerance: f64) -> Self {
        GeomLPropCLProps2d {
            p: c.point_at(u),
            tangent: c.tangent_at(u),
        }
    }

    /// OCCT GeomLProp_CLProps2d::Tangent(D).
    pub fn tangent(&self) -> DVec2 {
        self.tangent
    }
}

/// OCCT Extrema_POnCurv — the parameter/point pair.
#[derive(Debug, Clone, Copy, Default)]
pub struct ExtremaPOnCurv {
    parameter: f64,
    value: DVec3,
}

impl ExtremaPOnCurv {
    pub fn new(parameter: f64, value: DVec3) -> Self {
        ExtremaPOnCurv { parameter, value }
    }
    pub fn parameter(&self) -> f64 {
        self.parameter
    }
    pub fn value(&self) -> DVec3 {
        self.value
    }
}

/// OCCT Extrema_ExtPC — pending TKGeomAlgo translation; IsDone()=false
/// keeps the OCCT fallback branches of CurveHermite in effect.
pub struct ExtremaExtPC;

#[allow(dead_code)]
impl ExtremaExtPC {
    pub fn new() -> Self {
        ExtremaExtPC
    }
    pub fn is_done(&self) -> bool {
        false
    }
    pub fn nb_ext(&self) -> i32 {
        0
    }
    pub fn point(&self, _n: i32) -> ExtremaPOnCurv {
        ExtremaPOnCurv::default()
    }
}

/// OCCT Geom2dInt_GInter re-host (architecture difference #31, mirroring
/// brep_offset_inter2d.rs) — the TheIntPCurvePCurveOfGInter vehicle over the
/// translated IntRes2d/IntCurve machinery; the results surface (IsDone /
/// IsEmpty / NbPoints / Point / NbSegments) mirrors the OCCT GInter.  The
/// curve adaptors are the Geom2dAdaptorCurve loads above (curve + trimmed
/// domain).
pub struct Geom2dIntGInter {
    base: crate::geomalgo::int_res2d::IntersectionBase,
}

impl Geom2dIntGInter {
    pub fn new() -> Self {
        Geom2dIntGInter {
            base: crate::geomalgo::int_res2d::IntersectionBase::new(),
        }
    }

    fn engine(
        c1: &Curve2d,
        f1: f64,
        l1: f64,
        c2: &Curve2d,
        f2: f64,
        l2: f64,
        tol_conf: f64,
        tol: f64,
    ) -> crate::geomalgo::int_res2d::IntersectionBase {
        let mut inter = crate::geomalgo::geom2d_int::TheIntPCurvePCurveOfGInter::new();
        let d1 = crate::geomalgo::int_res2d::Domain::bounded(
            c1.point_at(f1),
            f1,
            tol_conf,
            c1.point_at(l1),
            l1,
            tol_conf,
        );
        let d2 = crate::geomalgo::int_res2d::Domain::bounded(
            c2.point_at(f2),
            f2,
            tol_conf,
            c2.point_at(l2),
            l2,
            tol_conf,
        );
        inter.perform(c1, &d1, c2, &d2, tol_conf, tol);
        inter.base
    }

    /// OCCT Geom2dInt_GInter::Perform(C, TolConf, Tol) — the single-curve
    /// (auto-intersection) form.
    pub fn perform(&mut self, c1: &Geom2dAdaptorCurve, tol_conf: f64, tol: f64) {
        let curve1 = match &c1.curve {
            Some(c) => c.clone(),
            None => return,
        };
        self.base = Self::engine(
            &curve1,
            c1.first,
            c1.last,
            &curve1,
            c1.first,
            c1.last,
            tol_conf,
            tol,
        );
    }

    /// OCCT Geom2dInt_GInter::Perform(C1, C2, TolConf, Tol).
    pub fn perform2(
        &mut self,
        c1: &Geom2dAdaptorCurve,
        c2: &Geom2dAdaptorCurve,
        tol_conf: f64,
        tol: f64,
    ) {
        let curve1 = match &c1.curve {
            Some(c) => c.clone(),
            None => return,
        };
        let curve2 = match &c2.curve {
            Some(c) => c.clone(),
            None => return,
        };
        self.base = Self::engine(
            &curve1,
            c1.first,
            c1.last,
            &curve2,
            c2.first,
            c2.last,
            tol_conf,
            tol,
        );
    }
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }
    pub fn is_empty(&self) -> bool {
        self.base.is_empty()
    }
    pub fn nb_segments(&self) -> i32 {
        self.base.nb_segments() as i32
    }
    pub fn nb_points(&self) -> i32 {
        self.base.nb_points() as i32
    }
    /// OCCT IntRes2d_Intersection::Point(i) — (Value, ParamOnFirst,
    /// ParamOnSecond).
    pub fn point(&self, n: i32) -> (DVec2, f64, f64) {
        let p = self.base.point(n as usize);
        (p.value(), p.param_on_first(), p.param_on_second())
    }
}

/// OCCT GeomPlate_PlateG0Criterion — pending TKGeomAlgo translation.
#[derive(Default)]
pub struct GeomPlatePlateG0Criterion;

/// OCCT GeomPlate_Surface — carrier over the surface payload pending the
/// GeomPlate module wiring used here.
pub struct GeomPlateSurfaceCarrier {
    #[allow(dead_code)]
    pub surface: rcad_kernel::geom::Surface3,
}

/// OCCT GeomPlate_MakeApprox — pending TKGeomAlgo translation; Surface()
/// reports null until the MakeApprox chain lands.
#[derive(Default)]
pub struct GeomPlateMakeApprox;

#[allow(dead_code)]
impl GeomPlateMakeApprox {
    /// OCCT GeomPlate_MakeApprox::GeomPlate_MakeApprox(SurfPlate,
    /// Criterion, Tol3d, Nbmax, Degmax).
    pub fn new(
        _plate: &GeomPlateSurfaceCarrier,
        _criterion: &GeomPlatePlateG0Criterion,
        _tol3d: f64,
        _nbmax: i32,
        _degmax: i32,
    ) -> Self {
        GeomPlateMakeApprox
    }
    pub fn surface(&self) -> Option<rcad_kernel::geom::Surface3> {
        None
    }
    pub fn approx_error(&self) -> f64 {
        0.0
    }
    pub fn criterion_error(&self) -> f64 {
        0.0
    }
}

/// OCCT GeomPlate_CurveConstraint — the curve constraint carrier (curve-on-
/// surface + constraint order + point count + tolerances).
pub struct GeomPlateCurveConstraint {
    /// OCCT: Handle(Adaptor3d_CurveOnSurface) myCurve.
    #[allow(dead_code)]
    pub my_curve: Adaptor3dCurveOnSurface,
    /// OCCT: int myOrder.
    #[allow(dead_code)]
    pub my_order: i32,
    /// OCCT: int myNbPoints.
    #[allow(dead_code)]
    pub my_nb_points: i32,
    /// OCCT: double myTolCurve / myTolAng / myTolCurv.
    #[allow(dead_code)]
    pub my_tol_curve: f64,
    #[allow(dead_code)]
    pub my_tol_ang: f64,
    #[allow(dead_code)]
    pub my_tol_curv: f64,
}

impl GeomPlateCurveConstraint {
    /// OCCT GeomPlate_CurveConstraint::CurveConstraint(ConstraintCurve,
    /// Order, NbPoints, TolDist, TolAng, TolCurv).
    pub fn new(
        constraint_curve: Adaptor3dCurveOnSurface,
        order: i32,
        nb_points: i32,
        tol_dist: f64,
        tol_ang: f64,
        tol_curv: f64,
    ) -> Self {
        GeomPlateCurveConstraint {
            my_curve: constraint_curve,
            my_order: order,
            my_nb_points: nb_points,
            my_tol_curve: tol_dist,
            my_tol_ang: tol_ang,
            my_tol_curv: tol_curv,
        }
    }
}

/// OCCT GeomPlate_BuildPlateSurface — the curve-constraint path used by
/// ChFi3d_Builder_CnCrn.cxx.  The rcad geomplate module covers the
/// anchor-out-of-scope point-constraint path (geomplate/build_plate_surface.rs
/// notes the curve path as not ported); this carrier keeps the OCCT call
/// surface (Add / Perform / IsDone / Surface / Curves2d / Disc2dContour /
/// Disc3dContour / G0Error) until the geomplate curve path lands.
pub struct GeomPlateBuildPlateSurface {
    #[allow(dead_code)]
    degree: i32,
    #[allow(dead_code)]
    nbcurvpnt: i32,
    #[allow(dead_code)]
    nbiter: i32,
    #[allow(dead_code)]
    tol2d: f64,
    #[allow(dead_code)]
    tolapp3d: f64,
    #[allow(dead_code)]
    angular: f64,
    /// OCCT: myLinCont — the added CurveConstraints.
    #[allow(dead_code)]
    pub my_lin_cont: Vec<GeomPlateCurveConstraint>,
    /// OCCT: myCurves2d — the 2d constraint curves after Perform.
    my_curves2d: Vec<rcad_kernel::geom::Curve2d>,
    done: bool,
    surface: Option<rcad_kernel::geom::Surface3>,
    g0_error: f64,
}

impl GeomPlateBuildPlateSurface {
    /// OCCT ctor with degree (GeomPlate_BuildPlateSurface.cxx L186-217;
    /// OCCT defaults TolCurv=0.1, Anisotropie=false).
    pub fn new(
        degree: i32,
        nbcurvpnt: i32,
        nbiter: i32,
        tol2d: f64,
        tol3d: f64,
        tolang: f64,
    ) -> Self {
        GeomPlateBuildPlateSurface {
            degree,
            nbcurvpnt,
            nbiter,
            tol2d,
            tolapp3d: tol3d,
            angular: tolang,
            my_lin_cont: Vec::new(),
            my_curves2d: Vec::new(),
            done: false,
            surface: None,
            g0_error: 0.0,
        }
    }

    /// OCCT GeomPlate_BuildPlateSurface::Add(Cont).
    #[allow(dead_code)]
    pub fn add(&mut self, cont: GeomPlateCurveConstraint) {
        self.my_lin_cont.push(cont);
    }

    /// OCCT GeomPlate_BuildPlateSurface::Perform() — pending the geomplate
    /// curve path; the rcad carrier marks not-done so the OCCT partial-
    /// result branch of PerformMoreThreeCorner (L3893-3926) runs.
    pub fn perform(&mut self) {
        self.done = false;
    }

    /// OCCT GeomPlate_BuildPlateSurface::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT GeomPlate_BuildPlateSurface::Surface().
    pub fn surface(&self) -> Option<GeomPlateSurfaceCarrier> {
        self.surface
            .as_ref()
            .map(|s| GeomPlateSurfaceCarrier { surface: s.clone() })
    }

    /// OCCT GeomPlate_BuildPlateSurface::Curves2d() — the HArray1 of the 2d
    /// constraint curves; OCCT indexing is 1-based.
    #[allow(dead_code)]
    pub fn curves2d_value(&self, i: i32) -> Option<rcad_kernel::geom::Curve2d> {
        self.my_curves2d.get((i - 1) as usize).cloned()
    }

    /// OCCT GeomPlate_BuildPlateSurface::Disc2dContour(NbIsos, Sequence2d)
    /// — pending.
    pub fn disc2d_contour(&self, _nb_isos: i32, seq: &mut Vec<DVec2>) {
        seq.clear();
    }

    /// OCCT GeomPlate_BuildPlateSurface::Disc3dContour(NbIsos, Order,
    /// Sequence3d) — pending.
    pub fn disc3d_contour(&self, _nb_isos: i32, _order: i32, seq: &mut Vec<DVec3>) {
        seq.clear();
    }

    /// OCCT GeomPlate_BuildPlateSurface::G0Error().
    pub fn g0_error(&self) -> f64 {
        self.g0_error
    }
}

/// OCCT GeomLib::BuildCurve3d(Tolerance, CurvOnSurf, First, Last, Curve3d,
/// Maxdev, AvDev) — pending TKTopAlgo translation; the carrier leaves the
/// 3d curve null (OCCT TopOpeBRepDS_Curve carries a nullable curve) with
/// zero deviation until the approximation lands.
#[allow(clippy::too_many_arguments)]
pub fn geom_lib_build_curve3d(
    _tolerance: f64,
    _curv_on_surf: &Adaptor3dCurveOnSurface,
    _first: f64,
    _last: f64,
    curve3d: &mut Option<rcad_kernel::geom::Curve3>,
    maxdev: &mut f64,
    av_dev: &mut f64,
) {
    *curve3d = None;
    *maxdev = 0.0;
    *av_dev = 0.0;
}

/// OCCT PLib::HermiteCoefficients(0, 1, 1, 1, MatCoefs) — the 4x4 power-
/// basis conversion matrix solved by math_Gauss in OCCT
/// (PLib.cxx L1404-1478).  For FirstParameter=0, LastParameter=1 and
/// derivative order 1 at both ends the constraint rows are
/// [F(0), F'(0), F(1), F'(1)] and the exact solution matrix is
/// [[1,0,0,0],[0,1,0,0],[-3,-2,3,-1],[2,1,-2,1]].
fn plib_hermite_coefficients_0_1_c1() -> [[f64; 4]; 4] {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [-3.0, -2.0, 3.0, -1.0],
        [2.0, 1.0, -2.0, 1.0],
    ]
}

/// OCCT PLib::CoefficientsPoles(Coefs, NoWeights, Poles, NoWeights) — the
/// power-basis -> Bezier poles conversion of a degree-3 vector polynomial
/// (PLib.cxx L1482-1493).  poles[0]=a0, poles[1]=a0+a1/3,
/// poles[2]=a0+(2a1)/3+a2/3, poles[3]=a0+a1+a2+a3.
fn plib_coefficients_poles_no_weights(coefs: &[DVec3; 4]) -> [DVec3; 4] {
    let a0 = coefs[0];
    let a1 = coefs[1];
    let a2 = coefs[2];
    let a3 = coefs[3];
    [
        a0,
        a0 + a1 / 3.0,
        a0 + a1 * (2.0 / 3.0) + a2 * (1.0 / 3.0),
        a0 + a1 + a2 + a3,
    ]
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L276-527 — CurveHermite: calculate a curve
// 3d using polynoms of Hermite.  the edge is a regular edge. Curve 3D is
// constructed between edges icmoins and icplus.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L283-527
#[allow(clippy::too_many_arguments)]
pub fn curve_hermite(
    brep: &topods::BRep,
    dstr: &TopOpeBRepDSHDataStructure,
    cd_icmoins: &super::chfi_ds::ChFiDSStripe,
    jficmoins: i32,
    icmoins: i32,
    picmoins: f64,
    sensicmoins: i32,
    sharpicmoins: bool,
    evive_icmoins: &Shape,
    cd_icplus: &super::chfi_ds::ChFiDSStripe,
    jficplus: i32,
    icplus: i32,
    picplus: f64,
    sensicplus: i32,
    sharpicplus: bool,
    evive_icplus: &Shape,
    nbface: i32,
    ecom: &mut Vec<Shape>,
    face: &[Shape],
    proj2d: &mut Vec<Option<rcad_kernel::geom::Curve2d>>,
    cproj: &mut Vec<Option<rcad_kernel::geom::Curve3>>,
    eproj: &mut Vec<Shape>,
    param: &mut Vec<f64>,
    error: &mut f64,
) {
    let mut p01;
    let mut p02;
    let mut d11;
    let mut d12;
    let mut up1 = 0.0f64;
    let mut up2 = 0.0f64;
    let mut ilin: i32;
    // OCCT: c1 = the 3d curve at the icmoins end.
    let c1: Option<rcad_kernel::geom::Curve3>;
    if sharpicmoins {
        // OCCT: c1 = BRep_Tool::Curve(Eviveicmoins, up1, up2);
        match brep.edge_curve_world(evive_icmoins) {
            Some((c, r)) => {
                up1 = r[0];
                up2 = r[1];
                c1 = Some(c);
            }
            None => {
                c1 = None;
            }
        }
    } else {
        if jficmoins == 1 {
            ilin = cd_icmoins.set_of_surf_data()[(icmoins - 1) as usize]
                .read()
                .expect("surfdata lock")
                .interference_on_s1()
                .line_index();
        } else {
            ilin = cd_icmoins.set_of_surf_data()[(icmoins - 1) as usize]
                .read()
                .expect("surfdata lock")
                .interference_on_s2()
                .line_index();
        }
        c1 = dstr.curve(ilin).curve().cloned();
    }
    // OCCT: c2 = the 3d curve at the icplus end.
    let c2: Option<rcad_kernel::geom::Curve3>;
    if sharpicplus {
        // OCCT: c2 = BRep_Tool::Curve(Eviveicplus, up1, up2);
        match brep.edge_curve_world(evive_icplus) {
            Some((c, r)) => {
                up1 = r[0];
                up2 = r[1];
                c2 = Some(c);
            }
            None => {
                c2 = None;
            }
        }
    } else {
        let jfp = 3 - jficplus;
        if jfp == 1 {
            ilin = cd_icplus.set_of_surf_data()[(icplus - 1) as usize]
                .read()
                .expect("surfdata lock")
                .interference_on_s1()
                .line_index();
        } else {
            ilin = cd_icplus.set_of_surf_data()[(icplus - 1) as usize]
                .read()
                .expect("surfdata lock")
                .interference_on_s2()
                .line_index();
        }
        c2 = dstr.curve(ilin).curve().cloned();
    }
    let c1 = match c1 {
        Some(c) => c,
        None => panic!("Standard_ConstructionError: Failed to get 3D curve of edge"),
    };
    let c2 = match c2 {
        Some(c) => c,
        None => panic!("Standard_ConstructionError: Failed to get 3D curve of edge"),
    };
    // OCCT: c1->D1(picmoins, p01, d11); c2->D1(picplus, p02, d12);
    p01 = c1.point_at(picmoins);
    d11 = c1.derivative_at(picmoins);
    p02 = c2.point_at(picplus);
    d12 = c2.derivative_at(picplus);
    let size: usize = 4;
    let mat_coefs = plib_hermite_coefficients_0_1_c1();
    let mut cont: [DVec3; 4] = [DVec3::ZERO; 4];
    let l1 = p01.distance(p02);
    let mut lambda = 1.0 / (d11.length() / l1).max(1.0e-6);
    cont[0] = p01;
    if sensicmoins == 1 {
        cont[1] = d11 * (-lambda);
    } else {
        cont[1] = d11 * lambda;
    }
    lambda = 1.0 / (d12.length() / l1).max(1.0e-6);
    cont[2] = p02;
    if sensicplus == 1 {
        cont[3] = d12 * lambda;
    } else {
        cont[3] = d12 * (-lambda);
    }
    // OCCT: ExtraCoeffs(jj).ChangeCoord() += MatCoefs(ii, jj) * Cont(ii);
    let mut extra_coeffs: [DVec3; 4] = [DVec3::ZERO; 4];
    for ii in 0..size {
        for jj in 0..size {
            extra_coeffs[jj] += mat_coefs[ii][jj] * cont[ii];
        }
    }
    // OCCT: PLib::CoefficientsPoles(ExtraCoeffs, NoWeights, ExtrapPoles,
    //       NoWeights); Bezier = new Geom_BezierCurve(ExtrapPoles);
    let extrap_poles = plib_coefficients_poles_no_weights(&extra_coeffs);
    let bezier = rcad_kernel::geom::Curve3::Bezier(rcad_kernel::geom::BezierCurve3 {
        control_points: extrap_poles.to_vec(),
        weights: vec![1.0; 4],
    });
    // OCCT: BRepLib_MakeEdge Bedge(Bezier); TopoDS_Edge edg = Bedge.Edge();
    // rcad: BRepLib_MakeEdge is the BRep_Builder edge vehicle
    // (BRep::add_tedge = the MakeEdge TEdgeData: curve + range, no vertices,
    // exactly the OCCT Bedge.Edge() state before it is added to a wire).
    let mut a_pool = topods::BRep::new();
    let a_bez_domain = rcad_kernel::geom::CurveEval::default_domain(&bezier);
    let edg = a_pool.add_tedge(
        Some(bezier.clone()),
        Shape::null(),
        Shape::null(),
        [a_bez_domain[0], a_bez_domain[1]],
    );
    let mut f_shape: Shape;
    *error = 1.0e-30;
    for nb in 1..=nbface {
        f_shape = face[(nb - 1) as usize].clone();
        // OCCT: NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>
        //       MapE1; TopoDS_Edge E1.
        let mut map_e1: Vec<Shape> = Vec::new();
        let mut e1 = Shape::null();
        let mut proj1: Option<rcad_kernel::geom::Curve2d> = None;
        let mut proj1c: Option<rcad_kernel::geom::Curve3> = None;
        // OCCT: BRepAlgo_NormalProjection OrtProj; OrtProj.Init(F);
        //       OrtProj.Add(edg); OrtProj.SetParams(1.e-4, 1.e-4, GeomAbs_C1,
        //       14, 16); OrtProj.Build();
        let mut ort_proj = BRepAlgoNormalProjection::new();
        ort_proj.init(&f_shape);
        ort_proj.add(&edg);
        ort_proj.set_params(1.0e-4, 1.0e-4, GeomAbsShape::C1, 14, 16);
        ort_proj.build();
        if ort_proj.is_done() {
            // OCCT: TopExp::MapShapes(OrtProj.Projection(), TopAbs_EDGE,
            //       MapE1) — the IndexedMap keeps the (TShape, Location)
            //       unique edges in explorer order (ShapeMapHasher).
            map_e1 = top_exp_map_shapes_edges(ort_proj.projection());
            if !map_e1.is_empty() {
                if map_e1.len() != 1 {
                    // OCCT: BRepLib_MakeFace Bface(BRep_Tool::Surface(F),
                    //       Precision::Confusion()); F = Bface.Face();
                    //       OrtProj.Init(F); OrtProj.Build(); MapE1.Clear();
                    // GAP: BRepLib_MakeFace (TKBRep/BRepLib) — the planar
                    // natural-bound face maker is not translated; F keeps the
                    // same surface (see the BRepLibMakeFace carrier note in
                    // offset/brep_offset_make_simple_offset.rs).
                    ort_proj.init(&f_shape);
                    ort_proj.build();
                    map_e1.clear();
                    if ort_proj.is_done() {
                        map_e1 = top_exp_map_shapes_edges(ort_proj.projection());
                    }
                }
                if !map_e1.is_empty() {
                    let mut trouve = false;
                    for ecur in &map_e1 {
                        if trouve {
                            break;
                        }
                        e1 = ecur.clone();
                        if !brep.is_edge_degenerated(&e1) {
                            trouve = true;
                        }
                    }
                    eproj.push(e1.clone());
                    // OCCT: proj1 = BRep_Tool::CurveOnSurface(E1, F, up1, up2);
                    match brep.curve_on_surface(&e1, &f_shape) {
                        Some((c, a, b)) => {
                            proj1 = Some(c);
                            up1 = a;
                            up2 = b;
                        }
                        None => {
                            proj1 = None;
                        }
                    }
                    if proj1.is_none() {
                        panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                    }
                    proj2d.push(Some(rcad_kernel::geom::Curve2d::Trimmed(
                        rcad_kernel::geom::TrimmedCurve2 {
                            curve: Box::new(proj1.clone().expect("proj1")),
                            t_min: up1,
                            t_max: up2,
                        },
                    )));
                    // OCCT: proj1c = BRep_Tool::Curve(E1, up1, up2);
                    match brep.edge_curve_world(&e1) {
                        Some((c, r)) => {
                            proj1c = Some(c);
                            up1 = r[0];
                            up2 = r[1];
                        }
                        None => {
                            proj1c = None;
                        }
                    }
                    if proj1c.is_none() {
                        panic!("Standard_ConstructionError: Failed to get 3D curve of edge");
                    }
                    // OCCT: cproj.Append(new Geom_TrimmedCurve(proj1c, up1, up2));
                    cproj.push(Some(rcad_kernel::geom::Curve3::Trimmed(
                        rcad_kernel::geom::TrimmedCurve3 {
                            curve: Box::new(proj1c.clone().expect("proj1c")),
                            first: up1,
                            last: up2,
                        },
                    )));
                    // OCCT: if (error > BRep_Tool::Tolerance(E1)) error = ...
                    let tol_e1 = e1.as_edge().map(|e| e.tolerance).unwrap_or(0.0);
                    if *error > tol_e1 {
                        *error = tol_e1;
                    }
                } else {
                    eproj.push(e1.clone());
                    proj2d.push(proj1.clone());
                    cproj.push(proj1c.clone());
                }
            } else {
                eproj.push(e1.clone());
                proj2d.push(proj1.clone());
                cproj.push(proj1c.clone());
            }
        }
    }
    // OCCT: for (nb = 1; nb <= nbface - 1; nb++) — update the parameters on
    //       the spine edges against the Hermite curve extremas.
    if nbface >= 2 {
        for nb in 1..=(nbface - 1) {
            let ecom_nb = ecom[(nb - 1) as usize].clone();
            // OCCT: BRepAdaptor_Curve C(TopoDS::Edge(Ecom.Value(nb)));
            //       C.D0(param.Value(nb), p02);
            let c = BRepAdaptorCurve::initialize(brep, &ecom_nb);
            p02 = c.value(param[(nb - 1) as usize]);
            // OCCT: GeomAdaptor_Curve L(Bezier);
            //       Extrema_ExtCC ext(C, L);
            let a_c_curve = match c.curve.clone() {
                Some(an_c) => an_c,
                None => panic!("Standard_ConstructionError: Failed to get 3D curve of edge"),
            };
            let a_c_adaptor = GeomCurveAdaptor::with_range(
                a_c_curve.clone(),
                c.first_parameter(),
                c.last_parameter(),
            );
            let a_l_adaptor = GeomCurveAdaptor::new(bezier.clone());
            let a_c_tool = CurveToolHandle::for_curve3(&a_c_curve, &a_c_adaptor, &a_c_adaptor);
            let a_l_tool = CurveToolHandle::for_curve3(&bezier, &a_l_adaptor, &a_l_adaptor);
            let ext = ExtremaExtCC::new_curves(&a_c_tool, &a_l_tool, 1.0e-10, 1.0e-10);
            if ext.is_done() {
                if !ext.is_parallel() && ext.nb_ext() != 0 {
                    // OCCT: Extrema_POnCurv POnC, POnL;
                    //       ext.Points(1, POnC, POnL);
                    let mut pon_c = POnCurve {
                        param: 0.0,
                        point: DVec3::ZERO,
                    };
                    let mut pon_l = POnCurve {
                        param: 0.0,
                        point: DVec3::ZERO,
                    };
                    ext.points(1, &mut pon_c, &mut pon_l);
                    if pon_c.point.distance(pon_l.point) < P_CONFUSION {
                        param[(nb - 1) as usize] = pon_c.param;
                    } else if let Some(cp) = cproj[(nb - 1) as usize].as_ref() {
                        p01 = cp.point_at(cp.default_domain()[1]);
                    } else if let Some(cp) = cproj[nb as usize].as_ref() {
                        p01 = cp.point_at(cp.default_domain()[0]);
                    }
                }
            }
            if !ext.is_done() || ext.nb_ext() == 0 {
                if let Some(cp) = cproj[(nb - 1) as usize].as_ref() {
                    p01 = cp.point_at(cp.default_domain()[1]);
                } else if let Some(cp) = cproj[nb as usize].as_ref() {
                    p01 = cp.point_at(cp.default_domain()[0]);
                }
                if p01.distance(p02) > 1.0e-4 {
                    // OCCT: Extrema_ExtPC ext1(p01, C);
                    let ext1 = ExtremaExtPC::new();
                    if ext1.is_done() && ext1.nb_ext() != 0 {
                        let pon_c = ext1.point(1);
                        param[(nb - 1) as usize] = pon_c.parameter();
                    }
                }
            }
        }
    }
    let _ = (p01, d11, d12);
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L530-543 — CalculDroite: calculate a 2D
// straight line passing through point p2d1 and direction xdir ydir.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L534-543
pub fn calcul_droite(p2d1: DVec2, xdir: f64, ydir: f64, pcurve: &mut Option<rcad_kernel::geom::Curve2d>) {
    let len = (xdir * xdir + ydir * ydir).sqrt();
    let dir = if len > 0.0 {
        DVec2::new(xdir / len, ydir / len)
    } else {
        DVec2::new(xdir, ydir)
    };
    let l = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
        origin: p2d1,
        direction: dir,
    });
    let l0 = len;
    *pcurve = Some(rcad_kernel::geom::Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
        curve: Box::new(l),
        t_min: 0.0,
        t_max: l0,
    }));
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L546-677 — CalculBatten: calcule a batten
// between curves 2d curv2d1 and curv2d2 at points p2d1 and p2d2.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L550-677
#[allow(clippy::too_many_arguments)]
pub fn calcul_batten(
    asurf: &GeomAdaptorSurface,
    face: &Shape,
    brep: &topods::BRep,
    xdir: f64,
    ydir: f64,
    p2d1: DVec2,
    p2d2: DVec2,
    contraint1: bool,
    contraint2: bool,
    curv2d1: &Option<rcad_kernel::geom::Curve2d>,
    curv2d2: &Option<rcad_kernel::geom::Curve2d>,
    picicplus: f64,
    picplusic: f64,
    inverseic: bool,
    inverseicplus: bool,
    pcurve: &mut Option<rcad_kernel::geom::Curve2d>,
) {
    let pi = std::f64::consts::PI;
    let mut isplane = asurf.get_type() == GeomAbsSurfaceType::Plane;
    let mut anglebig = false;
    let dir1 = DVec2::new(xdir, ydir);
    let cl1 = match curv2d1 {
        Some(c) => GeomLPropCLProps2d::new(c, picicplus, 1, 1.0e-4),
        None => GeomLPropCLProps2d {
            p: DVec2::ZERO,
            tangent: DVec2::ZERO,
        },
    };
    let cl2 = match curv2d2 {
        Some(c) => GeomLPropCLProps2d::new(c, picplusic, 1, 1.0e-4),
        None => GeomLPropCLProps2d {
            p: DVec2::ZERO,
            tangent: DVec2::ZERO,
        },
    };
    let mut dir3 = cl1.tangent();
    let mut dir4 = cl2.tangent();
    if inverseic {
        dir3 = -dir3;
    }
    if inverseicplus {
        dir4 = -dir4;
    }
    let h = p2d2.distance(p2d1) / 20.0;
    let mut bat = FairCurveBatten::new(p2d1, p2d2, h);
    bat.set_free_sliding(true);
    // OCCT: gp_Dir2d::Angle — the signed angle between two unit directions
    // in (-PI, PI].
    let angle_of = |a: DVec2, b: DVec2| -> f64 {
        let an = a.normalize_or_zero();
        let bn = b.normalize_or_zero();
        let mut ang = bn.y.atan2(bn.x) - an.y.atan2(an.x);
        while ang <= -pi {
            ang += 2.0 * pi;
        }
        while ang > pi {
            ang -= 2.0 * pi;
        }
        ang
    };
    let ang1 = angle_of(dir1, dir3);
    let dir1_angle_dir4 = angle_of(dir1, dir4);
    let ang2 = if dir1_angle_dir4 > 0.0 {
        pi - dir1_angle_dir4
    } else {
        -pi - dir1_angle_dir4
    };
    if contraint1 && contraint2 {
        anglebig = ang1.abs() > 1.2 || ang2.abs() > 1.2;
    } else if contraint1 {
        anglebig = ang1.abs() > 1.2;
    } else if contraint2 {
        anglebig = ang2.abs() > 1.2;
    }
    if isplane && (ang1.abs() > pi / 2.0 || ang2.abs() > pi / 2.0) {
        isplane = false;
    }
    if anglebig && !isplane {
        calcul_droite(p2d1, xdir, ydir, pcurve);
    } else {
        if contraint1 {
            bat.set_angle1(ang1);
        } else {
            bat.set_constraint_order1(0);
        }
        if contraint2 {
            bat.set_angle2(ang2);
        } else {
            bat.set_constraint_order2(0);
        }
        let mut iana = FairCurveAnalysisCode;
        let mut ok = bat.compute(&mut iana, 25, 1.0e-2);
        if !ok {
            bat.dump();
        }
        if ok {
            *pcurve = bat.curve();
            // OCCT: BRepTools::UVBounds(Face, umin, umax, vmin, vmax);
            let (umin, vmin, umax, vmax) = brep_tools_uv_bounds(brep, face);
            let mut bf = BndBox2d::new();
            let mut bc = BndBox2d::new();
            if let Some(pc) = pcurve.as_ref() {
                let acur = Geom2dAdaptorCurve::load(pc.clone(), pc.default_domain()[0], pc.default_domain()[1]);
                bnd_lib_add2d_curve(&acur, 0.0, &mut bc);
            }
            bf.update(umin, vmin, umax, vmax);
            let (uminc, vminc, umaxc, vmaxc) = bc.get();
            if uminc < umin - 1.0e-7 {
                ok = false;
            }
            if umaxc > umax + 1.0e-7 {
                ok = false;
            }
            if vminc < vmin - 1.0e-7 {
                ok = false;
            }
            if vmaxc > vmax + 1.0e-7 {
                ok = false;
            }
        }
        if !ok {
            calcul_droite(p2d1, xdir, ydir, pcurve);
        }
    }
}

/// OCCT BRepTools::UVBounds(F, U1, U2, V1, V2) — pending the wire pcurve
/// range scan; the unbounded domain is returned so the batten box check
/// (OCCT L647-671) passes through unchanged.
fn brep_tools_uv_bounds(brep: &topods::BRep, f: &Shape) -> (f64, f64, f64, f64) {
    match brep.face_surface(f) {
        Some(_) => (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::INFINITY),
        None => (0.0, 0.0, 0.0, 0.0),
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L680-715 — OrientationIcNonVive.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L685-715
pub fn orientation_ic_non_vive(
    cd_ic: &super::chfi_ds::ChFiDSStripe,
    jfic: i32,
    icicplus: i32,
    sensic: i32,
    orien: &mut Orientation,
) {
    let mut orinterf = Orientation::Forward;
    calcul_orientation(cd_ic, jfic, icicplus, &mut orinterf);
    if sensic != 1 {
        if orinterf == Orientation::Forward {
            *orien = Orientation::Forward;
        } else {
            *orien = Orientation::Reversed;
        }
    } else if orinterf == Orientation::Forward {
        *orien = Orientation::Reversed;
    } else {
        *orien = Orientation::Forward;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L717-754 — OrientationIcplusNonVive.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L723-754
pub fn orientation_icplus_non_vive(
    cd_icplus: &super::chfi_ds::ChFiDSStripe,
    jficplus: i32,
    icplusic: i32,
    sensicplus: i32,
    orien: &mut Orientation,
) {
    let mut orinterf = Orientation::Forward;
    let jfp = 3 - jficplus;
    calcul_orientation(cd_icplus, jfp, icplusic, &mut orinterf);
    if sensicplus == 1 {
        if orinterf == Orientation::Forward {
            *orien = Orientation::Forward;
        } else {
            *orien = Orientation::Reversed;
        }
    } else if orinterf == Orientation::Forward {
        *orien = Orientation::Reversed;
    } else {
        *orien = Orientation::Forward;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L756-806 — OrientationAreteViveConsecutive.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L762-806
pub fn orientation_arete_vive_consecutive(
    brep: &topods::BRep,
    fvive_icicplus: &Shape,
    evive_ic: &Shape,
    v1: &Shape,
    orien: &mut Orientation,
) {
    // orinterf is orientation of edge ic corresponding to face
    // Fviveicicplus taken FORWARD
    let mut orinterf = Orientation::Forward;
    let mut f = fvive_icicplus.clone();
    f.orientation = Orientation::Forward;
    let edges = topexp_face_edges(brep, &f);
    for ecur in &edges {
        if evive_ic.is_same(ecur) {
            orinterf = ecur.orientation;
            break;
        }
    }
    // if V1 is vertex REVERSED of edge ic the curve
    // has the same orientation as ic
    let vl = brep.last_vertex(evive_ic);
    if vl.is_same(v1) {
        if orinterf == Orientation::Forward {
            *orien = Orientation::Forward;
        } else {
            *orien = Orientation::Reversed;
        }
    } else if orinterf == Orientation::Forward {
        *orien = Orientation::Reversed;
    } else {
        *orien = Orientation::Forward;
    }
}

/// OCCT ChFiDS_SurfData::ChangeVertex(First, OnS) — OCCT returns a
/// reference; rcad copies (CommonPoint is a value payload here).
pub fn cp_change_vertex(fd: &ChFiDSSurfData, isfirst: bool, on_s: i32) -> ChFiDS_CommonPoint {
    fd.vertex(isfirst, on_s).clone()
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L809-951 — PerformTwoCornerSameExt:
// calculate intersection between two stripes stripe1 and stripe2.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L813-951
#[allow(clippy::too_many_arguments)]
pub fn perform_two_corner_same_ext(
    brep: &topods::BRep,
    dstr: &mut TopOpeBRepDSHDataStructure,
    stripe1: &SharedStripe,
    index1: i32,
    sens1: i32,
    stripe2: &SharedStripe,
    index2: i32,
    sens2: i32,
    trouve: &mut bool,
) {
    let isfirst1 = sens1 == 1;
    let com11 = {
        let st1 = stripe1.read().expect("stripe lock");
        let fd1 = st1.set_of_surf_data()[(index1 - 1) as usize].clone();
        let fdg = fd1.read().expect("surfdata lock");
        cp_change_vertex(&fdg, isfirst1, 1)
    };
    let com12 = {
        let st1 = stripe1.read().expect("stripe lock");
        let fd1 = st1.set_of_surf_data()[(index1 - 1) as usize].clone();
        let fdg = fd1.read().expect("surfdata lock");
        cp_change_vertex(&fdg, isfirst1, 2)
    };
    let isfirst2 = sens2 == 1;
    let com21 = {
        let st2 = stripe2.read().expect("stripe lock");
        let fd2 = st2.set_of_surf_data()[(index2 - 1) as usize].clone();
        let fdg = fd2.read().expect("surfdata lock");
        cp_change_vertex(&fdg, isfirst2, 1)
    };
    let indic1 = {
        let st1 = stripe1.read().expect("stripe lock");
        st1.set_of_surf_data()[(index1 - 1) as usize]
            .read()
            .expect("surfdata lock")
            .surf()
    };
    let indic2 = {
        let st2 = stripe2.read().expect("stripe lock");
        st2.set_of_surf_data()[(index2 - 1) as usize]
            .read()
            .expect("surfdata lock")
            .surf()
    };
    let fd1_shared = {
        let st1 = stripe1.read().expect("stripe lock");
        st1.set_of_surf_data()[(index1 - 1) as usize].clone()
    };
    let fd2_shared = {
        let st2 = stripe2.read().expect("stripe lock");
        st2.set_of_surf_data()[(index2 - 1) as usize].clone()
    };
    let f1 = fd1_shared.read().expect("surfdata lock");
    let f2 = fd2_shared.read().expect("surfdata lock");

    let fi11 = f1.interference_on_s1();
    let fi12 = f1.interference_on_s2();
    let fi21 = f2.interference_on_s1();
    let fi22 = f2.interference_on_s2();
    let isfirst = sens1 == 1;
    let pfi11 = fi11
        .pcurve_on_surf()
        .map(|pc| pc.point_at(fi11.parameter(isfirst)))
        .unwrap_or(DVec2::ZERO);
    let pfi12 = fi12
        .pcurve_on_surf()
        .map(|pc| pc.point_at(fi12.parameter(isfirst)))
        .unwrap_or(DVec2::ZERO);
    let isfirst = sens2 == 1;
    let (pfi21, pfi22) = if com11.point().distance(com21.point()) < 1.0e-4 {
        (
            fi21.pcurve_on_surf()
                .map(|pc| pc.point_at(fi21.parameter(isfirst)))
                .unwrap_or(DVec2::ZERO),
            fi22.pcurve_on_surf()
                .map(|pc| pc.point_at(fi22.parameter(isfirst)))
                .unwrap_or(DVec2::ZERO),
        )
    } else {
        (
            fi22.pcurve_on_surf()
                .map(|pc| pc.point_at(fi22.parameter(isfirst)))
                .unwrap_or(DVec2::ZERO),
            fi21.pcurve_on_surf()
                .map(|pc| pc.point_at(fi21.parameter(isfirst)))
                .unwrap_or(DVec2::ZERO),
        )
    };

    // OCCT: NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);
    let pardeb: [f64; 4] = [pfi11.x, pfi11.y, pfi21.x, pfi21.y];
    let parfin: [f64; 4] = [pfi12.x, pfi12.y, pfi22.x, pfi22.y];

    // OCCT: HS1 = ChFi3d_BoundSurf(DStr, Fd1, 1, 2); HS2 = ...;
    let hs1 = chfi3d_bound_surf(dstr, &f1, 1, 2);
    let hs2 = chfi3d_bound_surf(dstr, &f2, 1, 2);
    *trouve = false;
    let mut cint: Option<rcad_kernel::geom::Curve3> = None;
    let mut c2dint1: Option<rcad_kernel::geom::Curve2d> = None;
    let mut c2dint2: Option<rcad_kernel::geom::Curve2d> = None;
    let mut tol = 0.0f64;
    if let Some(cc) = chfi3d_compute_curves(&hs1, &hs2, pardeb, parfin, 1.0e-4, 1.0e-5, &mut tol) {
        cint = Some(cc.c3d);
        c2dint1 = Some(cc.pc1);
        c2dint2 = Some(cc.pc2);
        let cint_ref = cint.as_ref().expect("cint");
        let p1 = cint_ref.point_at(cint_ref.default_domain()[0]);
        let p2 = cint_ref.point_at(cint_ref.default_domain()[1]);
        *trouve = (com11.point().distance(p1) < 1.0e-4 || com11.point().distance(p2) < 1.0e-4)
            && (com12.point().distance(p1) < 1.0e-4 || com12.point().distance(p2) < 1.0e-4);
    }
    drop(f1);
    drop(f2);

    if *trouve {
        let isfirst = sens1 == 1;
        stripe1.write().expect("stripe lock").in_ds(isfirst, 1);
        let mut indpoint1 = super::chfi3d::chfi3d_index_point_in_ds(&com11, dstr);
        let mut indpoint2 = super::chfi3d::chfi3d_index_point_in_ds(&com12, dstr);
        {
            let mut w1 = stripe1.write().expect("stripe lock");
            w1.set_index_point(indpoint1, isfirst, 1);
            w1.set_index_point(indpoint2, isfirst, 2);
        }
        let isfirst = sens2 == 1;
        stripe2.write().expect("stripe lock").in_ds(isfirst, 1);
        if com11.point().distance(com21.point()) < 1.0e-4 {
            let mut w2 = stripe2.write().expect("stripe lock");
            w2.set_index_point(indpoint1, isfirst, 1);
            w2.set_index_point(indpoint2, isfirst, 2);
        } else {
            let mut w2 = stripe2.write().expect("stripe lock");
            w2.set_index_point(indpoint2, isfirst, 1);
            w2.set_index_point(indpoint1, isfirst, 2);
        }

        let orsurf1 = fd1_shared.read().expect("surfdata lock").orientation();
        let s1_index = fd1_shared.read().expect("surfdata lock").index_of_s1();
        let fi11_transition = fd1_shared
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .transition();
        let fd1_orientation = fd1_shared.read().expect("surfdata lock").orientation();
        let mut trafil1 = dstr.shape(s1_index).orientation;
        trafil1 = topabs_compose(trafil1, fd1_orientation);
        trafil1 = topabs_compose(topabs_reverse(fi11_transition), trafil1);
        // OCCT: orsurf2 = Fd2->Orientation();
        let orsurf2 = fd2_shared.read().expect("surfdata lock").orientation();
        // OCCT: TopOpeBRepDS_Curve tcurv3d(cint, tol); indcurve =
        //       DStr.AddCurve(tcurv3d);
        let indcurve = dstr.add_curve(super::chfi3d_ds::TopOpeBRepDSCurve::new(cint.clone(), tol));
        let cint_ref = cint.as_ref().expect("cint");
        let p1 = cint_ref.point_at(cint_ref.default_domain()[0]);
        let p2 = cint_ref.point_at(cint_ref.default_domain()[1]);
        // OCCT: Fi11.PCurveOnFace()->D0(...); Stemp->D0(p2d.X(), p2d.Y(), ...);
        let fi11_pc = fd1_shared
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .pcurve_on_face()
            .cloned();
        let stemp = brep
            .face_surface(dstr.shape(s1_index))
            .cloned();
        let lastp = fd1_shared
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .parameter_last();
        let firstp = fd1_shared
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .parameter_first();
        let p4 = match (&fi11_pc, &stemp) {
            (Some(pc), Some(s)) => {
                let p2d = pc.point_at(lastp);
                s.point_at(p2d.x, p2d.y)
            }
            _ => DVec3::ZERO,
        };
        let p3 = match (&fi11_pc, &stemp) {
            (Some(pc), Some(s)) => {
                let p2d = pc.point_at(firstp);
                s.point_at(p2d.x, p2d.y)
            }
            _ => DVec3::ZERO,
        };
        let mut orpcurve;
        if p1.distance(p4) < 1.0e-4 || p2.distance(p3) < 1.0e-4 {
            orpcurve = trafil1;
        } else {
            orpcurve = topabs_reverse(trafil1);
        }
        if com11.point().distance(p1) > 1.0e-4 {
            let ind = indpoint1;
            indpoint1 = indpoint2;
            indpoint2 = ind;
        }
        let interfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Forward,
            indcurve,
            indpoint1,
            cint_ref.default_domain()[0],
            false,
        );
        let interfp2 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
            Orientation::Reversed,
            indcurve,
            indpoint2,
            cint_ref.default_domain()[1],
            false,
        );
        dstr.change_curve_interferences(indcurve).push(interfp1);
        dstr.change_curve_interferences(indcurve).push(interfp2);
        let interfc =
            super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(indcurve, indic1, c2dint1.clone(), orpcurve);
        dstr.change_surface_interferences(indic1).push(interfc);
        if orsurf1 == orsurf2 {
            orpcurve = topabs_reverse(orpcurve);
        }
        let interfc =
            super::chfi3d_builder_0::chfi3d_fil_curve_in_ds(indcurve, indic2, c2dint2.clone(), orpcurve);
        dstr.change_surface_interferences(indic2).push(interfc);
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L954-984 — CpOnEdge.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L958-984
fn cp_on_edge(
    stripe: &super::chfi_ds::ChFiDSStripe,
    num: i32,
    isfirst: bool,
    eadj1: &Shape,
    eadj2: &Shape,
    compoint: &mut bool,
) {
    let fd = stripe.set_of_surf_data()[(num - 1) as usize]
        .read()
        .expect("surfdata lock");
    let cp1 = cp_change_vertex(&fd, isfirst, 1);
    let cp2 = cp_change_vertex(&fd, isfirst, 2);
    *compoint = false;
    if cp1.is_on_arc() && (cp1.arc().is_same(eadj1) || cp1.arc().is_same(eadj2)) {
        *compoint = true;
    }
    if cp2.is_on_arc() && (cp2.arc().is_same(eadj1) || cp2.arc().is_same(eadj2)) {
        *compoint = true;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L986-1080 — RemoveSurfData: for each stripe
// removal of unused surfdatas.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L991-1080
pub fn remove_surf_data(
    brep: &topods::BRep,
    my_vdata_map: &ChFiDSStripeMap,
    my_ef_map: &ChFiDSMap,
    edgecouture: &Shape,
    facecouture: &Shape,
    v1: &Shape,
) {
    for it in stripe_map_find(my_vdata_map, v1) {
        let nbsurf = it.read().expect("stripe lock").set_of_surf_data().len();
        let nbedge = it
            .read()
            .expect("stripe lock")
            .spine()
            .map(|sp| sp.base().nb_edges())
            .unwrap_or(0);
        if nbsurf != 1 {
            let mut sense = 0i32;
            let num = {
                let st = it.read().expect("stripe lock");
                chfi3d_index_of_surf_data(v1, &st, &mut sense)
            };
            let ecur = {
                let st = it.read().expect("stripe lock");
                let sp = st.spine().expect("spine");
                if sense == 1 {
                    sp.base().edges(1).clone()
                } else {
                    sp.base().edges(nbedge).clone()
                }
            };
            let (mut f1, mut f2) = (Shape::null(), Shape::null());
            chfi3d_edge_common_faces(my_ef_map.find(&ecur), &mut f1, &mut f2);
            let mut eadj1 = Shape::null();
            let mut eadj2 = Shape::null();
            if f1.is_same(facecouture) {
                eadj1 = edgecouture.clone();
            } else {
                let mut vbid = Shape::null();
                chfi3d_cherche_element(brep, v1, &ecur, &f1, &mut eadj1, &mut vbid);
            }
            let (mut fg, mut fd) = (Shape::null(), Shape::null());
            chfi3d_edge_common_faces(my_ef_map.find(&eadj1), &mut fg, &mut fd);
            if f2.is_same(facecouture) {
                eadj2 = edgecouture.clone();
            } else {
                let mut vbid = Shape::null();
                chfi3d_cherche_element(brep, v1, &ecur, &f2, &mut eadj2, &mut vbid);
            }
            chfi3d_edge_common_faces(my_ef_map.find(&eadj2), &mut fg, &mut fd);
            let mut compoint = false;
            let isfirst = sense == 1;
            if sense == 1 {
                let mut ind = 0i32;
                // among surfdatas find the greatest indice ind so that
                // surfdata could have one of commonpoint on Eadj1 and Eadj2
                // remove surfdata from 1 to ind-1
                for i in 1..=(nbsurf as i32) {
                    let st = it.read().expect("stripe lock");
                    cp_on_edge(&st, i, isfirst, &eadj1, &eadj2, &mut compoint);
                    if compoint {
                        ind = i;
                    }
                }
                if ind >= 2 {
                    let mut w = it.write().expect("stripe lock");
                    remove_sd(&mut w, 1, ind - 1);
                }
            } else {
                let mut ind = num;
                // among surfdatas find the smallest indice ind so that
                // surfdata could have one of commonpoint on Eadj1 and Eadj2
                // remove surfdata from ind+1 to num
                for i in (1..=num).rev() {
                    let st = it.read().expect("stripe lock");
                    cp_on_edge(&st, i, isfirst, &eadj1, &eadj2, &mut compoint);
                    if compoint {
                        ind = i;
                    }
                }
                if ind < num {
                    let mut w = it.write().expect("stripe lock");
                    remove_sd(&mut w, ind + 1, num);
                }
            }
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L1084-1109 — ParametrePlate.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L1084-1109
pub fn parametre_plate(
    n3d: i32,
    psurf: &GeomPlateBuildPlateSurface,
    surf: &rcad_kernel::geom::Surface3,
    point: DVec3,
    apperror: f64,
    uv: &mut DVec2,
) {
    let mut p1: DVec3;
    let mut par;
    let mut trouve = false;
    for ip in 1..=n3d {
        if trouve {
            break;
        }
        let Some(c) = psurf.curves2d_value(ip) else {
            continue;
        };
        par = c.default_domain()[0];
        *uv = c.point_at(par);
        p1 = surf.point_at(uv.x, uv.y);
        trouve = p1.abs_diff_eq(point, apperror);
        if !trouve {
            par = c.default_domain()[1];
            *uv = c.point_at(par);
            p1 = surf.point_at(uv.x, uv.y);
            trouve = p1.abs_diff_eq(point, apperror);
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L1113-1141 — SummarizeNormal.
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L1113-1141
pub fn summarize_normal(
    brep: &topods::BRep,
    v1: &Shape,
    fcur: &Shape,
    ecur: &Shape,
    sum_face_normal_at_v1: &mut DVec3,
) {
    let Some((c2d, fp, lp)) = brep.curve_on_surface(ecur, fcur) else {
        return;
    };
    // OCCT: BRep_Tool::UVPoints(Ecur, Fcur, uv1, uv2);
    let mut uv1 = c2d.point_at(fp);
    let uv2 = c2d.point_at(lp);
    let first_v = brep.first_vertex(ecur);
    if !v1.is_same(&first_v) {
        uv1 = uv2;
    }

    let Some(s) = brep.face_surface(fcur) else {
        return;
    };
    // OCCT: BRep_Tool::Surface(Fcur)->D1(uv1.X(), uv1.Y(), P, d1U, d1V);
    let (_p, d1u, d1v) = s.derivatives(uv1.x, uv1.y);
    let mut n = d1u.cross(d1v);
    if fcur.orientation == Orientation::Reversed {
        n = -n;
    }

    if n.length_squared() <= P_CONFUSION {
        return;
    }

    *sum_face_normal_at_v1 += n.normalize();
    *sum_face_normal_at_v1 = sum_face_normal_at_v1.normalize();
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L1143-1170 — ChFi3d_SurfType + SurfIndex.
// =========================================================================

/// OCCT enum ChFi3d_SurfType (L1143-1148) — for call SurfIndex(...).
#[derive(Debug, Clone, Copy)]
pub enum ChFi3dSurfType {
    ChFiSURFACE,
    FACE1,
    FACE2,
}

// OCCT ChFi3d_Builder_CnCrn.cxx L1152-1170
pub fn surf_index(
    stripe_array1: &[SharedStripe],
    stripe_index: usize,
    surf_data_index: i32,
    surf_type: ChFi3dSurfType,
) -> i32 {
    let st = stripe_array1[stripe_index].read().expect("stripe lock");
    let a_surf_data = st.set_of_surf_data()[(surf_data_index - 1) as usize]
        .read()
        .expect("surfdata lock");
    match surf_type {
        ChFi3dSurfType::ChFiSURFACE => a_surf_data.surf(),
        ChFi3dSurfType::FACE1 => a_surf_data.index_of_s1(),
        ChFi3dSurfType::FACE2 => a_surf_data.index_of_s2(),
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_CnCrn.cxx L1172-1230 — PlateOrientation: define
// Plate orientation compared to <theRefDir> previewing that Plate surface
// can have a sharp angle with adjacent filet (bug occ266: 2 chamfs, OnSame
// and OnDiff) and can be even twisted (grid tests cfi900 B1).
// =========================================================================

// OCCT ChFi3d_Builder_CnCrn.cxx L1180-1230
pub fn plate_orientation(
    the_plate_surf: &rcad_kernel::geom::Surface3,
    the_pc_arr: &[Option<rcad_kernel::geom::Curve2d>],
    the_ref_dir: DVec3,
) -> Orientation {
    let mut pp1 = DVec3::ZERO;
    let mut pp2;
    let mut pp3;
    let mut uv;
    let mut sum_scal1 = 0.0f64;
    let mut sum_scal2 = 0.0f64;

    let nb = the_pc_arr.len() as i32;
    if nb == 0 {
        return Orientation::Forward;
    }
    let a_pc0 = the_pc_arr[(nb - 1) as usize].clone();
    let Some(a_pc) = a_pc0 else {
        return Orientation::Forward;
    };
    let fpar = a_pc.default_domain()[0];
    let lpar = a_pc.default_domain()[1];
    uv = a_pc.point_at((fpar + lpar) / 2.0);
    pp1 = the_plate_surf.point_at(uv.x, uv.y);
    uv = a_pc.point_at(lpar);
    pp2 = the_plate_surf.point_at(uv.x, uv.y);

    for i in 1..=nb {
        let Some(pc_i) = the_pc_arr[(i - 1) as usize].clone() else {
            continue;
        };
        let fpar = pc_i.default_domain()[0];
        let lpar = pc_i.default_domain()[1];
        uv = pc_i.point_at(fpar);
        let (pp2v, du, dv) = the_plate_surf.derivatives(uv.x, uv.y);
        pp2 = pp2v;
        let n1 = du.cross(dv).normalize();

        uv = pc_i.point_at((fpar + lpar) / 2.0);
        pp3 = the_plate_surf.point_at(uv.x, uv.y);

        let vv1 = pp1 - pp2;
        let vv2 = pp3 - pp2;
        let n2 = vv2.cross(vv1).normalize();

        sum_scal1 += n1.dot(n2);
        sum_scal2 += n2.dot(the_ref_dir);

        pp1 = pp3;
    }
    if sum_scal2 * sum_scal1 > 0.0 {
        Orientation::Forward
    } else {
        Orientation::Reversed
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1222-1230 — recadre (pending relocation).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L1222-1230
pub fn recadre(p: f64, ref_p: f64, sens: i32, first: f64, last: f64) -> f64 {
    let pp = p + if sens > 0 { first - last } else { last - first };
    if (pp - ref_p).abs() < (p - ref_p).abs() {
        pp
    } else {
        p
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1234-1447 — ChFi3d_IntTraces (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L1234-1447
#[allow(clippy::too_many_arguments)]
pub fn chfi3d_int_traces(
    fd1: &ChFiDSSurfData,
    pref1: f64,
    p1: &mut f64,
    jf1: i32,
    sens1: i32,
    fd2: &ChFiDSSurfData,
    pref2: f64,
    p2: &mut f64,
    jf2: i32,
    sens2: i32,
    ref_p2d: DVec2,
    check2d_distance: bool,
    enlarge: bool,
) -> bool {
    // pcurves are enlarged to be sure that there is intersection
    // additionally all periodic curves are taken and points on
    // them are filtered using a specific criterion.
    let delta0 = 0.0f64;

    let fi1 = fd1.interference(jf1);
    let first = fi1.parameter_first();
    let last = fi1.parameter_last();
    if (last - first) < P_CONFUSION {
        return false;
    }
    let mut delta = delta0;
    if enlarge {
        delta = 0.1f64.min(0.05 * (last - first));
    }
    let Some(pcf1) = fi1.pcurve_on_face().cloned() else {
        return false;
    };
    let isper1 = pcf1.is_periodic();
    // OCCT: if the pcurve is periodic the basis curve of the trimmed curve
    // is loaded (occ::down_cast<Geom2d_TrimmedCurve>) — the rcad Curve2d
    // carries the full-range expression; the load keeps the curve as-is.
    let c1 = if isper1 {
        Geom2dAdaptorCurve::load(pcf1.clone(), pcf1.default_domain()[0], pcf1.default_domain()[1])
    } else {
        Geom2dAdaptorCurve::load(pcf1.clone(), first - delta, last + delta)
    };
    let first1 = pcf1.default_domain()[0];
    let last1 = pcf1.default_domain()[1];

    let fi2 = fd2.interference(jf2);
    let first = fi2.parameter_first();
    let last = fi2.parameter_last();
    if (last - first) < P_CONFUSION {
        return false;
    }
    if enlarge {
        delta = 0.1f64.min(0.05 * (last - first));
    }
    let Some(pcf2) = fi2.pcurve_on_face().cloned() else {
        return false;
    };
    let isper2 = pcf2.is_periodic();
    let c2 = if isper2 {
        Geom2dAdaptorCurve::load(pcf2.clone(), pcf2.default_domain()[0], pcf2.default_domain()[1])
    } else {
        // OCCT loads the original fd2 pcurve (not pcf2) here.
        Geom2dAdaptorCurve::load(pcf2.clone(), first - delta, last + delta)
    };
    let first2 = pcf2.default_domain()[0];
    let last2 = pcf2.default_domain()[1];

    // OCCT: handle identity PCurveOnFace() == PCurveOnFace() — rcad Curve2d
    // has no handle identity; the generic two-curve perform is used.
    let mut intersection = Geom2dIntGInter::new();
    intersection.perform2(&c1, &c2, P_CONFUSION, P_CONFUSION);
    if intersection.is_done() {
        if !intersection.is_empty() {
            let nbseg = intersection.nb_segments();
            if nbseg > 0 {
                // OCCT: no processing of the tangent segments.
            }
            let nbpt = intersection.nb_points();
            if nbpt >= 1 {
                // The criteria sets to filter the found points in a strict
                // way are missing. Two different criterions chosen somewhat
                // randomly are used :
                // - periodic curves : closest to the border.
                // - non-periodic curves : the closest to the left of 2
                //   curves modulo sens1 and sens2
                let (mut p2d, par1, par2) = intersection.point(1);
                *p1 = par1;
                *p2 = par2;
                if isper1 {
                    *p1 = recadre(*p1, pref1, sens1, first1, last1);
                }
                if isper2 {
                    *p2 = recadre(*p2, pref2, sens2, first2, last2);
                }
                for i in 2..=nbpt {
                    let (int2d_value, int2d_par1, int2d_par2) = intersection.point(i);
                    if isper1 {
                        let mut pp1 = int2d_par1;
                        pp1 = recadre(pp1, pref1, sens1, first1, last1);
                        if (pp1 - pref1).abs() < (*p1 - pref1).abs() {
                            *p1 = pp1;
                            *p2 = int2d_par2;
                            p2d = int2d_value;
                        } else if check2d_distance
                            && ref_p2d.distance(int2d_value) < ref_p2d.distance(p2d)
                        {
                            // Modified by skv - Mon Jun 16 15:51:21 2003 OCC615
                            let mut pp2 = int2d_par2;
                            if isper2 {
                                pp2 = recadre(pp2, pref2, sens2, first2, last2);
                            }
                            *p1 = pp1;
                            *p2 = pp2;
                            p2d = int2d_value;
                        }
                    } else if isper2 {
                        let mut pp2 = int2d_par2;
                        pp2 = recadre(pp2, pref2, sens2, first2, last2);
                        if (pp2 - pref2).abs() < (*p2 - pref2).abs() {
                            *p2 = pp2;
                            *p1 = int2d_par1;
                            p2d = int2d_value;
                        } else if check2d_distance
                            && ref_p2d.distance(int2d_value) < ref_p2d.distance(p2d)
                        {
                            // Modified by skv - Mon Jun 16 15:51:21 2003 OCC615
                            let mut pp1 = int2d_par1;
                            if isper1 {
                                pp1 = recadre(pp1, pref1, sens1, first1, last1);
                            }
                            *p1 = pp1;
                            *p2 = pp2;
                            p2d = int2d_value;
                        }
                    } else if ((int2d_par1 - *p1) * (sens1 as f64) < 0.0)
                        && ((int2d_par2 - *p2) * (sens2 as f64) < 0.0)
                    {
                        *p1 = int2d_par1;
                        *p2 = int2d_par2;
                        p2d = int2d_value;
                    } else if ((int2d_par1 - pref1).abs() < (*p1 - pref1).abs())
                        && ((int2d_par2 - pref2).abs() < (*p2 - pref2).abs())
                    {
                        *p1 = int2d_par1;
                        *p2 = int2d_par2;
                        p2d = int2d_value;
                    } else if check2d_distance
                        && ref_p2d.distance(int2d_value) < ref_p2d.distance(p2d)
                    {
                        *p1 = int2d_par1;
                        *p2 = int2d_par2;
                        p2d = int2d_value;
                    }
                }
                return true;
            }
            return false;
        }
        false
    } else {
        false
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L909-1218 — ChFi3d_IsInFront (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L909-1218
#[allow(clippy::too_many_arguments)]
pub fn chfi3d_is_in_front(
    brep: &topods::BRep,
    dstr: &TopOpeBRepDSHDataStructure,
    cd1: &SharedStripe,
    cd2: &SharedStripe,
    i1: i32,
    i2: i32,
    sens1: i32,
    sens2: i32,
    p1: &mut f64,
    p2: &mut f64,
    face: &mut Shape,
    sameside: &mut bool,
    jf1: &mut i32,
    jf2: &mut i32,
    visavis: &mut bool,
    vtx: &Shape,
    check2d_distance: bool,
    enlarge: bool,
) -> bool {
    let _ = (brep, vtx);
    let isf1 = sens1 == 1;
    let isf2 = sens2 == 1;
    let fd1 = cd1
        .read()
        .expect("stripe lock")
        .set_of_surf_data()[(i1 - 1) as usize]
        .clone();
    let fd2 = cd2
        .read()
        .expect("stripe lock")
        .set_of_surf_data()[(i2 - 1) as usize]
        .clone();

    let mut or;
    let mut or_face1;
    let mut or_face2;
    *visavis = false;
    let mut u1 = 0.0f64;
    let mut u2 = 0.0f64;
    let mut ss = false;
    let mut ok = false;
    let mut j1 = 0i32;
    let mut j2 = 0i32;
    let mut ff = Shape::null();
    let f1g = fd1.read().expect("surfdata lock");
    let f2g = fd2.read().expect("surfdata lock");
    let or1_save_slot = cd1.read().expect("stripe lock");
    let or2_save_slot = cd2.read().expect("stripe lock");
    // The four fd1/fd2 side combinations of OCCT L937-1216.  Each block:
    // face = DStr.Shape(fd->Index(jf)) with the NULL checks, the stripe and
    // face orientations, sameside, then ChFi3d_IntTraces with restore.
    for pass in 0..4 {
        let (jf1_new, jf2_new, index1, index2) = match pass {
            0 => (1, 1, f1g.index_of_s1(), f2g.index_of_s1()),
            1 => (2, 1, f1g.index_of_s2(), f2g.index_of_s1()),
            2 => (1, 2, f1g.index_of_s1(), f2g.index_of_s2()),
            _ => (2, 2, f1g.index_of_s2(), f2g.index_of_s2()),
        };
        if index1 != index2 {
            continue;
        }
        let restore_pass = pass != 0;
        *jf1 = jf1_new;
        *jf2 = jf2_new;
        *face = dstr.shape(f1g.index_of(*jf1)).clone();
        if face.is_null() {
            panic!("Standard_NullObject: ChFi3d_IsInFront : Trying to check orientation of NULL face");
        }
        let or_save1 = or1_save_slot.orientation_on_s(*jf1);
        or = face.orientation;
        or_face1 = or;
        let or_save2 = or2_save_slot.orientation_on_s(*jf2);
        let shape2 = dstr.shape(f2g.index_of(*jf2));
        if shape2.is_null() {
            panic!("Standard_NullObject: ChFi3d_IsInFront : Trying to check orientation of NULL shape");
        }
        or_face2 = shape2.orientation;
        *visavis = true;
        *sameside = chfi3d_same_side(or, or_save1, or_save2, or_face1, or_face2);
        // The parameters of the other side are not used for orientation.
        // This would raise problems
        let kf1 = *jf1;
        let kf2 = *jf2;
        let pref1 = f1g.interference(kf1).parameter(isf1);
        let pref2 = f2g.interference(kf2).parameter(isf2);
        // OCCT: P2d = BRep_Tool::Parameters(Vtx, face) — pending.
        let p2d = DVec2::ZERO;
        let _ = (or_save1, or_save2, or_face1, or_face2);
        if chfi3d_int_traces(
            &f1g, pref1, p1, *jf1, sens1, &f2g, pref2, p2, *jf2, sens2, p2d, check2d_distance,
            enlarge,
        ) {
            if restore_pass {
                let restore = ok
                    && ((j1 == *jf1 && (sens1 as f64) * (*p1 - u1) > 0.0)
                        || (j2 == *jf2 && (sens2 as f64) * (*p2 - u2) > 0.0));
                ok = true;
                if restore {
                    *p1 = u1;
                    *p2 = u2;
                    *sameside = ss;
                    *jf1 = j1;
                    *jf2 = j2;
                    *face = ff.clone();
                } else {
                    u1 = *p1;
                    u2 = *p2;
                    ss = *sameside;
                    j1 = *jf1;
                    j2 = *jf2;
                    ff = face.clone();
                }
            } else {
                u1 = *p1;
                u2 = *p2;
                ss = *sameside;
                j1 = *jf1;
                j2 = *jf2;
                ff = face.clone();
                ok = true;
            }
        } else if ok && restore_pass {
            // the re-initialization is added in case p1,... take wrong values
            *p1 = u1;
            *p2 = u2;
            *sameside = ss;
            *jf1 = j1;
            *jf2 = j2;
            *face = ff.clone();
        }
    }
    ok
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L4549-4663 — ChFi3d_SearchFD (pending
// relocation into chfi3d_builder_0.rs).
// =========================================================================

// OCCT ChFi3d_Builder_0.cxx L4549-4663
#[allow(clippy::too_many_arguments)]
pub fn chfi3d_search_fd(
    brep: &topods::BRep,
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
    let l1 = cd1.read().expect("stripe lock").set_of_surf_data().len() as i32;
    let l2 = cd2.read().expect("stripe lock").set_of_surf_data().len() as i32;
    let mut fini1 = false;
    let mut fini2 = false;
    let mut visavis = false;
    let vtx = Shape::null();
    while !found {
        let mut i = id1;
        while (i * sens1) <= (if1 * sens1) && !found && !fini2 {
            if chfi3d_is_in_front(
                brep, dstr, cd1, cd2, i, if2, sens1, sens2, p1, p2, face, sameside, jf1, jf2,
                &mut visavis, &vtx, false, false,
            ) {
                *i1 = i;
                *i2 = if2;
                found = true;
            } else if visavis {
                // OCCT: visavis && !visavisok — the ok latch.
                *i1 = i;
                *i2 = if2;
            }
            i += sens1;
        }
        if !fini1 {
            if1 += sens1;
            if if1 < 1 || if1 > l1 {
                if1 -= sens1;
                fini1 = true;
            }
        }

        let mut i = id2;
        while (i * sens2) <= (if2 * sens2) && !found && !fini1 {
            if chfi3d_is_in_front(
                brep, dstr, cd1, cd2, if1, i, sens1, sens2, p1, p2, face, sameside, jf1, jf2,
                &mut visavis, &vtx, false, false,
            ) {
                *i1 = if1;
                *i2 = i;
                found = true;
            } else if visavis {
                *i1 = if1;
                *i2 = i;
            }
            i += sens2;
        }
        if !fini2 {
            if2 += sens2;
            if if2 < 1 || if2 > l2 {
                if2 -= sens2;
                fini2 = true;
            }
        }
        if fini1 && fini2 {
            break;
        }
    }
    found
}
