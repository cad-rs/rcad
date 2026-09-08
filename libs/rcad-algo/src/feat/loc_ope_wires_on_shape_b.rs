// OCCT LocOpe_WiresOnShape.cxx L17-1623 — static helpers part (module b).
//
// This module carries the BRep_Tool/TopExp re-hosts shared with the class
// module, the pure-math re-hosts and the four cxx statics:
//   Project(V, p2d, F, theEdge, param)    cxx L491-675
//   Project(V, theEdge)                   cxx L679-697
//   Project(V, p2d, theEdge, theFace)     cxx L701-718
//   PutPCurve(Edg, Fac)                   cxx L722-909
//   PutPCurves(Efrom, Eto, myShape)       cxx L913-1287
//   FindInternalIntersections(...)        cxx L1291-1511
//
// Architecture differences (numbering continues from the class module):
// 6. ShapeAnalysis::AdjustToPeriod/AdjustByPeriod -> exact re-host below
//    (ShapeAnalysis.cxx L44-69).
// 7. Bnd_Box2d + BndLib_Add2dCurve::Add(PCurve, Tol, Box) -> BndBox2d filled
//    by base::bnd_lib::curve2d_bounding_box (the BRepTools::AddUVBounds
//    vehicle, already OCCT-aligned).
// 8. BRepTopAdaptor_FClass2d(face, Tol) -> topalgo::brep_top_adaptor::
//    fclass2d::FClass2d over FaceShapeSource (the rcad FClass2d is the
//    IntTools_FClass2d-equivalent classifier built from a bare face —
//    shape_source.rs; bop/int_tools/edge_face.rs precedent).
// 9. ShapeConstruct_ProjectCurveOnSurface (TKShHealing) is not translated
//    yet; re-hosted below with Perform failing (leaving the pcurve null) so
//    PutPCurve takes the OCCT null-curve early return (cxx L800-803).
//    GAP: closes with the TKShHealing batch.
// 10. GeomProjLib::Curve2d (TKTopAlgo) is not translated yet; the PutPCurves
//    call sites keep the OCCT structure and take the null-curve exit. Note:
//    OCCT L1227 dereferences the result without a null check on the seam
//    path; rcad returns early instead of crashing (marked at the call site).
//    GAP: closes with the GeomProjLib batch.
// 11. Extrema_ExtCC (TKGeomBase; the NbExt/IsParallel/TrimmedSquareDistances/
//    Points surface) is not translated yet; re-hosted below as ExtremaExtCC
//    with IsDone()=false (chfi3d_builder_cncrn.rs precedent) so
//    FindInternalIntersections takes the OCCT !IsDone() continue (cxx
//    L1347-1350) for every face edge. GAP: closes with the TKGeomAlgo
//    ExtCC batch.
// 12. ShapeAnalysis_Edge::CheckSameParameter -> reduced re-host
//    shape_analysis_edge_check_same_parameter below (pcurve/surface sampling
//    max-deviation, the same reduction as BRep::same_parameter in
//    rcad-kernel). GAP: the OCCT BRepLib_ValidateEdge walk.
// 13. Standard_Real Epsilon(theValue) (Standard_Real.hxx L238-247) ->
//    machine-epsilon re-host below (pure math).
// 14. BRepTools::UVBounds (BRepTools.cxx L137-181 via AddUVBounds) and
//    BRep_Tool::IsClosed(E, F) -> re-hosts below.

use crate::feat::brep_feat_builder::explorer;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    translate_curve2d, Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval,
    TrimmedCurve3,
};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use std::sync::Arc;

use super::loc_ope_wires_on_shape::{
    brep_tool_pnt, shape_key, top_exp_first_vertex, top_exp_last_vertex, top_exp_vertices,
};

pub(crate) type ShapeKey = (u64, u32);

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts (loc_ope_find_edges.rs / loc_ope_gluer.rs precedents).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Tolerance(shape) — vertex/edge/face tolerance.
pub(crate) fn brep_tool_tolerance(s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(E).
pub(crate) fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Curve(E, f, l) — the 3D curve and range (identity
/// location; loc_ope_find_edges.rs arch. diff. #1).
pub(crate) fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::Range(E, f, l) — the 3D range only.
pub(crate) fn brep_tool_range(edg: &Shape) -> (f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Surface(F) (loc_ope_gluer.rs arch. diff. #1).
pub(crate) fn brep_tool_surface(face: &Shape) -> Option<Surface3> {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(E, F, f, l) — the pcurve of the edge on
/// the face with its range (loc_ope_gluer.rs arch. diff. #1).
pub(crate) fn brep_tool_curve_on_surface(edg: &Shape, face: &Shape) -> Option<(Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .pcurves
            .get(&shape_key(face))
            .map(|(c, f, l)| (c.clone(), *f, *l)),
        _ => None,
    }
}

/// OCCT BRep_Tool::Parameter(V, E) — the stored vertex parameter on the edge.
pub(crate) fn brep_tool_parameter(vtx: &Shape, edg: &Shape) -> f64 {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .vertex_params
            .get(&vtx.ptr_id())
            .copied()
            .unwrap_or(0.0),
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Range(E, F, f, l) — the pcurve range of the edge on the
/// face.
pub(crate) fn brep_tool_range_on_face(edg: &Shape, face: &Shape) -> Option<(f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .pcurves
            .get(&shape_key(face))
            .map(|(_, f, l)| (*f, *l)),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Pure-math re-hosts (arch. diffs. #4/#6/#13).
// ---------------------------------------------------------------------------

/// OCCT Epsilon(theValue) (Standard_Real.hxx L238-247) — the distance to the
/// nearest representable double (pure-math re-host, arch. diff. #13).
fn epsilon(the_value: f64) -> f64 {
    if the_value == 0.0 {
        std::f64::EPSILON
    } else {
        std::f64::EPSILON * the_value.abs()
    }
}

/// OCCT ShapeAnalysis::AdjustByPeriod(Val, ToVal, Period)
/// (ShapeAnalysis.cxx L44-59).
fn shape_analysis_adjust_by_period(the_val: f64, to_val: f64, period: f64) -> f64 {
    let diff = the_val - to_val;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
}

/// OCCT ShapeAnalysis::AdjustToPeriod(Val, ValMin, ValMax)
/// (ShapeAnalysis.cxx L66-69).
fn shape_analysis_adjust_to_period(the_val: f64, val_min: f64, val_max: f64) -> f64 {
    shape_analysis_adjust_by_period(the_val, 0.5 * (val_min + val_max), val_max - val_min)
}

/// OCCT Geom_ConicSurface/Geom_SphericalSurface::UPeriod — the U period of
/// the elementary periodic surfaces (pure-math re-host, arch. diff. #4).
fn surface_u_period(surf: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match surf {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => TAU,
        _ => 0.0,
    }
}

/// OCCT Geom_SphericalSurface/ToroidalSurface::VPeriod (2*PI sphere,
/// 2*minorRadius torus; pure-math re-host, arch. diff. #4).
fn surface_v_period(surf: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match surf {
        Surface3::Sphere(_) => TAU,
        Surface3::Torus(t) => 2.0 * t.minor_radius,
        _ => 0.0,
    }
}

/// OCCT BndLib_Add2dCurve::Add(PCurve, Tol, Box) — fill a 2D box from the
/// pcurve range (architecture difference #7).
fn bnd_lib_add2d_curve(pc: &Curve2d, f: f64, l: f64, tol: f64, the_box: &mut BndBox2d) {
    // curve2d_bounding_box returns [u_min, u_max, v_min, v_max].
    let b = rcad_kernel::base::bnd_lib::curve2d_bounding_box(pc, f, l, tol);
    the_box.update(b[0], b[2], b[1], b[3]);
}

/// OCCT GeomAdaptor surface UPeriod wrapper (arch. diff. #4).
fn a_surf_u_period(s: &Surface3) -> f64 {
    surface_u_period(s)
}

/// OCCT GeomAdaptor surface VPeriod wrapper (arch. diff. #4).
fn a_surf_v_period(s: &Surface3) -> f64 {
    surface_v_period(s)
}

/// OCCT GeomAdaptor_Curve::Resolution(R3d) wrapper (arch. diff. #4).
fn curve_resolution_for(c: &Curve3, r3d: f64) -> f64 {
    rcad_kernel::topo::topods::curve_resolution(c, r3d)
}

/// OCCT Precision::IsInfinite.
fn is_infinite(v: f64) -> bool {
    v >= rcad_kernel::precision::INFINITE_VALUE || v <= -rcad_kernel::precision::INFINITE_VALUE
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
fn oriented(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT Geom2d_Curve::Translated(gp_Vec2d) — a translated copy.
fn translated_c2d(c: &Curve2d, dx: f64, dy: f64) -> Curve2d {
    translate_curve2d(c, DVec2::new(dx, dy))
}

/// OCCT BRepTools::UVBounds(F, Umin, Umax, Vmin, Vmax) — the union of the 2D
/// boxes of the face's edge pcurves (the AddUVBounds/BndLib_Add2dCurve
/// vehicle; returns [umin, umax, vmin, vmax]).
fn brep_tools_uv_bounds(face: &Shape) -> Option<[f64; 4]> {
    let mut a_box = BndBox2d::new();
    for edg in explorer(face, ShapeType::Edge, ShapeType::Shape) {
        if let Some((c2d, f, l)) = brep_tool_curve_on_surface(&edg, face) {
            bnd_lib_add2d_curve(&c2d, f, l, 0.0, &mut a_box);
        }
    }
    // -> (xmin, ymin, xmax, ymax) == (umin, vmin, umax, vmax).
    a_box.get().map(|g| [g.0, g.2, g.1, g.3])
}

// ---------------------------------------------------------------------------
// Re-hosted pending translations (arch. diffs. #9/#10/#11/#12).
// ---------------------------------------------------------------------------

/// OCCT ShapeConstruct_ProjectCurveOnSurface — pending TKShHealing
/// translation (architecture difference #9). Init stores the surface; Perform
/// fails (leaves the pcurve null) so the caller takes the OCCT null path.
pub(crate) struct ShapeConstructProjectCurveOnSurface {
    my_surf: Option<Surface3>, // OCCT: mySurf
    my_tol: f64,               // OCCT: myTol
}

impl Default for ShapeConstructProjectCurveOnSurface {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeConstructProjectCurveOnSurface {
    pub fn new() -> Self {
        ShapeConstructProjectCurveOnSurface {
            my_surf: None,
            my_tol: 0.0,
        }
    }

    /// OCCT ShapeConstruct_ProjectCurveOnSurface::Init(S, Tol2d).
    pub fn init(&mut self, the_s: Option<Surface3>, the_tol: f64) {
        self.my_surf = the_s;
        self.my_tol = the_tol;
    }

    /// OCCT Perform(C, First, Last, C2d, TolFirst, TolLast) — GAP (arch.
    /// diff. #9): the TKShHealing projection algorithm is not translated;
    /// C2d stays null (the OCCT failure output) and the caller returns.
    pub fn perform(
        &mut self,
        _the_c: &Curve3,
        _the_first: f64,
        _the_last: f64,
        the_c2d: &mut Option<Curve2d>,
        _the_tol_first: f64,
        _the_tol_last: f64,
    ) -> bool {
        *the_c2d = None;
        false
    }
}

/// OCCT GeomProjLib::Curve2d(C, S, Umin, Umax, Vmin, Vmax, Tol2d) — pending
/// TKTopAlgo translation (architecture difference #10): returns None (the
/// OCCT null handle).
fn geom_proj_lib_curve2d(
    _the_c: &Curve3,
    _the_s: &Surface3,
    _the_umin: f64,
    _the_umax: f64,
    _the_vmin: f64,
    _the_vmax: f64,
    _the_tol2d: f64,
) -> Option<Curve2d> {
    None
}

/// OCCT Extrema_ExtCC between two bounded 3D curves — pending TKGeomAlgo
/// translation (architecture difference #11). The OCCT accessor surface is
/// carried with "pending" outputs; IsDone() is false so the OCCT
/// !IsDone() continue path (cxx L1347-1350) is taken.
pub(crate) struct ExtremaExtCC;

impl ExtremaExtCC {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        _the_c1: &Curve3,
        _the_f1: f64,
        _the_l1: f64,
        _the_c2: &Curve3,
        _the_f2: f64,
        _the_l2: f64,
    ) -> Self {
        ExtremaExtCC
    }

    /// OCCT Extrema_ExtCC::IsDone().
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT Extrema_ExtCC::NbExt().
    pub fn nb_ext(&self) -> usize {
        0
    }

    /// OCCT Extrema_ExtCC::IsParallel().
    pub fn is_parallel(&self) -> bool {
        false
    }

    /// OCCT Extrema_ExtCC::SquareDistance(N).
    pub fn square_distance(&self, _the_n: usize) -> f64 {
        0.0
    }

    /// OCCT Extrema_ExtCC::TrimmedSquareDistances(D1, D2, D3, D4, P1, P2,
    /// P3, P4).
    pub fn trimmed_square_distances(
        &self,
        the_dists: &mut [f64; 4],
        _the_p11: &mut DVec3,
        _the_p12: &mut DVec3,
        _the_p21: &mut DVec3,
        _the_p22: &mut DVec3,
    ) {
        *the_dists = [0.0; 4];
    }

    /// OCCT Extrema_ExtCC::Points(N, P1, P2).
    pub fn points(&self, _the_n: usize, the_p1: &mut (f64, DVec3), the_p2: &mut (f64, DVec3)) {
        the_p1.0 = 0.0;
        the_p2.0 = 0.0;
    }
}

/// OCCT ShapeAnalysis_Edge::CheckSameParameter(edge, face, maxdev)
/// (ShapeAnalysis_Edge.cxx L704-833) — reduced re-host (architecture
/// difference #12): the deviation is the max sampled distance between the 3D
/// curve and the pcurve-lifted point (the BRep::same_parameter reduction of
/// BRepLib_ValidateEdge). Returns true when maxdev exceeds the edge
/// tolerance (the DONE1 status OCCT reports).
fn shape_analysis_edge_check_same_parameter(
    the_edge: &Shape,
    the_face: &Shape,
    the_maxdev: &mut f64,
) -> bool {
    // OCCT L710-713: degenerated edges are skipped.
    if brep_tool_degenerated(the_edge) {
        *the_maxdev = 0.0;
        return false;
    }
    *the_maxdev = 0.0;
    // OCCT L722-729: the 3D curve.
    let Some((c3d, f, l)) = brep_tool_curve(the_edge) else {
        return false;
    };
    // OCCT L744-...: the pcurve on the given face (the CurveOnSurface
    // representation matching it).
    let Some((pc, _, _)) = brep_tool_curve_on_surface(the_edge, the_face) else {
        return false;
    };
    let Some(surf) = brep_tool_surface(the_face) else {
        return false;
    };
    // Reduced BRepLib_ValidateEdge walk: sample both curves over the range.
    const NB_SAMPLES: usize = 23;
    for i in 0..NB_SAMPLES {
        let u = f + (l - f) * (i as f64) / (NB_SAMPLES - 1) as f64;
        let p3d = c3d.point_at(u);
        let uv = pc.point_at(u);
        let ps = surf.point_at(uv.x, uv.y);
        let dev = (p3d - ps).length();
        if dev > *the_maxdev {
            *the_maxdev = dev;
        }
    }
    // OCCT L833: status DONE1 when maxdev > tolerance.
    *the_maxdev > brep_tool_tolerance(the_edge)
}

// ---------------------------------------------------------------------------
// BRepAdaptor_Curve2d re-host (architecture difference #3).
// ---------------------------------------------------------------------------

/// OCCT BRepAdaptor_Curve2d(edg, fac) — the edge pcurve on the face with its
/// range; Value(u) evaluates the pcurve.
pub(crate) struct BRepAdaptorCurve2d {
    my_pc: Option<(Curve2d, f64, f64)>,
}

impl BRepAdaptorCurve2d {
    /// OCCT BRepAdaptor_Curve2d(E, F).
    pub fn new(the_edg: &Shape, the_face: &Shape) -> Self {
        BRepAdaptorCurve2d {
            my_pc: brep_tool_curve_on_surface(the_edg, the_face),
        }
    }

    /// OCCT Curve().IsNull().
    pub fn is_null(&self) -> bool {
        self.my_pc.is_none()
    }

    /// OCCT Value(u).
    pub fn value(&self, the_u: f64) -> DVec2 {
        match &self.my_pc {
            Some((c, _, _)) => c.point_at(the_u),
            None => DVec2::ZERO,
        }
    }
}

// ---------------------------------------------------------------------------
// BRep_Builder mutation vehicles (architecture difference #2): in-place
// Arc::make_mut edits on the Shape payload.
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::UpdateVertex(V, Tol) — BRep_TVertex::UpdateTolerance
/// keeps the max.
fn builder_update_vertex(the_v: &mut Shape, the_tol: f64) {
    if let TShape::Vertex(vd) = Arc::make_mut(&mut the_v.data) {
        vd.tolerance = vd.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, Tol).
fn builder_update_edge_tolerance(the_e: &mut Shape, the_tol: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, F, Tol) — bind the pcurve on the
/// face. The pcurve range is the edge 3D range (OCCT BRep_Tool::CurveOnSurface
/// returns the edge's First/Last for the pcurve bounds).
fn builder_update_edge_pcurve(the_e: &mut Shape, the_c2d: &Curve2d, the_f: &Shape, the_tol: f64) {
    let key = shape_key(the_f);
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        let (f0, l0) = (ed.range[0], ed.range[1]);
        ed.pcurves.insert(key, (the_c2d.clone(), f0, l0));
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C2df, C2dr, F, Tol) — the seam variant
/// (two pcurves on the closed surface).
fn builder_update_edge_pcurves_seam(
    the_e: &mut Shape,
    the_c2d_f: &Curve2d,
    the_c2d_r: &Curve2d,
    the_f: &Shape,
    the_tol: f64,
) {
    let key = shape_key(the_f);
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        let (f0, l0) = (ed.range[0], ed.range[1]);
        ed.pcurves.insert(key, (the_c2d_f.clone(), f0, l0));
        ed.representations
            .push(rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
                face: key,
                pcurve1: the_c2d_f.clone(),
                pcurve2: the_c2d_r.clone(),
                range: [f0, l0],
            });
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::Range(E, First, Last).
fn builder_range(the_e: &mut Shape, the_first: f64, the_last: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.range = [the_first, the_last];
    }
}

/// OCCT BRep_Builder::Add(E, V) — attach the vertex by orientation
/// (FORWARD -> first, REVERSED -> last).
fn builder_add_edge_vertex(the_e: &mut Shape, the_v: &Shape) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        match the_v.orientation {
            Orientation::Reversed => ed.last = the_v.clone(),
            _ => ed.first = the_v.clone(),
        }
    }
}

/// OCCT BRep_Builder::SameParameter(E, flag) — the TEdgeData flag write
/// (cxx L907).
fn builder_set_same_parameter(the_e: &mut Shape, flag: bool) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.same_parameter = flag;
    }
}

// ---------------------------------------------------------------------------
// Statics of LocOpe_WiresOnShape.cxx.
// ---------------------------------------------------------------------------
// NOTE: the "value assigned is never read" warnings inside the statics below
// are intentional — they mirror the OCCT reference flow (cxx L589 aPcur is
// consumed only on the periodic paths; the (f, l) pairs of L744-765 /
// L1107-1142 are overwritten by the following fetches).

/// OCCT static Project(V, p2d, F, theEdge, param) (cxx L491-675).
#[allow(unused_assignments)]
pub(crate) fn project_vertex_p2d_face_edge(
    the_v: &mut Shape,
    p2d: DVec2,
    the_f: &Shape,
    the_edge: &mut Shape,
    param: &mut f64,
) -> bool {
    // OCCT L497-498.
    let a_tol_v = brep_tool_tolerance(the_v);
    let mut dmin = if the_edge.is_null() {
        f64::MAX
    } else {
        a_tol_v * a_tol_v
    };

    let mut valret = false; // OCCT L500

    // OCCT L502.
    let Some(a_surf) = brep_tool_surface(the_f) else {
        return valret;
    };

    if the_edge.is_null() {
        // OCCT L506.
        let toproj = brep_tool_pnt(the_v).unwrap_or(DVec3::ZERO);
        // OCCT L507: exp(F.Oriented(TopAbs_FORWARD), TopAbs_EDGE).
        for edg in explorer(
            &oriented(the_f, Orientation::Forward),
            ShapeType::Edge,
            ShapeType::Shape,
        ) {
            // OCCT L510-513.
            let mut a_cur_dist = rcad_kernel::precision::INFINITE_VALUE;
            let mut a_cur_par = rcad_kernel::precision::INFINITE_VALUE;
            match brep_tool_curve(&edg) {
                Some((c, f, l)) => {
                    // OCCT L514-524 (the !C.IsNull() branch).
                    a_cur_par = project_vertex_edge(the_v, &edg);
                    if is_infinite(a_cur_par) {
                        continue;
                    }
                    let a_cur_p_bound = c.point_at(a_cur_par);
                    a_cur_dist = (a_cur_p_bound - toproj).length_squared();
                    let _ = (f, l);
                }
                None => {
                    // OCCT L526-546 (the C.IsNull() branch).
                    if !is_infinite(p2d.x) {
                        let Some((_, f2, l2)) = brep_tool_curve_on_surface(&edg, the_f) else {
                            continue;
                        };
                        let _ = (f2, l2);
                        a_cur_par = project_vertex_p2d_edge_face(the_v, p2d, &edg, the_f);
                        if is_infinite(a_cur_par) {
                            continue;
                        }
                        // OCCT L540-545 (PC re-fetched).
                        let Some((pc, _, _)) = brep_tool_curve_on_surface(&edg, the_f) else {
                            continue;
                        };
                        let a_p_proj = pc.point_at(a_cur_par);
                        let a_cur_p_bound = a_surf.point_at(a_p_proj.x, a_p_proj.y);
                        a_cur_dist = (a_cur_p_bound - toproj).length_squared();
                    }
                }
            }

            // OCCT L548-554.
            if a_cur_dist < dmin {
                *the_edge = edg.clone();
                the_edge.orientation = edg.orientation;
                dmin = a_cur_dist;
                *param = a_cur_par;
            }
        }
        if the_edge.is_null() {
            // OCCT L556-559.
            return false;
        }
    } else if is_infinite(*param) {
        // OCCT L561-566.
        match brep_tool_curve(the_edge) {
            Some(_) => *param = project_vertex_edge(the_v, the_edge),
            None => *param = project_vertex_p2d_edge_face(the_v, p2d, the_edge, the_f),
        }
    }

    // OCCT L568.
    let ttol = a_tol_v + brep_tool_tolerance(the_edge);
    if dmin <= ttol * ttol {
        // OCCT L571.
        valret = true;
        // OCCT L572-575: GeomAdaptor_Surface over aSurf (arch. diff. #4);
        // NOTE: OCCT L575 computes aVResolution with UResolution(1.) —
        // kept verbatim.
        let an_u_resolution = rcad_kernel::topo::topods::u_resolution_for_surface(&a_surf, 1.0);
        let a_v_resolution = rcad_kernel::topo::topods::u_resolution_for_surface(&a_surf, 1.0);

        // OCCT L578.
        let a_crv_bound = brep_tool_curve_on_surface(the_edge, the_f);
        if let Some((a_crv_bound_c, _, _)) = a_crv_bound {
            // OCCT L582.
            let mut a_p_bound2d = a_crv_bound_c.point_at(*param);

            // OCCT L584-587: distance in 2D space recomputed in the 3D space
            // in order to tolerance of vertex cover gap in 2D space. For
            // consistency with the check of the validity in the
            // BRepCheck_Wire.
            let dom = a_surf.default_domain();
            let mut dumax = 0.01 * (dom[1] - dom[0]);
            let mut dvmax = 0.01 * (dom[3] - dom[2]);

            let mut a_pcur = p2d; // OCCT L589
            let mut dumin = (a_pcur.x - a_p_bound2d.x).abs(); // OCCT L590
            let mut dvmin = (a_pcur.y - a_p_bound2d.y).abs(); // OCCT L591
            if dumin > dumax && a_surf.is_u_periodic() {
                // OCCT L592-610.
                let mut a_x1 = a_p_bound2d.x;
                let mut a_shift = shape_analysis_adjust_to_period(a_x1, dom[0], dom[1]);
                a_x1 += a_shift;
                a_p_bound2d.x = a_x1;
                let mut a_x2 = a_pcur.x;
                a_shift = shape_analysis_adjust_to_period(a_x2, dom[0], dom[1]);
                a_x2 += a_shift;
                dumin = (a_x2 - a_x1).abs();
                if dumin > dumax
                    && ((dumin - a_surf_u_period(&a_surf)).abs()
                        < rcad_kernel::precision::PCONFUSION)
                {
                    a_x2 = a_x1;
                    dumin = 0.0;
                }
                a_pcur.x = a_x2;
            }

            if dvmin > dvmax && a_surf.is_v_periodic() {
                // OCCT L612-630.
                let mut a_y1 = a_p_bound2d.y;
                let mut a_shift = shape_analysis_adjust_to_period(a_y1, dom[2], dom[3]);
                a_y1 += a_shift;
                a_p_bound2d.y = a_y1;
                let mut a_y2 = a_pcur.y;
                a_shift = shape_analysis_adjust_to_period(a_y2, dom[2], dom[3]);
                a_y2 += a_shift;
                dvmin = (a_y1 - a_y2).abs();
                if dvmin > dvmax
                    && ((dvmin - a_surf_v_period(&a_surf)).abs()
                        < rcad_kernel::precision::CONFUSION)
                {
                    a_y2 = a_y1;
                    dvmin = 0.0;
                }
                a_pcur.y = a_y2;
            }
            // OCCT L631.
            let mut a_dist3d = a_tol_v;
            if dumin > dumax || dvmin > dvmax {
                // OCCT L635-647.
                dumax = rcad_kernel::topo::topods::u_resolution_for_surface(&a_surf, a_tol_v);
                dvmax = rcad_kernel::topo::topods::v_resolution_for_surface(&a_surf, a_tol_v);
                let a_tol2d = 2.0 * dumax.max(dvmax);
                let a_dist2d = dumin.max(dvmin);

                if a_dist2d > a_tol2d {
                    let a_dist3d1 = a_dist2d / an_u_resolution.max(a_v_resolution);
                    if a_dist3d1 > a_dist3d {
                        a_dist3d = a_dist3d1;
                    }
                }
            }

            // OCCT L650-656: added check by 3D the same as in the
            // BRepCheck_Wire::SelfIntersect.
            let a_p_bound = a_surf.point_at(a_p_bound2d.x, a_p_bound2d.y);
            let a_pv2d = a_surf.point_at(p2d.x, p2d.y);
            let a_dist_points_3d = (a_pv2d - a_p_bound).length_squared();
            let a_max_dist = a_dist_points_3d.max(a_dist3d * a_dist3d);

            // OCCT L658-663.
            if a_tol_v * a_tol_v < a_max_dist {
                let a_new_tol = a_max_dist.sqrt();
                builder_update_vertex(the_v, a_new_tol);
            }
        }
    }
    // OCCT L666-673 (OCCT_DEBUG_MESH) not translated (arch. diff. #15).
    // OCCT L674.
    valret
}

/// OCCT static Project(V, theEdge) (cxx L679-697).
pub(crate) fn project_vertex_edge(the_v: &Shape, the_edge: &Shape) -> f64 {
    // OCCT L685-688.
    let toproj = brep_tool_pnt(the_v).unwrap_or(DVec3::ZERO);
    let Some((c, f, l)) = brep_tool_curve(the_edge) else {
        // OCCT inits the projector with a null curve; the rcad output is the
        // OCCT Precision::Infinite() branch value.
        return rcad_kernel::precision::INFINITE_VALUE;
    };
    // OCCT L689-694 (the Loc.IsIdentity() transform branch applies only for
    // identity locations — loc_ope_find_edges.rs arch. diff. #1).
    let _ = &c;
    // OCCT L696: proj.NbPoints() > 0 ? proj.LowerDistanceParameter() :
    // Precision::Infinite().
    let proj =
        rcad_kernel::base::geom_api::project::closest_point_on_curve_range(&c, toproj, f, l, 64);
    proj.param
}

/// OCCT static Project(V, p2d, theEdge, theFace) (cxx L701-718).
pub(crate) fn project_vertex_p2d_edge_face(
    the_v: &Shape,
    p2d: DVec2,
    the_edge: &Shape,
    the_face: &Shape,
) -> f64 {
    let _ = the_v;
    // OCCT L713: PC = BRep_Tool::CurveOnSurface(theEdge, theFace, f, l).
    let Some((pc, f, l)) = brep_tool_curve_on_surface(the_edge, the_face) else {
        // OCCT inits the 2D projector with a null curve; the rcad output is
        // the OCCT Precision::Infinite() branch value.
        return rcad_kernel::precision::INFINITE_VALUE;
    };
    // OCCT L715-717: Geom2dAPI_ProjectPointOnCurve -> ExtPC2d (arch.
    // diff. #5).
    let proj = rcad_kernel::base::extrema::ExtPC2d::new(
        p2d,
        &pc,
        rcad_kernel::precision::CONFUSION,
        f,
        l,
    );
    if proj.is_done() && proj.nb_ext() > 0 {
        proj.point(1).param
    } else {
        rcad_kernel::precision::INFINITE_VALUE
    }
}

/// OCCT static PutPCurve(Edg, Fac) (cxx L722-909).
#[allow(unused_assignments)]
pub(crate) fn put_pcurve(the_edg: &mut Shape, the_fac: &Shape) {
    // OCCT L727-734: S = BRep_Tool::Surface(Fac, LocFac); styp with the
    // rectangular-trimmed basis strip.
    let Some(mut s) = brep_tool_surface(the_fac) else {
        // OCCT would dereference a null surface; rcad returns (marked).
        return;
    };
    if matches!(s, Surface3::Trimmed(_)) {
        let (basis, _) = rcad_kernel::topo::topods::surface_adaptor_basis_and_bounds(&s);
        s = basis.clone();
    }

    // OCCT L736-739.
    if matches!(s, Surface3::Plane(_)) {
        return;
    }

    // OCCT L741-742: BRepTools::UVBounds(Fac, Umin, Umax, Vmin, Vmax).
    let Some(bounds) = brep_tools_uv_bounds(the_fac) else {
        return;
    };
    let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);

    // OCCT L744-765.
    let mut f;
    let mut l;
    if let Some((c2d, pf, pl)) = brep_tool_curve_on_surface(the_edg, the_fac) {
        f = pf;
        l = pl;
        // OCCT L749-764: p2d at mid-parameter; the IsIn check.
        let p2d = c2d.point_at((f + l) * 0.5);
        let mut is_in = true;
        if p2d.x < umin - rcad_kernel::precision::PCONFUSION
            || p2d.x > umax + rcad_kernel::precision::PCONFUSION
        {
            is_in = false;
        }
        if p2d.y < vmin - rcad_kernel::precision::PCONFUSION
            || p2d.y > vmax + rcad_kernel::precision::PCONFUSION
        {
            is_in = false;
        }
        if is_in {
            return;
        }
    }

    // OCCT L767-773: C = BRep_Tool::Curve(Edg, Loc, f, l) (+transform).
    let Some((mut c, cf, cl)) = brep_tool_curve(the_edg) else {
        // OCCT would dereference a null curve (C->DynamicType()); rcad
        // returns (marked).
        return;
    };
    f = cf;
    l = cl;

    // OCCT L775-778: wrap into a trimmed curve.
    if !matches!(c, Curve3::Trimmed(_)) {
        c = Curve3::Trimmed(TrimmedCurve3::new(c, f, l));
    }

    // OCCT L780: S = BRep_Tool::Surface(Fac).
    let Some(s) = brep_tool_surface(the_fac) else {
        return;
    };

    // OCCT L782-792: TolFirst/TolLast from the edge vertices.
    let mut tol_first = -1.0;
    let mut tol_last = -1.0;
    let (mut v1, mut v2) = top_exp_vertices(the_edg);
    if let Some(vv) = &v1 {
        tol_first = brep_tool_tolerance(vv);
    }
    if let Some(vv) = &v2 {
        tol_last = brep_tool_tolerance(vv);
    }

    // OCCT L794-799: the projection tool (architecture difference #9).
    let tol2d = rcad_kernel::precision::CONFUSION;
    let mut c2d_opt: Option<Curve2d> = None;
    let mut a_tool_proj = ShapeConstructProjectCurveOnSurface::new();
    a_tool_proj.init(Some(s.clone()), tol2d);
    a_tool_proj.perform(&c, f, l, &mut c2d_opt, tol_first, tol_last);
    let Some(mut c2d) = c2d_opt else {
        // OCCT L800-803: the null-curve return.
        return;
    };

    // OCCT L805-827: pf/pl, PF/PL and the vertex bookkeeping.
    let pf = c2d.point_at(f);
    let pl = c2d.point_at(l);
    let p_f = s.point_at(pf.x, pf.y);
    let p_l = s.point_at(pl.x, pl.y);
    if the_edg.orientation == Orientation::Reversed {
        v1 = top_exp_last_vertex(the_edg);
        if let Some(ref mut vv) = v1 {
            vv.orientation = Orientation::Reversed.compose(vv.orientation);
        }
    } else {
        v1 = top_exp_first_vertex(the_edg);
    }
    if the_edg.orientation == Orientation::Reversed {
        v2 = top_exp_first_vertex(the_edg);
        if let Some(ref mut vv) = v2 {
            vv.orientation = Orientation::Reversed.compose(vv.orientation);
        }
    } else {
        v2 = top_exp_last_vertex(the_edg);
    }

    // OCCT L829-840: handling of internal vertices (the OCCT L829 test reads
    // `!V1.IsNull() && V2.IsNull()` and then reads both V1 and V2 —
    // reproduced verbatim).
    if v1.is_some() && v2.is_none() {
        let old1 = v1.as_ref().map(brep_tool_tolerance).unwrap_or(0.0);
        let old2 = v2.as_ref().map(brep_tool_tolerance).unwrap_or(0.0);
        let pnt1 = v1.as_ref().and_then(brep_tool_pnt).unwrap_or(DVec3::ZERO);
        let pnt2 = v2.as_ref().and_then(brep_tool_pnt).unwrap_or(DVec3::ZERO);
        let tol1 = (pnt1 - p_f).length();
        let tol2 = (pnt2 - p_l).length();
        if let Some(vv) = &mut v1 {
            builder_update_vertex(vv, old1.max(tol1));
        }
        if let Some(vv) = &mut v2 {
            builder_update_vertex(vv, old2.max(tol2));
        }
    }

    // OCCT L842-871: the U-periodic shift.
    if s.is_u_periodic() {
        let up = a_surf_u_period(&s);
        let tolu = rcad_kernel::precision::PCONFUSION;
        let mut nbtra: i32 = 0;
        let mut the_umin = pf.x.min(pl.x);
        let mut the_umax = pf.x.max(pl.x);

        if the_umin < umin - tolu {
            while the_umin < umin - tolu {
                the_umin += up;
                nbtra += 1;
            }
        } else if the_umax > umax + tolu {
            while the_umax > umax + tolu {
                the_umax -= up;
                nbtra -= 1;
            }
        }

        if nbtra != 0 {
            c2d = translated_c2d(&c2d, nbtra as f64 * up, 0.0);
        }
    }

    // OCCT L873-904: the V-periodic shift.
    if s.is_v_periodic() {
        let vp = a_surf_v_period(&s);
        let tolv = rcad_kernel::precision::PCONFUSION;
        let mut nbtra: i32 = 0;
        let mut the_vmin = pf.y.min(pl.y);
        let mut the_vmax = pf.y.max(pl.y);

        if the_vmin < vmin - tolv {
            while the_vmin < vmin - tolv {
                the_vmin += vp;
                the_vmax += vp;
                nbtra += 1;
            }
        } else if the_vmax > vmax + tolv {
            while the_vmax > vmax + tolv {
                the_vmax -= vp;
                the_vmin -= vp;
                nbtra -= 1;
            }
        }

        if nbtra != 0 {
            c2d = translated_c2d(&c2d, 0.0, nbtra as f64 * vp);
        }
    }

    // OCCT L905-908.
    let tol_edg = brep_tool_tolerance(the_edg);
    builder_update_edge_pcurve(the_edg, &c2d, the_fac, tol_edg);
    // OCCT L907: B.SameParameter(Edg, false) — the TEdgeData flag.
    builder_set_same_parameter(the_edg, false);
    // OCCT L908: BRepLib::SameParameter(Edg, tol2d).
    crate::topalgo::brep_lib::brep_lib::BRepLib::same_parameter(the_edg, tol2d);
}

/// OCCT static PutPCurves(Efrom, Eto, myShape) (cxx L913-1287).
#[allow(unused_assignments)]
pub(crate) fn put_pcurves(the_efrom: &mut Shape, the_eto: &Shape, my_shape: &Shape) {
    // OCCT L916-928: collect the faces carrying Eto.
    let mut lfaces: Vec<Shape> = Vec::new();
    for exp in explorer(my_shape, ShapeType::Face, ShapeType::Shape) {
        for exp2 in explorer(&exp, ShapeType::Edge, ShapeType::Shape) {
            if exp2.is_same(the_eto) {
                lfaces.push(exp.clone());
            }
        }
    }

    // OCCT L930-933.
    if lfaces.len() != 1 && lfaces.len() != 2 {
        // OCCT: throw Standard_ConstructionError().
        panic!("Standard_ConstructionError");
    }

    // soit bord libre, soit connexite entre 2 faces, eventuellement edge
    // closed (OCCT L935).

    // OCCT L937-940.
    if lfaces.len() == 1 {
        return; // sera fait par PutPCurve.... on l`espere
    }

    // OCCT L950-1097: two DIFFERENT faces.
    let first_same_last = lfaces
        .first()
        .expect("Lfaces")
        .is_same(lfaces.last().expect("Lfaces"));
    if !first_same_last {
        for itl in lfaces.iter() {
            let fac = itl;
            // OCCT L957-960.
            if brep_tool_curve_on_surface(the_efrom, fac).is_some() {
                continue;
            }
            // OCCT L961-971: S/styp with the trimmed-surface strip.
            let Some(s0) = brep_tool_surface(fac) else {
                return;
            };
            let mut s = s0;
            if matches!(s, Surface3::Trimmed(_)) {
                let (basis, _) =
                    rcad_kernel::topo::topods::surface_adaptor_basis_and_bounds(&s);
                s = basis.clone();
            }
            if matches!(s, Surface3::Plane(_)) {
                continue;
            }

            // OCCT L973.
            let Some(bounds) = brep_tools_uv_bounds(fac) else {
                return;
            };
            let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);

            // OCCT L974-979: C with transform.
            let Some((mut c, cf, cl)) = brep_tool_curve(the_efrom) else {
                return;
            };
            let mut f = cf;
            let mut l = cl;

            // OCCT L981-984.
            if !matches!(c, Curve3::Trimmed(_)) {
                c = Curve3::Trimmed(TrimmedCurve3::new(c, f, l));
            }

            // OCCT L986.
            let Some(s) = brep_tool_surface(fac) else {
                return;
            };

            // OCCT L988-993: Compute the tol2d (arch. diff. #4).
            let tol3d = brep_tool_tolerance(the_efrom).max(brep_tool_tolerance(fac));
            let tol_u = rcad_kernel::topo::topods::u_resolution_for_surface(&s, tol3d);
            let tol_v = rcad_kernel::topo::topods::v_resolution_for_surface(&s, tol3d);
            let tol2d = tol_u.max(tol_v);

            // OCCT L995-999: GeomProjLib::Curve2d (architecture
            // difference #10 — pending; the OCCT null-curve return).
            let Some(mut c2d) = geom_proj_lib_curve2d(&c, &s, umin, umax, vmin, vmax, tol2d)
            else {
                return;
            };

            // OCCT L1001-1002.
            let pf = c2d.point_at(f);
            let pl = c2d.point_at(l);

            // OCCT L1004-1048: the U-periodic shift.
            if s.is_u_periodic() {
                let up = a_surf_u_period(&s);
                let tolu = rcad_kernel::precision::PCONFUSION;
                let mut nbtra: i32 = 0;
                let mut the_umin = pf.x.min(pl.x);
                let mut the_umax = pf.x.max(pl.x);

                if the_umin < umin - tolu {
                    while the_umin < umin - tolu {
                        the_umin += up;
                        the_umax += up;
                        nbtra += 1;
                    }
                } else if the_umax > umax + tolu {
                    while the_umax > umax + tolu {
                        the_umax -= up;
                        the_umin -= up;
                        nbtra -= 1;
                    }
                }
                // (the OCCT L1030-1043 commented block is not translated)
                if nbtra != 0 {
                    c2d = translated_c2d(&c2d, nbtra as f64 * up, 0.0);
                }
            }

            // OCCT L1050-1094: the V-periodic shift.
            if s.is_v_periodic() {
                let vp = a_surf_v_period(&s);
                let tolv = rcad_kernel::precision::PCONFUSION;
                let mut nbtra: i32 = 0;
                let mut the_vmin = pf.y.min(pl.y);
                let mut the_vmax = pf.y.max(pl.y);

                if the_vmin < vmin - tolv {
                    while the_vmin < vmin - tolv {
                        the_vmin += vp;
                        the_vmax += vp;
                        nbtra += 1;
                    }
                } else if the_vmax > vmax + tolv {
                    while the_vmax > vmax + tolv {
                        the_vmax -= vp;
                        the_vmin -= vp;
                        nbtra -= 1;
                    }
                }
                // (the OCCT L1076-1089 commented block is not translated)
                if nbtra != 0 {
                    c2d = translated_c2d(&c2d, 0.0, nbtra as f64 * vp);
                }
            }

            // OCCT L1095.
            let tol_efrom = brep_tool_tolerance(the_efrom);
            builder_update_edge_pcurve(the_efrom, &c2d, fac, tol_efrom);
            let _ = (&mut f, &mut l);
        }
    } else {
        // OCCT L1099-1286: the seam (closed-edge) branch.
        let fac = lfaces.first().expect("Lfaces").clone();
        // OCCT L1102-1105.
        if !brep_tool_is_closed_on_face(the_eto, &fac) {
            // OCCT: throw Standard_ConstructionError().
            panic!("Standard_ConstructionError");
        }

        // OCCT L1107-1142: the four pcurve fetches; every fetch overwrites
        // the shared (f, l) exactly like the OCCT reference-passing.
        let local_e_f = oriented(the_efrom, Orientation::Forward);
        let local_f_f = oriented(&fac, Orientation::Forward);
        let mut f = 0.0;
        let mut l = 0.0;

        let c2dff = brep_tool_curve_on_surface(&local_e_f, &local_f_f);
        let c2dfr = brep_tool_curve_on_surface(
            &oriented(the_efrom, Orientation::Reversed),
            &local_f_f,
        );
        let c2dtf =
            brep_tool_curve_on_surface(&oriented(the_eto, Orientation::Forward), &local_f_f);
        let c2dtr =
            brep_tool_curve_on_surface(&oriented(the_eto, Orientation::Reversed), &local_f_f);
        if let Some((_, ff, ll)) = &c2dtr {
            f = *ff;
            l = *ll;
        }

        // OCCT L1144-1145: ptf (sur courbe frw), ptr (sur courbe rev).
        let Some((c2dtf_c, _, _)) = &c2dtf else {
            // OCCT dereferences c2dtf without a check; rcad returns (marked).
            return;
        };
        let Some((c2dtr_c, _, _)) = &c2dtr else {
            return;
        };
        let ptf = c2dtf_c.point_at(f);
        let ptr = c2dtr_c.point_at(f);

        // OCCT L1147: bool isoU = (std::abs(ptf.Y() - ptr.Y()) <
        // Epsilon(ptf.X())) — meme V.
        let iso_u = (ptf.y - ptr.y).abs() < epsilon(ptf.x);

        // Efrom et Eto dans le meme sens??? (OCCT L1149)

        // OCCT L1151-1156: C from Efrom (+transform); C->D1(f, pt, d1f).
        let Some((c, cf, cl)) = brep_tool_curve(the_efrom) else {
            return;
        };
        f = cf;
        l = cl;
        let d1f = c.derivative_at(f);
        let _ = c.point_at(f);

        // OCCT L1163-1166.
        let first_vertex = top_exp_first_vertex(the_efrom);
        let vtx_param = first_vertex
            .as_ref()
            .map(|v| brep_tool_parameter(v, the_efrom))
            .unwrap_or(0.0);
        let ba_curve2d = BRepAdaptorCurve2d::new(the_efrom, &fac);
        let p2d = ba_curve2d.value(vtx_param);

        // OCCT L1168.
        let prmproj = match &first_vertex {
            Some(fv) => project_vertex_p2d_edge_face(fv, p2d, the_eto, &fac),
            None => rcad_kernel::precision::INFINITE_VALUE,
        };

        // OCCT L1170-1177: C from Eto (+transform); C->D1(prmproj, pt, d1t).
        let Some((c, _, _)) = brep_tool_curve(the_eto) else {
            return;
        };
        let d1t = c.derivative_at(prmproj);

        // OCCT L1179.
        let same_ori = d1t.dot(d1f) > 0.0;

        // OCCT L1181-1223: fill the missing pcurves.
        let pair: (Option<Curve2d>, Option<Curve2d>) = match (&c2dff, &c2dfr) {
            (None, None) => {
                // OCCT L1183-1189: S/styp.
                let Some(s0) = brep_tool_surface(&fac) else {
                    return;
                };
                let mut s = s0;
                if matches!(s, Surface3::Trimmed(_)) {
                    let (basis, _) =
                        rcad_kernel::topo::topods::surface_adaptor_basis_and_bounds(&s);
                    s = basis.clone();
                }

                // OCCT L1191-1196: C from Efrom (+transform, trimmed wrap).
                let Some((mut c, cf1, cl1)) = brep_tool_curve(the_efrom) else {
                    return;
                };
                f = cf1;
                l = cl1;
                if !matches!(c, Curve3::Trimmed(_)) {
                    c = Curve3::Trimmed(TrimmedCurve3::new(c, f, l));
                }

                // OCCT L1204-1210: Compute the tol2d.
                let Some(bounds) = brep_tools_uv_bounds(&fac) else {
                    return;
                };
                let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                let tol3d = brep_tool_tolerance(the_efrom).max(brep_tool_tolerance(&fac));
                let tol_u = rcad_kernel::topo::topods::u_resolution_for_surface(&s, tol3d);
                let tol_v = rcad_kernel::topo::topods::v_resolution_for_surface(&s, tol3d);
                let tol2d = tol_u.max(tol_v);

                // OCCT L1212-1214 (architecture difference #10).
                match geom_proj_lib_curve2d(&c, &s, umin, umax, vmin, vmax, tol2d) {
                    Some(c2d) => (Some(c2d.clone()), Some(c2d)),
                    None => {
                        // OCCT L1227 dereferences the null C2d; rcad returns
                        // instead of crashing (marked, arch. diff. #10).
                        return;
                    }
                }
            }
            (Some(ff), None) => {
                // OCCT L1216-1219: c2dfr = c2dff.
                (Some(ff.0.clone()), Some(ff.0.clone()))
            }
            (None, Some(fr)) => {
                // OCCT L1220-1223: c2dff = c2dfr.
                (Some(fr.0.clone()), Some(fr.0.clone()))
            }
            (Some(ff), Some(fr)) => (Some(ff.0.clone()), Some(fr.0.clone())),
        };
        let (Some(mut c2dff_c), Some(mut c2dfr_c)) = pair else {
            return;
        };

        // OCCT L1225: BRep_Tool::Range(Efrom, f, l).
        let (rf, rl) = brep_tool_range(the_efrom);
        f = rf;
        l = rl;

        // OCCT L1227-1228.
        let p2f = c2dff_c.point_at(f);
        let p2r = c2dfr_c.point_at(f);

        if iso_u {
            if same_ori {
                // OCCT L1232-1241.
                if (ptf.x - p2f.x).abs() > epsilon(ptf.x) {
                    c2dff_c = translated_c2d(&c2dff_c, ptf.x - p2f.x, 0.0);
                }
                if (ptr.x - p2r.x).abs() > epsilon(ptr.x) {
                    c2dfr_c = translated_c2d(&c2dfr_c, ptr.x - p2r.x, 0.0);
                }
            } else {
                // OCCT L1243-1254.
                if (ptr.x - p2f.x).abs() > epsilon(ptr.x) {
                    c2dff_c = translated_c2d(&c2dff_c, ptr.x - p2f.x, 0.0);
                }
                if (ptf.x - p2r.x).abs() > epsilon(ptf.x) {
                    c2dfr_c = translated_c2d(&c2dfr_c, ptf.x - p2r.x, 0.0);
                }
            }
            // on est bien en U, recalage si periodique en V a faire
        } else {
            // !isoU soit isoV
            if same_ori {
                // OCCT L1261-1270.
                if (ptf.y - p2f.y).abs() > epsilon(ptf.y) {
                    c2dff_c = translated_c2d(&c2dff_c, 0.0, ptf.y - p2f.y);
                }
                if (ptr.y - p2r.y).abs() > epsilon(ptr.y) {
                    c2dfr_c = translated_c2d(&c2dfr_c, 0.0, ptr.y - p2r.y);
                }
            } else {
                // OCCT L1272-1282.
                if (ptr.y - p2f.y).abs() > epsilon(ptr.y) {
                    c2dff_c = translated_c2d(&c2dff_c, 0.0, ptr.y - p2f.y);
                }
                if (ptf.y - p2r.y).abs() > epsilon(ptf.y) {
                    c2dfr_c = translated_c2d(&c2dfr_c, 0.0, ptf.y - p2r.y);
                }
            }
            // on est bien en V, recalage si periodique en U a faire
        }

        // OCCT L1285.
        let tol_efrom = brep_tool_tolerance(the_efrom);
        builder_update_edge_pcurves_seam(the_efrom, &c2dff_c, &c2dfr_c, &fac, tol_efrom);
        let _ = (&mut f, &mut l);
    }
}

/// OCCT BRep_Tool::IsClosed(E, F) — the edge is a seam (closed) on the face
/// (the BRep_CurveOnClosedSurface representation matching the face).
fn brep_tool_is_closed_on_face(edg: &Shape, face: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let fkey = shape_key(face);
            ed.representations.iter().any(|cr| match cr {
                rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
                    face,
                    ..
                } => *face == fkey,
                _ => false,
            })
        }
        _ => false,
    }
}

/// OCCT static FindInternalIntersections(theEdge, theFace, Splits,
/// isOverlapped) (cxx L1291-1511).
pub(crate) fn find_internal_intersections(
    the_edge: &Shape,
    the_face: &Shape,
    splits: &mut indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)>,
    is_overlapped: &mut bool,
) {
    // OCCT L1298.
    let tol_ext = rcad_kernel::precision::PCONFUSION;

    // OCCT L1301: BRepAdaptor_Surface anAdSurf(theFace, false) — unused in
    // the OCCT body; kept as a comment for form.
    let mut split_pars: Vec<f64> = Vec::new(); // OCCT L1302

    // OCCT L1304-1317.
    let (v_first, v_last) = top_exp_vertices(the_edge);
    let the_pnt = [
        v_first
            .as_ref()
            .and_then(brep_tool_pnt)
            .unwrap_or(DVec3::ZERO),
        v_last
            .as_ref()
            .and_then(brep_tool_pnt)
            .unwrap_or(DVec3::ZERO),
    ];
    let a_tol_v = [
        v_first.as_ref().map(brep_tool_tolerance).unwrap_or(0.0),
        v_last.as_ref().map(brep_tool_tolerance).unwrap_or(0.0),
    ];
    // OCCT L1313: ext = 16. — = 4 * 4 - to avoid creating microedges, area
    // around vertices is increased up to 4 vertex tolerance. Such approach
    // is usual for other topological algorithms, for example, Boolean
    // Operations.
    let ext = 16.0;
    let a_tol_v_ext = [
        ext * a_tol_v[0] * a_tol_v[0],
        ext * a_tol_v[1] * a_tol_v[1],
    ];

    // OCCT L1319-1321: the edge pcurve box.
    let Some((the_pcurve, pf0, pl0)) = brep_tool_curve_on_surface(the_edge, the_face) else {
        return;
    };
    let mut the_box = BndBox2d::new();
    bnd_lib_add2d_curve(
        &the_pcurve,
        pf0,
        pl0,
        brep_tool_tolerance(the_edge),
        &mut the_box,
    );

    // OCCT L1323-1329: the curve adaptor and vertex parameter tolerances.
    let Some((the_curve, the_par0, the_par1)) = brep_tool_curve(the_edge) else {
        return;
    };
    let the_par = [the_par0, the_par1];
    let a_tol_v2d = [
        curve_resolution_for(&the_curve, a_tol_v[0]).max(rcad_kernel::precision::PCONFUSION),
        curve_resolution_for(&the_curve, a_tol_v[1]).max(rcad_kernel::precision::PCONFUSION),
    ];
    // OCCT L1330.
    let mut a_dist_max =
        rcad_kernel::precision::CONFUSION * rcad_kernel::precision::CONFUSION;

    // OCCT L1331-1419: the face-edge loop.
    for an_edge in explorer(the_face, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L1335-1338: the pcurve box of the face edge.
        let Some((a_pcurve, af, al)) = brep_tool_curve_on_surface(&an_edge, the_face) else {
            continue;
        };
        let mut a_box = BndBox2d::new();
        bnd_lib_add2d_curve(&a_pcurve, af, al, brep_tool_tolerance(&an_edge), &mut a_box);
        if the_box.is_out_box(&a_box) {
            continue;
        }

        // OCCT L1343-1350: the extrema (architecture difference #11 —
        // pending; IsDone() is false so the OCCT continue is taken).
        let Some((a_curve, a_fpar, a_lpar)) = brep_tool_curve(&an_edge) else {
            continue;
        };
        let an_extrema = ExtremaExtCC::new(
            &the_curve,
            the_par[0],
            the_par[1],
            &a_curve,
            a_fpar,
            a_lpar,
        );

        if !an_extrema.is_done() || an_extrema.nb_ext() == 0 {
            continue;
        }

        // OCCT L1352-1418: unreachable while the pending ExtremaExtCC keeps
        // is_done()==false; translated for form against the pending
        // accessor surface.
        let a_nb_ext = an_extrema.nb_ext();
        let max_tol = brep_tool_tolerance(&an_edge);
        let a_max_tol2 = max_tol * max_tol;
        if an_extrema.is_parallel() && an_extrema.square_distance(1) <= a_max_tol2 {
            *is_overlapped = true;
            return;
        }
        // Check extremity distances (OCCT L1360-1379).
        let mut dists = [0.0f64; 4];
        let (mut a_p11, mut a_p12, mut a_p21, mut a_p22) =
            (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        an_extrema.trimmed_square_distances(
            &mut dists,
            &mut a_p11,
            &mut a_p12,
            &mut a_p21,
            &mut a_p22,
        );
        let _ = (&a_p11, &a_p12, &a_p21, &a_p22);
        for i in 0..4 {
            let j = if i < 2 { 0usize } else { 1usize };
            if dists[i] < a_tol_v_ext[j] / ext {
                return;
            }
        }

        // OCCT L1381-1418: the extrema loop.
        for i in 1..=a_nb_ext {
            let a_dist = an_extrema.square_distance(i);
            if a_dist > a_max_tol2 {
                continue;
            }

            // OCCT L1389-1392.
            let mut a_p_on_c1 = (0.0, DVec3::ZERO);
            let mut a_p_on_c2 = (0.0, DVec3::ZERO);
            an_extrema.points(i, &mut a_p_on_c1, &mut a_p_on_c2);
            let the_int_par = a_p_on_c1.0;
            let an_int_par = a_p_on_c2.0;
            // OCCT L1393-1399.
            let mut j = 2usize;
            for jj in 0..2 {
                if (the_int_par - the_par[jj]).abs() <= a_tol_v2d[jj] {
                    j = jj;
                    break;
                }
            }
            // intersection found in the middle of the edge (OCCT L1400).
            if j >= 2 {
                // intersection is inside "theEdge" => split
                let a_point = a_curve.point_at(an_int_par);
                let a_point_int = the_curve.point_at(the_int_par);

                if a_point_int.distance_squared(the_pnt[0]) > a_tol_v_ext[0]
                    && a_point_int.distance_squared(the_pnt[1]) > a_tol_v_ext[1]
                    && a_point.distance_squared(the_pnt[0]) > a_tol_v_ext[0]
                    && a_point.distance_squared(the_pnt[1]) > a_tol_v_ext[1]
                {
                    split_pars.push(the_int_par);
                    if a_dist > a_dist_max {
                        a_dist_max = a_dist;
                    }
                }
            }
        }
    }
    let _ = tol_ext;

    // OCCT L1421-1424.
    if split_pars.is_empty() {
        return;
    }

    // OCCT L1426-1438: Sort (the OCCT 1-based bubble sort).
    for i in 1..split_pars.len() {
        for j in (i + 1)..=split_pars.len() {
            if split_pars[i - 1] > split_pars[j - 1] {
                let tmp = split_pars[i - 1];
                split_pars[i - 1] = split_pars[j - 1];
                split_pars[j - 1] = tmp;
            }
        }
    }

    // OCCT L1439-1453: Remove repeating points.
    let mut i = 1usize;
    while i < split_pars.len() {
        let pnt1 = the_curve.point_at(split_pars[i - 1]);
        let pnt2 = the_curve.point_at(split_pars[i]);
        if pnt1.distance_squared(pnt2)
            <= rcad_kernel::precision::CONFUSION * rcad_kernel::precision::CONFUSION
        {
            split_pars.remove(i);
        } else {
            i += 1;
        }
    }

    // OCCT L1455-1505: Split.
    let mut new_edges: Vec<Shape> = Vec::new();
    let mut pool = rcad_kernel::topo::topods::BRep::new();

    let mut first_vertex = v_first.clone().unwrap_or_else(Shape::null);
    let mut last_vertex;
    let mut first_par = the_par[0];
    let mut last_par;
    for i in 1..=(split_pars.len() + 1) {
        first_vertex.orientation = Orientation::Forward; // OCCT L1463
        if i <= split_pars.len() {
            // OCCT L1465-1471.
            last_par = split_pars[i - 1];
            let last_point = the_curve.point_at(last_par);
            // OCCT L1468: BRepLib_MakeVertex(LastPoint).
            last_vertex = pool.add_tvertex(last_point);
            // OCCT L1469-1470: aB.UpdateVertex(LastVertex, sqrt(aDistMax)).
            builder_update_vertex(&mut last_vertex, a_dist_max.sqrt());
        } else {
            // OCCT L1473-1475.
            last_par = the_par[1];
            last_vertex = v_last.clone().unwrap_or_else(Shape::null);
        }
        last_vertex.orientation = Orientation::Reversed; // OCCT L1477

        // OCCT L1479-1486.
        let an_orient = the_edge.orientation;
        let mut new_edge = pool.empty_copied(the_edge);
        new_edge.orientation = Orientation::Forward;
        builder_range(&mut new_edge, first_par, last_par);
        builder_add_edge_vertex(&mut new_edge, &first_vertex);
        builder_add_edge_vertex(&mut new_edge, &last_vertex);
        new_edge.orientation = an_orient;

        // OCCT L1487-1493: ShapeAnalysis_Edge::CheckSameParameter
        // (architecture difference #12).
        let mut amaxdev = 0.0;
        if shape_analysis_edge_check_same_parameter(&new_edge, the_face, &mut amaxdev) {
            builder_update_edge_tolerance(&mut new_edge, amaxdev);
        }

        // OCCT L1495-1502.
        if an_orient == Orientation::Forward {
            new_edges.push(new_edge);
        } else {
            new_edges.insert(0, new_edge);
        }
        // OCCT L1503-1504.
        first_vertex = last_vertex;
        first_par = last_par;
    }

    // OCCT L1507-1510.
    if !new_edges.is_empty() {
        splits.insert(shape_key(the_edge), (the_edge.clone(), new_edges));
    }
}
