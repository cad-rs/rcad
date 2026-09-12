// OCCT BRepOffset_Inter2d.cxx L1-1187 + BRepOffset_Inter2d.hxx L37-119 —
// 1:1 translation (part A: the file statics through ExtendPCurve; the class
// entry points Compute / ConnexIntByInt / ConnexIntByIntInVert / FuseVertices
// / ExtentEdge / UpdateVertex / MakeChain live in brep_offset_inter2d_b.rs).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_Inter2d.cxx / .hxx
//
// OCCT inheritance chain (hxx L37): none — BRepOffset_Inter2d is a
// static-only class.
//
// Architecture differences (continuing the numbering of the offset package):
// 21. NCollection_DataMap / IndexedDataMap / IndexedMap / Map / List /
//     Sequence over TopoDS_Shape -> HashMap / IndexMap / HashSet / Vec keyed
//     by (TShape ptr, Location) (TopTools_ShapeMapHasher, orientation
//     ignored); the IndexedDataMap keeps the key shape (FindKey).
// 22. The OCCT global TShape arena: `const TopoDS_Edge&` parameters are
//     mutated in place through BRep_Builder (UpdateVertex); the rcad Shape is
//     a value handle — the callee mutates its own clone and the working copy
//     is annotated at each site (pool architecture pending, port plan §0.6).
// 23. BRepOffset_Tool::MapVertexEdges (BRepOffset_Tool.cxx) — consumed by
//     ConnexIntByInt; GAP leaf below (the parallel BRepOffset_Tool module
//     carries the 1:1 translation; switch to
//     `super::brep_offset_tool::BRepOffsetTool::map_vertex_edges` when it
//     lands).
// 24. BRepOffset_Analyse::EdgeReplacement (BRepOffset_Analyse.cxx) — the
//     Analyse is not translated yet; GAP carrier below.
// 25. BRepLib::BuildCurve3d — GAP leaf (the brep_offset_offset.rs carrier).
// 26. BRepLib::SameParameter — GAP no-op leaf (the brep_offset_offset_b.rs
//     precedent).
// 27. Geom2dConvert_CompCurveToBSplineCurve / GeomConvert_CompCurveToBSplineCurve
//     (the incremental Add + BSplineCurve forms) — GAP carriers below.
// 28. GeomProjLib::Curve2d(C, f, l, S) — GAP leaf (the
//     loc_ope_wires_on_shape_b.rs precedent).
// 29. GeomLib::BuildCurve3d — GAP leaf.
// 30. GeomAPI_ProjectPointOnCurve — the real class (1:1 over the kernel
//     Extrema_ExtPC translation) in its OCCT toolkit home,
//     `crate::geomalgo::geom_api_project_point_on_curve`.
// 31. Geom2dInt_GInter — the rcad TheIntPCurvePCurveOfGInter vehicle (the
//     brep_blend_surf_rst_line_builder_b.rs precedent); the OCCT
//     (GAC1, GAC2, TolConf, Tol) constructor form maps to the explicit-range
//     constructor below, the Perform(GAcurve, GAline, TolConf, Tol) form to
//     perform_natural (natural domains).
// 32. BRep_TEdge::ChangeCurves() / BRep_CurveRepresentation::Surface() —
//     the rcad CurveRepresentation carries the face KEY, not the surface
//     value; the Surface() read is a GAP leaf.  The Locations are identity in
//     this pipeline (the brep_algo/loop.rs header convention).
// 33. Message_ProgressRange / Message_ProgressScope — the rcad translation
//     carries the parameter as () (no progress surface yet).
// 34. Geom2d/Geom Line/TrimmedCurve constructions map to the rcad Curve2d /
//     Curve3 enum variants; gp_Dir2d normalization is the Line2d::new
//     invariant.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use glam::{DVec2, DVec3};
use indexmap::IndexMap;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3, TrimmedCurve2,
};
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_curve_on_surface, brep_tool_parameter, brep_tool_pnt,
    brep_tool_range, brep_tool_surface, brep_tool_tolerance, builder_make_vertex,
    builder_update_vertex_point_tol, builder_update_vertex_tol, explorer, oriented, shape_key,
    sub_shapes, top_exp_vertices_raw, ShapeKey,
};
use crate::feat::brep_feat_rib_slot::top_exp_last_vertex;
use crate::feat::loc_ope_generator_b::brep_tools_is_really_closed;
use crate::feat::loc_ope_wires_on_shape_b::{brep_tool_degenerated, BRepAdaptorCurve2d};
use crate::fillet::chfi3d_builder_0::topexp_common_vertex;
use crate::geomalgo::geom2d_int::TheIntPCurvePCurveOfGInter;
use crate::geomalgo::geom_api_project_point_on_curve::GeomAPIProjectPointOnCurve;
use crate::geomalgo::int_res2d::{Domain as Res2dDomain, IntersectionBase};

use super::brep_offset_offset::brep_lib_build_curve3d;

// ---------------------------------------------------------------------------
// Shared map types (architecture difference #21).
// ---------------------------------------------------------------------------

/// OCCT NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher> — the key shape is kept for FindKey.
pub type DmvvMap = IndexMap<ShapeKey, (Shape, Vec<Shape>)>;

/// OCCT NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>.
pub type IndexedShapeMap = IndexMap<ShapeKey, Shape>;

/// OCCT NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>
/// over vertices (the GetEdgesOrientedInFace aVEmap).
type VEMap = IndexMap<ShapeKey, (Shape, Vec<Shape>)>;

// ---------------------------------------------------------------------------
// GAP / architecture carriers (architecture differences #22-#33).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Analyse (BRepOffset_Analyse.hxx) — the real body lives in
/// super::brep_offset_analyse (the E0 carrier-switch list; the local
/// EdgeReplacement panic carrier is deleted).  Re-exported for the Inter2d
/// consumers (brep_offset_inter2d_b globs this module).
pub use super::brep_offset_analyse::BRepOffsetAnalyse;

/// OCCT BRepLib::SameParameter(E, Tol, DoAlsoMinmax) — GAP no-op leaf
/// (architecture difference #26; the brep_offset_offset_b.rs precedent).
pub fn brep_lib_same_parameter(_the_e: &Shape, _the_tol: f64) {}

/// OCCT GeomProjLib::Curve2d(C, f, l, S) — GAP leaf (architecture difference
/// #28).
pub fn geom_proj_lib_curve2d(_the_c: &Curve3, _the_f: f64, _the_l: f64, _the_s: &Surface3) -> Option<Curve2d> {
    panic!("GAP: GeomProjLib::Curve2d (TKTopAlgo/GeomProjLib not translated)");
}

/// OCCT GeomLib::BuildCurve3d(Tol, ConS, f, l, C3d, MaxDeviation,
/// AverageDeviation, Continuity, MaxDegree, MaxSegment) — GAP leaf
/// (architecture difference #29).
pub fn geom_lib_build_curve3d(
    _the_c3d: &mut Option<Curve3>,
    _the_max_deviation: &mut f64,
    _the_average_deviation: &mut f64,
) {
    panic!("GAP: GeomLib::BuildCurve3d (TKTopAlgo/GeomLib not translated)");
}

/// OCCT Geom2dConvert_CompCurveToBSplineCurve (TKGeomBase/Geom2dConvert) —
/// GAP carrier (architecture difference #27): the incremental
/// Add(Segment, Tol) + BSplineCurve() engine.
pub struct Geom2dConvertCompCurveToBSplineCurve;

impl Geom2dConvertCompCurveToBSplineCurve {
    /// OCCT Geom2dConvert_CompCurveToBSplineCurve(BasisCurve, Convert).
    pub fn new(_the_basis: &Curve2d) -> Self {
        panic!("GAP: Geom2dConvert_CompCurveToBSplineCurve (TKGeomBase/Geom2dConvert not translated)");
    }

    /// OCCT Add(NewCurve, Tol) — appends the segment with C1 continuity.
    pub fn add(&mut self, _the_new_curve: &Curve2d, _the_tol: f64) -> bool {
        panic!("GAP: Geom2dConvert_CompCurveToBSplineCurve::Add");
    }

    /// OCCT BSplineCurve().
    pub fn bspline_curve(&self) -> Curve2d {
        panic!("GAP: Geom2dConvert_CompCurveToBSplineCurve::BSplineCurve");
    }
}

/// OCCT GeomConvert_CompCurveToBSplineCurve (TKGeomBase/GeomConvert) — GAP
/// carrier (architecture difference #27).
pub struct GeomConvertCompCurveToBSplineCurve;

impl GeomConvertCompCurveToBSplineCurve {
    /// OCCT GeomConvert_CompCurveToBSplineCurve(BasisCurve, Convert).
    pub fn new(_the_basis: &Curve3) -> Self {
        panic!("GAP: GeomConvert_CompCurveToBSplineCurve (TKGeomBase/GeomConvert not translated)");
    }

    /// OCCT Add(NewCurve, Tol).
    pub fn add(&mut self, _the_new_curve: &Curve3, _the_tol: f64) -> bool {
        panic!("GAP: GeomConvert_CompCurveToBSplineCurve::Add");
    }

    /// OCCT BSplineCurve().
    pub fn bspline_curve(&self) -> Curve3 {
        panic!("GAP: GeomConvert_CompCurveToBSplineCurve::BSplineCurve");
    }
}

/// OCCT BRepAdaptor_Surface(F) re-host — the face surface value (architecture
/// difference: the rcad surface is the TFaceData payload; BAsurf.Value(u, v)
/// evaluates it, the face location is identity in this pipeline).
pub struct BRepAdaptorSurface {
    my_surface: Surface3,
}

impl BRepAdaptorSurface {
    /// OCCT BRepAdaptor_Surface(F).
    pub fn new(the_face: &Shape) -> Self {
        let surf = brep_tool_surface(the_face)
            .unwrap_or_else(|| panic!("BRepAdaptor_Surface: null face surface"));
        BRepAdaptorSurface { my_surface: surf }
    }

    /// OCCT BRepAdaptor_Surface::Value(u, v).
    pub fn value(&self, the_u: f64, the_v: f64) -> DVec3 {
        use rcad_kernel::geom::SurfaceEval;
        self.my_surface.point_at(the_u, the_v)
    }
}

/// OCCT BRepAdaptor_Curve(E, F) re-host — the edge 3D curve with its range
/// (architecture difference: the OCCT (E, F) form falls back to the curve on
/// surface F when the edge has no 3D curve; the rcad carrier GAP-panics on
/// that path — the pcurve-to-3D evaluation is not translated).
pub struct BRepAdaptorCurve {
    my_curve: Curve3,
    my_first: f64,
    my_last: f64,
}

impl BRepAdaptorCurve {
    /// OCCT BRepAdaptor_Curve(E, F).
    pub fn new(the_e: &Shape, _the_f: &Shape) -> Self {
        match brep_tool_curve(the_e) {
            Some((c, f, l)) => BRepAdaptorCurve {
                my_curve: c,
                my_first: f,
                my_last: l,
            },
            None => panic!("GAP: BRepAdaptor_Curve(E, F) — the curve-on-surface 3D fallback"),
        }
    }

    /// OCCT FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.my_last
    }

    /// OCCT Value(u).
    pub fn value(&self, the_u: f64) -> DVec3 {
        self.my_curve.point_at(the_u)
    }
}

// OCCT GeomAPI_ProjectPointOnCurve (architecture difference #30) lives in the
// geomalgo home of its OCCT toolkit (TKGeomAlgo/GeomAPI):
// `crate::geomalgo::geom_api_project_point_on_curve`.  The former local
// re-host is deleted; the consumers below use the real class.

/// OCCT Geom2dInt_GInter re-host (architecture difference #31) — the
/// TheIntPCurvePCurveOfGInter vehicle; the results surface (IsDone /
/// NbPoints / Point / NbSegments / Segment) mirrors the OCCT GInter.
pub struct Geom2dIntGInter {
    pub base: crate::geomalgo::int_res2d::IntersectionBase,
}

impl Geom2dIntGInter {
    fn engine(
        c1: &Curve2d,
        d1: &Res2dDomain,
        c2: &Curve2d,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> IntersectionBase {
        let mut inter = TheIntPCurvePCurveOfGInter::new();
        inter.perform(c1, d1, c2, d2, tol_conf, tol);
        inter.base
    }

    /// OCCT Geom2dInt_GInter(GAC1, GAC2, TolConf, Tol) — the constructor form
    /// over the bounded (pcurve, first, last) adaptor pair.
    pub fn new(
        gac1: &Curve2d,
        f1: f64,
        l1: f64,
        gac2: &Curve2d,
        f2: f64,
        l2: f64,
        the_tol_conf: f64,
        the_tol: f64,
    ) -> Self {
        let d1 = Res2dDomain::bounded(
            gac1.point_at(f1),
            f1,
            the_tol_conf,
            gac1.point_at(l1),
            l1,
            the_tol_conf,
        );
        let d2 = Res2dDomain::bounded(
            gac2.point_at(f2),
            f2,
            the_tol_conf,
            gac2.point_at(l2),
            l2,
            the_tol_conf,
        );
        Geom2dIntGInter {
            base: Self::engine(gac1, &d1, gac2, &d2, the_tol_conf, the_tol),
        }
    }

    /// OCCT IntCC.Perform(GAcurve, GAline, TolConf, Tol) — the Perform form
    /// over the adaptors' own natural domains (an infinite natural domain
    /// maps to the infinite IntRes2d_Domain).
    pub fn perform_natural(
        &mut self,
        ga_curve: &Curve2d,
        ga_line: &Curve2d,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        let natural = |c: &Curve2d, tol: f64| -> Res2dDomain {
            let dom = c.default_domain();
            // Bounded-domain test via Precision::IsInfinite
            // (Precision.hxx L350-353): unbounded curves carry 2e100.
            if !rcad_kernel::precision::is_infinite_value(dom[0])
                && !rcad_kernel::precision::is_infinite_value(dom[1])
            {
                Res2dDomain::bounded(
                    c.point_at(dom[0]),
                    dom[0],
                    tol,
                    c.point_at(dom[1]),
                    dom[1],
                    tol,
                )
            } else {
                Res2dDomain::infinite()
            }
        };
        let d1 = natural(ga_curve, the_tol_conf);
        let d2 = natural(ga_line, the_tol_conf);
        self.base = Self::engine(ga_curve, &d1, ga_line, &d2, the_tol_conf, the_tol);
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// OCCT Point(i) — 1-based.
    pub fn point(&self, the_i: usize) -> &crate::geomalgo::int_res2d::IntersectionPoint {
        self.base.point(the_i)
    }

    /// OCCT NbSegments().
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// OCCT Segment(i) — 1-based.
    pub fn segment(&self, the_i: usize) -> &crate::geomalgo::int_res2d::IntersectionSegment {
        self.base.segment(the_i)
    }
}

/// OCCT BRepTools_WireExplorer re-host (architecture difference #7 of the
/// MakeSimpleOffset module / the brep_feat_make_linear_form.rs precedent):
/// the wire edge list (the caller passes the FORWARD-oriented wire; the face
/// argument is the OCCT form and unused in the rcad reduction).
pub struct BRepToolsWireExplorer {
    my_edges: Vec<Shape>,
    my_index: usize,
}

impl BRepToolsWireExplorer {
    /// OCCT BRepTools_WireExplorer wexp; (the default constructor).
    pub fn new() -> Self {
        BRepToolsWireExplorer {
            my_edges: Vec::new(),
            my_index: 0,
        }
    }

    /// OCCT Init(W, F).
    pub fn init(&mut self, the_w: &Shape, _the_f: &Shape) {
        self.my_edges = sub_shapes(the_w);
        self.my_index = 0;
    }

    /// OCCT More().
    pub fn more(&self) -> bool {
        self.my_index < self.my_edges.len()
    }

    /// OCCT Current() — the cached myEdge (BRepTools_WireExplorer.hxx L103);
    /// on the exhausted explorer OCCT Current() returns the null edge
    /// (Next() nulls myEdge, BRepTools_WireExplorer.cxx L393-413).
    pub fn current(&self) -> Shape {
        match self.my_edges.get(self.my_index) {
            Some(e) => e.clone(),
            None => Shape::null(),
        }
    }

    /// OCCT Next().
    pub fn next(&mut self) {
        self.my_index += 1;
    }
}

impl Default for BRepToolsWireExplorer {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT Adaptor3d_CurveOnSurface re-host — the (pcurve, surface) pair carried
/// for evaluateMaxSegment / GeomLib::BuildCurve3d.
pub struct Adaptor3dCurveOnSurface {
    my_surface: Surface3,
    my_curve: Curve2d,
}

impl Adaptor3dCurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface(HC2d, HSurf).
    pub fn new(the_curve: &Curve2d, the_surface: &Surface3) -> Self {
        Adaptor3dCurveOnSurface {
            my_surface: the_surface.clone(),
            my_curve: the_curve.clone(),
        }
    }

    /// OCCT GetSurface().
    pub fn get_surface(&self) -> &Surface3 {
        &self.my_surface
    }

    /// OCCT GetCurve().
    pub fn get_curve(&self) -> &Curve2d {
        &self.my_curve
    }
}

/// OCCT BRep_CurveRepresentation::Surface() (architecture difference #32):
/// the rcad representation carries the face key (TShape pointer + location)
/// instead of the surface value; the surface resolves through the shared
/// BRep pool.
pub fn curve_rep_surface(
    the_brep: &rcad_kernel::topods::BRep,
    the_rep: &rcad_kernel::topods::CurveRepresentation,
) -> Surface3 {
    let face_key = match the_rep {
        rcad_kernel::topods::CurveRepresentation::CurveOnSurface { face, .. } => *face,
        rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
        _ => panic!("BRep_CurveRepresentation::Surface: not a curve-on-surface"),
    };
    let idx = the_brep
        .index_by_ptr(face_key.0)
        .unwrap_or_else(|| panic!("curve_rep_surface: face key not in the shared pool"));
    the_brep
        .shape_at(idx)
        .as_face()
        .and_then(|fd| fd.surface.clone())
        .expect("curve_rep_surface: null face surface")
}

/// OCCT BRep_Builder::UpdateVertex(V, P, E, Tol) — the vertex tolerance, the
/// stored parameter on the edge, and the edge tolerance
/// (BRep_Builder.cxx UpdateVertex Curve branch).
pub fn builder_update_vertex_on_edge(the_v: &mut Shape, the_u: f64, the_e: &mut Shape, the_tol: f64) {
    if let TShape::Vertex(vd) = Arc::make_mut(&mut the_v.data) {
        vd.tolerance = vd.tolerance.max(the_tol);
    }
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.vertex_params.insert(the_v.ptr_id(), the_u);
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C, Tol) — bind the 3D curve and raise the
/// tolerance (BRep_Builder.cxx UpdateEdge 3D curve branch).
pub fn builder_update_edge_curve(the_e: &mut Shape, the_c: Option<Curve3>, the_tol: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.curve = the_c;
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRepLib_MakeVertex(P) (BRepLib_MakeVertex.cxx L26-31):
/// B.MakeVertex(myShape, P, BRepLib::Precision()) — the BRepLib precision
/// default is Precision::Confusion.
pub fn brep_lib_make_vertex(the_p: DVec3) -> Shape {
    let mut v = builder_make_vertex();
    builder_update_vertex_point_tol(&mut v, the_p, rcad_kernel::precision::CONFUSION);
    v
}

/// OCCT BOPTools_AlgoTools::MakeVertex(LV, VNew) (BOPTools_AlgoTools.cxx
/// L1798-1813) — the single-vertex copy / BRepLib::BoundingVertex + MakeVertex
/// forms; the rcad BRepLib::BoundingVertex stand-in is the bop
/// centroid+max-tolerance helper (crate::bop::tools::algo_tools::make_vertex).
pub fn bop_tools_algo_tools_make_vertex(a_lv: &[Shape]) -> Shape {
    // OCCT L1800-1804: aNb == 1 -> aVnew = First().
    if a_lv.len() == 1 {
        return a_lv[0].clone();
    }
    // OCCT L1805-1812: BRepLib::BoundingVertex(aLV, aNC, aNTol);
    //                  aBB.MakeVertex(aVnew, aNC, aNTol).
    let pairs: Vec<(DVec3, f64)> = a_lv
        .iter()
        .map(|s| (brep_tool_pnt(s).unwrap_or(DVec3::ZERO), brep_tool_tolerance(s)))
        .collect();
    let (a_nc, a_ntol) = crate::bop::tools::algo_tools::make_vertex(&pairs);
    let mut a_v_new = builder_make_vertex();
    builder_update_vertex_point_tol(&mut a_v_new, a_nc, a_ntol);
    a_v_new
}

/// OCCT BRep_Tool::Parameter(V, E, T) — the Standard_Boolean overload
/// (BRep_Tool.cxx): true when the vertex carries a stored parameter on the
/// edge.
pub fn brep_tool_parameter_bool(the_v: &Shape, the_e: &Shape) -> Option<f64> {
    let ed = the_e.as_edge()?;
    ed.vertex_params.get(&the_v.ptr_id()).copied()
}

// ---------------------------------------------------------------------------
// gp / precision leaf vehicles.
// ---------------------------------------------------------------------------

/// OCCT gp_Vec2d operator^ — the 2D cross product (the Z component).
pub(super) fn cross2d(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// OCCT gp_Dir2d::Angle(Other) (gp_Dir2d.cxx L26-66) — the angle in [-PI, PI].
fn gp_dir2d_angle(a: DVec2, b: DVec2) -> f64 {
    let cosinus = a.dot(b);
    let sinus = cross2d(a, b);
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

/// OCCT gp_Dir2d::IsParallel(Other, AngularTolerance) — IsEqual or
/// IsOpposite: |cross| <= tolerance.
pub(super) fn gp_dir2d_is_parallel(a: DVec2, b: DVec2, the_angular_tolerance: f64) -> bool {
    cross2d(a, b).abs() <= the_angular_tolerance
}

/// OCCT IntTools_Tools::ComputeIntRange(Tol1, Tol2, Angle)
/// (IntTools_Tools.cxx L783-809).
fn int_tools_compute_int_range(the_tol1: f64, the_tol2: f64, the_angle: f64) -> f64 {
    if (std::f64::consts::FRAC_PI_2 - the_angle).abs() < rcad_kernel::precision::ANGULAR {
        the_tol2
    } else {
        let an_angle = if the_angle > std::f64::consts::FRAC_PI_2 {
            std::f64::consts::PI - the_angle
        } else {
            the_angle
        };
        let a1 = the_tol1 * (std::f64::consts::FRAC_PI_2 - an_angle).tan();
        let a2 = the_tol2 / an_angle.sin();
        a1 + a2
    }
}

/// OCCT Geom2d_Curve::Period() — the natural domain length (the rcad
/// vehicle; OCCT conics carry 2*PI).
pub(super) fn curve2d_period(the_curve: &Curve2d) -> f64 {
    let dom = the_curve.default_domain();
    dom[1] - dom[0]
}

/// OCCT Geom_Curve::Period() — the natural domain length (the rcad vehicle).
pub(super) fn curve3_period(the_curve: &Curve3) -> f64 {
    let dom = the_curve.default_domain();
    dom[1] - dom[0]
}

/// OCCT Geom2d_BSplineCurve::NbKnots() — the DISTINCT knot count (the rcad
/// knot vector carries the expanded multiplicities).
pub(super) fn distinct_knots_1d(knots: &[f64]) -> usize {
    let mut n = 0usize;
    let mut prev = f64::NAN;
    for k in knots {
        if *k != prev {
            n += 1;
            prev = *k;
        }
    }
    n
}

// ---------------------------------------------------------------------------
// OCCT TopExp::MapShapesAndAncestors(S, TS, TA, M) re-host
// (TopExp.cxx L87-119; the brep_algo/loop.rs model).
// ---------------------------------------------------------------------------

fn top_exp_map_shapes_and_ancestors(
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut VEMap,
) {
    for a_anc in explorer(s, ta, ShapeType::Shape) {
        for a_exs in explorer(&a_anc, ts, ShapeType::Shape) {
            let key = shape_key(&a_exs);
            let entry = m.entry(key).or_insert((a_exs.clone(), Vec::new()));
            entry.1.push(a_anc.clone());
        }
    }
    for a_ex in explorer(s, ts, ta) {
        let key = shape_key(&a_ex);
        m.entry(key).or_insert((a_ex.clone(), Vec::new()));
    }
}

// =========================================================================
// OCCT BRepOffset_Inter2d.cxx — the file statics.
// =========================================================================

/// OCCT CommonVertex(E1, E2) (BRepOffset_Inter2d.cxx L81-99) — the file-local
/// static: the common vertex checked with the cumulated orientation pair;
/// returns a null shape when there is none (the OCCT default TopoDS_Vertex).
pub(super) fn common_vertex(e1: &Shape, e2: &Shape) -> Shape {
    // OCCT L83-86: TopExp::Vertices(E, V1[2], V2[2], true) — CumOri=true.
    let (v1f, v1l) = top_exp_vertices(e1, true);
    let (v2f, v2l) = top_exp_vertices(e2, true);
    let v10 = v1f.unwrap_or_else(Shape::null);
    let v11 = v1l.unwrap_or_else(Shape::null);
    let v20 = v2f.unwrap_or_else(Shape::null);
    let v21 = v2l.unwrap_or_else(Shape::null);

    // OCCT L87-96: the last vertex of the first edge is checked first.
    if v11.is_same(&v20) || v11.is_same(&v21) {
        return v11;
    }
    if v10.is_same(&v20) || v10.is_same(&v21) {
        return v10;
    }
    // OCCT L98: return V; (the null vertex).
    Shape::null()
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri) (TopExp.cxx L214-253) —
/// with CumOri the iterator (TopoDS_Iterator(E, CumOri=true)) composes the
/// occurrence orientation into the children FIRST and the composed-FORWARD
/// child is Vfirst, the composed-REVERSED one is Vlast (the
/// brep_feat_rib_slot.rs top_exp_first_vertex convention).
fn top_exp_vertices(e: &Shape, cum_ori: bool) -> (Option<Shape>, Option<Shape>) {
    let (vf, vl) = top_exp_vertices_raw(e);
    if !cum_ori {
        return (vf, vl);
    }
    let compose = |mut v: Shape| {
        v.orientation = e.orientation.compose(v.orientation);
        v
    };
    if e.orientation == Orientation::Reversed {
        (vl.map(compose), vf.map(compose))
    } else {
        (vf.map(compose), vl.map(compose))
    }
}

/// OCCT DefineClosedness(theFace) (BRepOffset_Inter2d.cxx L101-126) —
/// 1 when a really-closed edge has a pcurve parallel to OY, 2 when parallel
/// to OX, 0 when none.
fn define_closedness(the_face: &Shape) -> i32 {
    // OCCT L103-106: TopExp_Explorer(theFace, TopAbs_EDGE) walk.
    for an_edge in explorer(the_face, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L107: BRepTools::IsReallyClosed(anEdge, theFace).
        if brep_tools_is_really_closed(&an_edge, the_face) {
            // OCCT L109-111: the pcurve + DN(fpar, 1) tangent.
            let Some((a_pcurve, fpar, _lpar)) = brep_tool_curve_on_surface(&an_edge, the_face)
            else {
                continue;
            };
            let a_tangent = a_pcurve.derivative_at(fpar);
            // OCCT L112-113: aTangent ^ gp::DX2d() / ^ gp::DY2d().
            let a_cross_prod1 = cross2d(a_tangent, DVec2::X);
            let a_cross_prod2 = cross2d(a_tangent, DVec2::Y);
            // OCCT L114-121: |cross2| < |cross1| — parallel to OY.
            if a_cross_prod2.abs() < a_cross_prod1.abs() {
                return 1;
            } else {
                return 2;
            }
        }
    }

    // OCCT L125.
    0
}

/// OCCT GetEdgesOrientedInFace(theShape, theFace, theAsDes, theSeqEdges)
/// (BRepOffset_Inter2d.cxx L128-266) — orders the descendants of the face
/// that occur in theShape into a connected wire-order sequence.
pub(super) fn get_edges_oriented_in_face(
    the_shape: &Shape,
    the_face: &Shape,
    the_as_des: &BRepAlgoAsDes,
    the_seq_edges: &mut Vec<Shape>,
) {
    // OCCT L133: the descendants of the face (the rcad borrow-split clone;
    // the list is not mutated inside this function).
    let a_edges = the_as_des.descendant(the_face).to_vec();

    // OCCT L135-150: keep the shape's edges that occur in aEdges, in the
    // explorer order, with the in-face representation.
    for an_edge in explorer(the_shape, ShapeType::Edge, ShapeType::Shape) {
        for an_edge_in_face in &a_edges {
            if an_edge_in_face.is_same(&an_edge) {
                the_seq_edges.push(an_edge_in_face.clone());
                break;
            }
        }
    }

    // OCCT L151-154.
    if the_seq_edges.len() == 1 {
        return;
    }

    // OCCT L156-161: MapShapesAndAncestors(seq(ii), VERTEX, EDGE, aVEmap).
    let mut a_ve_map: VEMap = IndexMap::new();
    for ii in 1..=the_seq_edges.len() {
        top_exp_map_shapes_and_ancestors(
            &the_seq_edges[ii - 1],
            ShapeType::Vertex,
            ShapeType::Edge,
            &mut a_ve_map,
        );
    }

    // OCCT L163-181: the first free-end vertex (with orientation).
    let mut a_first_vertex = Shape::null();
    let mut a_first_edge = Shape::null();
    for ii in 1..=a_ve_map.len() {
        let (a_vertex, a_elist) = &a_ve_map[ii - 1];
        if a_elist.len() == 1 {
            let an_edge = &a_elist[0];
            // OCCT L173: TopExp::Vertices(anEdge, aV1, aV2, true).
            let (a_v1, _a_v2) = top_exp_vertices(an_edge, true);
            if let Some(a_v1) = a_v1 {
                if a_v1.is_same(a_vertex) {
                    a_first_vertex = a_vertex.clone();
                    a_first_edge = an_edge.clone();
                    break;
                }
            }
        }
    }

    // OCCT L183-231: the closed set of edges — the largest coordinate gap
    // vertex starts the walk.
    if a_first_edge.is_null() {
        // OCCT L186: IndCoord = DefineClosedness(theFace).
        let ind_coord = define_closedness(the_face);

        if ind_coord != 0 {
            // OCCT L197-216: the max coordinate delta vertex.
            let mut a_max_delta = 0.0f64;
            for ii in 1..=a_ve_map.len() {
                let (a_vertex, a_elist) = &a_ve_map[ii - 1];
                let an_edge1 = &a_elist[0];
                let an_edge2 = &a_elist[a_elist.len() - 1];
                // OCCT L204-205: BRep_Tool::Parameter(aVertex, anEdge).
                let a_param1 = brep_tool_parameter(a_vertex, an_edge1);
                let a_param2 = brep_tool_parameter(a_vertex, an_edge2);
                // OCCT L206-209: BRepAdaptor_Curve2d Value(aParam).
                let a_ba_curve1 = BRepAdaptorCurve2d::new(an_edge1, the_face);
                let a_ba_curve2 = BRepAdaptorCurve2d::new(an_edge2, the_face);
                let a_pnt1 = a_ba_curve1.value(a_param1);
                let a_pnt2 = a_ba_curve2.value(a_param2);
                // OCCT L210: Coord(IndCoord) — 1 = X, 2 = Y.
                let a_delta = if ind_coord == 1 {
                    (a_pnt1.x - a_pnt2.x).abs()
                } else {
                    (a_pnt1.y - a_pnt2.y).abs()
                };
                if a_delta > a_max_delta {
                    a_max_delta = a_delta;
                    a_first_vertex = a_vertex.clone();
                }
            }
            // OCCT L217-229: the edge starting at aFirstVertex.
            if let Some((_, a_elist)) = a_ve_map.get(&shape_key(&a_first_vertex)) {
                for an_edge in a_elist {
                    let (a_v1, _a_v2) = top_exp_vertices(an_edge, true);
                    if let Some(a_v1) = a_v1 {
                        if a_v1.is_same(&a_first_vertex) {
                            a_first_edge = an_edge.clone();
                            break;
                        }
                    }
                }
            }
        }
    }

    // OCCT L233-265: the connected walk.
    let a_nb_edges = the_seq_edges.len();
    the_seq_edges.clear();
    the_seq_edges.push(a_first_edge.clone());
    let mut an_edge = a_first_edge;
    loop {
        // OCCT L239: TopExp::LastVertex(anEdge, true).
        let a_last_vertex = top_exp_last_vertex(&an_edge, true);
        if a_last_vertex.is_same(&a_first_vertex) {
            break;
        }

        let Some((_, a_elist)) = a_ve_map.get(&shape_key(&a_last_vertex)) else {
            break;
        };
        if a_elist.len() == 1 {
            break;
        }

        if a_elist[0].is_same(&an_edge) {
            an_edge = a_elist[a_elist.len() - 1].clone();
        } else {
            an_edge = a_elist[0].clone();
        }

        the_seq_edges.push(an_edge.clone());
        if the_seq_edges.len() == a_nb_edges {
            break;
        }
    }
}

/// OCCT Store(theEdge, theLV, theTol, IsToUpdate, theAsDes2d, theDMVV)
/// (BRepOffset_Inter2d.cxx L276-390) — stores the vertices into AsDes for the
/// edge despite coincidences; the coinciding chains are fused later by
/// FuseVertices().
pub(super) fn store(
    the_edge: &mut Shape,
    the_lv: &mut Vec<Shape>,
    the_tol: f64,
    is_to_update: bool,
    the_as_des2d: &mut BRepAlgoAsDes,
    the_dmvv: &mut DmvvMap,
) {
    // OCCT L285-291: update the vertices (BRep_Builder().UpdateVertex(aV,
    // theTol) — the shared-TShape mutation maps to the in-place rcad edit).
    for a_v in the_lv.iter_mut() {
        if a_v.shape_type() == ShapeType::Vertex {
            builder_update_vertex_tol(a_v, the_tol);
        }
    }

    // OCCT L294: the vertices already added to the edge (the rcad
    // borrow-split clone; the AsDes is mutated below).
    let a_lv_ex = the_as_des2d.descendant(the_edge).to_vec();
    // OCCT L295-302.
    if !is_to_update && a_lv_ex.is_empty() {
        if !the_lv.is_empty() {
            the_as_des2d.add_list(the_edge, the_lv);
        }
        return;
    }

    // OCCT L304-312: the projector for the vertex parameter on the edge.
    let mut a_proj_pc = GeomAPIProjectPointOnCurve::new();
    let mut a_tol_e = 0.0f64;
    if is_to_update {
        if let Some((a_c, a_t1, a_t2)) = brep_tool_curve(the_edge) {
            a_proj_pc.init_curve(&a_c, a_t1, a_t2);
        }
        a_tol_e = brep_tool_tolerance(the_edge);
    }

    // OCCT L314-389.
    let mut a_mv: HashSet<ShapeKey> = HashSet::new();
    for a_v in the_lv.iter_mut() {
        // OCCT L318-321: skip the duplicates within theLV.
        if !a_mv.insert(shape_key(a_v)) {
            continue;
        }

        // OCCT L323-324.
        let a_p = brep_tool_pnt(a_v).unwrap_or(DVec3::ZERO);
        let a_tol = brep_tool_tolerance(a_v);

        // OCCT L326-341: the already-added vertices within the combined
        // tolerance.
        let mut a_lvc: Vec<Shape> = Vec::new();
        let mut found_same = false;
        for a_v_ex in &a_lv_ex {
            if a_v.is_same(a_v_ex) {
                found_same = true;
                break;
            }
            let a_p_ex = brep_tool_pnt(a_v_ex).unwrap_or(DVec3::ZERO);
            let a_tol_v_ex = brep_tool_tolerance(a_v_ex);
            // OCCT L337: aP.IsEqual(aPEx, aTol + aTolVEx).
            if a_p.distance(a_p_ex) <= a_tol + a_tol_v_ex {
                a_lvc.push(a_v_ex.clone());
            }
        }
        // OCCT L343-346: the vertex is already stored on the edge.
        if found_same {
            continue;
        }

        // OCCT L348-365: the vertex parameter on the edge.
        if is_to_update {
            a_proj_pc.perform(a_p);
            if a_proj_pc.nb_points() == 0 {
                continue;
            }
            if a_proj_pc.lower_distance() > a_tol + a_tol_e {
                continue;
            }
            let a_t = a_proj_pc.lower_distance_parameter();
            // OCCT L363-364: UpdateVertex(V.Oriented(INTERNAL), aT,
            //             theEdge, aTol) — the orientation copy shares the
            // TShape; the rcad mutation targets a_v directly.
            builder_update_vertex_on_edge(a_v, a_t, the_edge, a_tol);
        }

        // OCCT L367-387: the fusion chains.
        if !a_lvc.is_empty() {
            for a_vc in &a_lvc {
                // OCCT L373-378: theDMVV.ChangeSeek/Add(aVC).Append(aV).
                match the_dmvv.get_mut(&shape_key(a_vc)) {
                    Some(entry) => entry.1.push(a_v.clone()),
                    None => {
                        the_dmvv.insert(shape_key(a_vc), (a_vc.clone(), vec![a_v.clone()]));
                    }
                }
            }
            // OCCT L381-386: theDMVV.ChangeSeek/Add(aV).Append(aLVC).
            match the_dmvv.get_mut(&shape_key(a_v)) {
                Some(entry) => entry.1.extend(a_lvc.iter().cloned()),
                None => {
                    the_dmvv.insert(shape_key(a_v), (a_v.clone(), a_lvc.clone()));
                }
            }
        }
        // OCCT L388.
        the_as_des2d.add(the_edge, a_v);
    }
}

/// OCCT Store(theE1, theE2, theLV1, theLV2, theTol, theAsDes2d, theDMVV)
/// (BRepOffset_Inter2d.cxx L396-412) — the two-edge overload (OCCT `!i ?
/// theE1 : theE2` ternary maps to the two arms).
fn store_pair(
    the_e1: &mut Shape,
    the_e2: &mut Shape,
    the_lv1: &mut Vec<Shape>,
    the_lv2: &mut Vec<Shape>,
    the_tol: f64,
    the_as_des2d: &mut BRepAlgoAsDes,
    the_dmvv: &mut DmvvMap,
) {
    for i in 0..2 {
        if i == 0 {
            store(the_e1, the_lv1, the_tol, false, the_as_des2d, the_dmvv);
        } else {
            store(the_e2, the_lv2, the_tol, false, the_as_des2d, the_dmvv);
        }
    }
}

/// OCCT EdgeInter(F, BAsurf, E1, E2, AsDes, Tol, WithOri, aDMVV)
/// (BRepOffset_Inter2d.cxx L416-699).
pub(super) fn edge_inter(
    f: &Shape,
    b_asurf: &BRepAdaptorSurface,
    e1: &Shape,
    e2: &Shape,
    as_des: &mut BRepAlgoAsDes,
    tol: f64,
    with_ori: bool,
    a_dmvv: &mut DmvvMap,
) {
    // OCCT L428-431.
    if e1.is_same(e2) {
        return;
    }

    // Architecture difference #22: the OCCT const& edge parameters are
    // mutated in place (BuildCurve3d / UpdateVertex); the rcad callee mutates
    // its own working copy.
    let mut e1 = e1.clone();
    let mut e2 = e2.clone();

    // OCCT L433-435: double f[3], l[3]; double TolDub = 1.e-7; int i;
    // (the OCCT `f`/`l` arrays; the face parameter `f` shadows the name).
    let mut f_range = [0.0f64; 3];
    let mut l_range = [0.0f64; 3];
    let tol_dub = 1.0e-7;

    // OCCT L437-438: BRep_Tool::Range(E1, f[1], l[1]) / Range(E2, f[2], l[2]).
    let (f1, l1) = brep_tool_range(&e1);
    f_range[1] = f1;
    l_range[1] = l1;
    let (f2, l2) = brep_tool_range(&e2);
    f_range[2] = f2;
    l_range[2] = l2;

    // OCCT L440-441: BRepAdaptor_Curve CE1(E1, F); CE2(E2, F).
    let ce1 = BRepAdaptorCurve::new(&e1, f);
    let ce2 = BRepAdaptorCurve::new(&e2, f);

    // OCCT L443-448: EI[1] = E1; EI[2] = E2; LV1; LV2; BRep_Builder B (the
    // rcad builder is the tool-leaf set).
    let ei0 = Shape::null();
    let ei: [&Shape; 3] = [&ei0, &e1, &e2];
    let mut lv1: Vec<Shape> = Vec::new();
    let mut lv2: Vec<Shape> = Vec::new();

    // OCCT L450-454: the common-vertex guard.
    if topexp_common_vertex(&e1, &e2).is_none() {
        // OCCT L453-454: BRepLib::BuildCurve3d(E1); BuildCurve3d(E2) — GAP
        // (architecture difference #25).
        brep_lib_build_curve3d(&e1, rcad_kernel::precision::CONFUSION);
        brep_lib_build_curve3d(&e2, rcad_kernel::precision::CONFUSION);

        // OCCT L456-457.
        let mut tol_sum = brep_tool_tolerance(&e1) + brep_tool_tolerance(&e2);
        tol_sum = tol_sum.max(1.0e-5);
        let _ = tol_sum;

        // OCCT L459-462: the solution accumulators; the degenerated point.
        let mut res_points: Vec<DVec3> = Vec::new();
        let mut res_params_on_e1: Vec<f64> = Vec::new();
        let mut res_params_on_e2: Vec<f64> = Vec::new();
        let mut deg_point = DVec3::ZERO;
        let with_degen = brep_tool_degenerated(&e1) || brep_tool_degenerated(&e2);

        // OCCT L464-478.
        if with_degen {
            let ideg = if brep_tool_degenerated(&e1) { 1usize } else { 2usize };
            // OCCT L467-472: TopoDS_Iterator iter(EI[ideg]) — the first
            // child vertex.
            let iter_children = sub_shapes(ei[ideg]);
            if !iter_children.is_empty() {
                let vdeg = &iter_children[0];
                deg_point = brep_tool_pnt(vdeg).unwrap_or(DVec3::ZERO);
            } else {
                // OCCT L476: CEdeg.Value(CEdeg.FirstParameter()).
                let cedeg = BRepAdaptorCurve::new(ei[ideg], f);
                deg_point = cedeg.value(cedeg.first_parameter());
            }
        }

        // OCCT L480-484: the pcurves (the out-parameters rewrite f/l) and
        // the 2D intersection.
        let Some((pcurve1, pf1, pl1)) = brep_tool_curve_on_surface(&e1, f) else {
            panic!("EdgeInter: E1 carries no pcurve on F (OCCT raises Standard_NoSuchObject)");
        };
        f_range[1] = pf1;
        l_range[1] = pl1;
        let Some((pcurve2, pf2, pl2)) = brep_tool_curve_on_surface(&e2, f) else {
            panic!("EdgeInter: E2 carries no pcurve on F (OCCT raises Standard_NoSuchObject)");
        };
        f_range[2] = pf2;
        l_range[2] = pl2;
        let inter2d = Geom2dIntGInter::new(&pcurve1, f_range[1], l_range[1], &pcurve2, f_range[2], l_range[2], tol_dub, tol_dub);

        // OCCT L485-500: collect the intersection points.
        for i in 1..=inter2d.nb_points() {
            let p3d = if with_degen {
                deg_point
            } else {
                let p2d = inter2d.point(i).value();
                b_asurf.value(p2d.x, p2d.y)
            };
            res_points.push(p3d);
            res_params_on_e1.push(inter2d.point(i).param_on_first());
            res_params_on_e2.push(inter2d.point(i).param_on_second());
        }

        // OCCT L502-589: the solution loop.
        for i in 1..=res_points.len() {
            let a_t1 = res_params_on_e1[i - 1];
            let a_t2 = res_params_on_e2[i - 1];
            // OCCT L506-512: the infinite-parameter rejection
            // (BRepOffset_Inter2d.cxx L506: Precision::IsInfinite(aT1/aT2)).
            if rcad_kernel::precision::is_infinite_value(a_t1)
                || rcad_kernel::precision::is_infinite_value(a_t2)
            {
                continue;
            }

            // OCCT L514-527: the new INTERNAL vertex with the max distance
            // tolerance.
            let p = res_points[i - 1];
            let mut a_new_vertex = brep_lib_make_vertex(p);
            a_new_vertex.orientation = Orientation::Internal;
            builder_update_vertex_on_edge(&mut a_new_vertex, a_t1, &mut e1, tol);
            builder_update_vertex_on_edge(&mut a_new_vertex, a_t2, &mut e2, tol);
            let p1 = ce1.value(a_t1);
            let p2 = ce2.value(a_t2);
            let dist1 = p1.distance(p);
            let dist2 = p2.distance(p);
            let dist3 = p1.distance(p2);
            let dist1 = dist1.max(dist2).max(dist3);
            builder_update_vertex_tol(&mut a_new_vertex, dist1);

            // OCCT L551-586: the orientation of the new vertex.
            let mut oo1 = Orientation::Reversed;
            let mut oo2 = Orientation::Reversed;
            if with_ori {
                // OCCT L556-561: the pcurve tangents (the BRepAdaptor_Curve2d
                // D1 — the rcad pcurve derivative).
                let v1 = pcurve1.derivative_at(a_t1);
                let v2 = pcurve2.derivative_at(a_t2);
                let mut v1or = v1;
                let mut v2or = v2;
                // OCCT L564-571.
                if e1.orientation == Orientation::Reversed {
                    v1or = -v1or;
                }
                if e2.orientation == Orientation::Reversed {
                    v2or = -v2or;
                }
                // OCCT L572-585.
                let mut cross_prod = cross2d(v2or, v1);
                if cross_prod > 0. {
                    oo1 = Orientation::Forward;
                }
                cross_prod = cross2d(v1or, v2);
                if cross_prod > 0. {
                    oo2 = Orientation::Forward;
                }
            }
            // OCCT L587-588.
            lv1.push(oriented(&a_new_vertex, oo1));
            lv2.push(oriented(&a_new_vertex, oo2));
        }
    }

    // OCCT L592-640: the end-vertex proximity test.
    let tol_conf = tol;
    let (v1_0, v1_1) = top_exp_vertices_raw(&e1);
    let (v2_0, v2_1) = top_exp_vertices_raw(&e2);
    let v1: [Option<&Shape>; 2] = [v1_0.as_ref(), v1_1.as_ref()];
    let v2: [Option<&Shape>; 2] = [v2_0.as_ref(), v2_1.as_ref()];

    for j in 0..2 {
        let Some(v1j) = v1[j] else { continue };
        for k in 0..2 {
            let Some(v2k) = v2[k] else { continue };
            // OCCT L614-620.
            if v1j.is_same(v2k) && as_des.has_ascendant(v1j) {
                continue;
            }
            // OCCT L622-625.
            let p1 = brep_tool_pnt(v1j).unwrap_or(DVec3::ZERO);
            let p2 = brep_tool_pnt(v2k).unwrap_or(DVec3::ZERO);
            let dist = p1.distance(p2);
            if dist < tol_conf {
                // OCCT L627-634.
                let a_tol = brep_tool_tolerance(v1j).max(brep_tool_tolerance(v2k));
                let mut v = brep_lib_make_vertex(p1);
                let u1 = if j == 0 { f_range[1] } else { l_range[1] };
                let u2 = if k == 0 { f_range[2] } else { l_range[2] };
                builder_update_vertex_on_edge(&mut v, u1, &mut e1, a_tol);
                builder_update_vertex_on_edge(&mut v, u2, &mut e2, a_tol);
                // OCCT L636-637: Prepend.
                lv1.insert(0, oriented(&v, v1j.orientation));
                lv2.insert(0, oriented(&v, v2k.orientation));
            }
        }
    }

    // OCCT L642-698: purge the doubles, then store.
    let affich_purge = false;

    if !lv1.is_empty() {
        // OCCT L650-691: remove all vertices — there can be doubles.
        let mut purge = true;
        while purge {
            let mut i = 1usize;
            purge = false;
            let mut pos = 0usize;
            while pos < lv1.len() {
                let a_v1 = &lv1[pos];
                let p1 = brep_tool_pnt(a_v1).unwrap_or(DVec3::ZERO);
                let a_tol1 = brep_tool_tolerance(a_v1);

                let mut j = 1usize;
                let mut jpos = 0usize;
                while j < i {
                    let a_v2 = &lv1[jpos];
                    let p2 = brep_tool_pnt(a_v2).unwrap_or(DVec3::ZERO);
                    let a_tol = a_tol1.max(brep_tool_tolerance(a_v2));
                    // OCCT L671: P1.IsEqual(P2, aTol).
                    if p1.distance(p2) <= a_tol {
                        lv1.remove(pos);
                        lv2.remove(pos);
                        if affich_purge {
                            // OCCT: std::cout << "Doubles removed in EdgeInter."
                        }
                        purge = true;
                        break;
                    }
                    j += 1;
                    jpos += 1;
                }
                if purge {
                    break;
                }
                i += 1;
                pos += 1;
            }
        }
        // OCCT L695-697: the storage tolerance and the two-edge Store.
        let tol_store = (brep_tool_tolerance(&e1) + brep_tool_tolerance(&e2)).max(tol);
        store_pair(
            &mut e1,
            &mut e2,
            &mut lv1,
            &mut lv2,
            tol_store,
            as_des,
            a_dmvv,
        );
    }
}

/// OCCT RefEdgeInter(F, BAsurf, E1, E2, theOr1, theOr2, AsDes, Tol, WithOri,
/// theVref, theImageVV, aDMVV, theCoincide) (BRepOffset_Inter2d.cxx
/// L703-1071).
pub(super) fn ref_edge_inter(
    f: &Shape,
    b_asurf: &BRepAdaptorSurface,
    e1: &Shape,
    e2: &Shape,
    the_or1: Orientation,
    the_or2: Orientation,
    as_des: &mut BRepAlgoAsDes,
    tol: f64,
    with_ori: bool,
    the_vref: &Shape,
    the_image_vv: &mut BRepAlgoImage,
    a_dmvv: &mut DmvvMap,
    the_coincide: &mut bool,
) {
    // OCCT L719-729.
    *the_coincide = false;
    if e1.is_same(e2) {
        return;
    }
    if e1.is_null() || e2.is_null() {
        return;
    }

    // Architecture difference #22: the working copies (see EdgeInter).
    let mut e1 = e1.clone();
    let mut e2 = e2.clone();

    // OCCT L731-733: double f[3], l[3]; TolDub = 1.e-7, TolLL = 0.0; int i;
    // (the OCCT `f`/`l` arrays; the face parameter `f` shadows the name).
    let mut f_range = [0.0f64; 3];
    let mut l_range = [0.0f64; 3];
    let tol_dub = 1.0e-7;
    let mut tol_ll = 0.0f64;

    // OCCT L735-740: the pcurves (the out-parameters rewrite f/l).
    let Some((pcurve1, pf1, pl1)) = brep_tool_curve_on_surface(&e1, f) else {
        return;
    };
    f_range[1] = pf1;
    l_range[1] = pl1;
    let Some((pcurve2, pf2, pl2)) = brep_tool_curve_on_surface(&e2, f) else {
        return;
    };
    f_range[2] = pf2;
    l_range[2] = pl2;

    // OCCT L742-750.
    let ce1 = BRepAdaptorCurve::new(&e1, f);
    let ce2 = BRepAdaptorCurve::new(&e2, f);
    let ei0 = Shape::null();
    let ei: [&Shape; 3] = [&ei0, &e1, &e2];
    let mut lv1: Vec<Shape> = Vec::new();
    let mut lv2: Vec<Shape> = Vec::new();

    // OCCT L752-753: BRepLib::BuildCurve3d(E1); BuildCurve3d(E2) — GAP
    // (architecture difference #25).
    brep_lib_build_curve3d(&e1, rcad_kernel::precision::CONFUSION);
    brep_lib_build_curve3d(&e2, rcad_kernel::precision::CONFUSION);

    // OCCT L755-774: the degenerated point.
    let mut res_points: Vec<DVec3> = Vec::new();
    let mut res_params_on_e1: Vec<f64> = Vec::new();
    let mut res_params_on_e2: Vec<f64> = Vec::new();
    let mut deg_point = DVec3::ZERO;
    let with_degen = brep_tool_degenerated(&e1) || brep_tool_degenerated(&e2);
    if with_degen {
        let ideg = if brep_tool_degenerated(&e1) { 1usize } else { 2usize };
        let iter_children = sub_shapes(ei[ideg]);
        if !iter_children.is_empty() {
            let vdeg = &iter_children[0];
            deg_point = brep_tool_pnt(vdeg).unwrap_or(DVec3::ZERO);
        } else {
            let cedeg = BRepAdaptorCurve::new(ei[ideg], f);
            deg_point = cedeg.value(cedeg.first_parameter());
        }
    }

    // OCCT L776-794: the coincident-line quick check.
    if matches!(pcurve1, Curve2d::Line(_)) && matches!(pcurve2, Curve2d::Line(_)) {
        let d1 = match &pcurve1 {
            Curve2d::Line(l) => l.direction,
            _ => unreachable!(),
        };
        let d2 = match &pcurve2 {
            Curve2d::Line(l) => l.direction,
            _ => unreachable!(),
        };
        // OCCT L781: |GAC1.Line().Direction().Angle(GAC2.Line().Direction())|.
        let an_angle = gp_dir2d_angle(d1, d2).abs();
        if an_angle <= 1.0e-8 || std::f64::consts::PI - an_angle <= 1.0e-8 {
            // OCCT L783-786.
            *the_coincide = true;
            return;
        } else {
            // OCCT L789-793: the line-line intersection range.
            tol_ll = int_tools_compute_int_range(tol_dub, tol_dub, an_angle);
            tol_ll = tol_ll.min(1.0e-5);
        }
    }

    // OCCT L796.
    let inter2d = Geom2dIntGInter::new(&pcurve1, f_range[1], l_range[1], &pcurve2, f_range[2], l_range[2], tol_dub, tol_dub);

    // OCCT L798-803.
    if !inter2d.is_done() || inter2d.nb_points() == 0 {
        *the_coincide = inter2d.nb_segments() > 0
            && matches!(pcurve1, Curve2d::Line(_))
            && matches!(pcurve2, Curve2d::Line(_));
        return;
    }

    // OCCT L805-820: collect the intersection points.
    for i in 1..=inter2d.nb_points() {
        let p3d = if with_degen {
            deg_point
        } else {
            let p2d = inter2d.point(i).value();
            b_asurf.value(p2d.x, p2d.y)
        };
        res_points.push(p3d);
        res_params_on_e1.push(inter2d.point(i).param_on_first());
        res_params_on_e2.push(inter2d.point(i).param_on_second());
    }

    // OCCT L822-919: the solution loop.
    for i in 1..=res_points.len() {
        let a_t1 = res_params_on_e1[i - 1];
        let a_t2 = res_params_on_e2[i - 1];
        // OCCT L826-832 (BRepOffset_Inter2d.cxx L826:
        // Precision::IsInfinite(aT1/aT2)).
        if rcad_kernel::precision::is_infinite_value(a_t1)
            || rcad_kernel::precision::is_infinite_value(a_t2)
        {
            continue;
        }

        // OCCT L834-847.
        let p = res_points[i - 1];
        let mut a_new_vertex = brep_lib_make_vertex(p);
        a_new_vertex.orientation = Orientation::Internal;
        builder_update_vertex_on_edge(&mut a_new_vertex, a_t1, &mut e1, tol);
        builder_update_vertex_on_edge(&mut a_new_vertex, a_t2, &mut e2, tol);
        let p1 = ce1.value(a_t1);
        let p2 = ce2.value(a_t2);
        let dist1 = p1.distance(p);
        let dist2 = p2.distance(p);
        let dist3 = p1.distance(p2);
        let dist1 = dist1.max(dist2).max(dist3);
        builder_update_vertex_tol(&mut a_new_vertex, dist1);

        // OCCT L871-906: the orientation of the new vertex.
        let mut oo1 = Orientation::Reversed;
        let mut oo2 = Orientation::Reversed;
        if with_ori {
            let v1 = pcurve1.derivative_at(a_t1);
            let v2 = pcurve2.derivative_at(a_t2);
            let mut v1or = v1;
            let mut v2or = v2;
            if e1.orientation == Orientation::Reversed {
                v1or = -v1or;
            }
            if e2.orientation == Orientation::Reversed {
                v2or = -v2or;
            }
            let mut cross_prod = cross2d(v2or, v1);
            if cross_prod > 0. {
                oo1 = Orientation::Forward;
            }
            cross_prod = cross2d(v1or, v2);
            if cross_prod > 0. {
                oo2 = Orientation::Forward;
            }
        }

        // OCCT L908-915: the forced orientations.
        if the_or1 != Orientation::External {
            oo1 = the_or1;
        }
        if the_or2 != Orientation::External {
            oo2 = the_or2;
        }

        // OCCT L917-918.
        lv1.push(oriented(&a_new_vertex, oo1));
        lv2.push(oriented(&a_new_vertex, oo2));
    }

    // OCCT L921-966: the end-vertex proximity test.
    let tol_conf = tol;
    let (v1_0, v1_1) = top_exp_vertices_raw(&e1);
    let (v2_0, v2_1) = top_exp_vertices_raw(&e2);
    let v1: [Option<&Shape>; 2] = [v1_0.as_ref(), v1_1.as_ref()];
    let v2: [Option<&Shape>; 2] = [v2_0.as_ref(), v2_1.as_ref()];

    for j in 0..2 {
        let Some(v1j) = v1[j] else { continue };
        for k in 0..2 {
            let Some(v2k) = v2[k] else { continue };
            if v1j.is_same(v2k) && as_des.has_ascendant(v1j) {
                continue;
            }
            let p1 = brep_tool_pnt(v1j).unwrap_or(DVec3::ZERO);
            let p2 = brep_tool_pnt(v2k).unwrap_or(DVec3::ZERO);
            let dist = p1.distance(p2);
            if dist < tol_conf {
                // OCCT L956-963: NOTE the plain Tol (no aTol computation in
                // the OCCT RefEdgeInter).
                let mut v = brep_lib_make_vertex(p1);
                let u1 = if j == 0 { f_range[1] } else { l_range[1] };
                let u2 = if k == 0 { f_range[2] } else { l_range[2] };
                builder_update_vertex_on_edge(&mut v, u1, &mut e1, tol);
                builder_update_vertex_on_edge(&mut v, u2, &mut e2, tol);
                lv1.insert(0, oriented(&v, v1j.orientation));
                lv2.insert(0, oriented(&v, v2k.orientation));
            }
        }
    }

    let affich_purge = false;

    if !lv1.is_empty() {
        // OCCT L976-1013: purge the doubles (the Tol comparison).
        let mut purge = true;
        while purge {
            let mut i = 1usize;
            purge = false;
            let mut pos = 0usize;
            while pos < lv1.len() {
                let mut j = 1usize;
                let mut jpos = 0usize;
                while j < i {
                    // OCCT L991-993: the points are read inside the inner
                    // loop.
                    let p1 = brep_tool_pnt(&lv1[pos]).unwrap_or(DVec3::ZERO);
                    let p2 = brep_tool_pnt(&lv1[jpos]).unwrap_or(DVec3::ZERO);
                    if p1.distance(p2) <= tol {
                        lv1.remove(pos);
                        lv2.remove(pos);
                        if affich_purge {
                            // OCCT: std::cout << "Doubles removed in EdgeInter."
                        }
                        purge = true;
                        break;
                    }
                    j += 1;
                    jpos += 1;
                }
                if purge {
                    break;
                }
                i += 1;
                pos += 1;
            }
        }

        // OCCT L1018-1047: keep only the vertex nearest theVref.
        if lv1.len() > 1 {
            let pref = brep_tool_pnt(the_vref).unwrap_or(DVec3::ZERO);
            let mut dmin = f64::MAX; // OCCT RealLast().
            let mut vmin: Option<Shape> = None;
            for it in &lv1 {
                let p = brep_tool_pnt(it).unwrap_or(DVec3::ZERO);
                let d = p.distance_squared(pref);
                if d < dmin {
                    dmin = d;
                    vmin = Some(it.clone());
                }
            }
            if let Some(vmin) = &vmin {
                let mut pos = 0usize;
                while pos < lv1.len() {
                    if !vmin.is_same(&lv1[pos]) {
                        lv1.remove(pos);
                        lv2.remove(pos);
                        if pos >= lv1.len() {
                            break;
                        }
                    } else {
                        pos += 1;
                    }
                }
            }
        }

        // OCCT L1049-1062: bind the surviving vertices as images of theVref.
        for it in &lv1 {
            let mut a_new_vertex = it.clone();
            a_new_vertex.orientation = Orientation::Forward;
            if the_image_vv.has_image(the_vref) {
                the_image_vv.add(&oriented(the_vref, Orientation::Forward), &a_new_vertex);
            } else {
                the_image_vv.bind(&oriented(the_vref, Orientation::Forward), &a_new_vertex);
            }
        }

        // OCCT L1065-1069: the storage tolerance (vs the Line-Line tolerance)
        // and the two-edge Store.
        let tol_store = (brep_tool_tolerance(&e1) + brep_tool_tolerance(&e2)).max(tol);
        let tol_store = tol_store.max(tol_ll);
        store_pair(
            &mut e1,
            &mut e2,
            &mut lv1,
            &mut lv2,
            tol_store,
            as_des,
            a_dmvv,
        );
    }
}

/// OCCT evaluateMaxSegment(aCurveOnSurface) (BRepOffset_Inter2d.cxx
/// L1078-1096) — the MaxSegment to pass in approximation.
pub(super) fn evaluate_max_segment(a_curve_on_surface: &Adaptor3dCurveOnSurface) -> i32 {
    let a_surf = a_curve_on_surface.get_surface();
    let a_curv2d = a_curve_on_surface.get_curve();

    // OCCT L1083-1089: the surface knot counts (the rcad expanded knot
    // vectors map to the distinct-knot count).
    let mut a_nb_s_knots = 0.0f64;
    if let Surface3::BSpline(a_bspline) = a_surf {
        let nb_u = distinct_knots_1d(&a_bspline.knots_u) as f64;
        let nb_v = distinct_knots_1d(&a_bspline.knots_v) as f64;
        a_nb_s_knots = nb_u.max(nb_v);
    }
    // OCCT L1090-1093: the pcurve knot count.
    let mut a_nb_c2d_knots = 0.0f64;
    if let Curve2d::BSpline(a_bspline) = a_curv2d {
        a_nb_c2d_knots = distinct_knots_1d(&a_bspline.knots) as f64;
    }
    // OCCT L1094-1095.
    (30.0 + a_nb_s_knots.max(a_nb_c2d_knots)) as i32
}

/// OCCT ExtendPCurve(aPCurve, anEf, anEl, a2Offset, NewPCurve)
/// (BRepOffset_Inter2d.cxx L1100-1187) — prolongs the bounded pcurve by
/// tangent segments; the result lands in `new_pcurve`.
pub(super) fn extend_pcurve(
    a_pcurve: &Curve2d,
    an_ef: f64,
    an_el: f64,
    a2_offset: f64,
    new_pcurve: &mut Curve2d,
) -> bool {
    // OCCT L1106.
    *new_pcurve = a_pcurve.clone();
    // OCCT L1107-1110: unwrap the trimmed curve to its basis.
    if let Curve2d::Trimmed(tc) = new_pcurve {
        let basis = tc.curve.clone(); // Box<Curve2d>
        *new_pcurve = *basis;
    }

    // OCCT L1112-1113.
    let mut first_par = new_pcurve.default_domain()[0];
    let mut last_par = new_pcurve.default_domain()[1];

    // OCCT L1115-1140: the two-pole line prolongation for bounded curves.
    let is_bounded = matches!(
        new_pcurve,
        Curve2d::BSpline(_) | Curve2d::Bezier(_) | Curve2d::Trimmed(_)
    );
    if is_bounded && (first_par > an_ef - a2_offset || last_par < an_el + a2_offset) {
        if let Curve2d::Bezier(a_bezier) = new_pcurve {
            if a_bezier.control_points.len() == 2 {
                // OCCT L1123-1126: the poles and the line through them.
                let p1 = a_bezier.control_points[0];
                let p2 = a_bezier.control_points[1];
                let a_vec = p2 - p1;
                *new_pcurve = Curve2d::Line(Line2d::new(p1, a_vec));
                return true;
            }
        } else if let Curve2d::BSpline(a_bspline) = new_pcurve {
            if distinct_knots_1d(&a_bspline.knots) == 2 && a_bspline.control_points.len() == 2 {
                // OCCT L1134-1137.
                let p1 = a_bspline.control_points[0];
                let p2 = a_bspline.control_points[1];
                let a_vec = p2 - p1;
                *new_pcurve = Curve2d::Line(Line2d::new(p1, a_vec));
                return true;
            }
        }
    }

    // OCCT L1142-1144: the trimmed copy of the ORIGINAL pcurve.
    first_par = a_pcurve.default_domain()[0];
    last_par = a_pcurve.default_domain()[1];
    let a_tr_curve = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(a_pcurve.clone()),
        t_min: first_par,
        t_max: last_par,
    });

    // OCCT L1146-1156: the comp-curve engine (GAP — architecture difference
    // #27).
    let mut a_comp_curve = Geom2dConvertCompCurveToBSplineCurve::new(&a_tr_curve);
    let a_tol = rcad_kernel::precision::CONFUSION;
    let a_delta = a2_offset.max(1.);

    // OCCT L1158-1170: the begin prolongation.
    if first_par > an_ef - a2_offset {
        let a_p_bnd = a_pcurve.point_at(first_par);
        let a_v_bnd = a_pcurve.derivative_at(first_par);
        let a_d_bnd = a_v_bnd.normalize_or_zero();
        let a_p_beg = a_p_bnd - a_delta * a_d_bnd;
        let a_lin = Curve2d::Line(Line2d::new(a_p_beg, a_d_bnd));
        let a_segment = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(a_lin),
            t_min: 0.0,
            t_max: a_delta,
        });

        if !a_comp_curve.add(&a_segment, a_tol) {
            return false;
        }
    }

    // OCCT L1172-1183: the end prolongation.
    if last_par < an_el + a2_offset {
        let a_p_beg = a_pcurve.point_at(last_par);
        let a_v_bnd = a_pcurve.derivative_at(last_par);
        let a_d_bnd = a_v_bnd.normalize_or_zero();
        let a_lin = Curve2d::Line(Line2d::new(a_p_beg, a_d_bnd));
        let a_segment = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(a_lin),
            t_min: 0.0,
            t_max: a_delta,
        });

        if !a_comp_curve.add(&a_segment, a_tol) {
            return false;
        }
    }

    // OCCT L1185-1186.
    *new_pcurve = a_comp_curve.bspline_curve();
    true
}
