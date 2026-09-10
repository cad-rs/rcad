// OCCT BRepOffset_MakeOffset.cxx L1-5659 + BRepOffset_MakeOffset.hxx L28-285 —
// 1:1 translation (split by OCCT order into brep_offset_make_offset.rs /
// brep_offset_make_offset_b.rs / brep_offset_make_offset_c.rs /
// brep_offset_make_offset_d.rs).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_MakeOffset.cxx / .hxx
//
// Module a carries cxx L143-L931 (the static helpers of the first half, the
// anonymous-namespace progress machinery, the class data + small methods) and
// the shared architecture-difference / GAP-carrier section.  Module b carries
// cxx L837-L2394 (MakeOffsetShape / MakeThickSolid / MakeOffsetFaces /
// BuildOffsetByInter / BuildOffsetByArc / SelfInter / ToContext /
// UpdateFaceOffset).  Module c carries cxx L2396-L4750 +
// L5231-L5659 class methods (CorrectConicalFaces, Intersection3D/2D,
// MakeLoops, MakeFaces, MakeMissingWalls, MakeShells, MakeSolid,
// SelectShells, EncodeRegularity, CheckInputData, RemoveInternalEdges,
// IntersectEdges, Generated/Modified/IsDeleted, analyzeProgress, IsPlanar).
// Module d carries the remaining file statics (UpdateInitOffset,
// ComputeMaxDist, UpdateTolerance, CorrectSolid, checkSinglePoint,
// RemoveShapes, UpdateHistory, TrimEdges, TrimEdge, GetEnlargedFaces,
// BuildShellsCompleteInter, GetSubShapes, RemoveSeamAndDegeneratedEdges,
// IsSolid, AppendToList).
//
// OCCT inheritance chain (hxx L65): none — BRepOffset_MakeOffset is a
// standalone class.  The OCCT 8.0 BRepOffset_MakeOffset_1.cxx is not a class
// header — it hosts the definitions of the private methods BuildSplitsOf*
// (cxx_1 L9494-9533) through the local BRepOffset_BuildOffsetFaces engine;
// the rcad carrier is super::brep_offset_make_offset_1::BRepOffsetMakeOffset1
// (field my_make_offset_1) holding that engine.
//
// Architecture differences (numbering continues from brep_offset_tool_d.rs):
// 38. NCollection_List<TopoDS_Shape> -> Vec<Shape>; NCollection_Map ->
//     OcctShapeSet; NCollection_IndexedMap -> OcctIndexedShapeMap;
//     NCollection_DataMap -> ShapeDataMap<V>; NCollection_IndexedDataMap ->
//     ShapeIndexedDataMap<V>; NCollection_Sequence -> Vec<Shape>;
//     NCollection_Array1<double>(0, Last-1) -> Vec<f64> (the OCCT array is
//     0-based here, the index mapping is preserved).
// 39. Handle(BRepAlgo_AsDes) -> the owned BRepAlgoAsDes.  The OCCT Handle
//     aliasing between myAsDes and the BRepOffset_Inter3d / Inter2d locals is
//     result-equal in rcad by a boundary resync: the Inter takes the AsDes by
//     value (its constructor form) and myAsDes is re-assigned from the Inter
//     at the end of the OCCT aliasing scope (annotated at the sites).
// 40. Message_ProgressScope / Message_ProgressRange -> the rcad NoopProgress
//     + ProgressScope forms; the OCCT child-scope chaining
//     (aPS.Next(aSteps(...))) is flattened to local scopes (the progress
//     machinery is result-inert; arch. diff. #37 of Inter3d).  The OCCT
//     `!aPS.More()` user-break checks map onto ProgressScope::user_break()
//     (the rcad NoopProgress never breaks — the OCCT no-indicator run has
//     the same behavior).  OSD_Chronometer + the OCCT_DEBUG blocks
//     (DEBVerticesControl cxx L143-191, AffichInt2d/ChronBuild tracing) are
//     debug-only and omitted (arch. diff. #31 of bi_tgte_blended.rs).
// 41. BRep_Builder -> BRepBuilder over the class-local BRep pool (arch.
//     diff. #4/#19); the file statics that construct containers take the
//     pool as an explicit `brep: &mut BRep` parameter (the OCCT local
//     `BRep_Builder B;` is stateless; the pool is the rcad arena stand-in).
// 42. BRepOffset_Analyse (TKOffset/BRepOffset, its own translation unit) ->
//     the BRepOffsetAnalyse of brep_offset_tool_d.rs (the accessor surface
//     consumed by Inter3d/Inter2d) plus the local GAP free functions below
//     for the member surface only MakeOffset consumes (Perform / Edges /
//     HasGenerated / AddFaces / SetOffsetValue / SetFaceOffsetMap).
// 43. BRepOffset_MakeLoops (TKOffset/BRepOffset) — GAP carrier (the local
//     BRepOffsetMakeLoops below; the Build / BuildOnContext / BuildFaces
//     leaves panic until the unit lands).
// 44. BRepTools_Quilt (TKBRep/BRepTools) — the real body lives in
//     crate::topalgo::brep_tools_quilt (the Glue form, re-exported below).
// 45. BRepBuilderAPI_Sewing — the BRepBuilderAPISewing GAP carrier of
//     bi_tgte_blended.rs (arch. diff. #24).
// 46. BRepLib::{BuildCurves3d, SameParameter, UpdateTolerances, SortFaces,
//     BuildCurve3d} — GAP static leaves (the bi_tgte #25 / inter2d #26
//     precedents; the TopoDS-form BRepLib::SameParameter(E, Tol, DoAlsoMinmax)
//     and BRepLib::BuildCurve3d(E, Tol) keep the OCCT failure-path structure).
// 47. BRepTools::{UVBounds, Update, IsReallyClosed} — GAP static leaves.
// 48. BRepTools_WireExplorer (TKTopAlgo/BRepTools) — GAP carrier.
// 49. BRepCheck_Analyzer / BRepCheck_Edge / BRepCheck_Vertex (TKTopAlgo/
//     BRepCheck) — GAP carriers (the CorrectSolid / MakeSolid /
//     UpdateTolerance validity probes keep the OCCT branch structure).
// 50. BRepClass3d_SolidClassifier -> the crate::topalgo::brep_class3d
//     SolidClassifier translation (perform / state).
// 51. BRepGProp::VolumeProperties + GProp_GProps::Mass — GAP leaf (the
//     brep_offset_offset.rs #17 precedent carries LinearProperties; the
//     volume re-host of the kernel base::gprop is keyed by the BRep pool,
//     not by a TopoDS_Shape, so the carrier keeps the OCCT zero-mass failure
//     path until the TopoDS form lands).
// 52. BRepLib_FindSurface / GeomFill_Generator /
//     IntTools_FClass2d / BOPAlgo_MakerVolume / BOPTools_AlgoTools::
//     MakeSplitEdge (the TopoDS form; the rcad bop::tools carrier is the
//     BOPDS-index form) / GeomAPI_ProjectPointOnCurve (the (P, C) form) /
//     GeomLib::BuildCurve3d (the Adaptor3d_CurveOnSurface form) /
//     GC_MakeCylindricalSurface / GC_MakeLine2d / gce_MakeCone / gce_MakeDir
//     / GeomFill_Generator / Adaptor3d_CurveOnSurface +
//     Geom2dAdaptor_Curve + GeomAdaptor_Surface — GAP carriers/leaves; each
//     keeps the OCCT call form and the failure path at its call site.
//     (GeomLib_IsPlanarSurface left this list: the real body of
//     crate::geomalgo::geom_lib_is_planar_surface is re-exported below; its
//     remaining iso-curve GAP leaves live in that unit.)
// 53. gp_Vec / gp_Pnt / gp_Dir arithmetic -> DVec3; gp_Circ / gp_Cone /
//     gp_Pln / gp_Sphere / gp_Ax1..3 -> the rcad_kernel::geom carriers
//     (Circle3 / ConicalSurface / Plane / SphericalSurface; the Ax3
//     position form folds into the surface struct fields).
// 54. BRepAdaptor_Surface / BRepAdaptor_Curve / BRepAdaptor_Curve2d -> the
//     local BRepAdaptorSurface / BRepAdaptorCurve / BRepAdaptorCurve2d
//     re-hosts below (the inter2d.rs precedents, extended by the consumed
//     GetType / Circle / Plane / D1 surface).
// 55. The OCCT `myFaceOffset` iteration order (NCollection_DataMap bucket
//     order) maps to the HashMap iteration order — the OCCT order is itself
//     unspecified (UpdateFaceOffset / SetFacesWithOffset / CheckInputData
//     are order-insensitive in OCCT: the Bind/UnBind results are
//     order-independent sets of the same final bindings).
// 56. TopoDS_Iterator -> bat::sub_shapes (the direct-children walk with the
//     composed orientation); TopExp::Vertices(E, V1, V2, CumOri = true) ->
//     the local top_exp_vertices_cum_ori re-host.

use std::collections::HashMap;

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, BRepTool, Orientation, ShapeType};
use rcad_kernel::topo::topods::curve_on_surface_pool_free as brep_tool_curve_on_surface_uv;
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_offset_b::BRepOffsetOffset;
use super::brep_offset_tool::{
    oriented, shape_data_map, top_exp_vertices, OcctIndexedShapeMap, OcctShapeSet, ShapeDataMap,
};
use super::brep_offset_tool_d::BRepOffsetAnalyse;

use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_tolerance,
};
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::tool as bat;
use crate::brep_fill::offset_wire::GeomAbsJoinType;
use crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity;

// ---------------------------------------------------------------------------
// Shared map forms (architecture difference #38).
// ---------------------------------------------------------------------------

/// OCCT NCollection_DataMap<TopoDS_Shape, BRepOffset_Offset> (theMapSF).
pub(crate) type MapSF = ShapeDataMap<BRepOffsetOffset>;

/// OCCT NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> (MES / Build /
/// Created / MEF / ShapeTgt / theETrimEInf / aFacesOrigins).
pub(crate) type DataMapOfShapeShape = ShapeDataMap<Shape>;

/// OCCT NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>
/// (anEdgesOrigins / myEdgeIntEdges form).
pub(crate) type DataMapOfShapeListOfShape = ShapeDataMap<Vec<Shape>>;

/// OCCT NCollection_DataMap<TopoDS_Shape, double> (myFaceOffset).
pub(crate) type DataMapOfShapeReal = ShapeDataMap<f64>;

/// OCCT NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>
/// (Contours / aDMVV / aDMELF / aDMEF forms).
pub(crate) type IndexedDataMapOfShapeListOfShape = indexmap::IndexMap<
    crate::feat::loc_ope_wires_on_shape_b::ShapeKey,
    (Shape, Vec<Shape>),
>;

/// OCCT Precision::Infinite().
pub(crate) const PRECISION_INFINITE: f64 = f64::INFINITY;
/// OCCT RealLast().
pub(crate) const REAL_LAST: f64 = f64::MAX;
/// OCCT RealFirst().
pub(crate) const REAL_FIRST: f64 = f64::MIN;
/// OCCT M_PI.
pub(crate) const OCCT_PI: f64 = std::f64::consts::PI;
/// OCCT M_PI_4.
pub(crate) const OCCT_PI_4: f64 = std::f64::consts::PI / 4.0;

// ===========================================================================
// GAP carriers / leaves (architecture differences #42-#54).
// ===========================================================================

/// OCCT BRepOffset_Analyse::Perform(S, Tol) — the R1 real body
/// (brep_offset_analyse.rs perform, BRepOffset_Analyse.cxx L200-380); the
/// free-function bridge keeps the MakeOffset call forms unchanged.
pub(crate) fn analyse_perform(a: &mut BRepOffsetAnalyse, s: &Shape, tol: f64) {
    a.perform(s, tol);
}

/// OCCT BRepOffset_Analyse::SetOffsetValue(Offset).
pub(crate) fn analyse_set_offset_value(a: &mut BRepOffsetAnalyse, offset: f64) {
    a.set_offset_value(offset);
}

/// OCCT BRepOffset_Analyse::SetFaceOffsetMap(Map).
pub(crate) fn analyse_set_face_offset_map(a: &mut BRepOffsetAnalyse, m: &DataMapOfShapeReal) {
    a.set_face_offset_map(m);
}

/// OCCT BRepOffset_Analyse::Edges(S, TC, L) — the tangent-edge collection;
/// dispatches on the shape kind (vertex -> EdgesOnVertex, else ->
/// EdgesOnFace), the E0 re-host split.
pub(crate) fn analyse_edges(
    a: &BRepOffsetAnalyse,
    s: &Shape,
    tc: ChFiDS_TypeOfConcavity,
    l: &mut Vec<Shape>,
) {
    if s.shape_type() == ShapeType::Vertex {
        a.edges_on_vertex(s, tc, l);
    } else {
        a.edges_on_face(s, tc, l);
    }
}

/// OCCT BRepOffset_Analyse::HasGenerated(S).
pub(crate) fn analyse_has_generated(a: &BRepOffsetAnalyse, s: &Shape) -> bool {
    a.has_generated(s)
}

/// OCCT BRepOffset_Analyse::AddFaces(F, Co, Dummy, TC) — the no-RT form.
pub(crate) fn analyse_add_faces(
    a: &BRepOffsetAnalyse,
    f: &Shape,
    co: &mut Shape,
    dummy: &mut OcctShapeSet,
    tc: ChFiDS_TypeOfConcavity,
) {
    a.add_faces(f, co, dummy, tc);
}

/// OCCT BRepOffset_Analyse::AddFaces(F, Co, Dummy, TC, RT) — the RT form.
pub(crate) fn analyse_add_faces_rt(
    a: &BRepOffsetAnalyse,
    f: &Shape,
    co: &mut Shape,
    dummy: &mut OcctShapeSet,
    tc: ChFiDS_TypeOfConcavity,
    rt: ChFiDS_TypeOfConcavity,
) {
    a.add_faces_2_types(f, co, dummy, tc, rt);
}

// OCCT BRepOffset_MakeLoops (TKOffset/BRepOffset/BRepOffset_MakeLoops.hxx /
// .cxx) — the real body lives in brep_offset_make_offset_loops.rs (the
// 2000-line module split; architecture difference #43); the class type
// keeps the OCCT import path through the re-export.
pub(crate) use super::brep_offset_make_offset_loops::BRepOffsetMakeLoops;

/// OCCT BRepTools_Quilt (TKBRep/BRepTools/BRepTools_Quilt.hxx / .cxx) — the
/// real body lives in crate::topalgo::brep_tools_quilt (architecture
/// difference #44); the Glue form used by IsConnectedShell / MakeThickSolid
/// / MakeShells keeps the OCCT import path through the re-export.
pub(crate) use crate::topalgo::brep_tools_quilt::BRepToolsQuilt;

/// OCCT BRepLib::SortFaces(S, LF) — the real body lives in
/// crate::topalgo::brep_lib::BRepLib::sort_faces (OCCT BRepLib.cxx
/// L2894-2953); this keeps the offset-module import path under the same
/// signature.
pub(crate) fn brep_lib_sort_faces(s: &Shape, lf: &mut Vec<Shape>) {
    crate::topalgo::brep_lib::brep_lib::BRepLib::sort_faces(s, lf);
}

/// OCCT BRepLib::BuildCurves3d(S, Tol) (TKTopAlgo/BRepLib, BRepLib.cxx
/// L468-489; BRepLib.hxx L99-103 defaults Continuity=GeomAbs_C1,
/// MaxDegree=14, MaxSegment=0) — the real body lives in
/// crate::topalgo::brep_lib::BRepLib::build_curves3d_tol; the brep
/// argument carries the rcad pool mutation arena (architecture
/// difference #4).
pub(crate) fn brep_lib_build_curves3d_tol(
    the_brep: &mut BRep,
    the_s: &Shape,
    the_tol: f64,
) -> bool {
    crate::topalgo::brep_lib::brep_lib::BRepLib::build_curves3d_tol(the_brep, the_s, the_tol)
}

/// OCCT BRepLib::BuildCurve3d(E, Tol) (BRepLib.cxx L301-456; BRepLib.hxx
/// L90-94 defaults Continuity=GeomAbs_C1, MaxDegree=14, MaxSegment=0) —
/// the same re-host for the edge form: rebuilds the 3D curve of the edge
/// from its pcurves.
pub(crate) fn brep_lib_build_curve3d_edge(
    the_brep: &mut BRep,
    the_e: &Shape,
    the_tol: f64,
) -> bool {
    crate::topalgo::brep_lib::brep_lib::BRepLib::build_curve3d(
        the_brep,
        the_e,
        the_tol,
        rcad_kernel::topo::topods::GeomAbsShape::C1,
        14,
        0,
    )
}

/// OCCT BRepLib::SameParameter(E, Tol, DoAlsoMinmax) — GAP static leaf
/// (architecture difference #46; the inter2d.rs no-op form is the (E, Tol)
/// pair; the MakeOffset call sites use the three-arg form whose OCCT
/// default DoAlsoMinmax = false, so the leaf is a no-op there too).
pub(crate) fn brep_lib_same_parameter_3(_e: &Shape, _tol: f64) {}

/// OCCT BRepLib::UpdateTolerances(S, Strict = false) — GAP static leaf
/// (architecture difference #46).
pub(crate) fn brep_lib_update_tolerances(_s: &mut Shape) {}

/// OCCT BRep_Tool::CurveOnSurface(E, F) (BRep_Tool.cxx L339-401) — the
/// pool-free offset form: the edge's representations are matched by
/// (surface value, L.Predivided(E.Location())) — the owning-face pointer
/// OCCT BRepTools::AddUVBounds(aF, aE, aB) (BRepTools.cxx L181-365).
fn brep_tools_add_uv_bounds_edge(
    the_brep: &BRep,
    a_f: &Shape,
    a_e: &Shape,
    a_b: &mut BndBox2d,
) {
    //
    // OCCT L186-193: aC2D = BRep_Tool::CurveOnSurface(aE, aF, aT1, aT2).
    let (a_c2d, a_t1, a_t2) = match brep_tool_curve_on_surface_uv(a_e, a_f) {
        Some(t) => t,
        None => return,
    };
    // OCCT L196: BndLib_Add2dCurve::Add(aC2D, aT1, aT2, 0., aBoxC).
    let a_box_c = rcad_kernel::curve2d_bounding_box(&a_c2d, a_t1, a_t2, 0.0);
    let mut a_xmin = 0.0f64;
    let mut a_ymin = 0.0f64;
    let mut a_xmax = 0.0f64;
    let mut a_ymax = 0.0f64;
    let a_box_c_void = !(a_box_c[0] <= a_box_c[1] && a_box_c[2] <= a_box_c[3]);
    if !a_box_c_void {
        a_xmin = a_box_c[0];
        a_xmax = a_box_c[1];
        a_ymin = a_box_c[2];
        a_ymax = a_box_c[3];
    }
    //
    // OCCT L199-201: aS = BRep_Tool::Surface(aF); aS->Bounds(aUmin, aUmax, aVmin, aVmax).
    let a_s = match bat::brep_tool_surface(a_f) {
        Some(s) => s,
        None => return,
    };
    let (a_umin, a_umax, a_vmin, a_vmax) = {
        let d = a_s.default_domain();
        (d[0], d[1], d[2], d[3])
    };
    //
    // OCCT L204-209: unwrap the RectangularTrimmedSurface basis.
    let a_s_basis = match &a_s {
        Surface3::Trimmed(t) => (*t.basis).clone(),
        _ => a_s.clone(),
    };
    //
    // OCCT L212+: if (!aS->IsUPeriodic()).
    if !a_s_basis.is_u_periodic() {
        let mut is_u_periodic = false;
        // OCCT L219-221: BSpline verification when the box exceeds bounds.
        if matches!(a_s_basis, Surface3::BSpline(_)) && (a_xmin < a_umin || a_xmax > a_umax) {
            let a_tol2 = 100.0
                * rcad_kernel::core::precision::CONFUSION
                * rcad_kernel::core::precision::CONFUSION;
            is_u_periodic = true;
            // 1. Verify that the surface is U-closed (L226-238).
            if !a_s_basis.is_u_closed() {
                let a_v_step = a_vmax - a_vmin;
                let mut a_v = a_vmin;
                while a_v <= a_vmax {
                    let p1 = a_s_basis.point_at(a_umin, a_v);
                    let p2 = a_s_basis.point_at(a_umax, a_v);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_u_periodic = false;
                        break;
                    }
                    a_v += a_v_step;
                }
            }
            // 2. Verify periodicity of surface inside UV-bounds (L240-306).
            if is_u_periodic {
                let a_v = (a_vmin + a_vmax) * 0.5;
                let mut a_u = [0.0f64; 6];
                let mut a_upp = [0.0f64; 6];
                let mut a_nb_pnt = 0usize;
                if a_xmin < a_umin {
                    a_u[0] = a_xmin;
                    a_u[1] = (a_xmin + a_umin) * 0.5;
                    a_u[2] = a_umin;
                    a_upp[0] = a_u[0] + a_umax - a_umin;
                    a_upp[1] = a_u[1] + a_umax - a_umin;
                    a_upp[2] = a_u[2] + a_umax - a_umin;
                    a_nb_pnt += 3;
                }
                if a_xmax > a_umax {
                    a_u[a_nb_pnt] = a_umax;
                    a_u[a_nb_pnt + 1] = (a_xmax + a_umax) * 0.5;
                    a_u[a_nb_pnt + 2] = a_xmax;
                    a_upp[a_nb_pnt] = a_u[a_nb_pnt] - a_umax + a_umin;
                    a_upp[a_nb_pnt + 1] = a_u[a_nb_pnt + 1] - a_umax + a_umin;
                    a_upp[a_nb_pnt + 2] = a_u[a_nb_pnt + 2] - a_umax + a_umin;
                    a_nb_pnt += 3;
                }
                for an_ind in 0..a_nb_pnt {
                    let p1 = a_s_basis.point_at(a_u[an_ind], a_v);
                    let p2 = a_s_basis.point_at(a_upp[an_ind], a_v);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_u_periodic = false;
                        break;
                    }
                }
            }
        }
        // OCCT L309-320: the clamp.
        if !is_u_periodic {
            if (a_xmin < a_umin) && (a_umin < a_xmax) {
                a_xmin = a_umin;
            }
            if (a_xmin < a_umax) && (a_umax < a_xmax) {
                a_xmax = a_umax;
            }
        }
    }
    //
    // OCCT L323+: if (!aS->IsVPeriodic()).
    if !a_s_basis.is_v_periodic() {
        let mut is_v_periodic = false;
        if matches!(a_s_basis, Surface3::BSpline(_)) && (a_ymin < a_vmin || a_ymax > a_vmax) {
            let a_tol2 = 100.0
                * rcad_kernel::core::precision::CONFUSION
                * rcad_kernel::core::precision::CONFUSION;
            is_v_periodic = true;
            if !a_s_basis.is_v_closed() {
                let a_u_step = a_umax - a_umin;
                let mut a_u = a_umin;
                while a_u <= a_umax {
                    let p1 = a_s_basis.point_at(a_u, a_vmin);
                    let p2 = a_s_basis.point_at(a_u, a_vmax);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_v_periodic = false;
                        break;
                    }
                    a_u += a_u_step;
                }
            }
            if is_v_periodic {
                let a_u = (a_umin + a_umax) * 0.5;
                let mut a_v = [0.0f64; 6];
                let mut a_vpp = [0.0f64; 6];
                let mut a_nb_pnt = 0usize;
                if a_ymin < a_vmin {
                    a_v[0] = a_ymin;
                    a_v[1] = (a_ymin + a_vmin) * 0.5;
                    a_v[2] = a_vmin;
                    a_vpp[0] = a_v[0] + a_vmax - a_vmin;
                    a_vpp[1] = a_v[1] + a_vmax - a_vmin;
                    a_vpp[2] = a_v[2] + a_vmax - a_vmin;
                    a_nb_pnt += 3;
                }
                if a_ymax > a_vmax {
                    a_v[a_nb_pnt] = a_vmax;
                    a_v[a_nb_pnt + 1] = (a_ymax + a_vmax) * 0.5;
                    a_v[a_nb_pnt + 2] = a_ymax;
                    a_vpp[a_nb_pnt] = a_v[a_nb_pnt] - a_vmax + a_vmin;
                    a_vpp[a_nb_pnt + 1] = a_v[a_nb_pnt + 1] - a_vmax + a_vmin;
                    a_vpp[a_nb_pnt + 2] = a_v[a_nb_pnt + 2] - a_vmax + a_vmin;
                    a_nb_pnt += 3;
                }
                for an_ind in 0..a_nb_pnt {
                    let p1 = a_s_basis.point_at(a_u, a_v[an_ind]);
                    let p2 = a_s_basis.point_at(a_u, a_vpp[an_ind]);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_v_periodic = false;
                        break;
                    }
                }
            }
        }
        if !is_v_periodic {
            if (a_ymin < a_vmin) && (a_vmin < a_ymax) {
                a_ymin = a_vmin;
            }
            if (a_ymin < a_vmax) && (a_vmax < a_ymax) {
                a_ymax = a_vmax;
            }
        }
    }
    // OCCT L362-366: aBoxS.Update(aXmin, aYmin, aXmax, aYmax); aB.Add(aBoxS).
    let mut a_box_s = BndBox2d::new();
    a_box_s.update(a_xmin, a_ymin, a_xmax, a_ymax);
    a_b.add_box(&a_box_s);
}

/// OCCT BRepTools::AddUVBounds(FF, B) (BRepTools.cxx L125-159) — the face
/// form.
fn brep_tools_add_uv_bounds_face(the_brep: &BRep, a_ff: &Shape, a_b: &mut BndBox2d) {
    // OCCT L128-129: F = FF oriented FORWARD.
    let a_f = bat::oriented(a_ff, Orientation::Forward);
    // OCCT L132-136: fill the box for the given face.
    let mut a_box = BndBox2d::new();
    for ex in bat::explorer(&a_f, ShapeType::Edge, ShapeType::Shape) {
        brep_tools_add_uv_bounds_edge(the_brep, &a_f, &ex, &mut a_box);
    }
    // OCCT L139-154: if the box is empty, get natural bounds.
    if a_box.is_void() {
        let a_surf = match bat::brep_tool_surface(&a_f) {
            Some(s) => s,
            None => return,
        };
        let (u_min, u_max, v_min, v_max) = {
            let d = a_surf.default_domain();
            // OCCT L150-151: aSurf->Bounds(UMin, UMax, VMin, VMax) — the OCCT
            // infinite surfaces (Geom_Plane::Bounds L181-184) return finite
            // Precision::Infinite() magnitudes (2e100); the rcad
            // default_domain carries true infinities, so map to the OCCT
            // finite form (architecture difference #61).
            fn to_occt_infinite(v: f64) -> f64 {
                if v == f64::INFINITY {
                    rcad_kernel::core::precision::INFINITE_VALUE
                } else if v == f64::NEG_INFINITY {
                    -rcad_kernel::core::precision::INFINITE_VALUE
                } else {
                    v
                }
            }
            (
                to_occt_infinite(d[0]),
                to_occt_infinite(d[1]),
                to_occt_infinite(d[2]),
                to_occt_infinite(d[3]),
            )
        };
        a_box.update(u_min, v_min, u_max, v_max);
    }
    // OCCT L157: add the face box to the result.
    a_b.add_box(&a_box);
}

/// OCCT BRepTools::UVBounds(F, Umin, Umax, Vmin, Vmax) (BRepTools.cxx L64-80).
/// (the_chfi3d_builder_cncrn.rs (brep, face) form; the pcurve resolution goes
/// through the BRepTool::CurveOnSurface surface-value matcher.)
pub(crate) fn brep_tools_uv_bounds(the_brep: &BRep, _f: &Shape) -> (f64, f64, f64, f64) {
    // OCCT L70-71: Bnd_Box2d B; AddUVBounds(F, B).
    let mut a_b = BndBox2d::new();
    brep_tools_add_uv_bounds_face(the_brep, _f, &mut a_b);
    // OCCT L72-77: if (!B.IsVoid()) B.Get(UMin, VMin, UMax, VMax); else all 0.
    if let Some((u_min, v_min, u_max, v_max)) = a_b.get() {
        (u_min, u_max, v_min, v_max)
    } else {
        (0.0, 0.0, 0.0, 0.0)
    }
}

/// OCCT BRepTools::Update(F) — GAP static leaf (architecture difference
/// #47): recomputes the pcurves/tolerances of the face.
pub(crate) fn brep_tools_update(_f: &mut Shape) {}

/// OCCT BRepTools::IsReallyClosed(E, F) — GAP static leaf (architecture
/// difference #47): true when E is a seam edge of F.
pub(crate) fn brep_tools_is_really_closed(_e: &Shape, _f: &Shape) -> bool {
    panic!("GAP: BRepTools::IsReallyClosed (TKTopAlgo/BRepTools not translated)");
}

/// OCCT BRepTools_WireExplorer (TKTopAlgo/BRepTools) — GAP carrier
/// (architecture difference #48): the ordered edge walk of a wire on a
/// face.
pub(crate) struct BRepToolsWireExplorer;

impl BRepToolsWireExplorer {
    /// OCCT BRepTools_WireExplorer::Init(W, F).
    pub fn init(&mut self, _w: &Shape, _f: &Shape) {
        panic!("GAP: BRepTools_WireExplorer::Init (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_WireExplorer::More().
    pub fn more(&self) -> bool {
        panic!("GAP: BRepTools_WireExplorer::More (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_WireExplorer::Current().
    pub fn current(&self) -> Shape {
        panic!("GAP: BRepTools_WireExplorer::Current (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_WireExplorer::Next().
    pub fn next(&mut self) {
        panic!("GAP: BRepTools_WireExplorer::Next (TKTopAlgo/BRepTools not translated)");
    }
}

/// OCCT BRepCheck_Analyzer(S, GeomChecks) (TKTopAlgo/BRepCheck) — GAP
/// carrier (architecture difference #49): the overall validity probe.
pub(crate) struct BRepCheckAnalyzer;

impl BRepCheckAnalyzer {
    /// OCCT BRepCheck_Analyzer::BRepCheck_Analyzer(S, GeomChecks = false).
    pub fn new(_s: &Shape, _geom_checks: bool) -> Self {
        BRepCheckAnalyzer
    }

    /// OCCT BRepCheck_Analyzer::IsValid().
    pub fn is_valid(&self) -> bool {
        panic!("GAP: BRepCheck_Analyzer::IsValid (TKTopAlgo/BRepCheck not translated)");
    }
}

/// OCCT BRepCheck_Edge(E)::Tolerance() — GAP leaf (architecture difference
/// #49): the edge validity tolerance.
pub(crate) fn brep_check_edge_tolerance(_e: &Shape) -> f64 {
    panic!("GAP: BRepCheck_Edge::Tolerance (TKTopAlgo/BRepCheck not translated)");
}

/// OCCT BRepCheck_Vertex(V)::Tolerance() — GAP leaf (architecture
/// difference #49): the vertex validity tolerance.
pub(crate) fn brep_check_vertex_tolerance(_v: &Shape) -> f64 {
    panic!("GAP: BRepCheck_Vertex::Tolerance (TKTopAlgo/BRepCheck not translated)");
}

/// OCCT BRepGProp::VolumeProperties(S, VProps, OnlyClosed) + GProp_GProps::
/// Mass() — GAP leaf (architecture difference #51).
pub(crate) fn brep_gprop_volume_properties(_s: &Shape) -> f64 {
    panic!("GAP: BRepGProp::VolumeProperties (TKTopAlgo/BRepGProp not translated)");
}

/// OCCT BRepLib_FindSurface (TKTopAlgo/BRepLib_FindSurface.hxx / .cxx) — the
/// real body lives in crate::topalgo::brep_lib_find_surface (architecture
/// difference #52); the plane finder of MakeMissingWalls keeps the OCCT
/// import path through the re-export.
pub(crate) use crate::topalgo::brep_lib_find_surface::BRepLibFindSurface;

/// OCCT GeomLib_IsPlanarSurface (TKGeomBase/GeomLib/GeomLib_IsPlanarSurface.hxx
/// / .cxx) — the real body lives in
/// crate::geomalgo::geom_lib_is_planar_surface (architecture difference
/// #52); the IsPlanar consumption keeps the OCCT import path through the
/// re-export.
pub(crate) use crate::geomalgo::geom_lib_is_planar_surface::GeomLibIsPlanarSurface;

/// OCCT GeomFill_Generator (TKGeomAlgo/GeomFill) — GAP carrier
/// (architecture difference #52): the two-section ruled-surface generator
/// of MakeMissingWalls (the OCCT Extrusion path).
pub(crate) struct GeomFillGenerator;

impl GeomFillGenerator {
    /// OCCT GeomFill_Generator::AddCurve(C).
    pub fn add_curve(&mut self, _c: &Curve3) {
        panic!("GAP: GeomFill_Generator::AddCurve (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Generator::Perform(Tol3d).
    pub fn perform(&mut self, _tol3d: f64) {
        panic!("GAP: GeomFill_Generator::Perform (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Generator::Surface().
    pub fn surface(&self) -> Surface3 {
        panic!("GAP: GeomFill_Generator::Surface (TKGeomAlgo/GeomFill not translated)");
    }
}

// OCCT IntTools_FClass2d — the real body lives in
// crate::bop::int_tools::int_tools_fclass2d (the C.3 carrier switch; the
// local panic carrier is deleted).

/// OCCT BOPAlgo_MakerVolume (TKBO/BOPAlgo) — GAP carrier (architecture
/// difference #52): the MakeVolume engine of BuildShellsCompleteInter.
pub(crate) struct BOPAlgoMakerVolume;

impl BOPAlgoMakerVolume {
    /// OCCT BOPAlgo_Options::SetArguments(theLS).
    pub fn set_arguments(&mut self, _ls: &[Shape]) {}

    /// OCCT BOPAlgo_Options::SetIntersect(flag).
    pub fn set_intersect(&mut self, _flag: bool) {}

    /// OCCT BOPAlgo_MakerVolume::SetAvoidInternalShapes(flag).
    pub fn set_avoid_internal_shapes(&mut self, _flag: bool) {}

    /// OCCT BOPAlgo_Options::Perform(theRange).
    pub fn perform(&mut self) {
        panic!("GAP: BOPAlgo_MakerVolume::Perform (TKBO/BOPAlgo not translated)");
    }

    /// OCCT BOPAlgo_Options::HasErrors().
    pub fn has_errors(&self) -> bool {
        panic!("GAP: BOPAlgo_MakerVolume::HasErrors (TKBO/BOPAlgo not translated)");
    }

    /// OCCT BOPAlgo_Builder::Modified(S) — the history hook of
    /// UpdateHistory.
    pub fn modified(&self, _s: &Shape) -> Vec<Shape> {
        panic!("GAP: BOPAlgo_Builder::Modified (TKBO/BOPAlgo not translated)");
    }

    /// OCCT BOPAlgo_Options::Shape().
    pub fn shape(&self) -> Shape {
        panic!("GAP: BOPAlgo_MakerVolume::Shape (TKBO/BOPAlgo not translated)");
    }
}

/// OCCT TopoDS_Shape::IsNull() for the myOffsetShape carrier — true when the
/// shape carries no real TShape.  The rcad pool-free builder products (arena
/// index = usize::MAX) still carry a REAL Compound/Shell/Solid TShape payload
/// (e.g. BRepTools_Quilt::Shells() — the deferred compound assembly); the
/// OCCT null handle is only the Shape::null() dummy payload.
pub(crate) fn offset_shape_is_null(the_s: &Shape) -> bool {
    if !the_s.is_null() {
        return false;
    }
    !matches!(
        the_s.data.as_ref(),
        rcad_kernel::topo::topods::TShape::Compound(_)
            | rcad_kernel::topo::topods::TShape::Shell(_)
            | rcad_kernel::topo::topods::TShape::Solid(_)
    )
}

/// OCCT BOPTools_AlgoTools::MakeSplitEdge(NE, V1, aT1, V2, aT2, aSourceEdge)
/// (BOPTools_AlgoTools_2.cxx L138-183) — the empty-copied FORWARD edge with
/// the oriented extremity vertices and the ordered range.
pub(crate) fn bop_algo_tools_make_split_edge(
    ne: &Shape,
    v1: &Shape,
    t1: f64,
    v2: &Shape,
    t2: f64,
    source_edge: &mut Shape,
) {
    // OCCT L143-144: E = TopoDS::Edge(aE.Oriented(TopAbs_FORWARD));
    // E.EmptyCopy().
    let mut e_fwd = ne.clone();
    e_fwd.orientation = Orientation::Forward;
    *source_edge = bat::empty_copied(&e_fwd);
    // OCCT L146-166: the oriented extremity Adds (V1 FORWARD when aP1 < aP2,
    // REVERSED otherwise; V2 REVERSED when aP1 < aP2, FORWARD otherwise).
    if let rcad_kernel::topo::topods::TShape::Edge(ed) =
        std::sync::Arc::make_mut(&mut source_edge.data)
    {
        if !v1.is_null() {
            let o1 = if t1 < t2 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            let ov1 = bat::oriented(v1, o1);
            ed.my_shapes.push(ov1.clone());
            if o1 == Orientation::Forward {
                ed.first = ov1;
            } else {
                ed.last = ov1;
            }
        }
        if !v2.is_null() {
            let o2 = if t1 < t2 {
                Orientation::Reversed
            } else {
                Orientation::Forward
            };
            let ov2 = bat::oriented(v2, o2);
            ed.my_shapes.push(ov2.clone());
            if o2 == Orientation::Forward {
                ed.first = ov2;
            } else {
                ed.last = ov2;
            }
        }
        // OCCT L168-174: BB.Range(E, aP1, aP2) / BB.Range(E, aP2, aP1).
        ed.range = if t1 < t2 { [t1, t2] } else { [t2, t1] };
    }
    // OCCT L176-177: aNewEdge = E; aNewEdge.Orientation(aE.Orientation()).
    source_edge.orientation = ne.orientation;
}

/// OCCT GeomAPI_ProjectPointOnCurve (TKTopAlgo/GeomAPI; the (P, C) form) —
/// GAP carrier (architecture difference #52; the brep_offset_offset.rs #12
/// precedent).
pub(crate) struct GeomAPIProjectPointOnCurve;

impl GeomAPIProjectPointOnCurve {
    /// OCCT GeomAPI_ProjectPointOnCurve::GeomAPI_ProjectPointOnCurve(P, C).
    pub fn new(_p: DVec3, _c: &Curve3) -> Self {
        GeomAPIProjectPointOnCurve
    }

    /// OCCT GeomAPI_ProjectPointOnCurve::NbPoints().
    pub fn nb_points(&self) -> i32 {
        panic!("GAP: GeomAPI_ProjectPointOnCurve::NbPoints (TKTopAlgo/GeomAPI not translated)");
    }

    /// OCCT GeomAPI_ProjectPointOnCurve::LowerDistanceParameter().
    pub fn lower_distance_parameter(&self) -> f64 {
        panic!(
            "GAP: GeomAPI_ProjectPointOnCurve::LowerDistanceParameter (TKTopAlgo/GeomAPI not translated)"
        );
    }
}

/// OCCT GeomLib::BuildCurve3d(Tol, ConS, FirstPar, LastPar, C3d,
/// MaxDeviation, AverageDeviation) (the Adaptor3d_CurveOnSurface form) —
/// GAP leaf (architecture difference #52; the inter2d.rs #29 form).
pub(crate) fn geom_lib_build_curve3d_cons(
    _the_c3d: &mut Option<Curve3>,
    _the_max_deviation: &mut f64,
    _the_average_deviation: &mut f64,
) {
    panic!("GAP: GeomLib::BuildCurve3d (TKTopAlgo/GeomLib not translated)");
}

/// OCCT Adaptor3d_CurveOnSurface(HC2d, HSurf) — GAP carrier (architecture
/// difference #52): the pcurve-on-surface adapter consumed only by
/// GeomLib::BuildCurve3d.
pub(crate) struct Adaptor3dCurveOnSurface;

impl Adaptor3dCurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface::Adaptor3d_CurveOnSurface(Curve2d,
    /// Surface).
    pub fn new(_c2d: &Curve2d, _s: &Surface3) -> Self {
        Adaptor3dCurveOnSurface
    }
}

/// OCCT GC_MakeCylindricalSurface(Circ).Value() — GAP leaf (architecture
/// difference #52): the cylinder of the circle.
pub(crate) fn gc_make_cylindrical_surface(_c: &rcad_kernel::geom::Circle3) -> Surface3 {
    panic!("GAP: GC_MakeCylindricalSurface (TKGeomBase/GC not translated)");
}

/// OCCT GC_MakeLine2d(P1, P2).Value() — GAP leaf (architecture difference
/// #52): the 2D line through two points.
pub(crate) fn gc_make_line2d(_p1: DVec2, _p2: DVec2) -> Curve2d {
    panic!("GAP: GC_MakeLine2d (TKGeomBase/GC not translated)");
}

/// OCCT gce_MakeDir(P1, P2).Value() — the unit direction P1 -> P2; the
/// OCCT gce error (null distance) is the caller's NotDone path (architecture
/// difference #52).
pub(crate) fn gce_make_dir(_p1: DVec3, _p2: DVec3) -> DVec3 {
    panic!("GAP: gce_MakeDir (TKMath/gce not translated)");
}

/// OCCT gce_MakeCone(P1, P2, R1, R2).Value() — GAP leaf (architecture
/// difference #52): the cone through two circles.
pub(crate) fn gce_make_cone(
    _p1: DVec3,
    _p2: DVec3,
    _r1: f64,
    _r2: f64,
) -> rcad_kernel::geom::ConicalSurface {
    panic!("GAP: gce_MakeCone (TKMath/gce not translated)");
}

/// OCCT BRepLib::SameParameter(S) — the face form — GAP static leaf
/// (architecture difference #46; the (E, Tol) edge form keeps the OCCT
/// failure-path structure).
pub(crate) fn brep_lib_same_parameter_face(_the_s: &mut Shape) {}

/// OCCT BRepLib_MakeFace(W, OnlyPlane = true) — GAP leaf (architecture
/// difference #52): the only-plane face maker of MakeMissingWalls.
pub(crate) fn brep_lib_make_face_wire_only_plane(_the_w: &Shape) -> Shape {
    panic!("GAP: BRepLib_MakeFace(W, OnlyPlane) (TKTopAlgo/BRepLib not translated)");
}

// ---------------------------------------------------------------------------
// TopExp / BRep_Tool re-hosts local to the MakeOffset cluster.
// ---------------------------------------------------------------------------

/// OCCT TopExp::MapShapes(S, T, M) (TopExp.cxx L481-511) — the unique-shape
/// IndexedMap collection (the rcad explorer order is the OCCT order).
pub(crate) fn top_exp_map_shapes(s: &Shape, t: ShapeType) -> OcctIndexedShapeMap {
    let mut m = OcctIndexedShapeMap::new();
    for c in bat::explorer(s, t, ShapeType::Shape) {
        m.add(&c);
    }
    m
}

/// OCCT TopExp::MapShapesAndAncestors(S, TS, TA, M) (TopExp.cxx L360-398)
/// — the ancestor map child -> parents (one entry per occurrence).
pub(crate) fn top_exp_map_shapes_and_ancestors(
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut IndexedDataMapOfShapeListOfShape,
) {
    for anc in bat::explorer(s, ta, ShapeType::Shape) {
        for child in bat::explorer(&anc, ts, ShapeType::Shape) {
            let key = crate::brep_algo::tool::shape_key(&child);
            let entry = m.entry(key).or_insert((child.clone(), Vec::new()));
            entry.1.push(anc.clone());
        }
    }
}

/// OCCT TopExp::MapShapesAndUniqueAncestors(S, TS, TA, M) (TopExp.cxx
/// L400-440) — the ancestor map where each ancestor is listed at most
/// once per child.
pub(crate) fn top_exp_map_shapes_and_unique_ancestors(
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut IndexedDataMapOfShapeListOfShape,
) {
    // OCCT L406-436: explore the ancestors, then the children of each.
    let mut built: IndexedDataMapOfShapeListOfShape = indexmap::IndexMap::new();
    for anc in bat::explorer(s, ta, ShapeType::Shape) {
        for child in bat::explorer(&anc, ts, ShapeType::Shape) {
            let key = crate::brep_algo::tool::shape_key(&child);
            let entry = built
                .entry(key)
                .or_insert((child.clone(), Vec::new()));
            if !entry.1.iter().any(|a| a.is_same(&anc)) {
                entry.1.push(anc.clone());
            }
        }
    }
    *m = built;
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri = true) (TopExp.cxx
/// L255-272): the extremities with the cumulated (edge) orientation —
/// Vfirst FORWARD, Vlast REVERSED, swapped when the edge is REVERSED.
pub(crate) fn top_exp_vertices_cum_ori(edg: &Shape) -> (Shape, Shape) {
    let (mut vfirst, mut vlast) = top_exp_vertices(edg);
    if edg.orientation == Orientation::Reversed {
        std::mem::swap(&mut vfirst, &mut vlast);
    }
    vfirst = oriented(&vfirst, Orientation::Forward);
    vlast = oriented(&vlast, Orientation::Reversed);
    (vfirst, vlast)
}

/// OCCT TopExp::FirstVertex(E) (TopExp.cxx L283-296) — the
/// CumOri-true first extremity.
pub(crate) fn top_exp_first_vertex(edg: &Shape) -> Shape {
    let (v, _) = top_exp_vertices_cum_ori(edg);
    v
}

// ===========================================================================
// Anonymous namespace (cxx L194-243).
// ===========================================================================

/// OCCT BRepOffset_PIOperation (cxx L194-206) — the operation list of the
/// Progress Indicator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub(crate) enum BRepOffset_PIOperation {
    PIOperation_CheckInputData = 0,
    PIOperation_Analyse,
    PIOperation_BuildOffsetBy,
    PIOperation_Intersection,
    PIOperation_MakeMissingWalls,
    PIOperation_MakeShells,
    PIOperation_MakeSolid,
    PIOperation_Sewing,
    PIOperation_Last,
}

/// OCCT normalizeSteps (cxx L227-243) — normalization of progress steps.
pub(crate) fn normalize_steps(the_whole: f64, the_steps: &mut [f64]) {
    let mut a_sum = 0.;
    for i in 0..the_steps.len() {
        a_sum += the_steps[i];
    }

    // Normalize steps
    for i in 0..the_steps.len() {
        the_steps[i] = the_whole * the_steps[i] / a_sum;
    }
}

// ===========================================================================
// File statics, first half.
// ===========================================================================

/// OCCT OrientationOfEdgeInFace (cxx L326-343).
pub(crate) fn orientation_of_edge_in_face(the_edge: &Shape, the_face: &Shape) -> Orientation {
    let mut an_or = Orientation::External;

    for an_edge in bat::explorer(the_face, ShapeType::Edge, ShapeType::Shape) {
        if an_edge.is_same(the_edge) {
            an_or = an_edge.orientation;
            break;
        }
    }

    an_or
}

/// OCCT FindParameter (cxx L346-512).
pub(crate) fn find_parameter(v: &Shape, e: &Shape, u: &mut f64) -> bool {
    // Search the vertex in the edge

    let mut rev = false;
    let mut vf = Shape::null();
    let mut orient = Orientation::Internal;

    let itv = bat::sub_shapes(&oriented(e, Orientation::Forward));

    // if the edge has no vertices
    // and is degenerated use the vertex orientation
    // RLE, june 94

    if itv.is_empty() && brep_tool_degenerated(e) {
        orient = v.orientation;
    }

    for vcur in &itv {
        if v.is_same(vcur) {
            if vf.is_null() {
                vf = vcur.clone();
            } else {
                rev = e.orientation == Orientation::Reversed;
                if vcur.orientation == v.orientation {
                    vf = vcur.clone();
                }
            }
        }
    }

    if !vf.is_null() {
        orient = vf.orientation;
    }

    if orient == Orientation::Forward {
        let (f, l) = bat::brep_tool_range(e);
        // return (rev) ? l : f;
        *u = if rev { l } else { f };
        return true;
    } else if orient == Orientation::Reversed {
        let (f, l) = bat::brep_tool_range(e);
        // return (rev) ? f : l;
        *u = if rev { f } else { l };
        return true;
    } else {
        // OCCT L419-421: L = L.Predivided(V.Location()) — the
        // identity-location convention (architecture difference #32 of
        // inter2d.rs).
        let c3d = bat::brep_tool_curve(e);
        if c3d.is_some() || brep_tool_degenerated(e) {
            // OCCT L425-490: the BRep_PointRepresentation walk of the
            // vertex (the rcad TVertexData::points list; the
            // IsPointOnCurve identity maps to the curve-index match of the
            // same-pool representation).
            let prs = match v.data.as_ref() {
                rcad_kernel::topo::topods::TShape::Vertex(vd) => vd.points.clone(),
                _ => Vec::new(),
            };
            for pr in &prs {
                match pr {
                    rcad_kernel::topo::topods::PointRepresentation::PointOnCurve {
                        curve,
                        parameter,
                        ..
                    } => {
                        let mut p = *parameter;
                        // Closed curves RLE 16 june 94
                        if let Some((c, f, l)) = &c3d {
                            let _ = curve; // identity: the same-pool index form
                            if f.is_infinite() && *f < 0. {
                                *u = *parameter;
                                return true;
                            }
                            if l.is_infinite() && *l > 0. {
                                *u = *parameter;
                                return true;
                            }
                            let pf = c.point_at(*f);
                            let pl = c.point_at(*l);
                            let tol = brep_tool_tolerance(v);
                            if pf.distance(pl) < tol {
                                if let Some(vp) = bat::brep_tool_pnt(v) {
                                    if pf.distance(vp) < tol {
                                        if v.orientation == Orientation::Forward {
                                            p = *f; // p = f;
                                        } else {
                                            p = *l; // p = l;
                                        }
                                    }
                                }
                            }
                        }
                        // return res;//p;
                        *u = p;
                        return true;
                    }
                    // OCCT L459-491: the first-pcurve fallback
                    // (pr->IsPointOnCurveOnSurface(PC, S, L)) — the rcad
                    // PointRepresentation carries no CurveOnSurface variant,
                    // so the OCCT non-match falls through to the next
                    // representation (architecture difference #56).
                    rcad_kernel::topo::topods::PointRepresentation::PointOnSurface { .. } => {}
                }
            }
        } else {
            // no 3d curve !!
            // let us try with the first pcurve
            // OCCT L462-491: the PC walk — the rcad PointRepresentation
            // carries no CurveOnSurface variant, so no representation can
            // match IsPointOnCurveOnSurface and the OCCT loop body never
            // fires (architecture difference #56).
            let _ = brep_tool_curve_on_surface(e, &Shape::null());
        }
    }

    // throw Standard_NoSuchObject("BRep_Tool:: no parameter on edge");
    false
}

/// OCCT GetEdgePoints (cxx L515-533).
pub(crate) fn get_edge_points(
    an_edge: &Shape,
    a_face: &Shape,
    fpnt: &mut DVec3,
    mpnt: &mut DVec3,
    lpnt: &mut DVec3,
) {
    let (the_curve, f, l) = brep_tool_curve_on_surface(an_edge, a_face)
        .expect("BRep_Tool::CurveOnSurface null");
    let fpnt2d = the_curve.point_at(f);
    let lpnt2d = the_curve.point_at(l);
    let mpnt2d = the_curve.point_at(0.5 * (f + l));
    let a_surf = bat::brep_tool_surface(a_face).expect("BRep_Tool::Surface null");
    *fpnt = a_surf.point_at(fpnt2d.x, fpnt2d.y);
    *lpnt = a_surf.point_at(lpnt2d.x, lpnt2d.y);
    *mpnt = a_surf.point_at(mpnt2d.x, mpnt2d.y);
}

/// OCCT FillContours (cxx L536-604).
pub(crate) fn fill_contours(
    a_shape: &Shape,
    analyser: &BRepOffsetAnalyse,
    contours: &mut IndexedDataMapOfShapeListOfShape,
    map_ef: &mut DataMapOfShapeShape,
) {
    let mut edges: Vec<Shape> = Vec::new();

    for explo in bat::explorer(a_shape, ShapeType::Face, ShapeType::Shape) {
        let a_face = explo;
        for itf in bat::sub_shapes(&a_face) {
            let a_wire = itf;
            // OCCT L552: BRepTools_WireExplorer Wexp (GAP carrier,
            // architecture difference #48).
            let mut wexp = BRepToolsWireExplorer;
            wexp.init(&a_wire, &a_face);
            while wexp.more() {
                let an_edge = wexp.current();
                if brep_tool_degenerated(&an_edge) {
                    wexp.next();
                    continue;
                }
                let lint = analyser.type_(&an_edge);
                if !lint.is_empty() && lint[0].type_of() == ChFiDS_TypeOfConcavity::FreeBound {
                    shape_data_map::bind(map_ef, &an_edge, a_face.clone());
                    edges.push(an_edge);
                }
                wexp.next();
            }
        }
    }

    while !edges.is_empty() {
        let start_edge = edges.remove(0);
        let (start_vertex, mut cur_vertex) = top_exp_vertices_cum_ori(&start_edge);
        let mut a_contour: Vec<Shape> = Vec::new();
        a_contour.push(start_edge);
        while !cur_vertex.is_same(&start_vertex) {
            for i in 0..edges.len() {
                let an_edge = edges[i].clone();
                let (v1, v2) = top_exp_vertices(&an_edge);
                if v1.is_same(&cur_vertex) || v2.is_same(&cur_vertex) {
                    a_contour.push(an_edge);
                    cur_vertex = if v1.is_same(&cur_vertex) { v2 } else { v1 };
                    edges.remove(i);
                    break;
                }
            }
        }
        shape_indexed_data_map_add(contours, &start_vertex, a_contour);
    }
}

/// OCCT IndexedDataMap::Add(K, V) on the
/// IndexedDataMapOfShapeListOfShape form.
pub(crate) fn shape_indexed_data_map_add(
    m: &mut IndexedDataMapOfShapeListOfShape,
    k: &Shape,
    v: Vec<Shape>,
) {
    let key = crate::brep_algo::tool::shape_key(k);
    if !m.contains_key(&key) {
        m.insert(key, (k.clone(), v));
    }
}

/// OCCT IndexedDataMap::FindKey(i) on the
/// IndexedDataMapOfShapeListOfShape form (1-based).
pub(crate) fn shape_indexed_data_map_find_key_1(
    m: &IndexedDataMapOfShapeListOfShape,
    i: usize,
) -> Shape {
    m.get_index(i - 1).expect("IndexedDataMap::FindKey out of range").1 .0.clone()
}

/// OCCT IndexedDataMap::operator(i) on the
/// IndexedDataMapOfShapeListOfShape form (1-based).
pub(crate) fn shape_indexed_data_map_value_1(
    m: &IndexedDataMapOfShapeListOfShape,
    i: usize,
) -> &Vec<Shape> {
    &m.get_index(i - 1).expect("IndexedDataMap::operator() out of range").1 .1
}

/// OCCT IndexedDataMap::Find(k) on the
/// IndexedDataMapOfShapeListOfShape form (asserts when unbound).
pub(crate) fn shape_indexed_data_map_find<'a>(
    m: &'a IndexedDataMapOfShapeListOfShape,
    k: &Shape,
) -> &'a Vec<Shape> {
    &m.get(&crate::brep_algo::tool::shape_key(k))
        .expect("IndexedDataMap::Find unbound")
        .1
}

/// OCCT RemoveCorks (cxx L713-740).
pub(crate) fn remove_corks(s: &mut Shape, faces: &mut OcctIndexedShapeMap) {
    let mut b = BRepBuilder::new();
    let mut brep = BRep::new();
    let ss = b.make_compound(&mut brep, vec![]);
    //-----------------------------------------------------
    // Construction of a shape without caps.
    // and Orientation of caps as in shape S.
    //-----------------------------------------------------
    for cork in bat::explorer(s, ShapeType::Face, ShapeType::Shape) {
        if !faces.contains(&cork) {
            b.add_to_compound(&mut brep, ss.clone(), cork);
        } else {
            // to reset it with proper orientation.
            super::brep_offset_make_offset::indexed_shape_map_remove_key(faces, &cork);
            faces.add(&cork);
        }
    }
    *s = ss;
}

/// OCCT IsConnectedShell (cxx L742-753).
pub(crate) fn is_connected_shell(s: &Shape) -> bool {
    let mut glue = BRepToolsQuilt::new();
    glue.add(s);

    let ss = glue.shells();
    let mut explo = bat::explorer(&ss, ShapeType::Shell, ShapeType::Shape).into_iter();
    explo.next();
    explo.next().is_none()
}

/// OCCT MakeList (cxx L755-777).
pub(crate) fn make_list(
    offset_faces: &mut Vec<Shape>,
    my_init_offset_face: &BRepAlgoImage,
    my_faces: &OcctIndexedShapeMap,
) {
    for root in my_init_offset_face.roots() {
        if !my_faces.contains(root) {
            if my_init_offset_face.has_image(root) {
                for a_it_ls in my_init_offset_face.image(root) {
                    offset_faces.push(a_it_ls);
                }
            }
        }
    }
}

/// OCCT EvalMax (cxx L779-791).
pub(crate) fn eval_max(s: &Shape, tol: &mut f64) {
    for exp in bat::explorer(s, ShapeType::Vertex, ShapeType::Shape) {
        let tolv = brep_tool_tolerance(&exp);
        if tolv > *tol {
            *tol = tolv;
        }
    }
}

/// OCCT BRep_Builder::Add(aCompound, aShape) — the shared-TShape mutation
/// form (BRep_Builder.cxx: the added child is visible through every
/// TopoDS_Shape handle of the compound).
///
/// Architecture difference #61: the kernel BRepBuilder::add_to_compound
/// mutates via Arc::make_mut, which clones the TShape whenever another Arc
/// handle exists — a compound handle retained by the caller then keeps the
/// pre-mutation (empty) TShape, while the BRep pool slot carries the child.
/// This wrapper re-binds the caller's handle onto the pool slot after the
/// add so the OCCT "all handles see the mutation" semantics hold.
pub(crate) fn bb_add_to_compound(brep: &mut BRep, compound: &mut Shape, shape: &Shape) {
    BRepBuilder::new().add_to_compound(brep, compound.clone(), shape.clone());
    compound.data = brep.tshapes[compound.index].clone();
}

// ===========================================================================
// OCCT class (BRepOffset_MakeOffset.hxx L65-285).
// ===========================================================================

/// OCCT BRepOffset_MakeOffset (BRepOffset_MakeOffset.hxx L65-285).
pub struct BRepOffsetMakeOffset {
    /// OCCT: myOffset.
    pub(crate) my_offset: f64,
    /// OCCT: myTol.
    pub(crate) my_tol: f64,
    /// OCCT: myInitialShape.
    pub(crate) my_initial_shape: Shape,
    /// OCCT: myShape.
    pub(crate) my_shape: Shape,
    /// OCCT: myFaceComp (TopoDS_Compound).
    pub(crate) my_face_comp: Shape,
    /// OCCT: myMode (BRepOffset_Mode).
    pub(crate) my_mode: i32,
    /// OCCT: myIsLinearizationAllowed.
    pub(crate) my_is_linearization_allowed: bool,
    /// OCCT: myInter.
    pub(crate) my_inter: bool,
    /// OCCT: mySelfInter.
    pub(crate) my_self_inter: bool,
    /// OCCT: myJoin (GeomAbs_JoinType).
    pub(crate) my_join: GeomAbsJoinType,
    /// OCCT: myThickening.
    pub(crate) my_thickening: bool,
    /// OCCT: myRemoveIntEdges.
    pub(crate) my_remove_int_edges: bool,
    /// OCCT: myFaceOffset.
    pub(crate) my_face_offset: DataMapOfShapeReal,
    /// OCCT: myFaces.
    pub(crate) my_faces: OcctIndexedShapeMap,
    /// OCCT: myOriginalFaces.
    pub(crate) my_original_faces: OcctIndexedShapeMap,
    /// OCCT: myAnalyse (BRepOffset_Analyse).
    pub(crate) my_analyse: BRepOffsetAnalyse,
    /// OCCT: myOffsetShape.
    pub(crate) my_offset_shape: Shape,
    /// OCCT: myInitOffsetFace.
    pub(crate) my_init_offset_face: BRepAlgoImage,
    /// OCCT: myInitOffsetEdge.
    pub(crate) my_init_offset_edge: BRepAlgoImage,
    /// OCCT: myImageOffset.
    pub(crate) my_image_offset: BRepAlgoImage,
    /// OCCT: myImageVV.
    pub(crate) my_image_vv: BRepAlgoImage,
    /// OCCT: myWalls.
    pub(crate) my_walls: Vec<Shape>,
    /// OCCT: myAsDes (Handle(BRepAlgo_AsDes)).
    pub(crate) my_as_des: BRepAlgoAsDes,
    /// OCCT: myEdgeIntEdges.
    pub(crate) my_edge_int_edges: HashMap<crate::feat::loc_ope_wires_on_shape_b::ShapeKey, Vec<Shape>>,
    /// OCCT: myDone.
    pub(crate) my_done: bool,
    /// OCCT: myError (BRepOffset_Error).
    pub(crate) my_error: BRepOffset_Error,
    /// OCCT: myMakeLoops (BRepOffset_MakeLoops).
    pub(crate) my_make_loops: BRepOffsetMakeLoops,
    /// OCCT: myIsPerformSewing — handle bad walls in thicksolid mode.
    pub(crate) my_is_perform_sewing: bool,
    /// OCCT: myIsPlanar.
    pub(crate) my_is_planar: bool,
    /// OCCT: myBadShape.
    pub(crate) my_bad_shape: Shape,
    /// OCCT: myFacePlanfaceMap.
    pub(crate) my_face_planface_map: DataMapOfShapeShape,
    /// OCCT: myGenerated.
    pub(crate) my_generated: Vec<Shape>,
    /// OCCT: myResMap.
    pub(crate) my_res_map: OcctShapeSet,
    /// The rcad class-local BRep pool (architecture difference #41) — the
    /// OCCT global TShape arena stand-in for the shapes built here.
    pub(crate) my_brep: BRep,
}
// ---------------------------------------------------------------------------
// OCCT BRepOffset_Error (BRepOffset_Error.hxx L24-49).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Error (BRepOffset_Error.hxx L20-33).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum BRepOffset_Error {
    NoError,
    UnknownError,
    BadNormalsOnGeometry,
    C0Geometry,
    NullOffset,
    NotConnectedShell,
    CannotTrimEdges,
    CannotFuseVertices,
    CannotExtentEdge,
    UserBreak,
    MixedConnectivity,
}

/// OCCT BRepOffset_Mode (BRepOffset_Mode.hxx L27-33).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum BRepOffset_Mode {
    Skin,
    Pipe,
    RectoVerso,
}

impl Default for BRepOffsetMakeOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetMakeOffset {
    /// OCCT BRepOffset_MakeOffset::BRepOffset_MakeOffset() (cxx L606-609).
    pub fn new() -> Self {
        BRepOffsetMakeOffset {
            my_offset: 0.,
            my_tol: 0.,
            my_initial_shape: Shape::null(),
            my_shape: Shape::null(),
            my_face_comp: Shape::null(),
            my_mode: 0,
            my_is_linearization_allowed: false,
            my_inter: false,
            my_self_inter: false,
            my_join: GeomAbsJoinType::Arc,
            my_thickening: false,
            my_remove_int_edges: false,
            my_face_offset: HashMap::new(),
            my_faces: OcctIndexedShapeMap::new(),
            my_original_faces: OcctIndexedShapeMap::new(),
            my_analyse: BRepOffsetAnalyse::new(),
            my_offset_shape: Shape::null(),
            my_init_offset_face: BRepAlgoImage::new(),
            my_init_offset_edge: BRepAlgoImage::new(),
            my_image_offset: BRepAlgoImage::new(),
            my_image_vv: BRepAlgoImage::new(),
            my_walls: Vec::new(),
            my_as_des: BRepAlgoAsDes::new(),
            my_edge_int_edges: HashMap::new(),
            my_done: false,
            my_error: BRepOffset_Error::NoError,
            my_make_loops: BRepOffsetMakeLoops::new(),
            my_is_perform_sewing: false,
            my_is_planar: false,
            my_bad_shape: Shape::null(),
            my_face_planface_map: HashMap::new(),
            my_generated: Vec::new(),
            my_res_map: HashMap::new(),
            my_brep: BRep::new(),
        }
    }

    /// OCCT BRepOffset_MakeOffset::BRepOffset_MakeOffset(S, Offset, Tol,
    /// Mode, Inter, SelfInter, Join, Thickening, RemoveIntEdges, theRange)
    /// (cxx L613-641).
    #[allow(clippy::too_many_arguments)]
    pub fn with_args(
        s: &Shape,
        offset: f64,
        tol: f64,
        mode: BRepOffset_Mode,
        inter: bool,
        self_inter: bool,
        join: GeomAbsJoinType,
        thickening: bool,
        remove_int_edges: bool,
    ) -> Self {
        let mut res = Self::new();
        res.my_offset = offset;
        res.my_tol = tol;
        res.my_initial_shape = s.clone();
        res.my_shape = s.clone();
        res.my_mode = mode as i32;
        res.my_inter = inter;
        res.my_self_inter = self_inter;
        res.my_join = join;
        res.my_thickening = thickening;
        res.my_remove_int_edges = remove_int_edges;
        res.my_done = false;

        res.my_is_linearization_allowed = true;

        res.make_offset_shape();
        res
    }

    /// OCCT BRepOffset_MakeOffset::Initialize (cxx L643-670).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize(
        &mut self,
        s: &Shape,
        offset: f64,
        tol: f64,
        mode: BRepOffset_Mode,
        inter: bool,
        self_inter: bool,
        join: GeomAbsJoinType,
        thickening: bool,
        remove_int_edges: bool,
    ) {
        self.my_offset = offset;
        self.my_initial_shape = s.clone();
        self.my_shape = s.clone();
        self.my_tol = tol;
        self.my_mode = mode as i32;
        self.my_inter = inter;
        self.my_self_inter = self_inter;
        self.my_join = join;
        self.my_thickening = thickening;
        self.my_remove_int_edges = remove_int_edges;
        self.my_is_linearization_allowed = true;
        self.my_done = false;
        self.my_is_perform_sewing = false;
        self.my_is_planar = false;
        self.clear();
    }

    /// OCCT BRepOffset_MakeOffset::Clear (cxx L672-689).
    pub fn clear(&mut self) {
        self.my_offset_shape = Shape::null();
        self.my_init_offset_face.clear();
        self.my_init_offset_edge.clear();
        self.my_image_offset.clear();
        self.my_image_vv.clear();
        self.my_faces.clear();
        self.my_original_faces.clear();
        self.my_face_offset.clear();
        self.my_edge_int_edges.clear();
        self.my_as_des.clear();
        self.my_done = false;
        self.my_generated.clear();
        self.my_res_map.clear();
    }

    /// OCCT BRepOffset_MakeOffset::AllowLinearization (cxx L691-694).
    pub fn allow_linearization(&mut self, the_is_allowed: bool) {
        self.my_is_linearization_allowed = the_is_allowed;
    }

    /// OCCT BRepOffset_MakeOffset::AddFace (cxx L698-702).
    pub fn add_face(&mut self, f: &Shape) {
        self.my_original_faces.add(f);
    }

    /// OCCT BRepOffset_MakeOffset::SetOffsetOnFace (cxx L706-709).
    pub fn set_offset_on_face(&mut self, f: &Shape, off: f64) {
        shape_data_map::bind(&mut self.my_face_offset, f, off);
    }

    /// OCCT BRepOffset_MakeOffset::SetFaces (cxx L795-816).
    pub(crate) fn set_faces(&mut self) {
        for ii in 1..=self.my_original_faces.extent() {
            let mut a_face = self.my_original_faces.at_1(ii).clone();
            if let Some(a_planface) = shape_data_map::seek(&self.my_face_planface_map, &a_face) {
                a_face = a_planface.clone();
            }

            self.my_faces.add(&a_face);
            //-------------
            // MAJ SD.
            //-------------
            self.my_init_offset_face.set_root(&a_face);
            self.my_init_offset_face.bind(&a_face, &a_face);
            self.my_image_offset.set_root(&a_face);
        }
    }

    /// OCCT BRepOffset_MakeOffset::SetFacesWithOffset (cxx L818-835).
    pub(crate) fn set_faces_with_offset(&mut self) {
        let copied_keys: Vec<Shape> = self
            .my_face_planface_map
            .values()
            .map(|(k, _)| k.clone())
            .collect();
        for a_face in &copied_keys {
            let a_planface = shape_data_map::value(&self.my_face_planface_map, a_face).clone();
            if shape_data_map::is_bound(&self.my_face_offset, a_face) {
                let an_offset = shape_data_map::find(&self.my_face_offset, a_face);
                shape_data_map::un_bind(&mut self.my_face_offset, a_face);
                shape_data_map::bind(&mut self.my_face_offset, &a_planface, an_offset);
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::IsDone (cxx L1181-1184).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepOffset_MakeOffset::Error (cxx L1188-1191).
    pub fn error(&self) -> BRepOffset_Error {
        self.my_error
    }

    /// OCCT BRepOffset_MakeOffset::Shape (cxx L1194-1197).
    pub fn shape(&self) -> &Shape {
        &self.my_offset_shape
    }

    /// OCCT BRepOffset_MakeOffset::InitShape (hxx L216).
    pub fn init_shape(&self) -> &Shape {
        &self.my_initial_shape
    }

    /// OCCT BRepOffset_MakeOffset::GetAnalyse (hxx L164).
    pub fn get_analyse(&self) -> &BRepOffsetAnalyse {
        &self.my_analyse
    }

    /// OCCT BRepOffset_MakeOffset::GetBadShape (cxx L4497-4500).
    pub fn get_bad_shape(&self) -> &Shape {
        &self.my_bad_shape
    }

    /// OCCT BRepOffset_MakeOffset::OffsetFacesFromShapes (cxx L3928-3931).
    pub fn offset_faces_from_shapes(&self) -> &BRepAlgoImage {
        &self.my_init_offset_face
    }

    /// OCCT BRepOffset_MakeOffset::GetJoinType (cxx L3938-3941).
    pub fn get_join_type(&self) -> GeomAbsJoinType {
        self.my_join
    }

    /// OCCT BRepOffset_MakeOffset::OffsetEdgesFromShapes (cxx
    /// L3944-3947).
    pub fn offset_edges_from_shapes(&self) -> &BRepAlgoImage {
        &self.my_init_offset_edge
    }

    /// OCCT BRepOffset_MakeOffset::ClosingFaces (cxx L3954-3958).
    pub fn closing_faces(&self) -> &OcctIndexedShapeMap {
        &self.my_original_faces
    }

    /// OCCT BRepOffset_MakeOffset::BuildFaceComp (cxx L1873-1891).
    /// OCCT BRepOffset_MakeOffset::BuildFaceComp (cxx L1869-1889).
    pub(crate) fn build_face_comp(&mut self) {
        let mut a_bb = BRepBuilder::new();
        self.my_face_comp = a_bb.make_compound(&mut self.my_brep, vec![]);
        let faces = bat::explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape);
        for explo in faces {
            let mut a_face = explo;
            let an_or = a_face.orientation;
            if let Some(a_planface) = shape_data_map::seek(&self.my_face_planface_map, &a_face) {
                a_face = a_planface.clone();
            }
            let a_face = oriented(&a_face, an_or);
            // OCCT L1881: aBB.Add(myFaceComp, aFace.Oriented(anOr)) — the
            // retained-handle form (architecture difference #61).
            bb_add_to_compound(&mut self.my_brep, &mut self.my_face_comp, &a_face);
        }
    }

    /// OCCT BRepOffset_MakeOffset::SelfInter (cxx L2154-2169).
    pub(crate) fn self_inter(&mut self, _modif: &mut OcctShapeSet) {
        // throw Standard_NotImplemented();
        panic!("Standard_NotImplemented: BRepOffset_MakeOffset::SelfInter");
    }
}

/// OCCT IndexedMap::RemoveKey — the rcad rebuild form (the
/// OcctIndexedShapeMap of brep_offset_tool.rs carries no removal; the
/// rebuild keeps the OCCT insertion order of the remaining entries).
pub(crate) fn indexed_shape_map_remove_key(m: &mut OcctIndexedShapeMap, k: &Shape) {
    let kept: Vec<Shape> = m.iter().filter(|s| !s.is_same(k)).cloned().collect();
    m.clear();
    for s in kept {
        m.add(&s);
    }
}

/// OCCT BRep_Tool::IsClosed(S) — the bat re-host alias (the MakeOffset
/// call sites keep the OCCT BRep_Tool::IsClosed form).
pub(crate) fn brep_tool_is_closed(the_s: &Shape) -> bool {
    bat::shape_is_closed(the_s)
}

/// OCCT BOPAlgo_MakerVolume (TKBO/BOPAlgo) constructor — the GAP carrier
/// default form.
impl BOPAlgoMakerVolume {
    pub fn new() -> Self {
        BOPAlgoMakerVolume
    }
}

impl Default for BOPAlgoMakerVolume {
    fn default() -> Self {
        Self::new()
    }
}
