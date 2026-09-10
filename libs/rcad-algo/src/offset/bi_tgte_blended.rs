// OCCT BiTgte_Blend.cxx L1-2664 + BiTgte_Blend.hxx L32-150 — 1:1
// translation (the class part split into bi_tgte_blended_b.rs).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BiTgte/
//         BiTgte_Blend.cxx / .hxx
//
// OCCT inheritance chain (hxx L32): none — BiTgte_Blend is a standalone
// class.
//
// Architecture differences (continuing the numbering of
// brep_offset_offset.rs):
// 21. NCollection_DataMap / NCollection_IndexedDataMap / NCollection_Map /
//     NCollection_IndexedMap of TopoDS_Shape -> HashMap<Shape, _> /
//     IndexMap<Shape, _> / HashSet<Shape> / IndexMap<Shape, _> (the rcad
//     Shape carries the TopTools_ShapeMapHasher identity). myAncestors
//     (TopExp::MapShapesAndAncestors) -> the shared re-host
//     feat::loc_ope_glued_shape::map_shapes_and_ancestors keyed by
//     (TShape ptr, Location). myIndices (HArray1<int> 1-based) -> Vec<i32>
//     0-based with an explicit -1 offset at the access sites.
// 22. BRepOffset_Analyse (TKOffset/BRepOffset) — GAP carrier (the Stage 2b
//     unit owns the translation; the leaves panic until it lands).
// 23. BRepOffset_Inter3d / BRepOffset_Inter2d / BRepOffset_MakeLoops /
//     BRepOffset_Tool::EnLargeFace / BRepOffset_Interval — GAP carriers
//     (same Stage 2b unit).
// 24. BRepBuilderAPI_Sewing (TKTopAlgo/BRepBuilderAPI) — GAP carrier.
// 25. BRepLib::SameParameter — GAP static leaf.  (BRepLib::BuildCurves3d
//     is translated — topalgo/brep_lib/build_curves3d.rs; the call sites
//     call the real body over my_brep.)
// 26. Approx_FitAndDivide + AppCont_Function + AppParCurves_MultiCurve +
//     Convert_CompBezierCurvesToBSplineCurve (TKGeomBase) — GAP carriers;
//     BSplCLib::Reparametrize is translated (bspl_lib::reparametrize).
// 27. GeomAPI_ProjectPointOnCurve / Geom2dAPI_ProjectPointOnCurve
//     (TKTopAlgo/GeomAPI) — GAP carriers (the Init/NearestPoint form is the
//     #12 carrier of brep_offset_offset.rs).
// 28. GeomAPI::To3d + ElSLib::*VIso/*UIso + gp_Circ::Rotate + gp_Lin::
//     Translate — GAP leaves inside KPartCurve3d (the TKMath/ElSLib iso
//     constructors are not translated; cf. #20 of brep_offset_offset.rs).
// 29. BRepLib_MakeEdge -> BRepBuilder::add_edge over the class-local BRep
//     pool (arch. diff. #4/#19); the (PC, S, V1, V2, f, l) pcurve form has
//     no face-less pcurve storage in rcad — GAP leaf.
// 30. GeomAdaptor_Surface / Geom2dAdaptor_Curve / GeomAdaptor_Curve
//     GetType() -> the rcad Surface3 / Curve3 variant matches (arch. diff.
//     #18 of brep_offset_offset.rs); Geom_TrimmedCurve flattening (the
//     GeomAdaptor load of the basis) is an explicit one-level strip.
// 31. Message_ProgressRange / OSD_Chronometer (OCCT_DEBUG) — omitted (rcad
//     carries no progress scope; the debug chronometers are debug-only).

use std::collections::{HashMap, HashSet};

use glam::{DVec2, DVec3};
use indexmap::IndexMap;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{Circle3, Curve2d, Curve3, Line3, Surface3, TrimmedCurve3, BSplineCurve3};
use rcad_kernel::math::bspl_lib;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::feat::loc_ope_wires_on_shape::brep_tool_pnt;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_tolerance, ShapeKey,
};

use super::brep_offset_make_simple_offset::{edge_curve_of, top_exp_vertices_shape};
use super::brep_offset_offset_b::BRepOffsetOffset;
use super::bi_tgte_contact::BiTgteContactType;
use super::bi_tgte_curve_on_edge::BiTgteCurveOnEdge;

use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity;
// OCCT BRepTools_Quilt (TKBRep/BRepTools/BRepTools_Quilt.hxx / .cxx) — the
// single 1:1 body lives in crate::topalgo::brep_tools_quilt (its correct
// location); the local GAP carrier of NbBranches was retired.
use crate::topalgo::brep_tools_quilt::BRepToolsQuilt;

use rcad_kernel::core::precision::{
    ANGULAR as PRECISION_ANGULAR, APPROXIMATION as PRECISION_APPROXIMATION,
    CONFUSION as PRECISION_CONFUSION,
};

type AncestorsMap = IndexMap<ShapeKey, (Shape, Vec<Shape>)>;

// ===========================================================================
// GAP carriers (architecture differences #22-#27, #29).
// ===========================================================================

// OCCT BRepOffset_Interval / BRepOffset_Analyse — the real bodies live in
// super::brep_offset_analyse (the E0 carrier-switch list; the local panic
// carriers are deleted).  Re-exported for the Blended consumers.
pub use super::brep_offset_analyse::{BRepOffsetAnalyse, BRepOffsetInterval};

/// OCCT BRepOffset_Inter3d (TKOffset/BRepOffset/BRepOffset_Inter3d.hxx /
/// .cxx) — GAP carrier (arch. diff. #23).
pub struct BRepOffsetInter3d;

impl BRepOffsetInter3d {
    /// OCCT BRepOffset_Inter3d::BRepOffset_Inter3d(AsDes, Side, Tol).
    pub fn new(_as_des: &BRepAlgoAsDes, _side: State, _tol: f64) -> Self {
        BRepOffsetInter3d
    }

    /// OCCT BRepOffset_Inter3d::IsDone(F1, F2).
    pub fn is_done(&self, _f1: &Shape, _f2: &Shape) -> bool {
        panic!("GAP: BRepOffset_Inter3d::IsDone (TKOffset/BRepOffset not translated)");
    }

    /// OCCT BRepOffset_Inter3d::FaceInter(F1, F2, InitOffset).
    pub fn face_inter(
        &mut self,
        _f1: &Shape,
        _f2: &Shape,
        _init_offset: &BRepAlgoImage,
    ) {
        panic!("GAP: BRepOffset_Inter3d::FaceInter (TKOffset/BRepOffset not translated)");
    }

    /// OCCT BRepOffset_Inter3d::NewEdges() -> const TopTools_MapOfShape&.
    pub fn new_edges(&self) -> Vec<Shape> {
        panic!("GAP: BRepOffset_Inter3d::NewEdges (TKOffset/BRepOffset not translated)");
    }
}

/// OCCT BRepOffset_Inter2d (TKOffset/BRepOffset/BRepOffset_Inter2d.hxx /
/// .cxx) — GAP statics (arch. diff. #23).
pub struct BRepOffsetInter2d;

impl BRepOffsetInter2d {
    /// OCCT BRepOffset_Inter2d::Compute(AsDes, F, NewEdges, Tol, EdgeInt,
    /// DMVV, theProgress).
    #[allow(clippy::too_many_arguments)]
    pub fn compute(
        _as_des: &BRepAlgoAsDes,
        _f: &Shape,
        _new_edges: &IndexMap<Shape, ()>,
        _tol: f64,
        _edge_int: &mut HashMap<Shape, Vec<Shape>>,
        _dmvv: &mut AncestorsMap,
    ) {
        panic!("GAP: BRepOffset_Inter2d::Compute (TKOffset/BRepOffset not translated)");
    }

    /// OCCT BRepOffset_Inter2d::FuseVertices(DMVV, AsDes, Image).
    pub fn fuse_vertices(
        _dmvv: &AncestorsMap,
        _as_des: &BRepAlgoAsDes,
        _image: &mut BRepAlgoImage,
    ) {
        panic!("GAP: BRepOffset_Inter2d::FuseVertices (TKOffset/BRepOffset not translated)");
    }
}

/// OCCT BRepOffset_MakeLoops (TKOffset/BRepOffset/BRepOffset_MakeLoops.hxx /
/// .cxx) — GAP carrier (arch. diff. #23).
#[derive(Default)]
pub struct BRepOffsetMakeLoops;

impl BRepOffsetMakeLoops {
    /// OCCT BRepOffset_MakeLoops::BRepOffset_MakeLoops().
    pub fn new() -> Self {
        BRepOffsetMakeLoops
    }

    /// OCCT BRepOffset_MakeLoops::Build(LOF, AsDes, ImageOffset, Image,
    /// theProgress).
    pub fn build(
        &mut self,
        _lof: &mut Vec<Shape>,
        _as_des: &BRepAlgoAsDes,
        _image_offset: &BRepAlgoImage,
        _image: &mut BRepAlgoImage,
    ) {
        panic!("GAP: BRepOffset_MakeLoops::Build (TKOffset/BRepOffset not translated)");
    }
}

/// OCCT BRepBuilderAPI_Sewing (TKTopAlgo/BRepBuilderAPI/
/// BRepBuilderAPI_Sewing.hxx / .cxx) — GAP carrier (arch. diff. #24).
#[derive(Default)]
pub struct BRepBuilderAPISewing;

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::BRepBuilderAPI_Sewing(Tolerance).
    pub fn new(_tolerance: f64) -> Self {
        BRepBuilderAPISewing
    }

    /// OCCT BRepBuilderAPI_Sewing::Add(shape).
    pub fn add(&mut self, _shape: &Shape) {
        panic!("GAP: BRepBuilderAPI_Sewing::Add (TKTopAlgo/BRepBuilderAPI not translated)");
    }

    /// OCCT BRepBuilderAPI_Sewing::Perform(theProgress).
    pub fn perform(&mut self) {
        panic!("GAP: BRepBuilderAPI_Sewing::Perform (TKTopAlgo/BRepBuilderAPI not translated)");
    }

    /// OCCT BRepBuilderAPI_Sewing::SewedShape().
    pub fn sewed_shape(&self) -> Shape {
        panic!("GAP: BRepBuilderAPI_Sewing::SewedShape (TKTopAlgo/BRepBuilderAPI not translated)");
    }

    /// OCCT BRepBuilderAPI_Sewing::IsModified(shape).
    pub fn is_modified(&self, _shape: &Shape) -> bool {
        panic!("GAP: BRepBuilderAPI_Sewing::IsModified (TKTopAlgo/BRepBuilderAPI not translated)");
    }

    /// OCCT BRepBuilderAPI_Sewing::Modified(shape).
    pub fn modified(&self, _shape: &Shape) -> Shape {
        panic!("GAP: BRepBuilderAPI_Sewing::Modified (TKTopAlgo/BRepBuilderAPI not translated)");
    }
}

/// OCCT BRepLib::SameParameter(E, Tol) — GAP static leaf (arch. diff. #25).
pub fn brep_lib_same_parameter(_the_e: &Shape, _the_tol: f64) {
    panic!("GAP: BRepLib::SameParameter (TKTopAlgo/BRepLib not translated)");
}

/// OCCT BRepTools::Update(S) — GAP no-op re-host (the same re-host as the
/// brep_offset_offset.rs one; the face storage it refreshes is carried by
/// the rcad TFaceData fields directly).
fn brep_tools_update(_the_s: &Shape) {}

/// OCCT BRepLib_MakeEdge(Curve, V1, V2) — the (done, edge) pair of the
/// builder (arch. diff. #29).
fn brep_lib_make_edge_3d(
    the_brep: &mut BRep,
    the_curve: &Curve3,
    the_v1: &Shape,
    the_v2: &Shape,
) -> (bool, Shape) {
    let mut b = BRepBuilder::new();
    let first = curve_first_parameter(the_curve);
    let last = curve_last_parameter(the_curve);
    let e = b.add_edge(the_brep, Some(the_curve.clone()), the_v1.clone(), the_v2.clone(), [
        first, last,
    ]);
    (true, e)
}

/// OCCT BRepLib_MakeEdge(PCurve, Surface, V1, V2, f, l) — the face-less
/// pcurve-on-surface edge form has no rcad storage (the pcurves are keyed
/// by the face TShape) — GAP leaf (arch. diff. #29).
fn brep_lib_make_edge_pcurve(
    _the_brep: &mut BRep,
    _the_pcurve: &Curve2d,
    _the_surface: &Surface3,
    _the_v1: &Shape,
    _the_v2: &Shape,
    _the_first: f64,
    _the_last: f64,
) -> (bool, Shape) {
    panic!("GAP: BRepLib_MakeEdge(PCurve, Surface, V1, V2, f, l) (arch. diff. #29)");
}

/// OCCT Geom2dAPI_ProjectPointOnCurve (TKTopAlgo/GeomAPI) — GAP carrier
/// (arch. diff. #27); the Init(P, PC, f, l) form of IsOnRestriction.
#[derive(Default)]
pub struct Geom2dApiProjectPointOnCurve;

impl Geom2dApiProjectPointOnCurve {
    /// OCCT Geom2dAPI_ProjectPointOnCurve::Init(P, Curve, Uf, Ul).
    pub fn init(&mut self, _the_p: DVec2, _the_curve: &Curve2d, _uf: f64, _ul: f64) {}

    /// OCCT Geom2dAPI_ProjectPointOnCurve::NbPoints().
    pub fn nb_points(&self) -> i32 {
        panic!("GAP: Geom2dAPI_ProjectPointOnCurve::NbPoints (TKTopAlgo/GeomAPI not translated)");
    }

    /// OCCT Geom2dAPI_ProjectPointOnCurve::LowerDistance().
    pub fn lower_distance(&self) -> f64 {
        panic!("GAP: Geom2dAPI_ProjectPointOnCurve::LowerDistance (TKTopAlgo/GeomAPI not translated)");
    }
}

// ---------------------------------------------------------------------------
// OCCT statics (BiTgte_Blend.cxx L88-703).
// ---------------------------------------------------------------------------

// OCCT BiTgte_Blend.cxx L98-136 — IsOnRestriction.
fn is_on_restriction(v: &Shape, cur_e: &Shape, f: &Shape, e: &mut Shape) -> bool {
    // find if Vertex V of CurE is on a restriction of F.
    // if yes, store this restriction in E.

    // dub - 03 01 97
    // Method somewhat brutal : possible to really optimize by a
    // direct call the SD of intersections -> See LBR

    // OCCT L111: CurC = BRep_Tool::CurveOnSurface(CurE, F, f, l).
    let Some((cur_c, _f_par, _l_par)) = brep_tool_curve_on_surface(cur_e, f) else {
        return false; // OCCT: the null pcurve handle would raise downstream.
    };
    // OCCT L112: U = BRep_Tool::Parameter(V, CurE, F) — the rcad re-host
    // carries the (V, E) form (the face-keyed parameter is the same edge
    // parameter here).
    let u = crate::feat::loc_ope_wires_on_shape_b::brep_tool_parameter(v, cur_e);
    // OCCT L113: P = CurC->Value(U).
    use rcad_kernel::geom::Curve2dEval;
    let p = cur_c.point_at(u);

    // OCCT L115: Geom2dAPI_ProjectPointOnCurve Proj.
    let mut proj = Geom2dApiProjectPointOnCurve::default();

    // The tolerance is exaggerated : it is better to construct too many
    // tubes than to miss intersections.
    // double Tol = 100 * BRep_Tool::Tolerance(V);
    // OCCT L120: Tol = BRep_Tool::Tolerance(V).
    let tol = brep_tool_tolerance(v);
    // OCCT L121-134.
    for a_local_edge in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        *e = a_local_edge;
        // OCCT L125: PC = BRep_Tool::CurveOnSurface(E, F, f, l).
        let Some((pc, pf, pl)) = brep_tool_curve_on_surface(e, f) else {
            continue;
        };
        proj.init(p, &pc, pf, pl);
        if proj.nb_points() > 0 {
            if proj.lower_distance() < tol {
                return true;
            }
        }
    }
    false
}

// OCCT BiTgte_Blend.cxx L140-195 — Add.
fn add(
    e: &Shape,
    map: &mut IndexMap<Shape, ()>,
    s: &Shape,
    of: &BRepOffsetOffset,
    analyse: &BRepOffsetAnalyse,
    warning_sur_bord_libre: bool,
) {
    // If WarningSurBordLibre = TRUE, no propagation if the edge is open.
    let type_ = s.shape_type();

    if type_ == ShapeType::Face {
        // OCCT L152-173.
        for ori_e in explorer(s, ShapeType::Edge, ShapeType::Shape) {
            // OCCT L156: aLocalShape = OF.Generated(OriE); IE =
            // TopoDS::Edge(aLocalShape).
            let ie = of.generated(&ori_e);
            if e.is_equal(&ie) {
                if warning_sur_bord_libre {
                    // It is checked that the border is not free.
                    // OCCT L164: L = Analyse.Ancestors(OriE).
                    let l = analyse.ancestors(&ori_e);
                    if l.len() == 1 {
                        break; // Nothing is done.
                    }
                }
                map.insert(ori_e, ());
                break;
            }
        }
    } else if type_ == ShapeType::Edge {
        // OCCT L175-194.
        for a_local_vertex in explorer(s, ShapeType::Vertex, ShapeType::Shape) {
            // OCCT L180: aLocalShape = OF.Generated(exp.Current()).
            let ie = of.generated(&a_local_vertex);
            if e.is_equal(&ie) {
                // OCCT L185: L = Analyse.Ancestors(exp.Current()).
                let l = analyse.ancestors(&a_local_vertex).clone();
                for it_value in l {
                    map.insert(it_value, ());
                }
                break;
            }
        }
    }
}

// OCCT BiTgte_Blend.cxx L199-210 — IsInFace.
fn is_in_face(e: &Shape, f: &Shape) -> bool {
    for a_local_edge in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        if e.is_same(&a_local_edge) {
            return true;
        }
    }
    false
}

// OCCT BiTgte_Blend.cxx L214-378 — KPartCurve3d.
fn k_part_curve_3d(edge: &Shape, curve: &Curve2d, surf: &Surface3, the_brep: &mut BRep) {
    // try to find the particular case
    // if not found call BRepLib::BuildCurve3d

    // OCCT L221: TopLoc_Location Loc — the rcad location index (0 =
    // identity).
    let loc: u32 = 0;
    let tol = PRECISION_CONFUSION;

    // Search only isos on analytical surfaces.
    // OCCT L225-228: Geom2dAdaptor_Curve C(Curve); GeomAdaptor_Surface
    // S(Surf); CTy/STy = GetType() — the rcad variant matches (arch.
    // diff. #30).
    let mut the_builder = BRepBuilder::new();

    if !matches!(surf, Surface3::Plane(_)) {
        // if plane buildcurve3d manage KPart
        if let Curve2d::Line(c_line) = curve {
            // OCCT L233: CTy == GeomAbs_Line.
            let d = c_line.direction;
            // OCCT L236: D.IsParallel(gp::DX2d(), Precision::Angular()) —
            // the cross-magnitude form of the brep_offset_offset.rs
            // compute_curve3d precedent.
            let is_parallel_dx = (d.x * 0.0 - d.y * 1.0).abs() < PRECISION_ANGULAR; // |d x DX2d|
            let is_parallel_dy = (d.x * 1.0 - d.y * 0.0).abs() < PRECISION_ANGULAR; // |d x DY2d|
            let p = c_line.origin;
            if is_parallel_dx {
                // Iso V.
                if let Surface3::Sphere(sph) = surf {
                    // OCCT L238-260.
                    if ((p.y.abs() - std::f64::consts::FRAC_PI_2).abs())
                        < rcad_kernel::core::precision::PCONFUSION
                    {
                        // OCCT L243: TheBuilder.Degenerated(Edge, true).
                        the_builder.set_edge_degenerated(the_brep, edge.clone(), true);
                    } else {
                        // OCCT L247-258: SphereVIso + Rotate + Reverse +
                        // UpdateEdge — the ElSLib iso/GP rotation leaves.
                        let circle = elslib_sphere_v_iso(sph, p.y);
                        let circle = circle_rotate(circle, sph.axis, p.x);
                        let mut circle = circle;
                        if d.x < 0.0 {
                            circle = circle_reversed(circle);
                        }
                        update_edge_curve3d_gap(edge, &Curve3::Circle(circle), loc, tol);
                    }
                } else if let Surface3::Cylinder(cyl) = surf {
                    // OCCT L261-276.
                    let circle = elslib_cylinder_v_iso(cyl, cyl.radius, p.y);
                    let circle = circle_rotate(circle, cyl.axis, p.x);
                    let mut circle = circle;
                    if d.x < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_curve3d_gap(edge, &Curve3::Circle(circle), loc, tol);
                } else if let Surface3::Cone(cone) = surf {
                    // OCCT L277-292.
                    let circle = elslib_cone_v_iso(cone, cone.radius, cone.half_angle_rad, p.y);
                    let circle = circle_rotate(circle, cone.axis, p.x);
                    let mut circle = circle;
                    if d.x < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_curve3d_gap(edge, &Curve3::Circle(circle), loc, tol);
                } else if let Surface3::Torus(tore) = surf {
                    // OCCT L293-308.
                    let circle = elslib_torus_v_iso(tore, tore.major_radius, tore.minor_radius, p.y);
                    let circle = circle_rotate(circle, tore.axis, p.x);
                    let mut circle = circle;
                    if d.x < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_curve3d_gap(edge, &Curve3::Circle(circle), loc, tol);
                }
            } else if is_parallel_dy {
                // Iso U.
                if let Surface3::Sphere(sph) = surf {
                    // OCCT L312-335: calculate iso 0; set to sameparameter
                    // (rotation of the circle - offset from Y);
                    // transformation by iso U (= P.X()).
                    let circle = elslib_sphere_u_iso(sph, sph.radius, 0.0);
                    let circle = circle_rotate(circle, sph.axis, p.y);
                    let circle = circle_rotate(circle, sph.axis, p.x);
                    let mut circle = circle;
                    if d.y < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_curve3d_gap(edge, &Curve3::Circle(circle), loc, tol);
                } else if let Surface3::Cylinder(cyl) = surf {
                    // OCCT L337-350.
                    let line = elslib_cylinder_u_iso(cyl, cyl.radius, p.x);
                    let tr = line.direction * p.y;
                    let line = line3_translate(line, tr);
                    let mut line = line;
                    if d.y < 0.0 {
                        line.direction = -line.direction;
                    }
                    update_edge_curve3d_gap(edge, &Curve3::Line(line), loc, tol);
                } else if let Surface3::Cone(cone) = surf {
                    // OCCT L352-366.
                    let line = elslib_cone_u_iso(cone, cone.radius, cone.half_angle_rad, p.x);
                    let tr = line.direction * p.y;
                    let line = line3_translate(line, tr);
                    let mut line = line;
                    if d.y < 0.0 {
                        line.direction = -line.direction;
                    }
                    update_edge_curve3d_gap(edge, &Curve3::Line(line), loc, tol);
                } else if let Surface3::Torus(_tore) = surf {
                    // OCCT L367-369: the empty Torus branch.
                }
            }
        }
    } else {
        // Case Plane
        // OCCT L375: C3d = GeomAPI::To3d(Curve, S.Plane()) — GAP leaf
        // (arch. diff. #28).
        if let Surface3::Plane(pl) = surf {
            let c3d = geom_api_to_3d(curve, pl);
            // OCCT L376: TheBuilder.UpdateEdge(Edge, C3d, Loc, Tol).
            update_edge_curve3d_gap(edge, &c3d, loc, tol);
        }
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C3d, L, Tol) — GAP no-op re-host (the
/// rcad 3d-curve representation is index-based; the same re-host as the
/// brep_offset_offset.rs update_edge_curve3d_gap).
fn update_edge_curve3d_gap(_the_e: &Shape, _the_c: &Curve3, _the_loc: u32, _the_tol: f64) {}

/// OCCT GeomAPI::To3d(Curve, Plane) — GAP leaf (arch. diff. #28).
fn geom_api_to_3d(_the_curve: &Curve2d, _the_plane: &rcad_kernel::geom::Plane) -> Curve3 {
    panic!("GAP: GeomAPI::To3d (TKTopAlgo/GeomAPI not translated)");
}

// --- KPartCurve3d GAP leaves (architecture difference #28; the same leaf
// --- family as the #20 stand-ins of brep_offset_offset.rs) ---

type GpCirc = Circle3;
type GpLin = Line3;

/// OCCT ElSLib::SphereVIso(Axis, Radius, V) — GAP.
fn elslib_sphere_v_iso(_sph: &rcad_kernel::geom::SphericalSurface, _v: f64) -> GpCirc {
    panic!("GAP: ElSLib::SphereVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::CylinderVIso(Axis, Radius, V) — GAP.
fn elslib_cylinder_v_iso(
    _cyl: &rcad_kernel::geom::CylindricalSurface,
    _radius: f64,
    _v: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::CylinderVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::ConeVIso(Axis, RefRadius, SemiAngle, V) — GAP.
fn elslib_cone_v_iso(
    _cone: &rcad_kernel::geom::ConicalSurface,
    _ref_radius: f64,
    _semi_angle: f64,
    _v: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::ConeVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::TorusVIso(Axis, MajorR, MinorR, V) — GAP.
fn elslib_torus_v_iso(
    _tore: &rcad_kernel::geom::ToroidalSurface,
    _major: f64,
    _minor: f64,
    _v: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::TorusVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::SphereUIso(Axis, Radius, U) — GAP.
fn elslib_sphere_u_iso(
    _sph: &rcad_kernel::geom::SphericalSurface,
    _radius: f64,
    _u: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::SphereUIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::CylinderUIso(Position, Radius, U) — GAP.
fn elslib_cylinder_u_iso(
    _cyl: &rcad_kernel::geom::CylindricalSurface,
    _radius: f64,
    _u: f64,
) -> GpLin {
    panic!("GAP: ElSLib::CylinderUIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::ConeUIso(Position, RefRadius, SemiAngle, U) — GAP.
fn elslib_cone_u_iso(
    _cone: &rcad_kernel::geom::ConicalSurface,
    _ref_radius: f64,
    _semi_angle: f64,
    _u: f64,
) -> GpLin {
    panic!("GAP: ElSLib::ConeUIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::TorusUIso(Axis, MajorR, MinorR, U) — GAP.
fn elslib_torus_u_iso(
    _tore: &rcad_kernel::geom::ToroidalSurface,
    _major: f64,
    _minor: f64,
    _u: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::TorusUIso (TKMath/ElSLib not translated)");
}

/// OCCT gp_Circ::Rotate(gp_Ax1, Angle) — GAP.
fn circle_rotate(_c: GpCirc, _axis: DVec3, _angle: f64) -> GpCirc {
    panic!("GAP: gp_Circ::Rotate (gp_Trsf rotation not translated)");
}

/// OCCT Geom_Circle::Reverse() — the parameter reversal (the rcad Circle3
/// carries no sense flag; the x/y frame swap is the direction reversal).
fn circle_reversed(mut c: GpCirc) -> GpCirc {
    std::mem::swap(&mut c.x_dir, &mut c.y_dir);
    c
}

/// OCCT gp_Lin::Translate(gp_Vec) — GAP.
fn line3_translate(_l: GpLin, _tr: DVec3) -> GpLin {
    panic!("GAP: gp_Lin::Translate (gp_Trsf translation not translated)");
}

/// OCCT Geom_Curve::FirstParameter() — the rcad curve-bound read (the
/// trimmed range / the BSpline domain / the conic period).
fn curve_first_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.first,
        Curve3::BSpline(b) => b.first_parameter(),
        _ => 0.0,
    }
}

/// OCCT Geom_Curve::LastParameter() — the rcad curve-bound read.
fn curve_last_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.last,
        Curve3::BSpline(b) => b.last_parameter(),
        _ => std::f64::consts::TAU,
    }
}

// ---------------------------------------------------------------------------
// OCCT MakeCurve_Function (BiTgte_Blend.cxx L382-412) — the AppCont_Function
// subclass; the base class is carried by the plain struct (the GAP carrier
// of Approx_FitAndDivide consumes it).
// ---------------------------------------------------------------------------

/// OCCT MakeCurve_Function (BiTgte_Blend.cxx L382-412).
struct MakeCurveFunction {
    my_curve: BiTgteCurveOnEdge, // OCCT: myCurve
}

impl MakeCurveFunction {
    /// OCCT MakeCurve_Function::MakeCurve_Function(C).
    fn new(c: &BiTgteCurveOnEdge) -> Self {
        // OCCT L390-391: myNbPnt = 1; myNbPnt2d = 0 — carried by the GAP
        // base (the approximation machinery is not translated).
        MakeCurveFunction { my_curve: c.clone() }
    }

    /// OCCT L394: FirstParameter().
    fn first_parameter(&self) -> f64 {
        self.my_curve.first_parameter()
    }

    /// OCCT L396: LastParameter().
    fn last_parameter(&self) -> f64 {
        self.my_curve.last_parameter()
    }

    /// OCCT L398-404: Value(theT, thePnt2d, thePnt).
    fn value(&self, the_t: f64, the_pnt: &mut [DVec3]) -> bool {
        the_pnt[0] = self.my_curve.eval_d0(the_t);
        true
    }

    /// OCCT L406-411: D1(theT, theVec2d, theVec).
    fn d1(&self, _the_t: f64, _the_vec: &mut [DVec3]) -> bool {
        false
    }
}

// OCCT Approx_FitAndDivide (TKGeomBase/Approx) — GAP carrier (arch. diff.
// #26).
struct ApproxFitAndDivide;

impl ApproxFitAndDivide {
    /// OCCT Approx_FitAndDivide(Func, DegMin, DegMax, Tol3D, Tol2D, Warning).
    #[allow(clippy::too_many_arguments)]
    fn new(
        _f: &MakeCurveFunction,
        _deg_min: i32,
        _deg_max: i32,
        _tol3d: f64,
        _tol2d: f64,
        _warning: bool,
    ) -> Self {
        panic!("GAP: Approx_FitAndDivide (TKGeomBase/Approx not translated)");
    }

    /// OCCT Approx_FitAndDivide::NbMultiCurves().
    fn nb_multi_curves(&self) -> i32 {
        panic!("GAP: Approx_FitAndDivide::NbMultiCurves (TKGeomBase/Approx not translated)");
    }

    /// OCCT Approx_FitAndDivide::Value(Index) -> AppParCurves_MultiCurve.
    fn value(&self, _index: i32) -> AppParCurvesMultiCurve {
        panic!("GAP: Approx_FitAndDivide::Value (TKGeomBase/Approx not translated)");
    }
}

/// OCCT AppParCurves_MultiCurve — GAP carrier (arch. diff. #26).
struct AppParCurvesMultiCurve;

impl AppParCurvesMultiCurve {
    /// OCCT AppParCurves_MultiCurve::Degree().
    fn degree(&self) -> i32 {
        panic!("GAP: AppParCurves_MultiCurve::Degree (TKGeomBase/AppParCurves not translated)");
    }

    /// OCCT AppParCurves_MultiCurve::Curve(Index, TPoints).
    fn curve(&self, _index: i32, _tpoints: &mut [DVec3]) {
        panic!("GAP: AppParCurves_MultiCurve::Curve (TKGeomBase/AppParCurves not translated)");
    }
}

/// OCCT Convert_CompBezierCurvesToBSplineCurve (TKGeomBase/Convert) — GAP
/// carrier (arch. diff. #26).
#[derive(Default)]
struct ConvertCompBezierCurvesToBSplineCurve;

impl ConvertCompBezierCurvesToBSplineCurve {
    /// OCCT ctor.
    fn new() -> Self {
        ConvertCompBezierCurvesToBSplineCurve
    }

    /// OCCT Convert_CompBezierCurvesToBSplineCurve::AddCurve(Poles).
    fn add_curve(&mut self, _poles: &[DVec3]) {
        panic!("GAP: Convert_CompBezierCurvesToBSplineCurve::AddCurve (TKGeomBase/Convert not translated)");
    }

    /// OCCT Convert_CompBezierCurvesToBSplineCurve::Perform().
    fn perform(&mut self) {
        panic!("GAP: Convert_CompBezierCurvesToBSplineCurve::Perform (TKGeomBase/Convert not translated)");
    }

    /// OCCT Convert_CompBezierCurvesToBSplineCurve::NbPoles().
    fn nb_poles(&self) -> i32 {
        panic!("GAP: Convert_CompBezierCurvesToBSplineCurve::NbPoles (TKGeomBase/Convert not translated)");
    }

    /// OCCT Convert_CompBezierCurvesToBSplineCurve::NbKnots().
    fn nb_knots(&self) -> i32 {
        panic!("GAP: Convert_CompBezierCurvesToBSplineCurve::NbKnots (TKGeomBase/Convert not translated)");
    }

    /// OCCT Convert_CompBezierCurvesToBSplineCurve::KnotsAndMults(SKnots,
    /// SMults).
    fn knots_and_mults(&self, _s_knots: &mut Vec<f64>, _s_mults: &mut Vec<i32>) {
        panic!("GAP: Convert_CompBezierCurvesToBSplineCurve::KnotsAndMults (TKGeomBase/Convert not translated)");
    }

    /// OCCT Convert_CompBezierCurvesToBSplineCurve::Poles(Poles).
    fn poles(&self, _poles: &mut Vec<DVec3>) {
        panic!("GAP: Convert_CompBezierCurvesToBSplineCurve::Poles (TKGeomBase/Convert not translated)");
    }

    /// OCCT Convert_CompBezierCurvesToBSplineCurve::Degree().
    fn degree(&self) -> i32 {
        panic!("GAP: Convert_CompBezierCurvesToBSplineCurve::Degree (TKGeomBase/Convert not translated)");
    }
}

// OCCT BiTgte_Blend.cxx L414-470 — MakeCurve.
fn make_curve(hc: &BiTgteCurveOnEdge) -> Option<Curve3> {
    let c: Option<Curve3>;

    if hc.get_type() == CurveType::Circle {
        // OCCT L425-426: C = new Geom_Circle(HC.Circle());
        // C = new Geom_TrimmedCurve(C, First, Last).
        let circle = Curve3::Circle(hc.circle());
        c = Some(Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(circle),
            first: hc.first_parameter(),
            last: hc.last_parameter(),
        }));
    } else {
        // the approximation is done
        // OCCT L430-436.
        let f = MakeCurveFunction::new(hc);
        let deg1: i32 = 8;
        let deg2: i32 = 8;
        let tol = PRECISION_APPROXIMATION;
        let fit = ApproxFitAndDivide::new(&f, deg1, deg2, tol, tol, true);
        let nb_curves = fit.nb_multi_curves();
        // it is attempted to make the curve at least C1
        let mut conv = ConvertCompBezierCurvesToBSplineCurve::new();

        for i in 1..=nb_curves {
            // OCCT L442-444: MC = Fit.Value(i); Poles(1, MC.Degree()+1);
            // MC.Curve(1, Poles).
            let mc = fit.value(i);
            let mut poles: Vec<DVec3> = vec![DVec3::ZERO; (mc.degree() + 1) as usize];
            mc.curve(1, &mut poles);

            // OCCT L446: Conv.AddCurve(Poles).
            conv.add_curve(&poles);
        }

        // OCCT L449: Conv.Perform().
        conv.perform();

        // OCCT L451-455.
        let nb_poles = conv.nb_poles();
        let nb_knots = conv.nb_knots();
        let mut new_poles: Vec<DVec3> = vec![DVec3::ZERO; nb_poles as usize];
        let mut new_knots: Vec<f64> = vec![0.0; nb_knots as usize];
        let mut new_mults: Vec<i32> = vec![0; nb_knots as usize];

        // OCCT L457-458.
        conv.knots_and_mults(&mut new_knots, &mut new_mults);
        conv.poles(&mut new_poles);

        // OCCT L460: BSplCLib::Reparametrize(First, Last, NewKnots).
        bspl_lib::reparametrize(hc.first_parameter(), hc.last_parameter(), &mut new_knots);

        // OCCT L462: C = new Geom_BSplineCurve(NewPoles, NewKnots, NewMults,
        // Conv.Degree()) — the rcad BSplineCurve3 carries the flat knot
        // vector, so the (knots, mults) pair is flattened here (arch.
        // diff. #26).
        let mut flat_knots: Vec<f64> = Vec::new();
        for (k, m) in new_knots.iter().zip(new_mults.iter()) {
            for _ in 0..*m {
                flat_knots.push(*k);
            }
        }
        let degree = conv.degree();
        c = Some(Curve3::BSpline(BSplineCurve3 {
            degree: degree as usize,
            knots: flat_knots,
            control_points: new_poles.clone(),
            weights: vec![1.0; new_poles.len()],
            is_periodic: false,
        }));
    }

    c
}

// OCCT BiTgte_Blend.cxx L472-496 — Touched.
// Only the faces connected with caps are given
fn touched(
    _analyse: &BRepOffsetAnalyse,
    _stop_faces: &HashSet<Shape>,
    _shape: &Shape,
    _touched_by_cork: &mut HashSet<Shape>,
) {
    // currently nothing is done !!
}

// OCCT BiTgte_Blend.cxx L500-543 — FindVertex.
fn find_vertex(
    p: DVec3,
    map: &HashSet<Shape>,
    tol: f64,
    the_brep: &mut BRep,
) -> Shape {
    // Find in <Map> a vertex which represents the point <P>.
    let mut b = BRepBuilder::new();
    let v = Shape::null();
    let mut vv = [Shape::null(), Shape::null()];
    let tol_carre = tol * tol;
    for it_key in map.iter() {
        let e = it_key;
        if !e.is_null() {
            // OCCT L515: TopExp::Vertices(E, VV[0], VV[1]).
            let (v1, v2) = top_exp_vertices_shape(e);
            vv[0] = v1;
            vv[1] = v2;

            for i in 0..2 {
                // if OK la Tolerance du Vertex
                // OCCT L520: Tol2 = BRep_Tool::Tolerance(VV[i]); Tol2 *=
                // Tol2.
                let mut tol2 = brep_tool_tolerance(&vv[i]);
                tol2 *= tol2;
                // OCCT L522-523: P1 = BRep_Tool::Pnt(VV[i]); Dist =
                // P.SquareDistance(P1).
                let p1 = brep_tool_pnt(&vv[i]).expect("FindVertex: null vertex point");
                let dist = p.distance_squared(p1);
                if dist <= tol2 {
                    return vv[i].clone();
                }
                // otherwise with the required tolerance.
                if tol_carre > tol2 {
                    if dist <= tol_carre {
                        // so it is necessary to update the tolerance of
                        // Vertex.
                        // OCCT L534: B.UpdateVertex(VV[i], Tol).
                        b.update_vertex_tolerance(the_brep, vv[i].clone(), tol);
                        return vv[i].clone();
                    }
                }
            }
        }
    }

    v
}

// OCCT BiTgte_Blend.cxx L547-583 — MakeDegeneratedEdge.
fn make_degenerated_edge(cc: &Curve3, vf_on_e: &Shape, the_brep: &mut BRep) -> Shape {
    let mut b = BRepBuilder::new();
    let tol = PRECISION_CONFUSION;
    // kill trimmed curves
    // OCCT L553-559: the while down-cast loop.
    let mut c: Curve3 = cc.clone();
    loop {
        if let Curve3::Trimmed(tc) = c {
            c = (*tc.curve).clone();
        } else {
            break;
        }
    }

    let v1: Shape;
    let v2: Shape;
    if vf_on_e.is_null() {
        // OCCT L564-566: P = C->Value(C->FirstParameter());
        // B.MakeVertex(V1, P, Tol); V2 = V1.
        use rcad_kernel::geom::CurveEval;
        let p = c.point_at(curve_first_parameter(&c));
        let made = b.add_vertex(the_brep, p, tol);
        v1 = made.clone();
        v2 = made;
    } else {
        // OCCT L570: V1 = V2 = VfOnE.
        v1 = vf_on_e.clone();
        v2 = vf_on_e.clone();
    }
    // OCCT L572-573.
    let mut v1 = v1;
    v1.orientation = Orientation::Forward;
    let mut v2 = v2;
    v2.orientation = Orientation::Reversed;

    // OCCT L575-578: B.MakeEdge(E, C, Tol); B.Add(E, V1); B.Add(E, V2) —
    // the rcad edge constructor binds the vertices together with the curve
    // (arch. diff. #29).
    let first = curve_first_parameter(&c);
    let last = curve_last_parameter(&c);
    let e = b.add_edge(the_brep, Some(c), v1, v2, [first, last]);
    // OCCT L581: B.Range(E, CC->FirstParameter(), CC->LastParameter()).
    b.set_edge_range(the_brep, e.clone(), first, last);
    e
}

/// OCCT TopAbs::Reverse (TopAbs.cxx L28-38) — the orientation reversal.
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}

/// OCCT TopoDS_Shape::Oriented(O) / Reversed() / Reverse() — the value-form
/// orientation setters (the rcad Shape carries a public orientation field).
fn oriented_shape(s: &Shape, o: Orientation) -> Shape {
    let mut out = s.clone();
    out.orientation = o;
    out
}

fn shape_reversed(s: &Shape) -> Shape {
    oriented_shape(s, top_abs_reverse(s.orientation))
}

// OCCT BiTgte_Blend.cxx L587-607 — Orientation (static).
fn orientation_of(e: &Shape, f: &Shape, l: &[Shape]) -> Orientation {
    let mut orien = Orientation::Forward;
    for itld in l {
        if itld.is_same(e) {
            orien = itld.orientation;
            break;
        }
    }
    if f.orientation == Orientation::Reversed {
        orien = top_abs_reverse(orien);
    }

    orien
}

// OCCT BiTgte_Blend.cxx L611-703 — FindCreatedEdge.
#[allow(clippy::too_many_arguments)]
fn find_created_edge(
    v1: &Shape,
    e: &Shape,
    map_sf: &HashMap<Shape, BRepOffsetOffset>,
    map_on_v: &mut HashSet<Shape>,
    center_analyse: &BRepOffsetAnalyse,
    radius: f64,
    tol: f64,
) -> Shape {
    let mut e1 = Shape::null();
    if !center_analyse.has_ancestor(v1) {
        return e1; // return a Null Shape.
    }

    // OCCT L626-627.
    let mut tang_e: Vec<Shape> = Vec::new();
    center_analyse.tangent_edges(e, v1, &mut tang_e);

    // OCCT L629-683.
    let mut find = false;
    for itl in &tang_e {
        if find {
            break;
        }
        let et = itl;
        // OCCT L634: MapSF.IsBound(ET).
        if let Some(of_et) = map_sf.get(et) {
            // OCCT L636-637: aLocalShape = MapSF(ET).Generated(V1); E1 =
            // TopoDS::Edge(aLocalShape).
            e1 = of_et.generated(v1);
            map_on_v.insert(e1.clone());
            find = true;
        } else {
            // Find the sharing of vertices in case of tangent consecutive
            // 3 edges the second of which is the edge that degenerates the
            // tube.
            // OCCT L646-648: CET = BRep_Tool::Curve(ET, CLoc, ff, ll) — the
            // location is carried by the rcad Shape (the CLoc out form is
            // not used further).
            let Some((mut cet, _ff, _ll)) = edge_curve_of(et) else {
                continue;
            };
            if let Curve3::Trimmed(tc) = cet {
                cet = (*tc.curve).clone();
            }
            // OCCT L653: Circ = down_cast<Geom_Circle>(CET).
            let circ = match &cet {
                Curve3::Circle(ci) => Some(*ci),
                _ => None,
            };
            let Some(circ) = circ else {
                continue;
            };
            if (circ.radius - radius.abs()).abs() > tol {
                continue;
            }

            // OCCT L663-664.
            let (mut u1, u2) = top_exp_vertices_shape(et);
            if u1.is_same(v1) {
                u1 = u2;
            }
            // OCCT L669-681.
            let mut tang2: Vec<Shape> = Vec::new();
            center_analyse.tangent_edges(et, &u1, &mut tang2);
            for it2 in &tang2 {
                let et2 = it2;
                if let Some(of_et2) = map_sf.get(et2) {
                    // OCCT L677-678.
                    let a_generated = of_et2.generated(&u1);
                    map_on_v.insert(a_generated);
                }
            }
        }
    }
    if !find {
        // OCCT L684-700.
        tang_e.clear();
        if center_analyse.has_ancestor(v1) {
            tang_e = center_analyse.ancestors(v1).clone();
            for itl in &tang_e {
                if find {
                    break;
                }
                if let Some(of_itl) = map_sf.get(itl) {
                    // OCCT L696: MapOnV.Add(MapSF(itl.Value()).Generated(V1)).
                    let a_generated = of_itl.generated(v1);
                    map_on_v.insert(a_generated);
                }
            }
        }
    }

    e1
}

// ===========================================================================
// OCCT class (BiTgte_Blend.hxx L32-150, cxx L709-2664).
// ===========================================================================

/// OCCT BiTgte_Blend (BiTgte_Blend.hxx L32-150).
pub struct BiTgteBlend {
    my_radius: f64,               // OCCT: myRadius (hxx L135)
    my_tol: f64,                  // OCCT: myTol (hxx L136)
    my_nubs: bool,                // OCCT: myNubs (hxx L137)
    my_shape: Shape,              // OCCT: myShape (hxx L138)
    my_result: Shape,             // OCCT: myResult (hxx L139)
    my_build_shape: bool,         // OCCT: myBuildShape (hxx L140)
    my_ancestors: AncestorsMap,   // OCCT: myAncestors (hxx L141)
    my_created: HashMap<Shape, HashMap<Shape, Vec<Shape>>>, // OCCT: myCreated (hxx L142-146)
    my_cut_edges: HashMap<Shape, Vec<Shape>>,               // OCCT: myCutEdges (hxx L147-149)
    my_faces: IndexMap<Shape, ()>,               // OCCT: myFaces (hxx L150)
    my_edges: IndexMap<Shape, ()>,               // OCCT: myEdges (hxx L151)
    my_stop_faces: HashSet<Shape>,               // OCCT: myStopFaces (hxx L152)
    my_analyse: BRepOffsetAnalyse,               // OCCT: myAnalyse (hxx L153)
    my_centers: IndexMap<Shape, ()>,             // OCCT: myCenters (hxx L154)
    my_map_sf: HashMap<Shape, BRepOffsetOffset>, // OCCT: myMapSF (hxx L155)
    my_init_offset_face: BRepAlgoImage,          // OCCT: myInitOffsetFace (hxx L156)
    my_image: BRepAlgoImage,                     // OCCT: myImage (hxx L157)
    my_image_offset: BRepAlgoImage,              // OCCT: myImageOffset (hxx L158)
    my_as_des: BRepAlgoAsDes,                    // OCCT: myAsDes (hxx L159, handle)
    my_nb_branches: i32,                         // OCCT: myNbBranches (hxx L160)
    my_indices: Vec<i32>,                        // OCCT: myIndices (hxx L161) — 1-based
                                                 // HArray1 carried by a 0-based Vec (arch. diff. #21)
    my_done: bool,                               // OCCT: myDone (hxx L162)
    my_brep: BRep,                               // rcad arena stand-in (arch. diff. #4/#19)
}

impl Default for BiTgteBlend {
    fn default() -> Self {
        Self::new()
    }
}

impl BiTgteBlend {
    /// OCCT BiTgte_Blend::BiTgte_Blend() (cxx L709-713).
    pub fn new() -> Self {
        BiTgteBlend {
            my_radius: 0.0,
            my_tol: 0.0,
            my_nubs: false,
            my_shape: Shape::null(),
            my_result: Shape::null(),
            my_build_shape: false,
            my_ancestors: IndexMap::new(),
            my_created: HashMap::new(),
            my_cut_edges: HashMap::new(),
            my_faces: IndexMap::new(),
            my_edges: IndexMap::new(),
            my_stop_faces: HashSet::new(),
            my_analyse: BRepOffsetAnalyse::new(),
            my_centers: IndexMap::new(),
            my_map_sf: HashMap::new(),
            my_init_offset_face: BRepAlgoImage::new(),
            my_image: BRepAlgoImage::new(),
            my_image_offset: BRepAlgoImage::new(),
            my_as_des: BRepAlgoAsDes::new(),
            my_nb_branches: -1,
            my_indices: Vec::new(),
            my_done: false,
            my_brep: BRep::new(),
        }
    }

    /// OCCT BiTgte_Blend::BiTgte_Blend(S, Radius, Tol, NUBS) (cxx
    /// L717-724).
    pub fn with_shape(s: &Shape, radius: f64, tol: f64, nubs: bool) -> Self {
        let mut res = BiTgteBlend::new();
        res.init(s, radius, tol, nubs);
        res
    }

    /// OCCT BiTgte_Blend::Init(S, Radius, Tol, NUBS) (cxx L728-740).
    pub fn init(&mut self, s: &Shape, radius: f64, tol: f64, nubs: bool) {
        self.clear();
        self.my_shape = s.clone();
        self.my_tol = tol;
        self.my_nubs = nubs;
        self.my_radius = radius;
        self.my_nb_branches = -1;
        //  TopExp::MapShapesAndAncestors(S,TopAbs_EDGE,TopAbs_FACE,myAncestors);
    }

    /// OCCT BiTgte_Blend::Clear() (cxx L744-754) — Clear all the Fields.
    pub fn clear(&mut self) {
        self.my_init_offset_face.clear();
        self.my_image.clear();
        self.my_image_offset.clear();
        self.my_stop_faces.clear();
        self.my_analyse.clear();
        self.my_as_des.clear();
        self.my_nb_branches = -1;
        self.my_done = false;
    }

    /// OCCT BiTgte_Blend::SetStoppingFace(Face) (cxx L758-767) — Set a face
    /// on which the fillet must stop.
    pub fn set_stopping_face(&mut self, face: &Shape) {
        self.my_stop_faces.insert(face.clone());
        //-------------
        // MAJ SD. -> To end loop, set faces of edges
        //-------------
        //  myInitOffsetFace.SetRoot(Face);
        //  myInitOffsetFace.Bind   (Face,Face);
        //  myImageOffset.SetRoot   (Face);
    }

    /// OCCT BiTgte_Blend::SetFaces(F1, F2) (cxx L771-775) — Set two faces
    /// of <myShape> on which the Sphere must roll.
    pub fn set_faces(&mut self, f1: &Shape, f2: &Shape) {
        self.my_faces.insert(f1.clone(), ());
        self.my_faces.insert(f2.clone(), ());
    }

    /// OCCT BiTgte_Blend::SetEdge(Edge) (cxx L779-782) — Set an edge of
    /// <myShape> to be rounded.
    pub fn set_edge(&mut self, edge: &Shape) {
        self.my_edges.insert(edge.clone(), ());
    }

    /// OCCT BiTgte_Blend::Perform(BuildShape) (cxx L786-940) — Compute the
    /// generated surfaces.
    pub fn perform(&mut self, build_shape: bool) {
        self.my_build_shape = build_shape;

        // Try cutting to avoid tubes on free borders
        // that are not actually free.
        // OCCT L792: Sew = new BRepBuilderAPI_Sewing(myTol).
        let mut sew = BRepBuilderAPISewing::new(self.my_tol);
        // OCCT L793: BRepLib::BuildCurves3d(myShape) (BRepLib.cxx L460-464).
        crate::topalgo::brep_lib::brep_lib::BRepLib::build_curves3d(
            &mut self.my_brep,
            &self.my_shape,
        );
        // OCCT L794-799.
        for a_face in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            sew.add(&a_face);
        }
        sew.perform();
        let mut sewed_shape = sew.sewed_shape();
        if sewed_shape.is_null() {
            // OCCT L801-804: throw Standard_Failure("Sewing aux fraises").
            panic!("Failure: Sewing aux fraises");
        }

        // Check if the sewing modified the orientation.
        // OCCT L807-823.
        let expf = explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape);
        let mut face_ref = expf
            .first()
            .cloned()
            .expect("BiTgte_Blend::Perform: no face in myShape");
        let ori_ref = face_ref.orientation;
        if sew.is_modified(&face_ref) {
            face_ref = sew.modified(&face_ref);
        }
        for ff in explorer(&sewed_shape, ShapeType::Face, ShapeType::Shape) {
            if face_ref.is_same(&ff) && ff.orientation != ori_ref {
                sewed_shape = shape_reversed(&sewed_shape);
                break;
            }
        }

        // Make SameParameter if Sew does not do it (Detect that edges
        // are not sameparameter but it does nothing.)
        // OCCT L827-832.
        for sec in explorer(&sewed_shape, ShapeType::Edge, ShapeType::Shape) {
            let tol = brep_tool_tolerance(&sec);
            brep_lib_same_parameter(&sec, tol);
        }

        // OCCT L834: TopExp::MapShapesAndAncestors(SewedShape, EDGE, FACE,
        // myAncestors) — the shared re-host.
        self.my_ancestors = IndexMap::new();
        map_shapes_and_ancestors(
            &sewed_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut self.my_ancestors,
        );

        // Extend myFaces with the faces of the sewed shape.
        // OCCT L836-846.
        for f in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            if self.my_faces.contains_key(&f) && sew.is_modified(&f) {
                self.my_faces.shift_remove(&f);
                let modified = sew.modified(&f);
                self.my_faces.insert(modified, ());
            }
        }

        self.my_shape = sewed_shape;
        // end Sewing for false free borders.

        // ----------------------------------------------------------------
        // place faces with the proper orientation in the initial shape
        // ----------------------------------------------------------------
        // OCCT L866-880.
        for f in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            if self.my_faces.contains_key(&f) {
                self.my_faces.shift_remove(&f);
                self.my_faces.insert(f, ());
            } else if self.my_stop_faces.contains(&f) {
                self.my_stop_faces.remove(&f);
                self.my_stop_faces.insert(f);
            }
        }

        // ----------------------------------------------
        // Calculate lines of centers and of surfaces
        // ----------------------------------------------

        self.compute_centers();

        // -----------------------------
        // Calculate connection Surfaces
        // -----------------------------

        self.compute_surfaces();

        // ----------------------------------
        // Calculate the generated shape if required
        // ----------------------------------

        if self.my_build_shape {
            self.compute_shape();
        }

        // Finally construct curves 3d from edges to be transferred
        // since the partition is provided ( A Priori);
        // OCCT L926: BRepLib::BuildCurves3d(myResult, Precision::Confusion())
        // (BRepLib.cxx L468-489).
        crate::topalgo::brep_lib::brep_lib::BRepLib::build_curves3d_tol(
            &mut self.my_brep,
            &self.my_result,
            PRECISION_CONFUSION,
        );

        self.my_done = true;
    }

    /// OCCT BiTgte_Blend::IsDone() (cxx L944-947).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BiTgte_Blend::Shape() (cxx L951-954) — returns the result.
    pub fn shape(&self) -> &Shape {
        &self.my_result
    }

    /// OCCT BiTgte_Blend::NbSurfaces() (cxx L958-961) — returns the Number
    /// of generated surfaces.
    pub fn nb_surfaces(&self) -> i32 {
        self.my_centers.len() as i32
    }

    /// OCCT BiTgte_Blend::Surface(Index) (cxx L965-968).
    pub fn surface(&self, index: i32) -> Option<Surface3> {
        self.surface_by_center_line(&self.center_at(index))
    }

    /// OCCT BiTgte_Blend::Face(Index) (cxx L975-978).
    pub fn face(&self, index: i32) -> Shape {
        self.face_by_center_line(&self.center_at(index))
    }

    /// OCCT myCenters(Index) — the 1-based IndexedMap read.
    fn center_at(&self, index: i32) -> Shape {
        self.my_centers
            .get_index((index - 1) as usize)
            .expect("BiTgte_Blend: center index out of range")
            .0
            .clone()
    }

    /// OCCT BiTgte_Blend::CenterLines(LC) (cxx L982-990) — set in <LC> all
    /// the center lines.
    pub fn center_lines(&self, lc: &mut Vec<Shape>) {
        lc.clear();
        let nb = self.nb_surfaces();
        for i in 1..=nb {
            lc.push(self.center_at(i));
        }
    }

    /// OCCT BiTgte_Blend::Surface(CenterLine) (cxx L994-998) — returns the
    /// surface generated by the centerline (Null if none).
    pub fn surface_by_center_line(&self, center_line: &Shape) -> Option<Surface3> {
        // OCCT L996: F = myMapSF(CenterLine).Face(); BRep_Tool::Surface(F).
        let f = self
            .my_map_sf
            .get(center_line)
            .expect("BiTgte_Blend::Surface: center line not in myMapSF")
            .face();
        Self::brep_tool_surface(&f)
    }

    /// OCCT BiTgte_Blend::Face(CenterLine) (cxx L1005-1013).
    pub fn face_by_center_line(&self, center_line: &Shape) -> Shape {
        if !self.my_map_sf.contains_key(center_line) {
            // OCCT L1009: throw Standard_DomainError("BiTgte_Blend::Face").
            panic!("DomainError: BiTgte_Blend::Face");
        }

        self.my_map_sf[center_line].face()
    }

    /// OCCT BiTgte_Blend::ContactType(Index) (cxx L1017-1082).
    pub fn contact_type(&self, index: i32) -> BiTgteContactType {
        let s1 = self.support_shape1(index);
        let s2 = self.support_shape2(index);

        let mut type1 = s1.shape_type();
        let mut type2 = s2.shape_type();

        // OCCT L1025-1030: the TopAbs ordinal compare — the rcad
        // ShapeType declaration order differs from the OCCT TopAbs
        // enum, so the OCCT rank (FACE=5 < EDGE=6 < VERTEX=7) is mapped
        // explicitly (arch. diff. #21).
        if top_abs_rank(type2) < top_abs_rank(type1) {
            std::mem::swap(&mut type1, &mut type2);
        }
        let mut type_ = BiTgteContactType::VertexVertex;

        match type1 {
            ShapeType::Vertex => match type2 {
                ShapeType::Vertex => type_ = BiTgteContactType::VertexVertex,
                ShapeType::Edge => type_ = BiTgteContactType::EdgeVertex,
                ShapeType::Face => type_ = BiTgteContactType::FaceVertex,
                _ => {}
            },

            ShapeType::Edge => match type2 {
                ShapeType::Edge => type_ = BiTgteContactType::EdgeEdge,
                ShapeType::Face => type_ = BiTgteContactType::FaceEdge,
                _ => {}
            },

            ShapeType::Face => match type2 {
                ShapeType::Face => type_ = BiTgteContactType::FaceEdge,
                _ => {}
            },

            _ => {}
        }

        type_
    }

    /// OCCT BiTgte_Blend::SupportShape1(Index) (cxx L1086-1098).
    pub fn support_shape1(&self, index: i32) -> &Shape {
        let cur_e = self.center_at(index);

        // OCCT L1090: L = myAsDes->Ascendant(CurE).
        let l = self.my_as_des.ascendant(&cur_e);

        // --------------------------------------------------------------
        // F1 and F2 = 2 parallel faces intersecting at CurE.
        // --------------------------------------------------------------
        let f1 = &l[0];
        // OCCT L1096: Or1 = myInitOffsetFace.ImageFrom(F1).
        self.my_init_offset_face.image_from(f1)
    }

    /// OCCT BiTgte_Blend::SupportShape2(Index) (cxx L1102-1114).
    pub fn support_shape2(&self, index: i32) -> &Shape {
        let cur_e = self.center_at(index);

        // OCCT L1106: L = myAsDes->Ascendant(CurE).
        let l = self.my_as_des.ascendant(&cur_e);

        // --------------------------------------------------------------
        // F1 and F2 = 2 parallel faces intersecting at CurE.
        // --------------------------------------------------------------
        let f2 = l.last().expect("BiTgte_Blend: empty ascendant list");
        // OCCT L1112: Or2 = myInitOffsetFace.ImageFrom(F2).
        self.my_init_offset_face.image_from(f2)
    }

    /// OCCT BiTgte_Blend::CurveOnShape1(Index) (cxx L1118-1136).
    pub fn curve_on_shape1(&self, index: i32) -> Option<Curve3> {
        let cur_e = self.center_at(index);
        let f = self
            .my_map_sf
            .get(&cur_e)
            .expect("BiTgte_Blend::CurveOnShape1: CurE not in myMapSF")
            .face();

        // somewhat brutal method based ONLY on the construction of the
        // fillet: the first edge of the tube is exactly the edge on Shape1.

        // OCCT L1126-1127: exp.Current() — the first edge.
        let e = explorer(&f, ShapeType::Edge, ShapeType::Shape)
            .into_iter()
            .next()
            .expect("BiTgte_Blend::CurveOnShape1: no edge in tube face");
        let mut c: Option<Curve3> = None;
        if !brep_tool_degenerated(&e) {
            // OCCT L1132-1133: C = BRep_Tool::Curve(E, f, l); C = new
            // Geom_TrimmedCurve(C, f, l).
            if let Some((ec, pf, pl)) = edge_curve_of(&e) {
                c = Some(Curve3::Trimmed(TrimmedCurve3 {
                    curve: Box::new(ec),
                    first: pf,
                    last: pl,
                }));
            }
        }
        c
    }

    /// OCCT BiTgte_Blend::CurveOnShape2(Index) (cxx L1140-1159).
    pub fn curve_on_shape2(&self, index: i32) -> Option<Curve3> {
        let cur_e = self.center_at(index);
        let f = self
            .my_map_sf
            .get(&cur_e)
            .expect("BiTgte_Blend::CurveOnShape2: CurE not in myMapSF")
            .face();

        // somewhat brutal method based ONLY on the construction of the
        // fillet: the first edge of the tube is exactly the edge on Shape2.

        // OCCT L1148-1150: exp.Next(); E = exp.Current() — the second edge.
        let mut edges = explorer(&f, ShapeType::Edge, ShapeType::Shape).into_iter();
        let _first = edges.next();
        let e = edges
            .next()
            .expect("BiTgte_Blend::CurveOnShape2: no second edge in tube face");
        let mut c: Option<Curve3> = None;
        if !brep_tool_degenerated(&e) {
            // OCCT L1155-1156.
            if let Some((ec, pf, pl)) = edge_curve_of(&e) {
                c = Some(Curve3::Trimmed(TrimmedCurve3 {
                    curve: Box::new(ec),
                    first: pf,
                    last: pl,
                }));
            }
        }
        c
    }

    /// OCCT BiTgte_Blend::PCurveOnFace1(Index) (cxx L1163-1167).
    pub fn pcurve_on_face1(&self, _index: i32) -> Option<Curve2d> {
        let c: Option<Curve2d> = None;
        c
    }

    /// OCCT BiTgte_Blend::PCurve1OnFillet(Index) (cxx L1171-1175).
    pub fn pcurve1_on_fillet(&self, _index: i32) -> Option<Curve2d> {
        let c: Option<Curve2d> = None;
        c
    }

    /// OCCT BiTgte_Blend::PCurveOnFace2(Index) (cxx L1179-1183).
    pub fn pcurve_on_face2(&self, _index: i32) -> Option<Curve2d> {
        let c: Option<Curve2d> = None;
        c
    }

    /// OCCT BiTgte_Blend::PCurve2OnFillet(Index) (cxx L1187-1191).
    pub fn pcurve2_on_fillet(&self, _index: i32) -> Option<Curve2d> {
        let c: Option<Curve2d> = None;
        c
    }

    /// OCCT BiTgte_Blend::NbBranches() (cxx L1195-1270).
    pub fn nb_branches(&mut self) -> i32 {
        if self.my_nb_branches != -1 {
            return self.my_nb_branches;
        }

        // else, compute the Branches.
        // OCCT L1203: BRepTools_Quilt Glue — GAP carrier (arch. diff. #24
        // family; the TKTopAlgo quilt is not translated).
        let mut glue = BRepToolsQuilt::new();

        let nb_faces = self.my_centers.len() as i32;

        if nb_faces == 0 {
            return 0;
        }

        // OCCT L1213-1217.
        for i in 1..=nb_faces {
            let center_line = self.center_at(i);
            let face = self.my_map_sf[&center_line].face();
            glue.add(&face);
        }

        // OCCT L1219: Shells = Glue.Shells().
        let shells = glue.shells();

        // Reorder Map myCenters.
        // The method is brutal and unpolished,
        // it is possible to refine it.
        self.my_nb_branches = 0;
        let mut tmp_map: IndexMap<Shape, ()> = IndexMap::new();

        // OCCT L1227-1231.
        for _cur_s in explorer(&shells, ShapeType::Shell, ShapeType::Shape) {
            self.my_nb_branches += 1;
        }

        // OCCT L1233: myIndices = new NCollection_HArray1<int>(1,
        // myNbBranches + 1) — the 1-based HArray1 carried by a 0-based Vec
        // of the same length (arch. diff. #21).
        self.my_indices = vec![0; (self.my_nb_branches + 1) as usize];

        // OCCT L1235: myIndices->SetValue(1, 0).
        self.my_indices[0] = 0;
        let mut count = 0;
        let mut index: usize = 2;

        for cur_s in explorer(&shells, ShapeType::Shell, ShapeType::Shape) {
            // CurS = the current Shell.

            // OCCT L1245-1263.
            for cur_f in explorer(&cur_s, ShapeType::Face, ShapeType::Shape) {
                // CurF = the current face of the current Shell.

                for i in 1..=nb_faces {
                    let center = self.center_at(i);
                    let rakk = self.my_map_sf[&center].face();
                    // Rakk = the ith generated connection face
                    if cur_f.is_equal(&rakk) {
                        tmp_map.insert(center, ());
                        count += 1;
                        break;
                    }
                }
            }
            // OCCT L1264: myIndices->SetValue(Index, Count).
            self.my_indices[index - 1] = count;
            index += 1;
        }

        self.my_centers = tmp_map;
        self.my_nb_branches
    }

    /// OCCT BiTgte_Blend::IndicesOfBranche(Index, From, To) (cxx
    /// L1274-1281) — Set in <From>,<To> the indices of the faces of the
    /// branche <Index>.
    ///
    /// i.e: Branche<Index> = Face(From) + Face(From+1) + ..+ Face(To)
    pub fn indices_of_branche(&self, index: i32, from: &mut i32, to: &mut i32) {
        // Attention to the ranking in myIndices:
        // If the branches are  1-4 5-9 10-12, it is ranked in myIndices:
        //                      0 4   9    12
        // OCCT L1279-1280 — the 1-based HArray1 reads via the 0-based Vec.
        *from = self.my_indices[(index - 1) as usize] + 1;
        *to = self.my_indices[index as usize];
    }

    /// OCCT BRep_Tool::Surface(F) — the rcad face surface read (the
    /// TFaceData::surface).
    fn brep_tool_surface(f: &Shape) -> Option<Surface3> {
        f.as_face().and_then(|fd| fd.surface.clone())
    }
}

/// OCCT TopAbs_ShapeEnum rank (COMPOUND=0, COMPSOLID=1, SOLID=2, SHELL=3,
/// WIRE=4, FACE=5, EDGE=6, VERTEX=7, SHAPE=8) — the rcad ShapeType
/// declaration order differs (arch. diff. #21).
fn top_abs_rank(t: ShapeType) -> i32 {
    match t {
        ShapeType::Compound => 0,
        ShapeType::CompSolid => 1,
        ShapeType::Solid => 2,
        ShapeType::Shell => 3,
        ShapeType::Wire => 4,
        ShapeType::Face => 5,
        ShapeType::Edge => 6,
        ShapeType::Vertex => 7,
        ShapeType::Shape => 8,
    }
}

// The member-function half (ComputeCenters / ComputeSurfaces / ComputeShape
// / Intersect) — the child-module form keeps the private class fields
// visible to the implementation split (the OCCT members are private to the
// class).
#[path = "bi_tgte_blended_b.rs"]
mod blended_b;
