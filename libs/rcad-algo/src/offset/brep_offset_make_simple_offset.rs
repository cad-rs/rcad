// OCCT BRepOffset_MakeSimpleOffset.cxx L1-703 + BRepOffset_MakeSimpleOffset.hxx
// L30-172 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_MakeSimpleOffset.cxx / .hxx
//
// OCCT inheritance chain (hxx L61): none — BRepOffset_MakeSimpleOffset is a
// standalone class.
//
// Architecture differences:
// 1. TopoDS_Shape -> rcad Shape (Arc<TShape> handle + location + orientation);
//    NCollection_DataMap<TopoDS_Vertex, TopoDS_Edge> (myMapVE) -> HashMap
//    keyed by (TShape ptr, Location) — the TopTools_ShapeMapHasher identity
//    (orientation ignored), the same map form as feat/loc_ope_wires_on_shape.
// 2. BRepTools_Modifier + BRepOffset_SimpleOffset (TKOffset/BRepOffset) — the
//    rebuild engine; the BRepToolsModifier / BRepOffsetSimpleOffset carriers
//    below keep the OCCT constructor/accessor surface with GAP panics (port
//    plan §0.6).  Everything around the engine calls is translated 1:1.
// 3. ShapeAnalysis_FreeBounds (TKShHealing), BRepTools_Quilt (TKTopAlgo),
//    ShapeFix_Edge::FixSameParameter (TKShHealing), GeomFill_Generator
//    (TKGeomAlgo), BRepLib::BuildCurves3d and the planar
//    BRepLib_MakeFace(W, OnlyPlane) constructor have no rcad translation yet
//    — GAP carriers below.
// 4. OCCT mutates TShapes in place through a global arena; rcad carries the
//    arena as the BRep pool.  The class holds my_brep (the pool stand-in for
//    the arena; consumed by BRepBuilder mutations and the rcad
//    ShapeBuildReShape API, which takes &mut BRep).  The pool-local stand-in
//    for `BRep_Builder aBB;` locals is a BRepBuilder over that pool.
// 5. BRepLib_MakeEdge(V1, V2) — the rcad BRepBuilder::add_edge construction
//    cannot fail, so the OCCT IsDone() guards (cxx L549/L570) have no rcad
//    counterpart (annotated at the call sites).
// 6. BRepAdaptor_Surface(F, false)::D1 — the rcad stand-in evaluates the face
//    surface directly through Surface3::dn (the GeomAdaptor_Surface::DN
//    vehicle, geom/eval.rs); the BRepAdaptor_Surface wrapper needs the BRep
//    pool and carries no additional state for analytic surfaces.
// 7. The rcad ShapeBuildReShape API takes &mut self and &mut BRep
//    (shhealing/shape_build/reshape.rs); the OCCT const Generated/Modified
//    and the myReShape->Apply calls take the class by &mut self.

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, Line2d, Surface3, SurfaceEval, TrimmedCurve3,
};
use rcad_kernel::topo::topods::{tshape_flags, BRep, BRepBuilder, GeomAbsShape, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::feat::loc_ope_wires_on_shape::{shape_key, top_exp_vertices};
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_range,
    brep_tool_surface, brep_tool_tolerance, ShapeKey,
};
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;

// ---------------------------------------------------------------------------
// OCCT BRepOffsetSimple_Status (hxx L30-38).
// ---------------------------------------------------------------------------

/// OCCT BRepOffsetSimple_Status (BRepOffset_MakeSimpleOffset.hxx L30-38).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepOffsetSimpleStatus {
    Ok,
    NullInputShape,
    ErrorOffsetComputation,
    ErrorWallFaceComputation,
    ErrorInvalidNbShells,
    ErrorNonClosedShell,
}

// ---------------------------------------------------------------------------
// GAP carriers (architecture differences #2/#3).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_SimpleOffset (BRepOffset_SimpleOffset.hxx L44-192;
/// BRepOffset_SimpleOffset.cxx L1-427) — the BRepTools_Modification mapper of
/// the simple offset algorithm (architecture difference #2; GAP: the
/// NewSurface/NewCurve/NewPoint/NewCurve2d/NewParameter computations have no
/// rcad translation yet — the GAP panics are the §0.6 annotation; the
/// constructor keeps the OCCT storage form).
pub struct BRepOffsetSimpleOffset {
    my_input_shape: Shape, // OCCT: myInputShape
    my_offset_value: f64,  // OCCT: myOffsetValue
    my_tolerance: f64,     // OCCT: myTolerance
}

impl BRepOffsetSimpleOffset {
    /// OCCT BRepOffset_SimpleOffset::BRepOffset_SimpleOffset(theInputShape,
    /// theOffsetValue, theTolerance) (BRepOffset_SimpleOffset.cxx L58-66).
    pub fn new(the_input_shape: &Shape, the_offset_value: f64, the_tolerance: f64) -> Self {
        BRepOffsetSimpleOffset {
            my_input_shape: the_input_shape.clone(),
            my_offset_value: the_offset_value,
            my_tolerance: the_tolerance,
        }
    }

    /// OCCT BRepOffset_SimpleOffset::NewSurface — GAP.
    pub fn new_surface(&self, the_f: &Shape) -> Option<(Surface3, f64, bool, bool)> {
        let _ = the_f;
        panic!("GAP: BRepOffset_SimpleOffset::NewSurface (BRepOffset_SimpleOffset not translated)");
    }

    /// OCCT BRepOffset_SimpleOffset::NewCurve — GAP.
    pub fn new_curve(&self, the_e: &Shape) -> Option<(Curve3, f64, bool)> {
        let _ = the_e;
        panic!("GAP: BRepOffset_SimpleOffset::NewCurve (BRepOffset_SimpleOffset not translated)");
    }

    /// OCCT BRepOffset_SimpleOffset::NewPoint — GAP.
    pub fn new_point(&self, the_v: &Shape) -> Option<(DVec3, f64)> {
        let _ = the_v;
        panic!("GAP: BRepOffset_SimpleOffset::NewPoint (BRepOffset_SimpleOffset not translated)");
    }

    /// OCCT BRepOffset_SimpleOffset::NewCurve2d — GAP.
    pub fn new_curve2d(
        &self,
        the_e: &Shape,
        the_f: &Shape,
        the_new_e: &Shape,
        the_new_f: &Shape,
    ) -> Option<(Curve2d, f64)> {
        let _ = (the_e, the_f, the_new_e, the_new_f);
        panic!("GAP: BRepOffset_SimpleOffset::NewCurve2d (BRepOffset_SimpleOffset not translated)");
    }

    /// OCCT BRepOffset_SimpleOffset::NewParameter — GAP.
    pub fn new_parameter(&self, the_v: &Shape, the_e: &Shape) -> Option<(f64, f64)> {
        let _ = (the_v, the_e);
        panic!(
            "GAP: BRepOffset_SimpleOffset::NewParameter (BRepOffset_SimpleOffset not translated)"
        );
    }

    /// OCCT BRepOffset_SimpleOffset::Continuity — GAP.
    pub fn continuity(&self, the_e: &Shape, the_f1: &Shape, the_f2: &Shape) -> GeomAbsShape {
        let _ = (the_e, the_f1, the_f2);
        panic!("GAP: BRepOffset_SimpleOffset::Continuity (BRepOffset_SimpleOffset not translated)");
    }
}

/// OCCT BRepTools_Modifier (TKTopAlgo/BRepTools, BRepTools_Modifier.hxx) —
/// the shape rebuild vehicle driven by a BRepTools_Modification (architecture
/// difference #2; GAP: the Perform/ModifiedShape engine has no rcad
/// translation yet — the GAP panics are the §0.6 annotation; Init keeps the
/// OCCT storage form; the loc_ope_prism.rs precedent carries the same GAP for
/// the BRepTools_TrsfModification specialization).
///
/// The carrier owns the BRep pool of the produced shapes (architecture
/// difference #4): OCCT BRep_Builder edits TShapes in the global arena, rcad
/// needs the pool handle for the same edits.
pub struct BRepToolsModifier {
    my_shape: Shape,  // OCCT: myShape
    my_brep: BRep,    // rcad arena stand-in (arch. diff. #4)
    my_is_done: bool, // OCCT: myIsDone
}

impl BRepToolsModifier {
    /// OCCT BRepTools_Modifier::BRepTools_Modifier().
    pub fn new() -> Self {
        BRepToolsModifier {
            my_shape: Shape::null(),
            my_brep: BRep::new(),
            my_is_done: false,
        }
    }

    /// OCCT BRepTools_Modifier::Init(S).
    pub fn init(&mut self, the_shape: &Shape) {
        self.my_shape = the_shape.clone();
    }

    /// OCCT BRepTools_Modifier::Perform(M) — GAP.
    pub fn perform(&mut self, the_m: &BRepOffsetSimpleOffset) {
        let _ = the_m;
        panic!("GAP: BRepTools_Modifier::Perform (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_Modifier::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepTools_Modifier::ModifiedShape(S) — GAP.
    pub fn modified_shape(&self, the_shape: &Shape) -> Shape {
        let _ = the_shape;
        panic!("GAP: BRepTools_Modifier::ModifiedShape (TKTopAlgo/BRepTools not translated)");
    }

    /// The rcad arena stand-in of the produced shapes (arch. diff. #4).
    pub fn brep_pool_mut(&mut self) -> &mut BRep {
        &mut self.my_brep
    }
}

/// OCCT ShapeAnalysis_FreeBounds (TKShHealing) — the free-bounds explorer of
/// BuildMissingWalls (architecture difference #3; GAP: no rcad translation
/// yet — the GAP panic is the §0.6 annotation; GetClosedWires keeps the OCCT
/// accessor surface).
pub struct ShapeAnalysisFreeBounds {
    my_closed_wires: Shape, // OCCT: myClosedWires
}

impl ShapeAnalysisFreeBounds {
    /// OCCT ShapeAnalysis_FreeBounds::ShapeAnalysis_FreeBounds(theShape)
    /// (defaults: theSewConnected = false, theShared = false, theSetProjPCur
    /// = false).  GAP: the free-bounds computation is not translated.
    pub fn new(_the_shape: &Shape) -> Self {
        panic!("GAP: ShapeAnalysis_FreeBounds (TKShHealing not translated)");
    }

    /// OCCT ShapeAnalysis_FreeBounds::GetClosedWires().
    pub fn get_closed_wires(&self) -> Shape {
        self.my_closed_wires.clone()
    }
}

/// OCCT BRepTools_Quilt (TKTopAlgo/BRepTools) — the sewing of
/// BuildMissingWalls (architecture difference #3; GAP: no rcad translation
/// yet — the GAP panics are the §0.6 annotation).
pub struct BRepToolsQuilt;

impl BRepToolsQuilt {
    /// OCCT BRepTools_Quilt::Add(S).
    pub fn add(&mut self, the_s: &Shape) {
        let _ = the_s;
        panic!("GAP: BRepTools_Quilt::Add (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_Quilt::Shells().
    pub fn shells(&self) -> Shape {
        panic!("GAP: BRepTools_Quilt::Shells (TKTopAlgo/BRepTools not translated)");
    }
}

/// OCCT ShapeFix_Edge (TKShHealing) — the FixSameParameter of
/// BuildMissingWalls (architecture difference #3; GAP: no rcad translation
/// yet; the context argument keeps the OCCT SetContext form).
fn shape_fix_edge_fix_same_parameter(the_context: &mut ShapeBuildReShape, the_e: &Shape) {
    let _ = (the_context, the_e);
    panic!("GAP: ShapeFix_Edge::FixSameParameter (TKShHealing not translated)");
}

/// OCCT BRepLib::BuildCurves3d(S) (TKTopAlgo/BRepTools, BRepLib.cxx) — GAP:
/// the 3d-curve rebuild walk is not translated yet.
fn brep_lib_build_curves3d(the_s: &Shape) {
    let _ = the_s;
    panic!("GAP: BRepLib::BuildCurves3d (TKTopAlgo/BRepTools not translated)");
}

/// OCCT GeomFill_Generator (TKGeomAlgo/GeomFill) — the thrusection generator
/// of BuildWallFace (architecture difference #3; GAP: no rcad translation yet
/// — the GAP panics are the §0.6 annotation).
pub struct GeomFillGenerator;

impl GeomFillGenerator {
    /// OCCT GeomFill_Generator::AddCurve(Curve).
    pub fn add_curve(&mut self, the_curve: &TrimmedCurve3) {
        let _ = the_curve;
        panic!("GAP: GeomFill_Generator::AddCurve (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Generator::Perform(Pres3d).
    pub fn perform(&mut self, the_pres3d: f64) {
        let _ = the_pres3d;
        panic!("GAP: GeomFill_Generator::Perform (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Generator::Surface().
    pub fn surface(&self) -> Surface3 {
        panic!("GAP: GeomFill_Generator::Surface (TKGeomAlgo/GeomFill not translated)");
    }
}

/// OCCT BRepLib_MakeFace(W, OnlyPlane) (TKTopAlgo/BRepLib_MakeFace) — the
/// planar face maker of BuildWallFace (architecture difference #3; GAP: the
/// planar-surface fitting is not translated yet; the carrier keeps the OCCT
/// constructor/IsDone/Face surface with IsDone() = false so the caller takes
/// the OCCT failure branch, the loc_ope_wires_on_shape_b.rs #9 precedent).
pub struct BRepLibMakeFace {
    my_face: Shape, // OCCT: myFace
}

impl BRepLibMakeFace {
    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const TopoDS_Wire& W, const
    /// bool OnlyPlane).
    pub fn from_wire(_the_w: &Shape, _the_only_plane: bool) -> Self {
        // GAP: the planar face maker (BRepLib_FindSurface vehicle) is not
        // translated; IsDone() = false reproduces the OCCT failure path.
        BRepLibMakeFace {
            my_face: Shape::null(),
        }
    }

    /// OCCT BRepLib_MakeFace::IsDone().
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT BRepLib_MakeFace::Face().
    pub fn face(&self) -> Shape {
        self.my_face.clone()
    }
}

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp re-hosts (module-private; the loc_ope_* precedents).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::MaxTolerance(theShape, theSubShape)
/// (BRep_Tool.cxx L1792-1830) — the max tolerance over the explored
/// sub-shapes (Face/Edge/Vertex branches).
fn brep_tool_max_tolerance(the_shape: &Shape, the_sub_shape: ShapeType) -> f64 {
    let mut a_tol: f64 = 0.0;

    // Explorer Shape-Subshape.
    let an_exp_ss = explorer(the_shape, the_sub_shape, ShapeType::Shape);
    if the_sub_shape == ShapeType::Face
        || the_sub_shape == ShapeType::Edge
        || the_sub_shape == ShapeType::Vertex
    {
        for a_current_sub_shape in &an_exp_ss {
            a_tol = a_tol.max(brep_tool_tolerance(a_current_sub_shape));
        }
    }

    a_tol
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) (BRep_Tool.cxx L1180-1188 -> the
/// representation scan) — the regularity stored on the (E, F1/F2) curve
/// representations; the default is GeomAbs_C0 when nothing matches
/// (hlr/brep/shape_to_hlr.rs precedent, reduced to the identity-location rcad
/// form).
fn brep_tool_continuity(the_e: &Shape, the_f1: &Shape, the_f2: &Shape) -> GeomAbsShape {
    use rcad_kernel::topo::topods::CurveRepresentation;

    // OCCT: S1 = Surface(F1); S2 = Surface(F2) — the local surfaces.
    let (Some(s1), Some(s2)) = (brep_tool_surface(the_f1), brep_tool_surface(the_f2)) else {
        // The OCCT null-surface case has no matching representation.
        return GeomAbsShape::C0;
    };
    let Some(ed) = the_e.as_edge() else {
        return GeomAbsShape::C0;
    };
    for cr in &ed.representations {
        // OCCT: cr->IsRegularity(S1, S2, l1, l2) — the rcad face locations
        // are identity in this pipeline (arch. diff.: compose_pcurve_location
        // reduced to 0).
        if cr.is_regularity_on(&s1, &s2, 0, 0) {
            // OCCT: return cr->Continuity();
            if let CurveRepresentation::CurveOn2Surfaces { continuity, .. } = cr {
                return *continuity;
            }
        }
    }
    // OCCT: return GeomAbs_C0;
    GeomAbsShape::C0
}

/// OCCT BRep_Tool::IsClosed(theShape) (BRep_Tool.cxx L1707-1771) — the
/// shell/wire/edge branches (the default branch returns theShape.Closed(),
/// read here from the TShape flags).
fn brep_tool_is_closed(the_shape: &Shape) -> bool {
    if the_shape.shape_type() == ShapeType::Shell {
        let mut a_map: HashMap<(u64, u32), ()> = HashMap::new();
        let mut has_bound = false;
        let mut oriented = the_shape.clone();
        oriented.orientation = Orientation::Forward;
        for e in explorer(&oriented, ShapeType::Edge, ShapeType::Shape) {
            if brep_tool_degenerated(&e)
                || e.orientation == Orientation::Internal
                || e.orientation == Orientation::External
            {
                continue;
            }
            has_bound = true;
            // OCCT: if (!aMap.Add(E)) { aMap.Remove(E); } — the parity toggle.
            if a_map.remove(&shape_key(&e)).is_none() {
                a_map.insert(shape_key(&e), ());
            }
        }
        has_bound && a_map.is_empty()
    } else if the_shape.shape_type() == ShapeType::Wire {
        let mut a_map: HashMap<(u64, u32), ()> = HashMap::new();
        let mut has_bound = false;
        let mut oriented = the_shape.clone();
        oriented.orientation = Orientation::Forward;
        for v in explorer(&oriented, ShapeType::Vertex, ShapeType::Shape) {
            if v.orientation == Orientation::Internal || v.orientation == Orientation::External {
                continue;
            }
            has_bound = true;
            if a_map.remove(&shape_key(&v)).is_none() {
                a_map.insert(shape_key(&v), ());
            }
        }
        has_bound && a_map.is_empty()
    } else if the_shape.shape_type() == ShapeType::Edge {
        let (a_v_first, a_v_last) = top_exp_vertices_shape(the_shape);
        !a_v_first.is_null() && a_v_first.is_same(&a_v_last)
    } else {
        // OCCT: return theShape.Closed() — the TShape CLOSED flag.
        tshape_is_closed(the_shape.data.as_ref())
    }
}

/// OCCT TopoDS_Shape::Closed() — the TShape CLOSED flag read across the rcad
/// TShape variants (the flag is a member of every OCCT TopoDS_TShape).
fn tshape_is_closed(the_tshape: &TShape) -> bool {
    match the_tshape {
        TShape::Vertex(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Edge(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Wire(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Face(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Shell(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Solid(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::CompSolid(_) => false,
        TShape::Compound(_) => false,
    }
}

/// OCCT TopExp::Vertices(E, VFirst, VLast) — the null Shape form of the
/// loc_ope_wires_on_shape re-host.
pub(crate) fn top_exp_vertices_shape(edg: &Shape) -> (Shape, Shape) {
    let (v1, v2) = top_exp_vertices(edg);
    (
        v1.unwrap_or_else(Shape::null),
        v2.unwrap_or_else(Shape::null),
    )
}

/// OCCT BRepAdaptor_Surface(F, false)::D1(U, V, P, D1U, D1V) — the rcad
/// stand-in evaluates the face surface directly (architecture difference #6):
/// the point through SurfaceEval::point_at, the derivatives through the
/// GeomAdaptor_Surface::DN vehicle.
fn brep_adaptor_surface_d1(the_f: &Shape, the_u: f64, the_v: f64) -> (DVec3, DVec3, DVec3) {
    let s = brep_tool_surface(the_f).expect("BRepAdaptor_Surface: the face carries no surface");
    let p = s.point_at(the_u, the_v);
    let d1u = s.dn(the_u, the_v, 1, 0);
    let d1v = s.dn(the_u, the_v, 0, 1);
    (p, d1u, d1v)
}

/// OCCT Geom_Surface::Bounds(U1, U2, V1, V2) — the rcad stand-in through the
/// SurfaceEval::default_domain vehicle.
fn surface_bounds(the_s: &Surface3) -> (f64, f64, f64, f64) {
    let [u1, u2, v1, v2] = the_s.default_domain();
    (u1, u2, v1, v2)
}

/// OCCT Geom_Surface::UIso(U) — GAP: the iso-curve construction has no rcad
/// translation yet (the kernel iso sampling is private to base/convert).
fn surface_u_iso(the_s: &Surface3, the_u: f64) -> Curve3 {
    let _ = (the_s, the_u);
    panic!("GAP: Geom_Surface::UIso (iso-curve construction not translated)");
}

// ---------------------------------------------------------------------------
// OCCT statics (BRepOffset_MakeSimpleOffset.cxx).
// ---------------------------------------------------------------------------

//=============================================================================
// function : tgtfaces
// purpose  : check the angle at the border between two squares.
//           Two shares should have a shared front edge.
//=============================================================================
// OCCT BRepOffset_MakeSimpleOffset.cxx L215-322.
fn tgtfaces(ed: &Shape, f1: &Shape, f2: &Shape, couture: bool, the_res_angle: &mut f64) {
    // Check that pcurves exist on both faces of edge.
    // OCCT L224-233: aCurve = BRep_Tool::CurveOnSurface(Ed, F1/F2, aFirst,
    // aLast); the handles are consumed by the null checks only.
    if brep_tool_curve_on_surface(ed, f1).is_none() {
        return;
    }
    if brep_tool_curve_on_surface(ed, f2).is_none() {
        return;
    }

    let mut e = ed.clone();
    // OCCT L237-249: BRepAdaptor_Surface aBAS1/aBAS2, HS1/HS2 (reduced to the
    // face handles — architecture difference #6).
    let hs1 = f1.clone();
    let hs2 = if couture { hs1.clone() } else { f2.clone() };
    // case when edge lies on the one face

    e.orientation = Orientation::Forward;
    // OCCT L253: BRepAdaptor_Curve2d C2d1(E, F1).
    let c2d1 = brep_tool_curve_on_surface(&e, f1).map(|(c, _, _)| c);
    if couture {
        e.orientation = Orientation::Reversed;
    }
    // OCCT L258: BRepAdaptor_Curve2d C2d2(E, F2).
    let c2d2 = brep_tool_curve_on_surface(&e, f2).map(|(c, _, _)| c);
    let (Some(c2d1), Some(c2d2)) = (c2d1, c2d2) else {
        // The OCCT adaptor cannot fail here (the pcurves were checked above);
        // the rcad Option form is unpacked defensively.
        return;
    };

    let rev1 = f1.orientation == Orientation::Reversed;
    let rev2 = f2.orientation == Orientation::Reversed;
    let (mut f, mut l) = brep_tool_range(&e);
    // OCCT L264: Extrema_LocateExtPC ext; — declared but never used in the
    // OCCT body; kept as the unused binding.
    let _ext: Option<u8> = None;

    let eps = (l - f) / 100.0;
    f += eps; // to avoid calculations on
    l -= eps; // points of pointed squares.

    const NBPNT: i32 = 23;
    for i in 0..=NBPNT {
        // First suppose that this is sameParameter
        let u = f + (l - f) * i as f64 / NBPNT as f64;

        // take derivatives of surfaces at the same u, and compute normals
        let p = c2d1.point_at(u);
        let (_pp1, du1, dv1) = brep_adaptor_surface_d1(&hs1, p.x, p.y);
        let mut d1 = du1.cross(dv1);
        let norm = d1.length();
        if norm > 1.0e-12 {
            d1 /= norm;
        } else {
            continue; // skip degenerated point
        }
        if rev1 {
            d1 = -d1;
        }

        let p = c2d2.point_at(u);
        let (_pp2, du2, dv2) = brep_adaptor_surface_d1(&hs2, p.x, p.y);
        let mut d2 = du2.cross(dv2);
        let norm = d2.length();
        if norm > 1.0e-12 {
            d2 /= norm;
        } else {
            continue; // skip degenerated point
        }
        if rev2 {
            d2 = -d2;
        }

        // Compute angle.
        let a_current_ang = d1.angle_between(d2);

        *the_res_angle = the_res_angle.max(a_current_ang);
    }
}

//=============================================================================
// function : ComputeMaxAngleOnShape
// purpose  : Code the regularities on all edges of the shape, boundary of
//            two faces that do not have it.
//=============================================================================
// OCCT BRepOffset_MakeSimpleOffset.cxx L329-388.
fn compute_max_angle_on_shape(s: &Shape, the_res_angle: &mut f64) {
    // OCCT L331-333: NCollection_IndexedDataMap<TopoDS_Shape,
    // NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> M;
    // TopExp::MapShapesAndAncestors(S, TopAbs_EDGE, TopAbs_FACE, M);
    let mut m: indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)> = indexmap::IndexMap::new();
    map_shapes_and_ancestors(s, ShapeType::Edge, ShapeType::Face, &mut m);
    for (_key, (e, ancestors)) in m.iter() {
        let e = e.clone();
        let mut found = false;
        let mut couture = false;
        let mut f1 = Shape::null();
        let mut f2 = Shape::null();
        // OCCT L344-358: for (It.Initialize(M.FindFromIndex(i));
        // It.More() && !found; It.Next()) { ... }
        for it in ancestors {
            if found {
                break;
            }
            if f1.is_null() {
                f1 = it.clone();
            } else if !f1.is_same(it) {
                found = true;
                f2 = it.clone();
            }
        }
        if !found && !f1.is_null() {
            // is it a sewing edge?
            let or_e = e.orientation;
            // OCCT L363-373: for (Ex.Init(F1, TopAbs_EDGE); Ex.More() &&
            // !found; Ex.Next()) { ... }
            for cur_e in explorer(&f1, ShapeType::Edge, ShapeType::Shape) {
                if found {
                    break;
                }
                if e.is_same(&cur_e) && or_e != cur_e.orientation {
                    found = true;
                    couture = true;
                    f2 = f1.clone();
                }
            }
        }
        if found {
            // OCCT L376: if (BRep_Tool::Continuity(E, F1, F2) <= GeomAbs_C0).
            if brep_tool_continuity(&e, &f1, &f2) <= GeomAbsShape::C0 {
                // OCCT L378-384: the Standard_Failure catch is a no-op; the
                // rcad tgtfaces does not panic on the same inputs.
                tgtfaces(&e, &f1, &f2, couture, the_res_angle);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT class (BRepOffset_MakeSimpleOffset.hxx L61-172, cxx L46-703).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_MakeSimpleOffset (BRepOffset_MakeSimpleOffset.hxx L61-172).
pub struct BRepOffsetMakeSimpleOffset {
    // Input data.
    my_input_shape: Shape,           // OCCT: myInputShape (hxx L134)
    my_offset_value: f64,            // OCCT: myOffsetValue (hxx L137)
    my_tolerance: f64,               // OCCT: myTolerance (hxx L140)
    my_is_build_solid: bool,         // OCCT: myIsBuildSolid (hxx L143)

    // Internal data.
    my_max_angle: f64,               // OCCT: myMaxAngle (hxx L148)
    my_error: BRepOffsetSimpleStatus, // OCCT: myError (hxx L151)
    my_is_done: bool,                // OCCT: myIsDone (hxx L154)
    my_map_ve: HashMap<ShapeKey, Shape>, // OCCT: myMapVE (hxx L158)
    my_builder: BRepToolsModifier,   // OCCT: myBuilder (hxx L161)
    my_re_shape: ShapeBuildReShape,  // OCCT: myReShape (hxx L164)
    my_brep: BRep,                   // rcad arena stand-in (arch. diff. #4)

    // Output data.
    my_res_shape: Shape,             // OCCT: myResShape (hxx L169)
}

impl Default for BRepOffsetMakeSimpleOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetMakeSimpleOffset {
    /// OCCT BRepOffset_MakeSimpleOffset::BRepOffset_MakeSimpleOffset()
    /// (cxx L48-57).
    pub fn new() -> Self {
        BRepOffsetMakeSimpleOffset {
            my_input_shape: Shape::null(),
            my_offset_value: 0.0,
            my_tolerance: rcad_kernel::core::precision::CONFUSION,
            my_is_build_solid: false,
            my_max_angle: 0.0,
            my_error: BRepOffsetSimpleStatus::Ok,
            my_is_done: false,
            my_map_ve: HashMap::new(),
            my_builder: BRepToolsModifier::new(),
            my_re_shape: ShapeBuildReShape::new(),
            my_brep: BRep::new(),
            my_res_shape: Shape::null(),
        }
    }

    /// OCCT BRepOffset_MakeSimpleOffset::BRepOffset_MakeSimpleOffset(
    /// theInputShape, theOffsetValue) (cxx L61-72) — Rust has no overloading.
    pub fn with_input(the_input_shape: &Shape, the_offset_value: f64) -> Self {
        let mut res = Self::new();
        res.my_input_shape = the_input_shape.clone();
        res.my_offset_value = the_offset_value;
        res
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Initialize (cxx L76-82) — Initialise
    /// shape for modifications.
    pub fn initialize(&mut self, the_input_shape: &Shape, the_offset_value: f64) {
        self.my_input_shape = the_input_shape.clone();
        self.my_offset_value = the_offset_value;
        self.clear();
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetErrorMessage (cxx L86-117) —
    /// Gets error message.
    pub fn get_error_message(&self) -> &'static str {
        if self.my_error == BRepOffsetSimpleStatus::NullInputShape {
            return "Null input shape";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorOffsetComputation {
            return "Error during offset construction";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorWallFaceComputation {
            return "Error during building wall face";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorInvalidNbShells {
            return "Result contains two or more shells";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorNonClosedShell {
            return "Result shell is not closed";
        }
        ""
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetError() (hxx L81).
    pub fn get_error(&self) -> BRepOffsetSimpleStatus {
        self.my_error
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetBuildSolidFlag() (hxx L85).
    pub fn get_build_solid_flag(&self) -> bool {
        self.my_is_build_solid
    }

    /// OCCT BRepOffset_MakeSimpleOffset::SetBuildSolidFlag(theBuildFlag)
    /// (hxx L88).
    pub fn set_build_solid_flag(&mut self, the_build_flag: bool) {
        self.my_is_build_solid = the_build_flag;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetOffsetValue() (hxx L91).
    pub fn get_offset_value(&self) -> f64 {
        self.my_offset_value
    }

    /// OCCT BRepOffset_MakeSimpleOffset::SetOffsetValue(theOffsetValue)
    /// (hxx L94).
    pub fn set_offset_value(&mut self, the_offset_value: f64) {
        self.my_offset_value = the_offset_value;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetTolerance() (hxx L97).
    pub fn get_tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT BRepOffset_MakeSimpleOffset::SetTolerance(theValue) (hxx L100).
    pub fn set_tolerance(&mut self, the_value: f64) {
        self.my_tolerance = the_value;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::IsDone() (hxx L103).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetResultShape() (hxx L106).
    pub fn get_result_shape(&self) -> Shape {
        self.my_res_shape.clone()
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Clear (cxx L121-128) — Clears
    /// previous result.
    fn clear(&mut self) {
        self.my_is_done = false;
        self.my_error = BRepOffsetSimpleStatus::Ok;
        self.my_max_angle = 0.0;
        self.my_map_ve.clear();
        self.my_re_shape.clear(); // Clear possible stored modifications.
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetSafeOffset (cxx L132-151) —
    /// Computes max safe offset value for the given tolerance.
    pub fn get_safe_offset(&mut self, the_expected_toler: f64) -> f64 {
        if self.my_input_shape.is_null() {
            return 0.0; // Input shape is null.
        }

        // Compute max angle in faces junctions.
        if self.my_max_angle == 0.0 {
            // Non-initialized.
            self.compute_max_angle();
        }

        // OCCT L146: aMaxTol = BRep_Tool::MaxTolerance(myInputShape,
        // TopAbs_VERTEX).
        let a_max_tol = brep_tool_max_tolerance(&self.my_input_shape, ShapeType::Vertex);

        // OCCT L148-149: std::max((theExpectedToler - aMaxTol) /
        // (2.0 * myMaxAngle), 0.0) — Minimal distance can't be lower than 0.0.
        let an_exp_offset = ((the_expected_toler - a_max_tol) / (2.0 * self.my_max_angle)).max(0.0);
        an_exp_offset
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Perform (cxx L155-208) — Computes
    /// offset shape.
    pub fn perform(&mut self) {
        // Clear result of previous computations.
        self.clear();

        // Check shape existence.
        if self.my_input_shape.is_null() {
            self.my_error = BRepOffsetSimpleStatus::NullInputShape;
            return;
        }

        if self.my_max_angle == 0.0 {
            // Non-initialized.
            self.compute_max_angle();
        }

        self.my_builder.init(&self.my_input_shape);
        let a_mapper = BRepOffsetSimpleOffset::new(
            &self.my_input_shape,
            self.my_offset_value,
            self.my_tolerance,
        );
        self.my_builder.perform(&a_mapper);

        if !self.my_builder.is_done() {
            self.my_error = BRepOffsetSimpleStatus::ErrorOffsetComputation;
            return;
        }

        self.my_res_shape = self.my_builder.modified_shape(&self.my_input_shape);

        // Fix degeneracy. Degenerated edge should be mapped to the degenerated.
        // OCCT L186-199: BRep_Builder aBB; the rcad edit goes through the
        // modifier arena (arch. diff. #4).
        let mut a_bb = BRepBuilder::new();
        let an_exp_se = explorer(&self.my_input_shape, ShapeType::Edge, ShapeType::Shape);
        for a_curr_edge in &an_exp_se {
            if !brep_tool_degenerated(a_curr_edge) {
                continue;
            }

            let an_edge = self.my_builder.modified_shape(a_curr_edge);
            // OCCT L198: aBB.Degenerated(anEdge, true).
            a_bb.set_edge_degenerated(self.my_builder.brep_pool_mut(), an_edge, true);
        }

        // Restore walls for solid.
        if self.my_is_build_solid && !self.build_missing_walls() {
            return;
        }

        self.my_is_done = true;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::ComputeMaxAngle (cxx L394-397) —
    /// Computes max angle in faces junction.
    fn compute_max_angle(&mut self) {
        let mut res_angle = self.my_max_angle;
        compute_max_angle_on_shape(&self.my_input_shape, &mut res_angle);
        self.my_max_angle = res_angle;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::BuildMissingWalls (cxx L403-509) —
    /// Builds walls to the result solid.
    fn build_missing_walls(&mut self) -> bool {
        // Internal list of new faces.
        let mut a_bb = BRepBuilder::new();
        let a_new_faces = a_bb.make_compound(&mut self.my_brep, Vec::new());

        // Compute outer bounds of original shape.
        // OCCT L411-412: ShapeAnalysis_FreeBounds aFB(myInputShape);
        // GetClosedWires.
        let a_fb = ShapeAnalysisFreeBounds::new(&self.my_input_shape);
        let a_free_wires = a_fb.get_closed_wires();

        // Build linear faces on each edge and its image.
        let an_exp_cw = explorer(&a_free_wires, ShapeType::Wire, ShapeType::Shape);
        for a_cur_wire in &an_exp_cw {
            // Iterate over outer edges in outer wires.
            let an_exp_we = explorer(a_cur_wire, ShapeType::Edge, ShapeType::Shape);
            for a_cur_edge in &an_exp_we {
                let a_new_face = self.build_wall_face(a_cur_edge);

                if a_new_face.is_null() {
                    self.my_error = BRepOffsetSimpleStatus::ErrorWallFaceComputation;
                    return false;
                }

                // OCCT L434: aBB.Add(aNewFaces, aNewFace).
                a_bb.add_to_compound(&mut self.my_brep, a_new_faces.clone(), a_new_face);
            }
        }

        // Update edges from wall faces.
        // OCCT L439-447: ShapeFix_Edge aSFE; aSFE.SetContext(myReShape);
        // aSFE.FixSameParameter(aCurrEdge).
        for a_curr_edge in explorer(&a_new_faces, ShapeType::Edge, ShapeType::Shape) {
            // Fix same parameter and same range flags.
            shape_fix_edge_fix_same_parameter(&mut self.my_re_shape, &a_curr_edge);
        }

        // Update result to be compound.
        // OCCT L450-451: TopoDS_Compound aResCompound; aBB.MakeCompound.
        let a_res_compound = a_bb.make_compound(&mut self.my_brep, Vec::new());

        // Add old faces the result.
        for a_f in explorer(&self.my_input_shape, ShapeType::Face, ShapeType::Shape) {
            a_bb.add_to_compound(&mut self.my_brep, a_res_compound.clone(), a_f);
        }

        // Add new faces the result.
        for a_f in explorer(&self.my_res_shape, ShapeType::Face, ShapeType::Shape) {
            a_bb.add_to_compound(&mut self.my_brep, a_res_compound.clone(), a_f);
        }

        // Add wall faces to the result.
        for a_f in explorer(&a_new_faces, ShapeType::Face, ShapeType::Shape) {
            a_bb.add_to_compound(&mut self.my_brep, a_res_compound.clone(), a_f);
        }

        // Apply stored modifications.
        // OCCT L476: aResCompound = TopoDS::Compound(myReShape->Apply(...)).
        let a_res_compound = self
            .my_re_shape
            .apply(&mut self.my_brep, &a_res_compound, ShapeType::Shape);

        // Create result shell.
        // OCCT L479-481: BRepTools_Quilt aQuilt; aQuilt.Add(aResCompound);
        // aShells = aQuilt.Shells().
        let mut a_quilt = BRepToolsQuilt;
        a_quilt.add(&a_res_compound);
        let a_shells = a_quilt.shells();

        let mut a_res_shell = Shape::null();
        for a_shell in explorer(&a_shells, ShapeType::Shell, ShapeType::Shape) {
            if !a_res_shell.is_null() {
                // Shell is not null -> explorer contains two or more shells.
                self.my_error = BRepOffsetSimpleStatus::ErrorInvalidNbShells;
                return false;
            }
            a_res_shell = a_shell;
        }

        // OCCT L496: if (!BRep_Tool::IsClosed(aResShell)).
        if !brep_tool_is_closed(&a_res_shell) {
            self.my_error = BRepOffsetSimpleStatus::ErrorNonClosedShell;
            return false;
        }

        // Create result solid.
        // OCCT L503-506: aBB.MakeSolid(aResSolid); aBB.Add(aResSolid,
        // aResShell); myResShape = aResSolid.
        let a_res_solid = a_bb.make_solid(&mut self.my_brep, vec![a_res_shell]);
        self.my_res_shape = a_res_solid;

        true
    }

    /// OCCT BRepOffset_MakeSimpleOffset::BuildWallFace (cxx L513-667) —
    /// Builds face on specified wall.
    fn build_wall_face(&mut self, the_orig_edge: &Shape) -> Shape {
        let a_res_face = Shape::null();

        // Get offset edge. offset edge is reversed to create correct wire.
        // OCCT L518: aNewEdge = TopoDS::Edge(myBuilder.ModifiedShape(...)).
        let mut a_new_edge = self.my_builder.modified_shape(the_orig_edge);
        a_new_edge.orientation = Orientation::Reversed;

        // OCCT L522: TopExp::Vertices(aNewEdge, aNewV1, aNewV2).
        let (a_new_v1, a_new_v2) = top_exp_vertices_shape(&a_new_edge);

        // Wire contour is:
        // theOrigEdge (forcible forward) -> wall1 -> aNewEdge (forcible reversed) -> wall2
        // Firstly it is necessary to create copy of original shape with forward direction.
        // This simplifies walls creation.
        let mut an_orig_copy = the_orig_edge.clone();
        an_orig_copy.orientation = Orientation::Forward;
        // OCCT L531: TopExp::Vertices(anOrigCopy, aV1, aV2).
        let (a_v1, a_v2) = top_exp_vertices_shape(&an_orig_copy);

        // To simplify work with map.
        // OCCT L534-535: TopoDS::Vertex(aV1.Oriented(TopAbs_FORWARD)).
        let a_forward_v1 = oriented_vertex(&a_v1, Orientation::Forward);
        let a_forward_v2 = oriented_vertex(&a_v2, Orientation::Forward);

        // Check existence of edges in stored map: Edge1
        let a_wall1 = if self.my_map_ve.contains_key(&shape_key(&a_forward_v2)) {
            // Edge exists - get it from map.
            self.my_map_ve[&shape_key(&a_forward_v2)].clone()
        } else {
            // Edge does not exist - create it and add to the map.
            // OCCT L547-552: BRepLib_MakeEdge aME1(...); if (!aME1.IsDone())
            // return aResFace; — the rcad BRepBuilder::add_edge cannot fail
            // (architecture difference #5).
            let mut a_bb = BRepBuilder::new();
            let a_me1 = a_bb.add_edge(
                &mut self.my_brep,
                None,
                oriented_vertex(&a_v2, Orientation::Forward),
                oriented_vertex(&a_new_v2, Orientation::Reversed),
                [0.0, 0.0],
            );
            let a_wall1 = a_me1;

            self.my_map_ve
                .insert(shape_key(&a_forward_v2), a_wall1.clone());
            a_wall1
        };

        // Check existence of edges in stored map: Edge2
        let a_wall2 = if self.my_map_ve.contains_key(&shape_key(&a_forward_v1)) {
            // Edge exists - get it from map.
            // OCCT L563: TopoDS::Edge(myMapVE(aForwardV1).Oriented(
            // TopAbs_REVERSED)).
            let mut w = self.my_map_ve[&shape_key(&a_forward_v1)].clone();
            w.orientation = Orientation::Reversed;
            w
        } else {
            // Edge does not exist - create it and add to the map.
            // OCCT L568-580: BRepLib_MakeEdge aME2(...) — the IsDone guard has
            // no rcad counterpart (architecture difference #5).
            let mut a_bb = BRepBuilder::new();
            let a_me2 = a_bb.add_edge(
                &mut self.my_brep,
                None,
                oriented_vertex(&a_v1, Orientation::Forward),
                oriented_vertex(&a_new_v1, Orientation::Reversed),
                [0.0, 0.0],
            );
            let mut a_wall2 = a_me2;

            self.my_map_ve
                .insert(shape_key(&a_forward_v1), a_wall2.clone());

            // Orient it in reversed direction.
            a_wall2.orientation = Orientation::Reversed;
            a_wall2
        };

        let mut a_bb = BRepBuilder::new();

        // OCCT L584-589: aBB.MakeWire(aWire); aBB.Add(aWire, ...) x4.
        let a_wire = a_bb.make_wire(&mut self.my_brep);
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), an_orig_copy.clone());
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), a_wall1.clone());
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), a_new_edge.clone());
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), a_wall2.clone());

        // Build 3d curves on wire
        // OCCT L592: BRepLib::BuildCurves3d(aWire).
        brep_lib_build_curves3d(&a_wire);

        // Try to build using simple planar approach.
        // OCCT L595-607: the face maker is wrapped by try/catch since it
        // generates exceptions sometimes; the rcad BRepLibMakeFace carrier
        // keeps IsDone() = false (architecture difference #3).
        let mut a_f = Shape::null();
        {
            // Call of face maker is wrapped by try/catch since it generates exceptions sometimes.
            let a_fm = BRepLibMakeFace::from_wire(&a_wire, true);
            if a_fm.is_done() {
                a_f = a_fm.face();
            }
        }

        if a_f.is_null() {
            // Exception in face maker or result is not computed.
            // Build using thrusections.
            // OCCT L612: bool ToReverse = false.
            let to_reverse = false;
            // OCCT L614-615: EdgeCurve = BRep_Tool::Curve(theOrigEdge, fpar,
            // lpar); TrEdgeCurve = new Geom_TrimmedCurve(EdgeCurve, fpar,
            // lpar).
            let (edge_curve, fpar, lpar) = match edge_curve_of(the_orig_edge) {
                Some(v) => v,
                None => return a_res_face,
            };
            let tr_edge_curve = TrimmedCurve3::new(edge_curve.clone(), fpar, lpar);
            // OCCT L616-618: OffsetCurve = BRep_Tool::Curve(aNewEdge, fparOE,
            // lparOE); TrOffsetCurve = new Geom_TrimmedCurve(...).
            let (offset_curve, fpar_oe, lpar_oe) = match edge_curve_of(&a_new_edge) {
                Some(v) => v,
                None => return a_res_face,
            };
            let tr_offset_curve = TrimmedCurve3::new(offset_curve.clone(), fpar_oe, lpar_oe);

            // OCCT L620-624: GeomFill_Generator ThrusecGenerator;
            // AddCurve x2; Perform(Precision::PConfusion()); theSurf =
            // Surface().
            let mut thrusec_generator = GeomFillGenerator;
            thrusec_generator.add_curve(&tr_edge_curve);
            thrusec_generator.add_curve(&tr_offset_curve);
            thrusec_generator.perform(rcad_kernel::core::precision::PCONFUSION);
            let the_surf = thrusec_generator.surface();
            // OCCT L626-627: theSurf->Bounds(Uf, Ul, Vf, Vl).
            let (uf, ul, vf, vl) = surface_bounds(&the_surf);
            // OCCT L628: TopLoc_Location Loc; — the rcad location index
            // (0 = identity).
            let loc: u32 = 0;
            // OCCT L629-633: Geom2d_Line pcurves bound at (0, Vf)/(0, Vl) in
            // the X direction.  GAP: the BRep_Builder::UpdateEdge(E, C2d, S,
            // L, Tol) surface-keyed pcurve storage has no rcad carrier — the
            // calls are kept as no-op re-hosts (annotated).
            let edge_line2d = line2d_x(0.0, vf);
            update_edge_pcurve_on_surface(&the_orig_edge, &edge_line2d, &the_surf, loc);
            let oe_line2d = line2d_x(0.0, vl);
            update_edge_pcurve_on_surface(&a_new_edge, &oe_line2d, &the_surf, loc);
            // OCCT L634-637.
            let u_on_v1 = if to_reverse { ul } else { uf };
            let u_on_v2 = if to_reverse { uf } else { ul };
            let a_line2d = line2d_y(u_on_v2, 0.0);
            let a_line2d2 = line2d_y(u_on_v1, 0.0);
            if a_wall1.is_same(&a_wall2) {
                // OCCT L640: aBB.UpdateEdge(aWall1, aLine2d, aLine2d2,
                // theSurf, Loc, Precision::Confusion()).
                update_edge_pcurves_on_surface(&a_wall1, &a_line2d, &a_line2d2, &the_surf, loc);
                // OCCT L641-642: BSplC34 = theSurf->UIso(Uf);
                // aBB.UpdateEdge(aWall1, BSplC34, Precision::Confusion()).
                let bspl_c34 = surface_u_iso(&the_surf, uf);
                update_edge_curve3d_gap(&a_wall1, &bspl_c34);
                // OCCT L643: aBB.Range(aWall1, Vf, Vl).
                a_bb.set_edge_range(&mut self.my_brep, a_wall1.clone(), vf, vl);
            } else {
                // OCCT L647-650: aBB.SameParameter/SameRange(aWall1/aWall2,
                // false).
                a_bb.set_edge_same_parameter(&mut self.my_brep, a_wall1.clone(), false);
                a_bb.set_edge_same_range(&mut self.my_brep, a_wall1.clone(), false);
                a_bb.set_edge_same_parameter(&mut self.my_brep, a_wall2.clone(), false);
                a_bb.set_edge_same_range(&mut self.my_brep, a_wall2.clone(), false);
                // OCCT L651-654.
                update_edge_pcurve_on_surface(&a_wall1, &a_line2d, &the_surf, loc);
                set_edge_range_on_surface_gap(&a_wall1, &the_surf, loc, vf, vl);
                update_edge_pcurve_on_surface(&a_wall2, &a_line2d2, &the_surf, loc);
                set_edge_range_on_surface_gap(&a_wall2, &the_surf, loc, vf, vl);
                // OCCT L655-660: BSplC3 = theSurf->UIso(UonV2);
                // aBB.UpdateEdge(aWall1, BSplC3, ...); aBB.Range(aWall1, Vf,
                // Vl, true); BSplC4 = theSurf->UIso(UonV1);
                // aBB.UpdateEdge(aWall2, BSplC4, ...); aBB.Range(aWall2, Vf,
                // Vl, true).
                let bspl_c3 = surface_u_iso(&the_surf, u_on_v2);
                update_edge_curve3d_gap(&a_wall1, &bspl_c3);
                set_edge_range3d_gap(&a_wall1, vf, vl); // only for 3d curve
                let bspl_c4 = surface_u_iso(&the_surf, u_on_v1);
                update_edge_curve3d_gap(&a_wall2, &bspl_c4);
                set_edge_range3d_gap(&a_wall2, vf, vl); // only for 3d curve
            }

            // OCCT L663: aF = BRepLib_MakeFace(theSurf, aWire) — the
            // rcad BRepBuilder::make_face is the same storage-only stand-in
            // for the BRepLib_MakeFace(S, W) constructor.
            a_f = a_bb.make_face(&mut self.my_brep, Some(the_surf), a_wire);
        }

        a_f
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Generated (cxx L671-686) — Returns
    /// result shape for the given one (if exists).  The OCCT const method
    /// takes &mut self (architecture difference #7: the rcad ReShape Apply).
    pub fn generated(&mut self, the_shape: &Shape) -> Shape {
        // Shape generated by modification.
        let mut a_res = self.my_builder.modified_shape(the_shape);

        if a_res.is_null() {
            return a_res;
        }

        // Shape modifications obtained in scope of shape healing.
        // OCCT L683: aRes = myReShape->Apply(aRes) (the until/until shape
        // argument defaults to TopAbs_SHAPE).
        a_res = self
            .my_re_shape
            .apply(&mut self.my_brep, &a_res, ShapeType::Shape);

        a_res
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Modified (cxx L690-703) — Returns
    /// modified shape for the given one (if exists).  The OCCT const method
    /// takes &mut self (architecture difference #7: the rcad ReShape Status).
    pub fn modified(&mut self, the_shape: &Shape) -> Shape {
        let an_empty_shape = Shape::null();

        // Get modification status and new shape.
        // OCCT L695: int aModStatus = myReShape->Status(theShape, aRes).
        let (a_mod_status, a_res) = self
            .my_re_shape
            .status(&mut self.my_brep, the_shape, false);

        if a_mod_status == 0 {
            return an_empty_shape; // No modifications are applied to the shape or its sub-shapes.
        }

        a_res
    }
}

// ---------------------------------------------------------------------------
// Small local helpers (OCCT expression stand-ins).
// ---------------------------------------------------------------------------

/// OCCT aV.Oriented(TopAbs_FORWARD) — TopoDS::Vertex cast included.
pub(crate) fn oriented_vertex(the_v: &Shape, the_orient: Orientation) -> Shape {
    let mut s = the_v.clone();
    s.orientation = the_orient;
    s
}

/// OCCT BRep_Tool::Curve(E, f, l) — the (curve, fpar, lpar) triple (the
/// loc_ope_wires_on_shape_b re-host in the Option form).
pub(crate) fn edge_curve_of(the_e: &Shape) -> Option<(Curve3, f64, f64)> {
    match the_e.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT new Geom2d_Line(gp_Pnt2d(x, y), gp_Dir2d(gp_Dir2d::D::X)) — the
/// X-direction 2d line stand-in (Geom2d_Line -> Curve2d::Line(Line2d)).
fn line2d_x(x: f64, y: f64) -> Curve2d {
    Curve2d::Line(Line2d::new(DVec2::new(x, y), DVec2::new(1.0, 0.0)))
}

/// OCCT new Geom2d_Line(gp_Pnt2d(x, y), gp_Dir2d(gp_Dir2d::D::Y)) — the
/// Y-direction 2d line stand-in.
fn line2d_y(x: f64, y: f64) -> Curve2d {
    Curve2d::Line(Line2d::new(DVec2::new(x, y), DVec2::new(0.0, 1.0)))
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, S, L, Tol) — GAP: the rcad
/// BRepBuilder carries only face-keyed pcurve storage; the surface-keyed
/// form is a no-op re-host (annotated at the call site).
fn update_edge_pcurve_on_surface(the_e: &Shape, the_c2d: &Curve2d, the_s: &Surface3, the_l: u32) {
    let _ = (the_e, the_c2d, the_s, the_l);
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol) — the two-pcurve
/// surface-keyed no-op re-host.
fn update_edge_pcurves_on_surface(
    the_e: &Shape,
    the_c1: &Curve2d,
    the_c2: &Curve2d,
    the_s: &Surface3,
    the_l: u32,
) {
    let _ = (the_e, the_c1, the_c2, the_s, the_l);
}

/// OCCT BRep_Builder::UpdateEdge(E, C3d, Tol) — GAP no-op re-host (the rcad
/// 3d-curve representation is index-based, not value-based).
fn update_edge_curve3d_gap(the_e: &Shape, the_c: &Curve3) {
    let _ = (the_e, the_c);
}

/// OCCT BRep_Builder::Range(E, S, L, f, l) — GAP no-op re-host
/// (surface-keyed range storage).
fn set_edge_range_on_surface_gap(
    the_e: &Shape,
    the_s: &Surface3,
    the_l: u32,
    the_f: f64,
    the_last: f64,
) {
    let _ = (the_e, the_s, the_l, the_f, the_last);
}

/// OCCT BRep_Builder::Range(E, f, l, only3d = true) — GAP no-op re-host (the
/// rcad BRepBuilder has no only-3d range form).
fn set_edge_range3d_gap(the_e: &Shape, the_f: f64, the_l: f64) {
    let _ = (the_e, the_f, the_l);
}
