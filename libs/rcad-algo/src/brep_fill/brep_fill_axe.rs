//! OCCT BRepFill::Axe (TKBool/BRepFill) — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/BRepFill.cxx
//! L672-878 (the static `BRepFill::Axe` member; there is no dedicated
//! BRepFill_Axe.cxx in this OCCT tree — the function lives in BRepFill.cxx).
//!
//! First consumer: BRepOffsetAPI_MakeEvolved (Stage 2e; the L69 call
//! `BRepFill::Axe(Spine, Profil, Axis, POS, max(Tol, Precision::Confusion()))`).
//!
//! Reported gaps (plan §0.6, annotated at the call sites):
//! - BRepLib_FindSurface (TKTopAlgo) — only used when the spine face is not
//!   a Geom_Plane; GAP carrier below (same gap as CompatibleWires::PlaneOfWire);
//! - BRepLib_MakeFace(Wire, onlyPlane = true) — planar face synthesis from a
//!   wire; GAP carrier below (needs BRepLib_FindSurface internally);
//! - BRepExtrema_ExtPC (TKTopAlgo) — GAP carrier below (same family as the
//!   BRepExtrema_DistShapeShape gap of CompatibleWires).
//!
//! Architecture notes:
//! - `TopExp::MapShapesAndAncestors(F, VERTEX, EDGE, Map)` maps to the local
//!   `map_shapes_and_ancestors_ve` helper built on the wire-ordered edge
//!   lists of `TWireData` (same mapping as offset_wire::map_shapes_and_ancestors_ve).
//! - `BRep_Tool::Curve(E, L, f, l)` maps to the `TEdgeData` curve + range
//!   fields; OCCT's `TopLoc_Location L` composition is implicit because the
//!   rcad `TEdgeData` curve is stored in world space.

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval, Surface3};
use rcad_kernel::math::gp::Ax3;
use rcad_kernel::topo::topods::{BRep, Orientation, Shape, ShapeType, TShape};

use super::offset_wire_b::{brep_tool_surface, edge_vertices, explored_children};

/// OCCT Precision::Infinite().
const INFINITE: f64 = f64::INFINITY;

// ---------------------------------------------------------------------------
// GAP carriers (plan §0.6)
// ---------------------------------------------------------------------------

/// OCCT BRepLib_FindSurface (TKTopAlgo) — GAP: not translated (plan §0.6;
/// same gap as `CompatibleWires::PlaneOfWire`).  The OCCT constructor surface
/// is kept so the call sites stay 1:1.
pub struct BRepLibFindSurface;

impl BRepLibFindSurface {
    /// OCCT BRepLib_FindSurface(S, Tol = -1, OnlyPlane = false).
    pub fn new(_brep: &BRep, _s: &Shape, _tol: f64, _only_plane: bool) -> Self {
        panic!(
            "GAP: BRepLib_FindSurface (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT Found().
    pub fn found(&self) -> bool {
        panic!(
            "GAP: BRepLib_FindSurface (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT Surface().
    pub fn surface(&self) -> Surface3 {
        panic!(
            "GAP: BRepLib_FindSurface (TKTopAlgo) is not translated — see file header"
        )
    }
}

/// OCCT BRepLib_MakeFace(Wire, OnlyPlane) — GAP: the planar-face synthesis
/// needs BRepLib_FindSurface internally (plan §0.6).
pub struct BRepLibMakeFaceWire;

impl BRepLibMakeFaceWire {
    /// OCCT BRepLib_MakeFace(W, OnlyPlane).
    pub fn new(_brep: &mut BRep, _w: &Shape, _only_plane: bool) -> Self {
        panic!(
            "GAP: BRepLib_MakeFace(Wire, OnlyPlane) needs BRepLib_FindSurface \
             (TKTopAlgo, not translated) — see file header"
        )
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        panic!(
            "GAP: BRepLib_MakeFace(Wire, OnlyPlane) needs BRepLib_FindSurface \
             (TKTopAlgo, not translated) — see file header"
        )
    }

    /// OCCT Face().
    pub fn face(&self) -> Shape {
        panic!(
            "GAP: BRepLib_MakeFace(Wire, OnlyPlane) needs BRepLib_FindSurface \
             (TKTopAlgo, not translated) — see file header"
        )
    }
}

/// OCCT BRepExtrema_ExtPC (TKTopAlgo) — GAP: not translated (plan §0.6; same
/// family as the BRepExtrema_DistShapeShape gap of CompatibleWires).
pub struct BRepExtremaExtPC;

impl BRepExtremaExtPC {
    /// OCCT Initialize(E).
    pub fn initialize(&mut self, _brep: &BRep, _e: &Shape) {
        panic!(
            "GAP: BRepExtrema_ExtPC (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT Perform(V).
    pub fn perform(&mut self, _brep: &BRep, _v: &Shape) {
        panic!(
            "GAP: BRepExtrema_ExtPC (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        panic!(
            "GAP: BRepExtrema_ExtPC (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT NbExt().
    pub fn nb_ext(&self) -> usize {
        panic!(
            "GAP: BRepExtrema_ExtPC (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT IsMin(n).
    pub fn is_min(&self, _n: usize) -> bool {
        panic!(
            "GAP: BRepExtrema_ExtPC (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT SquareDistance(n).
    pub fn square_distance(&self, _n: usize) -> f64 {
        panic!(
            "GAP: BRepExtrema_ExtPC (TKTopAlgo) is not translated — see file header"
        )
    }

    /// OCCT Parameter(n).
    pub fn parameter(&self, _n: usize) -> f64 {
        panic!(
            "GAP: BRepExtrema_ExtPC (TKTopAlgo) is not translated — see file header"
        )
    }
}

// ---------------------------------------------------------------------------
// Kernel-mapping helper
// ---------------------------------------------------------------------------

/// OCCT TopExp::MapShapesAndAncestors(F, TopAbs_VERTEX, TopAbs_EDGE, Map) —
/// the vertex -> incident-edges map under `face`, in exploration order
/// (NCollection_IndexedDataMap mapping: Vec of pairs, index = position + 1).
fn map_shapes_and_ancestors_ve(brep: &BRep, face: &Shape) -> Vec<(Shape, Vec<Shape>)> {
    let mut map: Vec<(Shape, Vec<Shape>)> = Vec::new();
    let find_or_insert = |map: &mut Vec<(Shape, Vec<Shape>)>, v: &Shape| -> usize {
        if let Some(pos) = map.iter().position(|(k, _)| k.ptr_id() == v.ptr_id()) {
            pos
        } else {
            map.push((v.clone(), Vec::new()));
            map.len() - 1
        }
    };
    for w in explored_children(brep, face, ShapeType::Wire) {
        for e in explored_children(brep, &w, ShapeType::Edge) {
            let (vf, vl) = edge_vertices(brep, &e);
            let i1 = find_or_insert(&mut map, &vf);
            let i2 = find_or_insert(&mut map, &vl);
            map[i1].1.push(e.clone());
            if vf.ptr_id() != vl.ptr_id() {
                map[i2].1.push(e);
            }
        }
    }
    map
}

fn map_find<'a>(map: &'a [(Shape, Vec<Shape>)], v: &Shape) -> &'a Vec<Shape> {
    for (k, list) in map {
        if k.ptr_id() == v.ptr_id() {
            return list;
        }
    }
    panic!("BRepFill::Axe: MapShapesAndAncestors key not found");
}

/// OCCT BRep_Tool::Curve(E, L, f, l) — the 3d curve and its range.
fn brep_tool_curve(brep: &BRep, e: &Shape) -> Option<(Curve3, f64, f64)> {
    let ed = brep.edge(e.clone());
    ed.curve
        .clone()
        .map(|c| (c, ed.range[0], ed.range[1]))
}

// ---------------------------------------------------------------------------
// OCCT static void BRepFill::Axe(...) — BRepFill.cxx L672-878
// ---------------------------------------------------------------------------

/// OCCT static BRepFill::Axe(Spine, Profile, AxeProf, ProfOnSpine, Tol)
/// (BRepFill.cxx L672-878).
pub fn brep_fill_axe(
    brep: &mut BRep,
    spine: &Shape,
    profile: &Shape,
    axe_prof: &mut Ax3,
    prof_on_spine: &mut bool,
    tol: f64,
) {
    // OCCT L678-684: Loc, Loc1, Loc2, Tang, Tang1, Tang2, Normal, S, L, aFace.
    let mut loc: DVec3 = DVec3::ZERO;
    let mut loc1: DVec3 = DVec3::ZERO;
    let mut loc2: DVec3 = DVec3::ZERO;
    let mut tang: DVec3 = DVec3::ZERO;
    let mut tang1: DVec3 = DVec3::ZERO;
    let mut tang2: DVec3 = DVec3::ZERO;
    let mut normal: DVec3 = DVec3::ZERO;

    let mut s: Option<Surface3> = None;

    let mut a_face = Shape::null();

    // OCCT L686-709: normal to the Spine.
    if spine.shape_type() == ShapeType::Face {
        a_face = shape_as_type(spine, ShapeType::Face);
        s = brep_tool_surface(brep, &a_face);
        // OCCT L691: `if (!S->IsKind(STANDARD_TYPE(Geom_Plane)))` — the
        // TopLoc_Location handling (L690 `BRep_Tool::Surface(Face, L)`) is
        // implicit: the rcad TFaceData surface is stored in world space.
        if !matches!(s, Some(Surface3::Plane(_))) {
            // OCCT L693-698.
            let fs = BRepLibFindSurface::new(brep, spine, -1.0, true);
            if fs.found() {
                s = Some(fs.surface());
            } else {
                // OCCT L701.
                panic!("BRepFill_Evolved : The Face is not planar");
            }
        }
    } else if spine.shape_type() == ShapeType::Wire {
        // OCCT L707-708.
        let mf = BRepLibMakeFaceWire::new(brep, spine, true);
        a_face = mf.face();
        s = brep_tool_surface(brep, &a_face);
    }

    // OCCT L711-714.
    let s = match s {
        Some(sv) => sv,
        None => panic!("BRepFill_Evolved::Axe"),
    };

    // OCCT L716-719: `if (!L.IsIdentity()) S = S->Transformed(...)` — implicit
    // (world-space storage, see the architecture note above).

    // OCCT L721: Normal = down_cast<Geom_Plane>(S)->Pln().Axis().Direction().
    normal = match &s {
        Surface3::Plane(p) => p.normal,
        _ => panic!("BRepFill::Axe: the surface is not a plane"),
    };

    // OCCT L723-732: Find vertex of the profile closest to the spine.
    let mut dist_min = INFINITE;
    let mut dist: f64;
    // OCCT L727: double Tol2 = 1.e-10;
    let tol2 = 1.0e-10;
    let mut be = BRepExtremaExtPC;
    let mut par: f64 = 0.0;
    let mut p1: DVec3;
    let mut p2: DVec3;

    // OCCT L734-760: First check if there is contact Vertex Vertex.
    let mut is_on_vertex = false;
    let face_fwd = shape_oriented(&a_face, Orientation::Forward);
    let se_vertices = explored_children(brep, &face_fwd, ShapeType::Vertex);
    let mut se_index = 0usize; // TopExp_Explorer SE cursor
    for v_of_face in &se_vertices {
        se_index += 1;
        p1 = brep.vertex(v_of_face.clone()).point;

        let pe_vertices = explored_children(brep, profile, ShapeType::Vertex);
        for v_of_prof in &pe_vertices {
            p2 = brep.vertex(v_of_prof.clone()).point;
            let dist_p1p2 = p1.distance_squared(p2);
            is_on_vertex = dist_p1p2 <= tol2;
            if is_on_vertex {
                break;
            }
        }
        // otherwise SE.Next() is done and VonF is wrong
        if is_on_vertex {
            break;
        }
        //  modified by NIZHNY-EAP Wed Jan 26 09:08:36 2000 ___END___
    }

    if is_on_vertex {
        // OCCT L762-776: try to find on which edge which shared this vertex,
        // the profile must be considered.  E1, E2 : those two edges.
        let map = map_shapes_and_ancestors_ve(brep, &face_fwd);

        let von_f = &se_vertices[se_index - 1];
        let list = map_find(&map, von_f);
        let e1 = list.first().cloned().expect("BRepFill::Axe: empty ancestor list");
        let e2 = list.last().cloned().expect("BRepFill::Axe: empty ancestor list");

        // OCCT L778-789.
        let (ce1, _f1, _l1) = brep_tool_curve(brep, &e1).expect("BRepFill::Axe: no 3d curve on E1");
        let par1 = brep_tool_parameter_on_face(brep, von_f, &e1);
        loc1 = ce1.point_at(par1);
        tang1 = ce1.tangent_at(par1);
        if e1.orientation == Orientation::Reversed {
            tang1 = -tang1;
        }

        // OCCT L791-802.
        let (ce2, _f2, _l2) = brep_tool_curve(brep, &e2).expect("BRepFill::Axe: no 3d curve on E2");
        let par2 = brep_tool_parameter_on_face(brep, von_f, &e2);
        loc2 = ce2.point_at(par2);
        tang2 = ce2.tangent_at(par2);
        if e2.orientation == Orientation::Reversed {
            tang2 = -tang2;
        }

        //  modified by NIZHNY-EAP Wed Feb  2 15:38:41 2000 ___BEGIN___
        // OCCT L805-819.
        tang1 = tang1.normalize_or_zero();
        tang2 = tang2.normalize_or_zero();
        let mut sca1 = 0.0f64;
        let mut sca2 = 0.0f64;
        let pe_edges = explored_children(brep, profile, ShapeType::Edge);
        for e in &pe_edges {
            let (v1, v2) = edge_vertices(brep, e);
            let p1v = brep.vertex(v1).point;
            let p2v = brep.vertex(v2).point;
            let vec = p2v - p1v;
            sca1 += tang1.dot(vec).abs();
            sca2 += tang2.dot(vec).abs();
        }
        //  modified by NIZHNY-EAP Wed Feb  2 15:38:44 2000 ___END___

        // OCCT L822-833.
        if sca1.abs() < sca2.abs() {
            loc = loc1;
            tang = tang1;
        } else {
            loc = loc2;
            tang = tang2;
        }
        dist_min = 0.0;
    } else {
        // OCCT L834-872.
        let se_edges = explored_children(brep, &face_fwd, ShapeType::Edge);
        for e in &se_edges {
            be.initialize(brep, e);
            let pe_vertices = explored_children(brep, profile, ShapeType::Vertex);
            for v in &pe_vertices {
                dist = INFINITE;
                be.perform(brep, v);
                if be.is_done() {
                    // extrema.
                    for i in 1..=be.nb_ext() {
                        if be.is_min(i) {
                            dist = be.square_distance(i).sqrt();
                            par = be.parameter(i);
                            break;
                        }
                    }
                }
                // save minimum.
                if dist < dist_min {
                    dist_min = dist;
                    // OCCT L863-864: BRepAdaptor_Curve BAC(E); BAC.D1(Par, Loc, Tang)
                    // — mapped to the stored 3d curve evaluation.
                    let (c, _f, _l) =
                        brep_tool_curve(brep, e).expect("BRepFill::Axe: no 3d curve");
                    loc = c.point_at(par);
                    tang = c.tangent_at(par);
                    if e.orientation == Orientation::Reversed {
                        tang = -tang;
                    }
                }
            }
        }
    }

    // OCCT L874-877.
    *prof_on_spine = dist_min < tol;
    // Construction AxeProf;
    let a3 = Ax3::from_pnt_n_vx(loc, normal, tang);
    *axe_prof = a3;
}

// ---------------------------------------------------------------------------
// Small kernel-mapping helpers (shared form with generator.rs)
// ---------------------------------------------------------------------------

/// OCCT TopoDS::Shape cast by type — the wrapper re-typed (rcad `Shape` keeps
/// one wrapper type; the cast asserts the underlying TShape kind).
fn shape_as_type(s: &Shape, t: ShapeType) -> Shape {
    assert_eq!(s.shape_type(), t, "TopoDS cast mismatch in BRepFill::Axe");
    s.clone()
}

/// OCCT shape with an explicit orientation.
fn shape_oriented(s: &Shape, o: Orientation) -> Shape {
    let mut r = s.clone();
    r.orientation = o;
    r
}

/// OCCT BRep_Tool::Parameter(V, E, F) — the stored vertex parameter on the
/// edge (the face argument only selects the pcurve in OCCT; the parameter
/// value comes from the edge's vertex record, cf. generator.rs).
fn brep_tool_parameter_on_face(brep: &BRep, v: &Shape, e: &Shape) -> f64 {
    let ed = brep.edge(e.clone());
    ed.vertex_params
        .get(&v.ptr_id())
        .copied()
        .unwrap_or_else(|| panic!("BRepFill::Axe: no stored parameter for the vertex"))
}

/// Marker to keep the `TShape` import referenced by the architecture notes.
#[allow(dead_code)]
fn _tshape_marker(_: &TShape) {}
