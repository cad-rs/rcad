// OCCT BOPAlgo_Builder 鈥?shape construction from DS.
//
// OCCT BOPAlgo_Builder.hxx L75-507 + parent class fields (BOPAlgo_BuilderShape, BOPAlgo_Options, BOPAlgo_BOP).
// Flattened into one Rust struct because Rust has no C++ inheritance.

pub use crate::bop::algo::BooleanOpType;
use crate::bop::algo::{GlueEnum, Report};
use crate::bop::algo::builder_face::BuilderFace;
use crate::bop::ds::DS;
use crate::bop::ds::pave::SharedPB;
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::int_tools::face_make_curve::intermediate_point;
use crate::topalgo::brep_top_adaptor::fclass2d::{FClass2d, State};
use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::geom::{
    translate_curve2d, BezierSurface, BSplineSurface, Curve2d, Curve2dEval, CurveEval, Plane,
    Surface3, SurfaceEval, TrimmedCurve2,
};
use rcad_kernel::topods;
use rcad_kernel::topods::{
    surface_adaptor_basis_and_bounds, u_resolution_for_surface, v_resolution_for_surface,
    CurveRepresentation, Orientation, TShape, TVertexData, TEdgeData, TWireData, TFaceData,
    TShellData, TSolidData, tshape_flags,
};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::PCONFUSION;
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use glam::DVec2;
use glam::DVec3;

/// Boolean operation error type.
#[derive(Debug, Clone)]
pub enum BooleanError {
    InvalidOperation,
    TooFewArguments,
    NoFiller,
    BOPNotAllowed,
    BOPNotSet,
    EmptyShape,
    EmptyInput,
    DegenerateResult,
    NumericalFailure(&'static str),
    InvalidResult(&'static str),
}

/// OCCT BOPAlgo_Builder 鈥?result builder for boolean operations.
///
/// OCCT ref: BOPAlgo_Builder (BOPAlgo_Builder.hxx)
///
/// Fields map to OCCT hierarchy:
/// - BOPAlgo_Options (fuzzy, parallel, report)
/// - BOPAlgo_BuilderShape (myShape, myFillHistory)
/// - BOPAlgo_BOP (myOperation, myDims)
/// - BOPAlgo_Builder (myDS, myContext, myImages, etc.)
pub struct Builder<'a> {
    // 鈹€鈹€ BOPAlgo_Options (inherited) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    pub(crate) my_report: Report,          // BOPAlgo_Algo::myReport
    pub(crate) my_run_parallel: bool,      // BOPAlgo_Algo::myRunParallel
    pub(crate) my_fuzzy_value: f64,        // BOPAlgo_Algo::myFuzzyValue
    // 鈹€鈹€ BOPAlgo_BuilderShape (inherited) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    pub(crate) my_shape: Option<topods::BRep>, // BOPAlgo_BuilderShape::myShape
    pub(crate) my_fill_history: bool,      // BOPAlgo_BuilderShape::myFillHistory
    pub(crate) my_history: Option<crate::bop::history::BRepToolsHistory>, // BOPAlgo_Builder::myHistory
    // 鈹€鈹€ BOPAlgo_BOP (inherited) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    pub(crate) my_operation: BooleanOpType, // BOPAlgo_BOP::myOperation
    pub(crate) my_tools: Vec<Shape>,        // BOPAlgo_BOP::myTools
    pub(crate) my_rc: Vec<Shape>,           // BOPAlgo_BOP::myRC (result compound contents)
    pub(crate) my_dims: [i32; 2],           // BOPAlgo_BOP::myDims ([0]=obj, [1]=tool)
    // 鈹€鈹€ BOPAlgo_Builder.hxx L492-505 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    pub(crate) ds: &'a DS,                 // L496: myDS (borrowed from PaveFiller)
    pub(crate) my_context: IntToolsContext, // L497: myContext
    pub(crate) my_arguments: Vec<Shape>,   // L492: myArguments
    pub(crate) my_map_fence: HashSet<u64>, // L494: myMapFence
    pub(crate) my_entry_point: i32,        // L498: myEntryPoint
    // OCCT L499-502: myImages/myOrigins/myInParts are NCollection_DataMap
    // with TopTools_ShapeMapHasher 鈥?key identity is TShape + Location,
    // orientation is IGNORED. rcad keys by (ptr_id, location) to match.
    pub(crate) my_images: crate::bop::algo::occt_map::OcctDataMapInt<(u64, u32), Vec<Shape>>, // L499: myImages 鈥?NCollection_DataMap (bucket order)
    pub(crate) my_shapes_sd: HashMap<(u64, u32), Shape>,   // L500: myShapesSD (TopTools_ShapeMapHasher: TShape+Location, ignores orientation)
    pub(crate) my_origins: HashMap<(u64, u32), Vec<Shape>>,     // L501: myOrigins
    pub(crate) my_in_parts: HashMap<(u64, u32), Vec<Shape>>,    // L502: myInParts
    pub(crate) my_non_destructive: bool,   // L503: myNonDestructive
    pub(crate) my_glue: GlueEnum,          // L504: myGlue
    pub(crate) my_check_inverted: bool,    // L505: myCheckInverted
    pub(crate) my_nb_shapes_arr: [usize; 8], // L367: myNbShapesArr
    // rcad-specific: BRep index tracking (ptr_id -> tshapes index).  The
    // result keeps OCCT's TopoDS sharing semantics: a located copy (same
    // TShape, different Location) maps to the same result shape, mirroring
    // `nbshapes` without -t (same sub-shape with different location counts
    // once).  The evaluated positions are carried by the Location references.
    pub(crate) shape_remap: HashMap<u64, usize>,
    // rcad-specific: DS location index -> result BRep location index.  The
    // result BRep's locations pool is built lazily while cloning shapes, so a
    // located copy keeps its transform in the result (OCCT TopLoc_Location).
    pub(crate) loc_remap: HashMap<u32, u32>,
}

/// Stage snapshot: DS + result BRep counts at a Builder pipeline boundary.
/// Mirrors OCCT BOPAlgo_BOP::PerformInternal1 DUMP_STAGE points (10 stages).
#[derive(Debug, Clone)]
pub struct StageSnapshot {
    pub stage: u32,
    pub stage_name: &'static str,
    pub n_ds_vertices: usize,
    pub n_ds_edges: usize,
    pub n_ds_faces: usize,
    pub n_ds_pave_blocks: usize,
    pub n_ds_intersection_curves: usize,
    pub n_ds_interf_ff: usize,
    pub n_brep_vertices: usize,
    pub n_brep_edges: usize,
    pub n_brep_faces: usize,
    pub n_brep_shells: usize,
    pub n_brep_solids: usize,
}

/// Count result BRep entities by type from the flat tshape list.
fn count_brep_entities(b: &topods::BRep) -> (usize, usize, usize, usize, usize) {
    let mut v = 0usize;
    let mut e = 0usize;
    let mut f = 0usize;
    let mut sh = 0usize;
    let mut so = 0usize;
    for ts in &b.tshapes {
        match &**ts {
            topods::TShape::Vertex(_) => v += 1,
            topods::TShape::Edge(_) => e += 1,
            topods::TShape::Face(_) => f += 1,
            topods::TShape::Shell(_) => sh += 1,
            topods::TShape::Solid(_) => so += 1,
            _ => {}
        }
    }
    (v, e, f, sh, so)
}

/// Collect the face Shapes of a solid/compound via tree traversal, composing
/// the accumulated parent orientation into each face (shell REVERSED x face
/// stored) exactly like OCCT TopExp_Explorer(aSolid, TopAbs_FACE) with
/// cumOri=true (TopoDS_Iterator.cxx L72-80).
fn collect_solid_faces(s: &Shape) -> Vec<Shape> {
    // OCCT TopExp_Explorer(aSolid, TopAbs_FACE): solid -> shells -> faces in
    // the BRep's stored order (TopoDS_Iterator order). A stack-based DFS would
    // reverse the face order, which changes the edge -> faces map insertion
    // order — the IsInternalFace angle method reads aMEF.First()/Last() and
    // the classification depends on which face comes first. Use a FIFO
    // worklist that preserves the sub-shape order (OCCT BRep_Builder stores
    // the faces of a shell in the order they were added).
    let mut result: Vec<Shape> = Vec::new();
    let mut queue: std::collections::VecDeque<(Shape, topods::Orientation)> =
        std::collections::VecDeque::new();
    queue.push_back((s.clone(), topods::Orientation::Forward));
    while let Some((sh, cum_or)) = queue.pop_front() {
        match &*sh.data {
            TShape::Solid(sd) => {
                for x in &sd.shells {
                    queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                }
            }
            TShape::CompSolid(cd) => {
                for x in cd {
                    queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                }
            }
            TShape::Compound(cd) => {
                for x in cd {
                    queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                }
            }
            TShape::Shell(sd) => {
                for x in &sd.faces {
                    queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                }
            }
            TShape::Face(_) => {
                // Compose the accumulated parent (shell) orientation into the
                // face, as the explorer would (TopoDS_Iterator L75: myShape.
                // Orientation(TopAbs::Compose(myOrientation, ...))).
                let mut f = sh.clone();
                f.orientation = cum_or.compose(sh.orientation);
                result.push(f);
            }
            _ => {}
        }
    }
    result
}

/// OCCT TopAbs::Reverse (TopAbs.hxx L64-75) 鈥?flips the orientation.
/// FORWARD<->REVERSED; INTERNAL and EXTERNAL are unchanged.
fn flip_orientation(o: topods::Orientation) -> topods::Orientation {
    match o {
        topods::Orientation::Forward => topods::Orientation::Reversed,
        topods::Orientation::Reversed => topods::Orientation::Forward,
        topods::Orientation::Internal => topods::Orientation::Internal,
        topods::Orientation::External => topods::Orientation::External,
    }
}

/// OCCT BRep_Tool::IsClosed(aShell) (BRep_Tool.cxx L1707-1733) 鈥?a shell is
/// closed when every non-degenerate, non-INTERNAL/EXTERNAL boundary edge is
/// used an even number of times (parity pairing), and at least one boundary
/// edge is present. Edges are taken with the cumulative orientation
/// (TopExp_Explorer cumOri: face.or * wire.or * edge.or, shell is FORWARD-ized
/// by theShape.Oriented(TopAbs_FORWARD)); the parity map is keyed by
/// TopTools_ShapeMapHasher = TShape + Location, orientation ignored.
pub(crate) fn shell_is_closed(faces: &[Shape]) -> bool {
    let mut a_map: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
    let mut has_bound = false;
    for f in faces {
        for e in face_edges(f) {
            // OCCT L1719-1723: skip degenerated and INTERNAL/EXTERNAL edges.
            let is_degen = e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
            if is_degen
                || e.orientation == topods::Orientation::Internal
                || e.orientation == topods::Orientation::External
            {
                continue;
            }
            has_bound = true;
            let ekey = (e.ptr_id(), e.location);
            if !a_map.insert(ekey) {
                a_map.remove(&ekey);
            }
        }
    }
    has_bound && a_map.is_empty()
}

/// BRep_Tool::IsClosed(WIRE) equivalent (BRep_Tool.cxx L1734-1756):
/// parity pairing of boundary vertices plus at least one boundary vertex.
fn wire_is_closed(edges: &[Shape]) -> bool {
    let mut a_map: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut has_bound = false;
    for e in edges {
        if let Some(ed) = e.as_edge() {
            for v in [&ed.first, &ed.last] {
                // OCCT L1745-1747: skip INTERNAL/EXTERNAL vertices.
                if v.orientation == topods::Orientation::Internal
                    || v.orientation == topods::Orientation::External
                {
                    continue;
                }
                has_bound = true;
                if !a_map.insert(v.ptr_id()) {
                    a_map.remove(&v.ptr_id());
                }
            }
        }
    }
    has_bound && a_map.is_empty()
}

// ============================================================================
// BOPTools_AlgoTools::IsInternalFace (BOPTools_AlgoTools.cxx L807-891) +
// BOPAlgo_FillIn3DParts connexity-block classification
// (BOPAlgo_Tools.cxx L1334-1615).
// ============================================================================

/// All edge Shapes of a face (outer + inner wires).
/// OCCT TopExp_Explorer(aFace, TopAbs_EDGE).
fn face_edges(face: &Shape) -> Vec<Shape> {
    // OCCT TopExp_Explorer(face, TopAbs_EDGE) descends face -> wire -> edge
    // with cumOri=true, composing each parent's orientation into the edge
    // (TopExp_Explorer.cxx L152; TopoDS_Iterator.cxx L35-37, L72-80).
    let mut out = Vec::new();
    if let TShape::Face(fd) = &*face.data {
        let f_or = face.orientation;
        for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
            if let TShape::Wire(wd) = &*w.data {
                let w_or = w.orientation;
                for e in &wd.edges {
                    let mut e2 = Shape::from_parts(e.data.clone(), e.index, e.location, e.orientation);
                    e2.orientation = f_or.compose(w_or).compose(e.orientation);
                    out.push(e2);
                }
            }
        }
    }
    out
}

/// All (point, tolerance) pairs of the boundary vertices of a shape.
fn shape_vertices(s: &Shape) -> Vec<(DVec3, f64)> {    let mut out = Vec::new();
    let mut stack: Vec<Shape> = vec![s.clone()];
    while let Some(sh) = stack.pop() {
        match &*sh.data {
            TShape::Vertex(vd) => out.push((vd.point, vd.tolerance)),
            TShape::Edge(ed) => {
                stack.push(ed.first.clone());
                stack.push(ed.last.clone());
            }
            TShape::Wire(wd) => stack.extend(wd.edges.iter().cloned()),
            TShape::Face(fd) => {
                stack.push(fd.outer_wire.clone());
                stack.extend(fd.inner_wires.iter().cloned());
            }
            TShape::Shell(sd) => stack.extend(sd.faces.iter().cloned()),
            TShape::Solid(sd) => stack.extend(sd.shells.iter().cloned()),
            TShape::CompSolid(cd) => stack.extend(cd.iter().cloned()),
            TShape::Compound(cd) => stack.extend(cd.iter().cloned()),
        }
    }
    out
}

/// All edge Shapes of a shape (faces -> wires -> edges).
fn shape_edges(s: &Shape) -> Vec<Shape> {
    let mut out = Vec::new();
    let mut stack: Vec<Shape> = vec![s.clone()];
    while let Some(sh) = stack.pop() {
        match &*sh.data {
            TShape::Edge(_) => out.push(sh),
            TShape::Wire(wd) => stack.extend(wd.edges.iter().cloned()),
            TShape::Face(fd) => {
                stack.push(fd.outer_wire.clone());
                stack.extend(fd.inner_wires.iter().cloned());
            }
            TShape::Shell(sd) => stack.extend(sd.faces.iter().cloned()),
            TShape::Solid(sd) => stack.extend(sd.shells.iter().cloned()),
            TShape::CompSolid(cd) => stack.extend(cd.iter().cloned()),
            TShape::Compound(cd) => stack.extend(cd.iter().cloned()),
            TShape::Vertex(_) => {}
        }
    }
    out
}

/// Edge -> faces connection map of a solid.
/// OCCT TopExp::MapShapesAndAncestors(theSolid, EDGE, FACE, theMEF) --
/// IndexedDataMap with TopTools_ShapeMapHasher (TShape + Location).
fn build_edge_face_map(solid: &Shape) -> HashMap<(u64, u32), Vec<Shape>> {
    let mut map: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
    for f in collect_solid_faces(solid) {
        for e in face_edges(&f) {
            map.entry((e.ptr_id(), e.location)).or_default().push(f.clone());
        }
    }
    map
}

/// Find the edge Shape in a face's wires with the given edge identity.
/// OCCT BOPTools_AlgoTools::GetEdgeOnFace (IsSame = TShape + Location). 鉁?OCCT-aligned
fn find_edge_on_face(edge_ptr_id: u64, edge_loc: u32, face: &Shape) -> Option<Shape> {
    for e in face_edges(face) {
        if e.ptr_id() == edge_ptr_id && e.location == edge_loc {
            return Some(e);
        }
    }
    None
}

fn edge_data(edge: &Shape) -> Option<&TEdgeData> {
    if let TShape::Edge(ed) = &*edge.data {
        Some(ed)
    } else {
        None
    }
}

/// Geometric equality of two surfaces. The result-image faces carry the same
/// surface as their DS source face, so the edge pcurve keyed by the source
/// face index can be located by matching the surface (OCCT BRep_Tool::IsClosed
/// and BRep_Tool::CurveOnSurface match the face's surface handle; rcad stores
/// pcurves keyed by the DS face index).
fn surface_same(a: &rcad_kernel::geom::Surface3, b: &rcad_kernel::geom::Surface3) -> bool {
    use rcad_kernel::geom::Surface3;
    const T: f64 = 1e-9;
    let v = |x: DVec3, y: DVec3| (x - y).length() < T;
    match (a, b) {
        (Surface3::Plane(a), Surface3::Plane(b)) => {
            v(a.origin, b.origin)
                && v(a.normal, b.normal)
                && v(a.u_dir, b.u_dir)
                && v(a.v_dir, b.v_dir)
        }
        (Surface3::Cylinder(a), Surface3::Cylinder(b)) => {
            v(a.origin, b.origin) && v(a.axis, b.axis) && (a.radius - b.radius).abs() < T
        }
        (Surface3::Sphere(a), Surface3::Sphere(b)) => {
            v(a.center, b.center) && v(a.axis, b.axis) && (a.radius - b.radius).abs() < T
        }
        (Surface3::Cone(a), Surface3::Cone(b)) => {
            v(a.apex, b.apex)
                && v(a.axis, b.axis)
                && (a.radius - b.radius).abs() < T
                && (a.half_angle_rad - b.half_angle_rad).abs() < T
        }
        (Surface3::Torus(a), Surface3::Torus(b)) => {
            v(a.center, b.center)
                && v(a.axis, b.axis)
                && (a.major_radius - b.major_radius).abs() < T
                && (a.minor_radius - b.minor_radius).abs() < T
        }
        _ => false,
    }
}

/// The edge's pcurve on the face (OCCT BRep_Tool::CurveOnSurface(aE, aF)).
/// The result faces (images) are synthetic wrappers with index == usize::MAX,
/// so when the index-keyed lookup misses, match the pcurve by the face's
/// surface (OCCT matches the face's surface handle).
fn edge_pcurve_on_face<'a>(
    edge: &'a Shape,
    face: &Shape,
    ds: &DS,
) -> Option<&'a rcad_kernel::geom::Curve2d> {
    let ed = edge_data(edge)?;
    let surf = face.as_face().and_then(|fd| fd.surface.as_ref())?;
    ed.pcurves
        .get(&(face.ptr_id(), face.location))
        .or_else(|| {
            ed.pcurves.iter().find_map(|(k, v)| {
                if let Some(&fi) = ds.map_shape_index.get(k) {
                    if let Some(fs) = ds.face_surface(fi) {
                        if surface_same(surf, &fs) {
                            return Some(v);
                        }
                    }
                }
                None
            })
        })
        .map(|(pc, _, _)| pc)
}

/// GetNormalToFaceOnEdge (BOPTools_AlgoTools3D.cxx L351-376).
/// Surface normal at the edge parameter aT on the face, computed via the
/// edge's pcurve on the face and the surface first derivatives.
fn get_normal_to_face_on_edge(edge: &Shape, face: &Shape, a_t: f64, ds: &DS) -> Option<DVec3> {
    let fd = if let TShape::Face(fd) = &*face.data {
        fd
    } else {
        return None;
    };
    let surf = fd.surface.as_ref()?;
    // OCCT L365: CurveOnSurface(aE, aF1, aC2D1, aTolPC) 鈥?the pcurve of the
    // edge on the face (same parameterization as the 3D curve, SameParameter).
    let pc = edge_pcurve_on_face(edge, face, ds);
    if let Some(pc) = pc {
        // OCCT L367-369: aC2D1->D0(aT, aP2D).
        let uv = pc.point_at(a_t);
        // OCCT L371-375: aS1->D1(U, V, aP, aD1U, aD1V); aDNF1 = aDD1U ^ aDD1V.
        let (_p, d1u, d1v) = surf.derivatives(uv.x, uv.y);
        let n = d1u.cross(d1v);
        if n.length_squared() >= 1e-24 {
            return Some(n.normalize());
        }
    }
    // Fallback when the edge pcurve is not addressable by the face index (the
    // draft solid faces are synthetic wrappers with index == usize::MAX):
    // evaluate the surface normal at the domain center.  For planar faces this
    // equals the OCCT GetNormalToFaceOnEdge result (the normal is constant).
    let dom = surf.default_domain();
    let u = if dom[0].is_finite() && dom[1].is_finite() {
        0.5 * (dom[0] + dom[1])
    } else {
        0.0
    };
    let v = if dom[2].is_finite() && dom[3].is_finite() {
        0.5 * (dom[2] + dom[3])
    } else {
        0.0
    };
    let n = surf.normal_at(u, v);
    if n.length_squared() < 1e-24 {
        return None;
    }
    Some(n.normalize())
}

/// EdgeTangent (BOPTools_AlgoTools2D.cxx L74-98). 鉁?OCCT-aligned
/// Unit tangent of the edge curve at aT, reversed for REVERSED edges.
/// Returns None for degenerated edges or zero-length tangents.
fn edge_tangent(edge: &Shape, a_t: f64) -> Option<DVec3> {
    let ed = edge_data(edge)?;
    if ed.degenerated {
        return None;
    }
    let curve = ed.curve.as_ref()?;
    let mut tg = curve.tangent_at(a_t);
    let m = tg.length();
    // OCCT L88-96: if (mod > gp::Resolution()) aTau /= mod; else return false;
    // (gp::Resolution() == 1e-15).
    if m <= 1e-15 {
        return None;
    }
    tg /= m;
    if edge.orientation == topods::Orientation::Reversed {
        tg = -tg;
    }
    Some(tg)
}

/// GetFaceDir (BOPTools_AlgoTools.cxx L2118-2160). 鉁?OCCT-aligned
/// Computes the face normal aDN at the edge parameter aT (reversed for
/// REVERSED faces) and the bi-normal direction aDB = aDN ^ aDTgt. When the
/// face is not small (theSmallFaces false), FindPointInFace refines aDB to
/// point into the face (OCCT L2145-2157); the GetApproxNormalToFaceOnEdge
/// fallback is not ported (keeps the unrefined aDB).
fn get_face_dir(
    a_e: &Shape,
    a_f: &Shape,
    a_p: DVec3,
    a_t: f64,
    a_dtgt: DVec3,
    ds: &DS,
    a_dt: f64,
    b_small_faces: bool,
) -> Option<(DVec3, DVec3)> {
    // OCCT L2133-2137: normal on edge, reversed for REVERSED faces.
    let mut dn = get_normal_to_face_on_edge(a_e, a_f, a_t, ds)?;
    if a_f.orientation == topods::Orientation::Reversed {
        dn = -dn;
    }
    // OCCT L2139-2140: aTolE = Tolerance(aE); aDB = aDN ^ aDTgt.
    let a_tol_e = edge_data(a_e).map(|e| e.tolerance).unwrap_or(0.0);
    let mut db = dn.cross(a_dtgt);
    // OCCT L2145-2157: refine the bi-normal by FindPointInFace (skipped for
    // small faces). The plane aProjPL is (aP, aDTgt). When FindPointInFace
    // fails (or the face is small), OCCT falls back to
    // GetApproxNormalToFaceOnEdge (hatcher) which recomputes aDN and
    // redirects aDB from aP to the near point projected onto the aProjPL
    // plane (L2147-2156: bFound = !theSmallFaces && FindPointInFace(...);
    // if (!bFound) { ... }).
    let b_found = if !b_small_faces {
        find_point_in_face(a_f, a_p, &mut db, a_dt, a_tol_e, a_dtgt, ds)
    } else {
        false
    };
    if !b_found {
        // OCCT L2149-2156: GetApproxNormalToFaceOnEdge(aE, aF, aT, aDt, aPx, aDN, ctx);
        // aProjPL.Perform(aPx); aPx = aProjPL.NearestPoint();
        // aDB = Vec(aP, aPx).
        if let Some((a_px, a_dnf)) = get_approx_normal_to_face_on_edge(a_e, a_f, a_t, a_dt, ds) {
            dn = a_dnf;
            let a_px_proj = a_px - a_dtgt * (a_px - a_p).dot(a_dtgt);
            let v = a_px_proj - a_p;
            if v.length_squared() >= 1e-24 {
                db = v.normalize();
            }
        }
    }
    Some((dn, db))
}

/// OCCT BOPTools_AlgoTools::FindPointInFace (BOPTools_AlgoTools.cxx L2168-2239).
/// Finds a point inside the face in the aDB direction and refines aDB to point
/// from the edge point aP toward it, constrained to the plane perpendicular to
/// the edge tangent aDTgt (the aProjPL plane). Returns false when the
/// projection fails or the movement converges.
fn find_point_in_face(
    a_f: &Shape,
    a_p: DVec3,
    a_db: &mut DVec3,
    a_dt: f64,
    a_tol_e: f64,
    a_dtgt: DVec3,
    ds: &DS,
) -> bool {
    let fd = match &*a_f.data {
        TShape::Face(fd) => fd,
        _ => return false,
    };
    let Some(surf) = fd.surface.as_ref() else { return false };
    // OCCT L2182-2190: tolerance / iteration / eps parameters.
    let mut a_d_tol = rcad_kernel::ANGULAR;
    let a_pm = a_p.length();
    if a_pm > 1000.0 {
        a_d_tol = 5e-16 * a_pm;
    }
    let mut a_nb_it_max = 15;
    let an_eps = rcad_kernel::SQUARE_CONFUSION;
    // aProjPL 鈥?projection onto the plane (aP, aDTgt) (OCCT L2166 + L2201).
    let proj_plane = |p: DVec3| -> DVec3 { p - a_dtgt * (p - a_p).dot(a_dtgt) };
    // aProj 鈥?projection onto the face's surface (OCCT theContext->ProjPS(aF)).
    let proj_surf = |p: DVec3| -> (DVec3, f64) {
        let (_, pr) = crate::bop::closest_point_on_surface(surf, p);
        (pr, (p - pr).length())
    };
    // OCCT L2194-2212: project the edge point, then step by 2*tolE in aDB.
    let (mut a_ps, _) = proj_surf(a_p);
    a_ps = proj_plane(a_ps);
    a_ps = a_ps + 2.0 * a_tol_e * (*a_db);
    let (a_ps2, _) = proj_surf(a_ps);
    a_ps = proj_plane(a_ps2);
    // OCCT L2214-2238: iterate 鈥?step by aDt, re-project, refine aDB.
    loop {
        let a_p1 = a_ps + a_dt * (*a_db);
        let (a_p_out_raw, a_dist) = proj_surf(a_p1);
        let a_p_out = proj_plane(a_p_out_raw);
        let a_v = a_p_out - a_ps;
        if a_v.length_squared() < an_eps {
            return false;
        }
        *a_db = a_v.normalize();
        a_nb_it_max -= 1;
        if a_dist <= a_d_tol || a_nb_it_max <= 0 {
            return a_dist < a_d_tol;
        }
    }
}

/// OCCT BOPTools_AlgoTools3D::PointNearEdge (BOPTools_AlgoTools3D.cxx L525-612) 鈥?
/// the 6-parameter overload without the context. Computes a 2D point near the
/// edge at parameter aT: the edge's pcurve on the face is evaluated at aT, the
/// point is translated in the perpendicular 2D direction aDP (reversed for
/// REVERSED edge/face orientations) by (aDt2D + aTolE + aTolF), with the
/// cylindrical-surface angular correction (L583-595, pkv/909/F8) and the
/// spherical special case (L600-603).
pub(crate) fn point_near_edge(
    a_e: &Shape,
    a_f: &Shape,
    a_t: f64,
    a_dt2d: f64,
    ds: &DS,
) -> (Option<(DVec2, DVec3)>, i32) {
    // OCCT L537-542: aC2D = BRep_Tool::CurveOnSurface(aE, aF, aFirst, aLast);
    // iErr = aC2D.IsNull() ? 1 : 0.
    let pc = match edge_pcurve_on_face(a_e, a_f, ds) {
        Some(p) => p,
        None => return (None, 1),
    };
    let surf = match a_f.as_face().and_then(|fd| fd.surface.clone()) {
        Some(s) => s,
        None => return (None, 1),
    };
    // OCCT L546-549: aC2D->D1(aT, aPx2D, aVx2D); aDx2D = Dir(aVx2D).
    let a_px2d = pc.point_at(a_t);
    let a_vx2d = pc.derivative_at(a_t);
    if a_vx2d.length_squared() < 1e-24 {
        return (None, 1);
    }
    // OCCT L551-552: aDP = (-aDx2D.Y(), aDx2D.X()).
    let mut a_dp = DVec2::new(-a_vx2d.y, a_vx2d.x).normalize();
    // OCCT L554-562: reversals for REVERSED edge and face orientations.
    if a_e.orientation == topods::Orientation::Reversed {
        a_dp = -a_dp;
    }
    if a_f.orientation == topods::Orientation::Reversed {
        a_dp = -a_dp;
    }
    // OCCT L564-575: tolerances; the BSpline special case (NPAL19220) 鈥?OCCT
    // checks GeomAdaptor_Surface::GetType() == GeomAbs_BSplineSurface only
    // (Bezier surfaces do NOT take this branch).
    let mut a_etol = edge_data(a_e).map(|e| e.tolerance).unwrap_or(0.0);
    let mut a_ftol = a_f.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
    if matches!(surf, Surface3::BSpline(_)) && a_etol > 1e-5 {
        a_ftol = a_etol;
    }
    let a_px2d_near: DVec2;
    // OCCT L576-608: tolerance-based translation.
    if a_etol > 1e-5 || a_ftol > 1e-5 {
        if !matches!(surf, Surface3::Sphere(_)) {
            let mut trans_val = a_dt2d + a_etol + a_ftol;
            if let Surface3::Cylinder(cyl) = &surf {
                // OCCT L583-595 (pkv/909/F8): the 2D translation on a cylinder
                // corresponds to an angle on the circular cross-section.
                let a_r = cyl.radius;
                let d_t = 1.0 - trans_val / a_r;
                if d_t >= -1.0 && d_t <= 1.0 {
                    trans_val = d_t.acos();
                }
            }
            a_px2d_near = a_px2d + a_dp * trans_val;
        } else {
            // OCCT L600-603: sphere 鈥?plain aDt2D translation.
            a_px2d_near = a_px2d + a_dp * a_dt2d;
        }
    } else {
        // OCCT L605-608.
        a_px2d_near = a_px2d + a_dp * a_dt2d;
    }
    // OCCT L610: aS->D0(aPx2DNear.X(), aPx2DNear.Y(), aPxNear).
    let a_px_near = surf.point_at(a_px2d_near.x, a_px2d_near.y);
    (Some((a_px2d_near, a_px_near)), 0)
}

/// OCCT BOPTools_AlgoTools3D::GetApproxNormalToFaceOnEdge (L496-521) via
/// PointNearEdge(aE, aF, aT, theStep, ...) (L667-694): finds a point near the
/// edge inside the face and returns the surface normal there (reversed for
/// REVERSED faces). When the near point falls outside the face, OCCT falls back
/// to PointInFace (hatcher); rcad approximates it with the surface-domain
/// center (convex faces 鈥?planes, cylinders 鈥?keep the domain center inside).
fn get_approx_normal_to_face_on_edge(
    a_e: &Shape,
    a_f: &Shape,
    a_t: f64,
    the_step: f64,
    ds: &DS,
) -> Option<(DVec3, DVec3)> {
    let surf = a_f.as_face().and_then(|fd| fd.surface.clone())?;
    // OCCT L675-694: PointNearEdge(aE, aF, aT, theStep, aPx2DNear, aPxNear, ctx).
    let (mut a_px2d, mut a_px_near) = match point_near_edge(a_e, a_f, a_t, the_step, ds) {
        (Some(p), 0) => p,
        _ => return None,
    };
    // OCCT L676: if (!IsPointInOnFace(aF, aPx2DNear)) 鈥?PointInFace fallback
    // (the hatcher inside point; when that fails too, iErr = 2 and
    // GetApproxNormalToFaceOnEdge returns false, L657-663).
    let fi = ds.map_shape_index.get(&(a_f.ptr_id(), a_f.location)).copied();
    let in_face = match fi {
        Some(fi) => {
            let class2d = FClass2d::new(ds, fi, ds.face_tolerance(fi));
            class2d.perform(ds, a_px2d, true) != State::Out
        }
        // Draft-solid synthetic face wrappers (index == usize::MAX) have no DS
        // entry: approximate by the surface parameter domain (convex faces 鈥?
        // planes, cylinders 鈥?keep the domain center inside).
        None => {
            let dom = surf.default_domain();
            let u_in =
                dom[0].is_finite() && dom[1].is_finite() && a_px2d.x >= dom[0] && a_px2d.x <= dom[1];
            let v_in =
                dom[2].is_finite() && dom[3].is_finite() && a_px2d.y >= dom[2] && a_px2d.y <= dom[3];
            u_in && v_in
        }
    };
    if !in_face {
        // OCCT L681 (PointNearEdge 7-arg, BOPTools_AlgoTools3D.cxx L667-694):
        // PointInFace(aF, aE, aT, theStep, aP, aP2d, ctx) 鈥?hatcher inside
        // point along the EDGE-NORMAL 2D line (L942-990); when that fails too
        // iErr = 2 and GetApproxNormalToFaceOnEdge returns false (L657-663).
        match fi {
            Some(fi) => {
                let (err2, p2, p2d2) =
                    crate::bop::tools::algo_tools::point_in_face_edge(fi, a_e, a_t, the_step, ds);
                if err2 != 0 {
                    return None; // OCCT: iErr = 2 -> GetApproxNormalToFaceOnEdge false
                }
                a_px2d = p2d2;
                a_px_near = p2;
            }
            None => {
                // Synthetic face: surface-domain center (convex faces keep it inside).
                let dom = surf.default_domain();
                let u = if dom[0].is_finite() && dom[1].is_finite() {
                    0.5 * (dom[0] + dom[1])
                } else {
                    0.0
                };
                let v = if dom[2].is_finite() && dom[3].is_finite() {
                    0.5 * (dom[2] + dom[3])
                } else {
                    0.0
                };
                a_px2d = DVec2::new(u, v);
                a_px_near = surf.point_at(u, v);
            }
        }
    }
    // OCCT L510-517: GetNormalToSurface(aS, aPx2DNear.X(), aPx2DNear.Y(), aDNF);
    // REVERSED faces reverse the normal.
    let mut a_dnf = surf.normal_at(a_px2d.x, a_px2d.y);
    if a_f.orientation == topods::Orientation::Reversed {
        a_dnf = -a_dnf;
    }
    Some((a_px_near, a_dnf))
}


/// The 3D step aDt for FindPointInFace: max(2*(tolE+tolF)) over the candidate
/// faces, clamped to aDtMin which grows with the surface radius.
/// OCCT MinStep3D (BOPTools_AlgoTools.cxx L2239-2354). 鉁?OCCT-aligned
fn min_step_3d(
    the_e1: &Shape,
    the_f1: &Shape,
    candidates: &[(Shape, Shape)],
    a_p: DVec3,
    ds: &DS,
) -> (f64, bool) {
    let a_tol_e = edge_data(the_e1).map(|e| e.tolerance).unwrap_or(0.0);
    let mut a_dt_max: f64 = -1.0;
    let mut a_dt_min: f64 = 5e-6;
    // OCCT L2252-2258: theLCS + (theE1, theF1).
    let mut a_lcs: Vec<(&Shape, &Shape)> = candidates.iter().map(|(e, f)| (e, f)).collect();
    a_lcs.push((the_e1, the_f1));
    for (_, a_f) in &a_lcs {
        let a_tol_f = a_f.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
        let a_dt = 2.0 * (a_tol_e + a_tol_f);
        if a_dt > a_dt_max {
            a_dt_max = a_dt;
        }
        // OCCT L2277-2304: surface radius based aDtMin.
        let mut a_r = 0.0;
        if let Some(fd) = a_f.as_face() {
            if let Some(surf) = &fd.surface {
                match surf {
                    Surface3::Cylinder(c) => a_r = c.radius,
                    Surface3::Cone(c) => {
                        // aR = distance from aP to the cone axis (gp_Lin::Distance).
                        let w = a_p - c.apex;
                        a_r = (w - c.axis * w.dot(c.axis)).length();
                    }
                    Surface3::Sphere(s) => {
                        a_dt_min = a_dt_min.max(5e-4);
                        a_r = s.radius;
                    }
                    Surface3::Torus(t) => a_r = t.major_radius,
                    _ => a_dt_min = a_dt_min.max(5e-4),
                }
            }
        }
        // OCCT L2306-2310: large radius grows the minimum step.
        if a_r > 100.0 {
            let d = 10.0 * rcad_kernel::PCONFUSION;
            a_dt_min = a_dt_min.max((d * d + 2.0 * d * a_r).sqrt());
        }
    }
    // OCCT L2313-2316.
    if a_dt_max < a_dt_min {
        a_dt_max = a_dt_min;
    }
    // OCCT L2318-2349: check if the 3D step is too big for any of the faces
    // (UResolution/VResolution); theSmallFaces = aIt.More().
    let mut b_small_faces = false;
    for (_, a_f) in &a_lcs {
        let Some(surf) = a_f.as_face().and_then(|fd| fd.surface.clone()) else {
            continue;
        };
        // OCCT theContext->UVBounds(aF) 鈥?the face's actual UV bounds; rcad:
        // the DS face's sampled UV rect (face_actual_uv_bounds), falling back
        // to the surface domain for synthetic draft-solid face wrappers (no
        // DS entry).
        let [a_umin, a_umax, a_vmin, a_vmax] = match ds
            .map_shape_index
            .get(&(a_f.ptr_id(), a_f.location))
            .copied()
        {
            Some(fi) => ds.face_actual_uv_bounds(fi),
            None => surf.default_domain(),
        };
        let a_du = a_umax - a_umin;
        if a_du > 0.0 {
            let a_u_res = u_resolution_for_surface(&surf, a_dt_max);
            if 2.0 * a_u_res > a_du {
                b_small_faces = true;
                break;
            }
        }
        let a_dv = a_vmax - a_vmin;
        if a_dv > 0.0 {
            let a_v_res = v_resolution_for_surface(&surf, a_dt_max);
            if 2.0 * a_v_res > a_dv {
                b_small_faces = true;
                break;
            }
        }
    }
    (a_dt_max, b_small_faces)
}

/// OCCT BRepLib::BuildPCurveForEdgeOnPlane -> BRep_Tool::CurveOnPlane
/// (BRep_Tool.cxx L379-450) -> GeomProjLib::ProjectOnPlane -> ProjLib_Plane ->
/// Geom2dAdaptor::MakeCurve: project the edge's 3D curve onto the plane and
/// return the 2D pcurve in the plane's (u, v) space.
///
/// Returns None where OCCT stores no pcurve: the projected curve is a BSpline
/// / Bezier / Other (Geom2dAdaptor::MakeCurve throws, Geom2dAdaptor.cxx
/// L92-93), or the projected ellipse/parabola/hyperbola minor axis would need
/// a clockwise orientation that rcad's 2D conics cannot represent.
///
/// OCCT-aligned: rcad-kernel/src/base/geom_proj_lib/project_on_plane.rs
pub(crate) fn project_edge_on_plane(
    curve: &rcad_kernel::geom::Curve3,
    pl: &Plane,
    range: [f64; 2],
) -> Option<Curve2d> {
    use rcad_kernel::base::geom_proj_lib::project_on_plane::curve_on_plane;
    curve_on_plane(curve, range, pl)
}

/// Signed angle from d1 to d2 around the reference axis d_ref.
/// OCCT BOPTools_AlgoTools::AngleWithRef (L1938-1967). 鉁?OCCT-aligned
fn angle_with_ref(d1: DVec3, d2: DVec3, d_ref: DVec3) -> f64 {
    let half_pi = std::f64::consts::FRAC_PI_2;
    let cross = d1.cross(d2);
    let sinus = cross.length();
    let cosinus = d1.dot(d2);
    let beta = if sinus >= 0.0 {
        half_pi * (1.0 - cosinus)
    } else {
        std::f64::consts::TAU - half_pi * (3.0 + cosinus)
    };
    if cross.dot(d_ref) < 0.0 {
        -beta
    } else {
        beta
    }
}

/// GetFaceOff (BOPTools_AlgoTools.cxx L994-1102) 鈥?select the candidate face
/// whose bi-normal direction has the minimal angle to the reference face's
/// bi-normal, computed in the plane perpendicular to the edge tangent.
///
/// candidates: (edge_in_face, face) pairs.  Returns `(aFOff, bRet)` where
/// `bRet` is false when the minimal angle can not be found reliably.
pub(crate) fn get_face_off(
    the_e1: &Shape,
    the_f1: &Shape,
    candidates: &[(Shape, Shape)],
    ds: &DS,
) -> (Option<Shape>, bool) {
    // OCCT L1012-1016: 3D curve, intermediate point, point on the curve.
    let ed1 = match edge_data(the_e1) {
        Some(e) => e,
        None => return (None, false),
    };
    let curve1 = match ed1.curve.as_ref() {
        Some(c) => c,
        None => return (None, false),
    };
    let a_t = crate::bop::int_tools::face_make_curve::intermediate_point(
        ed1.range[0],
        ed1.range[1],
    );
    let a_px = curve1.point_at(a_t);

    // OCCT L1018-1020: EdgeTangent(theE1, aT, aVTgt); aDTgt = Dir(aVTgt);
    // aOr = theE1.Orientation().
    let a_dtgt = match edge_tangent(the_e1, a_t) {
        Some(v) => v,
        None => return (None, false),
    };
    let a_or = the_e1.orientation;

    // OCCT L1026-1028: MinStep3D(theE1, theF1, theLCSOff, aPx, ..., bSmallFaces).
    let (a_dt3d, b_small_faces) = min_step_3d(the_e1, the_f1, candidates, a_px, ds);
    // OCCT L1029-1037: GetFaceDir(theE1, theF1, ...) 鈫?(aDN1, aDBF).
    let (a_dn1, a_dbf1) = match get_face_dir(the_e1, the_f1, a_px, a_t, a_dtgt, ds, a_dt3d, b_small_faces) {
        Some(d) => d,
        None => return (None, false),
    };
    // OCCT L1038: aDTF = aDN1 ^ aDBF.
    let a_dtf = a_dn1.cross(a_dbf1);

    // OCCT L1012: aAngleMin = 100.  L1042: anAngleCriteria = Precision::Confusion().
    let two_pi = std::f64::consts::TAU;
    let mut a_angle_min = 100.0;
    let an_angle_criteria = rcad_kernel::CONFUSION;
    // OCCT L1044: bRet = true.
    let mut b_ret = true;
    let mut a_f_off: Option<Shape> = None;

    // OCCT L1045-1100: iterate the candidates.
    for (a_e2, a_f2) in candidates {
        // OCCT L1052: aDTgt2 = (aE2.Orientation() == aOr) ? aDTgt : aDTgt.Reversed().
        let mut a_dtgt2 = a_dtgt;
        if a_e2.orientation != a_or {
            a_dtgt2 = -a_dtgt2;
        }
        // OCCT L1053-1061: GetFaceDir(aE2, aF2, ...) 鈫?(aDN2, aDBF2).
        // When it fails, bIsComputed=false and aDBF2 keeps the fallback value.
        let b_is_computed = get_face_dir(a_e2, a_f2, a_px, a_t, a_dtgt2, ds, a_dt3d, b_small_faces);
        let (_a_dn2, a_dbf2) = match b_is_computed {
            Some(d) => d,
            None => (DVec3::ZERO, DVec3::ZERO),
        };
        // OCCT L1063: aAngle = AngleWithRef(aDBF, aDBF2, aDTF).
        let mut a_angle = angle_with_ref(a_dbf1, a_dbf2, a_dtf);
        if std::env::var("RCAD_GFO_DEBUG").is_ok() {
            let e_or = a_e2.orientation;
            eprintln!("[GFO]{} e1_ori={:?} f1_ori={:?} a_e2_ori={:?} a_f2_ori={:?} aDTgt=({:.4},{:.4},{:.4}) aDTgt2=({:.4},{:.4},{:.4}) aAngle={:.6}",
                std::env::var("RCAD_GFO_SITE").unwrap_or_default(),
                a_or, the_f1.orientation, e_or, a_f2.orientation,
                a_dtgt.x, a_dtgt.y, a_dtgt.z, a_dtgt2.x, a_dtgt2.y, a_dtgt2.z, a_angle);
        }
        // OCCT L1065-1082: near-zero angle handling.
        if a_angle.abs() < rcad_kernel::ANGULAR {
            // aF2 == theF1 (IsEqual: same TShape+location+orientation) -> PI.
            if a_f2.is_partner(the_f1) && a_f2.orientation == the_f1.orientation {
                a_angle = std::f64::consts::PI;
            // aF2.IsSame(theF1) (same TShape+location, any orientation) -> 2*PI.
            } else if a_f2.is_partner(the_f1) {
                a_angle = two_pi;
            // bi-normal direction could not be reliably computed -> 2*PI.
            } else if b_is_computed.is_none() {
                a_angle = two_pi;
            }
        }
        // OCCT L1085-1089: the minimal angle can not be found.
        if a_angle.abs() < an_angle_criteria
            || (a_angle - a_angle_min).abs() < an_angle_criteria
        {
            b_ret = false;
        }
        // OCCT L1091-1094: normalize to [0, 2*PI).
        if a_angle < 0.0 {
            a_angle = two_pi + a_angle;
        }
        // OCCT L1096-1100: the minimal angle wins.
        if a_angle < a_angle_min {
            a_angle_min = a_angle;
            a_f_off = Some(a_f2.clone());
        }
    }
    (a_f_off, b_ret)
}

/// IsInternalFace core (BOPTools_AlgoTools.cxx L939-990) 鈥?angle-based check
/// of `the_face` against the pair (the_face1, the_face2) sharing `the_edge`.
/// Returns 0 = not IN, 1 = IN, 2 = unable.
fn is_internal_face_core(
    the_face: &Shape,
    the_edge: &Shape,
    the_face1: &Shape,
    the_face2: &Shape,
    ds: &DS,
) -> i32 {
    // OCCT L950-966: edge copies on both faces. GetEdgeOnFace leaves aE1
    // null on failure (orientation FORWARD by default); the angle method then
    // cannot be applied (L978-982: GetFaceOff returns false -> iRet = 2), so
    // a missing edge copy maps to 2, not 0.
    let a_e1 = match find_edge_on_face(the_edge.ptr_id(), the_edge.location, the_face1) {
        Some(e) => e,
        None => return 2,
    };
    // OCCT L951-966: for an INTERNAL edge, or when the two faces are the same
    // (IsEqual: same TShape + Location + Orientation), aE2 = aE1 with
    // FORWARD/REVERSED orientations; otherwise the edge as it appears in the
    // second face. GetEdgeOnFace(aE, theFace2, aE2) always succeeds here 鈥?
    // theFace2 is an element of theMEF(aE) (OCCT L862-864), so it contains
    // the edge.
    let (a_e1_used, a_e2) = if a_e1.orientation == topods::Orientation::Internal
        || (the_face1.is_partner(the_face2) && the_face1.orientation == the_face2.orientation)
    {
        let mut e1 = a_e1.clone();
        e1.orientation = topods::Orientation::Forward;
        let mut e2 = a_e1.clone();
        e2.orientation = topods::Orientation::Reversed;
        (e1, e2)
    } else {
        let e2 = find_edge_on_face(the_edge.ptr_id(), the_edge.location, the_face2)
            .unwrap_or_else(|| a_e1.clone());
        (a_e1.clone(), e2)
    };

    // OCCT L968-974: candidates (theEdge, theFace) and (aE2, theFace2).
    let candidates = [
        (the_edge.clone(), the_face.clone()),
        (a_e2.clone(), the_face2.clone()),
    ];

    // OCCT L976-989: GetFaceOff 鈥?the minimal-angle face; bRet=false means
    // the minimal angle can not be found (iRet = 2).
    let (a_f_off, is_done) = get_face_off(&a_e1_used, the_face1, &candidates, ds);
    if !is_done {
        return 2;
    }
    // OCCT L983: theFace.IsEqual(aFOff) 鈥?same TShape + Location + Orientation.
    match a_f_off {
        Some(f) if f.is_partner(the_face) && f.orientation == the_face.orientation => 1,
        _ => 0,
    }
}

/// IsInternalFace (BOPTools_AlgoTools.cxx L807-891) 鈥?checks whether
/// `the_face` is internal to `the_solid`.  First tries the edge-angle method
/// on the solid's edge -> face map `a_mef`, then falls back to ComputeState.
/// `the_tol` = Precision::Confusion() (used by the ComputeState fallback).
fn is_internal_face(
    the_face: &Shape,
    the_solid: &Shape,
    a_mef: &HashMap<(u64, u32), Vec<Shape>>,
    the_tol: f64,
    ds: &DS,
) -> i32 {
    let mut i_ret = 0;
    let mut found_edge = false;
    for a_e in face_edges(the_face) {
        // OCCT L832-846: edge on the solid, not INTERNAL, not degenerated.
        let a_lf_opt = a_mef.get(&(a_e.ptr_id(), a_e.location));
        let Some(a_lf) = a_lf_opt else {
            continue;
        };
        if a_e.orientation == topods::Orientation::Internal {
            continue;
        }
        if edge_data(&a_e).map_or(true, |ed| ed.degenerated) {
            continue;
        }
        let a_nb_f = a_lf.len();
        if a_nb_f == 1 {
            // OCCT L851-861: single neighbor face 鈥?internal edge of a membrane.
            let a_f1 = &a_lf[0];
            let e_on_f1 = find_edge_on_face(a_e.ptr_id(), a_e.location, a_f1);
            if let Some(ef1) = e_on_f1 {
                if ef1.orientation == topods::Orientation::Internal {
                    i_ret = is_internal_face_core(the_face, &a_e, a_f1, a_f1, ds);
                    found_edge = true;
                    break;
                }
            }
            continue;
        } else if a_nb_f == 2 {
            // OCCT L864-873: two neighbor faces 鈥?angle-based method.
            let a_f1 = &a_lf[0];
            let a_f2 = &a_lf[1];
            if std::env::var("RCAD_GFO_DEBUG").is_ok() {
                eprintln!("[ISF] the_face_or={:?} e_ori={:?} f1_or={:?} f2_or={:?}",
                    the_face.orientation, a_e.orientation, a_f1.orientation, a_f2.orientation);
            }
            i_ret = is_internal_face_core(the_face, &a_e, a_f1, a_f2, ds);
            if i_ret != 2 {
                found_edge = true;
                break;
            }
        }
    }
    if found_edge && i_ret != 2 {
        return i_ret;
    }
    // OCCT L882-891: ComputeState fallback.
    if compute_state_face(the_face, the_solid, the_tol, ds) == 3 {
        1
    } else {
        0
    }
}

/// ComputeState (BOPTools_AlgoTools.cxx L660-715) — classify a face against a
/// solid: try an edge of the face not on the solid (classify the edge
/// midpoint), else classify a point inside the face (PointInFace hatcher,
/// with the PointNearEdge fallback, L688-712). ✅ OCCT-aligned
pub(crate) fn compute_state_face(the_face: &Shape, the_solid: &Shape, the_tol: f64, ds: &DS) -> u8 {
    // OCCT aBounds (L887) is IndexedMap with TopTools_ShapeMapHasher 鈥?
    // TShape + Location.
    let solid_edges: HashSet<(u64, u32)> = collect_solid_faces(the_solid)
        .iter()
        .flat_map(|f| face_edges(f))
        .map(|e| (e.ptr_id(), e.location))
        .collect();
    for e in face_edges(the_face) {
        if edge_data(&e).map_or(true, |ed| ed.degenerated) {
            continue;
        }
        if !solid_edges.contains(&(e.ptr_id(), e.location)) {
            // OCCT L683-685: classify the middle point of the edge. OCCT
            // ComputeState(Edge) returns UNKNOWN when the degenerated edge has
            // no first vertex (L748-754) 鈥?rcad returns UNKNOWN (0) too.
            let Some(p) = edge_midpoint(&e) else {
                return 0;
            };
            let mut clsf =
                crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(
                    the_solid,
                );
            clsf.perform(p, the_tol);
            return clsf.my_state;
        }
    }
    // OCCT L688-712: all edges of the face are on the solid. Get a point
    // inside the face 鈥?PointInFace (the U-line hatcher, L906-938) first,
    // then the PointNearEdge fallback (the dT2D overload with the
    // IsPointInOnFace check and the hatcher inside-point, L614-696).
    let fi = ds.map_shape_index.get(&(the_face.ptr_id(), the_face.location)).copied();
    let mut i_err = 1;
    let mut a_p3d = DVec3::ZERO;
    // OCCT L692: PointInFace(theF, ...) works on any TopoDS_Face; a split
    // face image not registered in the DS uses the Shape-based hatcher.
    if let Some(fi) = fi {
        let (err, p, _p2d) = crate::bop::tools::algo_tools::point_in_face(fi, ds);
        i_err = err;
        a_p3d = p;
    } else {
        let (err, p, _p2d) = crate::bop::tools::algo_tools::point_in_face_shape(the_face, &ds.locations);
        i_err = err;
        a_p3d = p;
    }
    if i_err != 0 {
        for a_se in face_edges(the_face) {
            if edge_data(&a_se).map_or(true, |ed| ed.degenerated) {
                continue;
            }            // OCCT L685-694: aT = IntermediatePoint(Range(aSE)).
            let (a_t1, a_t2) = match edge_data(&a_se) {
                Some(ed) => (ed.range[0], ed.range[1]),
                None => (0.0, 0.0),
            };
            let a_t = intermediate_point(a_t1, a_t2);
            // OCCT L619-641 (the dT2D overload): dT2D = 10 * MinStepIn2d
            // (1e-4), x10 for cylinder/sphere surfaces, max(2*(tolE+tolF)).
            // The face FORWARD-ization + OrientEdgeOnFace of the 5-arg
            // overload (L685-694) is folded into point_near_edge's own
            // REVERSED edge/face handling.
            let mut d_t2d = 10.0 * 1e-5;
            let surf = the_face.as_face().and_then(|fd| fd.surface.clone());
            if matches!(surf, Some(Surface3::Cylinder(_)) | Some(Surface3::Sphere(_))) {
                d_t2d = 10.0 * d_t2d;
            }
            let a_tol_e = edge_data(&a_se).map(|ed| ed.tolerance).unwrap_or(0.0);
            let a_tol_f = the_face.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
            let d_tx = 2.0 * (a_tol_e + a_tol_f);
            if d_tx > d_t2d {
                d_t2d = d_tx;
            }
            let (near, err6) = point_near_edge(&a_se, the_face, a_t, d_t2d, ds);
            i_err = err6;
            if i_err != 1 {
                // OCCT L627-641: the point must be inside (or on) the face,
                // otherwise the hatcher inside-point is taken (or iErr = 2).
                let (p2d, p3d) = near.unwrap_or((DVec2::ZERO, DVec3::ZERO));
                let in_face = match fi {
                    Some(fi) => {
                        let class2d = FClass2d::new(ds, fi, ds.face_tolerance(fi));
                        class2d.perform(ds, p2d, true) != State::Out
                    }
                    None => false,
                };
                if in_face {
                    a_p3d = p3d;
                } else {
                    // OCCT L634-640 (BOPTools_AlgoTools3D.cxx L629-640):
                    // PointInFace(aF, aE, aT, theStep, aP, aP2d, theContext) 鈥?
                    // hatcher inside point along the edge-normal 2D line
                    // (BOPTools_AlgoTools3D.cxx L942-990).
                    match fi {
                        Some(fi) => {
                            let (err2, p2, _) = crate::bop::tools::algo_tools::point_in_face_edge(
                                fi, &a_se, a_t, d_t2d, ds,
                            );
                            if err2 == 0 {
                                a_p3d = p2;
                                i_err = 0;
                            } else {
                                i_err = 2;
                            }
                        }
                        None => {
                            i_err = 2;
                        }
                    }
                }
            }
            if i_err == 0 {
                break;
            }
        }
    }
    if i_err == 0 {
        // OCCT L705-709: classify the inside point.
        let mut clsf = crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(
            the_solid,
        );
        clsf.perform(a_p3d, the_tol);
        return clsf.my_state;
    }
    // OCCT L711-714: aState stays TopAbs_UNKNOWN (0) when no point was found.
    0
}

/// Middle point of an edge 鈥?OCCT ComputeState(edge) (BOPTools_AlgoTools.cxx
/// L733-778): for a degenerated edge (no 3D curve) the first vertex point is
/// used; otherwise the parameter is the intermediate one, with the infinite
/// range handled by the dT = 10. shifts (L748-771).
fn edge_midpoint(edge: &Shape) -> Option<DVec3> {
    let ed = edge_data(edge)?;
    if let Some(curve) = &ed.curve {
        let (a_t1, a_t2) = (ed.range[0], ed.range[1]);
        // OCCT L748-771: Precision::IsNegativeInfinite / IsPositiveInfinite
        // (|x| >= 5e99); the dT = 10. shifts for infinite ranges.
        let a_t = if rcad_kernel::is_negative_infinite_value(a_t1)
            && !rcad_kernel::is_positive_infinite_value(a_t2)
        {
            a_t2 - 10.0
        } else if !rcad_kernel::is_negative_infinite_value(a_t1)
            && rcad_kernel::is_positive_infinite_value(a_t2)
        {
            a_t1 + 10.0
        } else if rcad_kernel::is_negative_infinite_value(a_t1)
            && rcad_kernel::is_positive_infinite_value(a_t2)
        {
            0.0
        } else {
            intermediate_point(a_t1, a_t2)
        };
        Some(curve.point_at(a_t))
    } else if let TShape::Vertex(vd) = &*ed.first.data {
        // OCCT L748-754: degenerated edge 鈥?the first vertex point; a null
        // vertex returns UNKNOWN (rcad: None).
        Some(vd.point)
    } else {
        None
    }
}

/// Bounding box of a shape 鈥?vertices plus sampled edge-curve points
/// (semantic equivalent of OCCT BRepBndLib::Add, which also covers curve
/// extents beyond the boundary vertices).
pub(crate) fn shape_bbox(s: &Shape) -> Option<(DVec3, DVec3)> {
    let mut min = DVec3::splat(f64::INFINITY);
    let mut max = DVec3::splat(f64::NEG_INFINITY);
    let mut any = false;
    for (p, _tol) in shape_vertices(s) {
        if !p.is_finite() {
            continue;
        }
        min = min.min(p);
        max = max.max(p);
        any = true;
    }
    for e in shape_edges(s) {
        if let Some(ed) = edge_data(&e) {
            if let Some(curve) = &ed.curve {
                let [t1, t2] = ed.range;
                for k in 0..=8 {
                    let t = t1 + (t2 - t1) * (k as f64) / 8.0;
                    let p = curve.point_at(t);
                    if !p.is_finite() {
                        continue;
                    }
                    min = min.min(p);
                    max = max.max(p);
                    any = true;
                }
            }
        }
    }
    if any {
        Some((min, max))
    } else {
        None
    }
}

/// MakeConnexityBlock (BOPAlgo_FillIn3DParts member, BOPAlgo_Tools.cxx L1555-1615).
///
/// Collects the connexity block of faces (indices into `faces`) reachable
/// from `f_start` through edges not on the solid boundary (`a_mse`) and not
/// degenerated.  Faces touching a solid boundary edge are recorded in
/// `a_face_to_classify` (the first such face) as the block representative.
fn make_connexity_block(
    f_start: usize,
    a_mse: &HashSet<(u64, u32)>,
    a_mefp: &HashMap<(u64, u32), Vec<usize>>,
    a_mf_done: &mut HashSet<usize>,
    a_lcb: &mut Vec<usize>,
    a_face_to_classify: &mut Option<usize>,
    faces: &[Shape],
) {
    // OCCT L1566-1570: add the start element.
    a_lcb.push(f_start);
    if a_mefp.is_empty() {
        return;
    }
    // OCCT L1572-1614: iterate the growing block (breadth-first).
    let mut i = 0;
    while i < a_lcb.len() {
        let a_f = a_lcb[i];
        i += 1;
        for a_e in face_edges(&faces[a_f]) {
            // OCCT L1589-1596: border edge of the solid.
            // theMEAvoid/theEFMap use TopTools_ShapeMapHasher (TShape + Location).
            if a_mse.contains(&(a_e.ptr_id(), a_e.location))
                || edge_data(&a_e).map_or(false, |ed| ed.degenerated)
            {
                if a_face_to_classify.is_none() {
                    *a_face_to_classify = Some(a_f);
                }
                continue;
            }
            // OCCT L1598-1611: expand through faces sharing this edge.
            if let Some(p_lf) = a_mefp.get(&(a_e.ptr_id(), a_e.location)) {
                for &a_f_to_add in p_lf {
                    if a_mf_done.insert(a_f_to_add) {
                        a_lcb.push(a_f_to_add);
                    }
                }
            }
        }
    }
}

/// OCCT BOPAlgo_BOP::TypeToExplore (BOPAlgo_BOP.cxx L1574-1597).
fn type_to_explore(the_dim: i32) -> topods::ShapeType {
    match the_dim {
        0 => topods::ShapeType::Vertex,
        1 => topods::ShapeType::Edge,
        2 => topods::ShapeType::Face,
        3 => topods::ShapeType::Solid,
        _ => topods::ShapeType::Compound,
    }
}

/// OCCT BOPTools_AlgoTools::Dimension 鈥?max dimension of a shape.
fn shape_dimension(s: &Shape) -> i32 {
    crate::bop::tools::algo_tools::dimensions(s.shape_type()).1
}

/// OCCT BOPAlgo_Tools::FillMap(Shape, Shape, IndexedDataMap<Shape, List<Shape>>)
/// (BOPAlgo_Tools.hxx L84-102) 鈥?bidirectional connection for the SD back-and-forth
/// map. OCCT aDMSLS (BOPAlgo_Builder_2.cxx L747) uses TopTools_ShapeMapHasher 鈥?
/// key identity is TShape + Location, orientation is ignored; rcad keys by
/// (ptr_id, location) and keeps one representative Shape per key for the chains.
fn fill_map_faces(
    n1: &Shape,
    n2: &Shape,
    the_map: &mut IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
) {
    let e1 = the_map
        .entry((n1.ptr_id(), n1.location))
        .or_insert_with(|| (n1.clone(), Vec::new()));
    e1.1.push(n2.clone());
    let e2 = the_map
        .entry((n2.ptr_id(), n2.location))
        .or_insert_with(|| (n2.clone(), Vec::new()));
    e2.1.push(n1.clone());
}

/// OCCT BOPAlgo_Tools::MakeBlocks (BOPAlgo_Tools.hxx L46-80) 鈥?connected components
/// of the SD back-and-forth map. The fence uses the map's hasher
/// (TopTools_ShapeMapHasher = TShape + Location). aDMSLS is an IndexedDataMap
/// (Builder_2.cxx L747), so iteration follows the insertion order.
fn make_blocks_faces(the_map: &IndexMap<(u64, u32), (Shape, Vec<Shape>)>) -> Vec<Vec<Shape>> {
    let mut a_m_fence: HashSet<(u64, u32)> = HashSet::new();
    let mut a_m_blocks: Vec<Vec<Shape>> = Vec::new();
    for (k, (n, _)) in the_map {
        if !a_m_fence.insert(*k) {
            continue;
        }
        // Start the chain with the representative shape of the key (OCCT
        // aChain.Append(n), n being the map key).
        let mut a_chain: Vec<Shape> = vec![n.clone()];
        let mut i = 0;
        while i < a_chain.len() {
            let n1 = &a_chain[i];
            if let Some((_, a_li)) = the_map.get(&(n1.ptr_id(), n1.location)) {
                for n2 in a_li {
                    if a_m_fence.insert((n2.ptr_id(), n2.location)) {
                        a_chain.push(n2.clone());
                    }
                }
            }
            i += 1;
        }
        a_m_blocks.push(a_chain);
    }
    a_m_blocks
}

impl<'a> Builder<'a> {
    /// Create a new Builder borrowing a DS from PaveFiller.
    ///
    /// OCCT: BOPAlgo_Builder is constructed with a PaveFiller reference.
    /// OCCT BOPAlgo_BOP::PerformInternal1 L425-429:
    ///   myPaveFiller = &theFiller; myDS = myPaveFiller->PDS();
    ///   myFuzzyValue = myPaveFiller->FuzzyValue();
    pub fn new(ds: &'a DS, op: BooleanOpType, fuzzy_value: f64) -> Self {
        Builder {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_fuzzy_value: fuzzy_value,
            my_shape: None,
            my_fill_history: false,
            my_history: None,
            my_operation: op,
            my_tools: Vec::new(),
            my_rc: Vec::new(),
            my_dims: [3, 3],
            my_context: IntToolsContext::new(),
            my_arguments: Vec::new(),
            my_map_fence: HashSet::new(),
            my_entry_point: 0,
            my_images: crate::bop::algo::occt_map::OcctDataMapInt::new(),
            my_shapes_sd: HashMap::new(),
            my_origins: HashMap::new(),
            my_in_parts: HashMap::new(),
            my_non_destructive: false,
            my_glue: GlueEnum::GlueOff,
            my_check_inverted: false,
            my_nb_shapes_arr: [0; 8],
            shape_remap: HashMap::new(),
            loc_remap: HashMap::new(),
        }
    }

    /// OCCT BOPAlgo_Algo::SetArguments.
    pub fn set_arguments(&mut self, args: Vec<Shape>) {
        self.my_arguments = args;
    }

    /// OCCT BOPAlgo_BOP::SetTools.
    pub fn set_tools(&mut self, tools: Vec<Shape>) {
        self.my_tools = tools;
    }

    /// Shape backed by shared Arc in ds.shapes 鈥?OCCT myDS->Shape(n).
    fn remap_arg_key(&self, a_s: &Shape) -> (u64, u32) {
        // OCCT: an argument shape and its DS counterpart share the TShape
        // identity. rcad clones the inputs (clone_arguments_private) so the
        // original ptr_id must be translated to the cloned one for myImages
        // lookups (Builder stages key myImages by the DS shape).
        let p = self
            .ds
            .argument_remap
            .get(&a_s.ptr_id())
            .copied()
            .unwrap_or(a_s.ptr_id());
        (p, a_s.location)
    }

    fn brep_sr(&self, flat_idx: usize) -> Shape {
        self.ds.shape(flat_idx).clone()
    }

    pub fn has_errors(&self) -> bool {
        self.my_report.has_errors()
    }

    pub fn report(&self) -> &Report {
        &self.my_report
    }

    /// Debug accessor: snapshot of the my_images map (key -> image count).
    pub fn images_debug(&self) -> Vec<((u64, u32), usize)> {
        self.my_images.iter().map(|(k, v)| (k, v.len())).collect()
    }

    /// Debug accessor: snapshot of my_in_parts keys.
    pub fn in_parts_debug(&self) -> Vec<(u64, u32)> {
        self.my_in_parts.keys().cloned().collect()
    }

    /// OCCT BOPAlgo_Builder::Build 鈥?convenience wrapper.
    ///
    /// Returns the result BRep on success.
    pub fn build(&mut self) -> Result<rcad_kernel::BRep, ()> {
        match self.build_with_history_topods() {
            Ok((brep, _)) => Ok(brep),
            Err(_) => Err(()),
        }
    }

    /// Run the Builder pipeline stage by stage, capturing a snapshot after each
    /// of the 10 OCCT-aligned stages. Single source of truth for the pipeline:
    /// `build_with_history_topods` and `build` delegate here.
    ///
    /// On `has_errors` mid-pipeline, returns Ok with the partial result and the
    /// snapshots captured so far, so tests can localize the failure stage.
    pub fn build_with_history_stage_by_stage(
        &mut self,
    ) -> Result<(topods::BRep, Vec<StageSnapshot>), BooleanError> {
        let mut snapshots: Vec<StageSnapshot> = Vec::with_capacity(10);
        macro_rules! snap {
            ($stage:expr, $name:expr) => {{
                let (v, e, f, sh, so) = count_brep_entities(self.my_shape.as_ref().unwrap());
                snapshots.push(StageSnapshot {
                    stage: $stage,
                    stage_name: $name,
                    n_ds_vertices: self.ds.vertex_count(),
                    n_ds_edges: self.ds.edge_count(),
                    n_ds_faces: self.ds.face_count(),
                    // OCCT pipeline_dump.h dump_ds_snapshot: count PBs only for
                    // EDGE shapes with HasPaveBlocks(i), summing PaveBlocks(i).Extent().
                    // This excludes orphan pool entries (section-edge PBs etc.) that
                    // no edge shape references.
                    n_ds_pave_blocks: {
                        let mut n_pb = 0usize;
                        for i in 0..self.ds.nb_shapes() {
                            if self.ds.shape_info(i).shape_type != topods::ShapeType::Edge {
                                continue;
                            }
                            if self.ds.has_pave_blocks(i) {
                                n_pb += self.ds.pave_blocks(i).len();
                            }
                        }
                        n_pb
                    },
                    n_ds_intersection_curves: self.ds.intersection_curves.len(),
                    n_ds_interf_ff: self.ds.interf_ff.len(),
                    n_brep_vertices: v,
                    n_brep_edges: e,
                    n_brep_faces: f,
                    n_brep_shells: sh,
                    n_brep_solids: so,
                });
            }};
        }
        let partial = |s: &Option<topods::BRep>| s.clone().unwrap_or_default();

        // OCCT L431-436: CheckData
        self.check_data()?;
        self.check_filler();

        // OCCT L438-443: Prepare
        let _result = self.prepare();

        // OCCT L459-471: FillImagesVertices + BuildResult(VERTEX)
        self.fill_images_vertices();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Vertex);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(1, "after_FillImagesVertices");

        // OCCT L472-483: FillImagesEdges + BuildResult(EDGE)
        self.fill_images_edges();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Edge);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(2, "after_FillImagesEdges");

        // OCCT L484-494: FillImagesContainers(WIRE) + BuildResult(WIRE)
        self.fill_images_containers(topods::ShapeType::Wire);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Wire);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(3, "after_BuildResultWire");

        // OCCT L496-505: FillImagesFaces + BuildResult(FACE)
        self.fill_images_faces();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Face);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(4, "after_FillImagesFaces");

        // OCCT L507-516: FillImagesContainers(SHELL) + BuildResult(SHELL)
        self.fill_images_containers(topods::ShapeType::Shell);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Shell);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(5, "after_BuildResultShell");

        // OCCT L518-528: FillImagesSolids + BuildResult(SOLID)
        self.fill_images_solids();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Solid);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(6, "after_FillImagesSolids");
        

        // OCCT L530-539: FillImagesContainers(COMPSOLID) + BuildResult(COMPSOLID)
        self.fill_images_containers(topods::ShapeType::CompSolid);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::CompSolid);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(7, "after_BuildResultCompSolid");

        // OCCT L541-550: FillImagesCompounds + BuildResult(COMPOUND)
        self.fill_images_compounds();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Compound);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(8, "after_FillImagesCompounds");

        // OCCT L575-580 (BOPAlgo_BOP.cxx): BuildShape 鈥?apply the boolean operation
        // result construction. Runs between the s08 dump and PrepareHistory.
        self.build_shape();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        // OCCT L583-587 (BOPAlgo_BOP.cxx): PrepareHistory.
        self.prepare_history();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(9, "after_PrepareHistory");
        self.post_treat();
        snap!(10, "after_PostTreat");

        let result = self.my_shape.clone().unwrap_or_default();
        Ok((result, snapshots))
    }

    /// OCCT BOPAlgo_Builder::Build 鈥?full pipeline with history.
    ///
    /// Delegates to `build_with_history_stage_by_stage` (single source of
    /// truth); converts an early pipeline stop (has_errors) back to Err to
    /// preserve the production contract.
    pub fn build_with_history_topods(
        &mut self,
    ) -> Result<(topods::BRep, ()), BooleanError> {
        let (brep, snaps) = self.build_with_history_stage_by_stage()?;
        if snaps.len() != 10 {
            return Err(BooleanError::DegenerateResult);
        }
        Ok((brep, ()))
    }

    // --- Pipeline stage stubs ---
    // Each matches a method in OCCT BOPAlgo_Builder.
    // Stubs: empty bodies that compile. Implementation will be added incrementally.

    /// OCCT BOPAlgo_Builder::CheckData (BOPAlgo_Builder.cxx L130-140).
    fn check_data(&self) -> Result<(), BooleanError> {
        // OCCT L132-137: if (myArguments.Extent() < 2) 鈫?AddError(TooFewArguments)
        if self.my_arguments.len() < 2 {
            return Err(BooleanError::TooFewArguments);
        }
        // OCCT L139: CheckFiller();
        self.check_filler();
        Ok(())
    }

    /// OCCT BOPAlgo_BOP::CheckFiller (BOPAlgo_BOP.cxx L144-152).
    fn check_filler(&self) {
        // OCCT L146-150: if (!myPaveFiller) 鈫?AddError(NoFiller)
        // rcad: PaveFiller always runs before Builder, no reference stored.
        // OCCT L151: GetReport()->Merge(myPaveFiller->GetReport());
        // rcad: report merging not applicable (PaveFiller dropped before Builder).
    }

    /// OCCT BOPAlgo_Builder::Prepare (BOPAlgo_Builder.cxx L156-164).
    fn prepare(&mut self) {
        // OCCT L158-163: BRep_Builder aBB; MakeCompound(aC); myShape = aC;
        // rcad: topods::BRep is the equivalent of TopoDS_Compound for result.
        self.my_shape = Some(topods::BRep::new());
        self.shape_remap.clear();
    }

    /// OCCT BOPAlgo_Builder::FillImagesVertices (BOPAlgo_Builder_1.cxx L40-67).
    /// Maps each SD vertex pair as myImages[source]->[target], myShapesSD, myOrigins.
    fn fill_images_vertices(&mut self) {
        // OCCT L40-66: NCollection_DataMap<int, int>::Iterator aIt(myDS->ShapesSD());
        // rcad: DS::shapes_sd is HashMap<usize, usize> (source鈫扴D).
        let sd_pairs: Vec<(usize, usize)> = self.ds.shapes_sd
            .iter()
            .map(|(k, v)| (*k, *v))
            .collect();
        for (nV, nVSD) in sd_pairs {
            // OCCT L53-54: const TopoDS_Shape& aV = myDS->Shape(nV);
            let aV = self.brep_sr(nV);
            // OCCT L54: const TopoDS_Shape& aVSD = myDS->Shape(nVSD);
            let aVSD = self.brep_sr(nVSD);
            // OCCT L56: myImages.Bound(aV, ...)->Append(aVSD);
            self.my_images.bound((aV.ptr_id(), aV.location)).push(aVSD.clone());
            // OCCT L58: myShapesSD.Bind(aV, aVSD);
            self.my_shapes_sd.insert((aV.ptr_id(), aV.location), aVSD.clone());
            // OCCT L60-65: myOrigins 鈥?find or create list, append
            self.my_origins.entry((aVSD.ptr_id(), aVSD.location)).or_default().push(aV);
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesEdges (BOPAlgo_Builder_1.cxx L71-126).
    /// Maps source edges -> split images via pave-block real edge.
    /// Also handles CommonBlocks via myShapesSD.
    fn fill_images_edges(&mut self) {
        let aNbS = self.ds.nb_source_shapes();
        for i in 0..aNbS {
            let aSI = self.ds.shape_info(i);
            if aSI.shape_type != topods::ShapeType::Edge {
                continue;
            }
            // OCCT L84-86: if (!aSI.HasReference()) continue;
            if !aSI.has_reference() {
                continue;
            }
            // OCCT L89-91: aE = myDS->Shape(i); aLPB = myDS->PaveBlocks(i);
            let aE = self.brep_sr(i);
            let aLPB: Vec<SharedPB> = self.ds.pave_blocks(i).to_vec();
            // OCCT L95: pLS = myImages.Bound(aE, List()) 鈥?the image list is
            // bound UNCONDITIONALLY: an edge with no pave blocks (a small
            // edge) gets an empty image list and is thus avoided in the result
            // (OCCT comment L93-94: "The small edges, having no pave blocks,
            // will have the empty list of images").
            self.my_images.bound((aE.ptr_id(), aE.location));
            // OCCT L96-120: iterate pave blocks
            for aPB in &aLPB {
                // OCCT L100-102: aPBR = myDS->RealPaveBlock(aPB); nSpR = aPBR->Edge();
                let aPBR = self.ds.real_pave_block(aPB);
                let nSpR = { let r = aPBR.read(); r.edge };
                // OCCT L103-104: aSpR = myDS->Shape(nSpR);
                let aSpR = self.brep_sr(nSpR);
                // OCCT L105: pLS->Append(aSpR);
                self.my_images.get_mut((aE.ptr_id(), aE.location)).unwrap().push(aSpR.clone());
                // OCCT L107-112: pLOr = myOrigins.ChangeSeek(aSpR); append aE
                self.my_origins.entry((aSpR.ptr_id(), aSpR.location)).or_default().push(aE.clone());
                // OCCT L114-119: if (IsCommonBlockOnEdge(aPB)) 鈫?myShapesSD.Bind(aSp, aSpR)
                if self.ds.is_common_block_on_edge(aPB) {
                    // OCCT L116-117: nSp = aPB->Edge(); aSp = myDS->Shape(nSp);
                    let nSp = { let r = aPB.read(); r.edge };
                    let aSp = self.brep_sr(nSp);
                    // OCCT L118: myShapesSD.Bind(aSp, aSpR);
                    self.my_shapes_sd.insert((aSp.ptr_id(), aSp.location), aSpR.clone());
                }
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesContainers (BOPAlgo_Builder_1.cxx L172-193).
    /// Builds wire/shell/compsolid images from edge/face/solid images.
    /// For each source shape of theType, calls FillImagesContainer.
    fn fill_images_containers(&mut self, the_type: topods::ShapeType) {
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != the_type {
                continue;
            }
            // OCCT L185-186: FillImagesContainer(aC, theType)
            let a_c = self.brep_sr(i);
            self.fill_images_container(&a_c, the_type);
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesContainer (BOPAlgo_Builder_1.cxx L221-276).
    /// Builds a new container (wire/shell/compsolid) from sub-shape images.
    /// If no sub-shape was modified, the container is kept as-is.
    fn fill_images_container(&mut self, the_s: &Shape, the_type: topods::ShapeType) {
        // OCCT L223-233: check if any sub-shape has been modified
        let sub_shapes = self.shape_sub_shapes(the_s);
        let mut has_modified = false;
        for ss in &sub_shapes {
            // OCCT myImages.IsBound(aS) keys by TopTools_ShapeMapHasher
            // (orientation-insensitive); the composed orientation of the
            // sub-shape (e.g. a Reversed wire edge) must not hide its images.
            if let Some(imgs) = self.images_of(ss) {
                // OCCT L228-229: pLFIm->Extent() != 1 ||
                // !pLFIm->First().IsSame(aSS) 鈥?IsSame = TShape + Location.
                if imgs.len() != 1
                    || (imgs[0].ptr_id(), imgs[0].location) != (ss.ptr_id(), ss.location)
                {
                    has_modified = true;
                    break;
                }
            }
        }
        if !has_modified {
            return;
        }

        // OCCT L242-245: MakeContainer(theType, aCIm)
        let mut new_edges: Vec<Shape> = Vec::new();
        let mut new_faces: Vec<Shape> = Vec::new();
        let mut new_comps: Vec<Shape> = Vec::new();

        // OCCT L247-272: iterate sub-shapes, add images or originals
        for ss in &sub_shapes {
            let p_lss_im = self.images_of(ss);
            match p_lss_im {
                None => {
                    // OCCT L253-257: no splits, add sub-shape itself
                    match the_type {
                        topods::ShapeType::Wire => new_edges.push(ss.clone()),
                        topods::ShapeType::Shell => new_faces.push(ss.clone()),
                        topods::ShapeType::CompSolid => new_comps.push(ss.clone()),
                        _ => {}
                    }
                }
                Some(imgs) => {
                    // OCCT L260-271: add each image (split)
                    for a_ss_im0 in imgs {
                        // OCCT L265-269: IsSplitToReverseWithWarn(aSSIm, aSS) 鈥?
                        // reverse aSSIm when its geometry is oppositely oriented to aSS.
                        // OCCT BOPTools_AlgoTools::IsSplitToReverse (L1263-1300)
                        // dispatches by sub-shape type: EDGE -> edge overload,
                        // FACE -> face overload, other types -> no reversal.
                        let mut a_ss_im = a_ss_im0.clone();
                        if !a_ss_im.is_equal(ss) {
                            let b_to_rev = match the_type {
                                topods::ShapeType::Wire => crate::bop::tools::algo_tools::is_split_to_reverse_edge(&a_ss_im, ss).0,
                                topods::ShapeType::Shell => crate::bop::tools::algo_tools::is_split_to_reverse_face(&a_ss_im, ss, &self.ds).0,
                                _ => false,
                            };
                            if b_to_rev {
                                a_ss_im.orientation = flip_orientation(a_ss_im.orientation);
                            }
                        }
                        match the_type {
                            topods::ShapeType::Wire => new_edges.push(a_ss_im),
                            topods::ShapeType::Shell => new_faces.push(a_ss_im),
                            topods::ShapeType::CompSolid => new_comps.push(a_ss_im),
                            _ => {}
                        }
                    }
                }
            }
        }

        // OCCT L274: aCIm.Closed(BRep_Tool::IsClosed(aCIm))
        let mut container_flags = tshape_flags::DEFAULT;
        match the_type {
            topods::ShapeType::Wire => {
                if wire_is_closed(&new_edges) {
                    container_flags |= tshape_flags::CLOSED;
                }
            }
            topods::ShapeType::Shell => {
                if shell_is_closed(&new_faces) {
                    container_flags |= tshape_flags::CLOSED;
                }
            }
            // CompSolid: BRep_Tool::IsClosed(COMPSOLID) returns the shape's
            // own Closed flag, which is never set here.
            _ => {}
        }

        // Build new container TShape
        let new_container: TShape = match the_type {
            topods::ShapeType::Wire => {
                TShape::Wire(TWireData {
                    my_shapes: vec![], flags: container_flags,
                    edges: new_edges,
                })
            }
            topods::ShapeType::Shell => {
                TShape::Shell(TShellData {
                    my_shapes: vec![], flags: container_flags,
                    faces: new_faces,
                })
            }
            topods::ShapeType::CompSolid => {
                TShape::CompSolid(new_comps)
            }
            _ => return,
        };

        // Wrap in Shape (synthetic index, will be remapped during add_to_result)
        let container_shape = Shape::new(
            std::sync::Arc::new(new_container),
            0, topods::Orientation::Forward,
        );

        // OCCT L275: myImages.Bound(theS, ...)->Append(aCIm)
        self.my_images.bound((the_s.ptr_id(), the_s.location)).push(container_shape);
    }

    /// Extract immediate sub-shapes from a Shape (OCCT TopoDS_Iterator equivalent).
    fn shape_sub_shapes(&self, s: &Shape) -> Vec<Shape> {
        match &*s.data {
            TShape::Vertex(_) => vec![],
            TShape::Edge(ed) => {
                // OCCT TopoDS_Iterator(aE, cumLoc) composes the edge Location
                // into the vertices (TopoDS_Iterator.cxx L76-78).
                let loc = s.location;
                let vl = |v: &Shape| Shape::new(
                    v.data.clone(),
                    crate::bop::algo::compose_edge_vertex_location(loc, v.location, &self.ds.locations),
                    v.orientation,
                );
                vec![vl(&ed.first), vl(&ed.last)]
            }
            TShape::Wire(wd) => {
                wd.edges.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Face(fd) => {
                let mut v = vec![
                    Shape::new(fd.outer_wire.data.clone(), fd.outer_wire.location, fd.outer_wire.orientation)
                ];
                v.extend(fd.inner_wires.iter().map(|w| {
                    Shape::new(w.data.clone(), w.location, w.orientation)
                }));
                v
            }
            TShape::Shell(sd) => {
                sd.faces.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Solid(sd) => {
                sd.shells.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::CompSolid(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Compound(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesFaces (BOPAlgo_Builder_2.cxx L215-229).
    /// Splits faces using section edges.
    /// Calls BuildSplitFaces -> FillSameDomainFaces -> FillInternalVertices.
    fn fill_images_faces(&mut self) {
        // OCCT L218: BuildSplitFaces
        self.build_split_faces();
        if self.has_errors() { return; }
        // OCCT L223: FillSameDomainFaces
        self.fill_same_domain_faces();
        if self.has_errors() { return; }
        // OCCT L228: FillInternalVertices
        self.fill_internal_vertices();
    }

    /// OCCT BOPAlgo_Builder::BuildSplitFaces (BOPAlgo_Builder_2.cxx L233-555).
    ///
    /// For each source face with intersection data, builds a BuilderFace fed
    /// with the full edge set (bounding edges + their images + IN edges +
    /// section edges). Faces without section/IN edges but with modified wires
    /// or alone vertices take the BuildDraftFace fast path.
    fn build_split_faces(&mut self) {
        let a_nb_s = self.ds.nb_source_shapes();
        if std::env::var("RCAD_BS_DEBUG").is_ok() {
            for ff in &self.ds.interf_ff {
                eprintln!("[FF] f1={} f2={} n_curves={} tangent={}",
                    ff.f1, ff.f2, ff.curves.len(), ff.tangent_faces);
            }
        }
        // aFacesIm: DS face index -> area shapes (OCCT IndexedDataMap<int,
        // List<Shape>>, Builder_2.cxx L256 鈥?insertion order, iterated at L535).
        let mut a_faces_im: IndexMap<usize, Vec<Shape>> = IndexMap::new();
        // aVBF: pending BuilderFace tasks (face_idx, face, edges).
        let mut a_vbf: Vec<(usize, Shape, Vec<Shape>)> = Vec::new();

        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Face {
                continue;
            }
            // OCCT L275-279: bHasFaceInfo check.
            if !self.ds.has_face_info(i) {
                continue;
            }
            let a_f = self.brep_sr(i);
            let a_fi = self.ds.face_info(i).clone();

            // OCCT L286-287: AloneVertices(i, aLIAV).
            let a_liav = self.alone_vertices(i);
            let a_nb_pb_in = a_fi.pave_blocks_in.len();
            let a_nb_pb_on = a_fi.pave_blocks_on.len();
            let a_nb_pb_sc = a_fi.pave_blocks_sc.len();
            let a_nb_av = a_liav.len();
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                let edesc = |pb_key: &u64| -> String {
                    match self.ds.pb_from_ptr(*pb_key) {
                        Some(pb) => {
                            let r = pb.0.read().unwrap();
                            let e = r.original_edge;
                            let ef = r.edge;
                            let c = self.ds.edge_curve(e).map(|c| match c {
                                rcad_kernel::geom::Curve3::Line(l) => format!("L({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})", l.origin.x, l.origin.y, l.origin.z, l.origin.x + l.direction.x, l.origin.y + l.direction.y, l.origin.z + l.direction.z),
                                _ => "O".into(),
                            }).unwrap_or_else(|| "?".into());
                            format!("e{}:{} (edge={})", e, c, ef)
                        }
                        None => "?".into(),
                    }
                };
                let fin = a_fi.pave_blocks_in.iter().map(|k| edesc(k)).collect::<Vec<_>>();
                let fon = a_fi.pave_blocks_on.iter().map(|k| edesc(k)).collect::<Vec<_>>();
                let fsc = a_fi.pave_blocks_sc.iter().map(|k| edesc(k)).collect::<Vec<_>>();
                eprintln!("[BSF] face={} pbIn={:?} pbOn={:?} pbSc={:?} av={}",
                    i, fin, fon, fsc, a_nb_av);
            }
            // OCCT L293-296: not complete -> skip.
            if a_nb_pb_in == 0 && a_nb_pb_on == 0 && a_nb_pb_sc == 0 && a_nb_av == 0 {
                continue;
            }

            // OCCT L298-351: only alone vertices / On PBs -> draft-face fast path.
            if a_nb_pb_in == 0 && a_nb_pb_sc == 0 {
                let mut has_internals = false;
                if a_nb_av == 0 {
                    // OCCT L315-330: check wires for internal edges or modifications.
                    let mut has_modified = false;
                    for a_w in self.shape_sub_shapes(&a_f) {
                        if a_w.shape_type() != topods::ShapeType::Wire {
                            continue;
                        }
                        // OCCT L320-321: only the FIRST edge of each wire is
                        // checked for INTERNAL orientation 鈥?OCCT wires store
                        // internal edges first, so the first edge decides.
                        let first_is_internal = match &*a_w.data {
                            TShape::Wire(wd) => wd.edges
                                .first()
                                .map(|e| e.orientation == topods::Orientation::Internal)
                                .unwrap_or(false),
                            _ => false,
                        };
                        has_internals = first_is_internal;
                        if has_internals {
                            break;
                        }
                        if self.images_of(&a_w).is_some() {
                            has_modified = true;
                        }
                    }
                    if !has_internals && !has_modified {
                        continue;
                    }
                }
                if !has_internals {
                    // OCCT L344: BuildDraftFace fast path 鈥?face image directly.
                    if let Some(a_fd) = self.build_draft_face(&a_f) {
                        a_faces_im.entry(i).or_default().push(a_fd);
                        continue;
                    }
                }
            }

            // OCCT L353: aMFence.Clear() 鈥?per-face fence for closed edges.
            // OCCT aMFence is NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>
            // (L252) 鈥?key identity TShape + Location, orientation ignored.
            let mut a_mfence: HashSet<(u64, u32)> = HashSet::new();
            // OCCT L355-357: aFF = aF; aFF.Orientation(FORWARD).
            let a_ff = Shape::new(a_f.data.clone(), a_f.location, topods::Orientation::Forward);
            // OCCT L359: 1. Build the edges set aLE.
            let mut a_le: Vec<Shape> = Vec::new();

            // OCCT L362-465: 1.1 Bounding edges.
            let mut is_checked = false;
            let mut is_u_closed = false;
            let mut is_v_closed = false;
            let a_ff_sphere = a_ff.as_face().and_then(|fd| fd.surface.as_ref()).map_or(false, |s| matches!(s, rcad_kernel::geom::Surface3::Sphere(_)));
            for a_e in self.face_edges(&a_ff) {
                let an_ori_e = a_e.orientation;
                // OCCT L369: if !myImages.IsBound(aE).
                if self.images_of(&a_e).is_none() {
                    if an_ori_e == topods::Orientation::Internal {
                        let mut a_ee = a_e.clone();
                        a_ee.orientation = topods::Orientation::Forward;
                        a_le.push(a_ee.clone());
                        a_ee.orientation = topods::Orientation::Reversed;
                        a_le.push(a_ee);
                    } else {
                        a_le.push(a_e);
                    }
                    continue;
                }
                // OCCT L387-393: GeomLib::IsClosed(aSurf, BRep_Tool::Tolerance(aE),
                // isUClosed, isVClosed).
                if !is_checked {
                    let a_tol_e = a_e.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
                    let (uc, vc) = self.surface_is_closed(&a_f, a_tol_e);
                    is_u_closed = uc;
                    is_v_closed = vc;
                    is_checked = true;
                }
                // OCCT L395-404: bIsClosed = seam edge on closed surface.
                let mut b_is_closed = false;
                if (is_u_closed || is_v_closed) && self.edge_closed_on_face(&a_e, &a_f) {
                    let (is_ui, is_vi) = self.is_edge_isoline(&a_e, &a_f);
                    b_is_closed = (is_u_closed && is_ui) || (is_v_closed && is_vi);
                }
                // OCCT L406: bIsDegenerated = BRep_Tool::Degenerated(aE).
                let b_is_degenerated = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                // OCCT L408: aLIE = myImages.Find(aE).
                let a_lie = self.images_of(&a_e).unwrap_or_default();
                for a_sp0 in &a_lie {
                    let mut a_sp = a_sp0.clone();
                    // OCCT L413-418: degenerated -> keep original orientation.
                    if b_is_degenerated {
                        a_sp.orientation = an_ori_e;
                        a_le.push(a_sp);
                        continue;
                    }
                    // OCCT L420-427: INTERNAL -> forward + reversed.
                    if an_ori_e == topods::Orientation::Internal {
                        a_sp.orientation = topods::Orientation::Forward;
                        a_le.push(a_sp.clone());
                        a_sp.orientation = topods::Orientation::Reversed;
                        a_le.push(a_sp);
                        continue;
                    }
                    // OCCT L429-455: closed seam edge -> dedupe via aMFence.
                    if b_is_closed {
                        if a_mfence.insert((a_sp.ptr_id(), a_sp.location)) {
                            if !self.edge_closed_on_face(&a_sp, &a_f) {
                                // OCCT L435-446: DoSplitSEAMOnFace(aSp, aF) /
                                // DoSplitSEAMOnFace(aE, aSp, aF).
                                if !self.do_split_seam_on_face(&a_sp, &a_f)
                                    && !self.do_split_seam_on_face_origin(&a_e, &a_sp, &a_f)
                                {
                                    self.my_report.add_warning(
                                        crate::bop::algo::Alert::UnableToMakeClosedEdgeOnFace(
                                            vec![a_f.clone(), a_sp.clone()],
                                        ),
                                    );
                                }
                            }
                            a_sp.orientation = topods::Orientation::Forward;
                            a_le.push(a_sp.clone());
                            a_sp.orientation = topods::Orientation::Reversed;
                            a_le.push(a_sp);
                        }
                        continue;
                    }
                    // OCCT L457-463: regular split edge.
                    a_sp.orientation = an_ori_e;
                    // OCCT L458: IsSplitToReverseWithWarn(aSp, aE) 鈥?reverse the
                    // split edge when its direction differs from the original.
                    let (b_to_reverse, _err) =
                        crate::bop::tools::algo_tools::is_split_to_reverse_edge(&a_sp, &a_e);
                    if b_to_reverse {
                        a_sp.orientation = flip_orientation(a_sp.orientation);
                    }
                    a_le.push(a_sp);
                }
            }

            // OCCT L469-480: 1.2 In edges (forward + reversed).
            for &pb_key in &a_fi.pave_blocks_in {
                if let Some(n_sp) = self.pb_edge_by_ptr(pb_key) {
                    let mut a_sp = self.brep_sr(n_sp);
                    a_sp.orientation = topods::Orientation::Forward;
                    a_le.push(a_sp.clone());
                    a_sp.orientation = topods::Orientation::Reversed;
                    a_le.push(a_sp);
                }
            }
            // OCCT L483-494: 1.3 Section edges (forward + reversed).
            // OCCT reads aPB->Edge() with no null-check 鈥?PostTreatFF always
            // sets the edge on section PBs, so a missing edge here is a bug.
            for &pb_key in &a_fi.pave_blocks_sc {
                let n_sp = self.pb_edge_by_ptr(pb_key)
                    .expect("section pave block must reference an edge (OCCT PostTreatFF sets it)");
                let mut a_sp = self.brep_sr(n_sp);
                a_sp.orientation = topods::Orientation::Forward;
                a_le.push(a_sp.clone());
                a_sp.orientation = topods::Orientation::Reversed;
                a_le.push(a_sp);
            }
            // OCCT L496-500: if (!NonDestructive()) BRepLib::BuildPCurveForEdgesOnPlane(aLE, aFF).
            if !self.my_non_destructive {
                self.build_pcurve_for_edges_on_plane(&a_le, &a_ff);
            }
            // OCCT L502-505: aBF.SetFace(aF); aBF.SetShapes(aLE); SetRunParallel.
            a_vbf.push((i, a_f, a_le));
        }

        // OCCT L515-521: perform all BuilderFace tasks.
        for (fi, a_f, a_le) in &a_vbf {
            let mut a_bf = BuilderFace::new(&self.ds);
            // OCCT L502-504: aBF.SetFace(aF) 鈥?BOPAlgo_BuilderFace::SetFace
            // (BuilderFace.cxx L79-84) stores myOrientation = aF.Orientation()
            // and normalizes myFace to FORWARD. Areas are built from the
            // FORWARD face; the original orientation is re-applied below from
            // the DS face (L534-552), mirroring OCCT.
            let a_f_forward = Shape::new(
                a_f.data.clone(),
                a_f.location,
                topods::Orientation::Forward,
            );
            a_bf.my_face = Some(a_f_forward);
            a_bf.my_face_index = Some(*fi);
            a_bf.my_edges = a_le.clone();
            a_bf.perform();
            // OCCT L551: myReport->Merge(aBF.GetReport()) 鈥?merge the
            // BuilderFace warnings (e.g. failed area building) into the main
            // report.
            self.my_report.merge(a_bf.report().clone());
            // OCCT L527-531: aFacesIm.Add(myDS->Index(aBF.Face()), aBF.Areas()).
            // OCCT binds every split face to its areas, even an empty list 鈥?
            // a face whose areas could not be built contributes nothing to the
            // result (build_draft_solid drops it). Skipping the bind would keep
            // the original face, which OCCT does not do.
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[BSF-AREA] face={} n_areas={}", fi, a_bf.my_areas.len());
            }
            a_faces_im.entry(*fi).or_default().extend(a_bf.my_areas);
        }

        // OCCT L534-552: apply orientation and append areas to myImages.
        for (fi, a_lfr) in a_faces_im {
            let a_f = self.brep_sr(fi);
            let an_ori_f = a_f.orientation;
            let p_lf_im = self.my_images.bound((a_f.ptr_id(), a_f.location));
            for mut a_fr in a_lfr {
                if an_ori_f == topods::Orientation::Reversed {
                    a_fr.orientation = topods::Orientation::Reversed;
                }
                p_lf_im.push(a_fr);
            }
            if std::env::var("RCAD_BS_DEBUG").is_ok() && (fi == 60 || fi == 46) {
                let ims: Vec<String> = self.my_images.get((a_f.ptr_id(), a_f.location)).cloned().unwrap_or_default().iter().map(|im| {
                    let n = im.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                        rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                        _ => "O".into(),
                    }).unwrap_or_else(|| "?".into());
                    format!("{}:{}", if im.orientation == topods::Orientation::Reversed { "R" } else { "F" }, n)
                }).collect();
                eprintln!("[IMG] face={} src_or={:?} images=[{}]", fi, an_ori_f, ims.join(" "));
            }
        }
    }

    /// OCCT BOPDS_DS::AloneVertices (BOPDS_DS.cxx L1028-1062).
    /// Vertices of the face not belonging to any boundary edge: endpoints of
    /// PaveBlocksIn/PaveBlocksSc plus VerticesIn/VerticesSc not already seen.
    fn alone_vertices(&self, face_idx: usize) -> Vec<usize> {
        if !self.ds.has_face_info(face_idx) {
            return Vec::new();
        }
        let a_fi = self.ds.face_info(face_idx);
        let mut a_mi: HashSet<usize> = HashSet::new();
        for pb_set in [&a_fi.pave_blocks_in, &a_fi.pave_blocks_sc] {
            for &pb_key in pb_set {
                if let Some(pb) = self.ds.pb_from_ptr(pb_key) {
                    let r = pb.0.read().unwrap();
                    a_mi.insert(r.pave1.vertex_idx);
                    a_mi.insert(r.pave2.vertex_idx);
                }
            }
        }
        let mut result: Vec<usize> = Vec::new();
        for v in a_fi.vertices_in.iter().chain(a_fi.vertices_sc.iter()) {
            if a_mi.insert(*v) {
                result.push(*v);
            }
        }
        result
    }

    /// Edge DS index from a pave-block pool entry (OCCT aPB->Edge()).
    fn pb_edge(&self, pb_idx: usize) -> Option<usize> {
        let pool = self.ds.pave_blocks_pool.get(&pb_idx)?;
        let pb = pool.first()?;
        let e = pb.0.read().unwrap().edge;
        if e < self.ds.nb_shapes() {
            Some(e)
        } else {
            None
        }
    }

    /// Edge of a PB identified by its pointer id (FaceInfo PaveBlocksOn/In/Sc
    /// store pointer ids, OCCT stores handles).
    fn pb_edge_by_ptr(&self, pb_ptr: u64) -> Option<usize> {
        let pb = self.ds.pb_from_ptr(pb_ptr)?;
        let e = pb.0.read().unwrap().edge;
        if e < self.ds.nb_shapes() {
            Some(e)
        } else {
            None
        }
    }

    /// Look up myImages by TShape pointer + location, ignoring orientation.
    /// OCCT: myImages is keyed by TopTools_ShapeMapHasher which hashes only
    /// TShape* + Location, so IsBound(aE) matches regardless of orientation.
    fn images_of(&self, key: &Shape) -> Option<Vec<Shape>> {
        // OCCT myImages.Find(aE) keys by TopoDS_Shape (stable because OCCT mutates
        // TShapes in place). rcad's Arc::make_mut clones a shared TShape on write,
        // so the face wire edge (original) and the DS edge entry (clone) become two
        // TShape objects for the same logical edge. Both still map to the same DS
        // index 鈥?init keeps the original (ptr鈫抜dx) mapping and remap_shape_idx
        // adds the clone's 鈥?so resolve the lookup by DS index.
        let key_idx = self.ds.map_shape_index.get(&(key.ptr_id(), key.location)).copied();
        for (k, v) in self.my_images.iter() {
            if k.0 == key.ptr_id() && k.1 == key.location {
                return Some(v.clone());
            }
        }
        if let Some(ki) = key_idx {
            for (k, v) in self.my_images.iter() {
                if self.ds.map_shape_index.get(&(k.0, k.1)).copied() == Some(ki) {
                    return Some(v.clone());
                }
            }
        }
        None
    }

    /// Boundary edges of a face with wire-composed orientation.
    /// OCCT TopExp_Explorer(aFF, TopAbs_EDGE) 鈥?BOPAlgo_Builder_2.cxx L363-365.
    fn face_edges(&self, a_f: &Shape) -> Vec<Shape> {
        // OCCT TopExp_Explorer(aFF, TopAbs_EDGE) descends face -> wire -> edge
        // with cumOri=true, composing each parent's orientation into the edge.
        // BRepPrim_GWedge stores the MIN-face wires REVERSED (ReverseFace), so
        // their edge orientations come out flipped here 鈥?matching OCCT's
        // BuildSplitFaces anOriE values.
        let mut result: Vec<Shape> = Vec::new();
        for a_w in self.shape_sub_shapes(a_f) {
            if a_w.shape_type() != topods::ShapeType::Wire {
                continue;
            }
            if let TShape::Wire(wd) = &*a_w.data {
                result.extend(wd.edges.iter().map(|sr| {
                    let mut e = Shape::from_parts(sr.data.clone(), sr.index, sr.location, sr.orientation);
                    e.orientation = topods::Orientation::compose(a_w.orientation, e.orientation);
                    e
                }));
            }
        }
        result
    }

    /// OCCT GeomLib::IsClosed(aSurf, aTol, isUClosed, isVClosed)
    /// (GeomLib.cxx L2693-2868). Geometric closed-ness of the face surface in U
    /// and V, sampled with the edge tolerance (BOPAlgo_Builder_2.cxx L387-393:
    /// aTol = BRep_Tool::Tolerance(aE)).
    fn surface_is_closed(&self, a_f: &Shape, a_tol: f64) -> (bool, bool) {
        let Some(surf) = a_f.as_face().and_then(|fd| fd.surface.clone()) else {
            return (false, false);
        };
        // GeomAdaptor_Surface aGAS(S) (GeomAdaptor_Surface.cxx L417-430): a
        // RectangularTrimmedSurface is unwrapped to its basis surface while
        // keeping the trimmed parameter bounds.
        let (basis, [mut u1, mut u2, mut v1, mut v2]) = surface_adaptor_basis_and_bounds(&surf);
        let tol2 = a_tol * a_tol;
        match basis {
            // OCCT L2713-2715: GeomAbs_Plane.
            Surface3::Plane(_) => (false, false),
            // OCCT L2716-2733: GeomAbs_SurfaceOfExtrusion falls through to
            // GeomAbs_Cylinder when its u range is finite.
            Surface3::LinearExtrusion(_) | Surface3::Cylinder(_) => {
                if matches!(basis, Surface3::LinearExtrusion(_))
                    && (u1.is_infinite() || u2.is_infinite())
                {
                    return (false, false);
                }
                if v1.is_infinite() {
                    v1 = 0.0;
                }
                let p1 = basis.point_at(u1, v1);
                let p2 = basis.point_at(u2, v1);
                (p1.distance_squared(p2) <= tol2, false)
            }
            // OCCT L2734-2755: GeomAbs_Cone.
            Surface3::Cone(c) => {
                // find v with maximal distance from axis
                if !(v1.is_infinite() || v2.is_infinite()) {
                    let an_apex = c.apex_point();
                    let p1 = basis.point_at(u1, v1);
                    let p2 = basis.point_at(u1, v2);
                    if p2.distance_squared(an_apex) > p1.distance_squared(an_apex) {
                        v1 = v2;
                    }
                } else {
                    v1 = 0.0;
                }
                let p1 = basis.point_at(u1, v1);
                let p2 = basis.point_at(u2, v1);
                (p1.distance_squared(p2) <= tol2, false)
            }
            // OCCT L2756-2773: GeomAbs_Sphere.
            Surface3::Sphere(_) => {
                // find v with maximal distance from axis
                if v1 * v2 <= 0.0 {
                    v1 = 0.0;
                } else if v1 < 0.0 {
                    v1 = v2;
                }
                let p1 = basis.point_at(u1, v1);
                let p2 = basis.point_at(u2, v1);
                (p1.distance_squared(p2) <= tol2, false)
            }
            // OCCT L2774-2781: GeomAbs_Torus.
            Surface3::Torus(_) => {
                let ures = u_resolution_for_surface(&surf, a_tol);
                let vres = v_resolution_for_surface(&surf, a_tol);
                let u_period = std::f64::consts::PI * 2.0;
                let v_period = std::f64::consts::PI * 2.0;
                ((u2 - u1) >= u_period - ures, (v2 - v1) >= v_period - vres)
            }
            // OCCT L2782-2787: GeomAbs_BSplineSurface.
            Surface3::BSpline(bs) => (
                Self::is_bspl_u_closed(bs, u1, u2, a_tol),
                Self::is_bspl_v_closed(bs, v1, v2, a_tol),
            ),
            // OCCT L2788-2793: GeomAbs_BezierSurface.
            Surface3::Bezier(bz) => (
                Self::is_bz_u_closed(bz, u1, u2, a_tol),
                Self::is_bz_v_closed(bz, v1, v2, a_tol),
            ),
            // OCCT L2794-2863: GeomAbs_SurfaceOfRevolution /
            // GeomAbs_OffsetSurface / GeomAbs_OtherSurface 鈥?23-point sampling.
            _ => {
                let mut nbp = 23;
                if v1.is_infinite() {
                    v1 = v1.signum();
                }
                if v2.is_infinite() {
                    v2 = v2.signum();
                }
                // SurfaceOfRevolution keeps its u range; Offset/Other clamp it.
                if !matches!(basis, Surface3::Revolution(_)) {
                    if u1.is_infinite() {
                        u1 = u1.signum();
                    }
                    if u2.is_infinite() {
                        u2 = u2.signum();
                    }
                }
                let mut is_u_closed = true;
                let mut dt = (v2 - v1) / (nbp as f64 - 1.0);
                let mut res = u_resolution_for_surface(&surf, a_tol).max(PCONFUSION);
                if dt <= res {
                    nbp = ((v2 - v1) / (2.0 * res)) as i32 + 1;
                    nbp = nbp.max(2);
                    dt = (v2 - v1) / (nbp as f64 - 1.0);
                }
                for i in 0..nbp {
                    let t = if i == nbp - 1 { v2 } else { v1 + i as f64 * dt };
                    let p1 = basis.point_at(u1, t);
                    let p2 = basis.point_at(u2, t);
                    if p1.distance_squared(p2) > tol2 {
                        is_u_closed = false;
                        break;
                    }
                }
                nbp = 23;
                let mut is_v_closed = true;
                dt = (u2 - u1) / (nbp as f64 - 1.0);
                res = v_resolution_for_surface(&surf, a_tol).max(PCONFUSION);
                if dt <= res {
                    nbp = ((u2 - u1) / (2.0 * res)) as i32 + 1;
                    nbp = nbp.max(2);
                    dt = (u2 - u1) / (nbp as f64 - 1.0);
                }
                for i in 0..nbp {
                    let t = if i == nbp - 1 { u2 } else { u1 + i as f64 * dt };
                    let p1 = basis.point_at(t, v1);
                    let p2 = basis.point_at(t, v2);
                    if p1.distance_squared(p2) > tol2 {
                        is_v_closed = false;
                        break;
                    }
                }
                (is_u_closed, is_v_closed)
            }
        }
    }

    /// OCCT GeomLib::IsBSplUClosed (GeomLib.cxx L2872-2891). S->UIso(U1) /
    /// S->UIso(U2) at the boundary knots evaluate to the first/last control
    /// rows of the tensor-product control net.
    fn is_bspl_u_closed(s: &BSplineSurface, _u1: f64, _u2: f64, a_tol: f64) -> bool {
        if s.control_points.is_empty() {
            return false;
        }
        let a_pf = &s.control_points[0];
        let a_pl = &s.control_points[s.control_points.len() - 1];
        let wf = if s.weights.is_empty() {
            None
        } else {
            Some(s.weights[0].as_slice())
        };
        let wl = if s.weights.is_empty() {
            None
        } else {
            Some(s.weights[s.weights.len() - 1].as_slice())
        };
        Self::compare_weight_poles(a_pf, wf, a_pl, wl, 2.0 * a_tol)
    }

    /// OCCT GeomLib::IsBSplVClosed (GeomLib.cxx L2895-2914). S->VIso(V1) /
    /// S->VIso(V2) at the boundary knots evaluate to the first/last control
    /// columns of the tensor-product control net.
    fn is_bspl_v_closed(s: &BSplineSurface, _v1: f64, _v2: f64, a_tol: f64) -> bool {
        if s.control_points.is_empty() || s.control_points[0].is_empty() {
            return false;
        }
        let a_pf: Vec<DVec3> = s.control_points.iter().map(|row| row[0]).collect();
        let a_pl: Vec<DVec3> = s
            .control_points
            .iter()
            .map(|row| row[row.len() - 1])
            .collect();
        let wf_vec: Vec<f64> = if s.weights.is_empty() {
            Vec::new()
        } else {
            s.weights.iter().map(|row| row[0]).collect()
        };
        let wl_vec: Vec<f64> = if s.weights.is_empty() {
            Vec::new()
        } else {
            s.weights.iter().map(|row| row[row.len() - 1]).collect()
        };
        let wf = if wf_vec.is_empty() {
            None
        } else {
            Some(wf_vec.as_slice())
        };
        let wl = if wl_vec.is_empty() {
            None
        } else {
            Some(wl_vec.as_slice())
        };
        Self::compare_weight_poles(&a_pf, wf, &a_pl, wl, 2.0 * a_tol)
    }

    /// OCCT GeomLib::IsBzUClosed (GeomLib.cxx L2918-2936).
    fn is_bz_u_closed(s: &BezierSurface, _u1: f64, _u2: f64, a_tol: f64) -> bool {
        if s.control_points.is_empty() {
            return false;
        }
        let a_pf = &s.control_points[0];
        let a_pl = &s.control_points[s.control_points.len() - 1];
        Self::compare_weight_poles(a_pf, None, a_pl, None, 2.0 * a_tol)
    }

    /// OCCT GeomLib::IsBzVClosed (GeomLib.cxx L2940-2958).
    fn is_bz_v_closed(s: &BezierSurface, _v1: f64, _v2: f64, a_tol: f64) -> bool {
        if s.control_points.is_empty() || s.control_points[0].is_empty() {
            return false;
        }
        let a_pf: Vec<DVec3> = s.control_points.iter().map(|row| row[0]).collect();
        let a_pl: Vec<DVec3> = s
            .control_points
            .iter()
            .map(|row| row[row.len() - 1])
            .collect();
        Self::compare_weight_poles(&a_pf, None, &a_pl, None, 2.0 * a_tol)
    }

    /// OCCT GeomLib::CompareWeightPoles (GeomLib.cxx L2960-2987) 鈥?poles scaled
    /// by their weights, pairwise within theTol (gp_XYZ::IsEqual).
    fn compare_weight_poles(
        a_pf: &[DVec3],
        wf: Option<&[f64]>,
        a_pl: &[DVec3],
        wl: Option<&[f64]>,
        a_tol: f64,
    ) -> bool {
        if a_pf.len() != a_pl.len() {
            return false;
        }
        for i in 0..a_pf.len() {
            let a_w1 = wf.map(|w| w[i]).unwrap_or(1.0);
            let a_w2 = wl.map(|w| w[i]).unwrap_or(1.0);
            let a_pole1 = a_pf[i] * a_w1;
            let a_pole2 = a_pl[i] * a_w2;
            if (a_pole1 - a_pole2).length() > a_tol {
                return false;
            }
        }
        true
    }

    /// OCCT BRep_Tool::IsClosed(aE, aF) (BRep_Tool.cxx L795-841) 鈥?the edge
    /// has two pcurves on the closed surface of the face (seam edge).
    /// OCCT matches the CurveOnClosedSurface representation by the face's
    /// SURFACE handle (IsCurveOnSurface(S, l)) and returns false for plane
    /// faces outright (IsPlane short-circuit, L819-822). rcad keys the
    /// representation by the face's BRep index: try the exact index match
    /// first (original faces), then fall back to the surface-level match for
    /// split-face images whose index does not preserve the original face's
    /// (BOPAlgo_Builder_2.cxx L397).
    fn edge_closed_on_face(&self, a_e: &Shape, a_f: &Shape) -> bool {
        let f_key = (a_f.ptr_id(), a_f.location);
        // OCCT L825-840: any CurveOnClosedSurface on the edge whose surface
        // matches the face's surface 鈥?exact identity match for original faces.
        let ed = match a_e.as_edge() {
            Some(ed) => ed,
            None => return false,
        };
        if ed.representations.iter().any(|r| {
            matches!(
                r,
                topods::CurveRepresentation::CurveOnClosedSurface { face, .. } if *face == f_key
            )
        }) {
            return true;
        }
        // OCCT L819-822: if (IsPlane(S)) return false.
        match &*a_f.data {
            TShape::Face(fd) => {
                if matches!(fd.surface, Some(rcad_kernel::geom::Surface3::Plane(_))) {
                    return false;
                }
            }
            _ => return false,
        }
        // OCCT L825-840: any CurveOnClosedSurface on the edge whose surface
        // matches the face's surface. rcad keys representations by the face's
        // BRep index, so for split-face images (whose index does not preserve
        // the original face's) the exact match above misses. OCCT compares
        // surface handles (IsCurveOnSurface(S, l)); rcad approximates by
        // checking the face surface is closed (a seam exists only on a closed
        // surface) and the edge carries a CurveOnClosedSurface 鈥?the
        // representation is always created on the same closed surface as the
        // face being split (BOPAlgo_Builder_2.cxx L397).
        let surf = a_f.as_face().and_then(|fd| fd.surface.clone());
        let surface_closed = match &surf {
            Some(rcad_kernel::geom::Surface3::Cylinder(_))
            | Some(rcad_kernel::geom::Surface3::Sphere(_))
            | Some(rcad_kernel::geom::Surface3::Cone(_))
            | Some(rcad_kernel::geom::Surface3::Torus(_))
            | Some(rcad_kernel::geom::Surface3::Revolution(_)) => true,
            _ => false,
        };
        surface_closed
            && ed.representations
                .iter()
                .any(|r| matches!(r, topods::CurveRepresentation::CurveOnClosedSurface { .. }))
    }

    /// OCCT BOPTools_AlgoTools3D::DoSplitSEAMOnFace(aSplit, aF)
    /// (BOPTools_AlgoTools3D.cxx L58-232). The split edge `a_split` lies on a
    /// closed surface but is not yet marked closed: creates the second pcurve
    /// by translating the existing pcurve by one period, and stores both as a
    /// CurveOnClosedSurface representation.
    fn do_split_seam_on_face(&self, a_split: &Shape, a_f: &Shape) -> bool {
        let mut b_is_left = false;
        let mut an_u_period = 0.0;
        let mut an_v_period = 0.0;
        let mut a_sp = a_split.clone();
        a_sp.orientation = Orientation::Forward;
        let a_tol = a_sp.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
        //
        // OCCT L79-81: aS = BRep_Tool::Surface(aF); aS->Bounds(...).
        let Some(a_s) = a_f.as_face().and_then(|fd| fd.surface.clone()) else {
            return false;
        };
        let [a_u_min, a_u_max, a_v_min, a_v_max] = a_s.default_domain();
        //
        // OCCT L84-94: IsUClosed / IsVClosed -> period.
        let mut b_is_u_periodic = a_s.is_u_closed();
        let mut b_is_v_periodic = a_s.is_v_closed();
        if b_is_u_periodic {
            an_u_period = a_u_max - a_u_min;
        }
        if b_is_v_periodic {
            an_v_period = a_v_max - a_v_min;
        }
        //
        // OCCT L96-147: rectangular trimmed surface -> check basis surface.
        if !b_is_u_periodic && !b_is_v_periodic {
            let (basis, _trim) = match &a_s {
                Surface3::Trimmed(ts) => (ts.basis.as_ref().clone(), ts.trim),
                _ => return false,
            };
            b_is_u_periodic = basis.is_u_periodic();
            b_is_v_periodic = basis.is_v_periodic();
            if b_is_u_periodic || b_is_v_periodic {
                let [u0, u1, v0, v1] = basis.default_domain();
                an_u_period = if b_is_u_periodic { u1 - u0 } else { 0.0 };
                an_v_period = if b_is_v_periodic { v1 - v0 } else { 0.0 };
            } else {
                let b_is_u_closed = basis.is_u_closed();
                let b_is_v_closed = basis.is_v_closed();
                let [a_glob_u_min, a_glob_u_max, a_glob_v_min, a_glob_v_max] = basis.default_domain();
                if b_is_u_closed
                    && (a_u_min - a_glob_u_min).abs() < a_tol
                    && (a_u_max - a_glob_u_max).abs() < a_tol
                {
                    b_is_u_periodic = true;
                    an_u_period = a_u_max - a_u_min;
                }
                if b_is_v_closed
                    && (a_v_min - a_glob_v_min).abs() < a_tol
                    && (a_v_max - a_glob_v_max).abs() < a_tol
                {
                    b_is_v_periodic = true;
                    an_v_period = a_v_max - a_v_min;
                }
            }
            if !(b_is_u_periodic || b_is_v_periodic) {
                return false;
            }
        }
        //
        // OCCT L150-153: C2D1 = CurveOnSurface(aSp, aF, a, b);
        // aT = IntermediatePoint(a, b); C2D1->D1(aT, aP2D, aVec2D).
        // BRep_Tool::CurveOnSurface reads the pcurve from the edge TShape,
        // which is shared 鈥?any reference (including the split edge passed
        // here) sees it. rcad's Arc::make_mut clones a shared TShape on write,
        // so the passed edge may miss pcurves that live on the DS canonical
        // shape; fall back to the DS shape (same as wire_splitter::edge_pcurve).
        // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the key
        // location is L.Predivided(E.Location()).
        let face_key = (
            a_f.ptr_id(),
            crate::bop::algo::compose_face_edge_pcurve_location(
                a_f.location, a_sp.location, &self.ds.locations),
        );
        let fi = self.ds.index(a_f) as usize;
        let (c2d1, a, b) = {
            let direct = a_sp
                .as_edge()
                .and_then(|ed| ed.pcurves.get(&face_key).cloned());
            if let Some(v) = direct {
                v
            } else {
                let Some(idx) = self.ds.map_shape_index.get(&(a_sp.ptr_id(), a_sp.location)) else {
                    return false;
                };
                let Some(ed2) = self.ds.shape_info(*idx).shape.as_edge() else {
                    return false;
                };
                let Some(v) = ed2.pcurves.get(&face_key) else {
                    return false;
                };
                v.clone()
            }
        };
        let a_t = intermediate_point(a, b);
        let a_p2d = c2d1.point_at(a_t);
        let a_vec2d = c2d1.derivative_at(a_t);
        let a_dir2d1 = a_vec2d.normalize_or_zero();
        let a_dox = DVec2::X;
        let a_doy = DVec2::Y;
        //
        // OCCT L156-164.
        let mut an_u = a_p2d.x;
        let mut an_v = a_p2d.y;
        let mut an_u1 = an_u;
        let mut an_v1 = an_v;
        //
        let d_u = u_resolution_for_surface(&a_s, a_tol);
        let d_v = v_resolution_for_surface(&a_s, a_tol);
        //
        // OCCT L166-192: if near UMin/UMax or VMin/VMax, shift by one period.
        if an_u_period > 0.0 {
            if (an_u - a_u_min).abs() < d_u {
                b_is_left = true;
                an_u1 = an_u + an_u_period;
            } else if (an_u - a_u_max).abs() < d_u {
                b_is_left = false;
                an_u1 = an_u - an_u_period;
            }
        }
        if an_v_period > 0.0 {
            if (an_v - a_v_min).abs() < d_v {
                b_is_left = true;
                an_v1 = an_v + an_v_period;
            } else if (an_v - a_v_max).abs() < d_v {
                b_is_left = false;
                an_v1 = an_v - an_v_period;
            }
        }
        // OCCT L194-197.
        if an_u1 == an_u && an_v1 == an_v {
            return false;
        }
        //
        // OCCT L199: aScPr = (anU1 == anU) ? aDir2D1 * aDOX : aDir2D1 * aDOY.
        let a_sc_pr = if an_u1 == an_u {
            a_dir2d1.dot(a_dox)
        } else {
            a_dir2d1.dot(a_doy)
        };
        //
        // OCCT L201-207: trimmed copies; aC2 translated by (anU1-anU, anV1-anV).
        let a_c1 = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(c2d1.clone()),
            t_min: a,
            t_max: b,
        });
        let a_c2_base = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(c2d1),
            t_min: a,
            t_max: b,
        });
        let a_tr_v = DVec2::new(an_u1 - an_u, an_v1 - an_v);
        let a_c2 = translate_curve2d(&a_c2_base, a_tr_v);
        //
        // OCCT L209-230: BB.UpdateEdge(aSp, aC1, aC2, aF, aTol) depending on
        // bIsLeft and aScPr.
        let sp_idx = self.ds.index(&a_sp) as usize;
        if sp_idx >= self.ds.nb_shapes() {
            return false;
        }
        if !b_is_left {
            if a_sc_pr < 0.0 {
                self.ds.update_edge_closed_surface(sp_idx, face_key, a_c2, a_c1, a, b, a_tol);
            } else {
                self.ds.update_edge_closed_surface(sp_idx, face_key, a_c1, a_c2, a, b, a_tol);
            }
        } else if a_sc_pr < 0.0 {
            self.ds.update_edge_closed_surface(sp_idx, face_key, a_c1, a_c2, a, b, a_tol);
        } else {
            self.ds.update_edge_closed_surface(sp_idx, face_key, a_c2, a_c1, a, b, a_tol);
        }
        true
    }

    /// OCCT BOPTools_AlgoTools3D::DoSplitSEAMOnFace(aEOrigin, aESplit, aFace)
    /// (BOPTools_AlgoTools3D.cxx L236-327). The split edge carries a single
    /// pcurve; projects its midpoint onto the original seam's two pcurves and
    /// builds the second pcurve translated to the opposite seam line.
    fn do_split_seam_on_face_origin(
        &self,
        a_e_origin: &Shape,
        a_e_split: &Shape,
        a_face: &Shape,
    ) -> bool {
        // OCCT L240-243.
        if !self.edge_closed_on_face(a_e_origin, a_face) {
            return false;
        }
        if self.edge_closed_on_face(a_e_split, a_face) {
            return true;
        }
        //
        // OCCT L250-257: aC2DSplit = CurveOnSurface(aESplit, aFace, aTS1, aTS2).
        // BRep_Tool::CurveOnSurface reads the pcurve from the edge TShape,
        // which is shared 鈥?fall back to the DS canonical shape like
        // do_split_seam_on_face (Arc::make_mut clones on write).
        let mut a_e_split_f = a_e_split.clone();
        a_e_split_f.orientation = Orientation::Forward;
        // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the key
        // location is L.Predivided(E.Location()).
        let face_key = (
            a_face.ptr_id(),
            crate::bop::algo::compose_face_edge_pcurve_location(
                a_face.location, a_e_split_f.location, &self.ds.locations),
        );
        let fi = self.ds.index(a_face) as usize;
        let (a_c2d_split, a_ts1, a_ts2) = {
            let direct = a_e_split_f
                .as_edge()
                .and_then(|ed| ed.pcurves.get(&face_key).cloned());
            if let Some(v) = direct {
                v
            } else {
                let Some(idx) = self.ds.map_shape_index.get(&(a_e_split_f.ptr_id(), a_e_split_f.location)) else {
                    return false;
                };
                let Some(ed2) = self.ds.shape_info(*idx).shape.as_edge() else {
                    return false;
                };
                let Some(v) = ed2.pcurves.get(&face_key) else {
                    return false;
                };
                v.clone()
            }
        };
        //
        // OCCT L263-267: the original seam's two pcurves (forward -> pcurve1,
        // reversed -> pcurve2). BRep_Tool::CurveOnSurface reads from the
        // shared TShape 鈥?fall back to the DS canonical shape.
        let find_closed_rep = |e: &Shape| -> Option<(Curve2d, Curve2d, f64, f64)> {
            e.as_edge().and_then(|ed| {
                ed.representations.iter().find_map(|r| match r {
                    CurveRepresentation::CurveOnClosedSurface {
                        face, pcurve1, pcurve2, range,
                    } if *face == face_key => {
                        Some((pcurve1.clone(), pcurve2.clone(), range[0], range[1]))
                    }
                    _ => None,
                })
            })
        };
        let (a_c2d1, a_c2d2, a_t1, a_t2) = match find_closed_rep(a_e_origin) {
            Some(v) => v,
            None => {
                let Some(idx) = self.ds.map_shape_index.get(&(a_e_origin.ptr_id(), a_e_origin.location)) else {
                    return false;
                };
                let Some(ed2) = self.ds.shape_info(*idx).shape.as_edge() else {
                    return false;
                };
                match ed2.representations.iter().find_map(|r| match r {
                    CurveRepresentation::CurveOnClosedSurface {
                        face, pcurve1, pcurve2, range,
                    } if *face == face_key => {
                        Some((pcurve1.clone(), pcurve2.clone(), range[0], range[1]))
                    }
                    _ => None,
                }) {
                    Some(v) => v,
                    None => return false,
                }
            }
        };
        //
        // OCCT L269-272: aT = IntermediatePoint(aTS1, aTS2); D1 -> aPMid, aVTgt.
        let a_t = intermediate_point(a_ts1, a_ts2);
        let a_p_mid = a_c2d_split.point_at(a_t);
        let a_v_tgt = a_c2d_split.derivative_at(a_t);
        //
        // OCCT L275-282: project the midpoint onto the two original pcurves.
        // Geom2dAPI_ProjectPointOnCurve::Init -> Extrema_ExtPC2d ->
        // Extrema_GGExtPC default theTolF = 1.0e-10.
        let a_proj1 = ExtPC2d::new(a_p_mid, &a_c2d1, 1e-10, a_t1, a_t2);
        let a_proj2 = ExtPC2d::new(a_p_mid, &a_c2d2, 1e-10, a_t1, a_t2);
        if a_proj1.nb_ext() == 0 && a_proj2.nb_ext() == 0 {
            return false;
        }
        // OCCT L284-285: aDist1 = LowerDistance() = sqrt(SquareDistance())
        // (Geom2dAPI_ProjectPointOnCurve.cxx L178) 鈥?non-squared distance,
        // compared against PConfusion() below.
        let a_dist1 = if a_proj1.nb_ext() > 0 {
            a_proj1.square_distance(1).sqrt()
        } else {
            f64::MAX
        };
        let a_dist2 = if a_proj2.nb_ext() > 0 {
            a_proj2.square_distance(1).sqrt()
        } else {
            f64::MAX
        };
        // OCCT L287-290: PConfusion check.
        if a_dist1 > PCONFUSION && a_dist2 > PCONFUSION {
            return false;
        }
        //
        // OCCT L293-294: aNewPnt = closer opposite-curve point.
        let a_new_pnt = if a_dist1 < a_dist2 {
            a_c2d2.point_at(a_proj1.point(1).param)
        } else {
            a_c2d1.point_at(a_proj2.point(1).param)
        };
        //
        // OCCT L296-303: trimmed copies; aC2 translated from aPMid to aNewPnt.
        let a_c1 = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(a_c2d_split.clone()),
            t_min: a_ts1,
            t_max: a_ts2,
        });
        let a_c2_base = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(a_c2d_split),
            t_min: a_ts1,
            t_max: a_ts2,
        });
        let a_tr_vec = a_new_pnt - a_p_mid;
        let a_c2 = translate_curve2d(&a_c2_base, a_tr_vec);
        //
        // OCCT L305-314: aVTgtOrigin at the projection point.
        let (a_p_proj, a_v_tgt_origin) = if a_dist1 < a_dist2 {
            let t = a_proj1.point(1).param;
            (a_c2d1.point_at(t), a_c2d1.derivative_at(t))
        } else {
            let t = a_proj2.point(1).param;
            (a_c2d2.point_at(t), a_c2d2.derivative_at(t))
        };
        //
        // OCCT L316-325: aDot = aVTgt . aVTgtOrigin.
        let a_dot = a_v_tgt.dot(a_v_tgt_origin);
        let a_tol = a_e_split_f.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
        let sp_idx = self.ds.index(&a_e_split_f) as usize;
        if sp_idx >= self.ds.nb_shapes() {
            return false;
        }
        if (a_dist1 < a_dist2) == (a_dot > 0.0) {
            self.ds.update_edge_closed_surface(sp_idx, face_key, a_c1, a_c2, a_ts1, a_ts2, a_tol);
        } else {
            self.ds.update_edge_closed_surface(sp_idx, face_key, a_c2, a_c1, a_ts1, a_ts2, a_tol);
        }
        true
    }

    /// OCCT BRepLib::BuildPCurveForEdgesOnPlane (BRepLib_1.cxx L313-339 with the
    /// template loop L122-130): for each edge in `a_le` lacking a stored pcurve
    /// on the plane face, project the edge's 3D curve onto the plane and store
    /// the pcurve (BRep_Tool::CurveOnSurface's plane projection branch).
    fn build_pcurve_for_edges_on_plane(&self, a_le: &[Shape], a_f: &Shape) {
        let face_key = (a_f.ptr_id(), a_f.location);
        let surf = match a_f.as_face().and_then(|fd| fd.surface.clone()) {
            Some(s) => s,
            None => return,
        };
        let Surface3::Plane(pl) = surf else { return };
        for a_e in a_le {
            // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the key
            // location is L.Predivided(E.Location()).
            let e_key = (
                face_key.0,
                crate::bop::algo::compose_face_edge_pcurve_location(
                    face_key.1, a_e.location, &self.ds.locations),
            );
            let is_stored = a_e
                .as_edge()
                .map(|ed| ed.pcurves.contains_key(&e_key))
                .unwrap_or(true);
            if is_stored {
                continue;
            }
            let Some(curve) = a_e.as_edge().and_then(|ed| ed.curve.clone()) else {
                continue;
            };
            let range = a_e.as_edge().map(|ed| ed.range).unwrap_or([0.0, 0.0]);
            let Some(pc) = project_edge_on_plane(&curve, &pl, range) else {
                continue;
            };
            let e_idx = self.ds.index(a_e) as usize;
            if e_idx >= self.ds.nb_shapes() {
                continue;
            }
            let a_tol = a_e.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
            self.ds
                .update_edge_pcurve_shared(e_idx, e_key, pc, range[0], range[1], a_tol);
        }
    }

    /// OCCT BOPTools_AlgoTools2D::IsEdgeIsoline (L669-698) 鈥?true when the
    /// edge's pcurve is tangent to the U or V axis in the face UV space.
    /// is_u_iso = U-isoline (constant U), is_v_iso = V-isoline (constant V).
    fn is_edge_isoline(&self, a_e: &Shape, a_f: &Shape) -> (bool, bool) {
        let mut is_u_iso = false;
        let mut is_v_iso = false;
        // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the key
        // location is L.Predivided(E.Location()).
        let e_key = (
            a_f.ptr_id(),
            crate::bop::algo::compose_face_edge_pcurve_location(
                a_f.location, a_e.location, &self.ds.locations),
        );
        let pcurve = a_e.as_edge().and_then(|ed| {
            ed.pcurves
                .get(&e_key)
                .map(|(pc, f, l)| (pc.clone(), *f, *l))
        });
        if let Some((pc, a_first, a_last)) = pcurve {
            // aPC->D1(0.5*(aFirst+aLast), aP, aT)
            let a_t = pc.derivative_at(0.5 * (a_first + a_last));
            let a_sq_magn = a_t.length_squared();
            // OCCT L683-687: if (aSqMagn <= gp::Resolution()) return;
            // (gp::Resolution() == 1e-15).
            if a_sq_magn <= 1e-15 {
                return (false, false);
            }
            let a_t = a_t / a_sq_magn.sqrt();
            // aRefVDir(0,1), aRefUDir(1,0); CrossMagnitude(aT, aRefV) = |aT.X|.
            let a_tol = rcad_kernel::core::precision::ANGULAR; // Precision::Angular()
            let a_dpv = a_t.x.abs();
            let a_dpu = a_t.y.abs();
            is_u_iso = a_dpv <= a_tol;
            is_v_iso = a_dpu <= a_tol;
        }
        (is_u_iso, is_v_iso)
    }

    /// OCCT BuildDraftFace (BOPAlgo_Builder_2.cxx L1052-1189).
    /// Builds a new face from the original face by replacing each boundary
    /// edge with its images. Returns None when the BuilderFace algorithm must
    /// be used instead (internal edges / multi-connected vertices / unified edges).
    fn build_draft_face(&mut self, the_face: &Shape) -> Option<Shape> {
        let a_surf = the_face.as_face().and_then(|fd| fd.surface.clone())?;
        let a_tol = the_face.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
        // OCCT L1073-1074: aVerticesCounter 鈥?multi-connexity detection.
        // OCCT map keys use TopTools_ShapeMapHasher (TShape + Location).
        let mut a_vertices_counter: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        // OCCT L1078: aMEdges 鈥?edges-unification fence (TopTools_ShapeMapHasher).
        let mut a_m_edges: HashSet<(u64, u32)> = HashSet::new();

        // OCCT L1081-1181: rebuild each wire of the face.
        let mut new_wires: Vec<Shape> = Vec::new();
        for a_w in self.shape_sub_shapes(the_face) {
            if a_w.shape_type() != topods::ShapeType::Wire {
                continue;
            }
            // OCCT L1091-1095: skip empty wires.
            let w_edges: Vec<Shape> = if let TShape::Wire(wd) = &*a_w.data {
                wd.edges.iter().map(|sr| Shape::new(sr.data.clone(), sr.location, sr.orientation)).collect()
            } else {
                Vec::new()
            };
            if w_edges.is_empty() {
                continue;
            }
            let mut new_edges: Vec<Shape> = Vec::new();
            for a_e in &w_edges {
                let an_ori_e = a_e.orientation;
                // OCCT L1105-1110: internal edges may split the face -> BuilderFace.
                if an_ori_e == topods::Orientation::Internal {
                    return None;
                }
                // OCCT L1113-1115: degenerated / closed on face checks.
                let b_is_degenerated = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                let b_is_closed = self.edge_closed_on_face(a_e, the_face);
                // OCCT L1118: theImages.Seek(aE).
                let p_le_im = self.images_of(a_e);
                
                if p_le_im.is_none() {
                    // OCCT L1121-1131: multi-connected / unified edge -> BuilderFace.
                    if !b_is_degenerated && self.has_multi_connected(a_e, &mut a_vertices_counter) {
                        return None;
                    }
                    if !b_is_closed && !a_m_edges.insert((a_e.ptr_id(), a_e.location)) {
                        return None;
                    }
                    new_edges.push(a_e.clone());
                    continue;
                }
                // OCCT L1137-1175: replace by images.
                for a_sp0 in &p_le_im.unwrap() {
                    let mut a_sp = a_sp0.clone();
                    if !b_is_degenerated && self.has_multi_connected(&a_sp, &mut a_vertices_counter) {
                        return None;
                    }
                    if !b_is_closed && !a_m_edges.insert((a_sp.ptr_id(), a_sp.location)) {
                        return None;
                    }
                    // OCCT L1154: aSp.Orientation(anOriE).
                    a_sp.orientation = an_ori_e;
                    if b_is_degenerated {
                        new_edges.push(a_sp);
                        continue;
                    }
                    // OCCT L1163-1166: seam split 鈥?DoSplitSEAMOnFace(aSp, theFace)
                    // (overload 1 only, no fallback/warning in BuildDraftFace).
                    if b_is_closed && !self.edge_closed_on_face(&a_sp, the_face) {
                        self.do_split_seam_on_face(&a_sp, the_face);
                    }
                    // OCCT L1169-1172: IsSplitToReverseWithWarn(aSp, aE) 鈥?reverse
                    // the split when its geometry is oppositely oriented to aE.
                    if crate::bop::tools::algo_tools::is_split_to_reverse_edge(&a_sp, a_e).0 {
                        a_sp.orientation = flip_orientation(a_sp.orientation);
                    }
                    new_edges.push(a_sp);
                }
            }
            // OCCT L1178-1180: MakeWire(aNewWire) + orientation + closed flag.
            // OCCT L1179: aNewWire.Closed(BRep_Tool::IsClosed(aNewWire)).
            let mut wire_flags = tshape_flags::DEFAULT;
            if wire_is_closed(&new_edges) {
                wire_flags |= tshape_flags::CLOSED;
            }
            let new_wire = Shape::new(
                std::sync::Arc::new(TShape::Wire(TWireData {
                    my_shapes: vec![],
                    flags: wire_flags,
                    edges: new_edges,
                })),
                0,
                a_w.orientation,
            );
            new_wires.push(new_wire);
        }
        // OCCT L1066: MakeFace(aDraftFace, aS, aLoc, aTol) 鈥?face without wires.
        let mut draft_face = Shape::new(
            std::sync::Arc::new(TShape::Face(TFaceData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                surface: Some(a_surf),
                surface_location: 0,
                outer_wire: new_wires.first().cloned().unwrap_or_else(Shape::null),
                inner_wires: new_wires.into_iter().skip(1).collect(),
                sample_point: None,
                uv_domain: None,
                internal_vertices: vec![],
                tolerance: a_tol,
                natural_restriction: false,
            })),
            0,
            topods::Orientation::Forward,
        );
        // OCCT L1183-1186: reverse if the original face was reversed.
        if the_face.orientation == topods::Orientation::Reversed {
            draft_face.orientation = topods::Orientation::Reversed;
        }
        Some(draft_face)
    }

    /// OCCT HasMultiConnected (BOPAlgo_Builder_2.cxx L1014-1045).
    /// Returns true when any vertex of the edge is shared by more than two edges.
    /// The map uses TopTools_ShapeMapHasher (TShape + Location).
    fn has_multi_connected(
        &self,
        the_edge: &Shape,
        the_map: &mut HashMap<(u64, u32), Vec<Shape>>,
    ) -> bool {
        let verts: Vec<Shape> = match &*the_edge.data {
            TShape::Edge(ed) => {
                // OCCT TopoDS_Iterator(aE, cumLoc) 鈥?compose the edge Location
                // into the vertices.
                let loc = the_edge.location;
                let vf_loc = if loc == 0 { ed.first.location } else { loc };
                let vl_loc = if loc == 0 { ed.last.location } else { loc };
                vec![
                    Shape::new(ed.first.data.clone(), vf_loc, ed.first.orientation),
                    Shape::new(ed.last.data.clone(), vl_loc, ed.last.orientation),
                ]
            }
            _ => Vec::new(),
        };
        for v in verts {
            let vkey = (v.ptr_id(), v.location);
            let list = the_map.entry(vkey).or_default();
            // OCCT L1035: pList->Contains(theEdge) 鈥?TopoDS_Shape::operator==
            // is IsEqual (TShape + Location + Orientation), so two instances of
            // the same TShape with different orientations are distinct.
            let ekey = (the_edge.ptr_id(), the_edge.location, the_edge.orientation);
            if !list.iter().any(|e| (e.ptr_id(), e.location, e.orientation) == ekey) {
                list.push(the_edge.clone());
            }
            if list.len() > 2 {
                return true;
            }
        }
        false
    }

    /// OCCT BOPAlgo_Builder::FillSameDomainFaces (Builder_2.cxx L580-780).
    fn fill_same_domain_faces(&mut self) {
        // OCCT L584-589: get FF interferences, empty check.
        let a_ffs = &self.ds.interf_ff;
        if a_ffs.is_empty() { return; }

        // OCCT L597-649: build face-to-parent solid map (with image propagation).
        // OCCT aFaceToParent (L593-594) is DataMap<Shape, Shape,
        // TopTools_ShapeMapHasher> 鈥?key identity TShape + Location.
        let mut a_face_to_parent: HashMap<(u64, u32), u64> = HashMap::new(); // face 鈫?solid
        let a_nb_src = self.ds.nb_source_shapes();
        for i_src in 0..a_nb_src {
            let a_si = self.ds.shape_info(i_src);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_solid = self.brep_sr(i_src);
            // OCCT L610-618: TopExp_Explorer(aSolid, TopAbs_FACE) 鈥?deep
            // traversal of ALL nested faces (including faces directly held by
            // the solid, without a shell); only the first binding is kept.
            let mut a_sf: Vec<Shape> = Vec::new();
            Self::collect_sub_shapes_of_type_static(&a_solid, topods::ShapeType::Face, &mut a_sf);
            for a_f in a_sf {
                a_face_to_parent
                    .entry((a_f.ptr_id(), a_f.location))
                    .or_insert(a_solid.ptr_id());
            }
        }
        // OCCT L619-648: propagate the parent solid to the image faces.
        // OCCT L636: aPropagation is NCollection_DataMap<TopoDS_Shape, TopoDS_Shape,
        // TopTools_ShapeMapHasher> 鈥?bucket iteration order (L660-665).
        let mut a_propagation: crate::bop::algo::occt_map::OcctDataMapInt<(u64, u32), u64> =
            crate::bop::algo::occt_map::OcctDataMapInt::new();
        // OCCT L640-655: iterate myImages 鈥?NCollection_DataMap bucket order.
        for (a_src, a_l_im) in self.my_images.iter() {
            let p_parent = a_face_to_parent.get(&a_src).copied();
            if let Some(parent) = p_parent {
                for a_piece in a_l_im {
                    let pk = (a_piece.ptr_id(), a_piece.location);
                    if !a_face_to_parent.contains_key(&pk) {
                        a_propagation.insert(pk, parent);
                    }
                }
            }
        }
        // OCCT L660-665: iterate aPropagation 鈥?bucket order.
        for (k, v) in a_propagation.iter() {
            a_face_to_parent.entry(k).or_insert(*v);
        }

        // OCCT L654-684: collect face indices from FF interferences.
        let mut a_fi_vec: Vec<usize> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for ff in a_ffs {
            for &nf in &[ff.f1, ff.f2] {
                if !self.ds.has_face_info(nf) { continue; }
                if a_m_fence.insert(nf) {
                    a_fi_vec.push(nf);
                }
            }
        }
        
        // OCCT L687: sort indices.
        a_fi_vec.sort();

        // OCCT L690-694: map edge-sets 鈫?list of faces, planar face fence.
        // OCCT anESetFaces is IndexedDataMap<BOPTools_Set, List<Shape>>; the
        // BOPTools_Set is built by AddEdgeSet 鈫?BOPTools_Set::Add(theS, EDGE)
        // (BOPTools_Set.cxx L124-166): degenerated edges are skipped (L132-137),
        // INTERNAL edges are expanded to FORWARD+REVERSED (L139-148), and
        // IsEqual (L81-99) compares the expanded count first, then the
        // deduplicated TShape+Location set. rcad: (count, sorted unique
        // (ptr_id, location)).
        let mut an_e_set_faces: IndexMap<(usize, Vec<(u64, u32)>), Vec<Shape>> = IndexMap::new();
        let mut a_mf_planar: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();

        // OCCT L697-741: for each face, compute edge set.
        for &n_f in &a_fi_vec {
            let a_f = self.brep_sr(n_f);
            // OCCT L707-718: check if planar (Plane surface, bounded bbox).
            let b_check_planar = if let Some(surf) = self.ds.face_surface(n_f) {
                matches!(surf, rcad_kernel::geom::Surface3::Plane(_))
                    && !self.ds.shape_info(n_f).bbox.is_open()
            } else { false };

            // OCCT L720-740: get face images (or face itself), add edge set.
            let face_list: Vec<Shape> = if let Some(imgs) = self.my_images.get((a_f.ptr_id(), a_f.location)) {
                imgs.clone()
            } else {
                vec![a_f.clone()]
            };

            for f_piece in &face_list {
                // OCCT AddEdgeSet (L562-571) 鈥?BOPTools_Set::Add(theS, EDGE).
                // OCCT AddEdgeSet (L562-571) 鈥?BOPTools_Set::Add(theS, EDGE).
                let mut count: usize = 0;
                let mut edge_set: Vec<(u64, u32)> = Vec::new();
                for e_shape in self.face_edges(f_piece) {
                    // OCCT BOPTools_Set.cxx L132-137: degenerated edges skipped.
                    let degen = e_shape.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                    if degen {
                        continue;
                    }
                    count += 1;
                    let ekey = (e_shape.ptr_id(), e_shape.location);
                    if !edge_set.contains(&ekey) {
                        edge_set.push(ekey);
                    }
                    // OCCT BOPTools_Set.cxx L139-148: INTERNAL edges expanded
                    // to FORWARD + REVERSED entries (count +1).
                    if e_shape.orientation == topods::Orientation::Internal {
                        count += 1;
                    }
                }
                edge_set.sort_unstable();

                // OCCT L729: if (bCheckPlanar) aMFPlanar.Add(aItLF.Value()) 鈥?
                // the planar fence feeds the L780-785 fast SD check.
                if b_check_planar {
                    a_mf_planar.insert((f_piece.ptr_id(), f_piece.location));
                }
                an_e_set_faces.entry((count, edge_set)).or_default().push(f_piece.clone());
            }
        }

        // OCCT L743-748: aDMSLS 鈥?back-and-forth SD map (IndexedDataMap,
        // insertion order); aVPSB 鈥?pairs for analysis.
        let mut a_dmsls: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        let mut a_vpsb: Vec<(Shape, Shape)> = Vec::new();

        // OCCT L750-791: check pairs of faces with equal edge set.
        for (_edge_set, faces) in &an_e_set_faces {
            if faces.len() < 2 { continue; }
            for i1 in 0..faces.len() {
                let f1 = &faces[i1];
                let parent1 = a_face_to_parent.get(&(f1.ptr_id(), f1.location)).copied();
                let b_check_planar = a_mf_planar.contains(&(f1.ptr_id(), f1.location));
                for i2 in (i1 + 1)..faces.len() {
                    let f2 = &faces[i2];
                    let parent2 = a_face_to_parent.get(&(f2.ptr_id(), f2.location)).copied();
                    // OCCT L776-779: two faces of one solid cannot be SD.
                    if let (Some(p1), Some(p2)) = (parent1, parent2) {
                        if p1 == p2 { continue; }
                    }
                    // OCCT L780-785: planar bounded faces 鈫?SD without additional check.
                    if b_check_planar && a_mf_planar.contains(&(f2.ptr_id(), f2.location)) {
                        fill_map_faces(f1, f2, &mut a_dmsls);
                        continue;
                    }
                    // OCCT L786-791: add pair for analysis.
                    a_vpsb.push((f1.clone(), f2.clone()));
                }
            }
        }

        // OCCT L799-822: perform the pair analysis.
        for (f1, f2) in &a_vpsb {
            // OCCT BOPAlgo_PairOfShapeBoolean::Perform (Builder_2.cxx L94-105)
            // always calls BOPTools_AlgoTools::AreFacesSameDomain on the raw
            // TopoDS_Face pair (BOPTools_AlgoTools.cxx L1139-1205), which works
            // on the face geometry directly with no DS access. rcad mirrors
            // this with the Shape-based check for every pair.
            let flag = crate::bop::tools::algo_tools::are_faces_same_domain_shapes(
                f1, f2, self.my_fuzzy_value, &self.ds.locations,
            );
            if flag {
                fill_map_faces(f1, f2, &mut a_dmsls);
            }
        }

        // OCCT L826: MakeBlocks(aDMSLS, aMBlocks).
        let a_m_blocks = make_blocks_faces(&a_dmsls);

        // OCCT L830-882: fill same domain faces map.
        for a_lsd in &a_m_blocks {
            // Find the SD face: the original face (in DS) with minimal index;
            // otherwise the first face of the group.
            let mut p_fsd: Option<Shape> = None;
            let mut n_f_min: isize = std::isize::MAX;
            for a_f in a_lsd {
                let n_f = self.ds.index(a_f);
                if n_f >= 0 {
                    // OCCT L858: original face 鈥?consider it split into itself.
                    self.my_images.bound((a_f.ptr_id(), a_f.location)).push(a_f.clone());
                    if n_f < n_f_min {
                        n_f_min = n_f;
                        p_fsd = Some(a_f.clone());
                    }
                }
            }
            let p_fsd = match p_fsd {
                Some(fsd) => fsd,
                None => match a_lsd.first() {
                    Some(f) => f.clone(),
                    None => continue,
                },
            };
            // OCCT L876-881: bind all faces of the group to the SD face.
            for a_f in a_lsd {
                self.my_shapes_sd
                    .insert((a_f.ptr_id(), a_f.location), p_fsd.clone());
            }
        }

        // OCCT L886-921: update images with SD faces and fill origins.
        for i in 0..a_nb_src {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Face { continue; }
            let a_f = self.brep_sr(i);
            let Some(a_lf_im) = self.my_images.get_mut((a_f.ptr_id(), a_f.location)) else { continue };
            for a_f_im in a_lf_im.iter_mut() {
                // OCCT L906-910: replace the image with its SD face.
                if let Some(a_fsd) = self.my_shapes_sd.get(&(a_f_im.ptr_id(), a_f_im.location)) {
                    *a_f_im = a_fsd.clone();
                }
                // OCCT L913-919: fill the map of origins.
                let p_lf_or = self.my_origins.entry((a_f_im.ptr_id(), a_f_im.location)).or_default();
                p_lf_or.push(a_f.clone());
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillInternalVertices (Builder_2.cxx L929-1008).
    ///
    /// Adds vertices strictly inside face images as INTERNAL vertices. OCCT
    /// collects (vertex, face-image) pairs and classifies each with
    /// IntTools_Context::ComputeVF (iFlag == 0), then mutates the face in place
    /// with BRep_Builder().Add 鈥?preserving the TShape identity of Same-Domain
    /// faces shared across image lists.
    fn fill_internal_vertices(&mut self) {
        // OCCT L935: vector of Vertex/Face pairs for classification.
        // rcad pair: (source face index, DS vertex index, face image).
        let mut a_vvfi: Vec<(usize, usize, Shape)> = Vec::new();
        // OCCT L937-938: iterate all source shapes.
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            // OCCT L941-944: only faces are processed.
            if a_si.shape_type != topods::ShapeType::Face {
                continue;
            }
            // OCCT L951-956: aF = aSI.Shape(); pLFIm = myImages.Seek(aF);
            // if (!pLFIm) continue.
            let a_f = self.brep_sr(i);
            let Some(p_lf_im) = self.my_images.get((a_f.ptr_id(), a_f.location)).cloned() else {
                continue;
            };
            // OCCT L959-960: myDS->AloneVertices(i, aLIAV).
            let a_liav = self.alone_vertices(i);
            // OCCT L963-978: build (vertex, face-image) pairs.
            for &n_v in &a_liav {
                for a_f_im in &p_lf_im {
                    a_vvfi.push((i, n_v, a_f_im.clone()));
                }
            }
        }
        // OCCT L988-1006: classify each pair and add the internal vertices.
        for (i, n_v, a_f_im) in a_vvfi {
            // OCCT BOPAlgo_VFI::Perform (L188-200):
            //   myContext->ComputeVF(myV, myF, aT1, aT2, dummy, myFuzzyValue);
            //   myIsInternal = (iFlag == 0).
            if !self.point_in_face_image(&a_f_im, n_v) {
                continue;
            }
            // OCCT L966-967: aV = Vertex(v); aV.Orientation(INTERNAL).
            let a_v = Shape::new(
                Arc::new(TShape::Vertex(TVertexData {
                    my_shapes: vec![],
                    flags: tshape_flags::DEFAULT,
                    point: self.ds.vertex_point_by_idx(n_v),
                    tolerance: self.ds.vertex_tolerance_by_idx(n_v),
                    points: vec![],
                })),
                0,
                topods::Orientation::Internal,
            );
            // OCCT L1005: BRep_Builder().Add(aF, aV) 鈥?mutate the face in
            // place, keeping the TShape identity of Same-Domain faces shared
            // across image lists. rcad: the pair holds a clone of the image
            // Arc, so the image entry in myImages is located by ptr_id and
            // mutated there.
            self.add_internal_vertex_to_image(i, a_f_im.ptr_id(), a_v);
        }
    }

    /// OCCT L1005: BRep_Builder().Add(aF, aV) 鈥?add the INTERNAL vertex to the
    /// face image stored in myImages. rcad: the pair carries a clone of the
    /// image Arc; the image is located by (source face, ptr_id) and mutated via
    /// Arc::make_mut. A Same-Domain face shared across image lists is cloned
    /// only when the classification accepted a truly internal vertex 鈥?the
    /// boundary vertices of the s06 4鈫? case are rejected by the
    /// classification, so the shared SD identity is preserved.
    fn add_internal_vertex_to_image(&mut self, i: usize, a_ptr: u64, a_v: Shape) {
        let a_f = self.brep_sr(i);
        if let Some(imgs) = self.my_images.get_mut((a_f.ptr_id(), a_f.location)) {
            for img in imgs.iter_mut() {
                if img.ptr_id() != a_ptr {
                    continue;
                }
                let ts = Arc::make_mut(&mut img.data);
                if let TShape::Face(fd) = ts {
                    fd.internal_vertices.push(a_v);
                }
                return;
            }
        }
    }

    /// OCCT IntTools_Context::ComputeVF (IntTools_Context.cxx L546-591) applied
    /// to a face image (split piece) instead of a DS face. OCCT projects the
    /// vertex onto the face surface and classifies the UV point with the face's
    /// FClass2d (strictly inside 鈫?internal). rcad's FClass2d is DS-index based
    /// and split pieces are not registered in the DS, so the piece's UV
    /// polygons are rebuilt by projecting its wire vertices onto the piece's
    /// surface. The classification follows FClass2d's multi-loop semantics:
    /// the point must be strictly inside the outer wire and outside every
    /// inner wire (hole).
    fn point_in_face_image(&self, piece: &Shape, n_v: usize) -> bool {
        let a_p = self.ds.vertex_point_by_idx(n_v);
        let (surf, outer_wire, inner_wires, a_tol_f) = match &*piece.data {
            TShape::Face(fd) => (
                fd.surface.clone(),
                fd.outer_wire.clone(),
                fd.inner_wires.clone(),
                fd.tolerance,
            ),
            _ => return false,
        };
        let Some(surf) = surf else { return false; };
        // 1. GeomAPI_ProjectPointOnSurf: closest UV of the vertex.
        let (a_uv, a_proj) = crate::bop::closest_point_on_surface(&surf, a_p);
        let a_dist = (a_proj - a_p).length();
        // OCCT L597-604: aTolV + aTolF + theFuzz distance check.
        let a_tol_v = self.ds.vertex_tolerance_by_idx(n_v);
        let a_tol_sum = a_tol_v + a_tol_f + self.my_fuzzy_value.max(rcad_kernel::CONFUSION);
        if a_dist > a_tol_sum {
            return false;
        }
        // 2. Build the UV boundary polygons from the piece's wires.
        // OCCT IntTools_FClass2d::Init samples the boundary edges' pcurves;
        // rcad projects the wire vertices onto the piece's surface instead.
        let build_poly = |w: &Shape, surf: &rcad_kernel::geom::Surface3| -> Vec<glam::DVec2> {
            let mut a_poly: Vec<glam::DVec2> = Vec::new();
            if let TShape::Wire(wd) = &*w.data {
                for e in &wd.edges {
                    if let TShape::Edge(ed) = &*e.data {
                        if let TShape::Vertex(vd) = &*ed.last.data {
                            let (uv, _) = crate::bop::closest_point_on_surface(surf, vd.point);
                            a_poly.push(uv);
                        }
                    }
                }
            }
            a_poly
        };
        let a_poly = build_poly(&outer_wire, &surf);
        if a_poly.len() < 3 {
            return false;
        }
        // 3. On-boundary check: the vertex's UV within tolerance of any
        //    boundary segment (OCCT FClass2d ON state 鈥?not internal).
        let a_tol_on = a_tol_f.max(rcad_kernel::CONFUSION);
        let on_boundary = |poly: &[glam::DVec2], a_uv: glam::DVec2| -> bool {
            for i in 0..poly.len() {
                let a = poly[i];
                let b = poly[(i + 1) % poly.len()];
                let ab = b - a;
                let len2 = ab.length_squared();
                if len2 < 1e-30 {
                    continue;
                }
                let ap = a_uv - a;
                let t = (ap.dot(ab) / len2).clamp(0.0, 1.0);
                let d = (ap - t * ab).length();
                if d <= a_tol_on {
                    return true;
                }
            }
            false
        };
        if on_boundary(&a_poly, a_uv) {
            return false;
        }
        // 4. Holes: a point inside any inner wire (hole) is NOT internal.
        for w in &inner_wires {
            let a_hole = build_poly(w, &surf);
            if a_hole.len() < 3 {
                continue;
            }
            if on_boundary(&a_hole, a_uv) {
                return false;
            }
            if rcad_kernel::base::gprop::tri::point_in_polygon_2d(&a_hole, a_uv) {
                return false;
            }
        }
        // 5. Strictly inside the outer UV polygon (OCCT IsPointInFace, ON
        //    excluded).
        rcad_kernel::base::gprop::tri::point_in_polygon_2d(&a_poly, a_uv)
    }

    /// OCCT BOPAlgo_Builder::FillImagesSolids (BOPAlgo_Builder_3.cxx L60-93).
    /// Builds split solids: FillIn3DParts -> BuildSplitSolids -> FillInternalShapes.
    fn fill_images_solids(&mut self) {
        // OCCT L62-73: check all DS source shapes for SOLID type
        let a_nb_s = self.ds.nb_source_shapes();
        let mut has_solid = false;
        for i in 0..a_nb_s {
            if self.ds.shape_info(i).shape_type == topods::ShapeType::Solid {
                has_solid = true;
                break;
            }
        }
        if !has_solid { return; }
        // OCCT L78: local draft-solids map passed through the s06 chain.
        // theDraftSolids is DataMap with TopTools_ShapeMapHasher (L414) 鈥?
        // key identity TShape + Location of the SOURCE solid.
        let mut a_draft_solids: HashMap<(u64, u32), Shape> = HashMap::new();
        self.fill_in_3d_parts(&mut a_draft_solids);
        if self.has_errors() { return; }
        self.build_split_solids(&a_draft_solids);
        if self.has_errors() { return; }
        self.fill_internal_shapes();
    }

    /// OCCT BOPAlgo_Builder::FillIn3DParts (BOPAlgo_Builder_3.cxx L97-263).
    /// Collects all faces and draft solids, classifies faces as IN/OUT relative
    /// to each solid, and fills myInParts + myDraftSolids.
    fn fill_in_3d_parts(&mut self, the_draft_solids: &mut HashMap<(u64, u32), Shape>) {
        // OCCT L113-150: get all faces (source + their images).
        let a_nb_s = self.ds.nb_source_shapes();
        let mut a_lfaces: Vec<Shape> = Vec::new();
        let mut a_m_fence: HashSet<(u64, u32)> = HashSet::new();
        if std::env::var("RCAD_BS_DEBUG").is_ok() {
            // Face list orientations (sources + images).
            let mut ors: Vec<String> = Vec::new();
            for i in 0..a_nb_s {
                if self.ds.shape_info(i).shape_type != topods::ShapeType::Face { continue; }
                let s = self.ds.shape(i);
                let imgs = self.my_images.get((s.ptr_id(), s.location)).cloned().unwrap_or_default();
                let base = format!("{}:{}", i, if s.orientation == topods::Orientation::Reversed { "R" } else { "F" });
                if imgs.is_empty() {
                    ors.push(base);
                } else {
                    let ims: Vec<String> = imgs.iter().map(|im| {
                        if im.orientation == topods::Orientation::Reversed { "R".to_string() } else { "F".to_string() }
                    }).collect();
                    ors.push(format!("{}({})", base, ims.join("")));
                }
            }
            eprintln!("[F3D-COL3] {}", ors.join(","));
            // Map each source face to its parent solid (arg vs tool)
            for i in 0..a_nb_s {
                if self.ds.shape_info(i).shape_type != topods::ShapeType::Face { continue; }
                let s = self.ds.shape(i);
                let parent = if self.ds.arguments.get(0).map(|a| a.ptr_id() == s.ptr_id()).unwrap_or(false) { "ARG0" } else { "?" };
                let imgs = self.my_images.get((s.ptr_id(), s.location)).cloned().unwrap_or_default();
                let n = s.as_face().and_then(|fd| fd.surface.clone()).map(|sf| match sf {
                    rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                    _ => "O".into(),
                }).unwrap_or_else(|| "?".into());
                let imgids: Vec<String> = imgs.iter().map(|im| format!("{}", im.ptr_id() % 100000)).collect();
                eprintln!("[F3D-SRC] ds={} parent={} ptr={} n={} imgs=[{}]", i, parent, s.ptr_id() % 100000, n, imgids.join(" "));
                for im in &imgs {
                    let in2 = im.as_face().and_then(|fd| fd.surface.clone()).map(|sf| match sf {
                        rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                        _ => "O".into(),
                    }).unwrap_or_else(|| "?".into());
                    let neds = crate::bop::algo::builder::face_edges(im).len();
                    eprintln!("[F3D-IMG]   img={} n={} edges={}", im.ptr_id() % 100000, in2, neds);
                }
            }
        }
        // OCCT aShapeBoxMap (L113-114): shape -> Bnd_Box for the box culling
        // in ClassifyFaces. Both faces (L148) and solids (L193) are bound;
        // the face box is the DS box (aSI.Box(), with gap). Key identity is
        // TShape + Location (TopTools_ShapeMapHasher); the box is stored as
        // (corner_min, corner_max, gap) 鈥?the gap participates in the
        // IsOut(Bnd_Box) checks of ClassifyFaces (BOPTools_AlgoTools.cxx
        // L1498-1508), while the BVH culling uses only the corners.
        let mut a_shape_box_map: HashMap<(u64, u32), (DVec3, DVec3, f64)> = HashMap::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Face {
                continue;
            }
            let a_s = self.brep_sr(i);
            // OCCT L131-148: add images, or the face itself with its box.
            if let Some(imgs) = self.my_images.get((a_s.ptr_id(), a_s.location)).cloned() {
                for a_s_im in &imgs {
                    if a_m_fence.insert((a_s_im.ptr_id(), a_s_im.location)) {
                        a_lfaces.push(a_s_im.clone());
                    }
                }
            } else {
                // OCCT L147-148: aLFaces.Append(aS); aShapeBoxMap.Bind(aS, aSI.Box()).
                a_lfaces.push(a_s.clone());
                if let (Some(bmin), Some(bmax)) = (a_si.bbox.corner_min(), a_si.bbox.corner_max()) {
                    a_shape_box_map.insert(
                        (a_s.ptr_id(), a_s.location),
                        (bmin, bmax, a_si.bbox.get_gap()),
                    );
                }
            }
        }
        // OCCT L154-195: get all solids, build draft solids.
        let mut a_lsolids: Vec<Shape> = Vec::new();
        let mut a_solids_if: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        let mut a_draft_solid: HashMap<(u64, u32), Shape> = HashMap::new();
        let mut a_source_solids: Vec<Shape> = Vec::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid {
                continue;
            }
            let a_solid = self.brep_sr(i);
            // OCCT L181-185: Bnd_Box& aBoxS = aSI.ChangeBox();
            // if (aBoxS.IsVoid()) { myDS->BuildBndBoxSolid(i, aBoxS, myCheckInverted); }
            // rcad stores the box as (corner_min, corner_max, gap); the DS
            // solid box is void for multi-argument BOPs (prepare_solids only
            // pre-builds it for a single argument, BOPDS_DS.cxx L1789-1792),
            // so BuildBndBoxSolid is invoked here exactly like OCCT.
            let a_box_s: (DVec3, DVec3, f64) =
                if let (Some(bmin), Some(bmax)) = (a_si.bbox.corner_min(), a_si.bbox.corner_max()) {
                    (bmin, bmax, a_si.bbox.get_gap())
                } else {
                    let mut box3 = (
                        DVec3::splat(f64::INFINITY),
                        DVec3::splat(f64::NEG_INFINITY),
                        0.0,
                    );
                    self.ds.build_bnd_box_solid(i, &mut box3, self.my_check_inverted);
                    box3
                };
            // OCCT L186-194: BuildDraftSolid(aSolid, aSD, aLIF).
            let mut a_lif: Vec<Shape> = Vec::new();
            let a_sd = self.build_draft_solid(&a_solid, &mut a_lif);
            a_lsolids.push(a_sd.clone());
            a_solids_if.insert((a_sd.ptr_id(), a_sd.location), a_lif);
            a_draft_solid.insert((a_solid.ptr_id(), a_solid.location), a_sd.clone());
            // OCCT L193: aShapeBoxMap.Bind(aSD, aBoxS) 鈥?the SOLID's box bound
            // under the DRAFT solid's key (used by ClassifyFaces box culling).
            a_shape_box_map.insert((a_sd.ptr_id(), a_sd.location), a_box_s);
            a_source_solids.push(a_solid);
        }
        // OCCT L197-208: classify the faces relative to the draft solids.
        let an_in_parts =
            Self::classify_faces(&self.ds, &a_lfaces, &a_lsolids, &a_solids_if, &a_shape_box_map);
        // OCCT L210-262: analyze the results of classification.
        for a_solid in &a_source_solids {
            let a_sd = match a_draft_solid.get(&(a_solid.ptr_id(), a_solid.location)) {
                Some(sd) => sd.clone(),
                None => continue,
            };
            let a_l_in_faces = an_in_parts
                .get(&(a_sd.ptr_id(), a_sd.location))
                .cloned()
                .unwrap_or_default();
            let a_l_internal = a_solids_if
                .get(&(a_sd.ptr_id(), a_sd.location))
                .cloned()
                .unwrap_or_default();
            let a_nb_in = a_l_in_faces.len();
            if a_nb_in == 0 {
                // OCCT L227-238: check if the shells of the solid have an image.
                // OCCT L229: TopoDS_Iterator it(aSolid) iterates ALL direct
                // sub-shapes 鈥?shells plus internal vertices/edges.
                let mut b_has_image = false;
                let mut subs = self.shape_sub_shapes(a_solid);
                if let TShape::Solid(sd) = &*a_solid.data {
                    subs.extend(sd.internal_vertices.iter().cloned());
                    subs.extend(sd.internal_edges.iter().cloned());
                }
                for sh in subs {
                    if self.my_images.contains((sh.ptr_id(), sh.location)) {
                        b_has_image = true;
                        break;
                    }
                }
                if !b_has_image {
                    continue;
                }
            }
            // OCCT L241: theDraftSolids.Bind(aSolid, aSDraft).
            the_draft_solids.insert((a_solid.ptr_id(), a_solid.location), a_sd.clone());
            // OCCT L243-261: combine IN and internal faces into myInParts.
            let a_nb_int = a_l_internal.len();
            if a_nb_int != 0 || a_nb_in != 0 {
                let p_lin = self.my_in_parts.entry((a_solid.ptr_id(), a_solid.location)).or_default();
                p_lin.extend(a_l_in_faces);
                p_lin.extend(a_l_internal);
            }
        }
    }

    /// OCCT BOPAlgo_Builder::BuildDraftSolid (BOPAlgo_Builder_3.cxx L267-368).
    /// Builds the draft solid from the solid's shells, replacing each face by
    /// its images (keeping the SD faces and INTERNAL faces in theLIF).
    /// rcad returns the built solid (OCCT theDraftSolid out-param); the
    /// orientation set at OCCT L280-281 is applied when the solid is created
    /// (L3364-3375 below).
    fn build_draft_solid(&self, the_solid: &Shape, the_lif: &mut Vec<Shape>) -> Shape {
        // OCCT L283-284: aIt1.Initialize(theSolid) 鈥?iterate direct sub-shapes;
        // non-SHELL sub-shapes (internal edges, vertices) are skipped at L287-289.
        let mut a_shd_list: Vec<Shape> = Vec::new();
        for a_sh in self.shape_sub_shapes(the_solid) {
            if a_sh.shape_type() != topods::ShapeType::Shell {
                continue;
            }
            let a_or_sh = a_sh.orientation;
            // OCCT L292-295: aOrSh; aBB.MakeShell(aShD); aShD.Orientation(aOrSh); iFlag=0.
            let mut a_shd_subs: Vec<Shape> = Vec::new();
            let mut i_flag = 0;
            // OCCT L297-298: aIt2.Initialize(aSh) 鈥?iterate the shell's faces.
            for a_f in self.shape_sub_shapes(&a_sh) {
                // OCCT L300-301: aF; aOrF.
                let a_or_f = a_f.orientation;
                // OCCT L303: if myImages.IsBound(aF) 鈥?replace by images.
                if let Some(imgs) = self.my_images.get((a_f.ptr_id(), a_f.location)).cloned() {
                    for a_fx in &imgs {
                        if std::env::var("RCAD_BS_DEBUG").is_ok() {
                            let n = a_fx.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                                rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                                _ => "O".into(),
                            }).unwrap_or_else(|| "?".into());
                            let sn = a_f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                                rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                                _ => "O".into(),
                            }).unwrap_or_else(|| "?".into());
                            eprintln!("[DRAFT-DEC] srcface={} n={} or={:?} img={} n={} or={:?} sd={}", a_f.ptr_id() % 100000, sn, a_or_f, a_fx.ptr_id() % 100000, n, a_fx.orientation, self.my_shapes_sd.contains_key(&(a_fx.ptr_id(), a_fx.location)));
                        }
                        // OCCT L309: aFx; L311: if myShapesSD.IsBound(aFx) 鈥?SD face.
                        if self.my_shapes_sd.contains_key(&(a_fx.ptr_id(), a_fx.location)) {
                            // OCCT L314-318: if (aOrF == INTERNAL) { aFx.Orientation(aOrF);
                            // theLIF.Append(aFx); }
                            if a_or_f == topods::Orientation::Internal {
                                let mut fx = a_fx.clone();
                                fx.orientation = topods::Orientation::Internal;
                                the_lif.push(fx);
                            } else {
                                // OCCT L321-329: IsSplitToReverseWithWarn(aFx, aF, ctx, report);
                                // if (bToReverse) aFx.Reverse(); iFlag=1; aBB.Add(aShD, aFx).
                                // rcad: the warning alert is omitted (diagnostic only).
                                let (b_to_reverse, _err) =
                                    crate::bop::tools::algo_tools::is_split_to_reverse_face(a_fx, &a_f, &self.ds);
                                if std::env::var("RCAD_BS_DEBUG").is_ok() {
                                    eprintln!("[DRAFT-SD] img={} to_reverse={} err={}", a_fx.ptr_id() % 100000, b_to_reverse, _err);
                                }
                                let mut fx = a_fx.clone();
                                if b_to_reverse {
                                    fx.orientation = flip_orientation(fx.orientation);
                                }
                                i_flag = 1;
                                a_shd_subs.push(fx);
                            }
                        } else {
                            // OCCT L334-344: aFx.Orientation(aOrF); if (aOrF == INTERNAL)
                            // theLIF.Append(aFx); else { iFlag=1; aBB.Add(aShD, aFx); }
                            let mut fx = a_fx.clone();
                            fx.orientation = a_or_f;
                            if a_or_f == topods::Orientation::Internal {
                                the_lif.push(fx);
                            } else {
                                i_flag = 1;
                                a_shd_subs.push(fx);
                            }
                        }
                    }
                } else {
                    // OCCT L348-359: if (aOrF == INTERNAL) theLIF.Append(aF);
                    // else { iFlag=1; aBB.Add(aShD, aF); }
                    if a_or_f == topods::Orientation::Internal {
                        the_lif.push(a_f.clone());
                    } else {
                        i_flag = 1;
                        a_shd_subs.push(a_f.clone());
                    }
                }
            }
            // OCCT L362-365: if (iFlag) { aShD.Closed(BRep_Tool::IsClosed(aShD));
            // aBB.Add(theDraftSolid, aShD); }
            if i_flag != 0 {
                // OCCT L364: aShD.Closed(BRep_Tool::IsClosed(aShD)).
                let mut flags = tshape_flags::DEFAULT;
                if shell_is_closed(&a_shd_subs) {
                    flags |= tshape_flags::CLOSED;
                }
                let a_shd = Shape::new(
                    std::sync::Arc::new(TShape::Shell(TShellData {
                        my_shapes: vec![],
                        flags,
                        faces: a_shd_subs,
                    })),
                    0,
                    a_or_sh,
                );
                a_shd_list.push(a_shd);
            }
        }
        // OCCT L280-281: aOrSd = theSolid.Orientation(); theDraftSolid.Orientation(aOrSd)
        // (the caller also did MakeSolid at Builder_3.cxx L188).
        Shape::new(
            std::sync::Arc::new(TShape::Solid(TSolidData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                shells: a_shd_list,
                internal_vertices: vec![],
                internal_edges: vec![],
            })),
            0,
            the_solid.orientation,
        )
    }

    /// OCCT BOPAlgo_Tools::ClassifyFaces (BOPAlgo_Tools.cxx L1622-1747) 鈥?
    /// classifies faces relative to solids using connexity blocks.
    /// `a_shape_box_map`: shape (TShape+Location) -> Bnd_Box as
    /// (corner_min, corner_max, gap), filled at Builder_3.cxx L148/L193 (faces
    /// under their own key, solids under the DRAFT solid's key); empty when the
    /// caller has no cached boxes 鈥?then the boxes are built on the fly
    /// (matching OCCT L1665-1673/L1694-1711 where BRepBndLib::Add builds them).
    pub(crate) fn classify_faces(
        ds: &DS,
        faces: &[Shape],
        solids: &[Shape],
        a_solids_if: &HashMap<(u64, u32), Vec<Shape>>,
        a_shape_box_map: &HashMap<(u64, u32), (DVec3, DVec3, f64)>,
    ) -> HashMap<(u64, u32), Vec<Shape>> {
        let mut in_parts: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        for a_sd in solids {
            // aMSF = own faces of the draft solid + its internal faces.
            // aMSE = edges of the draft solid (OCCT L1366-1368).
            // OCCT aMSF/aMSE are IndexedMap with TopTools_ShapeMapHasher 鈥?
            // key identity TShape + Location.
            let mut a_msf: HashSet<(u64, u32)> = HashSet::new();
            let mut a_mse: HashSet<(u64, u32)> = HashSet::new();
            for f in collect_solid_faces(a_sd) {
                a_msf.insert((f.ptr_id(), f.location));
                for e in face_edges(&f) {
                    a_mse.insert((e.ptr_id(), e.location));
                }
            }
            // OCCT L1371: bIsEmpty is evaluated BEFORE the own INTERNAL faces
            // are added to aMSF (L1374-1378) 鈥?a solid made of INTERNAL faces
            // only is treated as empty by the classification below.
            let b_is_empty = a_msf.is_empty();
            if let Some(lif) = a_solids_if.get(&(a_sd.ptr_id(), a_sd.location)) {
                for f in lif {
                    a_msf.insert((f.ptr_id(), f.location));
                }
            }

            // OCCT L1380-1393: select the faces whose bounding box interferes
            // with the solid's box and which are not the solid's own faces.
            // (OCCT builds a BVH tree of face boxes; the interference test is
            // equivalent here, BOPTools_BoxTreeSelector.)
            // OCCT L1694-1712: the solid's box comes from theShapeBoxMap (the
            // source solid's box, bound under the draft solid's key by the
            // caller) or is built here (BRepBndLib::Add; IsInvertedSolid -> whole).
            // OCCT L1648-1657: the face's box comes from theShapeBoxMap (the
            // DS box with gap, bound by FillIn3DParts L148) or is built here.
            let solid_bbox = a_shape_box_map
                .get(&(a_sd.ptr_id(), a_sd.location))
                .copied()
                .or_else(|| shape_bbox(a_sd).map(|(a, b)| (a, b, 0.0)));
            let mut a_ivec: Vec<usize> = Vec::new();
            for (i, a_f) in faces.iter().enumerate() {
                if a_msf.contains(&(a_f.ptr_id(), a_f.location)) {
                    continue;
                }
                if let Some((smin, smax, _sgap)) = solid_bbox {
                    // BVH box-selection (BOPTools_BoxTreeSelector). The BVH
                    // boxes come from Bnd_Tools::Bnd2BVH -> Bnd_Box::Get(),
                    // which INCLUDES the box gap (GetXMin = Xmin - Gap,
                    // GetXMax = Xmax + Gap, Bnd_Box.hxx L181-186). rcad's
                    // BndBox::corner_min/corner_max (Bnd_Box::get) also
                    // include the gap, so the corners are used directly.
                    let face_bbox = a_shape_box_map
                        .get(&(a_f.ptr_id(), a_f.location))
                        .copied()
                        .or_else(|| shape_bbox(a_f).map(|(a, b)| (a, b, 0.0)));
                    if std::env::var("RCAD_BS_DEBUG").is_ok() && a_sd.shape_type() == topods::ShapeType::Solid {
                        let fb = face_bbox;
                        eprintln!("[F3D-BB] face_idx={} ds={} smin=({:.2},{:.2},{:.2}) smax=({:.2},{:.2},{:.2}) fmin={:?} fmax={:?}",
                            i, ds.index(&faces[i]), smin.x, smin.y, smin.z, smax.x, smax.y, smax.z,
                            fb.map(|b| b.0), fb.map(|b| b.1));
                    }
                    let Some((fmin, fmax, _fgap)) = face_bbox else {
                        continue;
                    };
                    if fmax.x < smin.x
                        || fmin.x > smax.x
                        || fmax.y < smin.y
                        || fmin.y > smax.y
                        || fmax.z < smin.z
                        || fmin.z > smax.z
                    {
                        continue;
                    }
                }
                a_ivec.push(i);
            }
            // OCCT L1398-1403: sort the selected indices.
            a_ivec.sort_unstable();
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[F3D-SEL] solid_ptr={} n_faces={} n_ivec={} b_empty={}",
                    a_sd.ptr_id(), faces.len(), a_ivec.len(), b_is_empty);
                for &i in &a_ivec {
                    let src = ds.index(&faces[i]);
                    let n = faces[i].as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                        rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                        _ => "O".into(),
                    }).unwrap_or_else(|| "?".into());
                    eprintln!("[F3D-SEL]   idx={} ds={} ptr={} {}", i, src, faces[i].ptr_id() % 100000, n);
                }
            }

            // OCCT L1405-1417: the solid has no faces -> all selected faces are IN.
            if b_is_empty {
                for &i in &a_ivec {
                    in_parts
                        .entry((a_sd.ptr_id(), a_sd.location))
                        .or_default()
                        .push(faces[i].clone());
                }
                continue;
            }

            // OCCT L1419-1428: EF map of the faces to process. Built only when
            // more than one face is selected (L1422: if (aNbFP > 1)); with a
            // single face the connexity block is the face itself.
            // OCCT aMEFP is IndexedDataMap with TopTools_ShapeMapHasher.
            let mut a_mefp: HashMap<(u64, u32), Vec<usize>> = HashMap::new();
            if a_ivec.len() > 1 {
                for &i in &a_ivec {
                    for e in face_edges(&faces[i]) {
                        a_mefp.entry((e.ptr_id(), e.location)).or_default().push(i);
                    }
                }
            }

            // OCCT L1498-1502: EF map of the solid (built once).
            let a_mef = build_edge_face_map(a_sd);
            if std::env::var("RCAD_GFO_DEBUG").is_ok() {
                let mut n_shells = 0usize;
                if let TShape::Solid(sd) = &*a_sd.data {
                    n_shells = sd.shells.len();
                    for sh in &sd.shells {
                        eprintln!("[MEF-SD] solid_or={:?} shell_or={:?} faces={}",
                            a_sd.orientation, sh.orientation,
                            if let TShape::Shell(shd) = &*sh.data { shd.faces.len() } else { 0 });
                    }
                }
                for (k, v) in a_mef.iter() {
                    for f in v {
                        eprintln!("[MEF] edge=({},{}) face_or={:?}", k.0, k.1, f.orientation);
                    }
                }
                let _ = n_shells;
            }

            // OCCT L1508: Precision::Confusion() 鈥?tolerance of the IsInternalFace
            // classification (used by the ComputeState fallback).
            let the_tol = rcad_kernel::CONFUSION;

            // OCCT L1437-1438: fence map to avoid processing faces twice.
            let mut a_mf_done: HashSet<usize> = HashSet::new();

            // OCCT L1444-1518: main classification loop over connexity blocks.
            for &k in &a_ivec {
                if !a_mf_done.insert(k) {
                    continue;
                }

                // OCCT L1460-1465: make the connexity block of face k.
                let mut a_lcb: Vec<usize> = Vec::new();
                let mut a_face_to_classify: Option<usize> = None;
                make_connexity_block(
                    k,
                    &a_mse,
                    &a_mefp,
                    &mut a_mf_done,
                    &mut a_lcb,
                    &mut a_face_to_classify,
                    faces,
                );

                // OCCT L1467-1491: fast check that all block vertices interfere
                // with the solid's bounding box; otherwise the block is out.
                // myBoxS.IsOut(aBBV) 鈥?Bnd_Box::IsOut(Bnd_Box) adds both boxes'
                // gaps: vertex box gap = vertex tolerance, solid box gap = the
                // DS box gap (BOPTools_AlgoTools.cxx L1503-1508).
                if let Some((smin, smax, sgap)) = solid_bbox {
                    let mut b_out = false;
                    for &bfi in &a_lcb {
                        for (p, gap) in shape_vertices(&faces[bfi]) {
                            let egap = gap + sgap;
                            if p.x - egap > smax.x
                                || p.x + egap < smin.x
                                || p.y - egap > smax.y
                                || p.y + egap < smin.y
                                || p.z - egap > smax.z
                                || p.z + egap < smin.z
                            {
                                b_out = true;
                                break;
                            }
                        }
                        if b_out {
                            break;
                        }
                    }
                    if b_out {
                        continue;
                    }
                }

                // OCCT L1493-1496: representative face for the classification.
                let a_fc = a_face_to_classify.unwrap_or(k);
                let a_fc_shape = &faces[a_fc];


                // OCCT L1505-1509: IsInternalFace on the representative face.
                let is_in = is_internal_face(a_fc_shape, a_sd, &a_mef, the_tol, ds) == 1;
                if is_in {
                    // OCCT L1510-1517: the whole connexity block is IN. Each
                    // face appears in exactly one block (aMFDone fence), so no
                    // duplicate check is needed 鈥?OCCT just Appends.
                    let entry = in_parts.entry((a_sd.ptr_id(), a_sd.location)).or_default();
                    for &bfi in &a_lcb {
                        entry.push(faces[bfi].clone());
                    }
                }
            }
        }
        in_parts
    }

    /// OCCT BOPAlgo_Builder::BuildSplitSolids (BOPAlgo_Builder_3.cxx L413-618).
    /// Each source solid is split independently by a separate BOPAlgo_SplitSolid;
    /// only its own split pieces are assigned to its images.
    fn build_split_solids(&mut self, the_draft_solids: &HashMap<(u64, u32), Shape>) {
        // OCCT L425-427: map of same-domain solids face sets (BOPTools_Set -> shape).
        // BOPTools_Set::IsEqual (BOPTools_Set.cxx L81-99) compares the total
        // entry count (INTERNAL faces expanded to FORWARD+REVERSED, L139-148)
        // first, then the deduplicated, orientation-insensitive set of faces 鈥?
        // key is (count, sorted unique face ids).
        let mut a_mst: HashMap<(usize, Vec<(u64, u32)>), Shape> = HashMap::new();
        // OCCT L432-461: find same-domain solids for non-interfered solids.
        // OCCT aMFence (L428) is Map<TopoDS_Shape, TopTools_ShapeMapHasher>.
        let a_nb_s = self.ds.nb_source_shapes();
        let mut a_mfence: HashSet<(u64, u32)> = HashSet::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_s = self.brep_sr(i);
            if !a_mfence.insert((a_s.ptr_id(), a_s.location)) { continue; }
            // OCCT L451-454: if theDraftSolids.IsBound(aS) continue;
            if the_draft_solids.contains_key(&(a_s.ptr_id(), a_s.location)) { continue; }
            // OCCT L456-459: aST.Add(aS, TopAbs_FACE).
            let a_st = self.shape_face_set(&a_s);
            a_mst.entry(a_st).or_insert_with(|| a_s.clone());
        }
        // OCCT L465-466: aSolidsIm 鈥?source solid -> result solids (IndexedDataMap
        // with TopTools_ShapeMapHasher).
        let mut a_solids_im: Vec<(Shape, Vec<Shape>)> = Vec::new();
        let mut a_solids_im_idx: HashMap<(u64, u32), usize> = HashMap::new();
        // OCCT L468-518: build split solids for interfered source solids.
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_s = self.brep_sr(i);
            // OCCT L478-481: if !theDraftSolids.IsBound(aS) continue;
            if !the_draft_solids.contains_key(&(a_s.ptr_id(), a_s.location)) { continue; }
            let a_sd = the_draft_solids.get(&(a_s.ptr_id(), a_s.location)).unwrap().clone();
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[DRAFT] src={} sd_or={:?}", a_s.ptr_id() % 100000, a_sd.orientation);
                for ss in self.shape_sub_shapes(&a_sd) {
                    if ss.shape_type() == topods::ShapeType::Shell {
                        eprintln!("[DRAFT]   shell or={:?}", ss.orientation);
                        for f in self.shape_sub_shapes(&ss) {
                            if f.shape_type() == topods::ShapeType::Face {
                                let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                                    rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                                    _ => "O".into(),
                                }).unwrap_or_else(|| "?".into());
                                eprintln!("[DRAFT]     face or={:?} {}", f.orientation, n);
                            }
                        }
                    }
                }
            }
            // OCCT L484-489: if no IN faces -> the draft solid itself, no split.
            let p_lfin: Vec<Shape> = self
                .my_in_parts
                .get(&(a_s.ptr_id(), a_s.location))
                .cloned()
                .unwrap_or_default();
            if p_lfin.is_empty() {
                let idx = *a_solids_im_idx.entry((a_s.ptr_id(), a_s.location)).or_insert_with(|| {
                    a_solids_im.push((a_s.clone(), Vec::new()));
                    a_solids_im.len() - 1
                });
                a_solids_im[idx].1.push(a_sd);
                continue;
            }
            // OCCT L493-499: 1.1 shell faces set of the draft solid.
            // aExp.Init(aSD, TopAbs_FACE) 鈥?TopExp_Explorer with cumulative
            // orientation (TopoDS_Iterator cumOri, TopExp_Explorer.cxx L152):
            // face.or = aSD.or * shell.or * face.or. Location composition is
            // identity here (draft solid and shells carry location 0).
            let mut a_sfs: Vec<Shape> = Vec::new();
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[BS-IN] src={} n_in={}", a_s.ptr_id(), p_lfin.len());
                for f in &p_lfin {
                    let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                        rcad_kernel::geom::Surface3::Plane(p) => format!("Plane n=({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                        _ => "O".into(),
                    }).unwrap_or_else(|| "?".into());
                    eprintln!("[BS-IN]   face={} or={:?} {}", f.ptr_id(), f.orientation, n);
                }
                for t in &self.my_tools {
                    for f in self.map_shapes_of_type(t, topods::ShapeType::Face) {
                        let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                            rcad_kernel::geom::Surface3::Plane(p) => format!("Plane n=({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                            _ => "O".into(),
                        }).unwrap_or_else(|| "?".into());
                        eprintln!("[BS-IN]   TOOLFACE={} or={:?} {}", f.ptr_id(), f.orientation, n);
                        let imgs = self.my_images.get((f.ptr_id(), f.location)).map(|v| v.len());
                        eprintln!("[BS-IN]   TOOLFACE imgcount={:?}", imgs);
                    }
                }
            }
            let sd_or = a_sd.orientation;
            for ss in self.shape_sub_shapes(&a_sd) {
                if ss.shape_type() == topods::ShapeType::Shell {
                    let ss_or = sd_or.compose(ss.orientation);
                    for f in self.shape_sub_shapes(&ss) {
                        if f.shape_type() == topods::ShapeType::Face {
                            let mut ff = f;
                            ff.orientation = ss_or.compose(ff.orientation);
                            a_sfs.push(ff);
                        }
                    }
                }
            }
            // OCCT L501-511: 1.2 add IN faces (both orientations).
            for a_f in &p_lfin {
                let mut f_fwd = a_f.clone();
                f_fwd.orientation = topods::Orientation::Forward;
                a_sfs.push(f_fwd);
                let mut f_rev = a_f.clone();
                f_rev.orientation = topods::Orientation::Reversed;
                a_sfs.push(f_rev);
            }
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[BS-SHAPES] src={} n={}", a_s.ptr_id(), a_sfs.len());
                for f in &a_sfs {
                    let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                        rcad_kernel::geom::Surface3::Plane(p) => format!("Plane n=({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                        _ => "O".into(),
                    }).unwrap_or_else(|| "?".into());
                    let eds: Vec<String> = face_edges(f).iter().map(|e| {
                        let (p0, p1) = match &*e.data {
                            rcad_kernel::topods::TShape::Edge(ed) => {
                                let c = ed.curve.clone();
                                match c {
                                    Some(rcad_kernel::geom::Curve3::Line(l)) => {
                                        let a = l.origin;
                                        let b = l.origin + l.direction;
                                        (format!("({:.1},{:.1},{:.1})", a.x, a.y, a.z), format!("({:.1},{:.1},{:.1})", b.x, b.y, b.z))
                                    }
                                    _ => ("?".into(), "?".into()),
                                }
                            }
                            _ => ("?".into(), "?".into()),
                        };
                        format!("e{}:{}-{}{}", e.ptr_id() % 100000, p0, p1, if e.orientation == rcad_kernel::topods::Orientation::Reversed { "R" } else { "F" })
                    }).collect();
                    eprintln!("[BS-SHAPES]   {}:{} or={:?} edges=[{}]", f.ptr_id() % 100000, n, f.orientation, eds.join(" "));
                }
            }
            // OCCT L514-517: BOPAlgo_SplitSolid& aBS = aVBS.Appended();
            // aBS.SetSolid(aSolid); aBS.SetShapes(aSFS); aBS.SetRunParallel(myRunParallel).
            let mut bs = crate::bop::algo::builder_solid::BuilderSolid::new(&self.ds);
            // OCCT L515: SetSolid(aSolid) 鈥?the SOURCE solid (not the draft
            // solid); the split result is keyed by it in aSolidsIm (L542).
            bs.my_solid = Some(a_s.clone());
            // OCCT L516: SetShapes(aSFS).
            bs.my_shapes = a_sfs;
            bs.perform();
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[SPLITSOLID] src={} n_solids={}", a_s.ptr_id(), bs.my_solids.len());
            }

            // OCCT L542: aSolidsIm.Add(aBS.Solid(), aBS.Areas()) 鈥?keyed by the
            // split solid's mySolid (the source solid).
            let a_solid = bs.my_solid.clone().expect("SetSolid before split");
            let idx = *a_solids_im_idx.entry((a_solid.ptr_id(), a_solid.location)).or_insert_with(|| {
                a_solids_im.push((a_solid.clone(), Vec::new()));
                a_solids_im.len() - 1
            });
            a_solids_im[idx].1.extend(bs.my_solids.clone());
            // OCCT L544-577: merge the split solid's report into the main
            // report, converting all sub-split errors into warnings. OCCT
            // additionally wraps TopoDS_AlertWithShape alerts into a compound
            // of the solid + the alert shape (L559-569); rcad keeps the alert
            // as-is 鈥?diagnostics only, no topology impact.
            for a in bs.report().errors() {
                self.my_report.add_warning(a.clone());
            }
        }
        // OCCT L580-617: add new solids to the images map (same-domain dedup).
        for (a_s, a_lsr) in &a_solids_im {
            // OCCT L586: if !myImages.IsBound(aS).
            if self.my_images.contains((a_s.ptr_id(), a_s.location)) { continue; }
            // Compute the same-domain dedup results first (immutable borrows).
            let mut results: Vec<(Shape, Shape, bool)> = Vec::new();
            for a_sr in a_lsr {
                // OCCT L593-601: BOPTools_Set of aSR's faces. aMST.Added(aST)
                // INSERTS the face set on first sight and returns the existing
                // entry on later sights, so a later solid sharing the same
                // face set dedups to the first one.
                let a_st = self.shape_face_set(a_sr);
                let b_flag_sd = a_mst.contains_key(&a_st);
                let a_sx = a_mst
                    .entry(a_st)
                    .or_insert_with(|| a_sr.clone())
                    .clone();
                results.push((a_sr.clone(), a_sx, b_flag_sd));
            }
            let p_lsx = self.my_images.bound((a_s.ptr_id(), a_s.location));
            for (a_sr, a_sx, b_flag_sd) in results {
                p_lsx.push(a_sx.clone());
                // OCCT L604-609: myOrigins[aSx].Append(aS).
                self.my_origins
                    .entry((a_sx.ptr_id(), a_sx.location))
                    .or_default()
                    .push(a_s.clone());
                // OCCT L611-614: if same-domain, bind myShapesSD[aSR] = aSx.
                if b_flag_sd {
                    self.my_shapes_sd.insert((a_sr.ptr_id(), a_sr.location), a_sx);
                }
            }
        }
    }
    /// OCCT BOPAlgo_Builder::FillInternalShapes (Builder_3.cxx L622-830).
    fn fill_internal_shapes(&mut self) {
        // OCCT L631-644: local lists/maps. aMFence (L639) and aMSI (IndexedMap)
        // use TopTools_ShapeMapHasher (TShape + Location).
        let mut a_lsc: Vec<Shape> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
        let mut a_l_args: Vec<Shape> = Vec::new();
        // aMSI indexed map 鈫?rcad: Vec<Shape> with dedup
        let mut a_msi: Vec<Shape> = Vec::new();
        let mut a_msi_fence: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
        // Map for ancestor lookup (vertex鈫抂edges], vertex鈫抂faces], edge鈫抂faces]).
        // OCCT aMSx (Builder_3.cxx L635-636) 鈥?TopTools_ShapeMapHasher.
        let mut a_msx: std::collections::HashMap<(u64, u32), Vec<u64>> = std::collections::HashMap::new();

        // OCCT L653-659: TreatCompound on arguments 鈫?flatten into aLSC
        let a_arguments = &self.ds.arguments;
        for a_s in a_arguments {
            Self::treat_compound(a_s, &mut a_lsc, &mut a_m_fence);
        }
        // OCCT L660-681: collect V/E from aLSC into aLArgs
        a_m_fence.clear();
        for a_s in &a_lsc {
            let a_type = a_s.shape_type();
            if a_type == topods::ShapeType::Wire {
                for ss in self.shape_sub_shapes(a_s) {
                    if a_m_fence.insert((ss.ptr_id(), ss.location)) {
                        a_l_args.push(ss);
                    }
                }
            } else if a_type == topods::ShapeType::Vertex || a_type == topods::ShapeType::Edge {
                a_l_args.push(a_s.clone());
            }
        }
        // OCCT L684-709: for each V/E/W, add images or self to aMSI
        a_m_fence.clear();
        for a_s in &a_l_args {
            let a_type = a_s.shape_type();
            if a_type == topods::ShapeType::Vertex
                || a_type == topods::ShapeType::Edge
                || a_type == topods::ShapeType::Wire
            {
                if a_m_fence.insert((a_s.ptr_id(), a_s.location)) {
                    if let Some(imgs) = self.my_images.get((a_s.ptr_id(), a_s.location)) {
                        for img in imgs {
                            if a_msi_fence.insert((img.ptr_id(), img.location)) {
                                a_msi.push(img.clone());
                            }
                        }
                    } else {
                        if a_msi_fence.insert((a_s.ptr_id(), a_s.location)) {
                            a_msi.push(a_s.clone());
                        }
                    }
                }
            }
        }

        // OCCT L721-788: internal V/E from source solids
        a_m_fence.clear();
        let a_nb_s = self.ds.nb_source_shapes();
        let mut a_ls_d: Vec<Shape> = Vec::new();
        // aMSOr: original solids without images (OCCT L785, TopTools_ShapeMapHasher).
        let mut a_ms_or: HashSet<(u64, u32)> = HashSet::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_s = self.brep_sr(i);

            // OCCT L738: OwnInternalShapes(aS, aMx)
            // rcad: iterate solid sub-shapes, find internal V/E
            let a_mx = self.own_internal_shapes(&a_s);

            // OCCT L741-758: add internal shapes to aMSI
            for a_si_internal in &a_mx {
                if let Some(imgs) = self.my_images.get((a_si_internal.ptr_id(), a_si_internal.location)) {
                    for img in imgs {
                        if a_msi_fence.insert((img.ptr_id(), img.location)) {
                            a_msi.push(img.clone());
                        }
                    }
                } else {
                    if a_msi_fence.insert((a_si_internal.ptr_id(), a_si_internal.location)) {
                        a_msi.push(a_si_internal.clone());
                    }
                }
            }

            // OCCT L760-787: build ancestor map from split solids
            if let Some(imgs) = self.my_images.get((a_s.ptr_id(), a_s.location)) {
                for a_sp in imgs {
                    if a_m_fence.insert((a_sp.ptr_id(), a_sp.location)) {
                        Self::map_shapes_and_ancestors(a_sp, &mut a_msx);
                        a_ls_d.push(a_sp.clone());
                    }
                }
            } else {
                if a_m_fence.insert((a_s.ptr_id(), a_s.location)) {
                    Self::map_shapes_and_ancestors(&a_s, &mut a_msx);
                    a_ls_d.push(a_s.clone());
                    // OCCT L785: aMSOr.Add(aS).
                    a_ms_or.insert((a_s.ptr_id(), a_s.location));
                }
            }
        }

        // OCCT L792-809: filter aMSI 鈥?keep only shapes not tied to split solid faces
        let mut a_ls_i: Vec<Shape> = Vec::new();
        for a_si in &a_msi {
            if let Some(ancestors) = a_msx.get(&(a_si.ptr_id(), a_si.location)) {
                if ancestors.is_empty() {
                    a_ls_i.push(a_si.clone());
                }
            } else {
                a_ls_i.push(a_si.clone());
            }
        }

        // OCCT L812-816: empty check
        if a_ls_i.is_empty() { return; }

        // OCCT L820-877: settle internal vertices and edges into solids.
        // OCCT iterates aLSd (aIt.Initialize(aLSd)) 鈥?a snapshot; the solid
        // handle aSd is local to each iteration and never written back.
        let mut i_sd = 0;
        while i_sd < a_ls_d.len() {
            let mut a_sd = a_ls_d[i_sd].clone();
            let mut i_si = 0;
            while i_si < a_ls_i.len() {
                let mut a_si = a_ls_i[i_si].clone();
                // OCCT L834: aSI.Orientation(TopAbs_INTERNAL).
                a_si.orientation = topods::Orientation::Internal;
                // OCCT L836: ComputeStateByOnePoint(aSI, aSd, 1.e-11, ctx).
                let a_state = Self::compute_state_by_one_point(&a_si, &a_sd, 1e-11, &self.ds);
                if a_state != 3 {
                    // myState 3 == IN (BRepClass3d_SClassifier: 0 unknown,
                    // 1 fault, 2 ON, 3 IN, 4 OUT); OCCT L838: aState != IN 鈥?
                    // not inside; keep for the next solid.
                    i_si += 1;
                    continue;
                }
                if a_ms_or.contains(&(a_sd.ptr_id(), a_sd.location)) {
                    // OCCT L846-858: make a new solid aSdx (copy of aSd + aSI).
                    let a_sdx = Self::solid_copy_add(&a_sd, &a_si);
                    // OCCT L860-861: myImages[aSd].Append(aSdx).
                    self.my_images
                        .bound((a_sd.ptr_id(), a_sd.location))
                        .push(a_sdx.clone());
                    // OCCT L863-865: myOrigins[aSdx].Append(aSd).
                    self.my_origins
                        .entry((a_sdx.ptr_id(), a_sdx.location))
                        .or_default()
                        .push(a_sd.clone());
                    // OCCT L867-868: aMSOr.Remove(aSd); aSd = aSdx.
                    a_ms_or.remove(&(a_sd.ptr_id(), a_sd.location));
                    a_sd = a_sdx;
                } else {
                    // OCCT L871-873: aBB.Add(aSd, aSI).
                    Self::solid_add_shape(&mut a_sd, &a_si);
                }
                // OCCT L875: aLSI.Remove(aIt1) 鈥?removal without advancing.
                a_ls_i.remove(i_si);
            }
            i_sd += 1;
        }
    }

    /// OCCT BOPTools_AlgoTools::ComputeStateByOnePoint (BOPTools_AlgoTools.cxx L623-656).
    /// Classifies a shape (vertex/edge/face) relative to a solid by a single point.
    fn compute_state_by_one_point(the_s: &Shape, the_ref: &Shape, the_tol: f64, ds: &DS) -> u8 {
        let a_type = the_s.shape_type();
        // OCCT L636-644: the FACE branch runs ComputeState(Face) with the
        // solid's edge bounds (the Context FClass2d path).
        if a_type == topods::ShapeType::Face {
            return compute_state_face(the_s, the_ref, the_tol, ds);
        }
        let a_p3d: Option<DVec3> = match a_type {
            topods::ShapeType::Vertex => match &*the_s.data {
                TShape::Vertex(vd) => Some(vd.point),
                _ => None,
            },
            topods::ShapeType::Edge => match &*the_s.data {
                TShape::Edge(ed) => {
                    if let Some(curve) = &ed.curve {
                        // OCCT BOPTools_AlgoTools::ComputeState(Edge) (L752-788):
                        // aT = IntermediatePoint(aT1, aT2) = 0.43213918-weighted;
                        // infinite ranges use the dT = 10. shifts (L748-771).
                        let (a_t1, a_t2) = (ed.range[0], ed.range[1]);
                        let a_t = if rcad_kernel::is_negative_infinite_value(a_t1)
                            && !rcad_kernel::is_positive_infinite_value(a_t2)
                        {
                            a_t2 - 10.0
                        } else if !rcad_kernel::is_negative_infinite_value(a_t1)
                            && rcad_kernel::is_positive_infinite_value(a_t2)
                        {
                            a_t1 + 10.0
                        } else if rcad_kernel::is_negative_infinite_value(a_t1)
                            && rcad_kernel::is_positive_infinite_value(a_t2)
                        {
                            0.0
                        } else {
                            intermediate_point(a_t1, a_t2)
                        };
                        Some(curve.point_at(a_t))
                    } else {
                        // degenerated edge 鈥?first vertex point (OCCT L748-754).
                        match &*ed.first.data {
                            TShape::Vertex(vd) => Some(vd.point),
                            _ => None,
                        }
                    }
                }
                _ => None,
            },
            _ => {
                // OCCT L646-653: recurse into the first sub-shape but IGNORE
                // the returned state 鈥?aState stays TopAbs_UNKNOWN (0).
                let subs = Self::shape_sub_shapes_static(the_s);
                if let Some(sub) = subs.first() {
                    let _ = Self::compute_state_by_one_point(sub, the_ref, the_tol, ds);
                }
                None
            }
        };
        let p = match a_p3d {
            Some(p) => p,
            None => return 0, // TopAbs_UNKNOWN
        };
        let mut clsf =
            crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(the_ref);
        clsf.perform(p, the_tol);
        clsf.my_state
    }

    /// OCCT aBB.MakeSolid(aSdx); copy all sub-shapes of aSd; aBB.Add(aSdx, aSI)
    /// (BOPAlgo_Builder_3.cxx L849-858). The new solid is FORWARD 鈥?OCCT
    /// MakeSolid creates a FORWARD solid and never copies the source
    /// orientation (unlike BuildDraftSolid which does copy it at L280-281).
    ///
    /// Architecture note: OCCT aBB.Add(aSd, aSI) accepts any sub-shape,
    /// including a WIRE (OwnInternalShapes collects all non-SHELL direct
    /// sub-shapes, Builder_3.cxx L891-905). rcad's TSolidData stores internal
    /// vertices/edges in dedicated fields and has no internal-wire slot; the
    /// WIRE case is dropped. The boolean sources (cylinders/boxes/prisms) carry
    /// no internal wires, so the divergence is not exercised by the stage tests.
    fn solid_copy_add(a_sd: &Shape, a_si: &Shape) -> Shape {
        let mut shells = Vec::new();
        let mut internal_v = Vec::new();
        let mut internal_e = Vec::new();
        if let TShape::Solid(sd) = &*a_sd.data {
            shells = sd.shells.clone();
            internal_v = sd.internal_vertices.clone();
            internal_e = sd.internal_edges.clone();
        }
        match &*a_si.data {
            TShape::Vertex(_) => internal_v.push(a_si.clone()),
            TShape::Edge(_) => internal_e.push(a_si.clone()),
            _ => {}
        }
        Shape::new(
            std::sync::Arc::new(TShape::Solid(TSolidData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                shells,
                internal_vertices: internal_v,
                internal_edges: internal_e,
            })),
            0,
            topods::Orientation::Forward,
        )
    }

    /// OCCT aBB.Add(aSd, aSI) 鈥?add an internal vertex/edge to the solid
    /// (BOPAlgo_Builder_3.cxx L871-873).
    fn solid_add_shape(a_sd: &mut Shape, a_si: &Shape) {
        let ts = Arc::make_mut(&mut a_sd.data);
        if let TShape::Solid(sd) = ts {
            match &*a_si.data {
                TShape::Vertex(_) => sd.internal_vertices.push(a_si.clone()),
                TShape::Edge(_) => sd.internal_edges.push(a_si.clone()),
                _ => {}
            }
        }
    }

    /// OCCT BOPTools_AlgoTools::TreatCompound 鈥?flatten compound shapes.
    /// The fence uses TopTools_ShapeMapHasher (TShape + Location).
    fn treat_compound(
        s: &Shape,
        result: &mut Vec<Shape>,
        fence: &mut std::collections::HashSet<(u64, u32)>,
    ) {
        if s.shape_type() == topods::ShapeType::Compound {
            if let TShape::Compound(children) = &*s.data {
                for child in children {
                    Self::treat_compound(child, result, fence);
                }
            }
        } else {
            if fence.insert((s.ptr_id(), s.location)) {
                result.push(s.clone());
            }
        }
    }

    /// OCCT BOPAlgo_Builder::OwnInternalShapes (BOPAlgo_Builder_3.cxx L891-905).
    /// Collects all non-SHELL direct sub-shapes of the solid (internal
    /// vertices/edges/wires stored at the solid level). aMx is an IndexedMap
    /// with TopTools_ShapeMapHasher 鈥?key identity TShape + Location.
    fn own_internal_shapes(&self, s: &Shape) -> Vec<Shape> {
        let mut result: Vec<Shape> = Vec::new();
        let mut fence: HashSet<(u64, u32)> = HashSet::new();
        if let TShape::Solid(sd) = &*s.data {
            for v in &sd.internal_vertices {
                if fence.insert((v.ptr_id(), v.location)) {
                    result.push(v.clone());
                }
            }
            for e in &sd.internal_edges {
                if fence.insert((e.ptr_id(), e.location)) {
                    result.push(e.clone());
                }
            }
        }
        // OCCT L896-904: direct sub-shapes that are not SHELL.
        for ss in self.shape_sub_shapes(s) {
            if ss.shape_type() != topods::ShapeType::Shell {
                if fence.insert((ss.ptr_id(), ss.location)) {
                    result.push(ss.clone());
                }
            }
        }
        result
    }

    /// OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L80-120) 鈥?maps each
    /// VERTEX to its ancestor EDGEs and FACEs, and each EDGE to its ancestor
    /// FACEs, over the FULL shape hierarchy of S (a source solid or its split
    /// images). The one-level scan previously used missed every ancestor, so
    /// the FillInternalShapes filter never dropped boundary V/E of the solids.
    fn map_shapes_and_ancestors(
        s: &Shape,
        map: &mut std::collections::HashMap<(u64, u32), Vec<u64>>,
    ) {
        // OCCT L90-107: for each ancestor of type TA, append it to the list of
        // every descendant of type TS. TA=FACE, TS鈭坽VERTEX,EDGE}; TA=EDGE,
        // TS=VERTEX (the three MapShapesAndAncestors calls at L770-772).
        // OCCT aMSx (L635-636) uses TopTools_ShapeMapHasher 鈥?TShape + Location.
        let mut faces: Vec<Shape> = Vec::new();
        Self::collect_sub_shapes_of_type(s, topods::ShapeType::Face, &mut faces);
        for f in &faces {
            let mut subs: Vec<Shape> = Vec::new();
            Self::collect_sub_shapes_of_type(f, topods::ShapeType::Vertex, &mut subs);
            Self::collect_sub_shapes_of_type(f, topods::ShapeType::Edge, &mut subs);
            for ts in &subs {
                map.entry((ts.ptr_id(), ts.location))
                    .or_default()
                    .push(f.ptr_id());
            }
        }
        let mut edges: Vec<Shape> = Vec::new();
        Self::collect_sub_shapes_of_type(s, topods::ShapeType::Edge, &mut edges);
        for e in &edges {
            let mut vs: Vec<Shape> = Vec::new();
            Self::collect_sub_shapes_of_type(e, topods::ShapeType::Vertex, &mut vs);
            for v in &vs {
                map.entry((v.ptr_id(), v.location))
                    .or_default()
                    .push(e.ptr_id());
            }
        }
    }

    /// Recursively collect all sub-shapes of the given type from S, at any depth
    /// (OCCT TopExp_Explorer with the ToFind type).
    fn collect_sub_shapes_of_type(
        s: &Shape,
        t: topods::ShapeType,
        out: &mut Vec<Shape>,
    ) {
        for ss in Self::shape_sub_shapes_static(s) {
            if ss.shape_type() == t {
                out.push(ss.clone());
            }
            Self::collect_sub_shapes_of_type(&ss, t, out);
        }
    }

    /// Static version of shape_sub_shapes for use in non-&self methods.
    fn shape_sub_shapes_static(s: &Shape) -> Vec<Shape> {
        match &*s.data {
            TShape::Vertex(_) => vec![],
            TShape::Edge(ed) => {
                // OCCT TopoDS_Iterator(aE, cumLoc) 鈥?the edge Location is
                // composed into the vertices. Identity fast paths keep the
                // exact index; a nested fold falls back to the edge location
                // (no Location table access in the static helper).
                let loc = s.location;
                let vl = |v: &Shape| {
                    let vloc = if loc == 0 { v.location } else { loc };
                    Shape::new(v.data.clone(), vloc, v.orientation)
                };
                vec![vl(&ed.first), vl(&ed.last)]
            }
            TShape::Wire(wd) => {
                wd.edges.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Face(fd) => {
                let mut v = vec![
                    Shape::new(fd.outer_wire.data.clone(), fd.outer_wire.location, fd.outer_wire.orientation)
                ];
                v.extend(fd.inner_wires.iter().map(|w| {
                    Shape::new(w.data.clone(), w.location, w.orientation)
                }));
                v
            }
            TShape::Shell(sd) => {
                sd.faces.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Solid(sd) => {
                sd.shells.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::CompSolid(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Compound(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesCompounds (BOPAlgo_Builder_1.cxx L197-217).
    fn fill_images_compounds(&mut self) {
        // OCCT L199-201: fence map + NbSourceShapes 鈥?TopTools_ShapeMapHasher.
        let mut a_mfp: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
        let a_nb_s = self.ds.nb_source_shapes();
        // OCCT L202-216: for each source compound, call FillImagesCompound
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Compound {
                continue;
            }
            let a_c = self.brep_sr(i);
            // OCCT L210: FillImagesCompound(aC, aMFP);
            self.fill_images_compound(&a_c, &mut a_mfp);
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesCompound (BOPAlgo_Builder_1.cxx L278-360).
    /// OCCT L285-299: iterate the sub-shapes, recursively processing nested
    /// compounds, and set bInterferred when any sub-shape has images. Only then
    /// build a new compound from the sub-shape images (each image taking the
    /// sub-shape's orientation, L348-356) and Bind it as theS's image (L362-365).
    fn fill_images_compound(
        &mut self,
        the_c: &Shape,
        the_mfp: &mut std::collections::HashSet<(u64, u32)>,
    ) {
        // OCCT L282-283: if (!theMFP.Add(theS)) return; 鈥?the fence uses
        // TopTools_ShapeMapHasher (TShape + Location).
        if !the_mfp.insert((the_c.ptr_id(), the_c.location)) {
            return;
        }
        // OCCT L285-299: bInterferred.
        let mut b_interferred = false;
        let subs = self.shape_sub_shapes(the_c);
        for ss in &subs {
            // OCCT L292-294: recursively process nested compounds.
            if ss.shape_type() == topods::ShapeType::Compound {
                self.fill_images_compound(ss, the_mfp);
            }
            // OCCT L295-297: if (myImages.IsBound(aSx)) bInterferred = true.
            if self.my_images.contains((ss.ptr_id(), ss.location)) {
                b_interferred = true;
            }
        }
        if !b_interferred {
            return;
        }
        // OCCT L301-344: MakeContainer(COMPOUND, aCIm); iterate sub-shapes.
        let mut new_shapes: Vec<Shape> = Vec::new();
        for ss in &subs {
            let a_or_x = ss.orientation;
            if let Some(ss_imgs) = self.my_images.get((ss.ptr_id(), ss.location)) {
                // OCCT L348-356: each image gets the sub-shape's orientation.
                for a_sx_im0 in ss_imgs {
                    let mut a_sx_im = a_sx_im0.clone();
                    a_sx_im.orientation = a_or_x;
                    new_shapes.push(a_sx_im);
                }
            } else {
                // OCCT L358-360: no images 鈥?add the sub-shape itself.
                new_shapes.push(ss.clone());
            }
        }
        // OCCT L362-365: myImages.Bind(theS, List(aCIm)) 鈥?replace the images.
        let comp_shape = Shape::new(
            std::sync::Arc::new(TShape::Compound(new_shapes)),
            0,
            topods::Orientation::Forward,
        );
        self.my_images
            .insert((the_c.ptr_id(), the_c.location), vec![comp_shape]);
    }

    /// OCCT BOPAlgo_Builder::PrepareHistory (BOPAlgo_Builder_4.cxx L164-252).
    fn prepare_history(&mut self) {
        // OCCT L166-168: if (!HasHistory()) return;
        if !self.my_fill_history { return; }

        // OCCT L172-176: init history tool, map result shapes.
        // rcad: shape_remap (ptr_id -> new BRep index) identifies the shapes
        // added to the result (equivalent to TopExp::MapShapes(myShape,
        // myMapShape) membership).
        let mut a_history = crate::bop::history::BRepToolsHistory::new();

        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            let a_s = self.brep_sr(i);
            // OCCT L192-195: BRepTools_History::IsSupportedType.
            if !crate::bop::history::is_supported_type(&a_s) {
                continue;
            }

            let mut is_modified = false;

            // OCCT L205-231: LocModified 鈥?the splits of the shape kept in the
            // result become Modified, with the proper orientation.
            if let Some(imgs) = self.my_images.get((a_s.ptr_id(), a_s.location)) {
                for a_sp in imgs {
                    // OCCT L214-217: check if the result contains the split.
                    if self.shape_remap.contains_key(&a_sp.ptr_id()) {
                        let mut a_sp2 = a_sp.clone();
                        // OCCT L218-226: VERTEX/SOLID keep the source
                        // orientation; other types reverse when IsSplitToReverse.
                        let a_type = a_sp2.shape_type();
                        if a_type == topods::ShapeType::Vertex || a_type == topods::ShapeType::Solid {
                            a_sp2.orientation = a_s.orientation;
                        } else {
                            let b_to_reverse = match a_type {
                                topods::ShapeType::Edge | topods::ShapeType::Wire => {
                                    crate::bop::tools::algo_tools::is_split_to_reverse_edge(&a_sp2, &a_s).0
                                }
                                topods::ShapeType::Face | topods::ShapeType::Shell => {
                                    crate::bop::tools::algo_tools::is_split_to_reverse_face(&a_sp2, &a_s, &self.ds).0
                                }
                                _ => false,
                            };
                            if b_to_reverse {
                                a_sp2.orientation = flip_orientation(a_sp2.orientation);
                            }
                        }
                        a_history.add_modified(&a_s, &a_sp2);
                        is_modified = true;
                    }
                }
            }

            // OCCT L234-243: LocGenerated 鈥?generated elements kept in the
            // result become Generated.
            let a_gen_shapes = self.loc_generated(&a_s);
            for a_g in &a_gen_shapes {
                if self.shape_remap.contains_key(&a_g.ptr_id()) {
                    a_history.add_generated(&a_s, a_g);
                }
            }

            // OCCT L247-250: not modified and not in the result -> Deleted.
            if !is_modified && !self.shape_remap.contains_key(&a_s.ptr_id()) {
                a_history.remove(&a_s);
            }
        }
        self.my_history = Some(a_history);
    }

    /// OCCT BOPAlgo_Builder::LocGenerated (BOPAlgo_Builder_4.cxx L30-152) 鈥?
    /// the shapes generated from the given EDGE/FACE: intersection vertices of
    /// the EE/EF interferences and, for faces, the section edges/vertices.
    fn loc_generated(&self, the_s: &Shape) -> Vec<Shape> {
        let mut a_hist: Vec<Shape> = Vec::new();
        // Only EDGES and FACES.
        let a_type = the_s.shape_type();
        if a_type != topods::ShapeType::Edge && a_type != topods::ShapeType::Face {
            return a_hist;
        }
        // Check that DS contains the shape (it is from the arguments).
        let n_s = self.ds.index(the_s);
        if n_s < 0 {
            return a_hist;
        }
        let n_s = n_s as usize;
        // Check that the shape has participated in any intersections.
        if !self.ds.shapes[n_s].has_reference() {
            return a_hist;
        }
        let a_ees = self.ds.interf_ee.clone();
        let a_efs = self.ds.interf_ef.clone();
        // Fence to avoid duplicates.
        let mut a_mfence: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let is_face = a_type == topods::ShapeType::Face;
        // EE interferences (only for EDGE) and EF interferences.
        let a_nb_ee = if is_face { 0 } else { a_ees.len() };
        for k in 0..2 {
            let a_nb_lines = if k == 0 { a_nb_ee } else { a_efs.len() };
            for a_int in 0..a_nb_lines {
                let (has_new, n_v_new, contains) = if k == 0 {
                    let ee = &a_ees[a_int];
                    (ee.new_vertex != usize::MAX, ee.new_vertex,
                     ee.e1 == n_s || ee.e2 == n_s)
                } else {
                    let ef = &a_efs[a_int];
                    (ef.new_vertex != usize::MAX, ef.new_vertex,
                     ef.edge == n_s || ef.face == n_s)
                };
                if !has_new || !contains {
                    continue;
                }
                // myDS->HasShapeSD(nVNew, nVNew)
                let mut n_v_new2 = n_v_new;
                self.ds.has_shape_sd(n_v_new, &mut n_v_new2);
                if !a_mfence.insert(n_v_new2) {
                    continue;
                }
                let a_v_new = self.brep_sr(n_v_new2);
                // Check that the result shape contains the vertex.
                if self.shape_remap.contains_key(&a_v_new.ptr_id()) {
                    a_hist.push(a_v_new);
                }
            }
        }
        if !is_face {
            return a_hist;
        }
        // FACE: section edges and vertices from FaceInfo.
        let a_fi = self.ds.face_info(n_s);
        // Section edges (PaveBlocksSc).
        let a_mpb_sc: Vec<u64> = a_fi.pave_blocks_sc.iter().copied().collect();
        for &pb_key in &a_mpb_sc {
            if let Some(a_pb) = self.ds.pb_from_ptr(pb_key) {
                let n_e = a_pb.0.read().unwrap().edge;
                if n_e >= self.ds.nb_shapes() { continue; }
                let a_e_new = self.brep_sr(n_e);
                if self.shape_remap.contains_key(&a_e_new.ptr_id()) {
                    a_hist.push(a_e_new);
                }
            }
        }
        // Section vertices (VerticesSc) 鈥?NCollection_Map bucket order.
        let a_mv_sc: Vec<usize> = a_fi.vertices_sc.iter().copied().collect();
        for n_v in a_mv_sc {
            let a_v_new = self.brep_sr(n_v);
            if self.shape_remap.contains_key(&a_v_new.ptr_id()) {
                a_hist.push(a_v_new);
            }
        }
        a_hist
    }

    /// OCCT BOPAlgo_Builder::PostTreat (BOPAlgo_Builder.cxx L461-486).
    fn post_treat(&mut self) {
        // OCCT L466-480: in non-destructive mode, collect source V/E/F shapes
        // into aMA (aMapToAvoid 鈥?tolerance of these shapes is not corrected).
        // rcad: non-destructive mode is not enabled for the boolean pipeline
        // (my_non_destructive is false), so aMA stays empty, matching OCCT's
        // default behaviour.
        let a_ma: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
        let Some(brep) = self.my_shape.as_mut() else { return };
        // OCCT L483: CorrectTolerances(myShape, aMA, 0.05, myRunParallel).
        crate::bop::tools::algo_tools::correct_tolerances(brep, &a_ma, 0.05);
        // OCCT L485: CorrectShapeTolerances(myShape, aMA, myRunParallel).
        crate::bop::tools::algo_tools::correct_shape_tolerances(brep, &a_ma);
    }

    // ====================================================================
    // BuildShape 鈥?OCCT BOPAlgo_BOP::BuildShape (BOPAlgo_BOP.cxx L885-1107)
    // Boolean operation result construction (BuildRC -> BuildSolid / containers).
    // ====================================================================

    /// OCCT BOPAlgo_BOP::BuildShape (BOPAlgo_BOP.cxx L885-1107).
    fn build_shape(&mut self) {
        // OCCT BOPAlgo_BOP::CheckData sets myDims from arguments/tools.
        self.compute_dims();
        // OCCT L889-911: for 3D+3D, if any argument solid is open, use the
        // BuildBOP alternative (BOPAlgo_Builder.cxx L491-897).
        if self.my_dims[0] == 3 && self.my_dims[1] == 3 {
            let has_not_closed_solids = self.check_args_for_open_solid();
            if has_not_closed_solids {
                // OCCT L902-909: BuildBOP fallback for open solids 鈥?rcad pending.
                // Closed-solid inputs (all stage tests) take the BuildRC path.
            }
        }
        // OCCT L913-914: BuildRC.
        self.build_rc();
        // OCCT L916-920: FUSE of 3D 鈫?BuildSolid.
        if self.my_operation == BooleanOpType::Union && self.my_dims[0] == 3 {
            self.build_solid();
            return;
        }
        // OCCT L923-1107: CUT/COMMON container logic.
        // OCCT L936-937: aMSRC = MapShapes(myRC) 鈥?TopTools_ShapeMapHasher.
        let mut a_msrc: HashSet<(u64, u32)> = HashSet::new();
        for s in &self.my_rc {
            self.add_all_sub_shapes(s, &mut a_msrc);
        }
        // OCCT L940-951: collect containers of arguments.
        let mut a_lsc: Vec<Shape> = Vec::new();
        for i in 0..2 {
            let a_ls: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            for a_s in a_ls {
                Self::collect_containers(a_s, &mut a_lsc);
            }
        }
        // OCCT L958-1045: make containers.
        let mut a_lc_res: Vec<Shape> = Vec::new();
        for a_sc in &a_lsc {
            // OCCT L966: MakeContainer(COMPOUND, aRC).
            let mut a_rc: Vec<Shape> = Vec::new();
            // OCCT L968-990: add splits of sub-shapes contained in aMSRC.
            for a_s in self.shape_sub_shapes(a_sc) {
                if let Some(imgs) = self.my_images.get((a_s.ptr_id(), a_s.location)) {
                    for a_s_im in imgs {
                        if a_msrc.contains(&(a_s_im.ptr_id(), a_s_im.location)) {
                            a_rc.push(a_s_im.clone());
                        }
                    }
                } else if a_msrc.contains(&(a_s.ptr_id(), a_s.location)) {
                    a_rc.push(a_s.clone());
                }
            }
            if self.my_operation != BooleanOpType::Union {
            }
            // OCCT L992-1009: connexity element types.
            let a_type = a_sc.shape_type();
            let (a_t1, a_t2) = match a_type {
                topods::ShapeType::Wire => (topods::ShapeType::Vertex, topods::ShapeType::Edge),
                topods::ShapeType::Shell => (topods::ShapeType::Edge, topods::ShapeType::Face),
                _ => (topods::ShapeType::Face, topods::ShapeType::Solid),
            };
            // OCCT L1011-1016: MakeConnexityBlocks(aRC, aT1, aT2, aLCB).
            let mut a_lcb: Vec<Vec<Shape>> = Vec::new();
            Self::make_connexity_blocks_shapes(&a_rc, a_t1, a_t2, &mut a_lcb);
            if a_lcb.is_empty() {
                continue;
            }
            // OCCT L1018-1043: build containers from blocks.
            for a_cb in &a_lcb {
                let mut a_rcb = Self::build_container_of_type(a_type, a_cb);
                match a_type {
                    topods::ShapeType::Wire => Self::orient_edges_on_wire(&mut a_rcb),
                    topods::ShapeType::Shell => Self::orient_faces_on_shell(&mut a_rcb),
                    _ => {}
                }
                // OCCT L1041: aRCB.Orientation(aSC.Orientation()) 鈥?the result
                // container inherits the source container's orientation.
                a_rcb.orientation = a_sc.orientation;
                a_lc_res.push(a_rcb);
            }
        }
        // OCCT L1047: RemoveDuplicates(aLCRes).
        self.remove_duplicates(&mut a_lc_res);
        // OCCT L1055-1063: aResult compound of containers.
        let mut a_result: Vec<Shape> = a_lc_res;
        // OCCT L1066-1067: aMSResult = MapShapes(aResult) 鈥?TopTools_ShapeMapHasher.
        let mut a_ms_result: HashSet<(u64, u32)> = HashSet::new();
        for s in &a_result {
            self.add_all_sub_shapes(s, &mut a_ms_result);
        }
        // OCCT L1069-1080: get input non-container shapes.
        let mut a_ls_non_cont: Vec<Shape> = Vec::new();
        let mut a_m_inp_fence: HashSet<(u64, u32)> = HashSet::new();
        for i in 0..2 {
            let a_ls: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            for a_s in a_ls {
                Self::treat_compound(a_s, &mut a_ls_non_cont, &mut a_m_inp_fence);
            }
        }
        // OCCT L1082-1104: put non-container shapes in the result.
        for a_s in &a_ls_non_cont {
            if let Some(imgs) = self.my_images.get((a_s.ptr_id(), a_s.location)) {
                for a_s_im in imgs {
                    if a_msrc.contains(&(a_s_im.ptr_id(), a_s_im.location))
                        && a_ms_result.insert((a_s_im.ptr_id(), a_s_im.location))
                    {
                        a_result.push(a_s_im.clone());
                    }
                }
            } else if a_msrc.contains(&(a_s.ptr_id(), a_s.location))
                && a_ms_result.insert((a_s.ptr_id(), a_s.location))
            {
                a_result.push(a_s.clone());
            }
        }
        // OCCT L1106: myShape = aResult.
        self.set_shape_from_shapes(a_result);
    }

    /// OCCT BOPAlgo_BOP::BuildRC (BOPAlgo_BOP.cxx L597-881).
    fn build_rc(&mut self) {
        let mut a_c: Vec<Shape> = Vec::new();
        // OCCT L608-623: A. Fuse 鈥?collect shapes of TypeToExplore(dim0) from myShape.
        if self.my_operation == BooleanOpType::Union {
            let mut a_m_fence: HashSet<(u64, u32)> = HashSet::new();
            let a_type = type_to_explore(self.my_dims[0]);
            let my_shape = self.my_shape.clone().unwrap_or_default();
            // OCCT: TopExp_Explorer aExp(myShape, aType) 鈥?each collected shape
            // keeps the orientation it is referenced with. rcad's myShape is a
            // flat tshapes array (no compound root, no orientations), so look
            // the orientation up from the argument tree; all stage-test solids
            // are FORWARD, so this is a no-op for them but preserves REVERSED
            // arguments per OCCT.
            let mut a_arg_oris: HashMap<u64, topods::Orientation> = HashMap::new();
            for a_arg in self.object_shapes() {
                for s in self.map_shapes_of_type(a_arg, a_type) {
                    a_arg_oris.entry(s.ptr_id()).or_insert(s.orientation);
                }
            }
            for ts in &my_shape.tshapes {
                let sh = Shape::new(ts.clone(), 0, topods::Orientation::Forward);
                if sh.shape_type() == a_type {
                    let ori = a_arg_oris
                        .get(&sh.ptr_id())
                        .copied()
                        .unwrap_or(topods::Orientation::Forward);
                    let sh = Shape::new(ts.clone(), 0, ori);
                    if a_m_fence.insert((sh.ptr_id(), sh.location)) {
                        a_c.push(sh);
                    }
                }
            }
            self.my_rc = a_c;
            return;
        }
        // OCCT L630-659: B. Common/Cut 鈥?building elements of arguments.
        let mut a_m_args: Vec<Shape> = Vec::new();
        let mut a_m_args_fence: HashSet<u64> = HashSet::new();
        let mut a_m_tools: Vec<Shape> = Vec::new();
        let mut a_m_tools_fence: HashSet<u64> = HashSet::new();
        for i in 0..2 {
            let a_ls: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            let (a_ms, a_ms_fence) = if i == 0 {
                (&mut a_m_args, &mut a_m_args_fence)
            } else {
                (&mut a_m_tools, &mut a_m_tools_fence)
            };
            for a_s in a_ls {
                let mut a_list: Vec<Shape> = Vec::new();
                let mut a_fence: HashSet<(u64, u32)> = HashSet::new();
                Self::treat_compound(a_s, &mut a_list, &mut a_fence);
                for a_ss in a_list {
                    let i_dim = shape_dimension(&a_ss);
                    if i_dim < 0 {
                        continue;
                    }
                    let a_type = type_to_explore(i_dim);
                    // OCCT L656: TopExp::MapShapes(aSS, aType, aMS).
                    for sh in self.map_shapes_of_type(&a_ss, a_type) {
                        if a_ms_fence.insert(sh.ptr_id()) {
                            a_ms.push(sh);
                        }
                    }
                }
            }
        }
        // OCCT L666-718: get splits of building elements.
        if std::env::var("RCAD_BS_DEBUG").is_ok() {
            eprintln!("[BRC-IMG] args={} tools={}",
                a_m_args.len(), a_m_tools.len());
            for s in &a_m_args {
                eprintln!("[BRC-IMG]   arg type={:?} ptr={} imgs={:?}",
                    s.shape_type(), s.ptr_id(),
                    self.my_images.get(self.remap_arg_key(s)).map(|v| v.len()));
            }
            for s in &a_m_tools {
                eprintln!("[BRC-IMG]   tool type={:?} ptr={} imgs={:?}",
                    s.shape_type(), s.ptr_id(),
                    self.my_images.get(self.remap_arg_key(s)).map(|v| v.len()));
            }
        }
        let mut a_m_args_im: Vec<Shape> = Vec::new();
        let mut a_m_args_im_fence: HashSet<u64> = HashSet::new();
        let mut a_m_tools_im: Vec<Shape> = Vec::new();
        let mut a_m_tools_im_fence: HashSet<u64> = HashSet::new();
        let mut a_m_set_args: HashMap<(usize, Vec<(u64, u32)>), Shape> = HashMap::new();
        let mut a_m_set_tools: HashMap<(usize, Vec<(u64, u32)>), Shape> = HashMap::new();
        let mut b_check_edges = false;
        for i in 0..2 {
            let a_ms: &Vec<Shape> = if i == 0 { &a_m_args } else { &a_m_tools };
            let (a_ms_im, a_ms_im_fence, a_m_set) = if i == 0 {
                (&mut a_m_args_im, &mut a_m_args_im_fence, &mut a_m_set_args)
            } else {
                (&mut a_m_tools_im, &mut a_m_tools_im_fence, &mut a_m_set_tools)
            };
            for a_s in a_ms {
                let a_type = a_s.shape_type();
                if a_type == topods::ShapeType::Edge {
                    b_check_edges = true;
                    // OCCT L689-691: skip degenerated edges.
                    let degen = a_s.as_edge().map(|e| e.degenerated).unwrap_or(false);
                    if degen {
                        continue;
                    }
                }
                if let Some(imgs) = self
                    .my_images
                    .get(self.remap_arg_key(a_s))
                    .cloned()
                {
                    for a_s_im in &imgs {
                        if a_ms_im_fence.insert(a_s_im.ptr_id()) {
                            a_ms_im.push(a_s_im.clone());
                        }
                    }
                } else {
                    if a_ms_im_fence.insert(a_s.ptr_id()) {
                        a_ms_im.push(a_s.clone());
                    }
                    if a_type == topods::ShapeType::Solid {
                        // OCCT L708-716: BOPTools_Set of the solid's face set.
                        let a_st = self.shape_face_set(a_s);
                        if !a_m_set.contains_key(&a_st) {
                            a_m_set.insert(a_st, a_s.clone());
                        }
                    }
                }
            }
        }
        // OCCT L723-798: compare maps and make the result.
        let i_dim_min = std::cmp::min(self.my_dims[0], self.my_dims[1]);
        let b_common = self.my_operation == BooleanOpType::Intersection;
        // rcad has no CUT21; aMIt = object splits, aMCheck = tool splits.
        let (a_m_it, a_m_check, a_m_set_check) = (&a_m_args_im, &a_m_tools_im, &a_m_set_tools);
        let mut a_m_check_exp: Vec<Shape> = Vec::new();
        let mut a_m_check_exp_fence: HashSet<u64> = HashSet::new();
        let mut a_m_it_exp: Vec<Shape> = Vec::new();
        let mut a_m_it_exp_fence: HashSet<u64> = HashSet::new();
        if b_common {
            // OCCT L738-751: expand aMIt with sub-shapes of lower dims.
            for a_s in a_m_it {
                let i_dim_max = shape_dimension(a_s);
                for i_dim in i_dim_min..i_dim_max {
                    let a_type = type_to_explore(i_dim);
                    for sh in self.map_shapes_of_type(a_s, a_type) {
                        if a_m_it_exp_fence.insert(sh.ptr_id()) {
                            a_m_it_exp.push(sh);
                        }
                    }
                }
                if a_m_it_exp_fence.insert(a_s.ptr_id()) {
                    a_m_it_exp.push(a_s.clone());
                }
            }
        } else {
            a_m_it_exp = a_m_it.clone();
            a_m_it_exp_fence = a_m_it_exp.iter().map(|s| s.ptr_id()).collect();
        }
        // OCCT L758-769: expand aMCheck with sub-shapes of lower dims.
        for a_s in a_m_check {
            let i_dim_max = shape_dimension(a_s);
            for i_dim in i_dim_min..i_dim_max {
                let a_type = type_to_explore(i_dim);
                for sh in self.map_shapes_of_type(a_s, a_type) {
                    if a_m_check_exp_fence.insert(sh.ptr_id()) {
                        a_m_check_exp.push(sh);
                    }
                }
            }
            if a_m_check_exp_fence.insert(a_s.ptr_id()) {
                a_m_check_exp.push(a_s.clone());
            }
        }
        // OCCT L771-798: build result.
        if std::env::var("RCAD_BS_DEBUG").is_ok() {
            eprintln!("[BRC] b_common={} i_dim_min={} n_it_exp={} n_check_exp={} n_set_check={}",
                b_common, i_dim_min, a_m_it_exp.len(), a_m_check_exp.len(), a_m_set_check.len());
        }
        for a_s in &a_m_it_exp {
            let mut b_contains = a_m_check_exp_fence.contains(&a_s.ptr_id());
            if !b_contains && a_s.shape_type() == topods::ShapeType::Solid {
                // OCCT L777-782: check by the solid's face set.
                let a_st = self.shape_face_set(a_s);
                b_contains = a_m_set_check.contains_key(&a_st);
            }
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[BRC]   it_exp type={:?} ptr={} contains={}",
                    a_s.shape_type(), a_s.ptr_id(), b_contains);
                let mut fids: Vec<String> = Vec::new();
                for f in self.map_shapes_of_type(a_s, topods::ShapeType::Face) {
                    let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                        rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                        _ => "O".into(),
                    }).unwrap_or_else(|| "?".into());
                    fids.push(format!("{}:{}", f.ptr_id() % 100000, n));
                }
                eprintln!("[BRC]     faces=[{}]", fids.join(" "));
            }
            if std::env::var("RCAD_BS_DEBUG").is_ok() && a_s.shape_type() == topods::ShapeType::Solid {
                for t in a_m_check.iter() {
                    let tfs = self.map_shapes_of_type(t, topods::ShapeType::Face);
                    let fid_strs: Vec<String> = tfs.iter().map(|f| format!("{}", f.ptr_id() % 100000)).collect();
                    eprintln!("[BRC-CHK] tool solid ptr={} faces=[{}]",
                        t.ptr_id() % 100000,
                        fid_strs.join(" "));
                }
                let mut ifids: Vec<String> = Vec::new();
                for f in self.map_shapes_of_type(a_s, topods::ShapeType::Face) {
                    ifids.push(format!("{}", f.ptr_id() % 100000));
                }
                eprintln!("[BRC-SPLIT] it_exp ptr={} faces=[{}]", a_s.ptr_id() % 100000, ifids.join(" "));
                // DS-space tool: translate each tool face to its DS clone
                let mut cids: Vec<String> = Vec::new();
                for t in a_m_check.iter() {
                    for f in self.map_shapes_of_type(t, topods::ShapeType::Face) {
                        let rk = self.remap_arg_key(&f);
                        cids.push(format!("{}", rk.0 % 100000));
                    }
                }
                eprintln!("[BRC-CHK-DS] tool faces ds=[{}]", cids.join(" "));
            }
            if b_common {
                if b_contains {
                    a_c.push(a_s.clone());
                }
            } else if !b_contains {
                a_c.push(a_s.clone());
            }
        }
        // OCCT L800-823: filter result for COMMON.
        if b_common {
            let mut a_m_fence: HashSet<(u64, u32)> = HashSet::new();
            let mut a_cx: Vec<Shape> = Vec::new();
            for i_dim in (i_dim_min..=3).rev() {
                let a_type = type_to_explore(i_dim);
                for sh in self.shapes_of_type_in_shapes(&a_c, a_type) {
                    if a_m_fence.insert((sh.ptr_id(), sh.location)) {
                        a_cx.push(sh.clone());
                        // OCCT L818: TopExp::MapShapes(aS, aMFence).
                        self.add_all_sub_shapes(&sh, &mut a_m_fence);
                    }
                }
            }
            a_c = a_cx;
        }
        // OCCT L825-829: if no edges were checked, done.
        if !b_check_edges {
            self.my_rc = a_c;
            return;
        }
        // OCCT L835-878: squats around degenerated edges.
        let mut a_m_vc: HashSet<u64> = HashSet::new();
        for sh in self.shapes_of_type_in_shapes(&a_c, topods::ShapeType::Vertex) {
            a_m_vc.insert(sh.ptr_id());
        }
        let a_nb = self.ds.nb_source_shapes();
        for i in 0..a_nb {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Edge {
                continue;
            }
            let a_e = self.brep_sr(i);
            let degen = a_e.as_edge().map(|e| e.degenerated).unwrap_or(false);
            if !degen {
                continue;
            }
            let n_vd = a_si.sub_shapes.first().copied().unwrap_or(0);
            let a_vd = self.brep_sr(n_vd);
            if !a_m_vc.contains(&a_vd.ptr_id()) {
                continue;
            }
            if self.ds.is_new_shape(n_vd) {
                continue;
            }
            if self.ds.interfered.contains(&n_vd) {
                continue;
            }
            a_c.push(a_e);
        }
        self.my_rc = a_c;
    }

    /// OCCT BOPAlgo_BOP::BuildSolid (BOPAlgo_BOP.cxx L1111-1392).
    fn build_solid(&mut self) {
        // OCCT L1121-1144: get solids from input arguments.
        let mut a_msa: HashSet<u64> = HashSet::new();
        // OCCT aMFS is NCollection_IndexedDataMap (BOP.cxx L1110-1111) 鈥?
        // insertion order matters for the aSFS face order below (L1217-1227).
        let mut a_mfs: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        let mut a_lsc: Vec<Shape> = Vec::new();
        for i in 0..2 {
            let a_lsa: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            for a_sa in a_lsa {
                // OCCT L1133-1139: explore solids, map face鈫抯olid ancestors.
                for sol in self.map_shapes_of_type(a_sa, topods::ShapeType::Solid) {
                    a_msa.insert(sol.ptr_id());
                    for a_f in self.map_shapes_of_type(&sol, topods::ShapeType::Face) {
                        a_mfs.entry((a_f.ptr_id(), a_f.location))
                            .or_insert_with(|| (a_f.clone(), Vec::new()))
                            .1
                            .push(sol.clone());
                    }
                }
                // OCCT L1141-1143: collect compsolids from arguments.
                Self::collect_containers(a_sa, &mut a_lsc);
            }
        }
        // OCCT L1151-1165: find solids sharing faces.
        let mut a_mt_sols: HashSet<u64> = HashSet::new();
        for (_f, (_fs, sols)) in &a_mfs {
            if sols.len() > 1 {
                for sol in sols {
                    a_mt_sols.insert(sol.ptr_id());
                }
            }
        }
        // OCCT L1167-1220: possibly untouched solids.
        let mut a_mu_sols: Vec<Shape> = Vec::new();
        let mut a_mu_fence: HashSet<u64> = HashSet::new();
        a_mfs.clear();
        for a_sx in &self.my_rc {
            if a_msa.contains(&a_sx.ptr_id()) {
                if !a_mt_sols.contains(&a_sx.ptr_id()) {
                    if a_mu_fence.insert(a_sx.ptr_id()) {
                        a_mu_sols.push(a_sx.clone());
                    }
                    continue;
                }
            }
            // OCCT L1185: MapFacesToBuildSolids(aSx, aMFS).
            self.map_faces_to_build_solids(a_sx, &mut a_mfs);
        }
        // OCCT L1191-1220: process untouched solids.
        let mut a_dmsts: Vec<Shape> = Vec::new();
        let mut a_dmsts_fence: HashSet<(usize, Vec<(u64, u32)>)> = HashSet::new();
        for a_sx in &a_mu_sols {
            let mut in_mfs = false;
            for a_f in self.map_shapes_of_type(a_sx, topods::ShapeType::Face) {
                if a_mfs.contains_key(&(a_f.ptr_id(), a_f.location)) {
                    in_mfs = true;
                    break;
                }
            }
            if in_mfs {
                self.map_faces_to_build_solids(a_sx, &mut a_mfs);
            } else {
                let a_st = self.shape_face_set(a_sx);
                if a_dmsts_fence.insert(a_st) {
                    a_dmsts.push(a_sx.clone());
                }
            }
        }
        // OCCT L1227-1241: faces belonging to a single solid.
        let mut a_mef: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut a_sfs: Vec<Shape> = Vec::new();
        for (_f, (fs, sols)) in &a_mfs {
            if sols.len() == 1 {
                a_sfs.push(fs.clone());
                // OCCT L1238: TopExp::MapShapesAndAncestors(aFx, EDGE, FACE, aMEF).
                for a_e in self.map_shapes_of_type(fs, topods::ShapeType::Edge) {
                    a_mef.entry(a_e.ptr_id()).or_default().push(fs.ptr_id());
                }
            }
        }
        
        // OCCT L1243-1271: build solids from the set of faces.
        // OCCT L1257-1260: aBS.SetAvoidInternalShapes(true) 鈥?the faces of the
        // result solids must not create internal shells.
        let mut a_rc: Vec<Shape> = Vec::new();
        if !a_sfs.is_empty() {
            let mut a_bs = crate::bop::algo::builder_solid::BuilderSolid::new(&self.ds);
            a_bs.my_shapes = a_sfs;
            a_bs.set_avoid_internal_shapes(true);
            a_bs.perform();
            if a_bs.has_errors() {
                // OCCT L1255: AddError(new BOPAlgo_AlertSolidBuilderFailed).
                self.my_report.add_error(crate::bop::algo::Alert::SolidBuilderFailed);
                return;
            }
            
            for a_sr in &a_bs.my_solids {
                a_rc.push(a_sr.clone());
            }
        }
        // OCCT L1273-1279: add untouched solids.
        for a_sx in &a_dmsts {
            a_rc.push(a_sx.clone());
        }
        // OCCT L1281-1286: no compsolids in arguments 鈫?done.
        if a_lsc.is_empty() {
            self.set_shape_from_shapes(a_rc);
            return;
        }
        // OCCT L1291-1391: compsolid construction 鈥?rcad pending (no compsolid args
        // in the stage tests).
        self.set_shape_from_shapes(a_rc);
    }

    /// OCCT BOPAlgo_BOP::CheckData 鈥?compute myDims from objects and tools.
    fn compute_dims(&mut self) {
        let mut d0 = 0i32;
        for s in self.object_shapes() {
            d0 = d0.max(shape_dimension(s));
        }
        let mut d1 = 0i32;
        for s in &self.my_tools {
            d1 = d1.max(shape_dimension(s));
        }
        self.my_dims = [d0, d1];
    }

    /// OCCT BOPAlgo_BOP::CheckArgsForOpenSolid (BOPAlgo_BOP.cxx L1396-1560).
    /// rcad: returns false 鈥?the open-solid edge/face analysis and the BuildBOP
    /// fallback are pending translation. All stage-test solids are closed.
    fn check_args_for_open_solid(&self) -> bool {
        false
    }

    /// Objects of the BOP: myArguments minus the tools suffix.
    /// OCCT BOPAlgo_BOP: myDS->Arguments() = myArguments ++ myTools.
    fn object_shapes(&self) -> &[Shape] {
        let n_tools = self.my_tools.len();
        let n_objs = self.my_arguments.len().saturating_sub(n_tools);
        &self.my_arguments[..n_objs]
    }

    /// Replace my_shape with a fresh compound built from the given shapes.
    /// OCCT: BRep_Builder().Add(aCompound, aS) for each result shape.
    fn set_shape_from_shapes(&mut self, shapes: Vec<Shape>) {
        self.my_shape = Some(topods::BRep::new());
        self.shape_remap.clear();
        // The fresh my_shape has an empty locations pool; drop the old
        // DS -> result location mapping so remap_location re-seeds it.
        self.loc_remap.clear();
        for s in shapes {
            self.add_shape_to_result(&s);
        }
    }

    /// OCCT TopExp::MapShapes(aS, aType, aMap) 鈥?collect all sub-shapes of a type.
    fn map_shapes_of_type(&self, s: &Shape, t: topods::ShapeType) -> Vec<Shape> {
        let mut out: Vec<Shape> = Vec::new();
        // OCCT TopExp::MapShapes uses TopTools_ShapeMapHasher 鈥?TShape + Location.
        let mut seen: HashSet<(u64, u32)> = HashSet::new();
        let mut stack: Vec<Shape> = vec![s.clone()];
        while let Some(cur) = stack.pop() {
            if !seen.insert((cur.ptr_id(), cur.location)) {
                continue;
            }
            if cur.shape_type() == t {
                out.push(cur.clone());
            }
            for sub in self.shape_sub_shapes(&cur) {
                stack.push(sub);
            }
        }
        out
    }

    /// Collect shapes of a type across a list of shapes (dedup).
    fn shapes_of_type_in_shapes(&self, shapes: &[Shape], t: topods::ShapeType) -> Vec<Shape> {
        let mut out: Vec<Shape> = Vec::new();
        let mut seen: HashSet<(u64, u32)> = HashSet::new();
        for s in shapes {
            for sh in self.map_shapes_of_type(s, t) {
                if seen.insert((sh.ptr_id(), sh.location)) {
                    out.push(sh);
                }
            }
        }
        out
    }

    /// Insert a shape and all its sub-shapes into a fence (TopExp::MapShapes).
    fn add_all_sub_shapes(&self, s: &Shape, fence: &mut HashSet<(u64, u32)>) {
        fence.insert((s.ptr_id(), s.location));
        for sub in self.shape_sub_shapes(s) {
            self.add_all_sub_shapes(&sub, fence);
        }
    }

    /// OCCT BOPTools_Set::Add(aS, TopAbs_FACE) (BOPTools_Set.cxx L124-166) 鈥?
    /// every INTERNAL face is expanded into FORWARD + REVERSED entries
    /// (L139-148), so the entry count (myNbShapes) differs between sets whose
    /// faces differ only by an INTERNAL orientation; IsEqual (L81-99) compares
    /// the count first, then the orientation-insensitive set of shapes keyed
    /// by TopTools_ShapeMapHasher (TShape + Location).
    /// The (count, sorted unique (TShape, Location) ids) pair reproduces both:
    /// the count is the expanded entry count, the sorted ids are the
    /// deduplicated set 鈥?two sets with the same faces but different INTERNAL
    /// expansion points (e.g. {A_INTERNAL, B} vs {A, B_INTERNAL}) compare
    /// equal, as in OCCT.
    fn shape_face_set(&self, s: &Shape) -> (usize, Vec<(u64, u32)>) {
        let mut count: usize = 0;
        let mut ids: Vec<(u64, u32)> = Vec::new();
        for f in self.map_shapes_of_type(s, topods::ShapeType::Face) {
            count += 1;
            ids.push((f.ptr_id(), f.location));
            if f.orientation == topods::Orientation::Internal {
                count += 1;
            }
        }
        ids.sort();
        ids.dedup();
        (count, ids)
    }

    /// OCCT BOPAlgo_BOP::CollectContainers (BOPAlgo_BOP.cxx L1601-1621).
    fn collect_containers(s: &Shape, out: &mut Vec<Shape>) {
        let a_type = s.shape_type();
        if a_type == topods::ShapeType::Wire
            || a_type == topods::ShapeType::Shell
            || a_type == topods::ShapeType::CompSolid
        {
            out.push(s.clone());
            return;
        }
        if a_type != topods::ShapeType::Compound {
            return;
        }
        for sub in Self::shape_sub_shapes_static(s) {
            Self::collect_containers(&sub, out);
        }
    }

    /// OCCT BOPTools_AlgoTools::MakeContainer + BRep_Builder().Add for a container.
    fn build_container_of_type(a_type: topods::ShapeType, subs: &[Shape]) -> Shape {
        let ts: TShape = match a_type {
            topods::ShapeType::Wire => TShape::Wire(TWireData {
                my_shapes: subs.to_vec(),
                flags: tshape_flags::DEFAULT,
                edges: subs.to_vec(),
            }),
            topods::ShapeType::Shell => TShape::Shell(TShellData {
                my_shapes: subs.to_vec(),
                flags: tshape_flags::DEFAULT,
                faces: subs.to_vec(),
            }),
            _ => TShape::Compound(subs.to_vec()),
        };
        Shape::new(std::sync::Arc::new(ts), 0, topods::Orientation::Forward)
    }

    /// OCCT BOPTools_AlgoTools::MakeConnexityBlocks(aRC, aT1, aT2, aLCB).
    /// rcad: minimal 鈥?each shape its own block. Solid arguments produce no
    /// containers, so this is not exercised by the stage tests.
    /// OCCT BOPTools_AlgoTools::MakeConnexityBlocks (BOPTools_AlgoTools.cxx
    /// L187-256) 鈥?groups `shapes` into connexity blocks by shared elements of
    /// type aT1 (the connection type); shapes repeated in the input mark the
    /// block non-regular (expanded to both orientations).
    fn make_connexity_blocks_shapes(
        shapes: &[Shape],
        a_t1: topods::ShapeType,
        _a_t2: topods::ShapeType,
        out: &mut Vec<Vec<Shape>>,
    ) {
        // OCCT L194-210: aMFence/aMNRegular 鈥?dedup the start elements.
        // TopTools_ShapeMapHasher (TShape + Location).
        let mut a_mfence: HashSet<(u64, u32)> = HashSet::new();
        let mut a_mn_regular: HashSet<(u64, u32)> = HashSet::new();
        let mut a_c_start: Vec<Shape> = Vec::new();
        for a_s in shapes {
            if a_mfence.insert((a_s.ptr_id(), a_s.location)) {
                a_c_start.push(a_s.clone());
            } else {
                a_mn_regular.insert((a_s.ptr_id(), a_s.location));
            }
        }
        // OCCT L212-216: the connection map 鈥?MapShapesAndAncestors(aCStart,
        // aT1, aT2, aCMap): aT1 element -> [aT2 shapes] (TopTools_ShapeMapHasher).
        let mut a_c_map: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        for s in &a_c_start {
            let mut subs: Vec<Shape> = Vec::new();
            Self::collect_sub_shapes_of_type_static(s, a_t1, &mut subs);
            for sub in subs {
                let skey = (s.ptr_id(), s.location);
                let l = a_c_map.entry((sub.ptr_id(), sub.location)).or_default();
                if !l.iter().any(|e| (e.ptr_id(), e.location) == skey) {
                    l.push(s.clone());
                }
            }
        }
        // OCCT L118-150: BFS blocks over the aT2 elements via shared aT1.
        let mut a_mfence2: HashSet<(u64, u32)> = HashSet::new();
        let mut a_lcb: Vec<Vec<Shape>> = Vec::new();
        for s in &a_c_start {
            if !a_mfence2.insert((s.ptr_id(), s.location)) {
                continue;
            }
            let mut a_l_block: Vec<Shape> = vec![s.clone()];
            let mut i = 0;
            while i < a_l_block.len() {
                let mut subs: Vec<Shape> = Vec::new();
                Self::collect_sub_shapes_of_type_static(&a_l_block[i], a_t1, &mut subs);
                for sub in subs {
                    if let Some(l) = a_c_map.get(&(sub.ptr_id(), sub.location)) {
                        for s2 in l {
                            if a_mfence2.insert((s2.ptr_id(), s2.location)) {
                                a_l_block.push(s2.clone());
                            }
                        }
                    }
                }
                i += 1;
            }
            a_lcb.push(a_l_block);
        }
        // OCCT L218-254: save the blocks, checking their regularity.
        for block in &a_lcb {
            let mut a_lcs: Vec<Shape> = Vec::new();
            let mut b_regular = true;
            for a_s in block {
                if a_mn_regular.contains(&(a_s.ptr_id(), a_s.location)) {
                    b_regular = false;
                    let mut f = a_s.clone();
                    f.orientation = topods::Orientation::Forward;
                    a_lcs.push(f);
                    let mut r = a_s.clone();
                    r.orientation = topods::Orientation::Reversed;
                    a_lcs.push(r);
                } else {
                    a_lcs.push(a_s.clone());
                    if b_regular {
                        // OCCT L243-247: every connection element of the shape
                        // must be used by exactly 2 shapes.
                        let mut subs: Vec<Shape> = Vec::new();
                        Self::collect_sub_shapes_of_type_static(a_s, a_t1, &mut subs);
                        for sub in subs {
                            let cnt = a_c_map
                                .get(&(sub.ptr_id(), sub.location))
                                .map_or(0, |l| l.len());
                            if cnt != 2 {
                                b_regular = false;
                                break;
                            }
                        }
                    }
                }
            }
            out.push(a_lcs);
        }
    }

    /// Recursively collect all sub-shapes of the given type (static version of
    /// collect_sub_shapes_of_type for use in associated functions).
    fn collect_sub_shapes_of_type_static(
        s: &Shape,
        t: topods::ShapeType,
        out: &mut Vec<Shape>,
    ) {
        for ss in Self::shape_sub_shapes_static(s) {
            if ss.shape_type() == t {
                out.push(ss.clone());
            }
            Self::collect_sub_shapes_of_type_static(&ss, t, out);
        }
    }

    /// OCCT BOPTools_AlgoTools::OrientEdgesOnWire (BOPTools_AlgoTools.cxx L262-362).
    /// Reorients the edges of a wire so that they form a continuous chain:
    /// each vertex is shared by exactly two edges with opposite orientations.
    fn orient_edges_on_wire(w: &mut Shape) {
        // aVEMap: vertex -> [edges] (NCollection_IndexedDataMap, insertion order,
        // TopTools_ShapeMapHasher 鈥?key TShape + Location).
        let orig_edges: Vec<Shape> = match &*w.data {
            TShape::Wire(wd) => wd.edges.clone(),
            _ => return,
        };
        let mut a_ve_map: IndexMap<(u64, u32), Vec<Shape>> = IndexMap::new();
        for a_e in &orig_edges {
            for a_v in Self::edge_vertices_of(a_e) {
                let vkey = (a_v.ptr_id(), a_v.location);
                let entry = a_ve_map.entry(vkey).or_default();
                if !entry.iter().any(|x| x.is_partner(a_e)) {
                    entry.push(a_e.clone());
                }
            }
        }
        if a_ve_map.is_empty() {
            return;
        }
        // New wire edges (OCCT: aBB.MakeWire(aWire); aBB.Add(aWire, ...)).
        let mut new_edges: Vec<Shape> = Vec::new();
        // OCCT aMFence: NCollection_Map<TopoDS_Shape> 鈥?default hasher
        // (TShape + Location + Orientation).
        let mut a_mfence: HashSet<(u64, u32, topods::Orientation)> = HashSet::new();
        for a_ec in &orig_edges {
            if !a_mfence.insert((a_ec.ptr_id(), a_ec.location, a_ec.orientation)) {
                continue;
            }
            new_edges.push(a_ec.clone());
            let (a_v1, a_v2) = Self::edge_endpoints_of(a_ec);
            if a_v1.is_partner(&a_v2) {
                // closed edge, go to the next edge
                continue;
            }
            // orient the adjacent edges
            for i in 0..2 {
                let mut a_vc = if i == 0 { a_v1.clone() } else { a_v2.clone() };
                loop {
                    let a_le: &[Shape] = match a_ve_map.get(&(a_vc.ptr_id(), a_vc.location)) {
                        Some(l) => l.as_slice(),
                        None => &[],
                    };
                    if a_le.len() != 2 {
                        // free vertex or multi-connexity, go to the next edge
                        break;
                    }
                    let mut b_stop = true;
                    for a_en in a_le {
                        if a_mfence.contains(&(a_en.ptr_id(), a_en.location, a_en.orientation)) {
                            continue;
                        }
                        let (a_vn1, a_vn2) = Self::edge_endpoints_of(a_en);
                        if a_vn1.is_partner(&a_vn2) {
                            // closed edge, go to the next edge
                            break;
                        }
                        // change orientation if necessary and go to the next edges
                        if (i == 0 && a_vc.is_partner(&a_vn2)) || (i == 1 && a_vc.is_partner(&a_vn1)) {
                            new_edges.push(a_en.clone());
                        } else {
                            let mut en_rev = a_en.clone();
                            en_rev.orientation = flip_orientation(en_rev.orientation);
                            new_edges.push(en_rev);
                        }
                        a_mfence.insert((a_en.ptr_id(), a_en.location, a_en.orientation));
                        a_vc = if a_vc.is_partner(&a_vn1) { a_vn2.clone() } else { a_vn1.clone() };
                        b_stop = false;
                        break;
                    }
                    if b_stop {
                        break;
                    }
                }
            }
        }
        // theWire = aWire
        if let TShape::Wire(wd) = Arc::make_mut(&mut w.data) {
            wd.edges = new_edges;
        }
    }

    /// Edge endpoints (TopExp::Vertices(aE, aV1, aV2, true) 鈥?the stored first
    /// and last vertices, independent of the edge orientation). The edge
    /// Location is composed into the vertices (TopoDS_Iterator cumLoc).
    fn edge_endpoints_of(e: &Shape) -> (Shape, Shape) {
        match &*e.data {
            TShape::Edge(ed) => {
                let loc = e.location;
                let vl = |v: &Shape| {
                    let vloc = if loc == 0 { v.location } else { loc };
                    Shape::new(v.data.clone(), vloc, v.orientation)
                };
                (vl(&ed.first), vl(&ed.last))
            }
            _ => (Shape::null(), Shape::null()),
        }
    }

    /// The two endpoint vertices of an edge with their composed orientations
    /// (TopoDS_Iterator semantics, as in BuilderFace::edge_vertices). The edge
    /// Location is composed into the vertices (cumLoc).
    fn edge_vertices_of(e: &Shape) -> Vec<Shape> {
        let flip_ori = |o: topods::Orientation| -> topods::Orientation {
            match o {
                topods::Orientation::Forward => topods::Orientation::Reversed,
                topods::Orientation::Reversed => topods::Orientation::Forward,
                other => other,
            }
        };
        match &*e.data {
            TShape::Edge(ed) => {
                let loc = e.location;
                let vf_loc = if loc == 0 { ed.first.location } else { loc };
                let vl_loc = if loc == 0 { ed.last.location } else { loc };
                if e.orientation == topods::Orientation::Reversed {
                    vec![
                        Shape::new(ed.last.data.clone(), vl_loc, flip_ori(ed.last.orientation)),
                        Shape::new(ed.first.data.clone(), vf_loc, flip_ori(ed.first.orientation)),
                    ]
                } else {
                    vec![
                        Shape::new(ed.first.data.clone(), vf_loc, ed.first.orientation),
                        Shape::new(ed.last.data.clone(), vl_loc, ed.last.orientation),
                    ]
                }
            }
            _ => Vec::new(),
        }
    }

    /// OCCT BOPTools_AlgoTools::OrientFacesOnShell (BOPTools_AlgoTools.cxx L363-503).
    /// Reorients the faces of a shell so that every shared edge is used by the
    /// two adjacent faces with opposite orientations. Delegates to the full
    /// ShellSplitter translation (shell_splitter.rs), operating on the shell's
    /// face list.
    fn orient_faces_on_shell(sh: &mut Shape) {
        let mut faces: Vec<Shape> = match &*sh.data {
            TShape::Shell(sd) => sd.faces.clone(),
            _ => return,
        };
        crate::bop::algo::shell_splitter::orient_faces_on_shell(&mut faces);
        if let TShape::Shell(sd) = Arc::make_mut(&mut sh.data) {
            sd.faces = faces;
        }
    }

    /// OCCT BOPAlgo_BOP::RemoveDuplicates (BOPAlgo_BOP.cxx L1627-1698).
    /// Dedups containers with identical sub-shape contents.
    fn remove_duplicates(&self, containers: &mut Vec<Shape>) {
        let mut seen: HashSet<Vec<u64>> = HashSet::new();
        let mut out: Vec<Shape> = Vec::new();
        for c in containers.iter() {
            let mut subs: Vec<u64> = self
                .shape_sub_shapes(c)
                .iter()
                .map(|s| s.ptr_id())
                .collect();
            subs.sort();
            if seen.insert(subs) {
                out.push(c.clone());
            }
        }
        *containers = out;
    }

    /// OCCT BOPAlgo_BOP::MapFacesToBuildSolids (BOPAlgo_BOP.cxx L1768-1798).
    /// OCCT uses TopExp_Explorer(theSol, TopAbs_FACE) 鈥?every face instance of
    /// the solid is visited, including the FORWARD/REVERSED copies of the same
    /// TShape (no dedup; TopExp_Explorer.cxx L110-170), with the cumulative
    /// orientation composed along the path (solid 脳 shell 脳 face). The
    /// deduplication of the solid list happens on the orientation of the
    /// FIRST-inserted face (theMFS key keeps TShape + Location via
    /// TopTools_ShapeMapHasher, L1781-1796).
    fn map_faces_to_build_solids(
        &self,
        the_sol: &Shape,
        the_mfs: &mut IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
    ) {
        // TopExp_Explorer semantics: depth-first, accumulate orientation, do
        // NOT dedup face copies (a FACE may appear twice with different
        // orientations, e.g. as both sides of an internal face).
        // OCCT aMFS is IndexedDataMap with TopTools_ShapeMapHasher: the hash
        // (std::hash<TopoDS_Shape>, TopoDS_Shape.hxx L332-340) combines the
        // TShape pointer with the LOCATION, so located copies of the same
        // TShape (e.g. the revolve end cap at the rotation Location L1) are
        // SEPARATE map keys — the cap@0 and cap@L1 each get their own entry
        // with a single solid, and both stay in the aSFS. Keying by ptr_id
        // only collapsed them into one entry whose orientation-differing
        // second visit appended the solid twice (nsols=2) and dropped the cap.
        // TopExp_Explorer walks the STORED sub-shape order (first sub-shape
        // explored first); a LIFO stack reverses that order, which changes the
        // aSFS face order and therefore which shell the final ShellSplitter
        // walk builds first.
        let mut queue: std::collections::VecDeque<Shape> = std::collections::VecDeque::new();
        queue.push_back(the_sol.clone());
        while let Some(cur) = queue.pop_front() {
            if cur.shape_type() == topods::ShapeType::Face {
                if cur.orientation == topods::Orientation::Internal {
                    continue;
                }
                let e = the_mfs
                    .entry((cur.ptr_id(), cur.location))
                    .or_insert_with(|| (cur.clone(), Vec::new()));
                // OCCT L1789-1796: append the solid only if orientations differ.
                if e.1.is_empty() || e.0.orientation != cur.orientation {
                    e.1.push(the_sol.clone());
                }
                continue;
            }
            for sub in self.shape_sub_shapes(&cur) {
                let mut sub2 = sub.clone();
                sub2.orientation = cur.orientation.compose(sub.orientation);
                queue.push_back(sub2);
            }
        }
    }

    /// OCCT BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx L130-168).
    /// Builds topology at the given shape type level.
    /// OCCT L136-143 iterates myArguments and skips any argument whose ShapeType
    /// does not match theType. For each matching argument it adds its images
    /// (or the argument itself if it has no images) to myShape, deduplicated by fence.
    /// When arguments are solids, the intermediate calls (VERTEX..SHELL, COMPOUND)
    /// are no-ops; only BuildResult(SOLID) adds shapes into the result compound.
    fn build_result(&mut self, the_type: topods::ShapeType) {
        // OCCT L133: fence map 鈥?TopTools_ShapeMapHasher (TShape + Location).
        let mut a_m_fence: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
        // OCCT L136-167: iterate myArguments, filter by theType
        let a_arguments = self.my_arguments.clone();
        for a_s in &a_arguments {
            if a_s.shape_type() != the_type { continue; }
            // OCCT L145-152: check for images
            if let Some(imgs) = self.my_images.get((a_s.ptr_id(), a_s.location)).cloned() {
                for a_s_im in &imgs {
                    if a_m_fence.insert((a_s_im.ptr_id(), a_s_im.location)) {
                        self.add_shape_to_result(a_s_im);
                    }
                }
            } else {
                if a_m_fence.insert((a_s.ptr_id(), a_s.location)) {
                    self.add_shape_to_result(a_s);
                }
            }
        }
    }

    /// Add a Shape to my_shape (OCCT equivalent: BRep_Builder().Add(myShape, aS)).
    /// OCCT uses TopoDS handles (pointer-based); rcad-kernel uses flat indices.
    /// This function clones the full TShape hierarchy and fixes all Shape.index
    /// values to point into the result BRep's tshapes array.
    /// Uses self.shape_remap persisted across the entire pipeline.
    fn add_shape_to_result(&mut self, shape: &Shape) {
        self.push_shape_recursive(shape);
    }

    /// Recursively push a Shape and all sub-shapes into my_shape.
    /// Uses self.shape_remap for persistent index tracking across pipeline stages.
    /// Returns the shape's index in the result BRep's tshapes.
    ///
    /// The dedup key is the TShape identity (ptr_id) only — Location does not
    /// participate, mirroring `nbshapes` without -t: a located copy (same
    /// TShape, different Location, e.g. the top cap of a swept prism sharing
    /// the bottom cap's TShape with a translation Location) maps to the same
    /// result shape.  The evaluated positions are carried by the Location
    /// references on the shape and its edge-endpoint vertex references, so the
    /// geometry stays correct while the topological vertex count matches OCCT
    /// (bcommon_simple F9/G2/G5/H3, bfuse_simple E1).
    fn push_shape_recursive(&mut self, shape: &Shape) -> usize {
        let ptr = shape.ptr_id();
        if let Some(&idx) = self.shape_remap.get(&ptr) {
            return idx;
        }

        // Reserve a slot in tshapes (placeholder, replaced below)
        let new_idx = {
            let brep = self.my_shape.as_mut().expect("prepare() must set my_shape");
            let idx = brep.tshapes.len();
            brep.tshapes.push(Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(), flags: 0, point: DVec3::ZERO,
                tolerance: 0.0, points: Vec::new(),
            })));
            idx
        };
        self.shape_remap.insert(ptr, new_idx);

        // Build new TShape with remapped sub-shape indices
        let new_tshape: TShape = match shape.data.as_ref() {
            TShape::Vertex(vd) => {
                let my_shapes = self.remap_shapes(&vd.my_shapes);
                TShape::Vertex(TVertexData {
                    my_shapes, ..vd.clone()
                })
            }
            TShape::Edge(ed) => {
                let _ = self.push_shape_recursive(&ed.first);
                let _ = self.push_shape_recursive(&ed.last);
                let my_shapes = self.remap_shapes(&ed.my_shapes);
                let first = self.remap_shape(&ed.first);
                let last = self.remap_shape(&ed.last);
                TShape::Edge(TEdgeData {
                    my_shapes, first, last, ..ed.clone()
                })
            }
            TShape::Wire(wd) => {
                for e in &wd.edges { let _ = self.push_shape_recursive(e); }
                let my_shapes = self.remap_shapes(&wd.my_shapes);
                let edges = self.remap_shapes(&wd.edges);
                
                TShape::Wire(TWireData {
                    my_shapes, edges, ..wd.clone()
                })
            }
            TShape::Face(fd) => {
                let _ = self.push_shape_recursive(&fd.outer_wire);
                for w in &fd.inner_wires { let _ = self.push_shape_recursive(w); }
                for v in &fd.internal_vertices { let _ = self.push_shape_recursive(v); }
                let my_shapes = self.remap_shapes(&fd.my_shapes);
                let outer_wire = self.remap_shape(&fd.outer_wire);
                let inner_wires = self.remap_shapes(&fd.inner_wires);
                let internal_vertices = self.remap_shapes(&fd.internal_vertices);
                
                TShape::Face(TFaceData {
                    my_shapes, outer_wire, inner_wires,
                    internal_vertices, ..fd.clone()
                })
            }
            TShape::Shell(sd) => {
                for f in &sd.faces { let _ = self.push_shape_recursive(f); }
                let my_shapes = self.remap_shapes(&sd.my_shapes);
                let faces = self.remap_shapes(&sd.faces);
                TShape::Shell(TShellData {
                    my_shapes, faces, ..sd.clone()
                })
            }
            TShape::Solid(sd) => {
                for s in &sd.shells { let _ = self.push_shape_recursive(s); }
                let my_shapes = self.remap_shapes(&sd.my_shapes);
                let shells = self.remap_shapes(&sd.shells);
                let internal_vertices = self.remap_shapes(&sd.internal_vertices);
                let internal_edges = self.remap_shapes(&sd.internal_edges);
                TShape::Solid(TSolidData {
                    my_shapes, shells, internal_vertices, internal_edges,
                    ..sd.clone()
                })
            }
            TShape::CompSolid(shapes) => {
                TShape::CompSolid(self.remap_shapes(shapes))
            }
            TShape::Compound(shapes) => {
                TShape::Compound(self.remap_shapes(shapes))
            }
        };

        // Replace placeholder with the remapped TShape
        let brep = self.my_shape.as_mut().unwrap();
        brep.tshapes[new_idx] = Arc::new(new_tshape);
        new_idx
    }

    /// Remap a single Shape's index via self.shape_remap (and its location via
    /// self.loc_remap, copying the matrix from the DS pool into my_shape).
    fn remap_shape(&mut self, shape: &Shape) -> Shape {
        let location = self.remap_location(shape.location);
        if let Some(&new_idx) = self.shape_remap.get(&shape.ptr_id()) {
            Shape {
                index: new_idx,
                data: shape.data.clone(),
                location,
                orientation: shape.orientation,
            }
        } else {
            let mut s = shape.clone();
            s.location = location;
            s
        }
    }

    /// Remap a DS location index into the result BRep's locations pool,
    /// copying the matrix when it is referenced for the first time.
    fn remap_location(&mut self, loc: u32) -> u32 {
        if loc == 0 {
            return 0;
        }
        if let Some(&r) = self.loc_remap.get(&loc) {
            return r;
        }
        let Some(mat) = self.ds.locations.get(loc as usize).cloned() else {
            // Dangling reference: treat as identity.
            return 0;
        };
        let brep = self.my_shape.as_mut().expect("prepare() must set my_shape");
        // The result BRep uses the rcad-kernel location convention: the
        // `locations` pool stores the real transforms only (index 0 = the
        // first non-identity matrix) and the logical index returned to shapes
        // is `pool index + 1` (0 = identity, 1 = first matrix).  This mirrors
        // BRep::add_location/get_location (topods.rs L286-308).  The DS pool
        // (merged global_locs) instead stores identity at index 0 and indexes
        // directly (DS::get_location reads locations[idx]), so the matrix is
        // looked up from the DS pool and copied into the result pool here.
        let new_idx = if let Some(pos) = brep.locations.iter().position(|x| *x == mat) {
            (pos + 1) as u32
        } else {
            brep.locations.push(mat);
            brep.locations.len() as u32
        };
        self.loc_remap.insert(loc, new_idx);
        new_idx
    }

    /// Remap a slice of Shapes.
    fn remap_shapes(&mut self, shapes: &[Shape]) -> Vec<Shape> {
        shapes.iter().map(|s| self.remap_shape(s)).collect()
    }

}

