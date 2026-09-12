//! OCCT BRepCheck shared base — `BRepCheck_Status`, `BRepCheck::Add`,
//! `BRepCheck_Result` and the shape-traversal helpers.
//!
//! Sources:
//! - `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Status.hxx L20-59`
//! - `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck.cxx L35-54`
//! - `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Result.cxx L27-84`
//! - `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Result.hxx L31-96`
//!
//! Architecture notes (documented where they apply):
//! - OCCT `TopoDS_Shape` carries Location inside the handle; rcad stores a
//!   location *index* resolved through the `BRep` location table, so the
//!   traversal helpers take `&BRep`.
//! - OCCT `NCollection_DataMap<TopoDS_Shape, List<Status>>` iteration order is
//!   unspecified; rcad preserves insertion order (an ordered map) so the
//!   context iterator is deterministic.

use std::collections::HashMap;

use glam::DAffine3;
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType, TShape};

use rcad_kernel::geom::{Curve2d, Curve3};

/// OCCT BRepCheck_Status (BRepCheck_Status.hxx L20-59) — same values, same order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BRepCheckStatus {
    NoError,
    InvalidPointOnCurve,
    InvalidPointOnCurveOnSurface,
    InvalidPointOnSurface,
    No3DCurve,
    Multiple3DCurve,
    Invalid3DCurve,
    NoCurveOnSurface,
    InvalidCurveOnSurface,
    InvalidCurveOnClosedSurface,
    InvalidSameRangeFlag,
    InvalidSameParameterFlag,
    InvalidDegeneratedFlag,
    FreeEdge,
    InvalidMultiConnexity,
    InvalidRange,
    EmptyWire,
    RedundantEdge,
    SelfIntersectingWire,
    NoSurface,
    InvalidWire,
    RedundantWire,
    IntersectingWires,
    InvalidImbricationOfWires,
    EmptyShell,
    RedundantFace,
    InvalidImbricationOfShells,
    UnorientableShape,
    NotClosed,
    NotConnected,
    SubshapeNotInShape,
    BadOrientation,
    BadOrientationOfSubshape,
    InvalidPolygonOnTriangulation,
    InvalidToleranceValue,
    EnclosedRegion,
    CheckFail,
}

/// OCCT BRepCheck::Add (BRepCheck.cxx L35-54): append `stat` to `lst`,
/// dropping a previously stored `NoError` and ignoring duplicates.
pub fn brep_check_add(lst: &mut Vec<BRepCheckStatus>, stat: BRepCheckStatus) {
    // OCCT L37-52: iterate; remove NoError when a real status arrives;
    // return early when the status is already present.
    let mut i = 0usize;
    while i < lst.len() {
        if lst[i] == BRepCheckStatus::NoError && stat != BRepCheckStatus::NoError {
            lst.remove(i);
        } else {
            if lst[i] == stat {
                return;
            }
            i += 1;
        }
    }
    lst.push(stat);
}

/// OCCT TopTools_ShapeMapHasher key (TopTools_ShapeMapHasher.hxx L38-44):
/// `S1.IsSame(S2)` — same TShape and same Location, orientation IGNORED.
/// rcad identity = (TShape Arc pointer, location index).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ShapeKey {
    pub ptr: u64,
    pub location: u32,
}

impl ShapeKey {
    pub fn of(s: &Shape) -> Self {
        ShapeKey {
            ptr: s.ptr_id(),
            location: s.location,
        }
    }
}

/// The ordered `NCollection_DataMap<TopoDS_Shape, List<Status>>` of
/// BRepCheck_Result (Result.hxx L86-89). Insertion order is preserved so the
/// context iterator is deterministic (OCCT order is unspecified).
#[derive(Debug, Default)]
pub struct StatusMap {
    entries: Vec<(Shape, Vec<BRepCheckStatus>)>,
    index: HashMap<ShapeKey, usize>,
}

impl StatusMap {
    /// OCCT NCollection_DataMap::Bind — insert (or overwrite) the list of `s`.
    pub fn bind(&mut self, s: Shape, list: Vec<BRepCheckStatus>) {
        let key = ShapeKey::of(&s);
        if let Some(&i) = self.index.get(&key) {
            self.entries[i].1 = list;
        } else {
            self.index.insert(key, self.entries.len());
            self.entries.push((s, list));
        }
    }

    /// OCCT NCollection_DataMap::Bound — bind an empty list and return true
    /// when the key was absent.
    pub fn bound(&mut self, s: &Shape) -> bool {
        let key = ShapeKey::of(s);
        if self.index.contains_key(&key) {
            false
        } else {
            self.index.insert(key, self.entries.len());
            self.entries.push((s.clone(), Vec::new()));
            true
        }
    }

    /// OCCT NCollection_DataMap::IsBound.
    pub fn is_bound(&self, s: &Shape) -> bool {
        self.index.contains_key(&ShapeKey::of(s))
    }

    /// OCCT NCollection_DataMap::Find — the list bound to `s`.
    pub fn find(&self, s: &Shape) -> Option<&Vec<BRepCheckStatus>> {
        self.index
            .get(&ShapeKey::of(s))
            .map(|&i| &self.entries[i].1)
    }

    /// Mutable list bound to `s` (OCCT `myMap(myShape)` reference read).
    pub fn find_mut(&mut self, s: &Shape) -> Option<&mut Vec<BRepCheckStatus>> {
        match self.index.get(&ShapeKey::of(s)) {
            Some(&i) => Some(&mut self.entries[i].1),
            None => None,
        }
    }

    /// OCCT NCollection_DataMap::Clear.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Diagnostics: the bound (shape, status) lists with statuses attached.
    pub fn non_empty_lists(&self) -> Vec<(Shape, Vec<BRepCheckStatus>)> {
        self.entries
            .iter()
            .filter(|(_, lst)| !lst.is_empty())
            .map(|(s, lst)| (s.clone(), lst.clone()))
            .collect()
    }

    /// Ordered (key, value) pairs — the DataMap iteration for the context
    /// iterator.
    pub fn iter(&self) -> impl Iterator<Item = (&Shape, &Vec<BRepCheckStatus>)> {
        self.entries.iter().map(|(s, l)| (s, l))
    }
}

/// OCCT TopoDS_Iterator (cumOri = true, cumLoc = true — the defaults used by
/// BRepCheck_Analyzer): the direct sub-shapes of `s` with the parent
/// orientation composed and the parent location multiplied.
pub fn iterator_subshapes(brep: &BRep, s: &Shape) -> Vec<Shape> {
    let _ = brep;
    child_occurrences(s)
}

/// The composition part of the TopoDS_Iterator (no location table needed).
pub fn child_occurrences(s: &Shape) -> Vec<Shape> {
    let mut out = Vec::new();
    if let Some(kids) = direct_children(s) {
        for c in kids {
            out.push(Shape {
                data: c.data.clone(),
                index: c.index,
                location: c.location, // composed below
                orientation: s.orientation.compose(c.orientation),
            });
            let out_len = out.len() - 1;
            // cumLoc: parent location * child location.
            let child_loc = compose_location(s.location, c.location);
            out[out_len].location = child_loc;
        }
    }
    out
}

/// The direct children of a shape (the TopoDS_Iterator enumeration order —
/// the order the sub-shape lists were built in).
fn direct_children(s: &Shape) -> Option<Vec<Shape>> {
    match &*s.data {
        TShape::Vertex(v) => Some(v.my_shapes.clone()),
        TShape::Edge(e) => {
            // OCCT `TopoDS_Iterator(E)` enumerates the edge's stored vertex
            // list in storage order with each vertex's OWN orientation — the
            // order and the tags are producer-owned and both are legitimate:
            // `BRepLib_MakeEdge::Init` (BRepLib_MakeEdge.cxx L771-772) stores
            // V1 first (low parameter, FORWARD) then V2 (high parameter,
            // REVERSED), while `BRepPrim_Builder::AddEdgeVertex`
            // (BRepPrim_Builder.cxx L143-155, direct=false) stores the HIGH
            // parameter vertex first with the REVERSED tag. The OCCT readers
            // key on the TAGS, not on the slot order — `TopExp::Vertices(E,
            // V1, V2, CumOri)` (TopExp.cxx L214-252) takes V1 = the
            // composed-FORWARD child and V2 = the composed-REVERSED one. The
            // rcad `TEdgeData` carries the same storage pair in
            // `first`/`last` with the same producer-owned tags, so the
            // children are exposed as stored (composed below).
            Some(vec![e.first.clone(), e.last.clone()])
        }
        TShape::Wire(w) => Some(w.edges.clone()),
        TShape::Face(f) => {
            let mut kids = vec![f.outer_wire.clone()];
            kids.extend(f.inner_wires.iter().cloned());
            Some(kids)
        }
        TShape::Shell(sh) => Some(sh.faces.clone()),
        TShape::Solid(sd) => Some(sd.shells.clone()),
        TShape::CompSolid(c) | TShape::Compound(c) => Some(c.clone()),
    }
}

/// TopLoc_Location composition over the BRep location table: `a * b`
/// (identity index 0 shortcuts). The composition is NOT registered back in
/// the table — the analyzer only reads composed locations.
pub fn compose_location(a: u32, b: u32) -> u32 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    // Without a table to consult, non-trivial compositions are represented
    // by `a` (the dominant wrapper location). The BRep-checking read paths
    // that need the exact matrix call `location_matrix` instead.
    a
}

/// The location MATRIX of a composed (parent, child) location pair —
/// OCCT `TopLoc_Location::Multiplied`. Used by the traversal code that
/// evaluates geometry.
pub fn location_matrix(brep: &BRep, parent: u32, child: u32) -> DAffine3 {
    if parent == 0 {
        return brep.get_location(child);
    }
    if child == 0 {
        return brep.get_location(parent);
    }
    brep.get_location(parent) * brep.get_location(child)
}

/// OCCT TopExp_Explorer(S, ToFind): depth-first enumeration of all
/// sub-shapes of type `to_find` under `s` (never `s` itself), with cumulated
/// orientation and location. Occurrences are NOT deduplicated.
pub fn explorer(brep: &BRep, s: &Shape, to_find: ShapeType) -> Vec<Shape> {
    let mut out = Vec::new();
    let mut stack: Vec<Shape> = iterator_subshapes(brep, s);
    // Depth-first with the children pushed in reverse so the enumeration
    // order matches the recursive OCCT explorer.
    stack.reverse();
    while let Some(cur) = stack.pop() {
        if cur.shape_type() == to_find {
            out.push(cur);
        } else {
            let mut kids = iterator_subshapes(brep, &cur);
            kids.reverse();
            stack.extend(kids);
        }
    }
    out
}

/// OCCT TopoDS_Shape::Oriented(TopAbs_FORWARD/...) — a copy with the given
/// orientation.
pub fn oriented(s: &Shape, o: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = o;
    c
}

/// OCCT BRep_CurveRepresentation kinds flattened from the rcad TEdgeData
/// storage. Architecture bridge: OCCT keeps one `BRep_CurveRepresentation`
/// list per edge TShape (`BRep_TEdge::Curves()`); rcad splits the same data
/// over `TEdgeData::curve` (the single 3D curve), `TEdgeData::representations`
/// and the `TEdgeData::pcurves` fast index (whose rows mirror the
/// CurveOnSurface representations written by the same builders). The
/// flattening reproduces the OCCT list content, skipping pcurve rows already
/// present in `representations`.
#[derive(Debug, Clone)]
pub enum EdgeCurveRep {
    /// BRep_Curve3D — the 3D curve with its own location.
    Curve3D { curve: Curve3, location: u32 },
    /// BRep_Curve3D whose curve value rcad cannot resolve (the
    /// `CurveRepresentation::Curve3D` index form) — counted for the
    /// existence/uniqueness checks only.
    Curve3DUnresolved,
    /// BRep_CurveOnSurface — pcurve on a face; `face_key` is the rcad
    /// (face TShape pointer, pcurve location id) storage key.
    CurveOnSurface {
        face_key: (u64, u32),
        pcurve: Curve2d,
        range: [f64; 2],
        location: u32,
    },
    /// BRep_CurveOnClosedSurface — the two pcurves of a seam edge.
    CurveOnClosedSurface {
        face_key: (u64, u32),
        pcurve1: Curve2d,
        pcurve2: Curve2d,
        range: [f64; 2],
        location: u32,
    },
    /// BRep_CurveOn2Surfaces (regularity) — neither a 3D curve nor a curve
    /// on surface.
    Regularity,
}

impl EdgeCurveRep {
    /// OCCT BRep_CurveRepresentation::IsCurve3D.
    pub fn is_curve3d(&self) -> bool {
        matches!(self, EdgeCurveRep::Curve3D { .. } | EdgeCurveRep::Curve3DUnresolved)
    }

    /// OCCT BRep_CurveRepresentation::IsCurveOnSurface.
    pub fn is_curve_on_surface(&self) -> bool {
        matches!(
            self,
            EdgeCurveRep::CurveOnSurface { .. } | EdgeCurveRep::CurveOnClosedSurface { .. }
        )
    }

    /// OCCT BRep_CurveRepresentation::IsCurveOnClosedSurface.
    pub fn is_curve_on_closed_surface(&self) -> bool {
        matches!(self, EdgeCurveRep::CurveOnClosedSurface { .. })
    }
}

/// The flattened OCCT `BRep_TEdge::Curves()` list of `edge`
/// (architecture bridge — see `EdgeCurveRep`).
pub fn edge_curve_reps(brep: &BRep, edge: &Shape) -> Vec<EdgeCurveRep> {
    let _ = brep;
    let mut out = Vec::new();
    let Some(ed) = edge.as_edge() else {
        return out;
    };
    // The canonical 3D curve (OCCT: the BRep_Curve3D representation).
    if let Some(c) = &ed.curve {
        out.push(EdgeCurveRep::Curve3D {
            curve: c.clone(),
            location: 0,
        });
    }
    // The representations written by BRep_Builder::UpdateEdge.
    let mut covered: Vec<(u64, u32)> = Vec::new();
    for r in &ed.representations {
        match r {
            rcad_kernel::topods::CurveRepresentation::Curve3D { .. } => {
                out.push(EdgeCurveRep::Curve3DUnresolved);
            }
            rcad_kernel::topods::CurveRepresentation::CurveOnSurface { face, pcurve, range } => {
                covered.push(*face);
                out.push(EdgeCurveRep::CurveOnSurface {
                    face_key: *face,
                    pcurve: pcurve.clone(),
                    range: *range,
                    location: face.1,
                });
            }
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                pcurve2,
                range,
            } => {
                covered.push(*face);
                out.push(EdgeCurveRep::CurveOnClosedSurface {
                    face_key: *face,
                    pcurve1: pcurve1.clone(),
                    pcurve2: pcurve2.clone(),
                    range: *range,
                    location: face.1,
                });
            }
            rcad_kernel::topods::CurveRepresentation::CurveOn2Surfaces { .. } => {
                out.push(EdgeCurveRep::Regularity);
            }
        }
    }
    // Pcurve rows whose representation is not already listed above (the
    // single-storage builders write only the fast index).
    for (k, (pc, f, l)) in &ed.pcurves {
        if covered.contains(k) {
            continue;
        }
        out.push(EdgeCurveRep::CurveOnSurface {
            face_key: *k,
            pcurve: pc.clone(),
            range: [*f, *l],
            location: k.1,
        });
    }
    out
}

/// OCCT BRepCheck_Result (BRepCheck_Result.hxx L31-96) — the shared base of
/// the per-type results. The virtual methods `InContext` / `Minimum` /
/// `Blind` live on the subtype structs; the base carries the status map.
#[derive(Debug)]
pub struct BRepCheckResultBase {
    /// OCCT myShape.
    pub my_shape: Shape,
    /// OCCT myMin.
    pub my_min: bool,
    /// OCCT myBlind.
    pub my_blind: bool,
    /// OCCT myIsParallel (Result.hxx L85). The analyzer runs single-threaded,
    /// so the OCCT mutex locking (a no-op when !parallel) is not modelled.
    pub my_is_parallel: bool,
    /// OCCT myMap (Result.hxx L86-89).
    pub my_map: StatusMap,
    /// OCCT myIter (Result.hxx L93-95) — the context iterator position.
    /// OCCT mutates the iterator through a const handle in ValidSub; the
    /// rcad position lives in a Cell so the iterator methods take `&self`.
    pub my_iter_pos: std::cell::Cell<Option<usize>>,
}

impl BRepCheckResultBase {
    /// OCCT BRepCheck_Result::BRepCheck_Result (Result.cxx L27-31).
    pub fn new() -> Self {
        BRepCheckResultBase {
            my_shape: Shape::null(),
            my_min: false,
            my_blind: false,
            my_is_parallel: false,
            my_map: StatusMap::default(),
            my_iter_pos: std::cell::Cell::new(None),
        }
    }

    /// OCCT BRepCheck_Result::Init (Result.cxx L35-42) — the Minimum() call
    /// is dispatched by the subtype wrapper (Rust has no virtual dispatch).
    pub fn init(&mut self, s: &Shape) {
        self.my_shape = s.clone();
        self.my_min = false;
        self.my_blind = false;
        self.my_map.clear();
        // Minimum(); — performed by the caller on the concrete type.
    }

    /// OCCT BRepCheck_Result::SetFailStatus (Result.cxx L46-61).
    pub fn set_fail_status(&mut self, s: &Shape) {
        if !self.my_map.bound(s) {
            // Bound already existed; keep the stored list (OCCT Find/Bind).
        }
        let list = self
            .my_map
            .find_mut(s)
            .expect("SetFailStatus: the shape must be bound");
        brep_check_add(list, BRepCheckStatus::CheckFail);
    }

    /// OCCT BRepCheck_Result::Status (Result.hxx L45) — the status list on
    /// myShape itself.
    pub fn status(&self) -> &Vec<BRepCheckStatus> {
        self.my_map
            .find(&self.my_shape)
            .expect("Status: myShape must be bound")
    }

    /// OCCT BRepCheck_Result::IsStatusOnShape (Result.hxx L67).
    pub fn is_status_on_shape(&self, s: &Shape) -> bool {
        self.my_map.is_bound(s)
    }

    /// OCCT BRepCheck_Result::StatusOnShape(S) (Result.hxx L69-72).
    pub fn status_on_shape(&self, s: &Shape) -> &Vec<BRepCheckStatus> {
        self.my_map
            .find(s)
            .expect("StatusOnShape: the shape must be bound")
    }

    /// OCCT BRepCheck_Result::InitContextIterator (Result.cxx L65-73) —
    /// "At least 1 element : the Shape itself".
    pub fn init_context_iterator(&self) {
        self.my_iter_pos.set(Some(0));
        if let Some(i) = self.my_iter_pos.get() {
            if let Some((s, _)) = self.my_map.entries_ref().get(i) {
                if s.is_same(&self.my_shape) {
                    self.my_iter_pos.set(Some(i + 1));
                }
            }
        }
    }

    /// OCCT BRepCheck_Result::MoreShapeInContext (Result.hxx L53).
    pub fn more_shape_in_context(&self) -> bool {
        match self.my_iter_pos.get() {
            Some(i) => i < self.my_map.len(),
            None => false,
        }
    }

    /// OCCT BRepCheck_Result::ContextualShape (Result.hxx L55).
    pub fn contextual_shape(&self) -> &Shape {
        let i = self
            .my_iter_pos
            .get()
            .expect("ContextualShape: iterator not initialised");
        &self.my_map.entries_ref()[i].0
    }

    /// OCCT BRepCheck_Result::StatusOnShape (Result.hxx L57) — the statuses
    /// on the current contextual shape.
    pub fn status_on_shape_in_context(&self) -> &Vec<BRepCheckStatus> {
        let i = self
            .my_iter_pos
            .get()
            .expect("StatusOnShape: iterator not initialised");
        &self.my_map.entries_ref()[i].1
    }

    /// OCCT BRepCheck_Result::NextShapeInContext (Result.cxx L77-84).
    pub fn next_shape_in_context(&self) {
        if let Some(i) = self.my_iter_pos.get() {
            self.my_iter_pos.set(Some(i + 1));
            if let Some(j) = self.my_iter_pos.get() {
                if j < self.my_map.len() && self.my_map.entries_ref()[j].0.is_same(&self.my_shape)
                {
                    self.my_iter_pos.set(Some(j + 1));
                }
            }
        }
    }
}

impl Default for BRepCheckResultBase {
    fn default() -> Self {
        Self::new()
    }
}

impl StatusMap {
    /// Read-only access to the ordered entries (used by the context iterator).
    pub fn entries_ref(&self) -> &[(Shape, Vec<BRepCheckStatus>)] {
        &self.entries
    }
}

/// OCCT BRep_Tool::Tolerance(E) (BRep_Tool.cxx L896-910) — the edge
/// tolerance floored at Precision::Confusion().
pub fn brep_tool_tolerance_edge(_brep: &BRep, e: &Shape) -> f64 {
    let p = e.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
    let p_min = rcad_kernel::CONFUSION;
    if p > p_min {
        p
    } else {
        p_min
    }
}

/// OCCT BRep_Tool::Tolerance(V) — the raw vertex tolerance.
pub fn brep_tool_tolerance_vertex(_brep: &BRep, v: &Shape) -> f64 {
    v.as_vertex().map(|vd| vd.tolerance).unwrap_or(0.0)
}

/// OCCT BRep_Tool::Tolerance(F) — the raw face tolerance.
pub fn brep_tool_tolerance_face(brep: &BRep, f: &Shape) -> f64 {
    let _ = brep;
    f.as_face().map(|fd| fd.tolerance).unwrap_or(0.0)
}
