//! OCCT TopOpeBRepBuild loop / classification machinery behind the
//! SplitFace1 face reconstruction (continuation file of
//! `fillet::hbuilder_face`, D6 carrier of the TKBool TopOpeBRepBuild
//! classes as components of the TopOpeBRepBuild_HBuilder translation):
//!
//! - `TopOpeBRepBuild_BlockIterator` (BlockIterator.cxx L21-34 + .lxx L17-56)
//! - `TopOpeBRepBuild_Loop` (Loop.cxx L23-61)
//! - `TopOpeBRepBuild_LoopSet` (LoopSet.cxx L23-66)
//! - `TopOpeBRepBuild_ShapeSet` (ShapeSet.cxx L36-502) + the subclass
//!   interface the BlockBuilder consumes
//! - `TopOpeBRepBuild_BlockBuilder` (BlockBuilder.cxx L23-285)
//! - `TopOpeBRepBuild_LoopClassifier` (LoopClassifier.hxx L27-41)
//! - `TopOpeBRepBuild_AreaBuilder` core (AreaBuilder.cxx L30-429) and the
//!   `TopOpeBRepBuild_Area2dBuilder::InitAreaBuilder` override
//!   (Area2dBuilder.cxx L41-262) — the operative virtual for the
//!   `TopOpeBRepBuild_FaceAreaBuilder` (FaceAreaBuilder.cxx L23-41)
//! - `TopOpeBRepBuild_CompositeClassifier` base state
//!   (CompositeClassifier.cxx L44-48)
//! - `TopOpeBRepBuild_WireEdgeClassifier` (WireEdgeClassifier.cxx L52-527)
//!
//! Architecture differences (Rust <-> C++), per the D6 carrier notes:
//! - `NCollection_List<TopoDS_Shape>` -> `Vec<Shape>`.
//! - `occ::handle<TopOpeBRepBuild_Loop>` -> `Arc<Loop>` (the handle is
//!   shared across the area / boundaryloops lists — the same aliasing
//!   semantics).
//! - `NCollection_IndexedDataMap<TopoDS_Shape, ...>` keyed by the
//!   `TopTools_ShapeMapHasher` identity (TShape + Location, orientation
//!   ignored) -> `HashMap<(u64, u32), ...>` (the same key the hbuilder
//!   split / merged tables use).
//! - The list iterators become index cursors over the owning `Vec`.
//! - OCCT classifiers read TShapes through handles; rcad's translated
//!   2D-classifier state machine (`FClass2dOfFClassifier`, the
//!   byte-identical twin of `BRepClass_FacePassiveClassifier`) is
//!   flat-index based over a `ShapeSource`, so the WireEdgeClassifier
//!   carries a `WecShapeSource` view (the classification face at index 0,
//!   the compared edges registered from 1 on).

// The carrier is unwired by design (the merge_solid wiring point belongs
// to the main agent — see the hbuilder_face module doc), so most of its
// surface is dead code until then.
#![allow(dead_code)]

use std::sync::Arc;

use glam::{DAffine3, DVec2, DVec3};

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Curve2d, Curve2dEval as _, CurveEval as _, Line2d};
use rcad_kernel::topods::{BRep, BRepTool as _, Orientation, Shape, ShapeType, State, TShape};

use crate::topalgo::brep_class::edge::ClassEdge;
use crate::topalgo::brep_class::face_classifier::{FClass2dOfFClassifier, FClassifier};
use crate::topalgo::shape_source::ShapeSource;

// =========================================================================
// Shared shape identity helpers (the TopTools_ShapeMapHasher semantics).
// =========================================================================

/// TopTools_ShapeMapHasher identity: TShape + Location (orientation
/// ignored) — the rcad `(ptr_id, location)` key (hbuilder.rs convention).
pub(crate) fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// TopoDS_Shape::IsSame — same TShape + Location.
pub(crate) fn shapes_same(a: &Shape, b: &Shape) -> bool {
    shape_key(a) == shape_key(b)
}

/// TopoDS_Shape::IsEqual — IsSame + same orientation.
pub(crate) fn shapes_equal(a: &Shape, b: &Shape) -> bool {
    shapes_same(a, b) && a.orientation == b.orientation
}

/// TopAbs::Complement (TopAbs.hxx L100-121): FORWARD<->REVERSED,
/// INTERNAL<->EXTERNAL.
pub(crate) fn top_abs_complement(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}

/// TopoDS_Shape::Oriented / the myBuildTool.Orientation(S, O) write — the
/// handle is value-copied with the new orientation (OCCT mutates the
/// handle's orientation field; rcad copies the handle).
pub(crate) fn shape_oriented(s: &Shape, o: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = o;
    c
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri = false) — the canonical
/// bounding vertices of the edge (BRep_Tool::FirstVertex / LastVertex).
pub(crate) fn top_exp_vertices(brep: &BRep, e: &Shape) -> (Shape, Shape) {
    (brep.first_vertex(e), brep.last_vertex(e))
}

/// The DS-convention location table (slot 0 = identity, kernel
/// transforms from slot 1) — the input the FaceShapeSource /
/// `edge_pcurve_on_face` family expects (shape_source.rs note: a kernel
/// BRep table must be converted first).
pub(crate) fn ds_locations(brep: &BRep) -> Vec<DAffine3> {
    let mut locs = Vec::with_capacity(brep.locations.len() + 1);
    locs.push(DAffine3::IDENTITY);
    locs.extend(brep.locations.iter().copied());
    locs
}

/// OCCT BRep_Tool::Surface(F) local-surface read (kernel free function).
pub(crate) fn face_surface_value(brep: &BRep, f: &Shape) -> Option<rcad_kernel::geom::Surface3> {
    rcad_kernel::topo::topods::face_surface_value(brep, f)
}

/// BRep_Tool::Tolerance(V) read off the TShape data.
pub(crate) fn vertex_tolerance_of(s: &Shape) -> f64 {
    match &*s.data {
        TShape::Vertex(vd) => vd.tolerance,
        _ => 0.0,
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_BlockIterator (BlockIterator.cxx L21-34 +
// BlockIterator.lxx L17-56).
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub(crate) struct BlockIterator {
    my_lower: i32,
    my_upper: i32,
    my_value: i32,
}

impl BlockIterator {
    /// BlockIterator.cxx L21-27.
    pub(crate) fn new() -> Self {
        BlockIterator { my_lower: 0, my_upper: 0, my_value: 1 }
    }

    /// BlockIterator.cxx L30-34.
    pub(crate) fn new2(lower: i32, upper: i32) -> Self {
        BlockIterator { my_lower: lower, my_upper: upper, my_value: lower }
    }

    /// BlockIterator.lxx L19-22.
    pub(crate) fn initialize(&mut self) {
        self.my_value = self.my_lower;
    }

    /// BlockIterator.lxx L26-30.
    pub(crate) fn more(&self) -> bool {
        let b = self.my_value <= self.my_upper;
        b
    }

    /// BlockIterator.lxx L34-36.
    pub(crate) fn next(&mut self) {
        self.my_value += 1;
    }

    /// BlockIterator.lxx L41-43.
    pub(crate) fn value(&self) -> i32 {
        self.my_value
    }

    /// BlockIterator.lxx L48-56.
    pub(crate) fn extent(&self) -> i32 {
        if self.my_lower != 0 {
            let n = self.my_upper - self.my_lower + 1;
            return n;
        }
        0
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Loop (Loop.cxx L23-61) — the member layout
// (myIsShape / myShape / myBlockIterator) carried as a payload selection.
// =========================================================================
#[derive(Debug, Clone)]
pub(crate) struct Loop {
    pub(crate) my_is_shape: bool,
    pub(crate) my_shape: Shape,
    pub(crate) my_block_iterator: BlockIterator,
}

impl Loop {
    /// Loop.cxx L23-28.
    pub(crate) fn new_shape(s: &Shape) -> Loop {
        Loop {
            my_is_shape: true,
            my_shape: s.clone(),
            my_block_iterator: BlockIterator::new2(0, 0),
        }
    }

    /// Loop.cxx L32-36.
    pub(crate) fn new_block(bi: BlockIterator) -> Loop {
        Loop {
            my_is_shape: false,
            my_shape: Shape::null(),
            my_block_iterator: bi,
        }
    }

    /// Loop.cxx L40-43.
    pub(crate) fn is_shape(&self) -> bool {
        self.my_is_shape
    }

    /// Loop.cxx L47-50.
    pub(crate) fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// Loop.cxx L54-57.
    pub(crate) fn block_iterator(&self) -> &BlockIterator {
        &self.my_block_iterator
    }
}

/// occ::handle<TopOpeBRepBuild_Loop>.
pub(crate) type LoopHandle = Arc<Loop>;

// =========================================================================
// OCCT TopOpeBRepBuild_LoopSet (LoopSet.cxx L23-66).
// =========================================================================
#[derive(Debug)]
pub(crate) struct LoopSet {
    my_list_of_loop: Vec<LoopHandle>,
    my_loop_index: i32,
    my_nb_loop: i32,
}

impl Default for LoopSet {
    fn default() -> Self {
        Self::new()
    }
}

impl LoopSet {
    /// LoopSet.cxx L25-29 (myLoopIndex = 1, myNbLoop = 0).
    pub(crate) fn new() -> Self {
        LoopSet { my_list_of_loop: Vec::new(), my_loop_index: 1, my_nb_loop: 0 }
    }

    /// LoopSet.cxx L36-39.
    pub(crate) fn change_list_of_loop(&mut self) -> &mut Vec<LoopHandle> {
        &mut self.my_list_of_loop
    }

    /// LoopSet.cxx L41-47.
    pub(crate) fn init_loop(&mut self) {
        self.my_loop_index = 1;
        self.my_nb_loop = self.my_list_of_loop.len() as i32;
    }

    /// LoopSet.cxx L49-52.
    pub(crate) fn more_loop(&self) -> bool {
        let b = self.my_loop_index >= 1
            && ((self.my_loop_index - 1) as usize) < self.my_list_of_loop.len();
        b
    }

    /// LoopSet.cxx L54-57.
    pub(crate) fn next_loop(&mut self) {
        self.my_loop_index += 1;
    }

    /// LoopSet.cxx L59-64.
    pub(crate) fn loop_(&self) -> LoopHandle {
        self.my_list_of_loop[(self.my_loop_index - 1) as usize].clone()
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_ShapeSet base data (ShapeSet.cxx L36-58 ctor +
// the member set) and the subclass interface (ShapeSet.hxx L68-118
// virtuals) the BlockBuilder / FaceBuilder consume.
// =========================================================================
#[derive(Debug)]
pub(crate) struct ShapeSetBase {
    /// OCCT: myShapeType / mySubShapeType (ctor L41-52).
    pub(crate) my_shape_type: ShapeType,
    pub(crate) my_sub_shape_type: ShapeType,
    /// OCCT: mySubShapeExplorer (a TopOpeBRepTool_ShapeExplorer) — the
    /// explorer content + cursor.
    pub(crate) my_sub_shape_explorer: Vec<Shape>,
    pub(crate) my_sub_shape_explorer_pos: usize,
    pub(crate) my_start_shapes: Vec<Shape>,
    pub(crate) my_start_shapes_iter: usize,
    /// OCCT: mySubShapeMap (IndexedDataMap keyed by the ShapeMapHasher
    /// identity).
    pub(crate) my_sub_shape_map: std::collections::HashMap<(u64, u32), Vec<Shape>>,
    pub(crate) my_incident_shapes: Vec<Shape>,
    pub(crate) my_incident_shapes_iter: usize,
    pub(crate) my_shapes: Vec<Shape>,
    pub(crate) my_shapes_iter: usize,
    pub(crate) my_current_shape: Shape,
    pub(crate) my_current_shape_neighbours: Vec<Shape>,
    pub(crate) my_deb_number: i32,
    /// OCCT: myOMSS / myOMES / myOMSH (IndexedMap — Contains only).
    pub(crate) my_omss: std::collections::HashSet<(u64, u32)>,
    pub(crate) my_omes: std::collections::HashSet<(u64, u32)>,
    pub(crate) my_omsh: std::collections::HashSet<(u64, u32)>,
    pub(crate) my_check_shape: bool,
}

impl ShapeSetBase {
    /// ShapeSet.cxx L36-56 (SubShapeType EDGE -> ShapeType FACE;
    /// VERTEX -> EDGE; myCheckShape = false "temporary NYI" L55).
    pub(crate) fn new(sub_shape_type: ShapeType) -> Self {
        let my_shape_type = if sub_shape_type == ShapeType::Edge {
            ShapeType::Face
        } else if sub_shape_type == ShapeType::Vertex {
            ShapeType::Edge
        } else {
            panic!("ShapeSet : bad ShapeType");
        };
        ShapeSetBase {
            my_shape_type,
            my_sub_shape_type: sub_shape_type,
            my_sub_shape_explorer: Vec::new(),
            my_sub_shape_explorer_pos: 0,
            my_start_shapes: Vec::new(),
            my_start_shapes_iter: 0,
            my_sub_shape_map: std::collections::HashMap::new(),
            my_incident_shapes: Vec::new(),
            my_incident_shapes_iter: 0,
            my_shapes: Vec::new(),
            my_shapes_iter: 0,
            my_current_shape: Shape::null(),
            my_current_shape_neighbours: Vec::new(),
            my_deb_number: 0,
            my_omss: std::collections::HashSet::new(),
            my_omes: std::collections::HashSet::new(),
            my_omsh: std::collections::HashSet::new(),
            my_check_shape: false,
        }
    }

    /// ShapeSet.cxx L110-117.
    pub(crate) fn process_add_shape(&mut self, s: &Shape) {
        if !self.my_omsh.contains(&shape_key(s)) {
            self.my_omsh.insert(shape_key(s));
            self.my_shapes.push(s.clone());
        }
    }

    /// ShapeSet.cxx L121-129.
    pub(crate) fn process_add_start_element(&mut self, brep: &BRep, s: &Shape) {
        if !self.my_omss.contains(&shape_key(s)) {
            self.my_omss.insert(shape_key(s));
            self.my_start_shapes.push(s.clone());
            self.process_add_element(brep, s);
        }
    }

    /// ShapeSet.cxx L133-151 (ProcessAddElement with the brep the
    /// subshape exploration reads).
    pub(crate) fn process_add_element(&mut self, brep: &BRep, s: &Shape) {
        if !self.my_omes.contains(&shape_key(s)) {
            self.my_omes.insert(shape_key(s));
            for subshape in shape_explore_children(brep, s, self.my_sub_shape_type) {
                let b = !self.my_sub_shape_map.contains_key(&shape_key(&subshape));
                if b {
                    self.my_sub_shape_map.insert(shape_key(&subshape), Vec::new());
                }
                self.my_sub_shape_map
                    .get_mut(&shape_key(&subshape))
                    .expect("bound above")
                    .push(s.clone());
            }
        }
    }

    /// ShapeSet.cxx L155-158.
    pub(crate) fn start_elements(&self) -> &Vec<Shape> {
        &self.my_start_shapes
    }

    /// ShapeSet.cxx L162-165.
    pub(crate) fn init_shapes(&mut self) {
        self.my_shapes_iter = 0;
    }

    /// ShapeSet.cxx L169-173.
    pub(crate) fn more_shapes(&self) -> bool {
        self.my_shapes_iter < self.my_shapes.len()
    }

    /// ShapeSet.cxx L177-180.
    pub(crate) fn next_shape(&mut self) {
        self.my_shapes_iter += 1;
    }

    /// ShapeSet.cxx L184-188.
    pub(crate) fn shape(&self) -> Shape {
        self.my_shapes[self.my_shapes_iter].clone()
    }

    /// ShapeSet.cxx L192-195.
    pub(crate) fn init_start_elements(&mut self) {
        self.my_start_shapes_iter = 0;
    }

    /// ShapeSet.cxx L199-203.
    pub(crate) fn more_start_elements(&self) -> bool {
        self.my_start_shapes_iter < self.my_start_shapes.len()
    }

    /// ShapeSet.cxx L207-210.
    pub(crate) fn next_start_element(&mut self) {
        self.my_start_shapes_iter += 1;
    }

    /// ShapeSet.cxx L214-218.
    pub(crate) fn start_element(&self) -> Shape {
        self.my_start_shapes[self.my_start_shapes_iter].clone()
    }

    /// ShapeSet.cxx L222-227 (base InitNeighbours).
    pub(crate) fn base_init_neighbours(&mut self, brep: &BRep, s: &Shape) {
        self.my_sub_shape_explorer = shape_explore_children(brep, s, self.my_sub_shape_type);
        self.my_sub_shape_explorer_pos = 0;
        self.my_current_shape = s.clone();
    }

    /// ShapeSet.cxx L231-235.
    pub(crate) fn more_neighbours(&self) -> bool {
        self.my_incident_shapes_iter < self.my_incident_shapes.len()
    }

    /// ShapeSet.cxx L256-260.
    pub(crate) fn neighbour(&self) -> Shape {
        self.my_incident_shapes[self.my_incident_shapes_iter].clone()
    }

    /// ShapeSet.cxx L264-267.
    pub(crate) fn change_start_shapes(&mut self) -> &mut Vec<Shape> {
        &mut self.my_start_shapes
    }

    /// ShapeSet.cxx L271-296 (base FindNeighbours).  The neighbours list
    /// `l` (an OCCT reference into mySubShapeMap /
    /// myCurrentShapeNeighbours) is stored into myIncidentShapes — the
    /// same data flow myIncidentShapesIter.Initialize(l) realizes.
    pub(crate) fn base_find_neighbours(&mut self, brep: &BRep) {
        while self.my_sub_shape_explorer_pos < self.my_sub_shape_explorer.len() {
            let v = self.my_sub_shape_explorer[self.my_sub_shape_explorer_pos].clone();
            let cur = self.my_current_shape.clone();
            let l = self.make_neighbours_list_impl(brep, &cur, &v);
            self.my_incident_shapes = l;
            self.my_incident_shapes_iter = 0;
            if self.more_neighbours() {
                break;
            }
            self.my_sub_shape_explorer_pos += 1;
        }
    }

    /// ShapeSet.cxx L302-308 (base MakeNeighboursList — FindFromKey).
    pub(crate) fn make_neighbours_list_impl(
        &mut self,
        _brep: &BRep,
        _earg: &Shape,
        varg: &Shape,
    ) -> Vec<Shape> {
        let l = self
            .my_sub_shape_map
            .get(&shape_key(varg))
            .cloned()
            .unwrap_or_default();
        l
    }

    /// ShapeSet.cxx L312-334.
    pub(crate) fn max_number_sub_shape(&mut self, brep: &BRep, shape: &Shape) -> i32 {
        let mut m: i32 = 0;
        for sub_shape in shape_explore_children(brep, shape, self.my_sub_shape_type) {
            if !self.my_sub_shape_map.contains_key(&shape_key(&sub_shape)) {
                continue;
            }
            let i = self
                .my_sub_shape_map
                .get(&shape_key(&sub_shape))
                .map(|l| l.len())
                .unwrap_or(0) as i32;
            m = m.max(i);
        }
        m
    }

    /// ShapeSet.cxx L360-370 — CheckShape(S) is a no-op while
    /// myCheckShape is false (L55 "temporary NYI"), so the check is
    /// always true.
    pub(crate) fn check_shape(&self, _s: &Shape) -> bool {
        if !self.my_check_shape {
            return true;
        }
        // BRepCheck_Analyzer path is unreachable on the D6 carrier
        // (myCheckShape is never set to true — ShapeSet.cxx L55).
        true
    }

    /// ShapeSet.cxx L338-349.
    pub(crate) fn set_check_shape(&mut self, checkshape: bool) {
        self.my_check_shape = checkshape;
    }
}

/// TopOpeBRepTool_ShapeExplorer over the direct subshapes this batch
/// explores (TopAbs_VERTEX of an EDGE; TopAbs_EDGE of a WIRE; TopAbs_WIRE
/// of a FACE; TopAbs_FACE of a SHELL).  OCCT TopExp_Explorer yields the
/// edge vertices with their in-edge orientation (FORWARD first vertex,
/// REVERSED last vertex); the wire / face walks yield the stored children
/// with their recorded orientations.
pub(crate) fn shape_explore_children(brep: &BRep, s: &Shape, t: ShapeType) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    match (&*s.data, t) {
        (TShape::Edge(ed), ShapeType::Vertex) => {
            out.push(shape_oriented(&ed.first, Orientation::Forward));
            out.push(shape_oriented(&ed.last, Orientation::Reversed));
        }
        (TShape::Wire(wd), ShapeType::Edge) => {
            out.extend(wd.edges.iter().cloned());
        }
        (TShape::Face(fd), ShapeType::Wire) => {
            if !fd.outer_wire.is_null() {
                out.push(fd.outer_wire.clone());
            }
            out.extend(fd.inner_wires.iter().cloned());
        }
        (TShape::Shell(sd), ShapeType::Face) => {
            out.extend(sd.faces.iter().cloned());
        }
        (TShape::Vertex(_), ShapeType::Vertex) => {
            out.push(s.clone());
        }
        // TopExp_Explorer recursion through compounds / compsolids (the
        // MapShapesAndAncestors walks over the DetectUnclosedWire /
        // DetectPseudoInternalEdge compounds).
        (TShape::Compound(cs), _) | (TShape::CompSolid(cs), _) => {
            for c in cs {
                out.extend(shape_explore_children(brep, c, t));
            }
        }
        _ => {}
    }
    let _ = brep;
    out
}

/// The OCCT ShapeSet virtual interface (ShapeSet.hxx L68-118: AddShape /
/// AddStartElement / AddElement / StartElements / the shape +
/// startelement iterators / InitNeighbours / MoreNeighbours /
/// NextNeighbour / Neighbour / FindNeighbours / MakeNeighboursList /
/// MaxNumberSubShape).
pub(crate) trait ShapeSet {
    fn add_shape(&mut self, brep: &BRep, s: &Shape);
    fn add_start_element(&mut self, brep: &BRep, s: &Shape);
    fn add_element(&mut self, brep: &BRep, s: &Shape);
    fn start_elements(&self) -> &Vec<Shape>;
    fn init_shapes(&mut self);
    fn more_shapes(&self) -> bool;
    fn next_shape(&mut self);
    fn shape(&self) -> Shape;
    fn init_start_elements(&mut self);
    fn more_start_elements(&self) -> bool;
    fn next_start_element(&mut self);
    fn start_element(&self) -> Shape;
    fn init_neighbours(&mut self, brep: &BRep, s: &Shape);
    fn more_neighbours(&self) -> bool;
    fn next_neighbour(&mut self, brep: &BRep);
    fn neighbour(&self) -> Shape;
    fn find_neighbours(&mut self, brep: &BRep);
    fn make_neighbours_list(&mut self, brep: &BRep, e: &Shape, v: &Shape) -> Vec<Shape>;
    fn max_number_sub_shape(&mut self, brep: &BRep, shape: &Shape) -> i32;
}

// =========================================================================
// OCCT TopOpeBRepBuild_BlockBuilder (BlockBuilder.cxx L23-285).
// =========================================================================
#[derive(Debug, Default)]
pub(crate) struct BlockBuilder {
    /// OCCT: myOrientedShapeMap (IndexedMap — Add returns the 1-based
    /// index, duplicates keep the existing one).
    pub(crate) my_oriented_shape_map: Vec<Shape>,
    pub(crate) my_oriented_shape_map_index: std::collections::HashMap<(u64, u32), i32>,
    /// OCCT: myOrientedShapeMapIsValid (IndexedDataMap<int, int>).
    pub(crate) my_oriented_shape_map_is_valid: std::collections::HashMap<i32, i32>,
    pub(crate) my_blocks: Vec<i32>,
    pub(crate) my_blocks_is_regular: Vec<i32>,
    pub(crate) my_is_done: bool,
    pub(crate) my_block_index: i32,
}

impl BlockBuilder {
    /// BlockBuilder.cxx L23-34.
    pub(crate) fn new() -> Self {
        BlockBuilder {
            my_oriented_shape_map: Vec::new(),
            my_oriented_shape_map_index: std::collections::HashMap::new(),
            my_oriented_shape_map_is_valid: std::collections::HashMap::new(),
            my_blocks: Vec::new(),
            my_blocks_is_regular: Vec::new(),
            my_is_done: false,
            my_block_index: 1,
        }
    }

    /// BlockBuilder.cxx L38-130.
    pub(crate) fn make_block(&mut self, brep: &BRep, ss: &mut dyn ShapeSet) {
        // Compute the set of connexity blocks of elements of element set SS.
        self.my_oriented_shape_map.clear();
        self.my_oriented_shape_map_index.clear();
        self.my_blocks.clear();
        self.my_blocks_is_regular.clear();

        let mut is_regular: bool;
        let mut cur_nei: i32;
        let mut m_extent: i32;
        let mut eindex: i32;

        ss.init_start_elements();
        while ss.more_start_elements() {
            let e = ss.start_element();
            m_extent = self.my_oriented_shape_map.len() as i32;
            eindex = self.add_element(&e);

            // E = current element of the element set SS
            // Eindex > Mextent => E is new in M
            let enewin_m = eindex > m_extent;
            if enewin_m {
                // make a new block starting at element Eindex
                self.my_blocks.push(eindex);
                is_regular = true;
                cur_nei = 0;
                // put in current block all the elements connex to E :
                // while an element E has been added to M
                //    - compute neighbours of E : N(E)
                //    - add each element N of N(E) to M
                m_extent = self.my_oriented_shape_map.len() as i32;
                let mut searchneighbours = eindex <= m_extent;
                while searchneighbours {
                    // E = element of M on which neighbours must be searched
                    let e1 = self.my_oriented_shape_map[(eindex - 1) as usize].clone();
                    cur_nei = ss.max_number_sub_shape(brep, &e1);
                    let condregu = cur_nei <= 2;
                    is_regular = is_regular && condregu;
                    // compute neighbours of E : add them to M to increase
                    // M.Extent().
                    ss.init_neighbours(brep, &e1);
                    while ss.more_neighbours() {
                        let n = ss.neighbour();
                        self.add_element(&n);
                        ss.next_neighbour(brep);
                    }

                    eindex += 1;
                    m_extent = self.my_oriented_shape_map.len() as i32;
                    searchneighbours = eindex <= m_extent;
                } // while (searchneighbours)
                let iiregu = if is_regular { 1 } else { 0 };
                self.my_blocks_is_regular.push(iiregu);
            } // if (EnewinM)
            ss.next_start_element();
        } // for ()

        // To value the l bound of the last connexity block created above,
        // we create an artificial block of value =
        // myOrientedShapeMap.Extent() + 1.  The real number of connexity
        // blocks is myBlocks.Length() - 1.
        m_extent = self.my_oriented_shape_map.len() as i32;
        self.my_blocks.push(m_extent + 1);
        self.my_is_done = true;
    }

    /// BlockBuilder.cxx L134-137.
    pub(crate) fn init_block(&mut self) {
        self.my_block_index = 1;
    }

    /// BlockBuilder.cxx L141-147.
    pub(crate) fn more_block(&self) -> bool {
        // the length of myBlocks is 1 + number of connexity blocks
        let l = self.my_blocks.len() as i32;
        let b = self.my_block_index < l;
        b
    }

    /// BlockBuilder.cxx L151-154.
    pub(crate) fn next_block(&mut self) {
        self.my_block_index += 1;
    }

    /// BlockBuilder.cxx L158-163.
    pub(crate) fn block_iterator(&self) -> BlockIterator {
        let lower = self.my_blocks[(self.my_block_index - 1) as usize];
        let upper = self.my_blocks[self.my_block_index as usize] - 1;
        BlockIterator::new2(lower, upper)
    }

    /// BlockBuilder.cxx L167-179.
    pub(crate) fn element_bi(&self, bi: &BlockIterator) -> Shape {
        let isbound = bi.more();
        if !isbound {
            panic!("OutOfRange");
        }
        let index = bi.value();
        self.my_oriented_shape_map[(index - 1) as usize].clone()
    }

    /// BlockBuilder.cxx L181-191.
    pub(crate) fn element_index(&self, index: i32) -> Shape {
        let isbound = self.my_oriented_shape_map_is_valid.contains_key(&index);
        if !isbound {
            panic!("OutOfRange");
        }
        self.my_oriented_shape_map[(index - 1) as usize].clone()
    }

    /// BlockBuilder.cxx L193-203.
    pub(crate) fn element_of(&self, e: &Shape) -> i32 {
        let isbound = self.my_oriented_shape_map_index.contains_key(&shape_key(e));
        if !isbound {
            panic!("OutOfRange");
        }
        self.my_oriented_shape_map_index[&shape_key(e)]
    }

    /// BlockBuilder.cxx L207-220.
    pub(crate) fn element_is_valid_bi(&self, bi: &BlockIterator) -> bool {
        let isbound = bi.more();
        if !isbound {
            return false;
        }
        let sindex = bi.value();
        self.element_is_valid_index(sindex)
    }

    /// BlockBuilder.cxx L222-234.
    pub(crate) fn element_is_valid_index(&self, sindex: i32) -> bool {
        let isbound = self.my_oriented_shape_map_is_valid.contains_key(&sindex);
        if !isbound {
            return false;
        }
        let isb = self.my_oriented_shape_map_is_valid[&sindex];
        let isvalid = isb == 1;
        isvalid
    }

    /// BlockBuilder.cxx L238-244.
    pub(crate) fn add_element(&mut self, s: &Shape) -> i32 {
        let sindex = match self.my_oriented_shape_map_index.get(&shape_key(s)) {
            Some(&i) => i,
            None => {
                self.my_oriented_shape_map.push(s.clone());
                let i = self.my_oriented_shape_map.len() as i32;
                self.my_oriented_shape_map_index.insert(shape_key(s), i);
                i
            }
        };
        self.my_oriented_shape_map_is_valid.insert(sindex, 1);
        sindex
    }

    /// BlockBuilder.cxx L248-260.
    pub(crate) fn set_valid_bi(&mut self, bi: &BlockIterator, isvalid: bool) {
        let isbound = bi.more();
        if !isbound {
            return;
        }
        let sindex = bi.value();
        self.set_valid_index(sindex, isvalid);
    }

    /// BlockBuilder.cxx L262-272.
    pub(crate) fn set_valid_index(&mut self, sindex: i32, isvalid: bool) {
        let isbound = self.my_oriented_shape_map_is_valid.contains_key(&sindex);
        if !isbound {
            return;
        }
        let i = if isvalid { 1 } else { 0 };
        self.my_oriented_shape_map_is_valid.insert(sindex, i);
    }

    /// BlockBuilder.cxx L276-285.
    pub(crate) fn current_block_is_regular(&self) -> bool {
        let mut b = false;
        let i = self.my_blocks_is_regular[(self.my_block_index - 1) as usize];
        if i == 1 {
            b = true;
        }
        b
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_LoopClassifier (LoopClassifier.hxx L27-41).
// =========================================================================
pub(crate) trait LoopClassifier {
    /// Returns the state of loop L1 compared with loop L2.  The `brep`
    /// argument is the D6 pool the classifier reads TShapes from (OCCT
    /// reads them through handles).
    fn compare(&mut self, brep: &mut BRep, l1: &LoopHandle, l2: &LoopHandle) -> State;
}

/// TopOpeBRepBuild_LoopEnum (AreaBuilder.hxx): ANYLOOP / BOUNDARY / BLOCK.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LoopEnum {
    AnyLoop,
    Boundary,
    Block,
}

// =========================================================================
// OCCT TopOpeBRepBuild_AreaBuilder core (AreaBuilder.cxx L30-121 core
// helpers + L334-429 iteration) + the Area2dBuilder::InitAreaBuilder
// override (Area2dBuilder.cxx L41-262) — the operative virtual for
// FaceAreaBuilder (FaceAreaBuilder.cxx L36-41 InitFaceAreaBuilder ->
// InitAreaBuilder, resolved to the 2d override).  The AreaBuilder base
// InitAreaBuilder body (AreaBuilder.cxx L125-330) is the overridden
// virtual — no call in the FaceBuilder chain resolves to it, so its body
// is not duplicated here.
// =========================================================================
#[derive(Debug, Default)]
pub(crate) struct FaceAreaBuilder {
    /// OCCT: myUNKNOWNRaise = false (ctor L30-33: "no raise if UNKNOWN
    /// state found").
    pub(crate) my_unknown_raise: bool,
    /// OCCT: myArea (list of list of loop handles).
    pub(crate) my_area: Vec<Vec<LoopHandle>>,
    pub(crate) my_area_iterator: usize,
    pub(crate) my_loop_iterator: usize,
    pub(crate) my_loop_active: bool,
}

impl FaceAreaBuilder {
    /// AreaBuilder.cxx L30-33.
    pub(crate) fn new() -> Self {
        FaceAreaBuilder {
            my_unknown_raise: false,
            my_area: Vec::new(),
            my_area_iterator: 0,
            my_loop_iterator: 0,
            my_loop_active: false,
        }
    }

    /// FaceAreaBuilder.cxx L36-41.
    pub(crate) fn init_face_area_builder(
        &mut self,
        brep: &mut BRep,
        ls: &mut LoopSet,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        self.init_area_builder_2d(brep, ls, lc, force_class);
    }

    /// AreaBuilder.cxx L59-103.
    pub(crate) fn compare_loop_with_list_of_loop(
        &self,
        lc: &mut dyn LoopClassifier,
        brep: &mut BRep,
        l: &LoopHandle,
        lol: &[LoopHandle],
        what: LoopEnum,
    ) -> State {
        let mut state = State::Unknown;
        let mut totest: bool;

        if lol.is_empty() {
            return State::Out;
        }

        for cur_l in lol {
            match what {
                LoopEnum::AnyLoop => totest = true,
                LoopEnum::Boundary => totest = cur_l.is_shape(),
                LoopEnum::Block => totest = !cur_l.is_shape(),
            }
            if totest {
                state = lc.compare(brep, l, cur_l);
                if state == State::Out {
                    // <L> is out of at least one Loop of <LOL> : stop to
                    // explore
                    break;
                }
            }
        }

        state
    }

    /// AreaBuilder.cxx L111-121.
    pub(crate) fn atomize(&self, state: &mut State, newstate: State) {
        if self.my_unknown_raise {
            assert!(state != &State::Unknown, "AreaBuilder : Position Unknown");
        } else {
            *state = newstate;
        }
    }

    /// Area2dBuilder.cxx L41-262 (the AreaBuilder::InitAreaBuilder
    /// override this batch consumes).
    pub(crate) fn init_area_builder_2d(
        &mut self,
        brep: &mut BRep,
        ls: &mut LoopSet,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        let mut state: State;
        let mut loopinside: bool;
        let mut loopoutside: bool;

        // boundaryloops : list of boundary loops out of the areas.
        let mut boundaryloops: Vec<LoopHandle> = Vec::new();

        self.my_area.clear(); // Clear the list of Area to be built

        ls.init_loop();
        while ls.more_loop() {
            // process a new loop : L is the new current Loop
            let l = ls.loop_();
            let boundary_l = l.is_shape();

            // L = shape et ForceClass  : on traite L comme un block
            // L = shape et !ForceClass : on traite L comme un pur shape
            // L = !shape               : on traite L comme un block
            let traitercommeblock = !boundary_l || force_class;
            if !traitercommeblock {
                // the loop L is a boundary loop :
                // - try to insert it in an existing area, such as L is
                //   inside all the block loops. Only block loops of the
                //   area are compared.
                // - if L could not be inserted, store it in list of
                //   boundary loops.
                loopinside = false;
                let mut area_iter = 0usize;
                while area_iter < self.my_area.len() {
                    let a_area = self.my_area[area_iter].clone();
                    if a_area.is_empty() {
                        area_iter += 1;
                        continue;
                    }
                    state =
                        self.compare_loop_with_list_of_loop(lc, brep, &l, &a_area, LoopEnum::Block);
                    if state == State::Unknown {
                        self.atomize(&mut state, State::In);
                    }
                    loopinside = state == State::In;
                    if loopinside {
                        break;
                    }
                    area_iter += 1;
                } // end of Area scan

                if loopinside {
                    // ADD_Loop_TO_LISTOFLoop(L, aArea, "IN, to current area")
                    self.my_area[area_iter].push(l.clone());
                } else if !loopinside {
                    // ADD_Loop_TO_LISTOFLoop(L, boundaryloops, "! IN, to
                    // boundaryloops")
                    boundaryloops.push(l.clone());
                }
            } // end of boundary loop
            else {
                // the loop L is a block loop
                // if L is IN theArea :
                //   - stop area scan, insert L in theArea.
                //   - remove from the area all the loops outside L
                //   - make a new area with them, unless they are all boundary
                //   - if they are all boundary put them back in boundaryLoops
                // else :
                //   - create a new area with L.
                //   - insert boundary loops that are IN the new area
                //     (and remove them from 'boundaryloops')
                loopinside = false;
                let mut area_iter = 0usize;
                while area_iter < self.my_area.len() {
                    let a_area = self.my_area[area_iter].clone();
                    if a_area.is_empty() {
                        area_iter += 1;
                        continue;
                    }
                    state = self.compare_loop_with_list_of_loop(
                        lc, brep, &l, &a_area, LoopEnum::AnyLoop,
                    );
                    if state == State::Unknown {
                        self.atomize(&mut state, State::In);
                    }
                    loopinside = state == State::In;
                    if loopinside {
                        break;
                    }
                    area_iter += 1;
                } // end of Area scan

                if loopinside {
                    let mut all_shape = true;
                    let mut removed_loops: Vec<LoopHandle> = Vec::new();
                    let mut loop_iter = 0usize;
                    while loop_iter < self.my_area[area_iter].len() {
                        let cur_l = self.my_area[area_iter][loop_iter].clone();
                        state = lc.compare(brep, &cur_l, &l);
                        if state == State::Unknown {
                            self.atomize(&mut state, State::In); // not OUT
                        }
                        loopoutside = state == State::Out;
                        if loopoutside {
                            // remove the loop from the area
                            // ADD_Loop_TO_LISTOFLoop(curL, removedLoops, ...)
                            removed_loops.push(cur_l.clone());
                            all_shape = all_shape && cur_l.is_shape();
                            // REM_Loop_FROM_LISTOFLoop(LoopIter, ...)
                            self.my_area[area_iter].remove(loop_iter);
                        } else {
                            loop_iter += 1;
                        }
                    }
                    // insert the loop in the area
                    // ADD_Loop_TO_LISTOFLoop(L, aArea, "area = current")
                    self.my_area[area_iter].push(l.clone());
                    if !removed_loops.is_empty() {
                        if all_shape {
                            // ADD_LISTOFLoop_TO_LISTOFLoop(removedLoops,
                            // boundaryloops, "allShape = 1")
                            boundaryloops.extend(removed_loops);
                        } else {
                            // make a new area with the removed loops
                            // ADD_LISTOFLoop_TO_LISTOFLoop(removedLoops,
                            // myArea.Last(), "allShape = 0")
                            self.my_area.push(removed_loops);
                        }
                    }
                } // Loopinside == True
                else {
                    let mut ashapeinside: bool;
                    let mut ablockinside: bool;
                    self.my_area.push(Vec::new());
                    let new_area0_idx = self.my_area.len() - 1;
                    // ADD_Loop_TO_LISTOFLoop(L, newArea0, "new area")
                    self.my_area[new_area0_idx].push(l.clone());

                    let mut loop_iter = 0usize;
                    while loop_iter < boundaryloops.len() {
                        ashapeinside = false;
                        ablockinside = false;
                        let lb = boundaryloops[loop_iter].clone();
                        state = lc.compare(brep, &lb, &l);
                        if state == State::Unknown {
                            self.atomize(&mut state, State::In);
                        }
                        ashapeinside = state == State::In;
                        if ashapeinside {
                            state = lc.compare(brep, &l, &lb);
                            if state == State::Unknown {
                                self.atomize(&mut state, State::In);
                            }
                            ablockinside = state == State::In;
                        }
                        if ashapeinside && ablockinside {
                            let cur_l = boundaryloops[loop_iter].clone();
                            // ADD_Loop_TO_LISTOFLoop(curL, newArea0,
                            // "ashapeinside && ablockinside, new area")
                            self.my_area[new_area0_idx].push(cur_l);
                            // REM_Loop_FROM_LISTOFLoop(LoopIter,
                            // boundaryloops, ...)
                            boundaryloops.remove(loop_iter);
                        } else {
                            loop_iter += 1;
                        }
                    } // end of boundaryloops scan
                } // Loopinside == False
            } // end of block loop
            ls.next_loop();
        } // end of LoopSet LS scan

        if !boundaryloops.is_empty() {
            if self.my_area.is_empty() {
                // Area2dBuilder.cxx L247-259: the purge area.
                let mut new_area3: Vec<LoopHandle> = Vec::new();
                new_area3.extend(boundaryloops.iter().cloned());
                self.my_area.push(new_area3);
            }
        }

        self.init_area();
    }

    /// AreaBuilder.cxx L334-340.
    pub(crate) fn init_area(&mut self) -> i32 {
        self.my_area_iterator = 0;
        self.init_loop();
        let n = self.my_area.len() as i32;
        n
    }

    /// AreaBuilder.cxx L344-348.
    pub(crate) fn more_area(&self) -> bool {
        self.my_area_iterator < self.my_area.len()
    }

    /// AreaBuilder.cxx L352-356.
    pub(crate) fn next_area(&mut self) {
        self.my_area_iterator += 1;
        self.init_loop();
    }

    /// AreaBuilder.cxx L360-374.
    pub(crate) fn init_loop(&mut self) -> i32 {
        let mut n = 0;
        if self.my_area_iterator < self.my_area.len() {
            n = self.my_area[self.my_area_iterator].len() as i32;
            self.my_loop_iterator = 0;
            self.my_loop_active = true;
        } else {
            // Create an empty iterator (AreaBuilder.cxx L370-371).
            self.my_loop_iterator = 0;
            self.my_loop_active = false;
        }
        n
    }

    /// AreaBuilder.cxx L378-382.
    pub(crate) fn more_loop(&self) -> bool {
        let b = self.my_loop_active
            && self.my_area_iterator < self.my_area.len()
            && self.my_loop_iterator < self.my_area[self.my_area_iterator].len();
        b
    }

    /// AreaBuilder.cxx L386-389.
    pub(crate) fn next_loop(&mut self) {
        self.my_loop_iterator += 1;
    }

    /// AreaBuilder.cxx L393-397.
    pub(crate) fn loop_(&self) -> LoopHandle {
        self.my_area[self.my_area_iterator][self.my_loop_iterator].clone()
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_CompositeClassifier base state
// (CompositeClassifier.cxx L44-48).  The WireEdgeClassifier overrides
// Compare (WireEdgeClassifier.cxx L62-163), so only the base carries the
// myBlockBuilder reference; the base Compare body
// (CompositeClassifier.cxx L52-132) is the overridden virtual — no call
// in the FaceBuilder chain resolves to it.
// =========================================================================
pub(crate) struct CompositeClassifierBase<'a> {
    pub(crate) my_block_builder: &'a BlockBuilder,
}

// =========================================================================
// OCCT BRepClass_FacePassiveClassifier carrier.
//
// OCCT BRepClass_FacePassiveClassifier.cxx L27-73 is byte-identical to
// BRepClass_FClass2dOfFClassifier.cxx (both delegate Reset / Compare to
// the same TopClass_Classifier2d state machine over identical members —
// verified by a modulo-class-name diff of the two .cxx files).  The rcad
// FClass2dOfFClassifier translation
// (topalgo/brep_class/face_classifier.rs) is therefore the 1:1 carrier of
// the passive classifier too; the WireEdgeClassifier drives it through
// my_fpc.reset / my_fpc.compare / my_fpc.state below.
// =========================================================================

// =========================================================================
// D6 flat-index view over the WireEdgeSet content for the index-based
// 2D classifier machinery (the OCCT WireEdgeClassifier works on TopoDS
// handles directly; rcad's translated FClass2dOfFClassifier /
// BRepClass_Intersector are flat-index based over a ShapeSource).
// Index 0 = the classification face, the compared edges are registered
// from 1 on (the FaceShapeSource convention, shape_source.rs L93+).
// =========================================================================
pub(crate) struct WecShapeSource {
    face: Shape,
    surf: Option<rcad_kernel::geom::Surface3>,
    edges: Vec<Shape>,
    edge_index: std::collections::HashMap<(u64, u32), usize>,
    locations: Vec<DAffine3>,
}

impl WecShapeSource {
    pub(crate) fn new(face: &Shape, brep: &BRep) -> Self {
        let surf = face_surface_value(brep, face);
        WecShapeSource {
            face: face.clone(),
            surf,
            edges: Vec::new(),
            edge_index: std::collections::HashMap::new(),
            locations: ds_locations(brep),
        }
    }

    /// Register (or find) the flat index of an edge to be classified.
    pub(crate) fn register_edge(&mut self, e: &Shape) -> usize {
        if let Some(&i) = self.edge_index.get(&shape_key(e)) {
            return i;
        }
        self.edges.push(e.clone());
        let i = self.edges.len();
        self.edge_index.insert(shape_key(e), i);
        i
    }
}

impl ShapeSource for WecShapeSource {
    fn nb_shapes(&self) -> usize {
        1 + self.edges.len()
    }
    fn shape_at(&self, i: usize) -> Shape {
        if i == 0 {
            self.face.clone()
        } else {
            self.edges.get(i - 1).cloned().unwrap_or_else(Shape::null)
        }
    }
    fn shape_type(&self, i: usize) -> ShapeType {
        self.shape_at(i).shape_type()
    }
    fn sub_shapes(&self, _i: usize) -> &[usize] {
        &[]
    }
    fn map_shape_index(&self, ptr_id: u64, location: u32) -> Option<usize> {
        if self.face.ptr_id() == ptr_id && self.face.location == location {
            Some(0)
        } else {
            self.edge_index.get(&(ptr_id, location)).copied()
        }
    }
    fn map_ve(&self, _vertex: usize) -> Option<&Vec<usize>> {
        None
    }
    fn face_surface(&self, i: usize) -> Option<rcad_kernel::geom::Surface3> {
        if i == 0 {
            self.surf.clone()
        } else {
            None
        }
    }
    fn vertex_tolerance(&self, i: usize) -> f64 {
        vertex_tolerance_of(&self.shape_at(i))
    }
    fn is_edge_degenerated(&self, i: usize) -> bool {
        match &*self.shape_at(i).data {
            TShape::Edge(ed) => ed.degenerated,
            _ => false,
        }
    }
    fn get_location(&self, idx: u32) -> DAffine3 {
        self.locations
            .get(idx as usize)
            .copied()
            .unwrap_or(DAffine3::IDENTITY)
    }
    fn locations(&self) -> &[DAffine3] {
        &self.locations
    }
}

// =========================================================================
// OCCT BRep_Tool / TopOpeBRepTool GAP carriers shared by the classifiers.
// =========================================================================

/// OCCT BRep_Tool::Parameter(V, E) (BRep_Tool.cxx L1631-1706, the
/// edge-representation read).  The rcad edge carries the vertex
/// parameters in TEdgeData::vertex_params keyed by vertex TShape.
pub(crate) fn brep_tool_parameter(brep: &BRep, v: &Shape, e: &Shape) -> f64 {
    if let TShape::Edge(ed) = &*e.data {
        if let Some(&p) = ed.vertex_params.get(&v.ptr_id()) {
            return p;
        }
    }
    let _ = brep;
    0.0
}

/// OCCT FC2D_HasCurveOnSurface (TopOpeBRepTool_2d.cxx L147-154) — the
/// real body is `crate::fillet::topopebrep_tool_2d::fc2d_has_curve_on_surface`;
/// this adapter keeps the historical (E, F) -> bool call-site form.
pub(crate) fn fc2d_has_curve_on_surface(brep: &BRep, e: &Shape, f: &Shape) -> bool {
    crate::fillet::topopebrep_tool_2d::fc2d_has_curve_on_surface(brep, e, f)
}

/// OCCT FC2D_CurveOnSurface(E, F, f, l, tol, trim3d) (TopOpeBRepTool_2d.cxx
/// L356-376) — the real body is
/// `crate::fillet::topopebrep_tool_2d::fc2d_curve_on_surface`; this
/// adapter repacks the OCCT out-parameters (f, l, tol) into the
/// historical tuple form the call sites consume.  `trim3d` follows the
/// OCCT call site (the no-EF overload defaults it to false in OCCT).
pub(crate) fn fc2d_curve_on_surface(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    trim3d: bool,
) -> Option<(Curve2d, f64, f64, f64)> {
    // OCCT call-site form: double f2, l2, tolpc; (uninitialized locals —
    // Rust requires initialization).
    let mut f2: f64 = 0.0;
    let mut l2: f64 = 0.0;
    let mut tolpc: f64 = 0.0;
    let c2d = crate::fillet::topopebrep_tool_2d::fc2d_curve_on_surface(
        brep, e, f, &mut f2, &mut l2, &mut tolpc, trim3d,
    );
    c2d.map(|c2d| (c2d, f2, l2, tolpc))
}

/// OCCT BB.UpdateEdge(E, C2D, F, tol) (BRep_Builder.lxx L92-98 ->
/// BRep_Builder.cxx L655-672 `UpdateEdge(E, C, S, L, Tol)` ->
/// BRep_Builder.cxx UpdateCurves L104-167) — the pcurve is stored under
/// the face-key the kernel pcurve map uses; the edge tolerance is raised
/// to `tol` (the UpdateEdge max-tolerance rule).
///
/// The stored range follows UpdateCurves' two-step rule: the representation
/// is seeded from the 2D curve's own range (L151-153, the
/// `new BRep_CurveOnSurface(C, S, L)` constructor -> `COS->Range(aFCur,
/// aLCur)`) and is then OVERWRITTEN by the range of the edge's
/// Curve3D representation whenever that range is finite (L116-129
/// `GC->Range(f, l)` on the IsCurve3D entry with the L112
/// `-/+Precision::Infinite()` seed, then L154-162
/// `if (!Precision::IsInfinite(f)) aFCur = f;`).  The OCCT overload carries
/// no f/l arguments — the `f2`/`l2` out-parameters of FC2D_CurveOnSurface
/// feed only the caller's own probes (WireEdgeClassifier.cxx ResetElement
/// L449-451 / CompareElement L490-492), never the UpdateEdge call.
pub(crate) fn bb_update_edge_pcurve(
    brep: &mut BRep,
    e: &Shape,
    f: &Shape,
    c2d: &Curve2d,
    tol: f64,
) {
    let key_loc = brep.compose_pcurve_location(f.location, e.location);
    let key = (f.ptr_id(), key_loc);
    let ed = brep.edge_mut_inplace(e.clone());
    let [mut a_f, mut a_l] = c2d.default_domain();
    if ed.curve.is_some() {
        if !is_infinite_value(ed.range[0]) {
            a_f = ed.range[0];
        }
        if !is_infinite_value(ed.range[1]) {
            a_l = ed.range[1];
        }
    }
    ed.pcurves.insert(key, (c2d.clone(), a_f, a_l));
    if ed.tolerance < tol {
        ed.tolerance = tol;
    }
}

/// OCCT BRep_Builder::Add(F, W) for wires (BRep_Builder.cxx L929-1158):
/// the first wire added to a face is its outer wire, the following ones
/// are inner wires.
pub(crate) fn bb_add_face_wire(brep: &mut BRep, f: &Shape, w: &Shape) {
    let fd = brep.face_mut(f.clone());
    if fd.outer_wire.is_null() {
        fd.outer_wire = w.clone();
    } else {
        fd.inner_wires.push(w.clone());
    }
    fd.my_shapes.push(w.clone());
}

/// WireEdgeClassifier.cxx L212-228 — FUN_tgINE: tg oriented INSIDE 1d(e);
/// vl : last vertex of e.
///
/// GAP carrier of TopOpeBRepTool_TOOL::TggeomE (TopOpeBRepTool_TOOL.cxx,
/// TKBool/TopOpeBRepTool — external untranslated dependency): the tangent
/// of the edge's 3D curve at parameter `par` (TggeomE reads the located
/// curve D1).  The OCCT !ok path returns the null (0,0,0) vector.
pub(crate) fn fun_tg_ine(brep: &BRep, v: &Shape, vl: &Shape, e: &Shape) -> DVec3 {
    let par = brep_tool_parameter(brep, v, e);
    let tg = match brep.edge_curve_world(e) {
        Some((curve, _range)) => {
            // TopOpeBRepTool_TOOL::TggeomE: the curve D1 at par.
            curve.derivative_at(par)
        }
        None => DVec3::ZERO, // NYIRAISE
    };
    let mut tg = tg;
    if shapes_same(v, vl) {
        tg = -tg;
    }
    tg
}

// =========================================================================
// OCCT TopOpeBRepBuild_WireEdgeClassifier (WireEdgeClassifier.cxx
// L52-527).
// =========================================================================

// WireEdgeClassifier.cxx L44-48 (the classification-result codes of the
// disabled FUN_tool_classiBnd2d tail — kept as the OCCT literal).
pub(crate) const SAME: i32 = -1;
#[allow(dead_code)]
pub(crate) const DIFF: i32 = -2;
#[allow(dead_code)]
pub(crate) const UNKNOWN: i32 = 0;
#[allow(dead_code)]
pub(crate) const ONE_INTWO: i32 = 1;
#[allow(dead_code)]
pub(crate) const TWO_INONE: i32 = 2;

pub(crate) struct WireEdgeClassifier<'a> {
    /// CompositeClassifier.cxx L44-48: myBlockBuilder.
    pub(crate) base: CompositeClassifierBase<'a>,
    my_first_compare: bool,
    my_point2d: DVec2,
    /// BRepClass_Edge myBCEdge — the (edge, face) pair
    /// (BRepClass_Edge.hxx L40-52).
    my_bcedge_face: Shape,
    my_bcedge_edge: Shape,
    /// BRepClass_FacePassiveClassifier myFPC — the FClass2dOfFClassifier
    /// carrier (see the carrier note above).
    my_fpc: FClass2dOfFClassifier,
    my_shape: Shape,
    /// D6: the flat-index view the index-based classifier reads.
    src: WecShapeSource,
}

impl<'a> WireEdgeClassifier<'a> {
    /// WireEdgeClassifier.cxx L52-58.
    pub(crate) fn new(brep: &BRep, f: &Shape, bb: &'a BlockBuilder, wes_shapes: &[Shape]) -> Self {
        let mut src = WecShapeSource::new(f, brep);
        // Register the edges the classifier will walk: the BlockBuilder
        // elements and the edges of the ShapeSet shapes (old wires).  D6:
        // the index space of the translated classifier machinery.
        for i in 1..=bb.my_oriented_shape_map.len() {
            let e = bb.my_oriented_shape_map[i - 1].clone();
            src.register_edge(&e);
        }
        for s in wes_shapes {
            for e in shape_explore_children(brep, s, ShapeType::Edge) {
                src.register_edge(&e);
            }
        }
        WireEdgeClassifier {
            base: CompositeClassifierBase { my_block_builder: bb },
            my_first_compare: true,
            my_point2d: DVec2::ZERO,
            my_bcedge_face: f.clone(),
            my_bcedge_edge: Shape::null(),
            my_fpc: FClass2dOfFClassifier::new(),
            my_shape: Shape::null(),
            src,
        }
    }

    /// WireEdgeClassifier.cxx L62-163.
    pub(crate) fn compare(&mut self, brep: &mut BRep, l1: &LoopHandle, l2: &LoopHandle) -> State {
        let mut state = State::Unknown;

        let isshape1 = l1.is_shape();
        let isshape2 = l2.is_shape();

        if isshape2 && isshape1 {
            // L1 is Shape , L2 is Shape
            let s1 = l1.shape().clone();
            let s2 = l2.shape().clone();
            state = self.compare_shapes(brep, &s1, &s2);
        } else if isshape2 && !isshape1 {
            // L1 is Block , L2 is Shape
            let mut bit1 = *l1.block_iterator();
            bit1.initialize();
            let mut yena1 = bit1.more();
            while yena1 {
                let s1 = self.base.my_block_builder.element_bi(&bit1);
                let s2 = l2.shape().clone();
                state = self.compare_element_to_shape(brep, &s1, &s2);
                yena1 = false;
                if state == State::Unknown {
                    if bit1.more() {
                        bit1.next();
                    }
                    yena1 = bit1.more();
                }
            }
        } else if !isshape2 && isshape1 {
            // L1 is Shape , L2 is Block
            let s1 = l1.shape().clone();
            self.reset_shape(brep, &s1);
            let mut bit2 = *l2.block_iterator();
            bit2.initialize();
            while bit2.more() {
                let s2 = self.base.my_block_builder.element_bi(&bit2);
                self.compare_element(brep, &s2);
                bit2.next();
            }
            state = self.state();
        } else if !isshape2 && !isshape1 {
            // L1 is Block , L2 is Block
            if state == State::Unknown {
                let mut bit1 = *l1.block_iterator();
                bit1.initialize();
                let mut yena1 = bit1.more();
                while yena1 {
                    let s1 = self.base.my_block_builder.element_bi(&bit1);
                    self.reset_element(brep, &s1);
                    let mut bit2 = *l2.block_iterator();
                    bit2.initialize();
                    while bit2.more() {
                        let s2 = self.base.my_block_builder.element_bi(&bit2);
                        self.compare_element(brep, &s2);
                        bit2.next();
                    }
                    state = self.state();
                    yena1 = false;
                    if state == State::Unknown {
                        if bit1.more() {
                            bit1.next();
                        }
                        yena1 = bit1.more();
                    }
                }
            } // UNKNOWN

            if state == State::Unknown {
                let s1 = self.loop_to_shape(brep, l1);
                if s1.is_null() {
                    return state;
                }
                let s2 = self.loop_to_shape(brep, l2);
                if s2.is_null() {
                    return state;
                }
                // WireEdgeClassifier.cxx L153-159: the
                // TopOpeBRepTool_ShapeClassifier fallback
                // (FSC_GetPSC().SetReference(s2) +
                // StateShapeReference(s1, s2) with SameDomain forced to 1)
                // — TKBool/TopOpeBRepTool external machinery, GAP carrier:
                // the state stays UNKNOWN, the OCCT State() value when the
                // fallback cannot decide is preserved.
                let _ = (&s1, &s2);
                state = State::Unknown;
            } // UNKNOWN
        }
        state
    }

    /// WireEdgeClassifier.cxx L167-210.
    pub(crate) fn loop_to_shape(&mut self, brep: &mut BRep, l: &LoopHandle) -> Shape {
        self.my_shape = Shape::null();
        let mut bit = *l.block_iterator();
        bit.initialize();
        if !bit.more() {
            return self.my_shape.clone();
        }

        let f1_shape = self.my_bcedge_face.clone();
        // F1.EmptyCopied() (WireEdgeClassifier.cxx L178-182).
        let f = brep.empty_copied(&f1_shape);
        let mut bb = rcad_kernel::topods::BRepBuilder::new();
        let w = bb.make_wire(brep);
        while bit.more() {
            let e = self.base.my_block_builder.element_bi(&bit);
            let tol_e = brep.tolerance(&e);
            let haspc = fc2d_has_curve_on_surface(brep, &e, &f1_shape);
            if !haspc {
                // WireEdgeClassifier.cxx L195-203: C2D =
                // FC2D_CurveOnSurface(E, F, f, l, tolpc); if (!C2D.IsNull())
                // BB.UpdateEdge(E, C2D, F, max(tolpc, tolE)).
                if let Some((c2d, _f2, _l2, tolpc)) = fc2d_curve_on_surface(brep, &e, &f1_shape, false)
                {
                    let tol = tolpc.max(tol_e);
                    bb_update_edge_pcurve(brep, &e, &f, &c2d, tol);
                }
            }
            // BB.Add(W, E).
            bb.add_to_wire(brep, w.clone(), e);
            bit.next();
        }
        // BB.Add(F, W).
        bb_add_face_wire(brep, &f, &w);

        self.my_shape = f.clone();
        self.my_shape.clone()
    }

    /// WireEdgeClassifier.cxx L232-385.
    pub(crate) fn compare_shapes(&mut self, brep: &mut BRep, b1: &Shape, b2: &Shape) -> State {
        // evolution xpu10198 : on classifie 1 wire / 1 des wires connexes.
        let mut state = State::Unknown;
        let ex1 = shape_explore_children(brep, b1, ShapeType::Edge);
        if ex1.is_empty() {
            return state;
        }
        for e1 in &ex1 {
            let (vf1, vl1) = top_exp_vertices(brep, e1); // xpu10198
            let e1clo = shapes_same(&vf1, &vl1);
            let mut mapv1: Vec<(u64, u32)> = Vec::new();
            mapv1.push(shape_key(&vf1));
            mapv1.push(shape_key(&vl1));

            self.reset_shape(brep, e1);
            let mut indy = false;
            let ex = shape_explore_children(brep, b2, ShapeType::Edge);
            for e in &ex {
                if shapes_same(e, e1) {
                    state = State::Unknown;
                    break;
                } // eap occ416
                let (vf, vl) = top_exp_vertices(brep, e); // xpu10198
                let eclo = shapes_same(&vf, &vl); // xpu10198
                let hasf = mapv1.contains(&shape_key(&vf)); // xpu10198
                let hasl = mapv1.contains(&shape_key(&vl)); // xpu10198
                let filter = (hasf || hasl) && (!e1clo) && (!eclo); // nyi : Eclo || e1clo
                if filter {
                    // xpu10198
                    let mut vshared = Shape::null();
                    if hasf {
                        vshared = vf.clone();
                    }
                    if hasl {
                        vshared = vl.clone();
                    }
                    let tg1 = fun_tg_ine(brep, &vshared, &vl1, e1);
                    let tg = fun_tg_ine(brep, &vshared, &vl, e);
                    let dot = tg1.dot(tg);
                    let tol = rcad_kernel::core::precision::ANGULAR * 1.0e4; // nyixpu
                    let undecided = (1.0 + dot).abs() < tol;
                    if undecided {
                        indy = true;
                    }
                } // xpu10198
                if indy {
                    state = State::Unknown;
                    break;
                }
                self.compare_element(brep, e);
                state = self.state();
            } // ex(B2,EDGE)
            if state != State::Unknown {
                break;
            }
        } // ex1

        let mut resta = state == State::Unknown;
        resta =
            resta && (b2.shape_type() == ShapeType::Wire) && (b1.shape_type() == ShapeType::Wire);
        if resta {
            let mut mape1: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
            for e in shape_explore_children(brep, b1, ShapeType::Edge) {
                mape1.insert(shape_key(&e));
            }
            // recall : avoid auto-intersection wires (ie B1 and B2 are
            // disjoint)
            let ex2 = shape_explore_children(brep, b2, ShapeType::Edge);
            for e2 in &ex2 {
                if mape1.contains(&shape_key(e2)) {
                    continue;
                }

                let the_face = self.my_bcedge_face.clone();

                // p2d on E2 of B2, E2 not shared by B1
                // (WireEdgeClassifier.cxx L337-343: ftmp = theFace
                // Oriented(FORWARD), F2 = ftmp.EmptyCopied(),
                // BB.Add(F2, TopoDS::Wire(B2))).
                let ftmp = shape_oriented(&the_face, Orientation::Forward);
                let f2 = brep.empty_copied(&ftmp);
                bb_add_face_wire(brep, &f2, b2);

                // WireEdgeClassifier.cxx L345-350: BRepAdaptor_Curve2d
                // BC2d(E2, F2); FUN_tool_bounds(E2, f, l);
                // p2 = (1 - 0.45678) * l + 0.45678 * f; p2d = BC2d.Value(p2).
                //
                // GAP carrier of BRepAdaptor_Curve2d::Value (the adaptor
                // evaluates BRep_Tool::CurveOnSurface(E2, F2)) and of
                // FUN_tool_bounds (TopOpeBRepTool_TOOL.cxx — the edge's
                // pcurve parameter range): the stored pcurve + its range
                // are read directly.
                let Some((c2d_e2, f_bnd, l_bnd)) = brep.curve_on_surface(e2, &f2) else {
                    continue;
                };
                let x = 0.45678;
                let p2 = (1.0 - x) * l_bnd + x * f_bnd;
                let p2d = c2d_e2.point_at(p2);

                // WireEdgeClassifier.cxx L352-360: F1 =
                // ftmp.EmptyCopied() + BB.Add(F1, TopoDS::Wire(B1));
                // tolF1 = BRep_Tool::Tolerance(F1);
                // BRepClass_FaceClassifier Fclass(F1, p2d, tolF1);
                // state = Fclass.State().
                let f1 = brep.empty_copied(&ftmp);
                bb_add_face_wire(brep, &f1, b1);
                let tol_f1 = brep.tolerance(&f1);
                let locations = ds_locations(brep);
                match face_surface_value(brep, &f1) {
                    Some(surf) => {
                        let fss = crate::topalgo::shape_source::FaceShapeSource::new(
                            &f1, surf, &locations,
                        );
                        let mut fclass = FClassifier::new();
                        fclass.perform(&fss, 0, p2d, tol_f1);
                        state = fclass.state();
                    }
                    None => {
                        state = State::Unknown;
                    }
                }
                return state;
            } // ex2
        }

        state
    }

    /// WireEdgeClassifier.cxx L389-402.
    pub(crate) fn compare_element_to_shape(
        &mut self,
        brep: &mut BRep,
        ee: &Shape,
        b: &Shape,
    ) -> State {
        // isEdge : edge E inits myPoint2d
        self.reset_element(brep, ee);
        for e in shape_explore_children(brep, b, ShapeType::Edge) {
            self.compare_element(brep, &e);
        }
        let state = self.state();
        state
    }

    /// WireEdgeClassifier.cxx L406-421.
    pub(crate) fn reset_shape(&mut self, brep: &mut BRep, b: &Shape) {
        if b.shape_type() == ShapeType::Edge {
            self.reset_element(brep, b);
        } else {
            let ex = shape_explore_children(brep, b, ShapeType::Edge);
            if let Some(e) = ex.first() {
                self.reset_element(brep, e);
            }
        }
    }

    /// WireEdgeClassifier.cxx L425-462.
    pub(crate) fn reset_element(&mut self, brep: &mut BRep, ee: &Shape) {
        let e = ee.clone();
        let f = self.my_bcedge_face.clone();
        let haspc = fc2d_has_curve_on_surface(brep, &e, &f); // jyl980406+
        if !haspc {
            // jyl980406+
            // bool trim3d = true; C2D = FC2D_CurveOnSurface(E,F,f2,l2,tolpc,trim3d);
            if let Some((c2d, _f2, _l2, tolpc)) = fc2d_curve_on_surface(brep, &e, &f, true) {
                let tol_e = brep.tolerance(&e); // jyl980406+
                let tol = tol_e.max(tolpc); // jyl980406+
                bb_update_edge_pcurve(brep, &e, &f, &c2d, tol); // jyl980406+
            }
        }

        // C2D = FC2D_CurveOnSurface(E, F, f2, l2, tolpc);
        match fc2d_curve_on_surface(brep, &e, &f, false) {
            Some((c2d, f2, l2, _tolpc)) => {
                let t = 0.397891143689;
                let par = (1.0 - t) * f2 + t * l2;
                self.my_point2d = c2d.point_at(par);
            }
            None => {
                // WireEdgeClassifier.cxx L444-447.
                panic!("WEC : ResetElement");
            }
        }

        self.my_first_compare = true;
    }

    /// WireEdgeClassifier.cxx L466-519.
    pub(crate) fn compare_element(&mut self, brep: &mut BRep, ee: &Shape) {
        let b_ret = true;
        let e = ee.clone();
        let f = self.my_bcedge_face.clone();

        let haspc = fc2d_has_curve_on_surface(brep, &e, &f); // jyl980402+
        if !haspc {
            // jyl980402+
            if let Some((c2d, _f2, _l2, tolpc)) = fc2d_curve_on_surface(brep, &e, &f, true) {
                let tol_e = brep.tolerance(&e); // jyl980402+
                let tol = tol_e.max(tolpc); // jyl980402+
                bb_update_edge_pcurve(brep, &e, &f, &c2d, tol); // jyl980402+
            }
        }

        if self.my_first_compare {
            let Some((c2d, f2, l2, _tolpc)) = fc2d_curve_on_surface(brep, &e, &f, false) else {
                return;
            };
            let t = 0.33334567;
            let par = (1.0 - t) * f2 + t * l2;
            let p2d = c2d.point_at(par);

            // NYI : p2d peut etre un point ou la courbe n'est pas C1.
            // WireEdgeClassifier.cxx L504-509: gp_Vec2d v2d(myPoint2d,
            // p2d); gp_Lin2d l2d(myPoint2d, v2d); dist =
            // myPoint2d.Distance(p2d); myFPC.Reset(l2d, dist,
            // Precision::PConfusion()).
            let v2d = p2d - self.my_point2d;
            let l2d = Line2d::new(self.my_point2d, v2d);
            let dist = self.my_point2d.distance(p2d);
            let tol2d = rcad_kernel::core::precision::p_confusion(); // NYI : a voir
            self.my_fpc.reset(&l2d, dist, tol2d);
            self.my_first_compare = false;
        }

        // WireEdgeClassifier.cxx L512-514: myBCEdge.Edge() = E;
        // myFPC.Compare(myBCEdge, Eori).
        self.my_bcedge_edge = e.clone();
        let e_ori = e.orientation;
        let e_idx = self.src.register_edge(&e);
        let class_edge = ClassEdge::new(e_idx, 0);
        self.my_fpc.compare(&class_edge, e_ori, &self.src);
        let _ = b_ret;
    }

    /// WireEdgeClassifier.cxx L523-527.
    pub(crate) fn state(&self) -> State {
        let state = self.my_fpc.state();
        state
    }
}

/// The LoopClassifier virtual dispatch (LoopClassifier.hxx L27-41):
/// FaceAreaBuilder::InitFaceAreaBuilder receives the classifier through
/// the base-class reference; the WireEdgeClassifier::Compare override is
/// the operative body.
impl LoopClassifier for WireEdgeClassifier<'_> {
    fn compare(&mut self, brep: &mut BRep, l1: &LoopHandle, l2: &LoopHandle) -> State {
        // Resolves to the inherent WireEdgeClassifier::compare
        // (WireEdgeClassifier.cxx L62-163).
        WireEdgeClassifier::compare(self, brep, l1, l2)
    }
}
