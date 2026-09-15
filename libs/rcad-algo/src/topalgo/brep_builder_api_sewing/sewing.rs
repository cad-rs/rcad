//! OCCT BRepBuilderAPI_Sewing class body — the struct, the map carriers,
//! the constructor / Init / Load / Add / Perform orchestration, the output
//! accessors (cxx L2353-2595), the .lxx inline setters and the
//! context accessors (cxx L5943-5946).
//!
//! OCCT anchors are given per item: `BRepBuilderAPI_Sewing.hxx L...` /
//! `BRepBuilderAPI_Sewing.cxx L...` / `BRepBuilderAPI_Sewing.lxx L...`.
//!
//! The owning `BRep` pool is the leading `brep` argument of every method
//! that reads/writes TShape data (the brep_tools_quilt.rs convention); the
//! class itself carries no pool (architecture difference #3 of mod.rs).

use std::collections::HashMap;

use indexmap::IndexMap;
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;

// The OCCT method groups, split across child files (each child carries the
// corresponding `impl BRepBuilderAPISewing` block + the OCCT static
// file-scope helpers of its group).
#[path = "sewing_selectors.rs"]
pub(crate) mod selectors;
#[path = "sewing_closed.rs"]
mod closed;
#[path = "sewing_same_param.rs"]
mod same_param;
#[path = "sewing_eval.rs"]
mod eval;
#[path = "sewing_analysis.rs"]
mod analysis;
#[path = "sewing_merging.rs"]
mod merging;
#[path = "sewing_cutting.rs"]
mod cutting;
#[path = "sewing_output.rs"]
mod output;

// ---------------------------------------------------------------------------
// NCollection carriers (architecture difference #1).
// ---------------------------------------------------------------------------

/// OCCT `NCollection_IndexedDataMap<TopoDS_Shape, TopoDS_Shape,
/// TopTools_ShapeMapHasher>`.
pub(crate) type IdxShapeMap = IndexMap<bat::ShapeKey, (Shape, Shape)>;
/// OCCT `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>`.
pub(crate) type IdxListMap = IndexMap<bat::ShapeKey, (Shape, Vec<Shape>)>;
/// OCCT `NCollection_IndexedMap` / `NCollection_Map` of shapes
/// (insertion-ordered set carrier).
pub(crate) type ShapeSet = IndexMap<bat::ShapeKey, Shape>;
/// OCCT `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, ...>`.
pub(crate) type DataShapeMap = HashMap<bat::ShapeKey, Shape>;
/// OCCT `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, ...>`.
pub(crate) type DataListMap = HashMap<bat::ShapeKey, Vec<Shape>>;

/// OCCT `NCollection_IndexedDataMap::Add(key, value)` — binds the key;
/// an already bound key keeps its existing value
/// (NCollection_IndexedDataMap.hxx L498-501).
pub(crate) fn idx_map_add(m: &mut IdxShapeMap, k: &Shape, v: Shape) {
    m.entry(bat::shape_key(k)).or_insert((k.clone(), v));
}

/// OCCT `NCollection_IndexedDataMap<Shape, List>::Add(key, value)`.
pub(crate) fn idx_list_add(m: &mut IdxListMap, k: &Shape, v: Vec<Shape>) {
    m.entry(bat::shape_key(k)).or_insert((k.clone(), v));
}

/// OCCT `NCollection_IndexedMap::Add` — inserts when absent, returns the
/// insertion rank (1-based); a repeated key keeps its rank.
pub(crate) fn set_add(set: &mut ShapeSet, sh: &Shape) -> usize {
    match set.entry(bat::shape_key(sh)) {
        indexmap::map::Entry::Occupied(e) => e.index() + 1,
        indexmap::map::Entry::Vacant(v) => {
            v.insert(sh.clone());
            set.len()
        }
    }
}

/// OCCT `NCollection_Map::Add` — the boolean form (true when newly added).
pub(crate) fn set_add_bool(set: &mut ShapeSet, sh: &Shape) -> bool {
    match set.entry(bat::shape_key(sh)) {
        indexmap::map::Entry::Occupied(_) => false,
        indexmap::map::Entry::Vacant(v) => {
            v.insert(sh.clone());
            true
        }
    }
}

/// OCCT `NCollection_Map::Contains`.
pub(crate) fn set_contains(set: &ShapeSet, sh: &Shape) -> bool {
    set.contains_key(&bat::shape_key(sh))
}

/// OCCT `BRep_Tool::Degenerated(E)` (BRep_Tool.cxx L1073-1086) — the edge
/// TShape's degenerated flag.
pub(crate) fn brep_tool_degenerated(edg: &Shape) -> bool {
    matches!(&*edg.data, rcad_kernel::topo::topods::TShape::Edge(ed) if ed.degenerated)
}

/// OCCT `TopoDS_Shape::Nullify()`.
pub(crate) fn nullify(s: &mut Shape) {
    *s = Shape::null();
}

// ---------------------------------------------------------------------------
// The class (hxx L68-420).
// ---------------------------------------------------------------------------

/// OCCT BRepBuilderAPI_Sewing (hxx L68) — identifies contiguous boundaries
/// and assembles contiguous shapes into one shape.
pub struct BRepBuilderAPISewing {
    // hxx L381-410 (protected members).
    pub(crate) my_tolerance: f64,      // hxx L381: myTolerance
    pub(crate) my_sewing: bool,        // hxx L382: mySewing
    pub(crate) my_analysis: bool,      // hxx L383: myAnalysis
    pub(crate) my_cutting: bool,       // hxx L384: myCutting
    pub(crate) my_nonmanifold: bool,   // hxx L385: myNonmanifold
    pub(crate) my_old_shapes: IdxShapeMap, // hxx L386: myOldShapes
    pub(crate) my_sewed_shape: Shape,  // hxx L387: mySewedShape
    pub(crate) my_degenerated: ShapeSet, // hxx L388: myDegenerated
    pub(crate) my_free_edges: ShapeSet, // hxx L389: myFreeEdges
    pub(crate) my_multiple_edges: ShapeSet, // hxx L390: myMultipleEdges
    pub(crate) my_contigous_edges: IdxListMap, // hxx L391-392: myContigousEdges
    pub(crate) my_contig_sec_bound: DataShapeMap, // hxx L393: myContigSecBound
    pub(crate) my_nb_shapes: i32,      // hxx L394: myNbShapes
    pub(crate) my_nb_vertices: i32,    // hxx L395: myNbVertices
    pub(crate) my_nb_edges: i32,       // hxx L396: myNbEdges
    pub(crate) my_bound_faces: IdxListMap, // hxx L397-398: myBoundFaces
    pub(crate) my_bound_sections: DataListMap, // hxx L399-400: myBoundSections
    pub(crate) my_section_bound: DataShapeMap, // hxx L401: mySectionBound
    pub(crate) my_vertex_node: IdxShapeMap, // hxx L402: myVertexNode
    pub(crate) my_vertex_node_free: IdxShapeMap, // hxx L403: myVertexNodeFree
    pub(crate) my_node_sections: DataListMap, // hxx L404-405: myNodeSections
    pub(crate) my_cutting_node: DataListMap, // hxx L406-407: myCuttingNode
    pub(crate) my_little_face: ShapeSet, // hxx L408: myLittleFace
    pub(crate) my_shape: Shape,        // hxx L409: myShape
    pub(crate) my_re_shape: crate::shhealing::shape_build::reshape::ShapeBuildReShape, // hxx L410: myReShape
    // hxx L413-419 (private members).
    pub(crate) my_face_mode: bool,             // hxx L413: myFaceMode
    pub(crate) my_floating_edges_mode: bool,   // hxx L414: myFloatingEdgesMode
    pub(crate) my_same_parameter_mode: bool,   // hxx L415: mySameParameterMode
    pub(crate) my_local_tolerance_mode: bool,  // hxx L416: myLocalToleranceMode
    pub(crate) my_min_tolerance: f64,          // hxx L417: myMinTolerance
    pub(crate) my_max_tolerance: f64,          // hxx L418: myMaxTolerance
    pub(crate) my_merged_edges: ShapeSet,      // hxx L419: myMergedEdges
}

impl Default for BRepBuilderAPISewing {
    fn default() -> Self {
        BRepBuilderAPISewing::new(1.0e-06)
    }
}

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::BRepBuilderAPI_Sewing(tolerance = 1.0e-06,
    /// option1 = true, option2 = true, option3 = true, option4 = false)
    /// (hxx L78-82, cxx L2090-2103).
    pub fn new(tolerance: f64) -> Self {
        BRepBuilderAPISewing::new_full(tolerance, true, true, true, false)
    }

    /// OCCT BRepBuilderAPI_Sewing::BRepBuilderAPI_Sewing — the explicit
    /// option form (Rust has no default arguments).
    pub fn new_full(
        tolerance: f64,
        option1: bool,
        option2: bool,
        option3: bool,
        option4: bool,
    ) -> Self {
        // OCCT L2092: myReShape = new BRepTools_ReShape;
        let mut this = BRepBuilderAPISewing {
            my_tolerance: 0.0,
            my_sewing: false,
            my_analysis: false,
            my_cutting: false,
            my_nonmanifold: false,
            my_old_shapes: IdxShapeMap::new(),
            my_sewed_shape: Shape::null(),
            my_degenerated: ShapeSet::new(),
            my_free_edges: ShapeSet::new(),
            my_multiple_edges: ShapeSet::new(),
            my_contigous_edges: IdxListMap::new(),
            my_contig_sec_bound: DataShapeMap::new(),
            my_nb_shapes: 0,
            my_nb_vertices: 0,
            my_nb_edges: 0,
            my_bound_faces: IdxListMap::new(),
            my_bound_sections: DataListMap::new(),
            my_section_bound: DataShapeMap::new(),
            my_vertex_node: IdxShapeMap::new(),
            my_vertex_node_free: IdxShapeMap::new(),
            my_node_sections: DataListMap::new(),
            my_cutting_node: DataListMap::new(),
            my_little_face: ShapeSet::new(),
            my_shape: Shape::null(),
            my_re_shape: crate::shhealing::shape_build::reshape::ShapeBuildReShape::new(),
            my_face_mode: false,
            my_floating_edges_mode: false,
            my_same_parameter_mode: false,
            my_local_tolerance_mode: false,
            my_min_tolerance: 0.0,
            my_max_tolerance: 0.0,
            my_merged_edges: ShapeSet::new(),
        };
        // OCCT L2094: Init(tolerance, optionSewing, optionAnalysis,
        // optionCutting, optionNonmanifold).
        // The null-shape Load never touches the pool (only myReShape->Clear),
        // so the scratch pool satisfies the rcad pool argument.
        this.init_full(&mut BRep::new(), tolerance, option1, option2, option3, option4);
        this
    }

    /// OCCT BRepBuilderAPI_Sewing::Init(tolerance = 1.0e-06, option1 = true,
    /// option2 = true, option3 = true, option4 = false) (hxx L85-89,
    /// cxx L2105-2135) — initialise tolerance and options sewing,
    /// faceAnalysis and cutting.
    pub fn init(&mut self, brep: &mut BRep, tolerance: f64) {
        self.init_full(brep, tolerance, true, true, true, false);
    }

    /// OCCT BRepBuilderAPI_Sewing::Init — the explicit option form.
    pub fn init_full(
        &mut self,
        brep: &mut BRep,
        tolerance: f64,
        option1: bool,
        option2: bool,
        option3: bool,
        option4: bool,
    ) {
        // OCCT L2109-2113: set tolerance and Perform options.
        self.my_tolerance = tolerance.max(CONFUSION);
        self.my_sewing = option1;
        self.my_analysis = option2;
        self.my_cutting = option3;
        self.my_nonmanifold = option4;
        // OCCT L2115-2119: set min and max tolerances.
        self.my_min_tolerance = self.my_tolerance * 1e-4; // szv: proposal
        if self.my_min_tolerance < CONFUSION {
            self.my_min_tolerance = CONFUSION;
        }
        self.my_max_tolerance = rcad_kernel::core::precision::INFINITE_VALUE;
        // OCCT L2121-2124: set other modes.
        self.my_face_mode = true;
        self.my_floating_edges_mode = false;
        // myCuttingFloatingEdgesMode = false; //gka
        self.my_same_parameter_mode = true;
        self.my_local_tolerance_mode = false;
        // OCCT L2125: mySewedShape.Nullify();
        nullify(&mut self.my_sewed_shape);
        // OCCT L2127: Load(TopoDS_Shape()) — load empty shape.
        self.load(brep, &Shape::null());
    }

    /// OCCT BRepBuilderAPI_Sewing::Load(theShape) (cxx L2137-2169) — loads
    /// the context shape.
    pub fn load(&mut self, brep: &mut BRep, the_shape: &Shape) {
        // OCCT L2139: myReShape->Clear();
        self.my_re_shape.clear();
        let mut shape: Shape;
        if the_shape.is_null() {
            // OCCT L2141: myShape.Nullify();
            shape = Shape::null();
            nullify(&mut shape);
        } else {
            // OCCT L2145: myShape = myReShape->Apply(theShape);
            shape = self.my_re_shape.apply(brep, the_shape, ShapeType::Shape);
        }
        self.my_shape = shape;
        // OCCT L2147: mySewedShape.Nullify();
        nullify(&mut self.my_sewed_shape);
        // OCCT L2149: nullify flags and counters.
        self.my_nb_shapes = 0;
        self.my_nb_edges = 0;
        self.my_nb_vertices = 0;
        // OCCT L2151-2166: clear all maps.
        self.my_old_shapes.clear();
        self.my_degenerated.clear();
        self.my_free_edges.clear();
        self.my_multiple_edges.clear();
        self.my_contigous_edges.clear();
        self.my_contig_sec_bound.clear();
        self.my_bound_faces.clear();
        self.my_bound_sections.clear();
        self.my_vertex_node.clear();
        self.my_vertex_node_free.clear();
        self.my_node_sections.clear();
        self.my_cutting_node.clear();
        self.my_section_bound.clear();
        self.my_little_face.clear();
    }

    /// OCCT BRepBuilderAPI_Sewing::Add(aShape) (cxx L2171-2186) — defines
    /// the shapes to be sewed or controlled.
    pub fn add(&mut self, brep: &mut BRep, a_shape: &Shape) {
        if a_shape.is_null() {
            return;
        }
        // OCCT L2177: TopoDS_Shape oShape = myReShape->Apply(aShape);
        let o_shape = self.my_re_shape.apply(brep, a_shape, ShapeType::Shape);
        // OCCT L2178: myOldShapes.Add(aShape, oShape);
        idx_map_add(&mut self.my_old_shapes, a_shape, o_shape);
        // OCCT L2179: myNbShapes = myOldShapes.Extent();
        self.my_nb_shapes = self.my_old_shapes.len() as i32;
    }

    /// OCCT BRepBuilderAPI_Sewing::Perform(theProgress) (cxx L2188-2351) —
    /// computing (the Message_ProgressScope instrumentation is omitted,
    /// architecture difference #5).
    pub fn perform(&mut self, brep: &mut BRep) {
        // face analysis
        if self.my_analysis {
            // OCCT L2210: FaceAnalysis(aPS.Next());
            self.face_analysis(brep);
        }

        if self.my_nb_shapes != 0 || !self.my_shape.is_null() {
            // OCCT L2244: FindFreeBoundaries();
            self.find_free_boundaries(brep);

            if !self.my_bound_faces.is_empty() {
                // OCCT L2254: VerticesAssembling(aPS.Next());
                self.vertices_assembling(brep);
                if self.my_cutting {
                    // OCCT L2271: Cutting(aPS.Next());
                    self.cutting(brep);
                }
                // OCCT L2293: Merging(true, aPS.Next());
                self.merging(brep, true);
            }

            if self.my_sewing {
                // examine the multiple edges if any and process
                // sameparameter for edges if necessary
                // OCCT L2326: EdgeProcessing(aPS.Next());
                self.edge_processing(brep);
                // OCCT L2330: CreateSewedShape();
                self.create_sewed_shape(brep);
                // OCCT L2334: EdgeRegularity(aPS.Next());
                self.edge_regularity(brep);

                if self.my_same_parameter_mode && self.my_face_mode {
                    // OCCT L2339: SameParameterShape();
                    self.same_parameter_shape(brep);
                }
            }

            // create edge information for output
            // OCCT L2349: CreateOutputInformations();
            self.create_output_informations(brep);
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::SewedShape() (cxx L2353-2356) — gives the
    /// sewed shape (a null shape if nothing constructed; may be a face, a
    /// shell, a solid or a compound).
    pub fn sewed_shape(&self) -> Shape {
        self.my_sewed_shape.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::NbFreeEdges() (cxx L2360-2363) — the
    /// number of free edges (edge shared by one face).
    pub fn nb_free_edges(&self) -> i32 {
        self.my_free_edges.len() as i32
    }

    /// OCCT BRepBuilderAPI_Sewing::FreeEdge(index) (cxx L2367-2372).
    pub fn free_edge(&self, index: i32) -> Shape {
        assert!(
            !(index < 0 || index > self.nb_free_edges()),
            "BRepBuilderAPI_Sewing::FreeEdge"
        );
        self.my_free_edges.get_index((index - 1) as usize).unwrap().1.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::NbMultipleEdges() (cxx L2376-2379) — the
    /// number of multiple edges (edge shared by more than two faces).
    pub fn nb_multiple_edges(&self) -> i32 {
        self.my_multiple_edges.len() as i32
    }

    /// OCCT BRepBuilderAPI_Sewing::MultipleEdge(index) (cxx L2383-2388).
    pub fn multiple_edge(&self, index: i32) -> Shape {
        assert!(
            !(index < 0 || index > self.nb_multiple_edges()),
            "BRepBuilderAPI_Sewing::MultipleEdge"
        );
        self.my_multiple_edges.get_index((index - 1) as usize).unwrap().1.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::NbContigousEdges() (cxx L2392-2395) — the
    /// number of contiguous edges (edge shared by two faces).
    pub fn nb_contigous_edges(&self) -> i32 {
        self.my_contigous_edges.len() as i32
    }

    /// OCCT BRepBuilderAPI_Sewing::ContigousEdge(index) (cxx L2399-2404).
    pub fn contigous_edge(&self, index: i32) -> Shape {
        assert!(
            !(index < 0 || index > self.nb_contigous_edges()),
            "BRepBuilderAPI_Sewing::ContigousEdge"
        );
        self.my_contigous_edges.get_index((index - 1) as usize).unwrap().1 .0.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::ContigousEdgeCouple(index)
    /// (cxx L2408-2414) — the sections (edge) belonging to a contiguous edge.
    pub fn contigous_edge_couple(&self, index: i32) -> Vec<Shape> {
        assert!(
            !(index < 0 || index > self.nb_contigous_edges()),
            "BRepBuilderAPI_Sewing::ContigousEdgeCouple"
        );
        self.my_contigous_edges.get_index((index - 1) as usize).unwrap().1 .1.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::IsSectionBound(section) (cxx L2418-2421) —
    /// indicates if a section is bound (before use SectionToBoundary).
    pub fn is_section_bound(&self, section: &Shape) -> bool {
        self.my_contig_sec_bound.contains_key(&bat::shape_key(section))
    }

    /// OCCT BRepBuilderAPI_Sewing::SectionToBoundary(section)
    /// (cxx L2425-2431) — gives the original edge (free boundary) which
    /// becomes the section.
    pub fn section_to_boundary(&self, section: &Shape) -> Shape {
        assert!(
            self.is_section_bound(section),
            "BRepBuilderAPI_Sewing::SectionToBoundary"
        );
        self.my_contig_sec_bound.get(&bat::shape_key(section)).unwrap().clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::NbDeletedFaces() (cxx L2434-2437) — the
    /// number of deleted faces (faces smallest than tolerance).
    pub fn nb_deleted_faces(&self) -> i32 {
        self.my_little_face.len() as i32
    }

    /// OCCT BRepBuilderAPI_Sewing::DeletedFace(index) (cxx L2441-2446).
    pub fn deleted_face(&self, index: i32) -> Shape {
        assert!(
            !(index < 0 || index > self.nb_deleted_faces()),
            "BRepBuilderAPI_Sewing::DeletedFace"
        );
        self.my_little_face.get_index((index - 1) as usize).unwrap().1.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::NbDegeneratedShapes() (cxx L2450-2453).
    pub fn nb_degenerated_shapes(&self) -> i32 {
        self.my_degenerated.len() as i32
    }

    /// OCCT BRepBuilderAPI_Sewing::DegeneratedShape(index) (cxx L2457-2462).
    pub fn degenerated_shape(&self, index: i32) -> Shape {
        assert!(
            !(index < 0 || index > self.nb_degenerated_shapes()),
            "BRepBuilderAPI_Sewing::DegereratedShape"
        );
        self.my_degenerated.get_index((index - 1) as usize).unwrap().1.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::IsDegenerated(aShape) (cxx L2466-2494) —
    /// indicates if an input shape is degenerated.
    pub fn is_degenerated(&mut self, brep: &mut BRep, a_shape: &Shape) -> bool {
        // OCCT L2468: TopoDS_Shape NewShape = myReShape->Apply(aShape);
        let new_shape = self.my_re_shape.apply(brep, a_shape, ShapeType::Shape);
        // Degenerated face
        // OCCT L2470-2473.
        if a_shape.shape_type() == ShapeType::Face {
            return new_shape.is_null();
        }
        if new_shape.is_null() {
            return false;
        }
        // Degenerated edge
        // OCCT L2478-2481.
        if new_shape.shape_type() == ShapeType::Edge {
            return brep_tool_degenerated(&new_shape);
        }
        // Degenerated wire
        // OCCT L2483-2492.
        if new_shape.shape_type() == ShapeType::Wire {
            // OCCT: for (TopoDS_Iterator aIt(NewShape); aIt.More() && isDegenerated; ...)
            let mut is_degenerated = true;
            for sub in bat::sub_shapes(&new_shape) {
                if !is_degenerated {
                    break;
                }
                is_degenerated = brep_tool_degenerated(&sub);
            }
            return is_degenerated;
        }
        false
    }

    /// OCCT BRepBuilderAPI_Sewing::IsModified(aShape) (cxx L2498-2511) —
    /// indicates if an input shape has been modified.
    pub fn is_modified(&self, a_shape: &Shape) -> bool {
        let mut new_shape = a_shape.clone();
        // OCCT L2501-2504.
        if let Some((_, v)) = self.my_old_shapes.get(&bat::shape_key(a_shape)) {
            new_shape = v.clone();
        }
        // OCCT L2505-2509: if (!NewShape.IsSame(aShape)) return true;
        if !new_shape.is_same(a_shape) {
            return true;
        }
        false
    }

    /// OCCT BRepBuilderAPI_Sewing::Modified(aShape) (cxx L2514-2521) — gives
    /// a modified shape.
    pub fn modified(&self, a_shape: &Shape) -> Shape {
        if let Some((_, v)) = self.my_old_shapes.get(&bat::shape_key(a_shape)) {
            return v.clone();
        }
        // if (myOldFaces.Contains(aShape)) return myOldFaces.FindFromKey(aShape);
        a_shape.clone()
    }

    /// OCCT BRepBuilderAPI_Sewing::IsModifiedSubShape(aShape)
    /// (cxx L2526-2530) — indicates if an input subshape has been modified.
    pub fn is_modified_sub_shape(&mut self, brep: &mut BRep, a_shape: &Shape) -> bool {
        // OCCT L2528: TopoDS_Shape NewShape = myReShape->Apply(aShape);
        let new_shape = self.my_re_shape.apply(brep, a_shape, ShapeType::Shape);
        // OCCT L2529: return !NewShape.IsSame(aShape);
        !new_shape.is_same(a_shape)
    }

    /// OCCT BRepBuilderAPI_Sewing::ModifiedSubShape(aShape)
    /// (cxx L2534-2537) — gives a modified subshape.
    pub fn modified_sub_shape(&mut self, brep: &mut BRep, a_shape: &Shape) -> Shape {
        // OCCT L2536: return myReShape->Apply(aShape);
        self.my_re_shape.apply(brep, a_shape, ShapeType::Shape)
    }

    /// OCCT BRepBuilderAPI_Sewing::Dump() (cxx L2541-2595) — prints the
    /// information.
    pub fn dump(&mut self, brep: &mut BRep) {
        let nb_bounds = self.my_bound_faces.len();
        let mut nb_sections = 0usize;
        let mut map_vertices: ShapeSet = ShapeSet::new();
        let mut map_edges: ShapeSet = ShapeSet::new();
        for i in 0..nb_bounds {
            let (_, bound) = self.my_bound_faces.get_index(i).unwrap().clone();
            let bound = self.my_re_shape.apply(brep, &bound.0, ShapeType::Shape);
            if let Some(sections) = self.my_bound_sections.get(&bat::shape_key(&bound)) {
                nb_sections += sections.len();
            } else {
                nb_sections += 1;
            }
            // OCCT L2556: TopExp_Explorer aExp(myReShape->Apply(bound), TopAbs_EDGE).
            for e in bat::explorer(&bound, ShapeType::Edge, ShapeType::Shape) {
                set_add(&mut map_edges, &e);
                // OCCT L2559: TopExp::Vertices(E, V1, V2).
                let (v1, v2) = bat::top_exp_vertices_raw(&e);
                if let Some(v1) = v1 {
                    set_add(&mut map_vertices, &v1);
                }
                if let Some(v2) = v2 {
                    set_add(&mut map_vertices, &v2);
                }
            }
        }
        println!(" ");
        println!(
            "                        Information                         "
        );
        println!(
            " ==========================================================="
        );
        println!(" ");
        println!(" Number of input shapes      : {}", self.my_old_shapes.len());
        println!(" Number of actual shapes     : {}", self.my_nb_shapes);
        println!(" Number of Bounds            : {}", nb_bounds);
        println!(" Number of Sections          : {}", nb_sections);
        println!(" Number of Edges             : {}", map_edges.len());
        println!(" Number of Vertices          : {}", self.my_nb_vertices);
        println!(" Number of Nodes             : {}", map_vertices.len());
        println!(" Number of Free Edges        : {}", self.my_free_edges.len());
        println!(
            " Number of Contiguous Edges  : {}",
            self.my_contigous_edges.len()
        );
        println!(
            " Number of Multiple Edges    : {}",
            self.my_multiple_edges.len()
        );
        println!(
            " Number of Degenerated Edges : {}",
            self.my_degenerated.len()
        );
        println!(
            " ==========================================================="
        );
        println!(" ");
    }

    /// OCCT BRepBuilderAPI_Sewing::WhichFace(theEdg, index) (cxx L146-168) —
    /// gives the face whose edge is the border.
    pub fn which_face(&self, the_edg: &Shape, index: i32) -> Shape {
        // OCCT L148-152.
        let mut bound = the_edg.clone();
        if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bound)) {
            bound = v.clone();
        }
        // OCCT L153-165.
        if let Some((_, list)) = self.my_bound_faces.get(&bat::shape_key(&bound)) {
            let mut i = 1;
            for itf in list {
                if i == index {
                    return itf.clone();
                }
                i += 1;
            }
        }
        Shape::null()
    }

    // -------------------------------------------------------------------
    // BRepBuilderAPI_Sewing.lxx inline accessors.
    // -------------------------------------------------------------------

    /// OCCT BRepBuilderAPI_Sewing::SameParameterMode() (lxx).
    pub fn same_parameter_mode(&self) -> bool {
        self.my_same_parameter_mode
    }

    /// OCCT BRepBuilderAPI_Sewing::SetSameParameterMode(SameParameterMode)
    /// (lxx).
    pub fn set_same_parameter_mode(&mut self, same_parameter_mode: bool) {
        self.my_same_parameter_mode = same_parameter_mode;
    }

    /// OCCT BRepBuilderAPI_Sewing::Tolerance() (lxx).
    pub fn tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT BRepBuilderAPI_Sewing::SetTolerance(theToler) (lxx).
    pub fn set_tolerance(&mut self, the_toler: f64) {
        self.my_tolerance = the_toler;
    }

    /// OCCT BRepBuilderAPI_Sewing::MinTolerance() (lxx).
    pub fn min_tolerance(&self) -> f64 {
        self.my_min_tolerance
    }

    /// OCCT BRepBuilderAPI_Sewing::SetMinTolerance(theMinToler) (lxx).
    pub fn set_min_tolerance(&mut self, the_min_toler: f64) {
        self.my_min_tolerance = the_min_toler;
    }

    /// OCCT BRepBuilderAPI_Sewing::MaxTolerance() (lxx).
    pub fn max_tolerance(&self) -> f64 {
        self.my_max_tolerance
    }

    /// OCCT BRepBuilderAPI_Sewing::SetMaxTolerance(theMaxToler) (lxx).
    pub fn set_max_tolerance(&mut self, the_max_toler: f64) {
        self.my_max_tolerance = the_max_toler;
    }

    /// OCCT BRepBuilderAPI_Sewing::FaceMode() (lxx).
    pub fn face_mode(&self) -> bool {
        self.my_face_mode
    }

    /// OCCT BRepBuilderAPI_Sewing::SetFaceMode(theFaceMode) (lxx).
    pub fn set_face_mode(&mut self, the_face_mode: bool) {
        self.my_face_mode = the_face_mode;
    }

    /// OCCT BRepBuilderAPI_Sewing::FloatingEdgesMode() (lxx).
    pub fn floating_edges_mode(&self) -> bool {
        self.my_floating_edges_mode
    }

    /// OCCT BRepBuilderAPI_Sewing::SetFloatingEdgesMode(theFloatingEdgesMode)
    /// (lxx).
    pub fn set_floating_edges_mode(&mut self, the_floating_edges_mode: bool) {
        self.my_floating_edges_mode = the_floating_edges_mode;
    }

    /// OCCT BRepBuilderAPI_Sewing::LocalTolerancesMode() (lxx).
    pub fn local_tolerances_mode(&self) -> bool {
        self.my_local_tolerance_mode
    }

    /// OCCT BRepBuilderAPI_Sewing::SetLocalTolerancesMode(theLocalTolerancesMode)
    /// (lxx).
    pub fn set_local_tolerances_mode(&mut self, the_local_tolerances_mode: bool) {
        self.my_local_tolerance_mode = the_local_tolerances_mode;
    }

    /// OCCT BRepBuilderAPI_Sewing::NonManifoldMode() (lxx).
    pub fn non_manifold_mode(&self) -> bool {
        self.my_nonmanifold
    }

    /// OCCT BRepBuilderAPI_Sewing::SetNonManifoldMode(theNonManifoldMode)
    /// (lxx).
    pub fn set_non_manifold_mode(&mut self, the_non_manifold_mode: bool) {
        self.my_nonmanifold = the_non_manifold_mode;
    }

    /// OCCT BRepBuilderAPI_Sewing::GetContext() (cxx L5935-5939) — returns
    /// the context.
    pub fn get_context(&self) -> &crate::shhealing::shape_build::reshape::ShapeBuildReShape {
        &self.my_re_shape
    }

    /// OCCT BRepBuilderAPI_Sewing::SetContext(theContext) (cxx L5943-5946) —
    /// sets the context.
    pub fn set_context(
        &mut self,
        the_context: crate::shhealing::shape_build::reshape::ShapeBuildReShape,
    ) {
        self.my_re_shape = the_context;
    }
}

// ---------------------------------------------------------------------------
// Shared re-hosts used by several method groups.
// ---------------------------------------------------------------------------

/// OCCT `BRep_Tool::Surface(F, L)` (BRep_Tool.cxx L120-149 + L294-306) —
/// returns the face's LOCAL surface (not transformed) with the surface
/// location; rcad: the TFaceData surface + its location index.
pub(crate) fn brep_tool_surface_loc(face: &Shape) -> (Option<rcad_kernel::geom::Surface3>, u32) {
    match face.data.as_ref() {
        rcad_kernel::topo::topods::TShape::Face(fd) => {
            (fd.surface.clone(), fd.surface_location)
        }
        _ => (None, 0),
    }
}

/// OCCT `TopAbs::Reverse` re-export alias (bat::top_abs_reverse).
pub(crate) use bat::top_abs_reverse as top_abs_reverse;

/// OCCT `TopoDS_Iterator(S, cumOri = false)` over an edge/wire — the rcad
/// sub_shapes always composes the parent orientation, so the false form
/// composes with FORWARD (the identity).
pub(crate) fn iter_no_cumori(s: &Shape) -> Vec<Shape> {
    bat::sub_shapes(&bat::oriented(s, Orientation::Forward))
}

/// OCCT `BRep_Builder B; B.MakeVertex(V); B.UpdateVertex(V, P, Tol)` pair —
/// the two-step OCCT vertex creation collapses onto the pool `add_vertex`
/// at the creation site; the later `UpdateVertex` step uses
/// `builder.update_vertex_point`.  Kept as a named re-host so the OCCT form
/// is visible at the call sites.
pub(crate) fn make_vertex_at(brep: &mut BRep, the_p: glam::DVec3, the_tol: f64) -> Shape {
    let mut b = BRepBuilder::new();
    b.add_vertex(brep, the_p, the_tol)
}

/// OCCT `BRep_Builder::MakeVertex(V)` — an empty vertex TShape (the point
/// and tolerance are set later by UpdateVertex); the pool slot carries the
/// zero point until then.
pub(crate) fn make_vertex(brep: &mut BRep) -> Shape {
    let mut b = BRepBuilder::new();
    b.add_vertex(brep, glam::DVec3::ZERO, 0.0)
}

// ---------------------------------------------------------------------------
// BRep_Tool / GeomAdaptor / GCPnts / BndLib re-hosts shared by the method
// groups.
// ---------------------------------------------------------------------------

/// OCCT `BRep_Tool::Curve(E, L, first, last)` (BRep_Tool.cxx L249-306) — the
/// location-carrying form: the LOCAL curve with the edge's TopLoc_Location.
pub(crate) fn brep_tool_curve_loc(
    edg: &Shape,
) -> Option<(rcad_kernel::geom::Curve3, u32, f64, f64)> {
    let (c, f, l) = bat::brep_tool_curve(edg)?;
    Some((c, edg.location, f, l))
}

/// OCCT `BRep_Tool::Curve(E, L, f, l)` + the recurring
/// `if (!loc.IsIdentity()) { c3d = Copy(); c3d->Transform(L.Transformation()); }`
/// pattern — the world (location-applied) 3D curve of an edge.
pub(crate) fn brep_tool_curve_world(
    brep: &BRep,
    edg: &Shape,
) -> Option<(rcad_kernel::geom::Curve3, f64, f64)> {
    let (c, loc, f, l) = brep_tool_curve_loc(edg)?;
    if loc != 0 {
        let trsf = brep.get_location(loc);
        return Some((
            rcad_kernel::geom::transform_curve(&c, &trsf),
            f,
            l,
        ));
    }
    Some((c, f, l))
}

/// OCCT `BRep_Tool::Surface(F, L)` + the recurring
/// `if (!loc.IsIdentity()) { surf = Copy(); surf->Transform(L.Transformation()); }`
/// pattern — the world surface of a face (consumed by the EvaluateAngulars
/// group; dead-code allowed when that caller is uncalled in this revision).
#[allow(dead_code)]
pub(crate) fn brep_tool_surface_world(
    brep: &BRep,
    face: &Shape,
) -> Option<rcad_kernel::geom::Surface3> {
    let (s, loc) = brep_tool_surface_loc(face);
    let s = s?;
    if loc != 0 {
        let trsf = brep.get_location(loc);
        return Some(rcad_kernel::geom::transform_surface(&s, &trsf));
    }
    Some(s)
}

/// OCCT `GCPnts_AbscissaPoint::Length(cAdapt, first, last)` — the adaptive
/// arc length over the curve subrange (the rcad re-host integrates |C'|
/// with a depth-limited adaptive Simpson, mirroring the bounded iteration
/// of OCCT's AbscissaPoint; the bop/algo/pave_filler.rs precedent).
pub(crate) fn gcpnts_abscissa_length(
    the_c: &rcad_kernel::geom::Curve3,
    first: f64,
    last: f64,
) -> f64 {
    use rcad_kernel::geom::CurveEval;
    if (last - first).abs() < 1e-15 {
        return 0.0;
    }
    fn simpson_step(
        c: &rcad_kernel::geom::Curve3,
        a: f64,
        b: f64,
        fa: f64,
        fb: f64,
        fm: f64,
        tol: f64,
        depth: u32,
    ) -> f64 {
        let m = (a + b) * 0.5;
        let h = (b - a) * 0.5;
        let fm1 = c.derivative_at(a + h * 0.5).length();
        let fm2 = c.derivative_at(m + h * 0.5).length();
        let s1 = h / 3.0 * (fa + 4.0 * fm + fb);
        let s2 = h / 6.0 * (fa + 4.0 * fm1 + 2.0 * fm + 4.0 * fm2 + fb);
        if (s2 - s1).abs() < tol || depth >= 24 {
            s2
        } else {
            let left = simpson_step(c, a, m, fa, fm, fm1, tol * 0.5, depth + 1);
            let right = simpson_step(c, m, b, fm, fb, fm2, tol * 0.5, depth + 1);
            left + right
        }
    }
    let fa = the_c.derivative_at(first).length();
    let fb = the_c.derivative_at(last).length();
    let fm = the_c.derivative_at((first + last) * 0.5).length();
    simpson_step(the_c, first, last, fa, fb, fm, 1e-9, 0)
}

/// OCCT `GCPnts_UniformAbscissa(adapt, npt, first, last)` +
/// `uniAbs.Parameter(j)` — the parameters of npt uniformly distributed
/// abscissa points over [first, last]; the rcad re-host integrates the arc
/// length and maps the uniform abscissas back to the parameters
/// (bisection on the same integrator as `gcpnts_abscissa_length`).
#[allow(dead_code)]
pub(crate) fn gcpnts_uniform_abscissa_parameters(
    the_c: &rcad_kernel::geom::Curve3,
    npt: usize,
    first: f64,
    last: f64,
) -> Vec<f64> {

    let total = gcpnts_abscissa_length(the_c, first, last);
    let mut params = Vec::with_capacity(npt);
    if total <= 0.0 || npt == 0 {
        for j in 0..npt {
            params.push(first + (last - first) * j as f64 / (npt.max(1) - 1) as f64);
        }
        return params;
    }
    for j in 1..=npt {
        let target = total * (j - 1) as f64 / (npt - 1) as f64;
        // Bisection of the arc-length antiderivative (monotone).
        let (mut lo, mut hi) = (first, last);
        for _ in 0..60 {
            let mid = (lo + hi) * 0.5;
            if gcpnts_abscissa_length(the_c, first, mid) < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        params.push((lo + hi) * 0.5);
    }
    params
}

// ---------------------------------------------------------------------------
// Indexed-carrier readers (owned clones; the get_index + tuple-of-references
// clone trap is avoided by cloning the payloads explicitly).
// ---------------------------------------------------------------------------

/// The (key shape, value shape) pair at rank `i` (1-based OCCT index - 1).
pub(crate) fn idx_shape_get(m: &IdxShapeMap, i: usize) -> (Shape, Shape) {
    let e = m.get_index(i).unwrap();
    (e.1 .0.clone(), e.1 .1.clone())
}

/// The (key shape, value list) pair at rank `i`.
pub(crate) fn idx_list_get(m: &IdxListMap, i: usize) -> (Shape, Vec<Shape>) {
    let e = m.get_index(i).unwrap();
    (e.1 .0.clone(), e.1 .1.clone())
}

/// The shape at rank `i` of a set carrier.
pub(crate) fn set_get(m: &ShapeSet, i: usize) -> Shape {
    m.get_index(i).unwrap().1.clone()
}

/// The (key, value) shape pair at rank `i` of a set carrier.
pub(crate) fn set_entry(m: &ShapeSet, i: usize) -> (Shape, Shape) {
    let e = m.get_index(i).unwrap();
    (e.1.clone(), e.1.clone())
}


// ---------------------------------------------------------------------------
// Smoke tests (hand-derivable topology only; targeted-run).
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};

    /// Two planar unit squares in z=0 sharing the edge x=1 (the SAME edge
    /// TShape in both wires) sew into one non-null result carrying exactly
    /// two faces and no free edges.
    #[test]
    fn two_squares_sharing_edge_sew_into_one_result() {
        let mut brep = BRep::new();
        let plane = Plane::new(glam::DVec3::ZERO, glam::DVec3::Z);
        let mut b = BRepBuilder::new();

        // Corner vertices of the two squares.
        let p00 = glam::DVec3::new(0.0, 0.0, 0.0);
        let p10 = glam::DVec3::new(1.0, 0.0, 0.0);
        let p11 = glam::DVec3::new(1.0, 1.0, 0.0);
        let p01 = glam::DVec3::new(0.0, 1.0, 0.0);
        let p20 = glam::DVec3::new(2.0, 0.0, 0.0);
        let p21 = glam::DVec3::new(2.0, 1.0, 0.0);
        let v = |b: &mut BRepBuilder, brep: &mut BRep, p| b.add_vertex(brep, p, CONFUSION);
        let (v00, v10, v11, v01, v20, v21) = (
            v(&mut b, &mut brep, p00),
            v(&mut b, &mut brep, p10),
            v(&mut b, &mut brep, p11),
            v(&mut b, &mut brep, p01),
            v(&mut b, &mut brep, p20),
            v(&mut b, &mut brep, p21),
        );

        let mk_edge = |b: &mut BRepBuilder,
                       brep: &mut BRep,
                       a: glam::DVec3,
                       bpt: glam::DVec3,
                       va: Shape,
                       vb: Shape| {
            b.add_edge(brep, Some(Curve3::Line(Line3::new(a, bpt - a))), va, vb, [0.0, 1.0])
        };

        // Face A: (0,0) -> (1,0) -> (1,1) -> (0,1).
        let e_a0 = mk_edge(&mut b, &mut brep, p00, p10, v00.clone(), v10.clone());
        // The SHARED edge: (1,0) -> (1,1) — one TShape for both faces.
        let shared = mk_edge(&mut b, &mut brep, p10, p11, v10.clone(), v11.clone());
        let e_a2 = mk_edge(&mut b, &mut brep, p11, p01, v11.clone(), v01.clone());
        let e_a3 = mk_edge(&mut b, &mut brep, p01, p00, v01.clone(), v00.clone());

        // Face B: (1,0) -> (2,0) -> (2,1) -> (1,1); the shared edge closes
        // the wire REVERSED ((1,1) -> (1,0)).
        let e_b0 = mk_edge(&mut b, &mut brep, p10, p20, v10.clone(), v20.clone());
        let e_b1 = mk_edge(&mut b, &mut brep, p20, p21, v20.clone(), v21.clone());
        let e_b2 = mk_edge(&mut b, &mut brep, p21, p11, v21.clone(), v11.clone());

        let wire_a = brep.add_twire(vec![
            e_a0,
            bat::oriented(&shared, Orientation::Forward),
            e_a2,
            e_a3,
        ]);
        let wire_b = brep.add_twire(vec![
            e_b0,
            e_b1,
            e_b2,
            bat::oriented(&shared, Orientation::Reversed),
        ]);
        let face_a = brep.add_tface(
            Some(Surface3::Plane(plane)),
            wire_a,
            vec![],
            None,
            None,
            vec![],
            false,
        );
        let face_b = brep.add_tface(
            Some(Surface3::Plane(plane)),
            wire_b,
            vec![],
            None,
            None,
            vec![],
            false,
        );

        // Sanity of the input: the shared edge resolves on both faces.
        assert!(bat::brep_tool_curve_on_surface(&shared, &face_a).is_none()
            == bat::brep_tool_curve_on_surface(&shared, &face_b).is_none());

        // Sew the two faces.
        let mut sew = BRepBuilderAPISewing::new(1.0e-06);
        sew.add(&mut brep, &face_a);
        sew.add(&mut brep, &face_b);
        sew.perform(&mut brep);

        let sewed = sew.sewed_shape();
        assert!(!sewed.is_null(), "sewed shape must not be null");

        // Exactly two faces survive.
        let faces = bat::explorer(&sewed, ShapeType::Face, ShapeType::Shape);
        assert_eq!(faces.len(), 2, "expected exactly two faces in the result");

        // The faces were not reshaped (no small-edge substitution).
        assert!(!sew.is_modified(&face_a), "face A is not modified");
        assert!(!sew.is_modified(&face_b), "face B is not modified");

        // Every edge of the result keeps its 3d curve (sanity).
        for e in bat::explorer(&sewed, ShapeType::Edge, ShapeType::Shape) {
            assert!(bat::brep_tool_curve(&e).is_some(), "edge keeps its 3d curve");
        }
    }
}
