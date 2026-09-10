// OCCT BRepFeat_SplitShape.hxx L17-125 + BRepFeat_SplitShape.cxx L17-90 +
// BRepFeat_SplitShape.lxx L17-87 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_SplitShape.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_SplitShape.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_SplitShape.lxx
//
// OCCT inheritance chain (BRepFeat_SplitShape.hxx L50):
//   BRepFeat_SplitShape : BRepBuilderAPI_MakeShape
//     BRepBuilderAPI_MakeShape : BRepBuilderAPI_Command
// Rust has no inheritance -> composition + delegation. The
// BRepBuilderAPI_MakeShape base sub-object is carried by the
// BRepFeatSplitShape fields listed below.
//
// Field mapping:
//   BRepBuilderAPI_Command::myDone          -> my_done
//   BRepBuilderAPI_MakeShape::myShape       -> my_shape (None = null)
//   BRepBuilderAPI_MakeShape::myGenerated   -> my_generated
//   BRepFeat_SplitShape::mySShape           -> my_sshape
//   BRepFeat_SplitShape::myWOnShape         -> my_w_on_shape (Option carries
//                                              the OCCT null handle)
//   BRepFeat_SplitShape::myRight (mutable)  -> my_right
//
// Architecture differences (referenced from the affected functions):
// 1. LocOpe_Spliter (LocOpe_Spliter.cxx L17-718) and LocOpe_WiresOnShape
//    (1623 lines) are deferred translations (Stage 3d; see the
//    loc_ope_spliter.rs header: Perform needs LocOpe_WiresOnShape,
//    LocOpe_BuildWires, BRepTools_Substitution, GeomAPI_ProjectPointOnCurve).
//    Both classes are carried below as API-surface carriers with the OCCT
//    member sets (LocOpe_Spliter.hxx L64-72, LocOpe_WiresOnShape.hxx
//    L109-114) and the methods BRepFeat_SplitShape consumes; their bodies
//    land with Stage 3d (same carrier model as BRepPrimCylinder in
//    brep_feat_make_cylindrical_hole.rs). Until then my_sshape.my_done
//    stays false and Build does not produce a result.

use crate::feat::brep_feat_builder::{explorer, OcctShapeMap};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
use std::collections::HashMap;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

// ---------------------------------------------------------------------------
// API-surface carriers for the deferred LocOpe classes (architecture
// difference #1).
// ---------------------------------------------------------------------------

/// OCCT LocOpe_Spliter (LocOpe_Spliter.hxx L30-72) — deferred body (Stage
/// 3d); the member set and the surface consumed by BRepFeat_SplitShape are
/// carried below.
pub struct LocOpeSpliter {
    my_shape: Shape,                          // OCCT: myShape
    my_done: bool,                            // OCCT: myDone
    my_res: Option<Shape>,                    // OCCT: myRes
    my_map: HashMap<(u64, u32), Vec<Shape>>,  // OCCT: myMap
    my_dleft: Vec<Shape>,                     // OCCT: myDLeft
    my_left: Vec<Shape>,                      // OCCT: myLeft
}

impl Default for LocOpeSpliter {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeSpliter {
    /// OCCT LocOpe_Spliter::LocOpe_Spliter() — empty constructor.
    pub fn new() -> Self {
        LocOpeSpliter {
            my_shape: Shape::null(),
            my_done: false,
            my_res: None,
            my_map: HashMap::new(),
            my_dleft: Vec::new(),
            my_left: Vec::new(),
        }
    }

    /// OCCT LocOpe_Spliter::LocOpe_Spliter(S).
    pub fn new_with_shape(the_s: &Shape) -> Self {
        let mut sp = LocOpeSpliter::new();
        sp.init(the_s);
        sp
    }

    /// OCCT LocOpe_Spliter::Init(S) (lxx) — myShape = S. Deferred: the body
    /// lands with Stage 3d.
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_done = false;
        self.my_res = None;
        self.my_map.clear();
        self.my_dleft.clear();
        self.my_left.clear();
    }

    /// OCCT LocOpe_Spliter::Perform(PW) (cxx L63-718) — deferred body
    /// (Stage 3d, see the header note).
    pub fn perform(&mut self, _p_w: &mut LocOpeWiresOnShape) {
        // deferred: myDone stays false until the Stage 3d body lands.
    }

    /// OCCT LocOpe_Spliter::IsDone() (lxx).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT LocOpe_Spliter::ResultingShape() (lxx).
    pub fn resulting_shape(&self) -> Option<&Shape> {
        self.my_res.as_ref()
    }

    /// OCCT LocOpe_Spliter::DirectLeft() (cxx) — deferred body.
    pub fn direct_left(&self) -> &Vec<Shape> {
        &self.my_dleft
    }

    /// OCCT LocOpe_Spliter::Left() (cxx) — deferred body.
    pub fn left(&self) -> &Vec<Shape> {
        &self.my_left
    }

    /// OCCT LocOpe_Spliter::DescendantShapes(S) — deferred body; the myMap
    /// lookup returns the bound list (an absent key carries the OCCT absent
    /// DataMap entry).
    pub fn descendant_shapes(&mut self, the_s: &Shape) -> Option<&Vec<Shape>> {
        self.my_map.get(&shape_key(the_s))
    }
}

/// OCCT LocOpe_WiresOnShape (LocOpe_WiresOnShape.hxx L38-114) — deferred
/// body (Stage 3d); the member set and the surface consumed by
/// BRepFeat_SplitShape are carried below.
pub struct LocOpeWiresOnShape {
    // OCCT: myShape.
    my_shape: Shape,
    // OCCT: myMapEF (NCollection_IndexedDataMap<Shape, Shape>); values are
    // Option to carry the Nullify() of the LocOpe body (same model as
    // LocOpeGluer::my_map_ef).
    my_map_ef: indexmap::IndexMap<(u64, u32), (Shape, Option<Shape>)>,
    // OCCT: myFacesWithSection (NCollection_Map<Shape>).
    my_faces_with_section: OcctShapeMap,
    // OCCT: myCheckInterior.
    my_check_interior: bool,
    // OCCT: myMap (NCollection_DataMap<Shape, Shape>).
    my_map: HashMap<(u64, u32), Shape>,
    // OCCT: myDone.
    my_done: bool,
}

impl LocOpeWiresOnShape {
    /// OCCT LocOpe_WiresOnShape::LocOpe_WiresOnShape(S) (cxx) — deferred
    /// body; the shape is carried.
    pub fn new(the_s: &Shape) -> Self {
        LocOpeWiresOnShape {
            my_shape: the_s.clone(),
            my_map_ef: indexmap::IndexMap::new(),
            my_faces_with_section: OcctShapeMap::new(),
            my_check_interior: true,
            my_map: HashMap::new(),
            my_done: false,
        }
    }

    /// OCCT LocOpe_WiresOnShape::Init(S) (cxx) — deferred body.
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_map_ef.clear();
        self.my_faces_with_section.clear();
        self.my_map.clear();
        self.my_done = false;
    }

    /// OCCT LocOpe_WiresOnShape::Add(theEdges) (cxx) — deferred body; the
    /// OCCT bool return carries false until Stage 3d.
    pub fn add(&mut self, _the_edges: &[Shape]) -> bool {
        false
    }

    /// OCCT LocOpe_WiresOnShape::SetCheckInterior(ToCheckInterior) (lxx).
    pub fn set_check_interior(&mut self, the_to_check_interior: bool) {
        self.my_check_interior = the_to_check_interior;
    }

    /// OCCT LocOpe_WiresOnShape::Bind(W, F) (cxx) — deferred body.
    pub fn bind_wire_face(&mut self, _the_w: &Shape, _the_f: &Shape) {}

    /// OCCT LocOpe_WiresOnShape::Bind(Comp, F) (cxx) — deferred body.
    pub fn bind_compound_face(&mut self, _the_comp: &Shape, _the_f: &Shape) {}

    /// OCCT LocOpe_WiresOnShape::Bind(E, F) (cxx) — deferred body.
    pub fn bind_edge_face(&mut self, _the_e: &Shape, _the_f: &Shape) {}

    /// OCCT LocOpe_WiresOnShape::Bind(EfromW, EonFace) (cxx) — deferred body.
    pub fn bind_edge_edge(&mut self, _the_efrom_w: &Shape, _the_eon_face: &Shape) {}
}

// ---------------------------------------------------------------------------
// BRepFeat_SplitShape proper.
// ---------------------------------------------------------------------------

/// OCCT BRepFeat_SplitShape — splitting a shape with wires or edges through
/// local operations (BRepFeat_SplitShape.hxx L37-50).
pub struct BRepFeatSplitShape {
    my_done: bool,              // BRepBuilderAPI_Command::myDone
    my_shape: Option<Shape>,    // BRepBuilderAPI_MakeShape::myShape (None = null)
    // myGenerated is the MakeShape base member; BRepFeat_SplitShape never
    // fills it (Modified reads mySShape.DescendantShapes instead).
    #[allow(dead_code)]
    my_generated: Vec<Shape>,
    my_sshape: LocOpeSpliter,   // BRepFeat_SplitShape::mySShape
    my_w_on_shape: Option<LocOpeWiresOnShape>, // BRepFeat_SplitShape::myWOnShape
    my_right: Vec<Shape>,       // BRepFeat_SplitShape::myRight (mutable)
}

impl Default for BRepFeatSplitShape {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFeatSplitShape {
    /// OCCT BRepBuilderAPI_Command::Done().
    fn done(&mut self) {
        self.my_done = true;
    }

    /// OCCT BRepBuilderAPI_Command::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepFeat_SplitShape::BRepFeat_SplitShape() (lxx L21) — = default.
    pub fn new() -> Self {
        BRepFeatSplitShape {
            my_done: false,
            my_shape: None,
            my_generated: Vec::new(),
            my_sshape: LocOpeSpliter::new(),
            my_w_on_shape: None,
            my_right: Vec::new(),
        }
    }

    /// OCCT BRepFeat_SplitShape::BRepFeat_SplitShape(S) (lxx L25-29).
    pub fn new_with_shape(the_s: &Shape) -> Self {
        let mut s = BRepFeatSplitShape {
            my_done: false,
            my_shape: None,
            my_generated: Vec::new(),
            my_sshape: LocOpeSpliter::new_with_shape(the_s),
            my_w_on_shape: None,
            my_right: Vec::new(),
        };
        s.my_w_on_shape = Some(LocOpeWiresOnShape::new(the_s));
        s
    }

    /// OCCT BRepFeat_SplitShape::Add(theEdges) (lxx L33-36) — add splitting
    /// edges or wires for the whole initial shape.
    pub fn add(&mut self, the_edges: &[Shape]) -> bool {
        match self.my_w_on_shape.as_mut() {
            Some(w) => w.add(the_edges),
            None => false,
        }
    }

    /// OCCT BRepFeat_SplitShape::Init(S) (lxx L40-51).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_sshape.init(the_s);
        match self.my_w_on_shape.as_mut() {
            None => {
                self.my_w_on_shape = Some(LocOpeWiresOnShape::new(the_s));
            }
            Some(w) => {
                w.init(the_s);
            }
        }
    }

    /// OCCT BRepFeat_SplitShape::SetCheckInterior(ToCheckInterior)
    /// (lxx L55-58).
    pub fn set_check_interior(&mut self, the_to_check_interior: bool) {
        if let Some(w) = self.my_w_on_shape.as_mut() {
            w.set_check_interior(the_to_check_interior);
        }
    }

    /// OCCT BRepFeat_SplitShape::Add(W, F) (lxx L62-65).
    pub fn add_wire_on_face(&mut self, the_w: &Shape, the_f: &Shape) {
        if let Some(w) = self.my_w_on_shape.as_mut() {
            w.bind_wire_face(the_w, the_f);
        }
    }

    /// OCCT BRepFeat_SplitShape::Add(E, F) (lxx L69-72).
    pub fn add_edge_on_face(&mut self, the_e: &Shape, the_f: &Shape) {
        if let Some(w) = self.my_w_on_shape.as_mut() {
            w.bind_edge_face(the_e, the_f);
        }
    }

    /// OCCT BRepFeat_SplitShape::Add(Comp, F) (lxx L76-79).
    pub fn add_compound_on_face(&mut self, the_comp: &Shape, the_f: &Shape) {
        if let Some(w) = self.my_w_on_shape.as_mut() {
            w.bind_compound_face(the_comp, the_f);
        }
    }

    /// OCCT BRepFeat_SplitShape::Add(E, EOn) (lxx L83-86).
    pub fn add_edge_on_edge(&mut self, the_e: &Shape, the_e_on: &Shape) {
        if let Some(w) = self.my_w_on_shape.as_mut() {
            w.bind_edge_edge(the_e, the_e_on);
        }
    }

    /// OCCT BRepFeat_SplitShape::Build (cxx L24-33).
    pub fn build(&mut self) {
        // OCCT cxx L26: mySShape.Perform(myWOnShape).
        let my_sshape = &mut self.my_sshape;
        match self.my_w_on_shape.as_mut() {
            Some(w) => my_sshape.perform(w),
            None => {}
        }
        // OCCT cxx L27: if (mySShape.IsDone()).
        if self.my_sshape.is_done() {
            self.done();
            // OCCT cxx L30: myShape = mySShape.ResultingShape().
            self.my_shape = self.my_sshape.resulting_shape().cloned();
            // OCCT cxx L31: myRight.Clear().
            self.my_right.clear();
        }
    }

    /// OCCT BRepFeat_SplitShape::DirectLeft() (cxx L37-40).
    pub fn direct_left(&self) -> &Vec<Shape> {
        self.my_sshape.direct_left()
    }

    /// OCCT BRepFeat_SplitShape::Left() (cxx L44-47).
    pub fn left(&self) -> &Vec<Shape> {
        self.my_sshape.left()
    }

    /// OCCT BRepFeat_SplitShape::Right() (cxx L51-72).
    pub fn right(&mut self) -> &Vec<Shape> {
        // OCCT cxx L53: if (myRight.IsEmpty()).
        if self.my_right.is_empty() {
            // OCCT cxx L55: NCollection_Map<TopoDS_Shape> aMapOfLeft.
            let mut a_map_of_left = OcctShapeMap::new();
            // OCCT cxx L57-60: add every face of mySShape.Left().
            for an_iterator in self.my_sshape.left() {
                a_map_of_left.add(shape_key(an_iterator), an_iterator.clone());
            }
            // OCCT cxx L61-69: for each FACE of myShape, keep the faces
            // outside aMapOfLeft.
            if let Some(my_shape) = self.my_shape.as_ref() {
                let my_shape = my_shape.clone();
                for a_face in explorer(&my_shape, ShapeType::Face, ShapeType::Shape) {
                    if !a_map_of_left.contains(shape_key(&a_face)) {
                        self.my_right.push(a_face);
                    }
                }
            }
        }
        &self.my_right
    }

    /// OCCT BRepFeat_SplitShape::IsDeleted(F) (cxx L76-82).
    pub fn is_deleted(&mut self, the_f: &Shape) -> bool {
        // OCCT cxx L78: itl over DescendantShapes(F) (the const swindle is
        // the &mut self receiver here).
        // OCCT cxx L81: return (!itl.More()).
        match self.my_sshape.descendant_shapes(the_f) {
            Some(l) => l.is_empty(),
            None => true,
        }
    }

    /// OCCT BRepFeat_SplitShape::Modified(F) (cxx L86-89).
    pub fn modified(&mut self, the_f: &Shape) -> Option<&Vec<Shape>> {
        self.my_sshape.descendant_shapes(the_f)
    }
}
