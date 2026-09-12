// OCCT BRepFeat_Builder.cxx L1-829 + BRepFeat_Builder.hxx L1-130 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Builder.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Builder.hxx
//
// OCCT inheritance chain (BRepFeat_Builder.hxx L47):
//   BRepFeat_Builder : BOPAlgo_BOP                 (BOPAlgo_BOP.hxx L63)
//     BOPAlgo_BOP : BOPAlgo_ToolsProvider
//       BOPAlgo_ToolsProvider : BOPAlgo_Builder
//         BOPAlgo_Builder : BOPAlgo_BuilderShape
//           BOPAlgo_BuilderShape : BOPAlgo_Algo
//             BOPAlgo_Algo : BOPAlgo_Options
//
// Rust has no inheritance -> composition + delegation. The BOPAlgo_BOP base
// sub-object is carried as the fields of BRepFeatBuilder listed below; the
// BOPAlgo_BOP *behaviour* is rcad's (PaveFiller, Builder) pair
// (crate::bop::algo::pave_filler::PaveFiller +
//  crate::bop::algo::builder::Builder), the same vehicle that
// bop/brep_algo_api/mod.rs uses for BRepAlgoAPI_BooleanOperation. The field
// mapping, inherited member by inherited member:
//
//   BOPAlgo_Options::myFuzzyValue          -> my_fuzzy_value
//   BOPAlgo_Algo::myReport                 -> my_report (crate::bop::algo::Report)
//   BOPAlgo_BuilderShape::myShape          -> my_shape (rcad: topods::BRep result
//                                             pool; OCCT is a TopoDS_Shape)
//   BOPAlgo_Builder::myArguments           -> my_arguments
//   BOPAlgo_Builder::myImages              -> my_images (OcctDataMapInt keyed by
//                                             (TShape ptr, Location) — OCCT keys
//                                             by TopTools_ShapeMapHasher, which
//                                             ignores orientation)
//   BOPAlgo_Builder::myOrigins             -> my_origins
//   BOPAlgo_Builder::myShapesSD            -> my_shapes_sd
//   BOPAlgo_Builder::myInParts             -> my_in_parts
//   BOPAlgo_Builder::myDS/myPaveFiller     -> my_filler (rcad: the DS lives
//                                             inside an owned PaveFiller; OCCT
//                                             binds myDS from the filler)
//   BOPAlgo_BOP::myOperation               -> my_operation (BooleanOpType =
//                                             BOPAlgo_Operation)
//   BOPAlgo_BOP::myTools                   -> my_tools
//   BOPAlgo_BOP::myDims                    -> my_dims
//   BRepFeat_Builder::myShapes             -> my_shapes (OcctShapeMap)
//   BRepFeat_Builder::myRemoved            -> my_removed (OcctShapeMap)
//   BRepFeat_Builder::myFuse               -> my_fuse
//
// Architecture differences (referenced from the affected functions):
// 1. Virtual overrides of the base class (Clear/Prepare/FillIn3DParts/
//    CheckArgsForOpenSolid) are methods on BRepFeatBuilder; where the base
//    behaviour runs inside rcad's closed Builder pipeline
//    (fill_images_solids -> FillIn3DParts, build_shape ->
//    CheckArgsForOpenSolid) the override cannot be dispatched — the override
//    body is translated here and its base-step coupling is marked at the call
//    site (see perform_result / fill_in_3d_parts).
// 2. OCCT BOPAlgo_BOP::Perform runs the PaveFiller + PerformInternal1 on ONE
//    object whose myDS stays alive across calls. rcad's Builder borrows the DS
//    from the PaveFiller, so it can only exist inside a method scope. The base
//    state (my_shape/my_images/my_origins/my_shapes_sd/my_in_parts/my_report)
//    is parked on BRepFeatBuilder between calls and injected into a
//    method-scoped Builder where base-class methods (FillImagesContainers,
//    FillImagesSolids, FillImagesCompounds, BuildResult) are invoked; the
//    state is parked back afterwards (park_base_state). This is pure Rust
//    scaffolding — the OCCT object has these as plain members.
// 3. OCCT BuildShape (BOPAlgo_BOP.cxx L885-1107) is rcad Builder::build_shape
//    (bop/algo/builder.rs L5522), which is PRIVATE to that module. Under the
//    file-ownership constraint of this stage it cannot be called from here;
//    the two call sites (perform_result L196/L267) are marked. The
//    CheckArgsForOpenSolid override (hxx L122, returns false) affects only
//    BuildShape's open-solid fallback for the same reason.
// 4. Message_ProgressScope / fillPISteps progress plumbing has no rcad
//    Builder counterpart (same reduction as bop/algo/section.rs); the
//    L200-228 progress-weight computation of PerformResult is reduced to a
//    comment at the translated spot.
// 5. NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> is modelled by
//    OcctShapeMap below (OcctDataMapInt keyed by (TShape ptr, Location),
//    orientation ignored) — the same key scheme builder.rs uses for
//    myImages/myShapesSD.
// 6. BOPTools_Set (BOPTools_Set.cxx L26-197) is re-hosted below as
//    BopToolsSet; NCollection_Map<BOPTools_Set> (aMST in RebuildFaces /
//    CheckSolidImages) only uses Add/Contains/Added (no bucket-order
//    iteration), so it is modelled as an insertion-ordered Vec with the
//    OCCT Added semantics (the first inserted element of an equal set is
//    kept and returned by Added).

use crate::bop::algo::builder::{shape_bbox, Builder, BooleanOpType};
use crate::bop::algo::builder_face::BuilderFace;
use crate::bop::algo::occt_map::{OcctDataMapInt, OcctMapInt};
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::algo::Report;
use crate::bop::ds::pave::{Pave, PaveBlock, SharedPB};
use crate::bop::ds::DS;
use crate::bop::tools::algo_tools;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, ShapeType, TShape};
use std::collections::HashMap;

/// OCCT NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> — a shape set
/// keyed by (TShape ptr, Location); orientation is ignored
/// (TopTools_ShapeMapHasher.hxx). Architecture note: OcctDataMapInt requires
/// V: Default (a value container), so the set is carried by a HashMap like
/// builder.rs does for myShapesSD/myOrigins; the OCCT bucket iteration order
/// is not reproduced — the only order-sensitive use (the aMST build in
/// RebuildFaces, which picks the first same-domain representative) iterates
/// myShapes; the representative selection across equal edge-sets follows the
/// HashMap order instead of the OCCT bucket order.
pub(crate) struct OcctShapeMap {
    items: HashMap<(u64, u32), Shape>,
}

impl OcctShapeMap {
    pub fn new() -> Self {
        OcctShapeMap { items: HashMap::new() }
    }
    /// OCCT NCollection_Map::Contains.
    pub fn contains(&self, key: (u64, u32)) -> bool {
        self.items.contains_key(&key)
    }
    /// OCCT NCollection_Map::Add — returns true when newly added.
    pub fn add(&mut self, key: (u64, u32), s: Shape) -> bool {
        if self.items.contains_key(&key) {
            return false;
        }
        self.items.insert(key, s);
        true
    }
    /// OCCT NCollection_Map iterator (key, value) pairs.
    pub fn iter(&self) -> impl Iterator<Item = ((u64, u32), &Shape)> {
        self.items.iter().map(|(k, v)| (*k, v))
    }
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.items.len()
    }
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    /// OCCT NCollection_Map::Clear.
    pub fn clear(&mut self) {
        self.items.clear();
    }
}

/// OCCT NCollection_Map::Add — returns true when the key was newly added.
fn map_add(m: &mut OcctShapeMap, s: &Shape) -> bool {
    m.add((s.ptr_id(), s.location), s.clone())
}

/// Key of a shape for the DS-anchored maps: OCCT holds the arguments by
/// reference, so myImages keys are the DS TShape identities. rcad's DS
/// deep-clones the inputs; an argument-derived shape must have its TShape ptr
/// translated through ds.argument_remap first (DS::argument_remap doc).
fn ds_shape_key(ds: Option<&DS>, s: &Shape) -> (u64, u32) {
    let ptr = ds
        .and_then(|d| d.argument_remap.get(&s.ptr_id()))
        .copied()
        .unwrap_or(s.ptr_id());
    (ptr, s.location)
}

// OCCT TopAbs::ShapeTypeValue order (TopAbs_ShapeEnum) — used by the explorer
// for the "ty > toFind" pruning (TopExp_Explorer.cxx L90-95).
fn topabs_value(t: ShapeType) -> i32 {
    match t {
        ShapeType::Compound => 0,
        ShapeType::CompSolid => 1,
        ShapeType::Solid => 2,
        ShapeType::Shell => 3,
        ShapeType::Face => 4,
        ShapeType::Wire => 5,
        ShapeType::Edge => 6,
        ShapeType::Vertex => 7,
        ShapeType::Shape => 8,
    }
}

/// OCCT TopoDS_Iterator (cumOri=true): the stored sub-shapes with
/// TopAbs::Compose(parent.Orientation(), child.Orientation()) applied
/// (TopoDS_Iterator.cxx L72-80). The child order is the rcad container field
/// order used by builder.rs push_shape_recursive (same model as
/// bop/algo/section.rs sub_shapes).
pub(crate) fn sub_shapes(sh: &Shape) -> Vec<Shape> {
    let composed = |c: &Shape| {
        let mut c = c.clone();
        c.orientation = sh.orientation.compose(c.orientation);
        c
    };
    match sh.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => vec![composed(&ed.first), composed(&ed.last)],
        TShape::Wire(wd) => wd.edges.iter().map(composed).collect(),
        TShape::Face(fd) => {
            let mut out =
                Vec::with_capacity(1 + fd.inner_wires.len() + fd.internal_vertices.len());
            out.push(composed(&fd.outer_wire));
            out.extend(fd.inner_wires.iter().map(composed));
            out.extend(fd.internal_vertices.iter().map(composed));
            out
        }
        TShape::Shell(sd) => sd.faces.iter().map(composed).collect(),
        TShape::Solid(sd) => {
            let mut out = Vec::new();
            out.extend(sd.shells.iter().map(composed));
            out.extend(sd.internal_vertices.iter().map(composed));
            out.extend(sd.internal_edges.iter().map(composed));
            out
        }
        TShape::CompSolid(cs) => cs.iter().map(composed).collect(),
        TShape::Compound(cd) => cd.iter().map(composed).collect(),
    }
}

/// OCCT TopExp_Explorer (TopExp_Explorer.cxx L68-171): yields the sub-shapes
/// of `s` of type `to_find`, descending into coarser container types, never
/// entering types of `to_avoid` (ShapeType::Shape is the "avoid nothing"
/// sentinel). Re-hosted from the bop/algo/section.rs model.
pub(crate) fn explorer(s: &Shape, to_find: ShapeType, to_avoid: ShapeType) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    // OCCT Init L84-87: if (toFind == TopAbs_SHAPE) { hasMore = false; }
    if to_find == ShapeType::Shape {
        return out;
    }
    // OCCT Init L90-95: if (ty > toFind) { hasMore = false; }
    if topabs_value(s.shape_type()) > topabs_value(to_find) {
        return out;
    }
    explore_rec(s, to_find, to_avoid, &mut out);
    out
}

/// The iterative walk of TopExp_Explorer::Next (TopExp_Explorer.cxx L110-171)
/// in recursive form (same model as bop/algo/section.rs explore_rec).
fn explore_rec(sh: &Shape, to_find: ShapeType, to_avoid: ShapeType, out: &mut Vec<Shape>) {
    let ty = sh.shape_type();
    if ty == to_find {
        out.push(sh.clone());
        return;
    }
    if topabs_value(ty) > topabs_value(to_find) {
        return;
    }
    // shouldAvoid (TopExp_Explorer.cxx L26-29): toAvoid != TopAbs_SHAPE.
    if to_avoid != ShapeType::Shape && ty == to_avoid {
        return;
    }
    for c in sub_shapes(sh) {
        explore_rec(&c, to_find, to_avoid, out);
    }
}

/// OCCT TopExp::MapShapes (TopExp.cxx L35-45): explore S for T and add every
/// found shape into M.
pub(crate) fn map_shapes(s: &Shape, t: ShapeType, m: &mut OcctShapeMap) {
    for x in explorer(s, t, ShapeType::Shape) {
        map_add(m, &x);
    }
}

/// OCCT TopExp::MapShapes over an rcad result pool: the pool is the
/// TopoDS_Compound root; its top-level (unreferenced) TShapes are the
/// compound children. Mirrors bop/brep_algo_api::brep_top_shapes.
fn pool_top_shapes(brep: &rcad_kernel::topods::BRep, t: ShapeType) -> Vec<Shape> {
    let mut referenced = vec![false; brep.tshapes.len()];
    fn mark(sr: &Shape, referenced: &mut Vec<bool>) {
        let i = sr.index;
        if i >= referenced.len() || referenced[i] {
            return;
        }
        referenced[i] = true;
        match &*sr.data {
            TShape::Solid(sd) => {
                for x in &sd.shells {
                    mark(x, referenced);
                }
                for x in &sd.internal_vertices {
                    mark(x, referenced);
                }
                for x in &sd.internal_edges {
                    mark(x, referenced);
                }
            }
            TShape::Shell(sd) => {
                for x in &sd.faces {
                    mark(x, referenced);
                }
            }
            TShape::Face(fd) => {
                mark(&fd.outer_wire, referenced);
                for w in &fd.inner_wires {
                    mark(w, referenced);
                }
                for v in &fd.internal_vertices {
                    mark(v, referenced);
                }
            }
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    mark(e, referenced);
                }
            }
            TShape::Edge(ed) => {
                mark(&ed.first, referenced);
                mark(&ed.last, referenced);
            }
            TShape::CompSolid(cs) => {
                for x in cs {
                    mark(x, referenced);
                }
            }
            TShape::Compound(cd) => {
                for x in cd {
                    mark(x, referenced);
                }
            }
            _ => {}
        }
    }
    for ts in &brep.tshapes {
        match ts.as_ref() {
            TShape::Solid(sd) => {
                for sr in &sd.shells {
                    mark(sr, &mut referenced);
                }
                for sr in &sd.internal_vertices {
                    mark(sr, &mut referenced);
                }
                for sr in &sd.internal_edges {
                    mark(sr, &mut referenced);
                }
            }
            TShape::Shell(sd) => {
                for sr in &sd.faces {
                    mark(sr, &mut referenced);
                }
            }
            TShape::Face(fd) => {
                mark(&fd.outer_wire, &mut referenced);
                for w in &fd.inner_wires {
                    mark(w, &mut referenced);
                }
                for v in &fd.internal_vertices {
                    mark(v, &mut referenced);
                }
            }
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    mark(e, &mut referenced);
                }
            }
            TShape::Edge(ed) => {
                mark(&ed.first, &mut referenced);
                mark(&ed.last, &mut referenced);
            }
            TShape::CompSolid(shapes) => {
                for sr in shapes {
                    mark(sr, &mut referenced);
                }
            }
            TShape::Compound(shapes) => {
                for sr in shapes {
                    mark(sr, &mut referenced);
                }
            }
            _ => {}
        }
    }
    brep.tshapes
        .iter()
        .enumerate()
        .filter(|(i, ts)| !referenced[*i] && tshape_is_type(ts, t))
        .map(|(i, ts)| {
            Shape::from_parts(ts.clone(), i, 0, rcad_kernel::topods::Orientation::Forward)
        })
        .collect()
}

/// Type test of a pooled TShape (the TShape variant -> ShapeType map).
fn tshape_is_type(ts: &TShape, t: ShapeType) -> bool {
    let st = match ts {
        TShape::Vertex(_) => ShapeType::Vertex,
        TShape::Edge(_) => ShapeType::Edge,
        TShape::Wire(_) => ShapeType::Wire,
        TShape::Face(_) => ShapeType::Face,
        TShape::Shell(_) => ShapeType::Shell,
        TShape::Solid(_) => ShapeType::Solid,
        TShape::CompSolid(_) => ShapeType::CompSolid,
        TShape::Compound(_) => ShapeType::Compound,
    };
    st == t
}

/// Faces of an rcad result pool (OCCT TopExp_Explorer(myShape, TopAbs_FACE)).
pub(crate) fn pool_faces(brep: &rcad_kernel::topods::BRep) -> Vec<Shape> {
    // The pool is flat: every Face TShape is a face of the result (the
    // compound-root explorer reaches them all through the containment tree).
    let mut out = Vec::new();
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if ts.shape_type() == ShapeType::Face {
            out.push(Shape::from_parts(
                ts.clone(),
                i,
                0,
                rcad_kernel::topods::Orientation::Forward,
            ));
        }
    }
    out
}

/// OCCT BRep_Tool::IsClosed(aE, aF) (BRep_Tool.cxx L1760-1790): an edge is
/// closed on a face when it carries a CurveOnClosedSurface (two pcurves) on
/// that face. Same model as bop/algo/shell_splitter.rs edge_closed_on_face.
fn edge_closed_on_face(a_e: &Shape, a_f: &Shape) -> bool {
    let f_key = (a_f.ptr_id(), a_f.location);
    a_e.as_edge()
        .map(|ed| {
            ed.representations.iter().any(|r| {
                matches!(
                    r,
                    topods::CurveRepresentation::CurveOnClosedSurface { face, .. }
                        if *face == f_key
                )
            })
        })
        .unwrap_or(false)
}

/// OCCT BOPTools_Set::NormalizedIds (BOPTools_Set.cxx L183-197).
fn normalized_ids(a_id: usize, a_div: usize) -> usize {
    if a_div == 0 {
        return a_id;
    }
    let a_max = usize::MAX;
    let a_tresh = a_max / a_div;
    if a_id > a_tresh {
        a_id % a_tresh
    } else {
        a_id
    }
}

/// OCCT BOPTools_Set (BOPTools_Set.cxx L26-197 + BOPTools_Set.hxx) — a
/// signature of the sub-shapes of a given type of one shape, used for the
/// same-domain detection among the kept parts.
pub(crate) struct BopToolsSet {
    my_shape: Option<Shape>, // OCCT: myShape (TopoDS_Shape)
    my_nb_shapes: i32,       // OCCT: myNbShapes
    my_sum: usize,           // OCCT: mySum
    my_upper: usize,         // OCCT: myUpper = 432123
    my_shapes: Vec<Shape>,   // OCCT: myShapes (NCollection_List)
}

impl BopToolsSet {
    /// OCCT BOPTools_Set::BOPTools_Set (Set.cxx L26-33).
    pub fn new() -> Self {
        BopToolsSet {
            my_shape: None,
            my_nb_shapes: 0,
            my_sum: 0,
            my_upper: 432123,
            my_shapes: Vec::new(),
        }
    }

    /// OCCT BOPTools_Set::Clear (Set.cxx L47-52).
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.my_nb_shapes = 0;
        self.my_sum = 0;
        self.my_shapes.clear();
    }

    /// OCCT BOPTools_Set::Shape (Set.cxx L66-69).
    pub fn shape(&self) -> Option<&Shape> {
        self.my_shape.as_ref()
    }

    /// OCCT BOPTools_Set::Add (Set.cxx L125-176).
    pub fn add(&mut self, the_s: &Shape, the_type: ShapeType) {
        self.my_shape = Some(the_s.clone());
        self.my_shapes.clear();
        self.my_nb_shapes = 0;
        self.my_sum = 0;
        //
        for a_sx in explorer(the_s, the_type, ShapeType::Shape) {
            if the_type == ShapeType::Edge {
                let is_degenerated = a_sx
                    .as_edge()
                    .map(|ed| ed.degenerated)
                    .unwrap_or(false);
                if is_degenerated {
                    continue;
                }
            }
            let a_or = a_sx.orientation;
            if a_or == rcad_kernel::topods::Orientation::Internal {
                let mut a_sy = a_sx.clone();
                a_sy.orientation = rcad_kernel::topods::Orientation::Forward;
                self.my_shapes.push(a_sy);
                let mut a_sy = a_sx.clone();
                a_sy.orientation = rcad_kernel::topods::Orientation::Reversed;
                self.my_shapes.push(a_sy);
            } else {
                self.my_shapes.push(a_sx);
            }
        }
        //
        self.my_nb_shapes = self.my_shapes.len() as i32;
        if self.my_nb_shapes == 0 {
            return;
        }
        //
        for a_sx in &self.my_shapes {
            // OCCT L171: aId = TopTools_ShapeMapHasher{}(aSx) % myUpper + 1 —
            // the hasher is the TShape handle hash (ptr), Location only in
            // equality.
            let a_id = a_sx.ptr_id() as usize % self.my_upper + 1;
            let a_id_n = normalized_ids(a_id, self.my_nb_shapes as usize);
            self.my_sum += a_id_n;
        }
    }

    /// OCCT BOPTools_Set::IsEqual (Set.cxx L73-98).
    pub fn is_equal(&self, the_other: &BopToolsSet) -> bool {
        if the_other.my_nb_shapes != self.my_nb_shapes {
            return false;
        }
        //
        let mut a_m1: OcctShapeMap = OcctShapeMap::new();
        //
        for s in &self.my_shapes {
            map_add(&mut a_m1, s);
        }
        //
        for s in &the_other.my_shapes {
            if !a_m1.contains((s.ptr_id(), s.location)) {
                return false;
            }
        }
        //
        true
    }
}

/// OCCT NCollection_Map<BOPTools_Set> — only Add/Contains/Added are used (no
/// bucket-order iteration), so an insertion-ordered Vec with the OCCT Added
/// semantics (the first inserted element is kept) is a faithful model.
pub(crate) struct BopToolsSetMap {
    items: Vec<BopToolsSet>,
}

impl BopToolsSetMap {
    pub fn new() -> Self {
        BopToolsSetMap { items: Vec::new() }
    }
    /// OCCT NCollection_Map::Contains.
    pub fn contains(&self, the_set: &BopToolsSet) -> bool {
        self.items.iter().any(|s| s.is_equal(the_set))
    }
    /// OCCT NCollection_Map::Add — returns true when newly added.
    pub fn add(&mut self, the_set: BopToolsSet) -> bool {
        if self.contains(&the_set) {
            return false;
        }
        self.items.push(the_set);
        true
    }
    /// OCCT NCollection_Map::Added — the element as it is in the map (the
    /// first inserted among equals).
    pub fn added(&self, the_set: &BopToolsSet) -> Option<&BopToolsSet> {
        self.items.iter().find(|s| s.is_equal(the_set))
    }
}

/// OCCT BRepFeat_Builder — provides a basic tool to implement features
/// topological operations (BRepFeat_Builder.hxx L32-46).
pub struct BRepFeatBuilder {
    // --- BOPAlgo_Algo (inherited) ---
    pub(crate) my_report: Report,      // BOPAlgo_Algo::myReport
    pub(crate) my_fuzzy_value: f64,    // BOPAlgo_Options::myFuzzyValue
    // --- BOPAlgo_BuilderShape (inherited) ---
    pub(crate) my_shape: Option<topods::BRep>, // BOPAlgo_BuilderShape::myShape
    // --- BOPAlgo_BOP (inherited) ---
    pub(crate) my_operation: BooleanOpType, // BOPAlgo_BOP::myOperation
    pub(crate) my_tools: Vec<Shape>,   // BOPAlgo_BOP::myTools
    pub(crate) my_dims: [i32; 2],      // BOPAlgo_BOP::myDims
    // --- BOPAlgo_Builder (inherited) ---
    pub(crate) my_arguments: Vec<Shape>, // BOPAlgo_Builder::myArguments
    pub(crate) my_images: OcctDataMapInt<(u64, u32), Vec<Shape>>, // BOPAlgo_Builder::myImages
    pub(crate) my_origins: HashMap<(u64, u32), Vec<Shape>>,       // BOPAlgo_Builder::myOrigins
    pub(crate) my_shapes_sd: HashMap<(u64, u32), Shape>,          // BOPAlgo_Builder::myShapesSD
    pub(crate) my_in_parts: HashMap<(u64, u32), Vec<Shape>>,      // BOPAlgo_Builder::myInParts
    pub(crate) my_filler: Option<PaveFiller>, // BOPAlgo_Builder::myPaveFiller/myDS carrier
    // --- BRepFeat_Builder (BRepFeat_Builder.hxx L124-126) ---
    pub(crate) my_shapes: OcctShapeMap, // L124: myShapes
    pub(crate) my_removed: OcctShapeMap, // L125: myRemoved
    pub(crate) my_fuse: i32,            // L126: myFuse
}

impl BRepFeatBuilder {
    /// OCCT BRepFeat_Builder::BRepFeat_Builder() (Builder.cxx L44-48) — calls
    /// Clear().
    pub fn new() -> Self {
        let mut b = BRepFeatBuilder {
            my_report: Report::new(),
            my_fuzzy_value: rcad_kernel::precision::CONFUSION,
            my_shape: None,
            my_operation: BooleanOpType::Unknown,
            my_tools: Vec::new(),
            my_dims: [-1, -1],
            my_arguments: Vec::new(),
            my_images: OcctDataMapInt::new(),
            my_origins: HashMap::new(),
            my_shapes_sd: HashMap::new(),
            my_in_parts: HashMap::new(),
            my_filler: None,
            my_shapes: OcctShapeMap::new(),
            my_removed: OcctShapeMap::new(),
            my_fuse: 0,
        };
        b.clear();
        b
    }

    /// OCCT BRepFeat_Builder::Clear (Builder.cxx L56-61): myShapes.Clear();
    /// myRemoved.Clear(); BOPAlgo_BOP::Clear().
    ///
    /// BOPAlgo_BOP::Clear (BOPAlgo_BOP.cxx L81-89) resets myOperation/myDims
    /// and calls BOPAlgo_ToolsProvider::Clear() -> BOPAlgo_Builder::Clear()
    /// -> ... which resets every inherited container. rcad: reset the fields.
    pub fn clear(&mut self) {
        self.my_shapes.clear();
        self.my_removed.clear();
        // BOPAlgo_BOP::Clear
        self.my_operation = BooleanOpType::Unknown;
        self.my_dims = [-1, -1];
        // BOPAlgo_ToolsProvider::Clear -> BOPAlgo_Builder::Clear ->
        // BOPAlgo_BuilderShape::Clear -> BOPAlgo_Algo::Clear
        self.my_arguments.clear();
        self.my_tools.clear();
        self.my_images.clear();
        self.my_origins.clear();
        self.my_shapes_sd.clear();
        self.my_in_parts.clear();
        self.my_shape = None;
        self.my_filler = None;
        self.my_report.clear();
        self.my_fuzzy_value = rcad_kernel::precision::CONFUSION;
    }

    /// OCCT BRepFeat_Builder::Init(theShape) (Builder.cxx L65-70) —
    /// initializes the object of local boolean operation.
    pub fn init(&mut self, the_shape: &Shape) {
        self.clear();
        //
        self.add_argument(the_shape.clone());
    }

    /// OCCT BRepFeat_Builder::Init(theShape, theTool) (Builder.cxx L74-80) —
    /// initializes the arguments of local boolean operation (the second Init
    /// overload; Rust has no overloading).
    pub fn init_with_tool(&mut self, the_shape: &Shape, the_tool: &Shape) {
        self.clear();
        //
        self.add_argument(the_shape.clone());
        self.add_tool(the_tool.clone());
    }

    /// OCCT BOPAlgo_Builder::AddArgument — appends to myArguments.
    pub fn add_argument(&mut self, the_shape: Shape) {
        self.my_arguments.push(the_shape);
    }

    /// OCCT BOPAlgo_BOP::AddTool — appends to myTools.
    pub fn add_tool(&mut self, the_shape: Shape) {
        self.my_tools.push(the_shape);
    }

    /// OCCT BRepFeat_Builder::SetOperation(theFuse) (Builder.cxx L84-88).
    /// If theFuse = 0 than the operation is CUT, otherwise FUSE.
    pub fn set_operation(&mut self, the_fuse: i32) {
        self.my_fuse = the_fuse;
        self.my_operation = if self.my_fuse != 0 {
            BooleanOpType::Union
        } else {
            BooleanOpType::Cut
        };
    }

    /// OCCT BRepFeat_Builder::SetOperation(theFuse, theFlag) (Builder.cxx
    /// L92-103) — the second SetOperation overload; Rust has no overloading.
    /// If theFlag = TRUE it means that no selection of parts of the tool is
    /// needed, t.e. no second part. In that case if theFuse = 0 than
    /// operation is COMMON, otherwise CUT21. If theFlag = FALSE
    /// SetOperation(theFuse) function is called.
    pub fn set_operation_with_flag(&mut self, the_fuse: i32, the_flag: bool) {
        self.my_fuse = the_fuse;
        if !the_flag {
            self.my_operation = if self.my_fuse != 0 {
                BooleanOpType::Union
            } else {
                BooleanOpType::Cut
            };
        } else {
            self.my_operation = if self.my_fuse != 0 {
                BooleanOpType::Cut21
            } else {
                BooleanOpType::Intersection
            };
        }
    }

    /// OCCT BRepAlgo_Algo::HasErrors — reads the parked report.
    pub fn has_errors(&self) -> bool {
        self.my_report.has_errors()
    }

    /// OCCT BOPAlgo_BuilderShape::Shape — the result of the operation. rcad:
    /// the result pool (architecture difference #2 of the header).
    pub fn shape_brep(&mut self) -> Option<topods::BRep> {
        self.my_shape.take()
    }

    /// OCCT BRepFeat_Builder::Perform — the base BOPAlgo_BOP::Perform
    /// (BOPAlgo_BOP.cxx L364-412), driven on rcad's (PaveFiller, Builder)
    /// vehicle (architecture difference #2 of the header).
    pub fn perform_bop(&mut self) {
        // OCCT L372: GetReport()->Clear();
        self.my_report.clear();
        // OCCT L375-379: myEntryPoint == 1 -> delete myPaveFiller — rcad builds
        // a fresh PaveFiller per Perform (a_filler below), same effect.
        // OCCT L381-399: aLS = myArguments + myTools.
        let mut a_ls = self.my_arguments.clone();
        a_ls.extend(self.my_tools.iter().cloned());
        // OCCT L401-408: pPF = new BOPAlgo_PaveFiller; SetArguments(aLS);
        // SetRunParallel; SetFuzzyValue; SetNonDestructive; SetGlue; SetUseOBB.
        let mut p_pf = PaveFiller::new();
        p_pf.set_arguments(a_ls);
        p_pf.set_fuzzy_value(self.my_fuzzy_value);
        // myRunParallel / myNonDestructive / myGlue / myUseOBB — rcad defaults
        // (false / false / GlueOff / false) match the OCCT defaults; the
        // Builder level knobs are not carried by this facade yet.
        // OCCT L404: Message_ProgressScope aPS(theRange, "Performing Boolean
        // operation", 10) — progress plumbing dropped (difference #4).
        // OCCT L409: pPF->Perform(aPS.Next(9));
        let a_prog = rcad_kernel::message::NoopProgress;
        let a_ps = rcad_kernel::message::ProgressScope::new(&a_prog, "Performing Boolean operation", 10);
        p_pf.perform(&a_ps);
        // OCCT L411-412: myEntryPoint = 1; PerformInternal(*pPF, ...) ->
        // PerformInternal1 — the rcad equivalent is the closed Builder
        // pipeline (Prepare + FillImages* + BuildResult + BuildShape).
        let mut a_builder = Builder::new(p_pf.ds(), self.my_operation, self.my_fuzzy_value);
        // rcad's DS deep-clones the inputs; the builder must carry the
        // DS-cloned shapes (same wiring as bop/brep_algo_api::run_build).
        a_builder.my_arguments = p_pf.ds().arguments.clone();
        let n_tools = self.my_tools.len();
        let n_objs = a_builder.my_arguments.len().saturating_sub(n_tools);
        a_builder.my_tools = a_builder.my_arguments[n_objs..].to_vec();
        let a_res = a_builder.build_with_history_topods();
        // Park the base state back (architecture difference #2).
        self.my_report = std::mem::replace(&mut a_builder.my_report, Report::new());
        self.my_images = std::mem::replace(&mut a_builder.my_images, OcctDataMapInt::new());
        self.my_origins = std::mem::take(&mut a_builder.my_origins);
        self.my_shapes_sd = std::mem::take(&mut a_builder.my_shapes_sd);
        self.my_in_parts = std::mem::take(&mut a_builder.my_in_parts);
        match a_res {
            Ok((brep, _)) => self.my_shape = Some(brep),
            Err(_) => self.my_shape = None,
        }
        self.my_filler = Some(p_pf);
    }

    /// OCCT BRepFeat_Builder::PartsOfTool (Builder.cxx L107-118) — collects
    /// parts of the tool. OCCT fills aLT in place; rcad returns the list.
    pub fn parts_of_tool(&self) -> Vec<Shape> {
        // OCCT L111: aLT.Clear();
        // OCCT L112-117: aExp.Init(myShape, TopAbs_SOLID) — the SOLIDs of the
        // result. rcad: the parked result pool's top-level SOLID TShapes
        // (architecture difference: OCCT myShape is a compound shape, rcad a
        // flat pool).
        let mut a_lt: Vec<Shape> = Vec::new();
        if let Some(brep) = self.my_shape.as_ref() {
            a_lt.extend(pool_top_shapes(brep, ShapeType::Solid));
        }
        a_lt
    }

    /// OCCT BRepFeat_Builder::KeepParts (Builder.cxx L122-131) — initializes
    /// parts of the tool for second step of algorithm.
    pub fn keep_parts(&mut self, the_im: &[Shape]) {
        for a_t_im in the_im {
            self.keep_part(a_t_im);
        }
    }

    /// OCCT BRepFeat_Builder::KeepPart (Builder.cxx L135-141) — adds shape
    /// theS and all its sub-shapes into myShapes map.
    pub fn keep_part(&mut self, the_part: &Shape) {
        // OCCT L138-140: TopExp::MapShapes(thePart, myShapes) (default type
        // TopAbs_SHAPE).
        map_shapes(the_part, ShapeType::Shape, &mut self.my_shapes);
    }

    /// OCCT BRepFeat_Builder::Prepare (Builder.cxx L145-155) — prepares
    /// builder of local operation (the base Prepare override).
    pub fn prepare(&mut self) {
        // OCCT L147: GetReport()->Clear();
        self.my_report.clear();
        //
        // OCCT L149-152: BRep_Builder aBB; TopoDS_Compound aC;
        // aBB.MakeCompound(aC); myShape = aC;
        // rcad: the result pool; the locations mirror of the DS (the
        // builder.rs Builder::prepare scaffolding) is replicated here.
        let mut a_c = topods::BRep::new();
        if let Some(a_filler) = self.my_filler.as_ref() {
            a_c.locations = a_filler.ds().locations.clone();
        }
        self.my_shape = Some(a_c);
        //
        // OCCT L154: FillRemoved();
        self.fill_removed();
    }

    /// OCCT BRepFeat_Builder::FillRemoved() (Builder.cxx L159-187) — collects
    /// the removed parts of the tool into myRemoved map.
    pub fn fill_removed(&mut self) {
        // OCCT L163-164: aArgs0 = myArguments.First(); aArgs1 =
        // myTools.First(); (OCCT assumes both lists are non-empty; rcad
        // guards the empty case.)
        let Some(a_args0) = self.my_arguments.first().cloned() else {
            return;
        };
        let Some(a_args1) = self.my_tools.first().cloned() else {
            return;
        };
        //
        // OCCT L166-171: for each SOLID of aArgs0: myImages.UnBind(aS).
        for a_s in explorer(&a_args0, ShapeType::Solid, ShapeType::Shape) {
            self.my_images.remove(ds_shape_key(self.ds_ref(), &a_s));
        }
        //
        // OCCT L173-176: if (!myImages.IsBound(aArgs1)) return;
        let a_key1 = ds_shape_key(self.ds_ref(), &a_args1);
        if !self.my_images.contains(a_key1) {
            return;
        }
        //
        // OCCT L178-186: for each image of aArgs1: FillRemoved(aS, myRemoved).
        // The list is cloned: the recursion only touches myShapes/myRemoved,
        // never myImages (same data flow as OCCT).
        let a_ls = self.my_images.get(a_key1).cloned().unwrap_or_default();
        for a_s in &a_ls {
            Self::fill_removed_in(&self.my_shapes, a_s, &mut self.my_removed);
        }
    }

    /// OCCT BRepFeat_Builder::PerformResult (Builder.cxx L191-268) — main
    /// function to build the result of the local operation required.
    pub fn perform_result(&mut self) {
        // OCCT L193: myOperation = myFuse ? BOPAlgo_FUSE : BOPAlgo_CUT;
        self.my_operation = if self.my_fuse != 0 {
            BooleanOpType::Union
        } else {
            BooleanOpType::Cut
        };
        // OCCT L194-198: if (myShapes.IsEmpty()) { BuildShape(theRange);
        // return; } — with no kept parts the result stays the base Perform
        // result; BuildShape (difference #3) is pending exposure, and the
        // parked my_shape already reflects that same base result.
        if self.my_shapes.is_empty() {
            return;
        }
        //
        // OCCT L200-228: Message_ProgressScope aPS("BRepFeat_Builder", 100)
        // and the aSteps weight computation from getNbShapes() — progress
        // plumbing dropped (difference #4).
        //
        // OCCT L230: Prepare();
        self.prepare();
        //
        // The DS carrier is taken out of self so rebuild_faces can mutate it
        // (architecture difference #2; restored at every exit below).
        let mut a_filler = match std::mem::take(&mut self.my_filler) {
            Some(f) => f,
            None => {
                // OCCT: myDS exists after Perform; rcad gap — nothing to
                // rebuild without the DS.
                return;
            }
        };
        //
        // OCCT L232: RebuildFaces();
        self.rebuild_faces(&mut a_filler);
        //
        // The base steps run on a method-scoped Builder over the DS
        // (architecture difference #2); the parked state is injected here.
        let mut a_builder = Builder::new(a_filler.ds(), self.my_operation, self.my_fuzzy_value);
        a_builder.my_arguments = a_filler.ds().arguments.clone();
        let n_tools = self.my_tools.len();
        let n_objs = a_builder.my_arguments.len().saturating_sub(n_tools);
        a_builder.my_tools = a_builder.my_arguments[n_objs..].to_vec();
        a_builder.my_report = std::mem::replace(&mut self.my_report, Report::new());
        a_builder.my_shape = self.my_shape.take();
        a_builder.my_images =
            std::mem::replace(&mut self.my_images, OcctDataMapInt::new());
        a_builder.my_origins = std::mem::take(&mut self.my_origins);
        a_builder.my_shapes_sd = std::mem::take(&mut self.my_shapes_sd);
        a_builder.my_in_parts = std::mem::take(&mut self.my_in_parts);
        a_builder.my_dims = self.my_dims;
        a_builder.shape_remap.clear();
        //
        'base_steps: {
            // OCCT L235-239: FillImagesContainers(TopAbs_SHELL); if
            // (HasErrors()) return;
            a_builder.fill_images_containers(ShapeType::Shell);
            if a_builder.has_errors() {
                break 'base_steps;
            }
            //
            // OCCT L241-245: FillImagesSolids(); if (HasErrors()) return;
            // (The BRepFeat FillIn3DParts override — difference #1 — is
            // applied by fill_in_3d_parts right after.)
            a_builder.fill_images_solids();
            if a_builder.has_errors() {
                break 'base_steps;
            }
            self.fill_in_3d_parts(&mut a_builder);
            //
            // OCCT L247: CheckSolidImages();
            self.check_solid_images(&mut a_builder);
            //
            // OCCT L249-253: BuildResult(TopAbs_SOLID); if (HasErrors())
            // return;
            a_builder.build_result(ShapeType::Solid);
            if a_builder.has_errors() {
                break 'base_steps;
            }
            //
            // OCCT L255-259: FillImagesCompounds(); if (HasErrors()) return;
            a_builder.fill_images_compounds();
            if a_builder.has_errors() {
                break 'base_steps;
            }
            //
            // OCCT L261-265: BuildResult(TopAbs_COMPOUND); if (HasErrors())
            // return;
            a_builder.build_result(ShapeType::Compound);
            if a_builder.has_errors() {
                break 'base_steps;
            }
            //
            // OCCT L267: BuildShape(aPS.Next(aBSPart)); — architecture gap
            // (difference #3): Builder::build_shape is private to
            // bop/algo/builder.rs; the final result-shape composition is
            // pending its exposure. The my_shape parked below carries the
            // state built up to BuildResult(COMPOUND).
        }
        //
        // Park the base state back and restore the DS carrier.
        self.park_base_state(&mut a_builder);
        drop(a_builder);
        self.my_filler = Some(a_filler);
    }

    /// Park the base state of a method-scoped Builder back onto self
    /// (architecture difference #2 scaffolding).
    fn park_base_state(&mut self, a_builder: &mut Builder<'_>) {
        self.my_report = std::mem::replace(&mut a_builder.my_report, Report::new());
        self.my_shape = a_builder.my_shape.take();
        self.my_images =
            std::mem::replace(&mut a_builder.my_images, OcctDataMapInt::new());
        self.my_origins = std::mem::take(&mut a_builder.my_origins);
        self.my_shapes_sd = std::mem::take(&mut a_builder.my_shapes_sd);
        self.my_in_parts = std::mem::take(&mut a_builder.my_in_parts);
    }

    /// Read-only access to the DS carrier (for the argument-remap of the
    /// FillRemoved lookups; the DS stays parked in my_filler).
    fn ds_ref(&self) -> Option<&DS> {
        self.my_filler.as_ref().map(|f| f.ds())
    }

    /// OCCT BRepFeat_Builder::RebuildFaces (Builder.cxx L272-558) — rebuilds
    /// faces in accordance with the kept parts of the tool.
    ///
    /// The DS is mutated (RebuildEdge appends split edges), so it is reached
    /// through the parked PaveFiller passed explicitly (architecture
    /// difference #2).
    pub fn rebuild_faces(&mut self, a_filler: &mut PaveFiller) {
        // OCCT L281: NCollection_Map<TopoDS_Shape> aME, aMESplit;
        // aME is cleared per edge (L364) and aMESplit fences the per-edge
        // RebuildEdge call; both live at the function scope in OCCT.
        // OCCT L285: NCollection_Map<BOPTools_Set> aMST;
        let mut a_mst = BopToolsSetMap::new();
        let mut a_me_split: OcctShapeMap = OcctShapeMap::new();
        //
        // OCCT L288-298: for each FACE in myShapes: BOPTools_Set of its
        // edges into aMST.
        for (_k, a_s) in self.my_shapes.iter() {
            if a_s.shape_type() == ShapeType::Face {
                let mut a_st = BopToolsSet::new();
                a_st.add(a_s, ShapeType::Edge);
                a_mst.add(a_st);
            }
        }
        //
        // OCCT L300: aNbS = myDS->NbSourceShapes();
        let a_nb_s = a_filler.ds().nb_source_shapes();
        for i in 0..a_nb_s {
            let a_si = a_filler.ds().shape_info(i).clone();
            //
            // OCCT L305: iRank = myDS->Rank(i);
            let i_rank = a_filler.ds().rank(i);
            if i_rank == 1 {
                // OCCT L308-326: the tool faces — drop the non-kept images.
                let a_s = a_si.shape().clone();
                let a_key = (a_s.ptr_id(), a_s.location);
                if self.my_images.contains(a_key) {
                    if let Some(a_l_im) = self.my_images.get_mut(a_key) {
                        // OCCT L313-323: remove while iterating (Remove at
                        // the iterator, no Next on removal).
                        let mut j = 0usize;
                        while j < a_l_im.len() {
                            let a_s_im = a_l_im[j].clone();
                            if !self.my_shapes.contains((a_s_im.ptr_id(), a_s_im.location)) {
                                a_l_im.remove(j);
                                continue;
                            }
                            j += 1;
                        }
                    }
                }
                continue;
            }
            //
            // OCCT L328-331: if (aSI.ShapeType() != TopAbs_FACE) continue;
            if a_si.shape_type() != ShapeType::Face {
                continue;
            }
            //
            // OCCT L333: const BOPDS_FaceInfo& aFI = myDS->FaceInfo(i);
            // rcad guards the absent FaceInfo (rcad face_info(i) panics on a
            // missing entry; OCCT only reaches faces with interference data
            // through the myImages bound check below).
            if !a_filler.ds().has_face_info(i) {
                continue;
            }
            let a_fi = a_filler.ds().face_info(i).clone();
            // OCCT L334: const TopoDS_Face& aF = aSI.Shape();
            let a_f = a_si.shape().clone();
            //
            // OCCT L336-339: if (!myImages.IsBound(aF)) continue;
            let a_f_key = (a_f.ptr_id(), a_f.location);
            if !self.my_images.contains(a_f_key) {
                continue;
            }
            //
            // OCCT L341-343: anOriF = aF.Orientation(); aFF = aF;
            // aFF.Orientation(TopAbs_FORWARD);
            let an_ori_f = a_f.orientation;
            let mut a_ff = a_f.clone();
            a_ff.orientation = rcad_kernel::topods::Orientation::Forward;
            //
            // OCCT L345-346: aMPBIn = aFI.PaveBlocksIn(); aMPBSc =
            // aFI.PaveBlocksSc(); (rcad FaceInfo keeps the pave block
            // pointers as u64 ptr ids, insertion-ordered like the OCCT
            // IndexedMap.)
            let a_m_pb_in: Vec<u64> = a_fi.pave_blocks_in.iter().copied().collect();
            let a_m_pb_sc: Vec<u64> = a_fi.pave_blocks_sc.iter().copied().collect();
            //
            // OCCT L348: aLE.Clear();
            let mut a_le: Vec<Shape> = Vec::new();
            //
            // OCCT L350-476: bounding edges.
            for a_e in explorer(&a_ff, ShapeType::Edge, ShapeType::Shape) {
                // OCCT L354-357: anOriE; bIsDegenerated; bIsClosed.
                let an_ori_e = a_e.orientation;
                let b_is_degenerated = a_e
                    .as_edge()
                    .map(|ed| ed.degenerated)
                    .unwrap_or(false);
                let b_is_closed = edge_closed_on_face(&a_e, &a_f);
                // OCCT L358: if (myImages.IsBound(aE))
                let a_e_key = (a_e.ptr_id(), a_e.location);
                if self.my_images.contains(a_e_key) {
                    // OCCT L360: aLEIm = myImages.ChangeFind(aE);
                    let mut a_me: OcctShapeMap = OcctShapeMap::new(); // L364: aME.Clear()
                    let mut a_le_im_new: Vec<Shape> = Vec::new(); // L365
                    //
                    // OCCT L367-401: classify the images by the kept parts.
                    let a_le_im = self.my_images.get(a_e_key).cloned().unwrap_or_default();
                    let mut b_rem = false; // L362
                    let mut b_im = false; // L363
                    for a_s in &a_le_im {
                        // OCCT L372-389: bVInShapes.
                        let mut b_v_in_shapes = false;
                        if self
                            .my_shapes
                            .contains((a_s.ptr_id(), a_s.location))
                        {
                            b_v_in_shapes = true;
                        } else {
                            for a_v in explorer(a_s, ShapeType::Vertex, ShapeType::Shape) {
                                if self
                                    .my_shapes
                                    .contains((a_v.ptr_id(), a_v.location))
                                {
                                    b_v_in_shapes = true;
                                    break;
                                }
                            }
                        }
                        //
                        if b_v_in_shapes {
                            b_im = true;
                            a_le_im_new.push(a_s.clone());
                        } else {
                            b_rem = true;
                            map_add(&mut a_me, a_s);
                        }
                    }
                    //
                    // OCCT L403-407: if (!bIm) { aLE.Append(aE); continue; }
                    if !b_im {
                        a_le.push(a_e.clone());
                        continue;
                    }
                    //
                    // OCCT L409-426: if (bRem && bIm)
                    if b_rem && b_im {
                        if a_le_im.len() == 2 {
                            a_le.push(a_e.clone());
                            continue;
                        }
                        if map_add(&mut a_me_split, &a_e) {
                            self.rebuild_edge(a_filler, &a_e, &a_ff, &a_me, &mut a_le_im_new);
                            // OCCT L419: aLEIm.Assign(aLEImNew);
                            if let Some(dst) = self.my_images.get_mut(a_e_key) {
                                *dst = a_le_im_new.clone();
                            }
                            if a_le_im_new.len() == 1 {
                                a_le.push(a_e.clone());
                                continue;
                            }
                        }
                    }
                    //
                    // OCCT L428-470: append the oriented split edges.
                    for a_s in self.my_images.get(a_e_key).cloned().unwrap_or_default() {
                        let mut a_sp = a_s; // OCCT L431: aSp = aItIm.Value()
                        //
                        if b_is_degenerated {
                            a_sp.orientation = an_ori_e;
                            a_le.push(a_sp);
                            continue;
                        }
                        //
                        if an_ori_e == rcad_kernel::topods::Orientation::Internal {
                            a_sp.orientation = rcad_kernel::topods::Orientation::Forward;
                            a_le.push(a_sp.clone());
                            a_sp.orientation = rcad_kernel::topods::Orientation::Reversed;
                            a_le.push(a_sp);
                            continue;
                        }
                        //
                        if b_is_closed {
                            if !edge_closed_on_face(&a_sp, &a_ff) {
                                // OCCT L453: BOPTools_AlgoTools3D::
                                // DoSplitSEAMOnFace(aSp, aFF);
                                // Architecture gap: hosted as a private method
                                // of Builder (bop/algo/builder.rs L3038);
                                // exposure pending (difference #3 class).
                            }
                            //
                            a_sp.orientation = rcad_kernel::topods::Orientation::Forward;
                            a_le.push(a_sp.clone());
                            a_sp.orientation = rcad_kernel::topods::Orientation::Reversed;
                            a_le.push(a_sp);
                            continue;
                        } // if (bIsClosed){
                        //
                        a_sp.orientation = an_ori_e;
                        // OCCT L464: bToReverse =
                        // BOPTools_AlgoTools::IsSplitToReverse(aSp, aE,
                        // myContext);
                        let (b_to_reverse, _err) =
                            algo_tools::is_split_to_reverse_edge(&a_sp, &a_e);
                        if b_to_reverse {
                            // OCCT L467: aSp.Reverse();
                            a_sp.orientation = match a_sp.orientation {
                                rcad_kernel::topods::Orientation::Forward => {
                                    rcad_kernel::topods::Orientation::Reversed
                                }
                                rcad_kernel::topods::Orientation::Reversed => {
                                    rcad_kernel::topods::Orientation::Forward
                                }
                                other => other,
                            };
                        }
                        a_le.push(a_sp);
                    }
                } else {
                    // OCCT L472-475: else { aLE.Append(aE); }
                    a_le.push(a_e);
                }
            }
            //
            // OCCT L478-480: aNbPBIn / aNbPBSc.
            let a_nb_pb_in = a_m_pb_in.len();
            let a_nb_pb_sc = a_m_pb_sc.len();
            //
            // OCCT L482-497: in edges.
            for j in 0..a_nb_pb_in {
                let a_pb = match a_filler.ds().pb_from_ptr(a_m_pb_in[j]) {
                    Some(pb) => pb,
                    None => continue,
                };
                let n_sp = a_pb.read().edge(); // OCCT L486: nSp = aPB->Edge();
                let a_sp = a_filler.ds().shape(n_sp).clone(); // OCCT L487
                if self
                    .my_removed
                    .contains((a_sp.ptr_id(), a_sp.location))
                {
                    continue;
                }
                //
                let mut a_sp_f = a_sp.clone();
                a_sp_f.orientation = rcad_kernel::topods::Orientation::Forward;
                a_le.push(a_sp_f);
                let mut a_sp_r = a_sp;
                a_sp_r.orientation = rcad_kernel::topods::Orientation::Reversed;
                a_le.push(a_sp_r);
            }
            // OCCT L498-513: section edges.
            for j in 0..a_nb_pb_sc {
                let a_pb = match a_filler.ds().pb_from_ptr(a_m_pb_sc[j]) {
                    Some(pb) => pb,
                    None => continue,
                };
                let n_sp = a_pb.read().edge();
                let a_sp = a_filler.ds().shape(n_sp).clone();
                if self
                    .my_removed
                    .contains((a_sp.ptr_id(), a_sp.location))
                {
                    continue;
                }
                //
                let mut a_sp_f = a_sp.clone();
                a_sp_f.orientation = rcad_kernel::topods::Orientation::Forward;
                a_le.push(a_sp_f);
                let mut a_sp_r = a_sp;
                a_sp_r.orientation = rcad_kernel::topods::Orientation::Reversed;
                a_le.push(a_sp_r);
            }
            //
            // OCCT L515-520: build new faces — BOPAlgo_BuilderFace aBF;
            // aBF.SetFace(aFF); aBF.SetShapes(aLE); aBF.Perform();
            let mut a_bf = BuilderFace::new(a_filler.ds());
            a_bf.my_face = Some(a_ff.clone());
            a_bf.my_edges = a_le;
            a_bf.perform();
            //
            // OCCT L522-523: aLFIm = myImages.ChangeFind(aF); aLFIm.Clear();
            let mut b_lf_im_empty = false;
            if let Some(a_lf_im) = self.my_images.get_mut(a_f_key) {
                a_lf_im.clear();
                //
                // OCCT L525-551: for each area of aBF.Areas(): the
                // BOPTools_Set same-domain bookkeeping.
                for a_fr in &a_bf.my_areas {
                    let mut a_st = BopToolsSet::new();
                    a_st.add(a_fr, ShapeType::Edge);
                    let b_flag_sd = a_mst.contains(&a_st);
                    //
                    let a_stx = a_mst.added(&a_st);
                    if let Some(a_stx) = a_stx {
                        if let Some(a_sx_shape) = a_stx.shape() {
                            let mut a_sx = a_sx_shape.clone();
                            a_sx.orientation = an_ori_f;
                            a_lf_im.push(a_sx.clone());
                            //
                            // OCCT L540-545: myOrigins bookkeeping.
                            let a_sx_key = (a_sx.ptr_id(), a_sx.location);
                            self.my_origins
                                .entry(a_sx_key)
                                .or_default()
                                .push(a_f.clone());
                            //
                            // OCCT L547-550: if (bFlagSD)
                            // myShapesSD.Bind(aFR, aSx);
                            if b_flag_sd {
                                self.my_shapes_sd
                                    .insert((a_fr.ptr_id(), a_fr.location), a_sx);
                            }
                        }
                    }
                }
                //
                // OCCT L553-556: if (aLFIm.Extent() == 0) myImages.UnBind(aF);
                b_lf_im_empty = a_lf_im.is_empty();
            }
            if b_lf_im_empty {
                self.my_images.remove(a_f_key);
            }
        }
    }

    /// OCCT BRepFeat_Builder::RebuildEdge (Builder.cxx L562-723) — rebuilds
    /// edges in accordance with the kept parts of the tool. OCCT mutates
    /// myDS in place; rcad reaches it through the parked PaveFiller
    /// (architecture difference #2).
    fn rebuild_edge(
        &mut self,
        a_filler: &mut PaveFiller,
        the_e: &Shape,
        the_f: &Shape,
        a_me: &OcctShapeMap,
        a_l_im: &mut Vec<Shape>,
    ) {
        // OCCT L574-583 declarations -> locals below.
        // OCCT L585: aSI.SetShapeType(TopAbs_EDGE) — folded into
        // DS::push_edge_inherit, which appends the split-edge ShapeInfo with
        // the Edge type directly (architecture note at L713-718).
        //
        // 1. collect origin vertices to aMV map (OCCT L587-602).
        let n_e_i = a_filler.ds().index(the_e);
        let n_e = if n_e_i < 0 { 0usize } else { n_e_i as usize };
        // OCCT L577: NCollection_Map<int> aMI, aMAdd, aMV, aMVOr; (aMI/aMAdd
        // are declared here and created at their use sites below.)
        let mut a_mv = OcctMapInt::new();
        let mut a_mv_or = OcctMapInt::new();
        {
            let a_ls: Vec<usize> = a_filler.ds().shape_info(n_e).sub_shapes().to_vec();
            for n_v in a_ls {
                let mut n_vx = n_v;
                let mut n_vsd = 0usize;
                if a_filler.ds().has_shape_sd(n_v, &mut n_vsd) {
                    n_vx = n_vsd;
                }
                a_mv.add(n_vx);
                a_mv_or.add(n_vx);
            }
        }
        //
        // 2. collect vertices that should be removed to aMI map (OCCT
        // L604-625).
        // OCCT L605: aPBNew = new BOPDS_PaveBlock (default ctor: unset edge).
        let mut a_pb_new = PaveBlock::new(
            usize::MAX,
            Pave::new(usize::MAX, 0.0),
            Pave::new(usize::MAX, 0.0),
        );
        // OCCT L606: aLPExt = aPBNew->ChangeExtPaves();
        // OCCT L607: aLPB = myDS->ChangePaveBlocks(nE);
        let a_lpb_snap: Vec<SharedPB> = a_filler.ds().pave_blocks(n_e).to_vec();
        //
        let mut a_mi: OcctMapInt = OcctMapInt::new();
        // OCCT L582: aMPB (NCollection_Map<handle<BOPDS_PaveBlock>>).
        let mut a_mpb: Vec<SharedPB> = Vec::new();
        let mut a_mpb_ptrs: std::collections::HashSet<u64> = std::collections::HashSet::new();
        for a_pb in &a_lpb_snap {
            // OCCT L611-613: nE1 = aPB->Edge(); aE1 = myDS->Shape(nE1);
            let a_e1 = {
                let pbr = a_pb.read();
                let n_e1 = pbr.edge();
                a_filler.ds().shape(n_e1).clone()
            };
            //
            if a_me.contains((a_e1.ptr_id(), a_e1.location)) {
                let (n_v1, n_v2) = a_pb.read().indices();
                a_mi.add(n_v1);
                a_mi.add(n_v2);
            } else {
                // OCCT: aMPB.Add(aPB) — the NCollection_Map<handle> identity
                // is the handle (Arc) pointer.
                let pb_ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
                if a_mpb_ptrs.insert(pb_ptr) {
                    a_mpb.push(a_pb.clone());
                }
            }
        }
        // 3. collect vertices that split the source shape (OCCT L626-640).
        for a_pb in &a_lpb_snap {
            let (n_v1, n_v2) = a_pb.read().indices();
            //
            if !a_mi.contains(n_v1) {
                a_mv.add(n_v1);
            }
            if !a_mi.contains(n_v2) {
                a_mv.add(n_v2);
            }
        }
        // 4. collect ext paves (OCCT L641-662).
        {
            let a_lp_ext = a_pb_new.change_ext_paves();
            let mut a_m_add: OcctMapInt = OcctMapInt::new();
            for a_pb in &a_lpb_snap {
                let (n_v1, n_v2) = a_pb.read().indices();
                //
                if a_mv.contains(n_v1) {
                    if a_m_add.add(n_v1) || a_mv_or.contains(n_v1) {
                        let p1 = a_pb.read().pave1().clone();
                        a_lp_ext.push(p1);
                    }
                }
                //
                if a_mv.contains(n_v2) {
                    if a_m_add.add(n_v2) || a_mv_or.contains(n_v2) {
                        let p2 = a_pb.read().pave2().clone();
                        a_lp_ext.push(p2);
                    }
                }
            }
        }
        //
        // OCCT L664-665: aE = theE; aE.Orientation(TopAbs_FORWARD);
        let mut a_e = the_e.clone();
        a_e.orientation = rcad_kernel::topods::Orientation::Forward;
        //
        // OCCT L667: aLIm.Clear();
        a_l_im.clear();
        //
        // 5. split edge by new set of vertices (OCCT L669-672).
        {
            let a_lpb = a_filler.ds_mut().change_pave_blocks(n_e);
            a_lpb.clear();
            a_pb_new.set_original_edge(n_e);
            a_pb_new.update(a_lpb, false);
        }
        let a_lpb: Vec<SharedPB> = a_filler.ds().pave_blocks(n_e).to_vec();
        for a_pb in &a_lpb {
            let (n_v1, a_t1, n_v2, a_t2) = {
                let pbr = a_pb.read();
                let (n_v1, n_v2) = pbr.indices();
                let (a_t1, a_t2) = pbr.range();
                (n_v1, a_t1, n_v2, a_t2)
            };
            //
            // OCCT L683-698: check if it is the old pave block.
            let mut b_old = false;
            for a_pb1 in &a_mpb {
                let (n_v11, n_v21) = a_pb1.read().indices();
                let (a_t11, a_t21) = a_pb1.read().range();
                if n_v1 == n_v11 && n_v2 == n_v21 && a_t1 == a_t11 && a_t2 == a_t21 {
                    let n_edge = a_pb1.read().edge();
                    let a_e_im = a_filler.ds().shape(n_edge).clone();
                    a_l_im.push(a_e_im);
                    b_old = true;
                    break;
                }
            }
            if b_old {
                continue;
            }
            //
            // OCCT L704-708: aV1 = myDS->Shape(nV1) FORWARD; aV2 =
            // myDS->Shape(nV2) REVERSED; BOPTools_AlgoTools::MakeSplitEdge
            // (aE, aV1, aT1, aV2, aT2, aSp) (BOPTools_AlgoTools_2.cxx
            // L138-183). rcad: the vertex orientations and the split-edge
            // production are folded into DS::push_edge_inherit (the
            // MakeSplitEdge host documented there).
            let a_curve_opt = a_filler.ds().edge_curve(n_e);
            let Some(a_curve) = a_curve_opt else {
                // OCCT asserts a curve on the source edge; rcad skips the
                // iteration when the DS has no curve (same fallback as
                // pave_filler.rs MakeSplitEdges).
                continue;
            };
            let (n_sp, a_sp) = {
                let ds = a_filler.ds_mut();
                let n_sp = ds.push_edge_inherit(a_curve, [a_t1, a_t2], n_v1, n_v2, Some(n_e));
                let a_sp = ds.shape(n_sp).clone();
                (n_sp, a_sp)
            };
            // OCCT L711: BOPTools_AlgoTools2D::BuildPCurveForEdgeOnFace(aSp,
            // theF, myContext) — rcad: the BOPTools_AlgoTools2D placeholder
            // (pcurve creation is handled by the PaveFiller MakePCurves step).
            let n_f_i = a_filler.ds().index(the_f);
            let n_f = if n_f_i < 0 { 0usize } else { n_f_i as usize };
            algo_tools::make_pcurve(n_sp, n_f, usize::MAX, a_filler.ds_mut());
            // OCCT L713-718: aSI.SetShape(aSp); Bnd_Box& aBox =
            // aSI.ChangeBox(); BRepBndLib::Add(aSp, aBox); nSp =
            // myDS->Append(aSI); — the ShapeInfo construction and append are
            // folded into push_edge_inherit; the box rebuild matches
            // BRepBndLib::Add + the PaveFiller box convention.
            a_filler.rebuild_edge_box(n_sp);
            //
            // OCCT L720: aPB->SetEdge(nSp);
            a_pb.write().set_edge(n_sp);
            // OCCT L721: aLIm.Append(aSp);
            a_l_im.push(a_sp);
        }
    }

    /// OCCT BRepFeat_Builder::CheckSolidImages (Builder.cxx L727-778) —
    /// collects the images of the object, that contains in the images of the
    /// tool. Runs on the injected base state (architecture difference #2).
    fn check_solid_images(&mut self, base: &mut Builder<'_>) {
        // OCCT L729: NCollection_Map<BOPTools_Set> aMST;
        let mut a_mst = BopToolsSetMap::new();
        // OCCT L730: NCollection_List<TopoDS_Shape> aLSImNew;
        let mut a_ls_im_new: Vec<Shape> = Vec::new();
        // OCCT L731: NCollection_Map<TopoDS_Shape> aMS — declared and unused
        // in OCCT.
        //
        // OCCT L736-737: aArgs0 = myArguments.First(); aArgs1 =
        // myTools.First();
        let Some(a_args0) = self.my_arguments.first().cloned() else {
            return;
        };
        let Some(a_args1) = self.my_tools.first().cloned() else {
            return;
        };
        //
        // OCCT L739-748: for each image of aArgs1: BOPTools_Set of its faces
        // into aMST.
        let a_key1 = ds_shape_key(Some(base.ds), &a_args1);
        let a_ls_im = base.my_images.get(a_key1).cloned().unwrap_or_default();
        for a_sol_im in &a_ls_im {
            let mut a_st = BopToolsSet::new();
            a_st.add(a_sol_im, ShapeType::Face);
            a_mst.add(a_st);
        }
        //
        // OCCT L750-777: for each SOLID of aArgs0 with images: keep only the
        // first same-domain representative.
        for a_solid in explorer(&a_args0, ShapeType::Solid, ShapeType::Shape) {
            let a_solid_key = ds_shape_key(Some(base.ds), &a_solid);
            if base.my_images.contains(a_solid_key) {
                let a_ls_im_sol = base.my_images.get(a_solid_key).cloned().unwrap_or_default();
                for a_sol_im in &a_ls_im_sol {
                    let mut a_st = BopToolsSet::new();
                    a_st.add(a_sol_im, ShapeType::Face);
                    let b_flag_sd = a_mst.contains(&a_st);
                    //
                    let a_stx = a_mst.added(&a_st);
                    if let Some(a_stx) = a_stx {
                        if let Some(a_sx) = a_stx.shape() {
                            a_ls_im_new.push(a_sx.clone());
                            //
                            if b_flag_sd {
                                base.my_shapes_sd
                                    .insert((a_sol_im.ptr_id(), a_sol_im.location), a_sx.clone());
                            }
                        }
                    }
                }
                // OCCT L775: aLSImSol.Assign(aLSImNew);
                if let Some(dst) = base.my_images.get_mut(a_solid_key) {
                    *dst = a_ls_im_new.clone();
                }
            }
        }
    }

    /// OCCT BRepFeat_Builder::FillRemoved(S, M) (Builder.cxx L782-797) —
    /// adds the shape S and its sub-shapes into myRemoved map (the recursive
    /// FillRemoved overload; Rust has no overloading). The myShapes/myRemoved
    /// members are passed explicitly: the whole-self receiver borrow would
    /// conflict with the myRemoved out-parameter (Rust scaffolding).
    fn fill_removed_in(my_shapes: &OcctShapeMap, s: &Shape, m: &mut OcctShapeMap) {
        if my_shapes.contains((s.ptr_id(), s.location)) {
            return;
        }
        //
        map_add(m, s);
        // OCCT L791: TopoDS_Iterator It(S);
        for child in sub_shapes(s) {
            Self::fill_removed_in(my_shapes, &child, m);
        }
    }

    /// OCCT BRepFeat_Builder::FillIn3DParts (Builder.cxx L801-828) — the base
    /// FillIn3DParts override; clears the report, runs the base fill, then
    /// removes the IN parts of the solids built from the removed faces.
    ///
    /// Architecture difference #1: the base call
    /// (BOPAlgo_Builder::FillIn3DParts) runs inside rcad's closed
    /// Builder::fill_images_solids; the override is applied right after that
    /// step in perform_result. Runs on the injected base state.
    fn fill_in_3d_parts(&mut self, base: &mut Builder<'_>) {
        // OCCT L805: GetReport()->Clear(); and L807:
        // BOPAlgo_Builder::FillIn3DParts(theDraftSolids, theRange); — both
        // covered by Builder::fill_images_solids (the report is cleared per
        // step in rcad).
        //
        // OCCT L809-827: clear the IN parts of the solids from the removed
        // faces.
        for (_k, a_list) in base.my_in_parts.iter_mut() {
            let mut j = 0usize;
            while j < a_list.len() {
                let v = &a_list[j];
                if self.my_removed.contains((v.ptr_id(), v.location)) {
                    a_list.remove(j);
                } else {
                    j += 1;
                }
            }
        }
    }

    /// Bounding box of a shape against a locations table
    /// (BRepBndLib::Add equivalent — bop/algo/builder.rs shape_bbox), used by
    /// the feature algorithms without a DS of their own.
    pub(crate) fn shape_box(s: &Shape, locations: &[glam::DAffine3]) -> Option<(glam::DVec3, glam::DVec3)> {
        shape_bbox(s, locations)
    }
}
