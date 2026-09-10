// OCCT LocOpe_GluedShape.hxx L32-66 + LocOpe_GluedShape.cxx L36-223 — 1:1
// translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_GluedShape.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_GluedShape.cxx
//
// OCCT inheritance chain (LocOpe_GluedShape.hxx L32):
//   LocOpe_GluedShape : LocOpe_GeneratedShape : Standard_Transient
// Rust has no inheritance: the LocOpe_GeneratedShape base maps to the
// LocOpeGeneratedShape trait (loc_ope_generated_shape.rs); the protected
// members myGEdges/myList of the base are carried as the fields
// my_g_edges/my_list below, and the four virtual overrides are the trait
// impl methods.
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> (myMap,
//    mapdone) — HashMap<(TShape ptr, Location), ()> plus an insertion-order
//    Vec for the iteration (the OCCT bucket iteration order is not
//    reproduced; myMap is iterated in MapEdgeAndVertices L92-96 and the
//    insertion order is the deterministic stand-in, the same reduction as
//    feat::brep_feat_builder::OcctShapeMap).
// 2. NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, ...> (myGShape) —
//    HashMap<(TShape ptr, Location), Option<Shape>>; the Option carries the
//    OCCT null TopoDS_Edge binding (cxx L153 myGShape.Bind(vtx,
//    TopoDS_Edge())).
// 3. TopExp::MapShapesAndAncestors (TopExp.cxx L80-120) is re-hosted below
//    (map_shapes_and_ancestors) over an IndexMap keyed by (TShape ptr,
//    Location) — the same model as bop/algo/section.rs.
// 4. TopExp_Explorer is feat::brep_feat_builder::explorer (same crate).
//
// first consumer: BRepFeat_Form family (3b) — LocOpe_Gluer::Perform (cxx
// L169) and BRepFeat_Form use LocOpe_GluedShape for the glued-face bookkeeping.

use crate::feat::brep_feat_builder::explorer;
use indexmap::IndexMap;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType };
use std::collections::HashMap;

/// OCCT TopAbs::Reverse (TopAbs.hxx) — FORWARD<->REVERSED, INTERNAL/
/// EXTERNAL unchanged.
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// OCCT TopoDS_Shape::Reversed() — a copy with reversed orientation.
fn shape_reversed(s: &Shape) -> Shape {
    let mut c = s.clone();
    c.orientation = top_abs_reverse(c.orientation);
    c
}

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L80-120) re-host
/// (architecture difference #3) — shared with loc_ope_build_shape.rs.
pub(crate) fn map_shapes_and_ancestors(
    the_s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
) {
    for a_anc in explorer(the_s, ta, ShapeType::Shape) {
        for a_exs in explorer(&a_anc, ts, ShapeType::Shape) {
            let key = shape_key(&a_exs);
            let entry = m.entry(key).or_insert((a_exs.clone(), Vec::new()));
            entry.1.push(a_anc.clone());
        }
    }
    for a_ex in explorer(the_s, ts, ta) {
        let key = shape_key(&a_ex);
        m.entry(key).or_insert((a_ex, Vec::new()));
    }
}

/// OCCT LocOpe_GluedShape (LocOpe_GluedShape.hxx L32-66).
pub struct LocOpeGluedShape {
    my_shape: Shape, // OCCT: myShape
    // OCCT: myMap (NCollection_Map) — arch. diff. #1
    my_map: HashMap<(u64, u32), ()>,
    my_map_order: Vec<Shape>,
    // OCCT: myGShape (NCollection_DataMap) — arch. diff. #2
    my_gshape: HashMap<(u64, u32), Option<Shape>>,
    // OCCT LocOpe_GeneratedShape protected members (hxx L50-51)
    my_g_edges: Vec<Shape>, // OCCT: myGEdges
    my_list: Vec<Shape>,    // OCCT: myList
}

impl Default for LocOpeGluedShape {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeGluedShape {
    /// OCCT LocOpe_GluedShape::LocOpe_GluedShape() (cxx L36).
    pub fn new() -> Self {
        LocOpeGluedShape {
            my_shape: Shape::null(),
            my_map: HashMap::new(),
            my_map_order: Vec::new(),
            my_gshape: HashMap::new(),
            my_g_edges: Vec::new(),
            my_list: Vec::new(),
        }
    }

    /// OCCT LocOpe_GluedShape::LocOpe_GluedShape(S) (cxx L40-43).
    pub fn with_shape(the_s: &Shape) -> Self {
        LocOpeGluedShape {
            my_shape: the_s.clone(),
            my_map: HashMap::new(),
            my_map_order: Vec::new(),
            my_gshape: HashMap::new(),
            my_g_edges: Vec::new(),
            my_list: Vec::new(),
        }
    }

    /// OCCT LocOpe_GluedShape::Init(S) (cxx L47-54).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_map.clear();
        self.my_map_order.clear();
        self.my_gshape.clear();
        self.my_list.clear();
        self.my_g_edges.clear();
    }

    /// OCCT LocOpe_GluedShape::GlueOnFace(F) (cxx L58-74).
    pub fn glue_on_face(&mut self, the_f: &Shape) {
        // OCCT cxx L61-68: find F among the faces of myShape.
        let mut found: Option<Shape> = None;
        for cur in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            if shape_key(&cur) == shape_key(the_f) {
                found = Some(cur);
                break;
            }
        }
        let Some(found) = found else {
            // OCCT cxx L69-72: throw Standard_ConstructionError().
            panic!("Standard_ConstructionError");
        };
        // OCCT cxx L73: myMap.Add(exp.Current()) — bonne orientation.
        if self.my_map.insert(shape_key(&found), ()).is_none() {
            self.my_map_order.push(found);
        }
    }

    /// OCCT LocOpe_GluedShape::MapEdgeAndVertices() (cxx L78-179).
    fn map_edge_and_vertices(&mut self) {
        // OCCT cxx L80-83: if (!myGShape.IsEmpty()) return.
        if !self.my_gshape.is_empty() {
            return;
        }

        // Edges et faces generes

        // OCCT cxx L87-89: theMapEF.
        let mut the_map_ef: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        map_shapes_and_ancestors(&self.my_shape, ShapeType::Edge, ShapeType::Face, &mut the_map_ef);

        // OCCT cxx L91-94: mapdone + the myMap iterator.
        let mut mapdone: HashMap<(u64, u32), ()> = HashMap::new();
        // The insertion-ordered keys of myMap stand in for the OCCT map
        // iterator (arch. diff. #1).
        let my_map_snapshot = self.my_map_order.clone();
        for fac in &my_map_snapshot {
            // OCCT cxx L98-99.
            for edg in explorer(fac, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L102-105: mapdone.Contains(edg) -> continue.
                if mapdone.contains_key(&shape_key(&edg)) {
                    continue;
                }
                // OCCT cxx L106-110: "Est-ce un edge de connexite entre les
                // faces collees" — FindFromKey(edg).Extent() != 2 -> throw.
                let Some((_, ancestors)) = the_map_ef.get(&shape_key(&edg)) else {
                    // OCCT FindFromKey on a missing key asserts; the edges of
                    // the glued faces are all under myShape, so the entry
                    // exists — panic keeps the OCCT assert surface.
                    panic!("Standard_NoSuchObject");
                };
                if ancestors.len() != 2 {
                    panic!("Standard_ConstructionError");
                }
                // OCCT cxx L111-117: first ancestor face not in myMap.
                let mut other: Option<Shape> = None;
                for f in ancestors {
                    if !self.my_map.contains_key(&shape_key(f)) {
                        other = Some(f.clone());
                        break;
                    }
                }

                // OCCT cxx L119-125.
                if let Some(other) = other {
                    // myGEdges.Append(edg.Reversed());
                    self.my_g_edges.push(shape_reversed(&edg));
                    // myGShape.Bind(edg, itl.Value());
                    self.my_gshape.insert(shape_key(&edg), Some(other));
                }

                // OCCT cxx L127.
                mapdone.insert(shape_key(&edg), ());
            }
        }

        // OCCT cxx L131-168: vertex bindings.
        let my_g_edges_snapshot = self.my_g_edges.clone();
        for edg in &my_g_edges_snapshot {
            for vtx in explorer(edg, ShapeType::Vertex, ShapeType::Shape) {
                // OCCT cxx L137-140: myGShape.IsBound(vtx) -> continue.
                if self.my_gshape.contains_key(&shape_key(&vtx)) {
                    continue;
                }
                // OCCT cxx L141-146: edges of the generated face of edg.
                let gen_face = self.my_gshape.get(&shape_key(edg)).cloned().flatten();
                let Some(gen_face) = gen_face else {
                    continue;
                };
                let mut done = false;
                for exp2 in explorer(&gen_face, ShapeType::Edge, ShapeType::Shape) {
                    if shape_key(&exp2) == shape_key(edg) {
                        continue;
                    }
                    // OCCT cxx L147-161: vertices of exp2.Current().
                    for exp3 in explorer(&exp2, ShapeType::Vertex, ShapeType::Shape) {
                        if shape_key(&exp3) == shape_key(&vtx) {
                            if self.my_gshape.contains_key(&shape_key(&exp2)) {
                                // OCCT cxx L153: myGShape.Bind(vtx, TopoDS_Edge())
                                // — a null edge binding.
                                self.my_gshape.insert(shape_key(&vtx), None);
                            } else {
                                // OCCT cxx L157: myGShape.Bind(vtx, exp2.Current())
                                self.my_gshape
                                    .insert(shape_key(&vtx), Some(exp2.clone()));
                            }
                            done = true;
                            break;
                        }
                    }
                    // OCCT cxx L162-165: if (exp3.More()) break.
                    if done {
                        break;
                    }
                }
            }
        }

        // liste de faces

        // OCCT cxx L172-178.
        for cur in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            if !self.my_map.contains_key(&shape_key(&cur)) {
                self.my_list.push(cur);
            }
        }
    }

    /// OCCT LocOpe_GluedShape::GeneratingEdges() (cxx L183-190) — the
    /// LocOpe_GeneratedShape override.
    fn lazy_guard(&mut self) {
        // OCCT cxx L185-188 / L196-199 / L207-210 / L218-221: the shared
        // "if (myGShape.IsEmpty()) MapEdgeAndVertices()" guard.
        if self.my_gshape.is_empty() {
            self.map_edge_and_vertices();
        }
    }
}

impl crate::feat::loc_ope_generated_shape::LocOpeGeneratedShape for LocOpeGluedShape {
    /// OCCT LocOpe_GluedShape::GeneratingEdges() (cxx L183-190).
    fn generating_edges(&mut self) -> &Vec<Shape> {
        self.lazy_guard();
        &self.my_g_edges
    }

    /// OCCT LocOpe_GluedShape::Generated(const TopoDS_Vertex& V)
    /// (cxx L194-201).
    fn generated_vertex(&mut self, v: &Shape) -> Shape {
        self.lazy_guard();
        // OCCT cxx L200: return TopoDS::Edge(myGShape(V)) — the bound value
        // (possibly the null edge of cxx L153).
        match self.my_gshape.get(&shape_key(v)) {
            Some(bound) => bound.clone().unwrap_or_else(Shape::null),
            None => {
                // OCCT NCollection_DataMap::operator() asserts IsBound.
                panic!("Standard_NoSuchObject");
            }
        }
    }

    /// OCCT LocOpe_GluedShape::Generated(const TopoDS_Edge& E)
    /// (cxx L205-212).
    fn generated_edge(&mut self, e: &Shape) -> Shape {
        self.lazy_guard();
        // OCCT cxx L211: return TopoDS::Face(myGShape(E)).
        match self.my_gshape.get(&shape_key(e)) {
            Some(bound) => bound.clone().unwrap_or_else(Shape::null),
            None => {
                panic!("Standard_NoSuchObject");
            }
        }
    }

    /// OCCT LocOpe_GluedShape::OrientedFaces() (cxx L216-223).
    fn oriented_faces(&mut self) -> &Vec<Shape> {
        self.lazy_guard();
        &self.my_list
    }
}
