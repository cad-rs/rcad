// OCCT HLRTopoBRep_Data (TKHLR/HLRTopoBRep/HLRTopoBRep_Data.hxx L39-147
// + .cxx L28-303 + .lxx L26-107) — the results of the OutLine and IsoLine
// processes: the per-edge split lists (mySplE), the per-face interference
// lists (myData), the new/old shape map (myOldS), the outline/internal
// vertex sets (myOutV / myIntV) and the per-edge vertex lists
// (myEdgesVertices) together with the (edge, vertex) iteration state.
//
// Map-key architecture note: OCCT keys every map with TopoDS_Shape under
// TopTools_ShapeMapHasher (TShape-pointer identity; IsSame ignores the
// location).  rcad encodes NCollection_DataMap/Map as insertion-ordered
// `Vec<(Shape, V)>` / `Vec<Shape>` keyed by `Shape::ptr_id()` (the Arc
// TShape pointer) — the same identity semantics and the same Vec mapping
// as geomalgo/hatch/elements.rs; the DataMap iteration order (unspecified
// in OCCT) becomes the deterministic insertion order.  The OCCT raw
// pointer `myVList` (into the map-owned list) becomes the map position
// `my_v_list`; the list iterators become position indices with
// `usize::MAX` as the default/exhausted state (More() == false).

use rcad_kernel::topods::Shape;

use super::face_data::FaceData;
use super::v_data::VData;

/// OCCT NCollection_DataMap::IsBound — the TopTools_ShapeMapHasher identity
/// (TShape pointer only).
fn map_is_bound<V>(map: &[(Shape, V)], s: &Shape) -> bool {
    let id = s.ptr_id();
    map.iter().any(|(k, _)| k.ptr_id() == id)
}

/// The map position of the shape — OCCT NCollection_DataMap::operator()(K)
/// (Find) raises Standard_NoSuchObject when unbound; the expect mirrors it.
fn map_pos<V>(map: &[(Shape, V)], s: &Shape) -> usize {
    let id = s.ptr_id();
    map.iter()
        .position(|(k, _)| k.ptr_id() == id)
        .expect("Standard_NoSuchObject: HLRTopoBRep_Data map key not bound")
}

/// OCCT HLRTopoBRep_Data.
pub struct Data {
    /// OCCT myOldS (hxx L139) — NCollection_DataMap<TopoDS_Shape,
    /// TopoDS_Shape, TopTools_ShapeMapHasher>.
    my_old_s: Vec<(Shape, Shape)>,
    /// OCCT mySplE (hxx L140) — NCollection_DataMap<TopoDS_Shape,
    /// NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>.
    my_spl_e: Vec<(Shape, Vec<Shape>)>,
    /// OCCT myData (hxx L141) — NCollection_DataMap<TopoDS_Shape,
    /// HLRTopoBRep_FaceData, TopTools_ShapeMapHasher>.
    my_data: Vec<(Shape, FaceData)>,
    /// OCCT myOutV (hxx L142) — NCollection_Map<TopoDS_Shape,
    /// TopTools_ShapeMapHasher>.
    my_out_v: Vec<Shape>,
    /// OCCT myIntV (hxx L143) — NCollection_Map<TopoDS_Shape,
    /// TopTools_ShapeMapHasher>.
    my_int_v: Vec<Shape>,
    /// OCCT myEdgesVertices (hxx L144-145) —
    /// NCollection_DataMap<TopoDS_Shape, NCollection_List<HLRTopoBRep_VData>,
    /// TopTools_ShapeMapHasher>.
    my_edges_vertices: Vec<(Shape, Vec<VData>)>,
    /// OCCT myEIterator (hxx L146) — the myEdgesVertices iterator position
    /// (usize::MAX = the default / exhausted iterator, More() == false).
    my_e_iterator: usize,
    /// OCCT myVIterator (hxx L147) — the vertex-list iterator position
    /// (usize::MAX = the default / exhausted iterator, More() == false).
    my_v_iterator: usize,
    /// OCCT myVList (hxx L148) — the NCollection_List<HLRTopoBRep_VData>*
    /// into myEdgesVertices; rcad keeps the map position (usize::MAX =
    /// unset).
    my_v_list: usize,
}

impl Data {
    /// OCCT HLRTopoBRep_Data() (cxx L28-32).
    pub fn new() -> Self {
        let mut d = Data {
            my_old_s: Vec::new(),
            my_spl_e: Vec::new(),
            my_data: Vec::new(),
            my_out_v: Vec::new(),
            my_int_v: Vec::new(),
            my_edges_vertices: Vec::new(),
            my_e_iterator: usize::MAX,
            my_v_iterator: usize::MAX,
            my_v_list: usize::MAX,
        };
        d.clear();
        d
    }

    /// OCCT Clear() (cxx L35-45) — clear of all the maps.
    pub fn clear(&mut self) {
        self.my_old_s.clear();
        self.my_spl_e.clear();
        self.my_data.clear();
        self.my_out_v.clear();
        self.my_int_v.clear();
        self.my_edges_vertices.clear();
    }

    /// OCCT Clean() (cxx L47-49) — clear of all the data not needed during
    /// and after the hiding process.
    pub fn clean(&mut self) {}

    /// OCCT EdgeHasSplE(const TopoDS_Edge& E) (cxx L51-58) — returns True if
    /// the Edge is split.
    pub fn edge_has_spl_e(&self, e: &Shape) -> bool {
        if !map_is_bound(&self.my_spl_e, e) {
            return false;
        }
        !self.edge_spl_e(e).is_empty()
    }

    /// OCCT FaceHasIntL(const TopoDS_Face& F) (cxx L62-69) — returns True if
    /// the Face has internal outline.
    pub fn face_has_int_l(&self, f: &Shape) -> bool {
        if !map_is_bound(&self.my_data, f) {
            return false;
        }
        let pos = map_pos(&self.my_data, f);
        !self.my_data[pos].1.face_int_l().is_empty()
    }

    /// OCCT FaceHasOutL(const TopoDS_Face& F) (cxx L73-80) — returns True if
    /// the Face has outlines on restriction.
    pub fn face_has_out_l(&self, f: &Shape) -> bool {
        if !map_is_bound(&self.my_data, f) {
            return false;
        }
        let pos = map_pos(&self.my_data, f);
        !self.my_data[pos].1.face_out_l().is_empty()
    }

    /// OCCT FaceHasIsoL(const TopoDS_Face& F) (cxx L84-91) — returns True if
    /// the Face has isolines.
    pub fn face_has_iso_l(&self, f: &Shape) -> bool {
        if !map_is_bound(&self.my_data, f) {
            return false;
        }
        let pos = map_pos(&self.my_data, f);
        !self.my_data[pos].1.face_iso_l().is_empty()
    }

    /// OCCT IsSplEEdgeEdge(const TopoDS_Edge& E1, const TopoDS_Edge& E2)
    /// (cxx L95-112).
    pub fn is_spl_e_edge_edge(&self, e1: &Shape, e2: &Shape) -> bool {
        let mut found = false;
        if self.edge_has_spl_e(e1) {
            // NCollection_List<TopoDS_Shape>::Iterator itS.
            let list = self.edge_spl_e(e1);
            let mut it_s = 0;
            while it_s < list.len() && !found {
                found = list[it_s].is_same(e2);
                it_s += 1;
            }
        } else {
            found = e1.is_same(e2);
        }
        found
    }

    /// OCCT IsIntLFaceEdge(const TopoDS_Face& F, const TopoDS_Edge& E)
    /// (cxx L116-131).
    pub fn is_int_l_face_edge(&self, f: &Shape, e: &Shape) -> bool {
        let mut found = false;
        if self.face_has_int_l(f) {
            let list = self.face_int_l(f);
            // NCollection_List<TopoDS_Shape>::Iterator itE.
            let mut it_e = 0;
            while it_e < list.len() && !found {
                // TopoDS::Edge(itE.Value()) — the list stores edges; rcad's
                // Shape is the untyped handle.
                found = self.is_spl_e_edge_edge(&list[it_e], e);
                it_e += 1;
            }
        }
        found
    }

    /// OCCT IsOutLFaceEdge(const TopoDS_Face& F, const TopoDS_Edge& E)
    /// (cxx L133-148).
    pub fn is_out_l_face_edge(&self, f: &Shape, e: &Shape) -> bool {
        let mut found = false;
        if self.face_has_out_l(f) {
            let list = self.face_out_l(f);
            // NCollection_List<TopoDS_Shape>::Iterator itE.
            let mut it_e = 0;
            while it_e < list.len() && !found {
                // TopoDS::Edge(itE.Value()).
                found = self.is_spl_e_edge_edge(&list[it_e], e);
                it_e += 1;
            }
        }
        found
    }

    /// OCCT IsIsoLFaceEdge(const TopoDS_Face& F, const TopoDS_Edge& E)
    /// (cxx L150-165).
    pub fn is_iso_l_face_edge(&self, f: &Shape, e: &Shape) -> bool {
        let mut found = false;
        if self.face_has_iso_l(f) {
            let list = self.face_iso_l(f);
            // NCollection_List<TopoDS_Shape>::Iterator itE.
            let mut it_e = 0;
            while it_e < list.len() && !found {
                // TopoDS::Edge(itE.Value()).
                found = self.is_spl_e_edge_edge(&list[it_e], e);
                it_e += 1;
            }
        }
        found
    }

    /// OCCT NewSOldS(const TopoDS_Shape& NewS) (cxx L167-176).
    pub fn new_s_old_s(&self, new_s: &Shape) -> Shape {
        if map_is_bound(&self.my_old_s, new_s) {
            let pos = map_pos(&self.my_old_s, new_s);
            self.my_old_s[pos].1.clone()
        } else {
            new_s.clone()
        }
    }

    /// OCCT EdgeSplE(const TopoDS_Edge& E) (lxx L26-31) — returns the list
    /// of the edges.
    pub fn edge_spl_e(&self, e: &Shape) -> &Vec<Shape> {
        let pos = map_pos(&self.my_spl_e, e);
        &self.my_spl_e[pos].1
    }

    /// OCCT FaceIntL(const TopoDS_Face& F) (lxx L33-38) — returns the list
    /// of the internal OutLines.
    pub fn face_int_l(&self, f: &Shape) -> &Vec<Shape> {
        let pos = map_pos(&self.my_data, f);
        self.my_data[pos].1.face_int_l()
    }

    /// OCCT FaceOutL(const TopoDS_Face& F) (lxx L40-45) — returns the list
    /// of the OutLines on restriction.
    pub fn face_out_l(&self, f: &Shape) -> &Vec<Shape> {
        let pos = map_pos(&self.my_data, f);
        self.my_data[pos].1.face_out_l()
    }

    /// OCCT FaceIsoL(const TopoDS_Face& F) (lxx L47-52) — returns the list
    /// of the IsoLines.
    pub fn face_iso_l(&self, f: &Shape) -> &Vec<Shape> {
        let pos = map_pos(&self.my_data, f);
        self.my_data[pos].1.face_iso_l()
    }

    /// OCCT IsOutV(const TopoDS_Vertex& V) (lxx L56-59) — returns True if V
    /// is an outline vertex on a restriction.
    pub fn is_out_v(&self, v: &Shape) -> bool {
        let id = v.ptr_id();
        self.my_out_v.iter().any(|s| s.ptr_id() == id)
    }

    /// OCCT IsIntV(const TopoDS_Vertex& V) (lxx L63-66) — returns True if V
    /// is an internal outline vertex.
    pub fn is_int_v(&self, v: &Shape) -> bool {
        let id = v.ptr_id();
        self.my_int_v.iter().any(|s| s.ptr_id() == id)
    }

    /// OCCT AddOldS(const TopoDS_Shape& NewS, const TopoDS_Shape& OldS)
    /// (cxx L181-187).
    pub fn add_old_s(&mut self, new_s: &Shape, old_s: &Shape) {
        if !map_is_bound(&self.my_old_s, new_s) {
            self.my_old_s.push((new_s.clone(), old_s.clone()));
        }
    }

    /// OCCT AddSplE(const TopoDS_Edge& E) (cxx L191-200).
    pub fn add_spl_e(&mut self, e: &Shape) -> &mut Vec<Shape> {
        if !map_is_bound(&self.my_spl_e, e) {
            self.my_spl_e.push((e.clone(), Vec::new()));
        }
        let pos = map_pos(&self.my_spl_e, e);
        &mut self.my_spl_e[pos].1
    }

    /// OCCT AddIntL(const TopoDS_Face& F) (cxx L203-212).
    pub fn add_int_l(&mut self, f: &Shape) -> &mut Vec<Shape> {
        if !map_is_bound(&self.my_data, f) {
            self.my_data.push((f.clone(), FaceData::new()));
        }
        let pos = map_pos(&self.my_data, f);
        self.my_data[pos].1.add_int_l()
    }

    /// OCCT AddOutL(const TopoDS_Face& F) (cxx L215-224).
    pub fn add_out_l(&mut self, f: &Shape) -> &mut Vec<Shape> {
        if !map_is_bound(&self.my_data, f) {
            self.my_data.push((f.clone(), FaceData::new()));
        }
        let pos = map_pos(&self.my_data, f);
        self.my_data[pos].1.add_out_l()
    }

    /// OCCT AddIsoL(const TopoDS_Face& F) (cxx L227-236).
    pub fn add_iso_l(&mut self, f: &Shape) -> &mut Vec<Shape> {
        if !map_is_bound(&self.my_data, f) {
            self.my_data.push((f.clone(), FaceData::new()));
        }
        let pos = map_pos(&self.my_data, f);
        self.my_data[pos].1.add_iso_l()
    }

    /// OCCT AddOutV(const TopoDS_Vertex& V) (lxx L71-74).
    pub fn add_out_v(&mut self, v: &Shape) {
        let id = v.ptr_id();
        if !self.my_out_v.iter().any(|s| s.ptr_id() == id) {
            self.my_out_v.push(v.clone());
        }
    }

    /// OCCT AddIntV(const TopoDS_Vertex& V) (lxx L77-80).
    pub fn add_int_v(&mut self, v: &Shape) {
        let id = v.ptr_id();
        if !self.my_int_v.iter().any(|s| s.ptr_id() == id) {
            self.my_int_v.push(v.clone());
        }
    }

    /// OCCT InitEdge() (cxx L239-248).
    pub fn init_edge(&mut self) {
        // myEIterator.Initialize(myEdgesVertices).
        self.my_e_iterator = 0;
        while self.more_edge() && self.my_edges_vertices[self.my_e_iterator].1.is_empty() {
            // myEIterator.Next().
            self.my_e_iterator += 1;
        }
    }

    /// OCCT MoreEdge() (lxx L84-87).
    pub fn more_edge(&self) -> bool {
        self.my_e_iterator < self.my_edges_vertices.len()
    }

    /// OCCT NextEdge() (cxx L251-258).
    pub fn next_edge(&mut self) {
        // myEIterator.Next().
        self.my_e_iterator += 1;
        while self.more_edge() && self.my_edges_vertices[self.my_e_iterator].1.is_empty() {
            self.my_e_iterator += 1;
        }
    }

    /// OCCT Edge() (lxx L89-92) — TopoDS::Edge(myEIterator.Key()).
    pub fn edge(&self) -> Shape {
        self.my_edges_vertices[self.my_e_iterator].0.clone()
    }

    /// OCCT InitVertex(const TopoDS_Edge& E) (cxx L263-273) — start an
    /// iteration on the vertices of E.
    pub fn init_vertex(&mut self, e: &Shape) {
        if !map_is_bound(&self.my_edges_vertices, e) {
            self.my_edges_vertices.push((e.clone(), Vec::new()));
        }
        let pos = map_pos(&self.my_edges_vertices, e);
        // OCCT myVList = &L — the map owns the list.
        self.my_v_list = pos;
        // myVIterator.Initialize(L).
        self.my_v_iterator = 0;
    }

    /// OCCT MoreVertex() (lxx L98-101).
    pub fn more_vertex(&self) -> bool {
        if self.my_v_list == usize::MAX {
            return false;
        }
        self.my_v_iterator < self.my_edges_vertices[self.my_v_list].1.len()
    }

    /// OCCT NextVertex() (lxx L103-106).
    pub fn next_vertex(&mut self) {
        self.my_v_iterator += 1;
    }

    /// OCCT Vertex() (cxx L277-280) — TopoDS::Vertex(myVIterator.Value().
    /// Vertex()).
    pub fn vertex(&self) -> Shape {
        self.my_edges_vertices[self.my_v_list].1[self.my_v_iterator]
            .vertex()
            .clone()
    }

    /// OCCT Parameter() (cxx L284-287).
    pub fn parameter(&self) -> f64 {
        self.my_edges_vertices[self.my_v_list].1[self.my_v_iterator].parameter()
    }

    /// OCCT InsertBefore(const TopoDS_Vertex& V, const double P) (cxx
    /// L291-296) — insert before the current position.
    pub fn insert_before(&mut self, v: &Shape, p: f64) {
        let vd = VData::new_with(p, v.clone());
        // myVList->InsertBefore(VD, myVIterator) — PInsertBefore links the
        // node before myCurrent and leaves the iterator pointing at the
        // same node; in the Vec encoding the current index shifts by one.
        self.my_edges_vertices[self.my_v_list]
            .1
            .insert(self.my_v_iterator, vd);
        self.my_v_iterator += 1;
    }

    /// OCCT Append(const TopoDS_Vertex& V, const double P) (cxx L299-303).
    pub fn append(&mut self, v: &Shape, p: f64) {
        let vd = VData::new_with(p, v.clone());
        self.my_edges_vertices[self.my_v_list].1.push(vd);
    }
}

impl Default for Data {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Curve3, Line3, Plane, Point3, Surface3, Vec3};
    use rcad_kernel::topods::{BRep, BRepBuilder, Orientation};

    /// The test rig: the unit square face of the topol_tool_brep fixture
    /// plus its two first wire edges and their end vertices, as flat lists.
    fn fixture() -> (BRep, Shape, Vec<Shape>, Vec<Shape>) {
        let mut brep = BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, Point3::new(0.0, 0.0, 0.0), 1e-7);
        let v2 = b.add_vertex(&mut brep, Point3::new(1.0, 0.0, 0.0), 1e-7);
        let v3 = b.add_vertex(&mut brep, Point3::new(1.0, 1.0, 0.0), 1e-7);
        // OCCT BRep_Builder::Add(E, V) stores the start vertex FORWARD and
        // the end vertex REVERSED on the edge.
        let e1 = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)))),
            v1.clone(),
            {
                let mut x = v2.clone();
                x.orientation = Orientation::Reversed;
                x
            },
            [0.0, 1.0],
        );
        let e2 = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(Point3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)))),
            v2.clone(),
            {
                let mut x = v3.clone();
                x.orientation = Orientation::Reversed;
                x
            },
            [0.0, 1.0],
        );
        let wire = brep.add_twire(vec![e1.clone(), e2.clone()]);
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_dir: Vec3::new(1.0, 0.0, 0.0),
                v_dir: Vec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );
        (brep, face, vec![v1, v2, v3], vec![e1, e2])
    }

    /// OCCT anchor: the SplE map round-trip — AddSplE binds an empty list
    /// (cxx L191-200), EdgeHasSplE checks bound + non-empty (cxx L51-58)
    /// and IsSplEEdgeEdge falls back to IsSame when E1 is not split
    /// (cxx L95-112).
    #[test]
    fn data_spl_e_round_trip() {
        let (_brep, _face, vs, es) = fixture();
        let (e1, e2) = (es[0].clone(), es[1].clone());
        let (v1, v2) = (vs[0].clone(), vs[1].clone());
        let mut ds = Data::new();

        assert!(!ds.edge_has_spl_e(&e1));
        ds.add_spl_e(&e1);
        // Bound but empty: not "split" yet.
        assert!(!ds.edge_has_spl_e(&e1));

        ds.add_spl_e(&e1).push(e2.clone());
        assert!(ds.edge_has_spl_e(&e1));
        assert!(!ds.edge_has_spl_e(&e2));
        let list = ds.edge_spl_e(&e1);
        assert_eq!(list.len(), 1);
        assert!(list[0].is_same(&e2));

        // E1 is split: the answer comes from the SplE list only.
        assert!(ds.is_spl_e_edge_edge(&e1, &e2));
        assert!(!ds.is_spl_e_edge_edge(&e1, &e1));
        // E2 is not split: the fallback compares E1.IsSame(E2).
        assert!(ds.is_spl_e_edge_edge(&e2, &e2));
        assert!(!ds.is_spl_e_edge_edge(&e2, &e1));
        // The vertices are unrelated.
        assert!(!ds.is_spl_e_edge_edge(&v1, &v2));
    }

    /// OCCT anchor: the per-face interference lists — AddIntL/AddOutL/
    /// AddIsoL bind a default FaceData on first touch (cxx L203-236) and
    /// the Is*LFaceEdge probes resolve through IsSplEEdgeEdge
    /// (cxx L116-165).
    #[test]
    fn data_face_lists() {
        let (_brep, face, _vs, es) = fixture();
        let (e1, e2) = (es[0].clone(), es[1].clone());
        let mut ds = Data::new();

        assert!(!ds.face_has_int_l(&face));
        assert!(!ds.face_has_out_l(&face));
        assert!(!ds.face_has_iso_l(&face));

        ds.add_spl_e(&e1).push(e2.clone());
        ds.add_int_l(&face).push(e1.clone());
        assert!(ds.face_has_int_l(&face));
        assert!(!ds.face_has_out_l(&face));
        assert!(!ds.face_has_iso_l(&face));
        // IntL holds e1; e1 is split into e2 → only e2 resolves.
        assert!(ds.is_int_l_face_edge(&face, &e2));
        assert!(!ds.is_int_l_face_edge(&face, &e1));

        ds.add_out_l(&face).push(e1.clone());
        assert!(ds.face_has_out_l(&face));
        assert!(ds.is_out_l_face_edge(&face, &e2));
        assert!(!ds.is_out_l_face_edge(&face, &e1));

        ds.add_iso_l(&face).push(e2.clone());
        assert!(ds.face_has_iso_l(&face));
        // e2 is not split: the fallback IsSame answers.
        assert!(ds.is_iso_l_face_edge(&face, &e2));
        assert!(!ds.is_iso_l_face_edge(&face, &e1));

        let int_l = ds.face_int_l(&face);
        assert_eq!(int_l.len(), 1);
        assert!(int_l[0].is_same(&e1));
    }

    /// OCCT anchor: AddOldS binds only once (cxx L181-187) and NewSOldS
    /// falls back to its argument (cxx L167-176).
    #[test]
    fn data_old_s() {
        let (_brep, face, _vs, es) = fixture();
        let (e1, e2) = (es[0].clone(), es[1].clone());
        let mut ds = Data::new();

        // Unbound: NewSOldS(New) == New.
        assert!(ds.new_s_old_s(&e2).is_same(&e2));
        ds.add_old_s(&e2, &e1);
        assert!(ds.new_s_old_s(&e2).is_same(&e1));
        // Already bound: the second AddOldS is ignored.
        ds.add_old_s(&e2, &face);
        assert!(ds.new_s_old_s(&e2).is_same(&e1));
        assert!(ds.new_s_old_s(&face).is_same(&face));
    }

    /// OCCT anchor: the outline/internal vertex sets (lxx L56-80).
    #[test]
    fn data_out_int_v() {
        let (_brep, _face, vs, _es) = fixture();
        let (v1, v2, v3) = (vs[0].clone(), vs[1].clone(), vs[2].clone());
        let mut ds = Data::new();

        assert!(!ds.is_out_v(&v1));
        assert!(!ds.is_int_v(&v2));
        ds.add_out_v(&v1);
        ds.add_int_v(&v2);
        ds.add_out_v(&v1); // NCollection_Map::Add is idempotent.
        assert!(ds.is_out_v(&v1));
        assert!(!ds.is_out_v(&v2));
        assert!(ds.is_int_v(&v2));
        assert!(!ds.is_int_v(&v3));
    }

    /// OCCT anchor: the (edge, vertex) iteration with the empty-entry skip
    /// (cxx L239-258), the vertex list mutation (InitVertex/Append/
    /// InsertBefore, cxx L263-303) and the PInsertBefore iterator semantics
    /// (the current node stays; the index shifts by one).
    #[test]
    fn data_edge_vertex_iteration() {
        let (_brep, _face, vs, es) = fixture();
        let (v1, v2, v3) = (vs[0].clone(), vs[1].clone(), vs[2].clone());
        let (e1, e2) = (es[0].clone(), es[1].clone());
        let mut ds = Data::new();

        // No bound edge yet: the iteration is over.
        ds.init_edge();
        assert!(!ds.more_edge());

        // A bound edge with an empty vertex list is skipped by InitEdge.
        ds.add_spl_e(&e1);
        ds.init_edge();
        assert!(!ds.more_edge());

        // The empty iteration on e1.
        ds.init_vertex(&e1);
        assert!(!ds.more_vertex());
        assert!(ds.my_edges_vertices.iter().any(|(k, _)| k.is_same(&e1)));

        // Append builds the vertex list; InitEdge now lands on e1.
        ds.append(&v1, 0.0);
        ds.init_edge();
        assert!(ds.more_edge());
        assert!(ds.edge().is_same(&e1));
        ds.init_vertex(&e1);
        assert!(ds.more_vertex());
        assert!(ds.vertex().is_same(&v1));
        assert_eq!(ds.parameter(), 0.0);

        // InsertBefore inserts before the current node and the iterator
        // still addresses v1 (PInsertBefore keeps myCurrent).
        ds.insert_before(&v2, 0.5);
        assert!(ds.more_vertex());
        assert!(ds.vertex().is_same(&v1));
        assert_eq!(ds.parameter(), 0.0);
        ds.next_vertex();
        assert!(!ds.more_vertex());

        // The list order is [v2@0.5, v1@0.0]; Append now goes at the end.
        ds.append(&v3, 1.0);
        ds.init_vertex(&e1);
        assert!(ds.vertex().is_same(&v2));
        assert_eq!(ds.parameter(), 0.5);
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v1));
        assert_eq!(ds.parameter(), 0.0);
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v3));
        assert_eq!(ds.parameter(), 1.0);
        ds.next_vertex();
        assert!(!ds.more_vertex());

        // The second edge joins the iteration; NextEdge skips nothing.
        ds.init_vertex(&e2);
        ds.append(&v2, 2.0);
        ds.init_edge();
        assert!(ds.more_edge());
        assert!(ds.edge().is_same(&e1));
        ds.next_edge();
        assert!(ds.more_edge());
        assert!(ds.edge().is_same(&e2));
        ds.next_edge();
        assert!(!ds.more_edge());
    }

    /// OCCT anchor: Clean is a literal no-op (cxx L47-49) and Clear resets
    /// every map (cxx L35-45).
    #[test]
    fn data_clean_clear() {
        let (_brep, face, vs, es) = fixture();
        let (v1, _v2, _v3) = (vs[0].clone(), vs[1].clone(), vs[2].clone());
        let (e1, e2) = (es[0].clone(), es[1].clone());
        let mut ds = Data::new();

        ds.add_spl_e(&e1).push(e2.clone());
        ds.add_int_l(&face).push(e1.clone());
        ds.add_old_s(&e2, &e1);
        ds.add_out_v(&v1);

        ds.clean();
        // Clean() {} — the data survives.
        assert!(ds.edge_has_spl_e(&e1));
        assert!(ds.face_has_int_l(&face));
        assert!(ds.is_out_v(&v1));
        assert!(ds.new_s_old_s(&e2).is_same(&e1));

        ds.clear();
        assert!(!ds.edge_has_spl_e(&e1));
        assert!(!ds.face_has_int_l(&face));
        assert!(!ds.is_out_v(&v1));
        assert!(ds.new_s_old_s(&e2).is_same(&e2));
        // After Clear the iteration state reports no more entries.
        ds.init_edge();
        assert!(!ds.more_edge());
    }
}
