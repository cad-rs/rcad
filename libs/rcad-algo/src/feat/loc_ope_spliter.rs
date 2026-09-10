// OCCT LocOpe_Spliter.hxx L30-77 + LocOpe_Spliter.cxx L17-718 +
// LocOpe_Spliter.lxx L19-65 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Spliter.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Spliter.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Spliter.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. TopTools_ShapeMapHasher identity -> the key (TShape ptr, Location);
//    NCollection_DataMap myMap -> HashMap<ShapeKey, Vec<Shape>>;
//    NCollection_IndexedMap Emap -> indexmap::IndexMap; the NCollection_Map
//    sets -> the local ShapeSet (insertion-ordered, Add/Contains/Remove/
//    First; the OCCT bucket iteration order is not reproduced — the same
//    reduction as OcctShapeMap in brep_feat_builder.rs).
// 2. BRepTools_Substitution (TKBRep, BRepTools_Substitution.cxx L17-166) is
//    not translated yet as a standalone module; re-hosted below as
//    BRepToolsSubstitution with the OCCT member set and the exact
//    Clear/Substitute/Build/IsCopied/Copy semantics (TopExp.cxx EmptyCopy
//    via the rcad BRep::empty_copy vehicle).
// 3. LocOpe_SplitShape (LocOpe_SplitShape.cxx L17-1776) is a separate,
//    still-deferred translation (its Add/Rebuild core needs the
//    BRepBuilderAPI/FaceRestrictor/AsDes tooling of BRepAlgo stage 2a).
//    Re-hosted below as LocOpeSplitShape with the OCCT member set
//    (LocOpe_SplitShape.hxx L82-86); the trivial members (ctor/Init/Shape/
//    Put) are translated, the splitting core (CanSplit/Add x3/AddOpenWire/
//    AddClosedWire/Rebuild/LeftOf and the Rebuild call inside
//    DescendantShapes) stays deferred — GAP, closes with the
//    LocOpe_SplitShape batch.
// 4. TopExp::Vertices(E, ...)/FirstVertex/LastVertex ->
//    loc_ope_wires_on_shape::top_exp_*; the Wire overload
//    (TopExp.cxx L255-312) -> top_exp_vertices_wire below.
// 5. TopExp::MapShapes(S, T, M) (TopExp.cxx L35-45) -> map_shapes below.
// 6. BRep_Builder UpdateVertex/Add/MakeWire -> Arc::make_mut edits on the
//    payload (bop/algo/builder.rs L6611 precedent); the incremental
//    MakeWire+Add sequence of the degenerate-wire build maps to pushes into
//    the TWireData edge list.
// 7. GeomAPI_ProjectPointOnCurve -> closest_point_on_curve_range (the
//    loc_ope_wires_on_shape.rs arch. diff. #5); the OCCT TopLoc_Location
//    transform of the 3D curves applies only for identity locations
//    (loc_ope_find_edges.rs arch. diff. #1).
// 8. StdFail_NotDone / Standard_ConstructionError / Standard_NullObject
//    raises -> panics with the same names (loc_ope_cs_intersector.rs
//    convention).

use crate::feat::brep_feat_builder::{explorer, sub_shapes};
use crate::feat::loc_ope_build_wires::LocOpeBuildWires;
use crate::feat::loc_ope_wires_on_shape::{top_exp_vertices, LocOpeWiresOnShape};
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use std::collections::HashMap;
use std::sync::Arc;

pub(crate) type ShapeKey = (u64, u32);

fn shape_key(s: &Shape) -> ShapeKey {
    (s.ptr_id(), s.location)
}

/// OCCT NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> — an
/// insertion-ordered set keyed by the shape identity (architecture
/// difference #1).
pub(crate) struct ShapeSet {
    keys: Vec<ShapeKey>,
    items: HashMap<ShapeKey, Shape>,
}

impl ShapeSet {
    pub fn new() -> Self {
        ShapeSet {
            keys: Vec::new(),
            items: HashMap::new(),
        }
    }

    /// OCCT NCollection_Map::Add — returns true when newly added.
    pub fn add(&mut self, the_s: &Shape) -> bool {
        let k = shape_key(the_s);
        if self.items.contains_key(&k) {
            return false;
        }
        self.keys.push(k);
        self.items.insert(k, the_s.clone());
        true
    }

    /// OCCT NCollection_Map::Contains.
    pub fn contains(&self, the_s: &Shape) -> bool {
        self.items.contains_key(&shape_key(the_s))
    }

    /// OCCT NCollection_Map::Remove.
    pub fn remove(&mut self, the_s: &Shape) {
        let k = shape_key(the_s);
        if self.items.remove(&k).is_some() {
            self.keys.retain(|&x| x != k);
        }
    }

    /// OCCT NCollection_Map::Extent.
    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    /// OCCT NCollection_Map::IsEmpty.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// OCCT the first key of the set (the map-iterator head).
    pub fn first(&self) -> Option<&Shape> {
        self.items.get(self.keys.first()?)
    }

    /// OCCT NCollection_Map::Clear.
    pub fn clear(&mut self) {
        self.keys.clear();
        self.items.clear();
    }

    /// OCCT NCollection_Map iterator (the values).
    pub fn iter(&self) -> impl Iterator<Item = &Shape> {
        self.keys
            .iter()
            .map(move |k| self.items.get(k).expect("entry"))
    }
}

// ---------------------------------------------------------------------------
// BRepTools_Substitution re-host (architecture difference #2).
// ---------------------------------------------------------------------------

/// OCCT BRepTools_Substitution (BRepTools_Substitution.hxx L24-40 +
/// BRepTools_Substitution.cxx L17-166).
pub(crate) struct BRepToolsSubstitution {
    my_map: HashMap<ShapeKey, Vec<Shape>>, // OCCT: myMap
}

impl Default for BRepToolsSubstitution {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepToolsSubstitution {
    /// OCCT BRepTools_Substitution::BRepTools_Substitution() (cxx L24).
    pub fn new() -> Self {
        BRepToolsSubstitution {
            my_map: HashMap::new(),
        }
    }

    /// OCCT BRepTools_Substitution::Clear() (cxx L29-33).
    pub fn clear(&mut self) {
        self.my_map.clear();
    }

    /// OCCT BRepTools_Substitution::Substitute(OS, NS) (cxx L36-43).
    pub fn substitute(&mut self, the_os: &Shape, the_ns: Vec<Shape>) {
        if self.is_copied(the_os) {
            // Standard_ConstructionError_Raise_if (cxx L38).
            panic!("Standard_ConstructionError");
        }
        self.my_map.insert(shape_key(the_os), the_ns);
    }

    /// OCCT BRepTools_Substitution::Build(S) (cxx L48-134).
    pub fn build(&mut self, the_s: &Shape) {
        if self.is_copied(the_s) {
            return;
        }

        // OCCT L52: B; L54: TopoDS_Iterator iteS(S.Oriented(TopAbs_FORWARD)).
        let mut is_modified = false; // OCCT L55

        //------------------------------------------
        // look S is modified and build subshapes.
        //------------------------------------------
        let subs = sub_shapes(&oriented(the_s, Orientation::Forward));
        for ss in &subs {
            self.build(ss);
            if self.is_copied(ss) {
                is_modified = true;
            }
        }

        // OCCT L62: TopoDS_Shape NewS = S.Oriented(TopAbs_FORWARD).
        let mut new_s = oriented(the_s, Orientation::Forward);
        if is_modified {
            //----------------------------------------
            // Rebuild S.
            //------------------------------------------
            // OCCT L68: NewS.EmptyCopy().
            new_s = empty_copy_of(the_s);

            // OCCT L70-74: the edge keeps the range.
            if new_s.shape_type() == ShapeType::Edge {
                let (f, l) = brep_tool_range(the_s);
                set_edge_range(&mut new_s, f, l);
            }

            //------------------------------------------
            // Add the copy of subshapes of S to NewS.
            //------------------------------------------
            let mut has_sub_shape = false; // OCCT L91
            for ss in &subs {
                let os = ss.orientation; // OCCT L95
                let l = self.my_map.get(&shape_key(ss)).cloned().unwrap_or_default();
                for nss in &l {
                    //------------------------------------------
                    // Rebuild NSS and add its copy to NewS.
                    //------------------------------------------
                    self.build(nss);

                    let nl = self
                        .my_map
                        .get(&shape_key(nss))
                        .cloned()
                        .unwrap_or_default();
                    let new_or = os.compose(nss.orientation); // OCCT L103
                    for e in &nl {
                        builder_add(&mut new_s, &oriented(e, new_or));
                        has_sub_shape = true;
                    }
                }
            }
            if !has_sub_shape {
                // OCCT L104-117: Wire,Solid,Shell,Compound must have
                // subshape else they disappear.
                if matches!(
                    new_s.shape_type(),
                    ShapeType::Wire
                        | ShapeType::Shell
                        | ShapeType::Solid
                        | ShapeType::Compound
                ) {
                    new_s = Shape::null();
                }
            }
        }

        //-------------------------------------------------------
        // NewS has the same orientation than S in its ancestors
        // so NewS is bound with orientation FORWARD.
        //-------------------------------------------------------
        let l: Vec<Shape> = if !new_s.is_null() {
            vec![oriented(&new_s, Orientation::Forward)]
        } else {
            Vec::new()
        };
        self.substitute(the_s, l);
    }

    /// OCCT BRepTools_Substitution::IsCopied(S) (cxx L139-154).
    pub fn is_copied(&self, the_s: &Shape) -> bool {
        if let Some(l) = self.my_map.get(&shape_key(the_s)) {
            if l.is_empty() {
                return true;
            }
            return !the_s.is_same(&l[0]);
        }
        false
    }

    /// OCCT BRepTools_Substitution::Copy(S) (cxx L159-165) — a clone of the
    /// bound list (the OCCT const-list&).
    pub fn copy(&self, the_s: &Shape) -> Vec<Shape> {
        if !self.is_copied(the_s) {
            // Standard_NoSuchObject_Raise_if (cxx L162).
            panic!("Standard_NoSuchObject");
        }
        self.my_map.get(&shape_key(the_s)).cloned().unwrap_or_default()
    }
}

// ---------------------------------------------------------------------------
// LocOpe_SplitShape re-host (architecture difference #3 — the splitting core
// is deferred; see the header).
// ---------------------------------------------------------------------------

/// OCCT LocOpe_SplitShape (LocOpe_SplitShape.hxx L36-86) — the member set is
/// carried 1:1; the trivial members are translated, the splitting core is
/// deferred (GAP, arch. diff. #3).
pub struct LocOpeSplitShape {
    my_done: bool,                         // OCCT: myDone
    my_shape: Shape,                       // OCCT: myShape
    my_map: HashMap<ShapeKey, Vec<Shape>>, // OCCT: myMap
    my_dbl_e: ShapeSet,                    // OCCT: myDblE
    // OCCT: myLeft — written only by the deferred LeftOf/AddOpenWire/
    // AddClosedWire bodies (arch. diff. #3).
    #[allow(dead_code)]
    my_left: Vec<Shape>,
}

impl LocOpeSplitShape {
    /// OCCT LocOpe_SplitShape::LocOpe_SplitShape() (lxx).
    pub fn new() -> Self {
        LocOpeSplitShape {
            my_done: false,
            my_shape: Shape::null(),
            my_map: HashMap::new(),
            my_dbl_e: ShapeSet::new(),
            my_left: Vec::new(),
        }
    }

    /// OCCT LocOpe_SplitShape::LocOpe_SplitShape(S) (lxx) — = Init(S).
    pub fn with_shape(the_s: &Shape) -> Self {
        let mut s = LocOpeSplitShape::new();
        s.init(the_s);
        s
    }

    /// OCCT LocOpe_SplitShape::Init(S) (cxx L94-101).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_done = false;
        self.my_shape = the_s.clone();
        self.my_dbl_e.clear();
        self.my_map.clear();
        let my_shape = self.my_shape.clone();
        self.put(&my_shape);
    }

    /// OCCT LocOpe_SplitShape::Shape() (lxx).
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT LocOpe_SplitShape::Put(S) (cxx L1355-1373).
    fn put(&mut self, the_s: &Shape) {
        if !self.my_map.contains_key(&shape_key(the_s)) {
            self.my_map.insert(shape_key(the_s), Vec::new());
            if the_s.shape_type() != ShapeType::Vertex {
                for it in sub_shapes(the_s) {
                    self.put(&it);
                }
            } else {
                self.my_map
                    .get_mut(&shape_key(the_s))
                    .expect("entry")
                    .push(the_s.clone());
            }
        }
    }

    /// OCCT LocOpe_SplitShape::DescendantShapes(S) (cxx L1337-1352) — the
    /// structure is kept: the !myDone gate runs Rebuild(myShape), whose body
    /// is the deferred splitting core (arch. diff. #3); until that lands the
    /// map carries the Put bindings only.
    pub fn descendant_shapes(&mut self, the_s: &Shape) -> Vec<Shape> {
        if !self.my_done {
            // OCCT L1341: Rebuild(myShape) — deferred (GAP).
            self.my_done = true;
        }
        // OCCT L1350: return myMap(S).
        self.my_map
            .get(&shape_key(the_s))
            .cloned()
            .unwrap_or_default()
    }

    /// OCCT LocOpe_SplitShape::Add(V, P, E) — DEFERRED (GAP, arch. diff. #3).
    pub fn add_vertex_on_edge(&mut self, _the_v: &Shape, _the_p: f64, _the_e: &Shape) {
        // deferred: the LocOpe_SplitShape splitting core (cxx L17-1776).
    }

    /// OCCT LocOpe_SplitShape::Add(W, F) — DEFERRED (GAP, arch. diff. #3);
    /// the OCCT bool return stays false (the failure output) until then.
    pub fn add_wire_on_face(&mut self, _the_w: &Shape, _the_f: &Shape) -> bool {
        false
    }

    /// OCCT LocOpe_SplitShape::Add(Lwires, F) — DEFERRED (GAP, arch.
    /// diff. #3); the OCCT bool return stays false until then.
    pub fn add_wires_on_face(&mut self, _the_lwires: &[Shape], _the_f: &Shape) -> bool {
        false
    }
}

// ---------------------------------------------------------------------------
// LocOpe_Spliter proper.
// ---------------------------------------------------------------------------

/// OCCT LocOpe_Spliter (LocOpe_Spliter.hxx L30-73).
pub struct LocOpeSpliter {
    my_shape: Shape, // OCCT: myShape
    my_done: bool,   // OCCT: myDone
    my_res: Option<Shape>, // OCCT: myRes (None = null)
    // OCCT: myMap (NCollection_DataMap<Shape, List<Shape>>); the key Shape
    // is carried with the list (OCCT reads the keys at cxx L122/L232).
    my_map: HashMap<ShapeKey, (Shape, Vec<Shape>)>,
    my_dleft: Vec<Shape>, // OCCT: myDLeft
    my_left: Vec<Shape>,  // OCCT: myLeft
}

impl Default for LocOpeSpliter {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeSpliter {
    /// OCCT LocOpe_Spliter::LocOpe_Spliter() (lxx L21-24).
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

    /// OCCT LocOpe_Spliter::LocOpe_Spliter(S) (lxx L28-32).
    pub fn with_shape(the_s: &Shape) -> Self {
        LocOpeSpliter {
            my_shape: the_s.clone(),
            my_done: false,
            my_res: None,
            my_map: HashMap::new(),
            my_dleft: Vec::new(),
            my_left: Vec::new(),
        }
    }

    /// OCCT LocOpe_Spliter::Init(S) (lxx L36-40).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_done = false;
    }

    /// OCCT LocOpe_Spliter::Perform(PW) (cxx L63-576).
    pub fn perform(&mut self, the_pw: &mut LocOpeWiresOnShape) {
        // OCCT L65-68.
        if self.my_shape.is_null() {
            panic!("Standard_NullObject");
        }
        self.my_done = false;
        self.my_map.clear();
        self.my_res = None;

        // OCCT L73: Put(myShape, myMap).
        put(&self.my_shape, &mut self.my_map);

        // OCCT L75-79.
        let mut map_v = ShapeSet::new();
        let mut map_e = ShapeSet::new();
        let mut edg_on_edg: HashMap<ShapeKey, (Shape, Shape)> = HashMap::new();
        let mut map_fe: IndexMap<ShapeKey, (Shape, Vec<Shape>)> = IndexMap::new();

        // 1ere etape : substitution des vertex (OCCT L81)

        let mut vb = Shape::null(); // OCCT L83
        let mut lsubs: Vec<Shape> = Vec::new(); // OCCT L84
        let mut the_subs = BRepToolsSubstitution::new(); // OCCT L85

        // OCCT L88-115.
        the_pw.init_edge_iterator();
        while the_pw.more_edge() {
            let edg = the_pw.edge();
            map_e.add(&edg);
            for exp in explorer(&edg, ShapeType::Vertex, ShapeType::Shape) {
                let vtx = exp;
                if !map_v.contains(&vtx) {
                    if the_pw.on_vertex(&vtx, &mut vb) {
                        map_v.add(&vtx);
                        if vtx.is_same(&vb) {
                            // OCCT L102: continue (the vertex loop).
                            continue;
                        }
                        lsubs.clear();
                        let vsub = oriented(&vtx, Orientation::Forward);
                        let p1 = brep_tool_pnt(&vsub);
                        let p2 = brep_tool_pnt(&vb);
                        let mut d = match (p1, p2) {
                            (Some(a), Some(b)) => (a - b).length(),
                            _ => 0.0,
                        };
                        d += brep_tool_tolerance(&vb);
                        builder_update_vertex(&vsub, d);
                        lsubs.push(vsub);
                        the_subs.substitute(&oriented(&vb, Orientation::Forward), lsubs.clone());
                    }
                }
            }
            the_pw.next_edge();
        }

        // OCCT L117-141.
        the_subs.build(&self.my_shape);
        if the_subs.is_copied(&self.my_shape) {
            // on n`a fait que des substitutions de vertex. Donc chaque
            // element est remplace par lui meme ou par un seul element du
            // meme type.
            let items: Vec<(Shape, Vec<Shape>)> = self
                .my_map
                .values()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect();
            for (key_shape, _v) in items {
                if the_subs.is_copied(&key_shape) {
                    let lsub = the_subs.copy(&key_shape);
                    // OCCT L136-138: myMap(key).Clear(); Append(lsub.First()).
                    let entry = self
                        .my_map
                        .get_mut(&shape_key(&key_shape))
                        .expect("entry");
                    entry.1.clear();
                    if let Some(first) = lsub.first() {
                        entry.1.push(first.clone());
                    }
                }
            }
        }

        // OCCT L143: myRes = myMap(myShape).First().
        let my_res_first = self
            .my_map
            .get(&shape_key(&self.my_shape))
            .and_then(|v| v.1.first().cloned());
        let Some(my_res_first) = my_res_first else {
            // OCCT First() on the empty list raises Standard_NoSuchObject.
            panic!("Standard_NoSuchObject");
        };
        let mut my_res = my_res_first;
        // OCCT L144: LocOpe_SplitShape theCFace(myRes).
        let mut the_cface = LocOpeSplitShape::with_shape(&my_res);

        // Adds every vertices lying on an edge of the shape, and prepares
        // work to rebuild wires on each face (OCCT L146-147).

        let mut ed = Shape::null(); // OCCT L148
        let mut prm = 0.0f64; // OCCT L149

        let mut the_faces_with_section = ShapeSet::new(); // OCCT L151

        // OCCT L152-199.
        the_pw.init_edge_iterator();
        while the_pw.more_edge() {
            let edg = the_pw.edge();
            for exp in explorer(&edg, ShapeType::Vertex, ShapeType::Shape) {
                let vtx = exp;
                if !map_v.contains(&vtx) {
                    map_v.add(&vtx);
                    // OCCT L161: PW->OnEdge(vtx, edg, Ed, prm).
                    if the_pw.on_edge_from(&vtx, &edg, &mut ed, &mut prm) {
                        // on devrait verifier que le vtx n`existe pas deja
                        // sur l`edge
                        let bound = self.my_map.contains_key(&shape_key(&ed));
                        if !bound {
                            // OCCT L185: continue (the vertex loop).
                            continue;
                        }
                        let first = self
                            .my_map
                            .get(&shape_key(&ed))
                            .and_then(|v| v.1.first().cloned())
                            .expect("myMap(Ed).First()");
                        ed = first;
                        the_cface.add_vertex_on_edge(&vtx, prm, &ed);
                    }
                }
            }
            let mut ebis = Shape::null();
            if the_pw.on_edge_iter(&mut ebis) {
                //	EBis = TopoDS::Edge(myMap(Ebis).First()); (OCCT L176)
                edg_on_edg.insert(shape_key(&edg), (edg.clone(), ebis));
            } else {
                let mut fac = the_pw.on_face();
                let bound = self.my_map.contains_key(&shape_key(&fac));
                if !bound {
                    the_pw.next_edge();
                    continue;
                }
                let is_face_with_sec = the_pw.is_face_with_section(&fac);
                let first = self
                    .my_map
                    .get(&shape_key(&fac))
                    .and_then(|v| v.1.first().cloned())
                    .expect("myMap(fac).First()");
                fac = first;
                if is_face_with_sec {
                    the_faces_with_section.add(&fac);
                }
                if !map_fe.contains_key(&shape_key(&fac)) {
                    map_fe.insert(shape_key(&fac), (fac.clone(), Vec::new()));
                }
                map_fe
                    .get_mut(&shape_key(&fac))
                    .expect("entry")
                    .1
                    .push(edg.clone());
            }
            the_pw.next_edge();
        }

        // Rebuilds wires on each face of the shape (OCCT L201).

        // OCCT L203-223.
        for i in 1..=map_fe.len() {
            let (fe_key, fe_entry) = {
                let (k, v) = map_fe.get_index(i - 1).expect("mapFE entry");
                (*k, v.1.clone())
            };
            let fac = map_fe.get(&fe_key).expect("entry").0.clone();
            let mut ledges = fe_entry;
            // Modified by skv - Mon May 31 12:32:54 2004 OCC5865 (OCCT
            // L208-211): RebuildWires(ledges, PW).
            rebuild_wires(&mut ledges, the_pw);
            map_fe.get_mut(&fe_key).expect("entry").1 = ledges.clone();
            if the_faces_with_section.contains(&fac) {
                // OCCT L214: theCFace.Add(ledges, fac).
                the_cface.add_wires_on_face(&ledges, &fac);
            } else {
                // OCCT L218-221.
                for itl in &ledges {
                    the_cface.add_wire_on_face(itl, &fac);
                }
            }
        }

        // Mise a jour des descendants (OCCT L225).

        // OCCT L227-235.
        let items: Vec<(Shape, Vec<Shape>)> = self
            .my_map
            .values()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (sori, v) in items {
            let Some(scib) = v.first().cloned() else {
                continue;
            };
            // OCCT L234: myMap(sori) = theCFace.DescendantShapes(scib).
            let desc = the_cface.descendant_shapes(&scib);
            self.my_map.insert(shape_key(&sori), (sori, desc));
        }

        // OCCT L237: lres = myMap(myShape).
        let lres = self
            .my_map
            .get(&shape_key(&self.my_shape))
            .map(|v| v.1.clone())
            .unwrap_or_default();

        // OCCT L239-269.
        let typ_s = self.my_shape.shape_type();
        if typ_s == ShapeType::Face && lres.len() >= 2 {
            // OCCT L242-250: B.MakeShell(TopoDS::Shell(myRes)) — a NEW shell
            // TShape (the rcad vehicle is a method-scoped pool, arch.
            // diff. #6).
            let mut pool = BRep::new();
            let mut shell = pool.add_tshell(Vec::new());
            shell.orientation = Orientation::Forward;
            for itl in &lres {
                builder_add(&mut shell, &oriented(itl, self.my_shape.orientation));
            }
            my_res = shell;
        } else if typ_s == ShapeType::Edge && lres.len() >= 2 {
            // OCCT L251-261: a wire from the edges.
            let mut pool = BRep::new();
            let mut wire = pool.add_twire(Vec::new());
            wire.orientation = Orientation::Forward;
            for itl in &lres {
                builder_add(&mut wire, &oriented(itl, self.my_shape.orientation));
            }
            my_res = wire;
        } else {
            if lres.len() != 1 {
                // OCCT L266: return (myDone stays false).
                return;
            }
            my_res = lres[0].clone();
        }

        // OCCT L271-344: the EdgOnEdg substitution.
        the_subs.clear();
        let edg_items: Vec<(Shape, Shape)> = edg_on_edg.values().cloned().collect();
        for (e1, e2_base) in &edg_items {
            let e1 = e1.clone();
            // on recherche dans les descendants de e2 l`edge qui correspont
            // a e1 (OCCT L278)

            let (vf1, vl1) = top_exp_vertices(&e1);
            lsubs.clear();
            let descendants = self.bound_of(&e2_base.clone());
            for itl in &descendants {
                let e2 = itl.clone();
                let (vf2, vl2) = top_exp_vertices(&e2);

                let vl1 = vl1.clone();
                let vf1 = vf1.clone();
                let vl1_not_vf1 = !match (&vl1, &vf1) {
                    (Some(a), Some(b)) => a.is_same(b),
                    _ => false,
                };
                if vl1_not_vf1 {
                    let f1f2 = match (&vf1, &vf2) {
                        (Some(a), Some(b)) => a.is_same(b),
                        _ => false,
                    };
                    let l1l2 = match (&vl1, &vl2) {
                        (Some(a), Some(b)) => a.is_same(b),
                        _ => false,
                    };
                    let f1l2 = match (&vf1, &vl2) {
                        (Some(a), Some(b)) => a.is_same(b),
                        _ => false,
                    };
                    let l1f2 = match (&vl1, &vf2) {
                        (Some(a), Some(b)) => a.is_same(b),
                        _ => false,
                    };
                    if f1f2 && l1l2 {
                        lsubs.push(oriented(&e2, Orientation::Forward));
                        //	break;
                    } else if f1l2 && l1f2 {
                        lsubs.push(oriented(&e2, Orientation::Reversed));
                        //	break;
                    }
                } else {
                    // discrimination sur les tangentes (OCCT L302)
                    let f2l2 = match (&vf2, &vl2) {
                        (Some(a), Some(b)) => a.is_same(b),
                        _ => false,
                    };
                    let l2l1 = match (&vl2, &vl1) {
                        (Some(a), Some(b)) => a.is_same(b),
                        _ => false,
                    };
                    if f2l2 && l2l1 {
                        // tout au meme point (OCCT L304)
                        let (c1, f_1, l_1) = match brep_tool_curve(&e1) {
                            Some(v) => v,
                            None => continue,
                        };
                        let (c2, f_2, l_2) = match brep_tool_curve(&e2) {
                            Some(v) => v,
                            None => continue,
                        };
                        let v1 = c1.derivative_at(f_1);
                        let v2 = c2.derivative_at(f_2);
                        let _ = (l_1, l_2);
                        if v1.dot(v2) > 0.0 {
                            lsubs.push(oriented(&e2, Orientation::Forward));
                        } else {
                            lsubs.push(oriented(&e2, Orientation::Reversed));
                        }
                    }
                }
            }
            if lsubs.len() >= 2 {
                // il faut faire un choix (OCCT L328)
                select(&e1, &mut lsubs);
            }
            if lsubs.len() == 1 {
                let ebase = lsubs[0].clone();
                lsubs.clear();
                lsubs.push(oriented(&e1, ebase.orientation));
                the_subs.substitute(&ebase, lsubs.clone());
            } else {
                // OCCT L340-342 (OCCT_DEBUG): "Pb pour substitution".
            }
        }

        // OCCT L346: theSubs.Build(myRes).
        the_subs.build(&my_res);

        // OCCT L348-374.
        let items: Vec<(Shape, Vec<Shape>)> = self
            .my_map
            .values()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (key_shape, ldesc) in items {
            let mut newdesc: Vec<Shape> = Vec::new();
            for itl in &ldesc {
                if the_subs.is_copied(itl) {
                    let lsub = the_subs.copy(itl);
                    if let Some(first) = lsub.first() {
                        newdesc.push(first.clone());
                    }
                } else {
                    newdesc.push(itl.clone());
                }
            }
            self.my_map
                .insert(shape_key(&key_shape), (key_shape, newdesc));
        }

        // OCCT L376-379.
        if the_subs.is_copied(&my_res) {
            let copy = the_subs.copy(&my_res);
            if let Some(first) = copy.first() {
                my_res = first.clone();
            }
        }

        ////remove superfluous vertices on degenerated edges (OCCT L381)
        the_subs.clear();
        let mut emap: IndexMap<ShapeKey, Shape> = IndexMap::new();
        map_shapes(&my_res, ShapeType::Edge, &mut emap);
        let mut deg_edges: Vec<Shape> = Vec::new();
        for i in 1..=emap.len() {
            let an_edge = emap.get_index(i - 1).expect("Emap entry").1.clone();
            if brep_tool_degenerated(&an_edge) {
                deg_edges.push(an_edge);
            }
        }

        let mut pool = BRep::new();
        let mut deg_wires: Vec<Shape> = Vec::new();
        loop {
            // OCCT L397-431.
            if deg_edges.is_empty() {
                break;
            }
            let mut a_deg_wire = pool.add_twire(Vec::new());
            builder_add(&mut a_deg_wire, &deg_edges[0]);
            deg_edges.remove(0);
            loop {
                // OCCT L410: TopExp::Vertices(aDegWire, Vfirst, Vlast).
                let (v_first, v_last) = top_exp_vertices_wire(&a_deg_wire);
                let mut found = false;
                for i in 0..deg_edges.len() {
                    let an_edge = deg_edges[i].clone();
                    let (v1, v2) = top_exp_vertices(&an_edge);
                    let hits = |a: &Option<Shape>, b: &Option<Shape>| match (a, b) {
                        (Some(x), Some(y)) => x.is_same(y),
                        _ => false,
                    };
                    if hits(&v1, &v_first)
                        || hits(&v1, &v_last)
                        || hits(&v2, &v_first)
                        || hits(&v2, &v_last)
                    {
                        builder_add(&mut a_deg_wire, &an_edge);
                        deg_edges.remove(i);
                        found = true;
                        break;
                    }
                }
                if !found {
                    break;
                }
            }
            deg_wires.push(a_deg_wire);
        }

        for i in 1..=deg_wires.len() {
            // OCCT L433-446.
            let mut vmap: IndexMap<ShapeKey, Shape> = IndexMap::new();
            map_shapes(&deg_wires[i - 1], ShapeType::Vertex, &mut vmap);
            let mut lv: Vec<Shape> = Vec::new();
            if let Some(first) = vmap.get_index(0) {
                lv.push(oriented(first.1, Orientation::Forward));
            }
            for j in 2..=vmap.len() {
                let entry = vmap.get_index(j - 1).expect("Vmap entry").1.clone();
                let first_key = vmap.get_index(0).expect("Vmap entry").1.clone();
                if !entry.is_same(&first_key) {
                    the_subs.substitute(&entry, lv.clone());
                }
            }
        }
        // OCCT L447-451.
        the_subs.build(&my_res);
        if the_subs.is_copied(&my_res) {
            let copy = the_subs.copy(&my_res);
            if let Some(first) = copy.first() {
                my_res = first.clone();
            }
        }
        ////

        // OCCT L454-457.
        self.my_dleft.clear();
        self.my_left.clear();
        map_v.clear();

        // OCCT L460-487.
        for exp in explorer(&my_res, ShapeType::Face, ShapeType::Shape) {
            let fac = exp;
            let mut found = false;
            for exp2 in explorer(&fac, ShapeType::Edge, ShapeType::Shape) {
                let edg = exp2;
                for itms in map_e.iter() {
                    if itms.is_same(&edg) && edg.orientation == itms.orientation {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }
            if found {
                self.my_dleft.push(fac.clone());
                self.my_left.push(fac.clone());
            } else {
                map_v.add(&fac);
            }
        }

        // (the OCCT L489-518 commented JAG block is not translated)

        // Map des edges ou les connexions sont possibles (OCCT L520).
        let mut mapebord = ShapeSet::new();
        for itl in &self.my_left {
            for exp in explorer(itl, ShapeType::Edge, ShapeType::Shape) {
                if !map_e.contains(&exp) {
                    if !mapebord.add(&exp) {
                        mapebord.remove(&exp);
                    }
                }
            }
        }

        // OCCT L536-573.
        while mapebord.extent() != 0 {
            let edg = mapebord.first().expect("Mapebord first").clone();

            let mut hit: Option<Shape> = None;
            for itms in map_v.iter() {
                let fac = itms.clone();
                let mut found = false;
                for exp in explorer(&fac, ShapeType::Edge, ShapeType::Shape) {
                    if exp.is_same(&edg) {
                        found = true;
                        break;
                    }
                }
                if found {
                    hit = Some(fac);
                    break; // face a gauche
                }
            }
            if let Some(fac) = hit {
                for exp in explorer(&fac, ShapeType::Edge, ShapeType::Shape) {
                    if !mapebord.add(&exp) {
                        mapebord.remove(&exp);
                    }
                }
                map_v.remove(&fac);
                self.my_left.push(fac);
            } else {
                mapebord.remove(&edg);
            }
        }

        // OCCT L575.
        self.my_res = Some(my_res);
        self.my_done = true;
    }

    /// OCCT LocOpe_Spliter::IsDone() (lxx L44-47).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT LocOpe_Spliter::ResultingShape() (lxx L58-65).
    pub fn resulting_shape(&self) -> Option<&Shape> {
        if !self.my_done {
            // OCCT lxx L61-64: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_res.as_ref()
    }

    /// OCCT LocOpe_Spliter::Shape() (lxx L51-54).
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT LocOpe_Spliter::DirectLeft() (cxx L599-606).
    pub fn direct_left(&self) -> &Vec<Shape> {
        if !self.my_done {
            // OCCT L601-604: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        &self.my_dleft
    }

    /// OCCT LocOpe_Spliter::Left() (cxx L610-617).
    pub fn left(&self) -> &Vec<Shape> {
        if !self.my_done {
            // OCCT L612-615: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        &self.my_left
    }

    /// OCCT LocOpe_Spliter::DescendantShapes(F) (cxx L580-595).
    pub fn descendant_shapes(&mut self, the_f: &Shape) -> Option<&Vec<Shape>> {
        if !self.my_done {
            // OCCT L582-585: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        // OCCT L586-594: myMap(F) or the static empty list.
        self.my_map.get(&shape_key(the_f)).map(|v| &v.1)
    }

    /// myMap(S) access for the perform body (the bound list).
    fn bound_of(&self, the_s: &Shape) -> Vec<Shape> {
        self.my_map
            .get(&shape_key(the_s))
            .map(|v| v.1.clone())
            .unwrap_or_default()
    }
}

// ---------------------------------------------------------------------------
// Static helpers (cxx L51-59 declarations; bodies L623-718).
// ---------------------------------------------------------------------------

/// OCCT static RebuildWires(ledge, PW) (cxx L623-633) — Modified by skv -
/// Mon May 31 12:31:39 2004 OCC5865.
fn rebuild_wires(ledge: &mut Vec<Shape>, the_pw: &mut LocOpeWiresOnShape) {
    let mut the_build = LocOpeBuildWires::new();
    the_build.perform(ledge, the_pw);
    if !the_build.is_done() {
        // OCCT L630-631: throw Standard_ConstructionError().
        panic!("Standard_ConstructionError");
    }
    // OCCT L632: ledge = theBuild.Result().
    *ledge = the_build.result().clone();
}

/// OCCT static Put(S, theMap) (cxx L637-653).
fn put(the_s: &Shape, the_map: &mut HashMap<ShapeKey, (Shape, Vec<Shape>)>) {
    if the_map.contains_key(&shape_key(the_s)) {
        return;
    }
    the_map.insert(shape_key(the_s), (the_s.clone(), Vec::new()));
    the_map
        .get_mut(&shape_key(the_s))
        .expect("entry")
        .1
        .push(the_s.clone());
    for it in sub_shapes(the_s) {
        put(&it, the_map);
    }
}

/// OCCT static Select(Ebase, lsubs) (cxx L657-718).
fn select(the_ebase: &Shape, lsubs: &mut Vec<Shape>) {
    // Choix d`un point (OCCT L660)

    // OCCT L662-666.
    let dmin_init = f64::MAX;
    let mut dmin = dmin_init;
    let mut i = 0usize;
    let mut imin = 0usize;

    // OCCT L668-674: C from Ebase (+transform); Pt at mid-range.
    let Some((c, f, l)) = brep_tool_curve(the_ebase) else {
        // OCCT would dereference a null curve; rcad clears the list (the
        // OCCT imin==0 branch output).
        lsubs.clear();
        return;
    };
    let pt = c.point_at((f + l) / 2.0);

    // OCCT L676-698.
    for itl in lsubs.iter() {
        i += 1;
        let edg = itl.clone();
        let Some((c, f, l)) = brep_tool_curve(&edg) else {
            // OCCT would dereference a null curve inside proj.Init; rcad
            // skips this candidate (NbPoints()==0 branch).
            continue;
        };
        // OCCT L689: proj.Init(Pt, C, f, l).
        let proj = rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
            &c, pt, f, l, 64,
        );
        // OCCT L690-697: proj.NbPoints() > 0.
        {
            if proj.distance < dmin {
                imin = i;
                dmin = proj.distance;
            }
        }
    }

    // OCCT L699-717.
    if imin == 0 {
        lsubs.clear();
    } else {
        // keep only the imin-th element.
        let kept = lsubs[imin - 1].clone();
        lsubs.clear();
        lsubs.push(kept);
    }
}

// ---------------------------------------------------------------------------
// Local re-hosts (arch. diffs. #4/#5/#6).
// ---------------------------------------------------------------------------

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
fn oriented(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT BRep_Tool re-hosts (loc_ope_wires_on_shape_b.rs vehicles).
fn brep_tool_pnt(vtx: &Shape) -> Option<DVec3> {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => Some(vd.point),
        _ => None,
    }
}

fn brep_tool_tolerance(s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

fn brep_tool_range(edg: &Shape) -> (f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT TopExp::Vertices(W, Vfirst, Vlast) (TopExp.cxx L255-312) — the wire
/// overload: add/remove the edge endpoints in a vertex map; closed wires
/// answer with V2, open ones with the FORWARD/REVERSED survivors.
fn top_exp_vertices_wire(the_w: &Shape) -> (Option<Shape>, Option<Shape>) {
    let mut vmap = ShapeSet::new();
    let mut v2_last: Option<Shape> = None;
    let edges = match the_w.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => return (None, None),
    };
    for e in &edges {
        let (mut v1, mut v2) = top_exp_vertices(e);
        if e.orientation == Orientation::Reversed {
            std::mem::swap(&mut v1, &mut v2);
        }
        // OCCT L275-283: add or remove in the vertex map.
        if let Some(a) = &v1 {
            let mut af = a.clone();
            af.orientation = Orientation::Forward;
            if !vmap.add(&af) {
                vmap.remove(&af);
            }
        }
        if let Some(b) = &v2 {
            let mut br = b.clone();
            br.orientation = Orientation::Reversed;
            if !vmap.add(&br) {
                vmap.remove(&br);
            }
            v2_last = v2;
        }
    }
    if vmap.is_empty() {
        // closed (OCCT L288-298)
        let v2 = v2_last.unwrap_or_else(Shape::null);
        return (
            Some(oriented(&v2, Orientation::Forward)),
            Some(oriented(&v2, Orientation::Reversed)),
        );
    }
    if vmap.extent() == 2 {
        // open (OCCT L299-312)
        let mut v_first = Shape::null();
        let mut v_last = Shape::null();
        for k in vmap.iter() {
            if k.orientation == Orientation::Forward {
                v_first = k.clone();
                break;
            }
        }
        for k in vmap.iter() {
            if k.orientation == Orientation::Reversed {
                v_last = k.clone();
                break;
            }
        }
        return (Some(v_first), Some(v_last));
    }
    (None, None)
}

/// OCCT TopExp::MapShapes(S, T, M) (TopExp.cxx L35-45).
fn map_shapes(the_s: &Shape, the_t: ShapeType, the_m: &mut IndexMap<ShapeKey, Shape>) {
    for ex in explorer(the_s, the_t, ShapeType::Shape) {
        the_m.entry(shape_key(&ex)).or_insert_with(|| ex.clone());
    }
}

/// OCCT TopoDS_Shape::EmptyCopied via the rcad BRep::empty_copy vehicle —
/// the source Arc is registered in a method-scoped pool first so the kernel
/// copy logic sees a valid index (architecture difference #2).
pub(crate) fn empty_copy_of(the_s: &Shape) -> Shape {
    let mut pool = BRep::new();
    let idx = pool.tshapes.len();
    pool.tshapes.push(Arc::clone(&the_s.data));
    let probe = Shape::from_parts(Arc::clone(&the_s.data), idx, the_s.location, the_s.orientation);
    pool.empty_copied(&probe)
}

/// OCCT BRep_Builder::Add(parent, child) — append the sub-shape to the
/// container payload (architecture difference #6).
fn builder_add(the_parent: &mut Shape, the_child: &Shape) {
    if let TShape::Wire(wd) = Arc::make_mut(&mut the_parent.data) {
        wd.edges.push(the_child.clone());
        return;
    }
    if let TShape::Shell(sd) = Arc::make_mut(&mut the_parent.data) {
        sd.faces.push(the_child.clone());
        return;
    }
    if let TShape::Solid(sd) = Arc::make_mut(&mut the_parent.data) {
        sd.shells.push(the_child.clone());
        return;
    }
    if let TShape::Compound(cd) = Arc::make_mut(&mut the_parent.data) {
        cd.push(the_child.clone());
        return;
    }
    if let TShape::CompSolid(cs) = Arc::make_mut(&mut the_parent.data) {
        cs.push(the_child.clone());
        return;
    }
    if let TShape::Face(fd) = Arc::make_mut(&mut the_parent.data) {
        if fd.outer_wire.is_null() && the_child.is_wire() {
            fd.outer_wire = the_child.clone();
        } else {
            fd.inner_wires.push(the_child.clone());
        }
        return;
    }
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_parent.data) {
        match the_child.orientation {
            Orientation::Reversed => ed.last = the_child.clone(),
            _ => ed.first = the_child.clone(),
        }
    }
}

/// OCCT BRep_Builder::Range(E, f, l).
fn set_edge_range(the_e: &mut Shape, the_f: f64, the_l: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.range = [the_f, the_l];
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, Tol) — the max-tolerance update.
fn builder_update_vertex(the_v: &Shape, the_tol: f64) {
    let mut data = Arc::clone(&the_v.data);
    if let TShape::Vertex(vd) = Arc::make_mut(&mut data) {
        vd.tolerance = vd.tolerance.max(the_tol);
    }
}
