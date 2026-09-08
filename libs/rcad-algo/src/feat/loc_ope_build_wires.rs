// OCCT LocOpe_BuildWires.hxx L28-50 + LocOpe_BuildWires.cxx L17-281 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_BuildWires.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_BuildWires.cxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. TopTools_ShapeMapHasher identity -> the key (TShape ptr, Location);
//    NCollection_IndexedDataMap theMapVE -> indexmap::IndexMap keyed by the
//    shape key with the (key, ancestors) pair; NCollection_IndexedMap mapE ->
//    indexmap::IndexMap; NCollection_Map theMap/Bords/mapV -> the local
//    ShapeSet below (insertion-ordered set with Add/Contains/Remove —
//    OcctShapeMap has no Remove; the OCCT bucket iteration order is not
//    reproduced, the same reduction as OcctShapeMap in
//    brep_feat_builder.rs).
// 2. BRep_Builder MakeCompound/MakeWire/Add -> a method-scoped rcad BRep
//    pool (brep_feat_form.rs precedent); the OCCT incremental
//    MakeWire+Add(edge) sequence maps to the accumulated TWireData edge
//    list in the same order.
// 3. TopExp::MapShapesAndAncestors(C, VERTEX, EDGE, theMapVE) (TopExp.cxx
//    L80-120) -> local re-host below.
// 4. TopExp::Vertices/FirstVertex/LastVertex(E, CumOri=true) ->
//    loc_ope_wires_on_shape::top_exp_*.
// 5. StdFail_NotDone and Standard_ConstructionError raises -> panic with
//    the same names (loc_ope_cs_intersector.rs convention).

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape::{
    top_exp_first_vertex, top_exp_last_vertex, top_exp_vertices, LocOpeWiresOnShape,
};
use indexmap::IndexMap;
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{tshape_flags, Orientation, ShapeType, TShape};
use std::collections::HashMap;
use std::sync::Arc;

pub(crate) type ShapeKey = (u64, u32);

fn shape_key(s: &Shape) -> ShapeKey {
    (s.ptr_id(), s.location)
}

/// OCCT NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> — an
/// insertion-ordered set keyed by the shape identity (architecture
/// difference #1; Remove is required by cxx L187/L503/L530/L566).
pub(crate) struct ShapeSet {
    items: HashMap<ShapeKey, Shape>,
}

impl ShapeSet {
    pub fn new() -> Self {
        ShapeSet {
            items: HashMap::new(),
        }
    }

    /// OCCT NCollection_Map::Add — returns true when newly added.
    pub fn add(&mut self, the_s: &Shape) -> bool {
        if self.items.contains_key(&shape_key(the_s)) {
            return false;
        }
        self.items.insert(shape_key(the_s), the_s.clone());
        true
    }

    /// OCCT NCollection_Map::Contains.
    pub fn contains(&self, the_s: &Shape) -> bool {
        self.items.contains_key(&shape_key(the_s))
    }

    /// OCCT NCollection_Map::Remove.
    pub fn remove(&mut self, the_s: &Shape) {
        self.items.remove(&shape_key(the_s));
    }

    /// OCCT NCollection_Map iterator (the values).
    pub fn iter(&self) -> impl Iterator<Item = &Shape> {
        self.items.values()
    }
}

/// OCCT TopExp::MapShapesAndAncestors(S, TS, TA, M) (TopExp.cxx L80-120) —
/// local re-host over the rcad explorer (architecture difference #3).
fn map_shapes_and_ancestors(
    the_s: &Shape,
    the_ts: ShapeType,
    the_ta: ShapeType,
    the_m: &mut IndexMap<ShapeKey, (Shape, Vec<Shape>)>,
) {
    // visit ancestors (OCCT L89-107).
    for anc in explorer(the_s, the_ta, ShapeType::Shape) {
        for exs in explorer(&anc, the_ts, ShapeType::Shape) {
            if the_m.get_index_of(&shape_key(&exs)).is_none() {
                the_m.insert(shape_key(&exs), (exs.clone(), Vec::new()));
            }
            let idx = the_m.get_index_of(&shape_key(&exs)).expect("entry");
            the_m
                .get_index_mut(idx)
                .expect("entry")
                .1
                 .1
                .push(anc.clone());
        }
    }

    // visit shapes not under ancestors (OCCT L109-119).
    for ex in explorer(the_s, the_ts, the_ta) {
        if the_m.get_index_of(&shape_key(&ex)).is_none() {
            the_m.insert(shape_key(&ex), (ex.clone(), Vec::new()));
        }
    }
}

/// OCCT LocOpe_BuildWires (LocOpe_BuildWires.hxx L28-48).
pub struct LocOpeBuildWires {
    my_done: bool,      // OCCT: myDone
    my_res: Vec<Shape>, // OCCT: myRes
}

impl Default for LocOpeBuildWires {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeBuildWires {
    /// OCCT LocOpe_BuildWires::LocOpe_BuildWires() (cxx L44-47).
    pub fn new() -> Self {
        LocOpeBuildWires {
            my_done: false,
            my_res: Vec::new(),
        }
    }

    /// OCCT LocOpe_BuildWires::LocOpe_BuildWires(L, PW) (cxx L52-56) —
    /// Modified by skv - Mon May 31 12:58:27 2004 OCC5865.
    pub fn with_edges_and_wires(the_l: &[Shape], the_pw: &mut LocOpeWiresOnShape) -> Self {
        let mut w = LocOpeBuildWires::new();
        w.perform(the_l, the_pw);
        w
    }

    /// OCCT LocOpe_BuildWires::Perform(L, PW) (cxx L63-228) — Modified by
    /// skv - Mon May 31 12:59:09 2004 OCC5865.
    pub fn perform(&mut self, the_l: &[Shape], the_pw: &mut LocOpeWiresOnShape) {
        // OCCT L67-68.
        self.my_done = false;
        self.my_res.clear();

        // OCCT L70-72: B; TopoDS_Compound C; B.MakeCompound(C).
        let mut pool = BRep::new();
        let mut c = pool.add_tcompound(Vec::new());

        // OCCT L74-84.
        let mut the_map = ShapeSet::new();
        for itl in the_l {
            let edg = itl.clone();
            if the_map.add(&edg) && edg.shape_type() == ShapeType::Edge {
                // orientation importante pour appel a TopExp::Vertices
                // (OCCT L81-82).
                let mut child = edg.clone();
                child.orientation = Orientation::Forward;
                if let TShape::Compound(cd) = Arc::make_mut(&mut c.data) {
                    cd.push(child);
                }
            }
        }

        // OCCT L86-88.
        let mut the_map_ve: IndexMap<ShapeKey, (Shape, Vec<Shape>)> = IndexMap::new();
        map_shapes_and_ancestors(&c, ShapeType::Vertex, ShapeType::Edge, &mut the_map_ve);

        // OCCT L90-108.
        let mut bords = ShapeSet::new();
        for i in 1..=the_map_ve.len() {
            // Modified by skv - Mon May 31 13:07:50 2004 OCC5865 (OCCT
            // L96-106).
            let entry = the_map_ve.get_index(i - 1).expect("theMapVE entry");
            let vtx = entry.1 .0.clone();
            let mut etmp = Shape::null();
            let mut a_v_border = Shape::null();
            let mut partmp = 0.0;
            if entry.1 .1.is_empty()
                || (the_pw.on_vertex(&vtx, &mut a_v_border)
                    || the_pw.on_edge(&vtx, &mut etmp, &mut partmp))
            {
                bords.add(&vtx);
            }
        }

        // OCCT L110-225.
        loop {
            let i = find_first_edge(&the_map_ve, &bords);
            if i > the_map_ve.len() {
                break;
            }
            let mut map_e: IndexMap<ShapeKey, Shape> = IndexMap::new();
            let mut map_v = ShapeSet::new();
            // OCCT L114: edgf = TopoDS::Edge(theMapVE(i).First()).
            let edgf = {
                let entry = the_map_ve.get_index(i - 1).expect("theMapVE entry");
                entry.1 .1.first().expect("list first").clone()
            };

            // OCCT L116-117.
            let (mut v_f, mut v_l) = top_exp_vertices(&edgf);

            // OCCT L119-129.
            let vl_in_bords = v_l.as_ref().map(|v| bords.contains(v)).unwrap_or(false);
            let vf_in_bords = v_f.as_ref().map(|v| bords.contains(v)).unwrap_or(false);
            if vl_in_bords && !vf_in_bords {
                // OCCT L121-124.
                let mut rev = edgf.clone();
                rev.orientation = Orientation::Reversed;
                map_e.insert(shape_key(&rev), rev);
                let temp = v_f.clone();
                v_f = v_l.clone();
                v_l = temp;
            } else {
                // OCCT L128.
                let mut fwd = edgf.clone();
                fwd.orientation = Orientation::Forward;
                map_e.insert(shape_key(&fwd), fwd);
            }
            if let Some(vf) = &v_f {
                map_v.add(vf);
            }

            // OCCT L132-162.
            loop {
                let vl_in_map_v = v_l.as_ref().map(|v| map_v.contains(v)).unwrap_or(false);
                let vl_in_bords = v_l.as_ref().map(|v| bords.contains(v)).unwrap_or(false);
                if vl_in_map_v || vl_in_bords {
                    break;
                }
                let vl = v_l.clone().expect("VL defined");
                // OCCT L134.
                let ind = the_map_ve
                    .get_index_of(&shape_key(&vl))
                    .expect("theMapVE.FindIndex(VL)");
                // OCCT L135-142.
                let mut found_idx: Option<usize> = None;
                for (k, e) in the_map_ve.get_index(ind).expect("entry").1 .1.iter().enumerate() {
                    if !map_e.contains_key(&shape_key(e)) {
                        found_idx = Some(k);
                        break;
                    }
                }
                // OCCT L144-147.
                let Some(found_idx) = found_idx else {
                    // OCCT: throw Standard_ConstructionError().
                    panic!("Standard_ConstructionError");
                };
                // OCCT L148-150.
                let the_edge = the_map_ve.get_index(ind).expect("entry").1 .1[found_idx].clone();
                let (v_f2, v_l2) = top_exp_vertices(&the_edge);
                // OCCT L151: mapV.Add(VL).
                map_v.add(&vl);
                // OCCT L152-161.
                let vf_same = v_f2.as_ref().map(|v| v.is_same(&vl)).unwrap_or(false);
                if vf_same {
                    let mut fwd = the_edge.clone();
                    fwd.orientation = Orientation::Forward;
                    map_e.insert(shape_key(&fwd), fwd);
                    v_l = v_l2;
                } else {
                    // on doit avoir Vl == VL
                    let mut rev = the_edge.clone();
                    rev.orientation = Orientation::Reversed;
                    map_e.insert(shape_key(&rev), rev);
                    v_l = v_f2;
                }
            }

            // OCCT L164-165: TopoDS_Wire newWire; B.MakeWire(newWire) — the
            // rcad vehicle builds the wire from the accumulated edge list
            // after the OCCT mapE walk below (architecture difference #2);
            // the Adds happen in the same order.
            let mut new_wire_edges: Vec<Shape> = Vec::new();

            // OCCT L167-194.
            let vl_in_map_v = v_l.as_ref().map(|v| map_v.contains(v)).unwrap_or(false);
            if vl_in_map_v {
                // on sort avec une boucle a recreer
                let vl = v_l.clone().expect("VL defined");
                // OCCT L170-172.
                let mut j = 1usize;
                while j <= map_e.len() {
                    // OCCT L174-182.
                    let edg = map_e.get_index(j - 1).expect("mapE entry").1.clone();
                    let v_f = if edg.orientation == Orientation::Forward {
                        top_exp_first_vertex(&edg)
                    } else {
                        top_exp_last_vertex(&edg)
                    };
                    let vf_same_vl = v_f.as_ref().map(|v| v.is_same(&vl)).unwrap_or(false);
                    if vf_same_vl {
                        break;
                    }
                    // OCCT L187: mapV.Remove(Vf).
                    if let Some(vf) = &v_f {
                        map_v.remove(vf);
                    }
                    j += 1;
                }
                // OCCT L189-192: B.Add(newWire, mapE(j)).
                for jj in j..=map_e.len() {
                    let edg = map_e.get_index(jj - 1).expect("mapE entry").1.clone();
                    new_wire_edges.push(edg);
                }
                // OCCT L193: newWire.Closed(true).
                let mut new_wire = pool.add_twire(std::mem::take(&mut new_wire_edges));
                set_closed_flag(&mut new_wire, true);
                // OCCT L205: myRes.Append(newWire).
                self.my_res.push(new_wire);
            } else {
                // on sort sur un bord : wire ouvert... (OCCT L195-203)
                if let Some(vl) = &v_l {
                    // OCCT L197: mapV.Add(VL).
                    map_v.add(vl);
                }
                // OCCT L198-201: B.Add(newWire, mapE(j)).
                for jj in 1..=map_e.len() {
                    let edg = map_e.get_index(jj - 1).expect("mapE entry").1.clone();
                    new_wire_edges.push(edg);
                }
                // OCCT L202: newWire.Closed(false).
                let mut new_wire = pool.add_twire(std::mem::take(&mut new_wire_edges));
                set_closed_flag(&mut new_wire, false);
                // OCCT L205: myRes.Append(newWire).
                self.my_res.push(new_wire);
            }

            // OCCT L206-224.
            let map_v_items: Vec<Shape> = map_v.iter().cloned().collect();
            for vtx in &map_v_items {
                // OCCT L210: Bords.Add(vtx).
                bords.add(vtx);
                // OCCT L211.
                let ind = the_map_ve
                    .get_index_of(&shape_key(vtx))
                    .expect("theMapVE.FindIndex(vtx)");
                // OCCT L212-223: remove the consumed edges from the vertex's
                // ancestor list.
                let consumed: Vec<Shape> =
                    the_map_ve.get_index(ind).expect("entry").1 .1.clone();
                let mut kept: Vec<Shape> = Vec::new();
                for e in consumed {
                    if !map_e.contains_key(&shape_key(&e)) {
                        kept.push(e);
                    }
                }
                the_map_ve.get_index_mut(ind).expect("entry").1 .1 = kept;
            }
        }

        // OCCT L227.
        self.my_done = true;
    }

    /// OCCT LocOpe_BuildWires::IsDone() (cxx L232-235).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT LocOpe_BuildWires::Result() (cxx L239-246).
    pub fn result(&self) -> &Vec<Shape> {
        if !self.my_done {
            // OCCT L241-244: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        &self.my_res
    }
}

/// OCCT static FindFirstEdge(theMapVE, theBord) (cxx L250-281) — the 1-based
/// index of the first non-empty vertex entry, preferring one on the border;
/// Extent()+1 when the map is exhausted.
fn find_first_edge(
    the_map_ve: &IndexMap<ShapeKey, (Shape, Vec<Shape>)>,
    the_bord: &ShapeSet,
) -> usize {
    // OCCT L255-263: the first index with a non-empty list.
    let mut i = 1usize;
    while i <= the_map_ve.len() {
        if !the_map_ve.get_index(i - 1).expect("entry").1 .1.is_empty() {
            break;
        }
        i += 1;
    }
    if i > the_map_ve.len() {
        return i;
    }

    // OCCT L270-279: the first index at or after i whose vertex is on the
    // border.
    let mut goodi = i;
    for j in i..=the_map_ve.len() {
        let entry = the_map_ve.get_index(j - 1).expect("entry");
        if !entry.1 .1.is_empty() && the_bord.contains(&entry.1 .0) {
            goodi = j;
            break;
        }
    }
    goodi
}

/// OCCT TopoDS_Shape::Closed(theFlag) on the wire TShape (cxx L193/L202).
fn set_closed_flag(wire: &mut Shape, flag: bool) {
    if let TShape::Wire(wd) = Arc::make_mut(&mut wire.data) {
        if flag {
            wd.flags |= tshape_flags::CLOSED;
        } else {
            wd.flags &= !tshape_flags::CLOSED;
        }
    }
}
