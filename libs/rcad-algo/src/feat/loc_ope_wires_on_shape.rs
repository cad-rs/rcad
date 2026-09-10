// OCCT LocOpe_WiresOnShape.hxx L38-118 + LocOpe_WiresOnShape.cxx L17-1623 +
// LocOpe_WiresOnShape.lxx L19-36 — 1:1 translation (class part).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_WiresOnShape.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_WiresOnShape.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_WiresOnShape.lxx
//
// The static helpers of the cxx (Project x3, PutPCurve, PutPCurves,
// FindInternalIntersections) and their re-host vehicles live in the sibling
// module loc_ope_wires_on_shape_b (file kept under 2000 lines).
//
// OCCT inheritance chain (hxx L38):
//   LocOpe_WiresOnShape : Standard_Transient (handle-managed)
// Rust has no handles/inheritance -> a plain owned struct; the OCCT
// occ::handle<LocOpe_WiresOnShape> of the consumers (LocOpe_Spliter::Perform,
// LocOpe_BuildWires::Perform) becomes a plain borrow.
//
// Architecture differences (referenced from the affected functions):
// 1. TopTools_ShapeMapHasher identity -> the key (TShape ptr, Location),
//    orientation ignored. NCollection_IndexedDataMap myMapEF ->
//    indexmap::IndexMap; Add (returns the index of an already-bound key, the
//    item is ignored) -> entry().or_insert(); RemoveKey (swap-with-last,
//    NCollection_IndexedMap.hxx L541-575) -> swap_remove.
//    NCollection_DataMap myMap -> HashMap; the key Shape is carried with the
//    value where OCCT reads the key's ShapeType (cxx L355).
//    NCollection_Map -> OcctShapeMap (feat::brep_feat_builder). The OCCT
//    bucket iteration order is not reproduced (HashMap order) — the same
//    reduction as OcctShapeMap in brep_feat_builder.rs.
// 2. BRep_Tool::Curve/Surface/Pnt/Tolerance/Parameter/Degenerated/Range ->
//    TShape data access (loc_ope_find_edges.rs arch. diff. #1);
//    CurveOnSurface -> TEdgeData.pcurves keyed by the face key
//    (loc_ope_gluer.rs arch. diff. #1). BRep_Builder::UpdateVertex/
//    UpdateEdge/Range/Add -> in-place Arc::make_mut edits on the Shape
//    payload (bop/algo/builder.rs L6611 precedent); the OCCT const-ref +
//    shared-TShape mutation becomes &mut / owned-payload threading.
// 3. BRepAdaptor_Curve2d(edg, fac) -> re-host BRepAdaptorCurve2d (module b).
// 4. GeomAdaptor_Surface First/LastU/VParameter -> the surface default
//    domain; U/VResolution -> u/v_resolution_for_surface (loc_ope_gluer.rs
//    arch. diff. #2); UPeriod/VPeriod -> surface_u_period/surface_v_period
//    helpers (module b); IsUPeriodic/IsVPeriodic -> SurfaceEval.
// 5. GeomAPI_ProjectPointOnCurve -> closest_point_on_curve_range
//    (base::geom_api::project; the rcad vehicle always returns one point, so
//    NbPoints()>0 is true by construction);
//    Geom2dAPI_ProjectPointOnCurve -> base::extrema::ExtPC2d;
//    Extrema_ExtPS -> base::extrema::ExtPS (loc_ope_gluer.rs arch. diff. #3;
//    the Initialize(...)+SetFlag(MIN) pair collapses into the constructor).
// 6-14. ShapeAnalysis/ExtCC/ShapeConstruct/GeomProjLib/FClass2d/BRepTools
//    re-hosts — see the module-b header (loc_ope_wires_on_shape_b.rs).
// 15. The OCCT_DEBUG_MESH / OCCT_DEBUG compile-time branches are not
//    translated (not compiled in the reference build).

use crate::feat::brep_feat_builder::{explorer, OcctShapeMap};
use indexmap::IndexMap;
use rcad_kernel::geom::{CurveEval, SurfaceEval};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use std::collections::HashMap;

use super::loc_ope_wires_on_shape_b::{
    brep_tool_curve, brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_parameter,
    brep_tool_range_on_face, brep_tool_surface, brep_tool_tolerance, find_internal_intersections,
    project_vertex_edge, project_vertex_p2d_edge_face, project_vertex_p2d_face_edge,
    put_pcurve, put_pcurves, BRepAdaptorCurve2d, ShapeKey,
};

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
pub(crate) fn shape_key(s: &Shape) -> ShapeKey {
    (s.ptr_id(), s.location)
}

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp / BRepTools re-hosts (arch. diffs. #1/#2/#14).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Pnt(V).
pub(crate) fn brep_tool_pnt(vtx: &Shape) -> Option<glam::DVec3> {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => Some(vd.point),
        _ => None,
    }
}

/// OCCT TopExp::FirstVertex(E, CumOri=true) (TopExp.cxx L182-194) — the
/// vertex sub-shape whose composed orientation is FORWARD (None = the OCCT
/// null vertex).
pub(crate) fn top_exp_first_vertex(edg: &Shape) -> Option<Shape> {
    let ed = match edg.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    if edg.orientation.compose(ed.first.orientation) == Orientation::Forward {
        return Some(ed.first.clone());
    }
    if edg.orientation.compose(ed.last.orientation) == Orientation::Forward {
        return Some(ed.last.clone());
    }
    None
}

/// OCCT TopExp::LastVertex(E, CumOri=true) — the REVERSED-composed vertex.
pub(crate) fn top_exp_last_vertex(edg: &Shape) -> Option<Shape> {
    let ed = match edg.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    if edg.orientation.compose(ed.last.orientation) == Orientation::Reversed {
        return Some(ed.last.clone());
    }
    if edg.orientation.compose(ed.first.orientation) == Orientation::Reversed {
        return Some(ed.first.clone());
    }
    None
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri=true) — FORWARD-composed
/// sub-shape -> Vfirst, REVERSED-composed -> Vlast (TopExp.cxx L540+).
pub(crate) fn top_exp_vertices(edg: &Shape) -> (Option<Shape>, Option<Shape>) {
    (top_exp_first_vertex(edg), top_exp_last_vertex(edg))
}

/// OCCT BRepTools::Compare(V1, V2) (BRepTools.cxx L527-546).
pub(crate) fn brep_tools_compare(v1: &Shape, v2: &Shape) -> bool {
    if v1.is_same(v2) {
        return true;
    }
    let (Some(p1), Some(p2)) = (brep_tool_pnt(v1), brep_tool_pnt(v2)) else {
        return false;
    };
    let l = (p1 - p2).length();
    if l <= brep_tool_tolerance(v1) {
        return true;
    }
    if l <= brep_tool_tolerance(v2) {
        return true;
    }
    false
}

// ---------------------------------------------------------------------------
// OCCT NCollection_DataMap myMap: Bind / IsBound / value access.
// ---------------------------------------------------------------------------

pub(crate) mod my_map {
    use super::*;

    pub type Map = HashMap<ShapeKey, (Shape, Shape)>;

    /// OCCT NCollection_DataMap::Bind.
    pub fn bind(the_map: &mut Map, the_key: &Shape, the_value: &Shape) {
        the_map.insert(shape_key(the_key), (the_key.clone(), the_value.clone()));
    }

    /// OCCT NCollection_DataMap::IsBound.
    pub fn is_bound(the_map: &Map, the_key: &Shape) -> bool {
        the_map.contains_key(&shape_key(the_key))
    }

    /// OCCT NCollection_DataMap::operator() — the value.
    pub fn value<'a>(the_map: &'a Map, the_key: &Shape) -> &'a Shape {
        &the_map.get(&shape_key(the_key)).expect("myMap bound").1
    }

    /// OCCT NCollection_DataMap::Find — a clone of the value.
    pub fn find(the_map: &Map, the_key: &Shape) -> Shape {
        the_map
            .get(&shape_key(the_key))
            .expect("myMap bound")
            .1
            .clone()
    }
}

/// OCCT LocOpe_WiresOnShape (LocOpe_WiresOnShape.hxx L38-116).
pub struct LocOpeWiresOnShape {
    my_shape: Shape, // OCCT: myShape
    // OCCT: myMapEF (NCollection_IndexedDataMap<TopoDS_Shape, TopoDS_Shape>).
    my_map_ef: IndexMap<ShapeKey, (Shape, Shape)>,
    // OCCT: myFacesWithSection (NCollection_Map<TopoDS_Shape>).
    my_faces_with_section: OcctShapeMap,
    my_check_interior: bool, // OCCT: myCheckInterior
    // OCCT: myMap (NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>); the key
    // Shape is carried with the value for the cxx L355 key-type check.
    my_map: my_map::Map,
    my_done: bool, // OCCT: myDone
    my_index: i32, // OCCT: myIndex
}

impl LocOpeWiresOnShape {
    /// OCCT LocOpe_WiresOnShape::LocOpe_WiresOnShape(S) (cxx L98-104).
    pub fn new(the_s: &Shape) -> Self {
        LocOpeWiresOnShape {
            my_shape: the_s.clone(),
            my_map_ef: IndexMap::new(),
            my_faces_with_section: OcctShapeMap::new(),
            my_check_interior: true,
            my_map: HashMap::new(),
            my_done: false,
            my_index: -1,
        }
    }

    /// OCCT LocOpe_WiresOnShape::Init(S) (cxx L108-115).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_check_interior = true;
        self.my_done = false;
        self.my_map.clear();
        self.my_map_ef.clear();
    }

    /// OCCT LocOpe_WiresOnShape::Bind(W, F) (cxx L119-125).
    pub fn bind_wire_face(&mut self, the_w: &Shape, the_f: &Shape) {
        for exp in explorer(the_w, ShapeType::Edge, ShapeType::Shape) {
            self.bind_edge_face(&exp, the_f);
        }
    }

    /// OCCT LocOpe_WiresOnShape::Bind(Comp, F) (cxx L129-136).
    pub fn bind_compound_face(&mut self, the_comp: &Shape, the_f: &Shape) {
        for exp in explorer(the_comp, ShapeType::Edge, ShapeType::Shape) {
            self.bind_edge_face(&exp, the_f);
        }
        // OCCT L135: myFacesWithSection.Add(F).
        self.my_faces_with_section
            .add(shape_key(the_f), the_f.clone());
    }

    /// OCCT LocOpe_WiresOnShape::Bind(E, F) (cxx L140-164).
    pub fn bind_edge_face(&mut self, the_e: &Shape, the_f: &Shape) {
        // OCCT L143: if (!myMapEF.Contains(E)).
        if !self.my_map_ef.contains_key(&shape_key(the_e)) {
            // OCCT L146-153: look for E among the edges of F.
            let mut found = false;
            for exp in explorer(the_f, ShapeType::Edge, ShapeType::Shape) {
                if exp.is_same(the_e) {
                    found = true;
                    break;
                }
            }
            // OCCT L154-158: if (!exp.More()) myMapEF.Add(E, F).
            if !found {
                self.my_map_ef
                    .insert(shape_key(the_e), (the_e.clone(), the_f.clone()));
            }
        } else {
            // OCCT L162: throw Standard_ConstructionError().
            panic!("Standard_ConstructionError");
        }
    }

    /// OCCT LocOpe_WiresOnShape::Bind(Ewir, Efac) (cxx L168-171).
    pub fn bind_edge_edge(&mut self, the_efrom_w: &Shape, the_eon_face: &Shape) {
        my_map::bind(&mut self.my_map, the_efrom_w, the_eon_face);
    }

    /// OCCT LocOpe_WiresOnShape::BindAll() (cxx L175-363).
    pub fn bind_all(&mut self) {
        // OCCT L177-180.
        if self.my_done {
            return;
        }
        let mut the_map = OcctShapeMap::new(); // OCCT L181

        // Detection des vertex a projeter ou a "binder" avec des vertex
        // existants (OCCT L183-219).
        let mut map_v: my_map::Map = HashMap::new();
        let ite_pairs: Vec<(Shape, Shape)> = self.my_map.values().cloned().collect();
        for pair in &ite_pairs {
            let mut eref = pair.0.clone();
            let eimg = pair.1.clone();
            // OCCT L192.
            put_pcurves(&mut eref, &eimg, &self.my_shape);

            for exp in explorer(&eref, ShapeType::Vertex, ShapeType::Shape) {
                let vtx = exp;
                if !the_map.contains(shape_key(&vtx)) {
                    // pas deja traite
                    let mut matched = false;
                    for exp2 in explorer(&eimg, ShapeType::Vertex, ShapeType::Shape) {
                        let vtx2 = exp2;
                        if vtx2.is_same(&vtx) {
                            matched = true;
                            break;
                        } else if brep_tools_compare(&vtx, &vtx2) {
                            my_map::bind(&mut map_v, &vtx, &vtx2);
                            matched = true;
                            break;
                        }
                    }
                    if !matched {
                        // OCCT L212-215: mapV.Bind(vtx, eimg).
                        my_map::bind(&mut map_v, &vtx, &eimg);
                    }
                    // OCCT L216: theMap.Add(vtx).
                    the_map.add(shape_key(&vtx), vtx.clone());
                }
            }
        }

        // OCCT L221-224: mapV -> myMap.
        for (_k, (key, value)) in map_v.iter() {
            my_map::bind(&mut self.my_map, key, value);
        }

        // OCCT L226-252: the Splits / overlapped-edge detection.
        let mut splits: IndexMap<ShapeKey, (Shape, Vec<Shape>)> = IndexMap::new();
        let mut an_overlapped_edges = OcctShapeMap::new();
        for ind in 1..=self.my_map_ef.len() {
            let entry = self
                .my_map_ef
                .get_index(ind - 1)
                .expect("myMapEF entry")
                .1
                .clone();
            let mut edg = entry.0.clone();
            let fac = entry.1.clone();

            // OCCT L235: PutPCurve(edg, fac).
            put_pcurve(&mut edg, &fac);
            // OCCT L237: aPCurve = BRep_Tool::CurveOnSurface(edg, fac, pf, pl).
            let a_pcurve = brep_tool_curve_on_surface(&edg, &fac);
            if a_pcurve.is_none() {
                // OCCT L238-241.
                continue;
            }

            // OCCT L243-251.
            if self.my_check_interior {
                let mut is_overlapped = false;
                find_internal_intersections(&edg, &fac, &mut splits, &mut is_overlapped);
                if is_overlapped {
                    an_overlapped_edges.add(shape_key(&edg), edg.clone());
                }
            }
        }

        // OCCT L254-269: apply the splits (remove the original edge, add the
        // pieces bound to the same face).
        for ind in 1..=splits.len() {
            let (ek, an_edge, new_edges) = {
                let (k, v) = splits.get_index(ind - 1).expect("Splits entry");
                (*k, v.0.clone(), v.1.clone())
            };
            if an_overlapped_edges.contains(ek) {
                continue;
            }
            // OCCT L261: aFace = myMapEF.FindFromKey(anEdge).
            let a_face = match self.my_map_ef.get(&shape_key(&an_edge)) {
                Some(v) => v.1.clone(),
                None => continue,
            };
            // Remove "anEdge" from "myMapEF" (OCCT L263).
            self.my_map_ef.swap_remove(&shape_key(&an_edge));
            // OCCT L264-268.
            for itl in new_edges.iter() {
                self.my_map_ef
                    .entry(shape_key(itl))
                    .or_insert((itl.clone(), a_face.clone()));
            }
        }

        // OCCT L271: aVertParam.
        let mut a_vert_param: HashMap<ShapeKey, f64> = HashMap::new();

        // OCCT L273-350.
        for ind in 1..=self.my_map_ef.len() {
            let entry = self
                .my_map_ef
                .get_index(ind - 1)
                .expect("myMapEF entry")
                .1
                .clone();
            let edg = entry.0.clone();
            let fac = entry.1.clone();
            // JAG 02.02.96 : On verifie les pcurves...
            // (the OCCT L279 PutPCurve call is commented out in the source).

            for exp in explorer(&edg, ShapeType::Vertex, ShapeType::Shape) {
                let mut vtx = exp;
                if the_map.contains(shape_key(&vtx)) {
                    continue;
                }

                // OCCT L289.
                let vtx_param = brep_tool_parameter(&vtx, &edg);
                // OCCT L290.
                let ba_curve2d = BRepAdaptorCurve2d::new(&edg, &fac);

                // OCCT L292-294.
                let p2d = if !ba_curve2d.is_null() {
                    ba_curve2d.value(vtx_param)
                } else {
                    glam::DVec2::new(
                        rcad_kernel::precision::INFINITE_VALUE,
                        rcad_kernel::precision::INFINITE_VALUE,
                    )
                };

                let mut e_pro = Shape::null(); // OCCT L296
                let mut prm = rcad_kernel::precision::INFINITE_VALUE; // OCCT L297
                let is_projected = my_map::is_bound(&self.my_map, &vtx); // OCCT L298

                // OCCT L300-315: if the vertex was already projected on the
                // current edge on the previous face it is necessary to check
                // tolerance of the vertex in the 2D space on the current face
                // without projection and update tolerance of vertex if it is
                // necessary.
                if is_projected {
                    let a_sh = my_map::find(&self.my_map, &vtx);
                    if a_sh.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    e_pro = a_sh;
                    if let Some(p) = a_vert_param.get(&shape_key(&vtx)) {
                        prm = *p;
                    }
                }
                // OCCT L316.
                let ok = project_vertex_p2d_face_edge(&mut vtx, p2d, &fac, &mut e_pro, &mut prm);
                if ok && !is_projected {
                    // OCCT L320-342.
                    let mut matched = false;
                    for exp2 in explorer(&e_pro, ShapeType::Vertex, ShapeType::Shape) {
                        let vtx2 = exp2;
                        if vtx2.is_same(&vtx) {
                            my_map::bind(&mut self.my_map, &vtx, &vtx2);
                            the_map.add(shape_key(&vtx), vtx.clone());
                            matched = true;
                            break;
                        } else if brep_tools_compare(&vtx, &vtx2) {
                            // OCCT L331-332.
                            if let Some((a_f1, a_l1)) = brep_tool_range_on_face(&e_pro, &fac) {
                                if !brep_tool_degenerated(&e_pro)
                                    && ((prm - a_f1).abs() <= rcad_kernel::precision::PCONFUSION
                                        || (prm - a_l1).abs()
                                            <= rcad_kernel::precision::PCONFUSION)
                                {
                                    my_map::bind(&mut self.my_map, &vtx, &vtx2);
                                    the_map.add(shape_key(&vtx), vtx.clone());
                                    matched = true;
                                    break;
                                }
                            }
                        }
                    }
                    // OCCT L343-347.
                    if !matched {
                        my_map::bind(&mut self.my_map, &vtx, &e_pro);
                        a_vert_param.insert(shape_key(&vtx), prm);
                    }
                }
            }
        }

        // Modified by Sergey KHROMOV (OCCT L352-360): add the edge-keyed
        // myMap entries into myMapEF.
        let ite_pairs: Vec<(Shape, Shape)> = self.my_map.values().cloned().collect();
        for (key, value) in &ite_pairs {
            if key.shape_type() == ShapeType::Edge {
                self.my_map_ef
                    .entry(shape_key(key))
                    .or_insert((key.clone(), value.clone()));
            }
        }

        // OCCT L362.
        self.my_done = true;
    }

    /// OCCT LocOpe_WiresOnShape::InitEdgeIterator() (cxx L367-372).
    pub fn init_edge_iterator(&mut self) {
        self.bind_all();
        self.my_index = 1;
    }

    /// OCCT LocOpe_WiresOnShape::MoreEdge() (cxx L376-380).
    pub fn more_edge(&self) -> bool {
        self.my_index <= self.my_map_ef.len() as i32
    }

    /// OCCT LocOpe_WiresOnShape::Edge() (cxx L384-388).
    pub fn edge(&self) -> Shape {
        self.my_map_ef
            .get_index((self.my_index - 1) as usize)
            .expect("myMapEF entry")
            .1
            .0
            .clone()
    }

    /// OCCT LocOpe_WiresOnShape::OnFace() (cxx L392-396) — the face of the
    /// shape on which the current edge is projected.
    pub fn on_face(&self) -> Shape {
        self.my_map_ef
            .get_index((self.my_index - 1) as usize)
            .expect("myMapEF entry")
            .1
            .1
            .clone()
    }

    /// OCCT LocOpe_WiresOnShape::OnEdge(E) (cxx L400-410) — if the current
    /// edge is projected on an edge, returns true and sets E.
    pub fn on_edge_iter(&self, the_e: &mut Shape) -> bool {
        let key = *self
            .my_map_ef
            .get_index((self.my_index - 1) as usize)
            .expect("myMapEF entry")
            .0;
        if let Some(entry) = self.my_map.get(&key) {
            *the_e = entry.1.clone();
            return true;
        }
        false
    }

    /// OCCT LocOpe_WiresOnShape::NextEdge() (cxx L414-418).
    pub fn next_edge(&mut self) {
        self.my_index += 1;
    }

    /// OCCT LocOpe_WiresOnShape::OnVertex(Vw, Vs) (cxx L422-434).
    pub fn on_vertex(&self, the_vw: &Shape, the_vs: &mut Shape) -> bool {
        if my_map::is_bound(&self.my_map, the_vw) {
            let v = my_map::value(&self.my_map, the_vw);
            if v.shape_type() == ShapeType::Vertex {
                *the_vs = v.clone();
                return true;
            }
            return false;
        }
        false
    }

    /// OCCT LocOpe_WiresOnShape::OnEdge(V, Ed, prm) (cxx L438-448).
    pub fn on_edge(&self, the_v: &Shape, the_ed: &mut Shape, the_prm: &mut f64) -> bool {
        if !my_map::is_bound(&self.my_map, the_v)
            || my_map::value(&self.my_map, the_v).shape_type() == ShapeType::Vertex
        {
            return false;
        }
        let ed = my_map::find(&self.my_map, the_v);
        *the_ed = ed.clone();
        *the_prm = project_vertex_edge(the_v, &ed);
        true
    }

    /// OCCT LocOpe_WiresOnShape::OnEdge(V, EdgeFrom, Ed, prm) (cxx L452-487).
    pub fn on_edge_from(
        &self,
        the_v: &Shape,
        the_edge_from: &Shape,
        the_ed: &mut Shape,
        the_prm: &mut f64,
    ) -> bool {
        if !my_map::is_bound(&self.my_map, the_v)
            || my_map::value(&self.my_map, the_v).shape_type() == ShapeType::Vertex
        {
            return false;
        }
        let ed = my_map::find(&self.my_map, the_v);
        *the_ed = ed.clone();
        if !self.my_map_ef.contains_key(&shape_key(the_edge_from)) {
            return false;
        }

        let a_shape = self
            .my_map_ef
            .get(&shape_key(the_edge_from))
            .expect("myMapEF entry")
            .1
            .clone();
        let a_c = brep_tool_curve(&ed);
        if a_c.is_none() && a_shape.shape_type() == ShapeType::Face {
            // OCCT L474-477.
            let a_face = a_shape;
            let vtx_param = brep_tool_parameter(the_v, the_edge_from);
            let ba_curve2d = BRepAdaptorCurve2d::new(the_edge_from, &a_face);
            let p2d = ba_curve2d.value(vtx_param);
            // OCCT L479.
            *the_prm = project_vertex_p2d_edge_face(the_v, p2d, &ed, &a_face);
        } else {
            // OCCT L483.
            *the_prm = project_vertex_edge(the_v, &ed);
        }
        true
    }

    /// OCCT LocOpe_WiresOnShape::IsFaceWithSection(aFace) (lxx L33-36).
    pub fn is_face_with_section(&self, the_a_face: &Shape) -> bool {
        self.my_faces_with_section.contains(shape_key(the_a_face))
    }

    /// OCCT LocOpe_WiresOnShape::Add(theEdges) (cxx L1515-1623) — add
    /// splitting edges or wires for the whole initial shape.
    pub fn add(&mut self, the_edges: &[Shape]) -> bool {
        use glam::{DAffine3, DVec2};
        use rcad_kernel::math::bnd::BndBox;
        use rcad_kernel::topods::State;

        // OCCT L1517-1519.
        let mut an_edges: Vec<Shape> = Vec::new();
        let nb = the_edges.len() as i32;
        let mut an_edge_boxes: Vec<BndBox> = (0..the_edges.len()).map(|_| BndBox::new()).collect();
        for i in 1..=nb {
            // OCCT L1522: aCurSplit = theEdges(i).
            let a_cur_split = &the_edges[(i - 1) as usize];
            for a_cur_e in explorer(a_cur_split, ShapeType::Edge, ShapeType::Shape) {
                // OCCT L1528-1530.
                let mut a_box_e = BndBox::new();
                crate::topalgo::brep_bnd_lib::BRepBndLib::add_close(&a_cur_e, &mut a_box_e);
                if a_box_e.is_void() {
                    continue;
                }
                // OCCT L1534-1537.
                let a_tol_e = brep_tool_tolerance(&a_cur_e);
                a_box_e.set_gap(a_tol_e);
                an_edge_boxes[(i - 1) as usize] = a_box_e;
                an_edges.push(a_cur_e);
            }
        }
        // OCCT L1540-1543.
        let mut an_used_edges: std::collections::HashSet<i32> = std::collections::HashSet::new();
        for a_cur_f in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            // OCCT L1546-1551.
            let mut a_box_f = BndBox::new();
            crate::topalgo::brep_bnd_lib::BRepBndLib::add(&a_cur_f, &mut a_box_f, false);
            if a_box_f.is_void() {
                continue;
            }
            // OCCT L1552: BRepAdaptor_Surface anAdF(aCurF, false) — its UV
            // domain feeds the extrema projector (arch. diff. #5).
            let Some(an_adf_surf) = brep_tool_surface(&a_cur_f) else {
                continue;
            };
            let dom = an_adf_surf.default_domain();
            let (fu, lu, fv, lv) = (dom[0], dom[1], dom[2], dom[3]);
            let mut a_check_state_tool: Option<
                crate::topalgo::brep_top_adaptor::fclass2d::FClass2d,
            > = None;

            // OCCT L1556-1564: the extrema projector initialized once per
            // face (the rcad constructor carries Initialize+SetFlag(MIN),
            // loc_ope_gluer.rs arch. diff. #3).
            let tol_u = rcad_kernel::topo::topods::u_resolution_for_surface(
                &an_adf_surf,
                rcad_kernel::precision::CONFUSION,
            );
            let tol_v = rcad_kernel::topo::topods::v_resolution_for_surface(
                &an_adf_surf,
                rcad_kernel::precision::CONFUSION,
            );

            // OCCT L1566-1620.
            let nb2 = an_edge_boxes.len() as i32;
            for i in 1..=nb2 {
                if an_used_edges.contains(&i) {
                    continue;
                }
                // OCCT L1575.
                if a_box_f.is_out_box(&an_edge_boxes[(i - 1) as usize]) {
                    continue;
                }
                // OCCT L1580.
                let a_cur_e = an_edges[(i - 1) as usize].clone();
                // OCCT L1583: aC = BRep_Tool::Curve(aCurE, aF, aL).
                let Some((a_c, a_f, a_l)) = brep_tool_curve(&a_cur_e) else {
                    // OCCT L1585-1587.
                    an_used_edges.insert(i);
                    continue;
                };
                // OCCT L1589-1590.
                let a_p = a_c.point_at((a_f + a_l) * 0.5);
                let an_extr = rcad_kernel::base::extrema::ExtPS::with_domain(
                    a_p, &an_adf_surf, fu, lu, fv, lv, tol_u, tol_v,
                );

                // OCCT L1592-1594.
                if !an_extr.is_done() || an_extr.nb_ext() == 0 {
                    continue;
                }
                // OCCT L1596-1597.
                let a_tol_e = brep_tool_tolerance(&a_cur_e);
                let a_tol2 = (a_tol_e + rcad_kernel::precision::CONFUSION)
                    * (a_tol_e + rcad_kernel::precision::CONFUSION);
                // OCCT L1598-1619.
                for n in 1..=an_extr.nb_ext() {
                    let a_dist2 = an_extr.square_distance(n);
                    if a_dist2 > a_tol2 {
                        continue;
                    }
                    let a_ps = an_extr.point(n);
                    let (a_u, a_v) = (a_ps.u, a_ps.v);

                    if a_check_state_tool.is_none() {
                        // OCCT L1612: new BRepTopAdaptor_FClass2d(aCurF,
                        // Precision::PConfusion()) (arch. diff. #8).
                        let src = crate::topalgo::shape_source::FaceShapeSource::new(
                            &a_cur_f,
                            an_adf_surf.clone(),
                            &[DAffine3::IDENTITY],
                        );
                        a_check_state_tool = Some(
                            crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
                                &src,
                                0,
                                rcad_kernel::precision::PCONFUSION,
                            ),
                        );
                    }
                    // OCCT L1614: Perform(gp_Pnt2d(aU, aV)) == TopAbs_IN.
                    let src = crate::topalgo::shape_source::FaceShapeSource::new(
                        &a_cur_f,
                        an_adf_surf.clone(),
                        &[DAffine3::IDENTITY],
                    );
                    let in_face = a_check_state_tool
                        .as_ref()
                        .expect("FClass2d")
                        .perform(&src, DVec2::new(a_u, a_v), false)
                        == State::In;
                    if in_face {
                        // OCCT L1616: Bind(aCurE, aCurF).
                        self.bind_edge_face(&a_cur_e, &a_cur_f);
                        an_used_edges.insert(i);
                    }
                }
            }
        }
        // OCCT L1622.
        !an_used_edges.is_empty()
    }

    /// OCCT LocOpe_WiresOnShape::SetCheckInterior(ToCheckInterior) (lxx
    /// L19-22).
    pub fn set_check_interior(&mut self, the_to_check_interior: bool) {
        self.my_check_interior = the_to_check_interior;
    }

    /// OCCT LocOpe_WiresOnShape::IsDone() (lxx L26-29).
    pub fn is_done(&self) -> bool {
        self.my_done
    }
}
