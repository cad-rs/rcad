// OCCT LocOpe_SplitShape.hxx L36-92 + LocOpe_SplitShape.cxx L17-1776 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_SplitShape.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_SplitShape.cxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. TopTools_ShapeMapHasher identity -> the key (TShape ptr, Location),
//    orientation ignored. NCollection_DataMap myMap -> HashMap<ShapeKey,
//    (Shape, Vec<Shape>)> (the key Shape travels with the list because the
//    OCCT iteration reads itm.Key(), cxx L125-127). NCollection_Map myDblE ->
//    feat::loc_ope_spliter::ShapeSet. NCollection_List myLeft -> Vec<Shape>.
//    NCollection_IndexedMap PossE / NCollection_Map MapE -> ShapeSet, whose
//    first()/extent()/add()/remove() carry the OCCT carriers 1:1.
// 2. BRep_Builder mutations -> the crate::brep_algo::tool in-place re-hosts
//    (feat loc_ope_spliter.rs arch. diff. #6).
// 3. TopoDS_Shape::EmptyCopied / BRep_Builder::MakeWire allocate through the
//    owned BRep pool (myPool) so the copies carry a VALID TShape index; the
//    pool-free (index == usize::MAX) builder shapes produced by the
//    brep_algo::tool vehicles are read as null by Shape::is_null and would
//    make builder_add_face_wire overwrite the outer wire. This is the
//    documented pool-free trap of the feat package.
// 4. BRepTools::Update(F) (BRepTools.cxx L383-390) refreshes the cached UV
//    points of the face edges (UpdateFaceUVPoints -> BRep_GCurve::Update /
//    BRep_CurveOnSurface::Update, both of which mutate a UV-point cache). The
//    rcad TEdgeData has no such cache: BRep_Tool::UVPoints is computed from
//    the pcurves on demand (brep_algo::tool::brep_tool_uv_points), so the call
//    is a no-op (brep_tools_update below).
// 5. BRepTopAdaptor_FClass2d(face, Tol) -> topalgo::brep_top_adaptor::
//    fclass2d::FClass2d over FaceShapeSource (feat loc_ope_wires_on_shape_b.rs
//    arch. diff. #8).
// 6. BRepAdaptor_Surface -> the local BRepAdaptorSurface re-host below
//    (IsUPeriodic/IsVPeriodic/UPeriod/VPeriod/D0/GetType; feat
//    loc_ope_gluer.rs arch. diffs. #2/#4).
// 7. Geom2dAPI_ProjectPointOnCurve -> rcad_kernel::base::extrema::ExtPC2d
//    (feat loc_ope_wires_on_shape_b.rs arch. diff. #5).
// 8. BRepLib_MakeWire -> topalgo::brep_lib_make_wire::MakeWire (landed with
//    this batch; it was a GAP carrier before).
// 9. BRepTools_WireExplorer -> topalgo::brep_tools_wire_explorer::WireExplorer
//    (landed with this batch; it was a GAP carrier before).
// 10. BRepTools::OuterWire / BRepTools::UVBounds(F, W, ...) -> the re-hosts
//    below (BRepTools.cxx L556-590 / L137-181).
// 11. StdFail_NotDone / Standard_NoSuchObject raises -> panics with the same
//    names (feat loc_ope_spliter.rs arch. diff. #8).
// 12. The OCCT_DEBUG / OCCT_DEBUG_MESH compile-time branches are not
//    translated (not compiled in the reference build).

use crate::brep_algo::tool::{
    brep_tool_curve_on_surface, brep_tool_parameter, brep_tool_range, brep_tool_tolerance,
    builder_add_edge_vertex, builder_add_face_wire, builder_add_wire_edge, builder_set_closed,
    empty_copied, oriented, sub_shapes, top_exp_vertices_raw, top_exp_vertices_wire,
};
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope::closed_wire as loc_ope_closed_wire;
use crate::feat::loc_ope_spliter::ShapeSet;
use crate::feat::loc_ope_wires_on_shape::{top_exp_first_vertex, top_exp_last_vertex};
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_degenerated, shape_analysis_edge_check_same_parameter,
};
use crate::topalgo::brep_lib_make_wire::{
    builder_update_vertex_parameter, MakeWire,
};
use crate::topalgo::brep_tools_wire_explorer::{
    shape_key, vec2_angle, ShapeKeySet, WireExplorer,
};
use crate::feat::loc_ope_split_shape_b::{
    add_vertex_and_update, brep_tool_is_closed_on_surface, brep_tool_is_closed_shape,
    brep_tools_outer_wire, brep_tools_update, builder_add_shape, builder_continuity,
    builder_update_edge_tol, classif_perform_infinite_point, classif_perform_p,
    face_classifier, has_continuity, is_negative_infinite, is_positive_infinite,
    make_wire_closed, reverse_orientation, sub_shapes_no_cum_ori, BRepAdaptorSurface,
};
use glam::DVec2;
use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::geom::Curve2dEval;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use std::collections::HashMap;
use std::sync::Arc;

/// Shape identity key — re-exported for the sibling modules of the feat
/// package (TopTools_ShapeMapHasher: TShape + Location; arch. diff. #1).
pub(crate) use crate::topalgo::brep_tools_wire_explorer::ShapeKey;

/// OCCT LocOpe_SplitShape (LocOpe_SplitShape.hxx L38-92).
pub struct LocOpeSplitShape {
    my_done: bool,                                       // OCCT: myDone
    my_shape: Shape,                                     // OCCT: myShape
    my_map: HashMap<ShapeKey, (Shape, Vec<Shape>)>,      // OCCT: myMap
    my_dbl_e: ShapeSet,                                  // OCCT: myDblE
    my_left: Vec<Shape>,                                 // OCCT: myLeft
    /// OCCT: the TShape allocator of TopoDS_Shape::EmptyCopied /
    /// BRep_Builder::MakeWire (architecture difference #3).
    my_pool: BRep,
}

impl Default for LocOpeSplitShape {
    fn default() -> Self {
        Self::new()
    }
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
            my_pool: BRep::new(),
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

    // -----------------------------------------------------------------------
    // myMap / myLeft carriers.
    // -----------------------------------------------------------------------

    /// OCCT NCollection_DataMap::operator()(S) on myMap (raises
    /// Standard_NoSuchObject when the key is not bound — NCollection_DataMap.hxx).
    fn map_list(&self, the_s: &Shape) -> &Vec<Shape> {
        &self
            .my_map
            .get(&shape_key(the_s))
            .unwrap_or_else(|| panic!("Standard_NoSuchObject"))
            .1
    }

    /// OCCT NCollection_DataMap::ChangeFind(S) on myMap.
    fn map_list_mut(&mut self, the_s: &Shape) -> &mut Vec<Shape> {
        &mut self
            .my_map
            .get_mut(&shape_key(the_s))
            .unwrap_or_else(|| panic!("Standard_NoSuchObject"))
            .1
    }

    /// OCCT `myMap(theS).Append(theV)`.
    fn map_append(&mut self, the_s: &Shape, the_v: &Shape) {
        self.map_list_mut(the_s).push(the_v.clone());
    }

    /// OCCT `myMap(theS).Clear()`.
    fn map_clear_list(&mut self, the_s: &Shape) {
        self.map_list_mut(the_s).clear();
    }

    /// OCCT `myMap(theS) = theList`.
    fn map_assign(&mut self, the_s: &Shape, the_list: Vec<Shape>) {
        *self.map_list_mut(the_s) = the_list;
    }

    // -----------------------------------------------------------------------
    // Pool allocators (architecture difference #3).
    // -----------------------------------------------------------------------

    /// OCCT TopoDS_Shape::EmptyCopied — the copy is registered in myPool so
    /// its TShape index is valid.
    fn empty_copied(&mut self, the_s: &Shape) -> Shape {
        let idx = self.my_pool.tshapes.len();
        self.my_pool.tshapes.push(Arc::clone(&the_s.data));
        let probe = Shape::from_parts(
            Arc::clone(&the_s.data),
            idx,
            the_s.location,
            the_s.orientation,
        );
        self.my_pool.empty_copied(&probe)
    }

    /// OCCT BRep_Builder::MakeWire(W) — an empty wire in myPool.
    fn make_wire(&mut self) -> Shape {
        self.my_pool.add_twire(Vec::new())
    }

    // -----------------------------------------------------------------------
    // CanSplit (cxx L105-139).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::CanSplit(E) (cxx L105-139).
    fn can_split(&self, the_e: &Shape) -> bool {
        // OCCT L107-110.
        if self.my_done {
            return false;
        }
        if self.my_map.is_empty() {
            return false;
        }

        // OCCT L116-119.
        if !self.my_map.contains_key(&shape_key(the_e)) {
            return false;
        }

        // On verifie que l`edge n`appartient pas a un wire deja reconstruit
        // (OCCT L121-137).
        for (_k, (key, value)) in self.my_map.iter() {
            if key.shape_type() == ShapeType::Wire && !value.is_empty() {
                for exp in explorer(key, ShapeType::Edge, ShapeType::Shape) {
                    if exp.is_same(the_e) {
                        return false;
                    }
                }
            }
        }
        true
    }

    // -----------------------------------------------------------------------
    // Add(V, P, E) (cxx L143-249).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::Add(const TopoDS_Vertex& V, const double P,
    /// const TopoDS_Edge& E) (cxx L143-249).
    pub fn add_vertex_on_edge(&mut self, the_v: &Shape, the_p: f64, the_e: &Shape) {
        // OCCT L145-148.
        if !self.can_split(the_e) {
            return;
        }

        // OCCT L151-155: NCollection_List<TopoDS_Shape>& le = myMap(E);
        // if (le.IsEmpty()) le.Append(E);
        if self.map_list(the_e).is_empty() {
            self.map_append(the_e, the_e);
        }

        // OCCT L156-172: walk `le` until an edge holding P strictly inside is
        // found; the skipped edges accumulate in aNewList.
        let mut le: Vec<Shape> = self.map_list(the_e).clone();
        let mut a_new_list: Vec<Shape> = Vec::new();
        let mut found: Option<usize> = None;
        for (i, itl) in le.iter().enumerate() {
            let edg = itl.clone();
            let (f, l) = brep_tool_range(&edg);
            if the_p > f + PCONFUSION && the_p < l - PCONFUSION {
                found = Some(i);
                break;
            }
            a_new_list.push(edg);
        }
        let Some(idx) = found else {
            // OCCT L169-172: `if (!itl.More()) return;` — le keeps the
            // appended E.
            self.map_assign(the_e, le);
            return;
        };
        let mut edg = le.remove(idx); // OCCT L174: le.Remove(itl)

        if the_v.orientation == Orientation::Forward
            || the_v.orientation == Orientation::Reversed
        {
            // OCCT L178.
            let edg = oriented(&edg, Orientation::Forward);
            let (a_cur_v1, a_cur_v2) = top_exp_vertices_raw(&edg);
            let mut a_cur_v1 = a_cur_v1.unwrap_or_else(Shape::null);
            let mut a_cur_v2 = a_cur_v2.unwrap_or_else(Shape::null);
            let a_par1 = brep_tool_parameter(&a_cur_v1, &edg);
            let a_par2 = brep_tool_parameter(&a_cur_v2, &edg);

            // OCCT L185-192.
            let mut e1 = self.empty_copied(&edg);
            let mut e2 = self.empty_copied(&edg);
            e1.orientation = Orientation::Forward;
            e2.orientation = Orientation::Forward;
            let mut new_vtx = the_v.clone();
            let a_tol_split_v = brep_tool_tolerance(the_v);

            a_cur_v1 = oriented(&a_cur_v1, Orientation::Forward);

            // OCCT L198-205.
            let a_tol_v1 = if brep_tool_degenerated(&edg) {
                brep_tool_tolerance(&a_cur_v1).max(a_tol_split_v)
            } else {
                brep_tool_tolerance(&a_cur_v1)
            };
            add_vertex_and_update(&mut e1, &mut a_cur_v1, a_par1, a_tol_v1);

            // OCCT L206-211.
            new_vtx = oriented(&new_vtx, Orientation::Reversed);
            add_vertex_and_update(&mut e1, &mut new_vtx, the_p, brep_tool_tolerance(the_v));
            new_vtx = oriented(&new_vtx, Orientation::Forward);
            add_vertex_and_update(&mut e2, &mut new_vtx, the_p, brep_tool_tolerance(the_v));

            // OCCT L213-219.
            a_cur_v2 = oriented(&a_cur_v2, Orientation::Reversed);
            let a_tol_v2 = if brep_tool_degenerated(&edg) {
                a_tol_v1
            } else {
                brep_tool_tolerance(&a_cur_v2)
            };
            builder_add_edge_vertex(&mut e2, &a_cur_v2);
            builder_update_vertex_parameter(&mut a_cur_v2, a_par2, &mut e2, a_tol_v2);
            builder_add_edge_vertex(&mut e2, &a_cur_v2);

            // OCCT L221-227.
            a_new_list.push(e1);
            a_new_list.push(e2);
            for edg1 in le.iter() {
                a_new_list.push(edg1.clone());
            }
            // OCCT L228-229: myMap.UnBind(E); myMap.Bind(E, aNewList).
            self.map_assign(the_e, a_new_list);
        } else {
            // OCCT L233-247.
            let mut e1 = self.empty_copied(&edg);

            for vtx in explorer(&edg, ShapeType::Vertex, ShapeType::Shape) {
                let mut vtx = vtx;
                let f = brep_tool_parameter(&vtx, &edg);
                builder_add_edge_vertex(&mut e1, &vtx);
                let vtx_tol = brep_tool_tolerance(&vtx);
                builder_update_vertex_parameter(&mut vtx, f, &mut e1, vtx_tol);
                builder_add_edge_vertex(&mut e1, &vtx);
            }
            let mut v = the_v.clone();
            builder_add_edge_vertex(&mut e1, &v);
            builder_update_vertex_parameter(&mut v, the_p, &mut e1, brep_tool_tolerance(the_v));
            builder_add_edge_vertex(&mut e1, &v);
            le.push(e1);
            // OCCT L247: le.Append(E1).
            self.map_assign(the_e, le);
        }
    }

    // -----------------------------------------------------------------------
    // Add(W, F) (cxx L602-654).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::Add(const TopoDS_Wire& W, const TopoDS_Face& F)
    /// (cxx L602-654).
    pub fn add_wire_on_face(&mut self, the_w: &Shape, the_f: &Shape) -> bool {
        // OCCT L605-608.
        if self.my_done {
            return false;
        }

        // OCCT L610-615.
        if self.map_list(the_f).is_empty() {
            self.rebuild(the_f);
        }

        // OCCT L616-642: the try/catch (Standard_Failure) around the two
        // internal calls. The rcad branch keeps the OCCT control flow: a
        // failure inside the internal calls returns false.
        if !loc_ope_closed_wire(the_w, the_f) {
            // OCCT L620-623.
            if !self.add_open_wire(the_w, the_f) {
                return false;
            }
        } else {
            // OCCT L626-631.
            if !self.add_closed_wire(the_w, the_f) {
                return false;
            }
        }

        // JAG 10.11.95 Codage des regularites (OCCT L643-652).
        for edg in explorer(the_w, ShapeType::Edge, ShapeType::Shape) {
            if !has_continuity(&edg, the_f, the_f) {
                let mut edg = edg;
                builder_continuity(&mut edg, the_f, the_f);
            }
        }
        true
    }

    // -----------------------------------------------------------------------
    // Add(Lwires, F) (cxx L256-598).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::Add(const NCollection_List<TopoDS_Shape>&
    /// Lwires, const TopoDS_Face& F) (cxx L256-598).
    pub fn add_wires_on_face(&mut self, the_lwires: &[Shape], the_f: &Shape) -> bool {
        // OCCT L259-262.
        if self.my_done {
            return false;
        }

        // OCCT L264-268.
        if self.map_list(the_f).is_empty() {
            self.rebuild(the_f);
        }

        // On cherche la face descendante de F qui contient le wire
        // (OCCT L270-298).
        let lf = self.map_list(the_f).clone();
        let mut a_l_inside: Vec<Shape> = Vec::new();
        let mut itl_face: Option<Shape> = None;
        for fac in lf.iter() {
            let fac = fac.clone();
            let mut local_inside: Vec<Shape> = Vec::new();
            for a_wire in the_lwires.iter() {
                if is_inside_wire(&fac, a_wire) {
                    local_inside.push(a_wire.clone());
                }
            }
            if !local_inside.is_empty() {
                a_l_inside = local_inside;
                itl_face = Some(fac);
                break;
            }
        }
        // OCCT L295-298.
        if a_l_inside.is_empty() || itl_face.is_none() {
            return false;
        }

        // OCCT L300-302.
        let face_ref_src = itl_face.expect("itl");
        let face_ref = oriented(&face_ref_src, Orientation::Forward);
        // lf.Remove(itl).
        let mut lf = lf;
        if let Some(pos) = lf.iter().position(|x| x.is_same(&face_ref_src)) {
            lf.remove(pos);
        }

        // OCCT L304.
        let mut new_wires: Vec<Shape> = Vec::new();

        // OCCT L306-310: SectionsTimes Bind(aLInside(i), 2).
        let mut sections_times: HashMap<ShapeKey, i32> = HashMap::new();
        for w in a_l_inside.iter() {
            sections_times.insert(shape_key(w), 2);
        }

        // OCCT L312-313.
        let mut break_vertices: Vec<Shape> = Vec::new();
        let mut break_on_wires: Vec<Shape> = Vec::new();

        // OCCT L315: VerWireMap.
        let mut ver_wire_map: HashMap<ShapeKey, Shape> = HashMap::new();

        // OCCT L317-348.
        for itl in a_l_inside.iter() {
            let a_section = itl.clone();
            // TopExp::Vertices(aSection, Ver[0], Ver[1]).
            let (v0, v1) = top_exp_vertices_wire(&a_section);
            let ver = [v0.unwrap_or_else(Shape::null), v1.unwrap_or_else(Shape::null)];
            for i in 0..2 {
                if ver_wire_map.contains_key(&shape_key(&ver[i])) {
                    continue;
                }
                for a_wire in explorer(&face_ref, ShapeType::Wire, ShapeType::Shape) {
                    let mut a_ver = Shape::null();
                    for v in explorer(&a_wire, ShapeType::Vertex, ShapeType::Shape) {
                        a_ver = v;
                        if a_ver.is_same(&ver[i]) {
                            break;
                        }
                    }
                    if a_ver.is_same(&ver[i]) {
                        ver_wire_map.insert(shape_key(&a_ver), a_wire.clone());
                        break;
                    }
                }
            }
        }

        // OCCT L350-368: VerSecMap.
        let mut ver_sec_map: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();
        for itl in a_l_inside.iter() {
            let a_wire = itl.clone();
            let (v1, v2) = top_exp_vertices_wire(&a_wire);
            let v1 = v1.unwrap_or_else(Shape::null);
            let v2 = v2.unwrap_or_else(Shape::null);
            // OCCT L358-362.
            ver_sec_map
                .entry(shape_key(&v1))
                .or_insert_with(Vec::new)
                .push(a_wire.clone());
            // OCCT L363-367.
            ver_sec_map
                .entry(shape_key(&v2))
                .or_insert_with(Vec::new)
                .push(a_wire.clone());
        }

        // OCCT L373-375.
        let outer_w = brep_tools_outer_wire(&face_ref);
        let mut cur_wire = outer_w.clone();
        let mut mw = MakeWire::new();
        let mut wexp = WireExplorer::with_wire_face(&cur_wire, &face_ref);

        // OCCT L378-517: the global for (;;).
        loop {
            // OCCT L380-383.
            let mut the_start_vertex = wexp.current_vertex().clone();
            let mut cur_edge = wexp.current().clone();
            let mut last_edge = cur_edge.clone();
            mw.add_edge(&cur_edge);
            let mut a_section_wire = Shape::null();
            let mut a_break_vertex;
            wexp.next();
            if !wexp.more() {
                // OCCT L387-390: wexp.Init(CurWire, FaceRef).
                wexp.init_wire_face(&cur_wire, &face_ref);
            }

            // OCCT L391-413.
            loop {
                // OCCT L393-396.
                if make_wire_closed(&mw) {
                    break;
                }
                let cur_vertex = wexp.current_vertex().clone();
                if let Some(l) = ver_sec_map.get(&shape_key(&cur_vertex)) {
                    // OCCT L400-402: ChooseDirection(LastEdge, CurVertex,
                    // FaceRef, VerSecMap(CurVertex)).
                    let l = l.clone();
                    a_section_wire = choose_direction(&last_edge, &cur_vertex, &face_ref, &l);
                    break;
                }
                cur_edge = wexp.current().clone();
                mw.add_edge(&cur_edge);
                last_edge = cur_edge;
                wexp.next();
                if !wexp.more() {
                    wexp.init_wire_face(&cur_wire, &face_ref);
                }
            }

            // OCCT L414-428.
            if make_wire_closed(&mw) {
                new_wires.push(mw.wire().clone());
                the_start_vertex = break_vertices.first().expect("BreakVertices").clone();
                break_vertices.remove(0);
                cur_wire = break_on_wires.first().expect("BreakOnWires").clone();
                break_on_wires.remove(0);
                wexp.init_wire_face(&cur_wire, &face_ref);
                while !wexp.current_vertex().is_same(&the_start_vertex) {
                    wexp.next();
                }
                mw = MakeWire::new();
                continue;
            }
            a_break_vertex = wexp.current_vertex().clone();
            break_vertices.push(a_break_vertex.clone());
            break_on_wires.push(cur_wire.clone());

            // OCCT L432-512: the section-wire loop.
            loop {
                // OCCT L434-439.
                mw.add_edge(&a_section_wire);
                if let Some(t) = sections_times.get_mut(&shape_key(&a_section_wire)) {
                    *t -= 1;
                }
                if sections_times.get(&shape_key(&a_section_wire)) == Some(&0) {
                    sections_times.remove(&shape_key(&a_section_wire));
                }
                // OCCT L440-458.
                if make_wire_closed(&mw) {
                    new_wires.push(mw.wire().clone());
                    if sections_times.is_empty() {
                        break;
                    }
                    the_start_vertex = break_vertices.first().expect("BreakVertices").clone();
                    break_vertices.remove(0);
                    cur_wire = break_on_wires.first().expect("BreakOnWires").clone();
                    break_on_wires.remove(0);
                    wexp.init_wire_face(&cur_wire, &face_ref);
                    while !wexp.current_vertex().is_same(&the_start_vertex) {
                        wexp.next();
                    }
                    mw = MakeWire::new();
                    break;
                } else {
                    // OCCT L461-511.
                    let (v1, v2) = top_exp_vertices_wire(&a_section_wire);
                    let v1 = v1.unwrap_or_else(Shape::null);
                    let v2 = v2.unwrap_or_else(Shape::null);
                    let a_start_vertex = if v1.is_same(&a_break_vertex) { v2 } else { v1 };
                    cur_wire = ver_wire_map
                        .get(&shape_key(&a_start_vertex))
                        .cloned()
                        .unwrap_or_else(Shape::null);

                    wexp.init_wire_face(&cur_wire, &face_ref);
                    while !wexp.current_vertex().is_same(&a_start_vertex) {
                        wexp.next();
                    }

                    // OCCT L472-476.
                    let l_sections = ver_sec_map
                        .get(&shape_key(&a_start_vertex))
                        .cloned()
                        .unwrap_or_default();
                    if l_sections.len() == 1 {
                        break;
                    }

                    // else: choose the way (OCCT L478-508).
                    let next_section_wire = if a_section_wire.is_same(&l_sections[0]) {
                        l_sections[l_sections.len() - 1].clone()
                    } else {
                        l_sections[0].clone()
                    };

                    // OCCT L482-491: how many VerWireMap entries point at
                    // CurWire.
                    let mut times = 0usize;
                    for (_k, v) in ver_wire_map.iter() {
                        if v.is_same(&cur_wire) {
                            times += 1;
                        }
                    }
                    if times == 1 {
                        // it is inner touching wire (OCCT L492-495).
                    } else {
                        // we have to choose the direction (OCCT L496-508).
                        let a_start_edge = wexp.current().clone();
                        let mut l_dirs: Vec<Shape> = Vec::new();
                        l_dirs.push(a_start_edge.clone());
                        l_dirs.push(next_section_wire.clone());
                        let the_direction =
                            choose_direction(&a_section_wire, &a_start_vertex, &face_ref, &l_dirs);
                        if the_direction.is_same(&a_start_edge) {
                            break;
                        }
                    }
                    a_section_wire = next_section_wire;
                    a_break_vertex = a_start_vertex;
                }
            }
            // OCCT L513-516.
            if sections_times.is_empty() {
                break;
            }
        }

        // OCCT L519-528: the new faces.
        let mut new_faces: Vec<Shape> = Vec::new();
        for itl in new_wires.iter() {
            let mut a_new_face = self.empty_copied(&face_ref);
            a_new_face.orientation = Orientation::Forward;
            builder_add_face_wire(&mut a_new_face, itl);
            new_faces.push(a_new_face);
        }

        // Inserting holes (OCCT L530-572).
        let mut holes: Vec<Shape> = Vec::new();
        for a_wire in explorer(&face_ref, ShapeType::Wire, ShapeType::Shape) {
            let mut found = false;
            let wire_edges: Vec<Shape> = explorer(&a_wire, ShapeType::Edge, ShapeType::Shape);
            let an_edge = wire_edges.first().cloned().unwrap_or_else(Shape::null);
            for a_new_wire in new_wires.iter() {
                for e in explorer(a_new_wire, ShapeType::Edge, ShapeType::Shape) {
                    if an_edge.is_same(&e) {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }
            if !found {
                holes.push(a_wire.clone());
            }
        }
        for a_hole in holes.iter() {
            for i in 0..new_faces.len() {
                if is_inside_wire(&new_faces[i], a_hole) {
                    builder_add_face_wire(&mut new_faces[i], a_hole);
                    break;
                }
            }
        }

        // Update "myMap" (OCCT L574-582).
        lf.extend(new_faces);
        self.map_assign(the_f, lf);
        for expf in explorer(the_f, ShapeType::Wire, ShapeType::Shape) {
            self.map_clear_list(&expf);
        }

        // JAG 10.11.95 Codage des regularites (OCCT L585-596).
        for itl in a_l_inside.iter() {
            for edg in explorer(itl, ShapeType::Edge, ShapeType::Shape) {
                if !has_continuity(&edg, the_f, the_f) {
                    let mut edg = edg;
                    builder_continuity(&mut edg, the_f, the_f);
                }
            }
        }
        true
    }

    // -----------------------------------------------------------------------
    // AddClosedWire (cxx L658-731).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::AddClosedWire(W, F) (cxx L658-731).
    fn add_closed_wire(&mut self, the_w: &Shape, the_f: &Shape) -> bool {
        // On cherche la face descendante de F qui contient le wire
        // (OCCT L662-678).
        let lf = self.map_list(the_f).clone();
        let mut itl_face: Option<Shape> = None;
        for fac in lf.iter() {
            if is_inside_wire(fac, the_w) {
                itl_face = Some(fac.clone());
                break;
            }
        }
        let Some(face_found) = itl_face else {
            return false;
        };

        // OCCT L682-701.
        let mut or_wire = the_w.orientation;
        let mut new_face = self.empty_copied(the_f);
        new_face.orientation = Orientation::Forward;
        builder_add_face_wire(&mut new_face, the_w);
        let classif = face_classifier(&new_face);
        if classif_perform_infinite_point(&new_face, &classif) == rcad_kernel::topods::State::In {
            // le wire donne defini un trou (OCCT L694-700).
            new_face = self.empty_copied(the_f);
            new_face.orientation = Orientation::Forward;
            or_wire = reverse_orientation(or_wire);
            let w_rev = oriented(the_w, or_wire);
            builder_add_face_wire(&mut new_face, &w_rev);
        }

        // OCCT L703-710.
        let face_ref = oriented(&face_found, Orientation::Forward);
        let mut lf = lf;
        if let Some(pos) = lf.iter().position(|x| x.is_same(&face_found)) {
            lf.remove(pos);
        }
        let mut new_ref = self.empty_copied(&face_ref);
        new_ref.orientation = Orientation::Forward;

        // On suppose que les edges du wire ont des courbes 2d. Comme on ne
        // change pas de surface de base, pas besoin d`UpdateEdge
        // (OCCT L712-726).
        for wir in explorer(&oriented(&face_ref, Orientation::Forward), ShapeType::Wire, ShapeType::Shape)
        {
            if is_inside_two_wires(the_f, &wir, the_w) {
                builder_add_face_wire(&mut new_face, &wir);
            } else {
                builder_add_face_wire(&mut new_ref, &wir);
            }
        }
        // OCCT L727-730.
        let w_or = oriented(the_w, reverse_orientation(or_wire));
        builder_add_face_wire(&mut new_ref, &w_or);
        lf.push(new_ref);
        lf.push(new_face);
        self.map_assign(the_f, lf);
        true
    }

    // -----------------------------------------------------------------------
    // AddOpenWire (cxx L781-1269).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::AddOpenWire(W, F) (cxx L781-1269).
    fn add_open_wire(&mut self, the_w: &Shape, the_f: &Shape) -> bool {
        // On cherche la face descendante de F qui contient le wire
        // (OCCT L784-846).
        let lf = self.map_list(the_f).clone();

        // OCCT L788: BRepTools::Update(F) — architecture difference #4.
        brep_tools_update(the_f);

        // OCCT L792-797.
        let w_forward = oriented(the_w, Orientation::Forward);
        let (vfirst_opt, vlast_opt) = top_exp_vertices_wire(&w_forward);
        let mut vfirst = vfirst_opt.unwrap_or_else(Shape::null);
        let mut vlast = vlast_opt.unwrap_or_else(Shape::null);

        let mut tolf = brep_tool_tolerance(&vfirst);
        let mut toll = brep_tool_tolerance(&vlast);
        let mut tol1 = tolf.max(toll);

        // OCCT L801-846.
        let mut wfirst = Shape::null();
        let mut wlast = Shape::null();
        let mut face_found: Option<Shape> = None;
        for itl in lf.iter() {
            let mut fac = itl.clone();
            if !is_inside_wire(&fac, the_w) {
                continue;
            }

            fac = oriented(&fac, Orientation::Forward);
            let mut ffound = false;
            let mut lfound = false;
            let wires: Vec<Shape> = explorer(&fac, ShapeType::Wire, ShapeType::Shape);
            for wir in wires.iter() {
                let wir = wir.clone();
                for vtx in explorer(&wir, ShapeType::Vertex, ShapeType::Shape) {
                    if !ffound && vtx.is_same(&vfirst) {
                        ffound = true;
                        wfirst = wir.clone();
                    } else if !lfound && vtx.is_same(&vlast) {
                        lfound = true;
                        wlast = wir.clone();
                    }
                    if ffound && lfound {
                        break;
                    }
                }
                if ffound && lfound {
                    break;
                }
            }
            if ffound && lfound {
                face_found = Some(fac);
                break;
            }
        }
        let Some(face_ref_src) = face_found else {
            // OCCT L843-846: `if (!itl.More()) return false;`.
            return false;
        };

        // OCCT L848-851.
        let mut face_ref = oriented(&face_ref_src, Orientation::Forward);
        let mut lf = lf;
        if let Some(pos) = lf.iter().position(|x| x.is_same(&face_ref_src)) {
            lf.remove(pos);
        }

        // OCCT L853-857.
        let bas = BRepAdaptorSurface::new(&face_ref);
        let is_periodic = bas.is_u_periodic() || bas.is_v_periodic();
        tol1 = bas.u_resolution(tol1).max(bas.v_resolution(tol1));

        if wfirst.is_same(&wlast) {
            // on cree 2 faces en remplacement de itl.Value()
            // Essai JAG (OCCT L859-1174).
            let mut wires_first: Vec<Shape> = Vec::new();
            for e in explorer(&wfirst, ShapeType::Edge, ShapeType::Shape) {
                if brep_tool_is_closed_on_surface(&e, &face_ref) {
                    self.my_dbl_e.add(&e);
                }
                wires_first.push(e);
            }

            // OCCT L873-878.
            let mut new_w1 = self.make_wire();
            new_w1.orientation = Orientation::Forward;
            let mut new_w2 = self.make_wire();
            new_w2.orientation = Orientation::Forward;

            // OCCT L880-889.
            let mut nb_e = 0usize;
            for e in explorer(&oriented(the_w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape)
            {
                nb_e += 1;
                let orient = e.orientation;
                wires_first.push(e.clone());
                wires_first.push(oriented(&e, reverse_orientation(orient)));
                self.my_dbl_e.add(&e);
            }

            // OCCT L891-892.
            let mut poss_e = ShapeKeySet::new();
            let mut map_e = ShapeKeySet::new();

            // On recherche l`edge contenant Vlast (OCCT L896-919).
            let mut last_edge = Shape::null();
            for e in explorer(&oriented(the_w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape)
            {
                let mut hit = false;
                for vtx in explorer(&e, ShapeType::Vertex, ShapeType::Shape) {
                    if vtx.is_same(&vfirst) {
                        hit = true;
                        break;
                    }
                }
                if hit {
                    last_edge = oriented(&e, e.orientation);
                    break;
                }
            }

            // OCCT L921-922.
            let a_local_face = oriented(&face_ref, wfirst.orientation);
            let Some((c2d, cf, cl)) = brep_tool_curve_on_surface(&last_edge, &a_local_face) else {
                // OCCT dereferences the null handle; the rcad re-host reports
                // the Standard_NoSuchObject surface.
                panic!("Standard_NoSuchObject");
            };

            // OCCT L924-931.
            let mut pfirst = if last_edge.orientation == Orientation::Forward {
                c2d.point_at(cf)
            } else {
                c2d.point_at(cl)
            };

            // OCCT L933-953.
            for e in explorer(&oriented(the_w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape)
            {
                if nb_e > 1 && e.is_same(&last_edge) {
                    continue;
                }
                let mut hit = false;
                for vtx in explorer(&e, ShapeType::Vertex, ShapeType::Shape) {
                    if vtx.is_same(&vlast) {
                        hit = true;
                        break;
                    }
                }
                if hit {
                    last_edge = oriented(&e, e.orientation);
                    break;
                }
            }
            // OCCT L954-955.
            let a_local_face = oriented(&face_ref, wfirst.orientation);
            let (mut plast, mut dlast) = get_direction(&last_edge, &a_local_face, false);

            // OCCT L957-969.
            let cond = if is_periodic {
                !(vfirst.is_same(&vlast) && same_uv(pfirst, plast, &bas))
            } else {
                !vfirst.is_same(&vlast)
            };

            // OCCT L969-1049.
            while cond {
                // OCCT L971: PossE.Clear().
                poss_e.clear();

                // On enchaine par la fin (OCCT L973-987).
                for edg in wires_first.iter() {
                    let edg = edg.clone();
                    let orient = edg.orientation;
                    let (vdeb, vfin) = top_exp_vertices_raw(&edg);
                    let vdeb = vdeb.unwrap_or_else(Shape::null);
                    let vfin = vfin.unwrap_or_else(Shape::null);

                    if (orient == Orientation::Forward && vlast.is_same(&vdeb))
                        || (orient == Orientation::Reversed && vlast.is_same(&vfin))
                    {
                        poss_e.add(&edg);
                    }
                }
                let nb_poss = poss_e.extent();
                if nb_poss == 0 {
                    break;
                }

                // OCCT L994-1010.
                let mut a_next_edge = Shape::null();
                if nb_poss == 1 {
                    a_next_edge = poss_e.first().expect("PossE").clone();
                    let a_local_face_temp = oriented(&face_ref, wfirst.orientation);
                    let (p, d) = get_direction(&a_next_edge, &a_local_face_temp, false);
                    plast = p;
                    dlast = d;
                } else if nb_poss > 1 {
                    // Faire choix en U,V... (OCCT L1001-1010).
                    let a_local_face_temp = oriented(&face_ref, wfirst.orientation);
                    if !choix_uv(
                        &last_edge,
                        &a_local_face_temp,
                        &poss_e,
                        &mut a_next_edge,
                        &mut plast,
                        &mut dlast,
                    ) {
                        return false;
                    }
                }

                // OCCT L1012-1045.
                if nb_poss >= 1 {
                    if a_next_edge.is_null() {
                        // loop is not closed. Split is not possible
                        return false;
                    }
                    if map_e.contains(&a_next_edge) {
                        break;
                    }
                    builder_add_wire_edge(&mut new_w1, &a_next_edge);
                    map_e.add(&a_next_edge);
                    last_edge = a_next_edge;

                    // OCCT L1028-1035.
                    vlast = if last_edge.orientation == Orientation::Forward {
                        top_exp_last_vertex(&last_edge)
                            .unwrap_or_else(Shape::null)
                    } else {
                        top_exp_first_vertex(&last_edge)
                            .unwrap_or_else(Shape::null)
                    };

                    toll = brep_tool_tolerance(&vlast);
                    tol1 = tolf.max(toll);
                } else {
                    // MODIFICATION PIERRE SMEYERS : si pas de possibilite, on
                    // sort avec erreur (OCCT L1040-1045).
                    return false;
                }

                tol1 = bas.u_resolution(tol1).max(bas.v_resolution(tol1));
            }

            // OCCT L1051-1071: the second face boundary.
            let mut nb_add_bound = 0usize;
            let mut an_e1 = Shape::null();
            let mut an_e2 = Shape::null();
            for edg in wires_first.iter() {
                let edg = edg.clone();
                if !map_e.contains(&edg) {
                    builder_add_wire_edge(&mut new_w2, &edg);
                    map_e.add(&edg);
                    nb_add_bound += 1;
                    if an_e1.is_null() {
                        an_e1 = edg;
                    } else {
                        an_e2 = edg;
                    }
                }
            }
            // check overlapping edges for second face (OCCT L1072-1083).
            if nb_add_bound < 2 {
                return false;
            }
            if nb_add_bound == 2 && !an_e1.is_null() && !an_e2.is_null() {
                if check_overlapping(&an_e1, &an_e2, &face_ref) {
                    return false;
                }
            }

            // OCCT L1085-1116.
            nb_add_bound = 0;
            let mut an_e11 = Shape::null();
            let mut an_e12 = Shape::null();
            for e in sub_shapes_no_cum_ori(&new_w1) {
                if e.shape_type() != ShapeType::Edge {
                    continue;
                }
                nb_add_bound += 1;
                if an_e11.is_null() {
                    an_e11 = e;
                } else {
                    an_e12 = e;
                }
            }
            // check overlapping edges for first face.
            if nb_add_bound < 2 {
                return false;
            }
            if nb_add_bound == 2 && !an_e11.is_null() && !an_e12.is_null() {
                if check_overlapping(&an_e11, &an_e12, &face_ref) {
                    return false;
                }
            }

            // OCCT L1118-1129.
            let mut new_f1 = self.empty_copied(&face_ref);
            new_f1.orientation = Orientation::Forward;
            let mut new_f2 = self.empty_copied(&face_ref);
            new_f2.orientation = Orientation::Forward;
            builder_add_face_wire(&mut new_f1, &new_w1);
            builder_add_face_wire(&mut new_f2, &new_w2);

            // modifs JAG 97.05.28 (OCCT L1130-1151).
            for wir in explorer(&oriented(&face_ref, Orientation::Forward), ShapeType::Wire, ShapeType::Shape)
            {
                if !wir.is_same(&wfirst) {
                    if is_inside_wire(&new_f1, &wir) {
                        builder_add_face_wire(&mut new_f1, &wir);
                    } else if is_inside_wire(&new_f2, &wir) {
                        builder_add_face_wire(&mut new_f2, &wir);
                    }
                    // else: OCCT prints the "ce wire est ni dans newF2 ni dans
                    // newF1" warning (OCCT_DEBUG text, not compiled).
                }
            }
            lf.push(new_f1);
            lf.push(new_f2);

            // Mise a jour des descendants des wires (OCCT L1155-1173).
            for w in explorer(the_f, ShapeType::Wire, ShapeType::Shape) {
                let mut ls = self.map_list(&w).clone();
                if let Some(pos) = ls.iter().position(|x| x.is_same(&wfirst)) {
                    ls.remove(pos);
                    ls.push(new_w1.clone());
                    ls.push(new_w2.clone());
                    self.map_assign(&w, ls);
                }
            }
            self.map_assign(the_f, lf);
        } else {
            // on ne cree qu`une seule face (OCCT L1175-1267).
            let outer_w = brep_tools_outer_wire(&face_ref);
            let mut new_wire = self.make_wire();
            new_wire.orientation = Orientation::Forward;

            // OCCT L1183-1192.
            let or_relat = if wfirst.orientation == wlast.orientation {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };

            // OCCT L1194-1198.
            if wlast.is_same(&outer_w) {
                wlast = wfirst.clone();
                wfirst = outer_w.clone();
            }

            // Edges de wfirst (OCCT L1201-1204).
            for e in explorer(&oriented(&wfirst, Orientation::Forward), ShapeType::Edge, ShapeType::Shape)
            {
                builder_add_wire_edge(&mut new_wire, &e);
            }

            // Edges de wlast (OCCT L1206-1212).
            for edg in explorer(&oriented(&wlast, Orientation::Forward), ShapeType::Edge, ShapeType::Shape)
            {
                let orient = edg.orientation.compose(or_relat);
                let e = oriented(&edg, orient);
                builder_add_wire_edge(&mut new_wire, &e);
            }

            // Edges du wire ajoute, et dans les 2 sens (OCCT L1214-1222).
            for e in explorer(&oriented(the_w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape)
            {
                let orient = e.orientation;
                builder_add_wire_edge(&mut new_wire, &oriented(&e, orient));
                builder_add_wire_edge(&mut new_wire, &oriented(&e, reverse_orientation(orient)));
                self.my_dbl_e.add(&oriented(&e, orient));
            }

            // on refait une face (OCCT L1226-1241).
            let mut new_face = self.empty_copied(&face_ref);
            face_ref = oriented(&face_ref, Orientation::Forward);
            for wir in explorer(&oriented(&face_ref, Orientation::Forward), ShapeType::Wire, ShapeType::Shape)
            {
                if wir.is_same(&wfirst) {
                    let w = oriented(&new_wire, wir.orientation);
                    builder_add_face_wire(&mut new_face, &w);
                } else if !wir.is_same(&wlast) {
                    builder_add_face_wire(&mut new_face, &wir);
                }
            }
            lf.push(new_face);

            // Mise a jour des descendants des wires (OCCT L1244-1266).
            for w in explorer(the_f, ShapeType::Wire, ShapeType::Shape) {
                let mut ls = self.map_list(&w).clone();
                let mut touch = false;
                let mut i = 0;
                while i < ls.len() {
                    if ls[i].is_same(&wfirst) || ls[i].is_same(&wlast) {
                        ls.remove(i);
                        touch = true;
                    } else {
                        i += 1;
                    }
                }
                if touch {
                    ls.push(new_wire.clone());
                    self.map_assign(&w, ls);
                }
            }
            self.map_assign(the_f, lf);
        }
        true
    }

    // -----------------------------------------------------------------------
    // LeftOf (cxx L1273-1333).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::LeftOf(W, F) (cxx L1273-1333).
    pub fn left_of(&mut self, the_w: &Shape, the_f: &Shape) -> Vec<Shape> {
        // OCCT L1276-1279.
        if self.my_shape.is_null() {
            panic!("Standard_NoSuchObject");
        }

        // OCCT L1281-1293.
        let mut the_face = Shape::null();
        for exp in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            if exp.is_same(the_f) {
                the_face = exp;
                break;
            }
        }
        if the_face.is_null() {
            panic!("Standard_NoSuchObject");
        }
        self.my_left.clear();

        // OCCT L1296-1297.
        let the_face = the_face;
        let or_face = the_face.orientation;

        // OCCT L1300-1331.
        for edg in explorer(the_w, ShapeType::Edge, ShapeType::Shape) {
            let map_list = self.map_list(&the_face).clone();
            for itl in map_list.iter() {
                let fac = oriented(itl, or_face);
                let mut face_found = false;
                for edgbis in explorer(&fac, ShapeType::Edge, ShapeType::Shape) {
                    if edgbis.is_same(&edg) && edgbis.orientation == edg.orientation {
                        // la face n`est pas deja presente (OCCT L1312-1322).
                        if !self.my_left.iter().any(|x| x.is_same(&fac)) {
                            self.my_left.push(fac.clone());
                        }
                        face_found = true;
                        break;
                    }
                }
                if face_found {
                    break;
                }
            }
        }
        self.my_left.clone()
    }

    // -----------------------------------------------------------------------
    // DescendantShapes (cxx L1337-1351) / Put (cxx L1355-1373).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::DescendantShapes(S) (cxx L1337-1351).
    pub fn descendant_shapes(&mut self, the_s: &Shape) -> Vec<Shape> {
        if !self.my_done {
            let my_shape = self.my_shape.clone();
            self.rebuild(&my_shape);
            self.my_done = true;
        }
        self.map_list(the_s).clone()
    }

    /// OCCT LocOpe_SplitShape::Put(S) (cxx L1355-1373).
    fn put(&mut self, the_s: &Shape) {
        if !self.my_map.contains_key(&shape_key(the_s)) {
            self.my_map
                .insert(shape_key(the_s), (the_s.clone(), Vec::new()));
            if the_s.shape_type() != ShapeType::Vertex {
                for the_iterator in sub_shapes(the_s) {
                    self.put(&the_iterator);
                }
            } else {
                self.map_append(the_s, the_s);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Rebuild (cxx L1411-1455).
    // -----------------------------------------------------------------------

    /// OCCT LocOpe_SplitShape::Rebuild(S) (cxx L1411-1455).
    fn rebuild(&mut self, the_s: &Shape) -> bool {
        // OCCT L1414-1417.
        if the_s.shape_type() == ShapeType::Face {
            self.update_toleraces(the_s);
        }

        // OCCT L1418-1422.
        if !self.map_list(the_s).is_empty() {
            return !self.map_list(the_s)[0].is_same(the_s);
        }

        // OCCT L1423-1428.
        let mut rebuild = false;
        for it in sub_shapes(the_s) {
            rebuild = self.rebuild(&it) || rebuild;
        }

        // OCCT L1430-1454.
        if rebuild {
            let mut result = self.empty_copied(the_s);
            for it in sub_shapes(the_s) {
                let orient = it.orientation;
                let l = self.map_list(&it).clone();
                for itr in l.iter() {
                    let e = oriented(itr, orient);
                    builder_add_shape(&mut result, &e);
                }
            }
            // Assign "Closed" flag for Wires and Shells only (OCCT L1443-1447).
            if result.shape_type() == ShapeType::Wire || result.shape_type() == ShapeType::Shell {
                let is_closed = brep_tool_is_closed_shape(&result);
                builder_set_closed(&mut result, is_closed);
            }
            self.map_append(the_s, &result);
        } else {
            self.map_append(the_s, the_s);
        }
        rebuild
    }

    /// OCCT static updateToleraces(theFace, theMap) (cxx L1375-1407).
    fn update_toleraces(&mut self, the_face: &Shape) {
        for e in explorer(the_face, ShapeType::Edge, ShapeType::Shape) {
            // OCCT L1383-1385.
            if !self.my_map.contains_key(&shape_key(&e)) {
                continue;
            }
            // OCCT L1387-1391.
            let l_edges = self.map_list(&e).clone();
            if l_edges.len() <= 1 {
                continue;
            }

            // OCCT L1393-1405.
            for itr_e in l_edges.iter() {
                let a_cur_e = itr_e.clone();
                let mut amaxdev = 0.0;
                if shape_analysis_edge_check_same_parameter(&a_cur_e, the_face, &mut amaxdev) {
                    let mut e = a_cur_e;
                    builder_update_edge_tol(&mut e, amaxdev);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Static helpers (cxx L735-777 / L1459-1776).
// ---------------------------------------------------------------------------

/// OCCT static checkOverlapping(theEdge1, theEdge2, theFace)
/// (LocOpe_SplitShape.cxx L735-777).
fn check_overlapping(the_edge1: &Shape, the_edge2: &Shape, the_face: &Shape) -> bool {
    // OCCT L740: BRepAdaptor_Surface anAdS(theFace, false).
    let an_ad_s = BRepAdaptorSurface::new(the_face);

    // OCCT L742-744.
    let max_tol = brep_tool_tolerance(the_edge1) + brep_tool_tolerance(the_edge2);
    let a_max_tol2d = an_ad_s.u_resolution(max_tol).max(an_ad_s.v_resolution(max_tol));
    let a_tol_ang = std::f64::consts::PI / 180.0;

    // OCCT L749-750.
    let Some((a_crv1, a_f1, a_l1)) = brep_tool_curve_on_surface(the_edge1, the_face) else {
        return false;
    };
    let Some((a_crv2, a_f2, a_l2)) = brep_tool_curve_on_surface(the_edge2, the_face) else {
        return false;
    };

    // OCCT L751-775.
    let nb_p = 4;
    let a_step = (a_l2 - a_f2) / (nb_p as f64);
    for j in 1..nb_p {
        let par2 = a_f2 + a_step * (j as f64);
        let a_p2d = a_crv2.point_at(par2);
        let a_v2 = a_crv2.derivative_at(par2);

        // proj.Init(aP2d, aCrv1, aF1, aL1) (arch. diff. #7).
        let proj = ExtPC2d::new(a_p2d, &a_crv1, CONFUSION, a_f1, a_l1);
        // check intermediate points (OCCT L762-765).
        if proj.nb_ext() == 0 || proj.square_distance(1).sqrt() > a_max_tol2d {
            return false;
        }
        let par1 = proj.point(1).param;
        let a_p2d1 = a_crv1.point_at(par1);
        let _ = a_p2d1;
        let a_v1 = a_crv1.derivative_at(par1);

        // OCCT L771-774: aV1.IsParallel(aV2, aTolAng).
        if !vec_is_parallel(a_v1, a_v2, a_tol_ang) {
            return false;
        }
    }
    true
}

/// OCCT gp_Vec2d::IsParallel(V, AngularTolerance).
fn vec_is_parallel(the_a: DVec2, the_b: DVec2, the_ang_tol: f64) -> bool {
    let la = the_a.length();
    let lb = the_b.length();
    if la <= rcad_kernel::math::gp::GP_RESOLUTION || lb <= rcad_kernel::math::gp::GP_RESOLUTION {
        return false;
    }
    (the_a.x * the_b.y - the_a.y * the_b.x).abs() <= the_ang_tol * la * lb
}

/// OCCT static IsInside(F, W1, W2) (cxx L1459-1512).
fn is_inside_two_wires(the_f: &Shape, the_w1: &Shape, the_w2: &Shape) -> bool {
    // Attention, c`est tres boeuf !!!! (OCCT L1460-1484).
    let mut new_face = empty_copied(the_f);
    new_face.orientation = Orientation::Forward;
    builder_add_face_wire(&mut new_face, the_w2);
    let classif = face_classifier(&new_face);
    let mut reversed = false;
    if classif_perform_infinite_point(&new_face, &classif) == rcad_kernel::topods::State::In {
        // le wire donne defini un trou
        reversed = true;
    }

    // OCCT L1487-1501.
    let wire_edges = explorer(the_w1, ShapeType::Edge, ShapeType::Shape);
    let edg = wire_edges.first().cloned().unwrap_or_else(Shape::null);
    let vtx = explorer(&edg, ShapeType::Vertex, ShapeType::Shape)
        .first()
        .cloned()
        .unwrap_or_else(Shape::null);
    let prm = brep_tool_parameter(&vtx, &edg);
    let Some((c2d, _f, _l)) = brep_tool_curve_on_surface(&edg, the_f) else {
        // OCCT L1494-1500: "Edge is not on surface".
        return false;
    };
    let pt2d = c2d.point_at(prm);
    // OCCT L1504-1511.
    if !reversed {
        classif_perform_p(&new_face, pt2d) == rcad_kernel::topods::State::In
    } else {
        classif_perform_p(&new_face, pt2d) == rcad_kernel::topods::State::Out
    }
}

/// OCCT static IsInside(F, W) (cxx L1516-1591).
fn is_inside_wire(the_f: &Shape, the_w: &Shape) -> bool {
    // Attention, c`est tres boeuf !!!! (OCCT L1518-1519).
    for edg in explorer(the_w, ShapeType::Edge, ShapeType::Shape) {
        let Some((c2d, f, l)) = brep_tool_curve_on_surface(&edg, the_f) else {
            continue;
        };
        // OCCT L1528-1546.
        let mut prm;
        if !is_negative_infinite(f) && !is_positive_infinite(l) {
            prm = (f + l) / 2.0;
        } else if is_negative_infinite(f) && is_positive_infinite(l) {
            prm = 0.0;
        } else if is_negative_infinite(f) {
            prm = l - 1.0;
        } else {
            prm = f + 1.0;
        }

        let classif = face_classifier(the_f);
        let mut pt2d = c2d.point_at(prm);
        let mut stat = classif_perform_p(the_f, pt2d);
        // OCCT L1554-1557.
        if stat == rcad_kernel::topods::State::Out {
            return false;
        }

        // OCCT L1559-1588.
        if stat == rcad_kernel::topods::State::On {
            let nb_pnt = 10;
            let mut nb_out = 0;
            let mut nb_in = 0;
            let mut nb_on = 0;
            for j in 1..=nb_pnt {
                // check neighbouring point
                prm = f + (l - f) / (nb_pnt as f64) * ((j - 1) as f64);
                pt2d = c2d.point_at(prm);
                stat = classif_perform_p(the_f, pt2d);
                if stat == rcad_kernel::topods::State::Out {
                    nb_out += 1;
                } else if stat == rcad_kernel::topods::State::In {
                    nb_in += 1;
                } else {
                    nb_on += 1;
                }
            }
            if nb_out > nb_in + nb_on {
                return false;
            }
        }
    }
    true
}

/// OCCT static GetDirection(theEdge, theFace, thePnt, theDir, isFirstEnd)
/// (cxx L1595-1627).
fn get_direction(
    the_edge: &Shape,
    the_face: &Shape,
    is_first_end: bool,
) -> (DVec2, DVec2) {
    let Some((a_c2d, a_first, a_last)) = brep_tool_curve_on_surface(the_edge, the_face) else {
        panic!("Standard_NoSuchObject");
    };

    // OCCT L1604-1607.
    let an_or = the_edge.orientation;
    let take_first = (an_or == Orientation::Forward && is_first_end)
        || (an_or == Orientation::Reversed && !is_first_end);

    // OCCT L1609-1626.
    let dpar = (a_last - a_first) * 0.01;
    let (the_pnt, mut the_dir);
    if take_first {
        the_pnt = a_c2d.point_at(a_first);
        let a_next_pnt = a_c2d.point_at(a_first + dpar);
        the_dir = a_next_pnt - the_pnt;
    } else {
        the_pnt = a_c2d.point_at(a_last);
        let a_prev_pnt = a_c2d.point_at(a_last - dpar);
        the_dir = the_pnt - a_prev_pnt;
    }
    if an_or == Orientation::Reversed {
        the_dir = -the_dir;
    }
    (the_pnt, the_dir)
}

/// OCCT ChoixUV(Last, F, Poss, theResEdge, plst, dlst) (cxx L1631-1687).
fn choix_uv(
    last: &Shape,
    the_f: &Shape,
    poss: &ShapeKeySet,
    the_res_edge: &mut Shape,
    plst: &mut DVec2,
    dlst: &mut DVec2,
) -> bool {
    // OCCT L1638-1645.
    let surf = BRepAdaptorSurface::new(the_f);
    let a_plst = surf.d0(plst.x, plst.y);
    let _ = a_plst;
    let ref2d = *dlst;

    // OCCT L1649-1650.
    let mut imin = 0usize;
    let mut angmax = -std::f64::consts::PI;

    // OCCT L1652-1678.
    for index in 1..=poss.extent() {
        let an_edge = poss.find_key(index).expect("Poss").clone();

        let (p2d, v2d) = get_direction(&an_edge, the_f, true);
        if !same_uv(*plst, p2d, &surf) {
            continue;
        }

        let a_pcur = surf.d0(p2d.x, p2d.y);
        let _ = a_pcur;

        let ang = if !last.is_same(&an_edge) {
            vec2_angle(ref2d, v2d)
        } else {
            -std::f64::consts::PI
        };

        if ang > angmax {
            imin = index;
            angmax = ang;
        }
    }

    // OCCT L1680-1686.
    if imin != 0 {
        *the_res_edge = poss.find_key(imin).expect("Poss").clone();
        let (p, d) = get_direction(the_res_edge, the_f, false);
        *plst = p;
        *dlst = d;
    }
    imin != 0
}

/// OCCT static ChooseDirection(RefDir, RefVertex, theFace, Ldirs)
/// (cxx L1691-1776).
fn choose_direction(
    ref_dir: &Shape,
    ref_vertex: &Shape,
    the_face: &Shape,
    ldirs: &[Shape],
) -> Shape {
    // OCCT L1696-1714.
    let mut ref_edge = Shape::null();
    let mut an_or = Orientation::Forward;
    for e in explorer(ref_dir, ShapeType::Edge, ShapeType::Shape) {
        ref_edge = e.clone();
        let (v1, v2) = top_exp_vertices_raw(&ref_edge);
        let v1 = v1.unwrap_or_else(Shape::null);
        let v2 = v2.unwrap_or_else(Shape::null);
        if v1.is_same(ref_vertex) {
            an_or = Orientation::Reversed;
            break;
        } else if v2.is_same(ref_vertex) {
            an_or = Orientation::Forward;
            break;
        }
    }

    // OCCT L1716-1728.
    let Some((ref_curve, ref_first, ref_last)) = brep_tool_curve_on_surface(&ref_edge, the_face)
    else {
        panic!("Standard_NoSuchObject");
    };
    let ref_par = if an_or == Orientation::Forward {
        ref_last
    } else {
        ref_first
    };
    let (ref_pnt, ref_vec0) = (ref_curve.point_at(ref_par), ref_curve.derivative_at(ref_par));
    let mut ref_vec = ref_vec0;
    if an_or == Orientation::Forward {
        ref_vec = -ref_vec;
    }

    // OCCT L1730-1773.
    let mut min_angle = rcad_kernel::precision::REAL_LAST;
    let mut target_dir = Shape::null();
    for a_shape in ldirs.iter() {
        let mut an_edge = Shape::null();
        let mut or_now = an_or;
        for e in explorer(a_shape, ShapeType::Edge, ShapeType::Shape) {
            an_edge = e.clone();
            let (v1, v2) = top_exp_vertices_raw(&an_edge);
            let v1 = v1.unwrap_or_else(Shape::null);
            let v2 = v2.unwrap_or_else(Shape::null);
            if v1.is_same(ref_vertex) {
                or_now = Orientation::Forward;
                break;
            } else if v2.is_same(ref_vertex) {
                or_now = Orientation::Reversed;
                break;
            }
        }
        an_or = or_now;
        let Some((a_curve, a_first, a_last)) = brep_tool_curve_on_surface(&an_edge, the_face)
        else {
            continue;
        };
        let a_par = if an_or == Orientation::Forward {
            a_first
        } else {
            a_last
        };
        let (_, a_vec0) = (a_curve.point_at(a_par), a_curve.derivative_at(a_par));
        let mut a_vec = a_vec0;
        if an_or == Orientation::Reversed {
            a_vec = -a_vec;
        }
        let mut an_angle = angle_between(a_vec, ref_vec);
        if an_angle < 0.0 {
            an_angle += 2.0 * std::f64::consts::PI;
        }

        if an_angle < min_angle {
            min_angle = an_angle;
            target_dir = a_shape.clone();
        }
    }

    target_dir
}

/// OCCT gp_Vec2d::Angle(gp_Vec2d) (gp_Vec2d.cxx) — the signed angle in
/// [-PI, PI] (returns 0 when either vector is null, the OCCT
/// `ElCLib::InPeriod`-protected path).
fn angle_between(the_a: DVec2, the_b: DVec2) -> f64 {
    if the_a.length() <= rcad_kernel::math::gp::GP_RESOLUTION
        || the_b.length() <= rcad_kernel::math::gp::GP_RESOLUTION
    {
        return 0.0;
    }
    vec2_angle(the_a, the_b)
}

/// OCCT inline SameUV(P1, P2, theBAS) (LocOpe_SplitShape.hxx L75-90).
fn same_uv(p1: DVec2, p2: DVec2, the_bas: &BRepAdaptorSurface) -> bool {
    let mut is_same = true;
    if the_bas.is_u_periodic() {
        is_same = (p1.x - p2.x).abs() < the_bas.u_period() * 0.5;
    }
    if the_bas.is_v_periodic() {
        is_same = is_same && ((p1.y - p2.y).abs() < the_bas.v_period() * 0.5);
    }
    is_same
}
