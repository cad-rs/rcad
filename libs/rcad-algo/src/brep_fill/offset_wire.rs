//! OCCT BRepFill_OffsetWire (TKBool/BRepFill) — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/BRepFill_OffsetWire.cxx
//! (L88-2940) + BRepFill_OffsetWire.hxx (L48-133).
//!
//! First consumers: BRepOffset_MakeOffset (Stage 2b, the wall / offset-wire
//! path); BRepOffsetAPI_MakeOffset (2e).
//!
//! Reported gaps (plan §0.6, annotated at each call site).  OffsetWire sits
//! on top of the TKMath bisecting-locus stack which is not translated yet;
//! the classes are represented by the placeholder types below so the class
//! structure, member layout and control flow stay 1:1 with OCCT:
//! - BRepMAT2d_Explorer / BRepMAT2d_BisectingLocus / BRepMAT2d_LinkTopoBilo
//!   (TKMath/MAT2d) and MAT_Node / MAT_Arc / MAT_Graph (TKMath/MAT);
//! - Bisector_Bisec (TKMath/Bisector);
//! - BRepFill_TrimEdgeTool (TKBool/BRepFill — own file, not in this stage);
//! - BRepTools_Substitution (TKTopAlgo/BRepTools);
//! - MAT2d_CutCurve (used by CutEdge) and GeomLProp_CLProps2d (CheckBadEdges);
//! - BRepLib_MakeEdge(pcurve, plane) — pcurve-only edge creation used by
//!   KPartCircle / MakeOffset / MakeCircle, and BRepLib::BuildCurve3d(s).
//!
//! Architecture notes:
//! - `TopoDS_Shape` identity in maps (TopTools_ShapeMapHasher == IsSame) maps
//!   to `Shape::ptr_id()`; NCollection_IndexedDataMap maps to a Vec of pairs
//!   (insertion order preserved, OCCT index = position + 1).
//! - `BRep_Tool::Curve / CurveOnSurface` map to the `TEdgeData` fields.
//! - `BRep_Builder::Remove` maps to rebuilding the parent's child list.

use glam::{DVec2, DVec3};
use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Plane, Surface3};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape, ShapeType, TShape};

pub use super::offset_wire_b::{
    BRepFillTrimEdgeTool, BRepMAT2dBisectingLocus, BRepMAT2dExplorer, BRepMAT2dLinkTopoBilo,
    BisectorBisec, MatArc, MatGraph, MatNode, BRepToolsSubstitution,
};
use crate::brep_fill::generator::{shape_key, shape_oriented, shape_reversed, ShapeKey};
use super::offset_wire_b::{
    perform_curve, brep_tool_curve_on_surface, brep_tool_surface, brep_tool_tolerance, builder_remove_from_wire,
    check_bad_edges, check_small_param_on_edge, cut_edge, edge_vertices, explored_children,
    face_wires, is_closed_wire, is_small_closed_edge, k_part_circle, make_edge_vertices,
    make_vertex, set_wire_closed, tshape_type, update_vertex_point, update_vertex_tolerance,
    vertex_point, wire_edges,
};

/// OCCT Precision::Confusion().
pub(super) const TOL_CONFUSION: f64 = CONFUSION;
/// OCCT Precision::PConfusion().
pub(super) const TOL_PCONFUSION: f64 = PCONFUSION;
/// OCCT gp::Resolution().
const GP_RESOLUTION: f64 = 2.2250738585072014e-308;
/// OCCT RealLast().
const REAL_LAST: f64 = f64::MAX;

/// OCCT GeomAbs_JoinType (GeomAbs_JoinType.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsJoinType {
    Arc,
    Intersection,
}

/// OCCT MAT_Side (MAT_Side.hxx) — only MAT_Left is used here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatSide {
    Left,
    Right,
}


// BRepFill_OffsetWire class (BRepFill_OffsetWire.hxx L48-133)
// =============================================================================

/// OCCT BRepFill_OffsetWire — constructs an offset wire to a spine (wire or
/// face).
#[derive(Debug)]
pub struct BRepFillOffsetWire {
    my_spine: Shape,
    my_work_spine: Shape,
    my_offset: f64,
    my_is_open_result: bool,
    my_shape: Shape,
    my_is_done: bool,
    my_join_type: GeomAbsJoinType,
    /// NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>
    /// myMap — insertion-ordered.
    my_map: Vec<(Shape, Vec<Shape>)>,
    my_bilo: BRepMAT2dBisectingLocus,
    my_link: BRepMAT2dLinkTopoBilo,
    /// NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> myMapSpine.
    my_map_spine: HashMap<ShapeKey, Shape>,
    my_call_gen: bool,
}

impl Default for BRepFillOffsetWire {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFillOffsetWire {
    /// OCCT BRepFill_OffsetWire::BRepFill_OffsetWire (L310-314).
    pub fn new() -> Self {
        BRepFillOffsetWire {
            my_spine: Shape::null(),
            my_work_spine: Shape::null(),
            my_offset: 0.0,
            my_is_open_result: false,
            my_shape: Shape::null(),
            my_is_done: false,
            my_join_type: GeomAbsJoinType::Arc,
            my_map: Vec::new(),
            my_bilo: BRepMAT2dBisectingLocus,
            my_link: BRepMAT2dLinkTopoBilo,
            my_map_spine: HashMap::new(),
            my_call_gen: false,
        }
    }

    /// OCCT BRepFill_OffsetWire::BRepFill_OffsetWire(Spine, Join, IsOpenResult)
    /// (L318-323).
    pub fn new_with_spine(brep: &mut BRep, spine: &Shape, join: GeomAbsJoinType, is_open_result: bool) -> Self {
        let mut me = Self::new();
        me.init(brep, spine, join, is_open_result);
        me
    }

    /// OCCT BRepFill_OffsetWire::Init (L327-365).
    pub fn init(&mut self, brep: &mut BRep, spine: &Shape, join: GeomAbsJoinType, is_open_result: bool) {
        self.my_is_done = false;
        self.my_spine = shape_oriented(spine, Orientation::Forward);
        self.my_join_type = join;
        self.my_is_open_result = is_open_result;

        self.my_map.clear();
        self.my_map_spine.clear();

        //------------------------------------------------------------------
        // cut the spine for bissectors.
        //------------------------------------------------------------------
        let mut exp = BRepMAT2dExplorer;
        exp.perform(&self.my_spine);
        self.my_spine = exp.modified_shape(&self.my_spine);
        self.prepare_spine(brep);

        let mut a_shape = Shape::null();
        let mut a_map: Vec<(Shape, Vec<Shape>)> = Vec::new();
        let mut done = false;
        if k_part_circle(
            brep,
            &self.my_work_spine,
            1.0,
            self.my_is_open_result,
            0.0,
            &mut a_shape,
            &mut a_map,
            &mut done,
        ) {
            return;
        }

        //-----------------------------------------------------
        // Calculate the map of bissectors to the left.
        // and Links Topology -> base elements of the map.
        //-----------------------------------------------------
        exp.perform(&self.my_work_spine);
        self.my_bilo
            .compute(&exp, 1, MatSide::Left, self.my_join_type, self.my_is_open_result);
        self.my_link.perform(&exp, &self.my_bilo);
    }

    /// OCCT BRepFill_OffsetWire::IsDone (L369-372).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepFill_OffsetWire::Spine (L376-379).
    pub fn spine(&self) -> &Shape {
        &self.my_spine
    }

    /// OCCT BRepFill_OffsetWire::Shape (L383-386).
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT BRepFill_OffsetWire::GeneratedShapes (L390-443).
    pub fn generated_shapes(&mut self, brep: &BRep, spine_shape: &Shape) -> Vec<Shape> {
        if !self.my_call_gen {
            if !self.my_map_spine.is_empty() {
                // myMapSpine can be empty if passed by PerformWithBilo.
                let pairs: Vec<(ShapeKey, Shape)> = self.my_map_spine.iter().map(|(k, v)| (*k, v.clone())).collect();
                for (k, v) in pairs {
                    let key_shape = find_shape_by_key(brep, k);
                    if self.map_contains(&key_shape) {
                        if !self.map_contains(&v) {
                            self.my_map.push((v.clone(), Vec::new()));
                        }
                        if !v.is_same(&key_shape) {
                            // myMap.ChangeFromKey(it.Value())
                            //   .Append(myMap.ChangeFromKey(it.Key()));
                            // myMap.RemoveKey(it.Key());
                            let from_key = self.map_find_from_key(&key_shape);
                            self.map_change_from_key(&v).extend(from_key);
                            self.map_remove_key(&key_shape);
                        }
                    }
                    if self.map_contains(&shape_reversed(&key_shape)) {
                        if !self.map_contains(&shape_reversed(&v)) {
                            self.my_map.push((shape_reversed(&v), Vec::new()));
                        }
                        if !v.is_same(&key_shape) {
                            let from_key_rev = self.map_find_from_key(&shape_reversed(&key_shape));
                            self.map_change_from_key(&shape_reversed(&v)).extend(from_key_rev);
                            self.map_remove_key(&shape_reversed(&key_shape));
                        }
                    }
                }
            }
            self.my_call_gen = true;
        }

        if self.map_contains(spine_shape) {
            return self.map_find_from_key(spine_shape);
        }
        Vec::new()
    }

    /// OCCT BRepFill_OffsetWire::JoinType (L447-450).
    pub fn join_type(&self) -> GeomAbsJoinType {
        self.my_join_type
    }

    /// OCCT BRepFill_OffsetWire::Generated (L1187-1191).
    pub fn generated(&mut self) -> &mut Vec<(Shape, Vec<Shape>)> {
        &mut self.my_map
    }

    fn map_contains(&self, s: &Shape) -> bool {
        self.my_map.iter().any(|(k, _)| k.is_same(s))
    }

    fn map_find_from_key(&self, s: &Shape) -> Vec<Shape> {
        self.my_map
            .iter()
            .find(|(k, _)| k.is_same(s))
            .map(|(_, v)| v.clone())
            .unwrap_or_default()
    }

    fn map_change_from_key(&mut self, s: &Shape) -> &mut Vec<Shape> {
        let idx = self.my_map.iter().position(|(k, _)| k.is_same(s)).unwrap();
        &mut self.my_map[idx].1
    }

    fn map_remove_key(&mut self, s: &Shape) {
        self.my_map.retain(|(k, _)| !k.is_same(s));
    }

    fn map_add(&mut self, s: &Shape, list: Vec<Shape>) {
        if !self.map_contains(s) {
            self.my_map.push((s.clone(), list));
        }
    }

    /// OCCT BRepFill_OffsetWire::Perform (L454-639).
    pub fn perform(&mut self, brep: &mut BRep, offset: f64, alt: f64) {
        // OCCT wraps the body in try/catch(Standard_Failure); rcad panics
        // propagate (no SEH translation — architecture note).
        self.my_call_gen = false;
        let mut my_is_done = self.my_is_done;
        let mut my_shape = self.my_shape.clone();
        let mut my_map = std::mem::take(&mut self.my_map);
        if k_part_circle(
            brep,
            &self.my_work_spine,
            offset,
            self.my_is_open_result,
            alt,
            &mut my_shape,
            &mut my_map,
            &mut my_is_done,
        ) {
            self.my_shape = my_shape;
            self.my_map = my_map;
            self.my_is_done = my_is_done;
            return;
        }

        let old_work_spain = self.my_work_spine.clone();

        let mut bad_edges: Vec<Shape> = Vec::new();
        // CheckBadEdges needs the locus/link — GAP panics inside.
        check_bad_edges(
            brep,
            &self.my_work_spine,
            offset,
            &self.my_bilo,
            &self.my_link,
            &mut bad_edges,
        );

        if !bad_edges.is_empty() {
            // Modification of myWorkSpine (OCCT L472-595): the BRepTools_
            // Substitution workflow over the bad edges, then a fresh
            // bisecting locus on the substituted spine.
            let mut a_subst = BRepToolsSubstitution;
            let mut a_l: Vec<Shape> = Vec::new();
            let a_defl = 0.01 * offset.abs();
            let mut parameters: Vec<f64> = Vec::new();
            let mut points: Vec<DVec3> = Vec::new();

            for an_e in &bad_edges {
                a_l.clear();
                parameters.clear();
                points.clear();

                let (vf, vl) = crate::brep_fill::generator::top_exp_wire_vertices(brep, an_e);
                let _ = (&vf, &vl);

                // OCCT L493-496: GeomAdaptor_Curve over the edge 3D curve,
                // PerformCurve at the deflection.
                let ed = brep.edge(an_e.clone());
                if let Some(g3d) = &ed.curve {
                    let (f, l) = (ed.range[0], ed.range[1]);
                    perform_curve(
                        &mut parameters,
                        &mut points,
                        g3d,
                        a_defl,
                        f,
                        l,
                        TOL_CONFUSION,
                        2,
                    );
                }

                let npnts = points.len();
                let mut fv = vf.clone();
                if npnts > 2 {
                    for np in 1..npnts - 1 {
                        let lp = points[np];
                        let lv = make_vertex(brep, lp);
                        let new_e = make_edge_vertices(brep, &fv, &lv);
                        a_l.push(new_e);
                        fv = lv;
                    }
                    let new_e = make_edge_vertices(brep, &fv, &vl);
                    a_l.push(new_e);
                } else {
                    let new_e = make_edge_vertices(brep, &vf, &vl);
                    a_l.push(new_e);
                }
                // Update myMapSpine (OCCT L525-544)
                if self.my_map_spine.contains_key(&shape_key(an_e)) {
                    let spine_edge = self.my_map_spine.get(&shape_key(an_e)).cloned().unwrap();
                    for new_e in &a_l {
                        self.my_map_spine.insert(shape_key(new_e), spine_edge.clone());
                        let (nv1, nv2) = edge_vertices(brep, new_e);
                        if !self.my_map_spine.contains_key(&shape_key(&nv1)) {
                            self.my_map_spine.insert(shape_key(&nv1), spine_edge.clone());
                        }
                        if !self.my_map_spine.contains_key(&shape_key(&nv2)) {
                            self.my_map_spine.insert(shape_key(&nv2), spine_edge.clone());
                        }
                    }
                    self.my_map_spine.remove(&shape_key(an_e));
                }
                ///////////////////
                a_subst.substitute(an_e, &a_l);
            }

            // OCCT L549-573: rebuild the wires through the substitution.
            let mut wwmap: Vec<(Shape, Vec<Shape>)> = Vec::new();
            for a_wire in face_wires(brep, &self.my_work_spine) {
                a_subst.build(brep, &a_wire);
                if a_subst.is_copied(&a_wire) {
                    let mut new_wire = a_subst.copy(&a_wire)[0].clone();
                    set_wire_closed(brep, &new_wire, is_closed_wire(brep, &a_wire));
                    wwmap.push((a_wire.clone(), vec![new_wire.clone()]));
                    let _ = &mut new_wire;
                }
            }
            a_subst.clear();
            for (k, v) in &wwmap {
                a_subst.substitute(k, v);
            }

            a_subst.build(brep, &self.my_work_spine);

            if a_subst.is_copied(&self.my_work_spine) {
                self.my_work_spine = a_subst.copy(&self.my_work_spine)[0].clone();

                let mut new_exp = BRepMAT2dExplorer;
                new_exp.perform(&self.my_work_spine);
                let mut new_bilo = BRepMAT2dBisectingLocus;
                let mut new_link = BRepMAT2dLinkTopoBilo;
                new_bilo.compute(&new_exp, 1, MatSide::Left, self.my_join_type, self.my_is_open_result);

                if !new_bilo.is_done() {
                    self.my_shape = Shape::null();
                    self.my_is_done = false;
                    self.my_map = my_map;
                    return;
                }

                new_link.perform(&new_exp, &new_bilo);
                let work = self.my_work_spine.clone();
                let join = self.my_join_type;
                self.perform_with_bilo(brep, &work, offset, &new_bilo, &mut new_link, join, alt);
                self.my_work_spine = old_work_spain;
            } else {
                let work = self.my_work_spine.clone();
                let join = self.my_join_type;
                let mut link = self.my_link.clone();
                let bilo = self.my_bilo.clone();
                self.perform_with_bilo(brep, &work, offset, &bilo, &mut link, join, alt);
            }
        } else {
            let work = self.my_work_spine.clone();
            let join = self.my_join_type;
            let mut link = self.my_link.clone();
            let bilo = self.my_bilo.clone();
            self.perform_with_bilo(brep, &work, offset, &bilo, &mut link, join, alt);
        }
        self.my_shape = my_shape;
        self.my_map = my_map;
        self.my_is_done = my_is_done;
    }

    /// OCCT BRepFill_OffsetWire::PerformWithBiLo (L680-1183) — the offset
    /// body.  The MAT2d-dependent branches hit the GAP panics of the
    /// placeholder types.
    pub fn perform_with_bilo(
        &mut self,
        _brep: &mut BRep,
        _spine: &Shape,
        _offset: f64,
        _locus: &BRepMAT2dBisectingLocus,
        _link: &mut BRepMAT2dLinkTopoBilo,
        _join: GeomAbsJoinType,
        _alt: f64,
    ) {
        panic!("GAP: PerformWithBiLo requires BRepMAT2d_* / Bisector_Bisec / BRepFill_TrimEdgeTool — see file header")
    }

    /// OCCT BRepFill_OffsetWire::PrepareSpine (L1195-1290).
    fn prepare_spine(&mut self, brep: &mut BRep) {
        self.my_work_spine = Shape::null();
        self.my_map_spine.clear();

        // const Geom_Surface& S = BRep_Tool::Surface(mySpine, L);
        // double TolF = BRep_Tool::Tolerance(mySpine);
        // B.MakeFace(myWorkSpine, S, L, TolF);
        let s = brep_tool_surface(brep, &self.my_spine);
        let tol_f = brep_tool_tolerance(brep, &self.my_spine);
        self.my_work_spine = brep.add_tface_tol(
            s,
            Shape::null(),
            vec![],
            None,
            None,
            vec![],
            false,
            tol_f,
        );

        for w in face_wires(brep, &self.my_spine) {
            let nw = brep.add_twire(Vec::new());

            // Modified by Sergey KHROMOV - Thu Nov 16 17:29:55 2000 Begin
            let mut forced_cut = 0i32;
            let mut nb_res_edges: i32 = -1;
            // TopExp::MapShapes(IteF.Value(), TopAbs_EDGE, EdgeMap)
            let edge_map = wire_edges(brep, &w);
            let nb_edges = edge_map.len();

            if nb_edges == 1 && !self.my_is_open_result {
                // in case of open wire there's no need to do it
                forced_cut = 2;
            }
            // Modified by Sergey KHROMOV - Thu Nov 16 17:29:48 2000 End

            for e in &edge_map {
                let (v1, v2) = edge_vertices(brep, e);
                self.my_map_spine.insert(shape_key(&v1), v1.clone());
                self.my_map_spine.insert(shape_key(&v2), v2.clone());
                // Cuts.Clear()
                let mut cuts: Vec<Shape> = Vec::new();

                // TopoDS_Shape aLocalShape = E.Oriented(TopAbs_FORWARD);
                let a_local_shape = shape_oriented(e, Orientation::Forward);
                // Modified by Sergey KHROMOV - Thu Nov 16 17:29:29 2000 Begin
                let a_num_curves_in_edge = {
                    let ed = brep.edge(e.clone());
                    ed.pcurves.len() + ed.representations.len()
                };
                if nb_edges == 2 && nb_res_edges == 0 && a_num_curves_in_edge > 1 {
                    forced_cut = 1;
                }
                // Modified by Sergey KHROMOV - Thu Nov 16 17:29:33 2000 End
                nb_res_edges = cut_edge(brep, &a_local_shape, &self.my_spine, forced_cut, &mut cuts) as i32;

                if cuts.is_empty() {
                    // B.Add(NW, E); myMapSpine.Bind(E, E);
                    append_edge_to_wire(brep, &nw, e);
                    self.my_map_spine.insert(shape_key(e), e.clone());
                } else {
                    for ne in &cuts {
                        let ne = shape_oriented(ne, e.orientation);
                        append_edge_to_wire(brep, &nw, &ne);
                        self.my_map_spine.insert(shape_key(&ne), e.clone());
                        let (nv1, nv2) = edge_vertices(brep, &ne);
                        if !self.my_map_spine.contains_key(&shape_key(&nv1)) {
                            self.my_map_spine.insert(shape_key(&nv1), e.clone());
                        }
                        if !self.my_map_spine.contains_key(&shape_key(&nv2)) {
                            self.map_spine_bind(brep, &nv2, e);
                        }
                    }
                }
            }
            // Modified by Sergey KHROMOV - Thu Mar  7 09:17:41 2002 Begin
            let (a_v1, a_v2) = crate::brep_fill::generator::top_exp_wire_vertices(brep, &nw);
            set_wire_closed(brep, &nw, a_v1.is_same(&a_v2));
            // Modified by Sergey KHROMOV - Thu Mar  7 09:17:43 2002 End
            // B.Add(myWorkSpine, NW)
            let idx = self.my_work_spine.index;
            if let TShape::Face(fd) = &mut *Arc::make_mut(&mut brep.tshapes[idx]) {
                if fd.outer_wire.is_null() {
                    fd.outer_wire = nw.clone();
                } else {
                    fd.inner_wires.push(nw.clone());
                }
                fd.my_shapes.push(nw.clone());
            }
        }
    }

    fn map_spine_bind(&mut self, brep: &BRep, v: &Shape, e: &Shape) {
        let _ = brep;
        self.my_map_spine.insert(shape_key(v), e.clone());
    }

    /// OCCT BRepFill_OffsetWire::UpdateDetromp (L104-113 decl, L1302-1402).
    #[allow(clippy::too_many_arguments)]
    fn update_detromp(
        &self,
        brep: &BRep,
        detromp: &mut HashMap<ShapeKey, Vec<Shape>>,
        detromp_order: &mut Vec<Shape>,
        shape1: &Shape,
        shape2: &Shape,
        vertices: &[Shape],
        params: &[DVec3],
        bisec: &BisectorBisec,
        s_on_e: bool,
        e_on_e: bool,
        trim: &BRepFillTrimEdgeTool,
    ) {
        let mut ii = 1usize;

        fn bind_detromp<'a>(
            brep: &BRep,
            detromp: &'a mut HashMap<ShapeKey, Vec<Shape>>,
            order: &mut Vec<Shape>,
            s: &Shape,
        ) -> &'a mut Vec<Shape> {
            let k = shape_key(s);
            if !detromp.contains_key(&k) {
                detromp.insert(k, Vec::new());
                order.push(s.clone());
            }
            let _ = brep;
            detromp.get_mut(&k).unwrap()
        }
        let _ = bind_detromp;

        if self.my_join_type == GeomAbsJoinType::Intersection {
            while ii <= vertices.len() {
                let a_vertex = vertices[ii - 1].clone();
                bind_detromp(brep, detromp, detromp_order, shape1).push(a_vertex.clone());
                bind_detromp(brep, detromp, detromp_order, shape2).push(a_vertex);
                ii += 1;
            }
        } else {
            // myJoinType == GeomAbs_Arc
            let mut u1;
            let mut v1 = Shape::null();
            let mut v2;

            let bis = bisec.value();
            let mut force_add = false;
            // OCCT: aTC = down_cast<Geom2d_TrimmedCurve>(Bis);
            if let Curve2d::Trimmed(t) = &bis {
                if is_curve2d_periodic(&t.curve) {
                    let pf = bis.point_at(t.t_min);
                    let pl = bis.point_at(t.t_max);
                    force_add = pf.distance(pl) <= TOL_CONFUSION;
                }
            }

            u1 = curve2d_first_parameter(&bis);

            if s_on_e {
                // the first point of the bissectrice is on the offset
                v1 = vertices[ii - 1].clone();
                ii += 1;
            }

            while ii <= vertices.len() && ii <= params.len() {
                let u2 = params[ii - 1].x;
                v2 = vertices[ii - 1].clone();

                let p = bis.point_at((u2 + u1) * 0.5);
                if !trim.is_inside(p) || force_add {
                    if !v1.is_null() {
                        bind_detromp(brep, detromp, detromp_order, shape1).push(v1.clone());
                        bind_detromp(brep, detromp, detromp_order, shape2).push(v1.clone());
                    }
                    bind_detromp(brep, detromp, detromp_order, shape1).push(v2.clone());
                    bind_detromp(brep, detromp, detromp_order, shape2).push(v2.clone());
                }
                u1 = u2;
                v1 = v2;
                ii += 1;
            }

            // test medium point between the last parameter and the end of the bissectrice.
            let u2 = curve2d_last_parameter(&bis);
            if !e_on_e {
                if !is_infinite_value(u2) {
                    let p = bis.point_at((u2 + u1) * 0.5);
                    if !trim.is_inside(p) || force_add {
                        if !v1.is_null() {
                            bind_detromp(brep, detromp, detromp_order, shape1).push(v1.clone());
                            bind_detromp(brep, detromp, detromp_order, shape2).push(v1.clone());
                        }
                    }
                } else if !v1.is_null() {
                    bind_detromp(brep, detromp, detromp_order, shape1).push(v1.clone());
                    bind_detromp(brep, detromp, detromp_order, shape2).push(v1.clone());
                }
            }
        }
    }

    /// OCCT BRepFill_OffsetWire::MakeWires (L1406-1579) — the trimmed offset
    /// edges are chained into wires.
    fn make_wires(&mut self, brep: &mut BRep) {
        //--------------------------------------------------------
        // creation of a single list of created parallel edges.
        //--------------------------------------------------------
        // MVE: NCollection_IndexedDataMap<Vertex, List<Edge>> — ordered.
        let mut mve: Vec<(Shape, Vec<Shape>)> = Vec::new();
        let map = std::mem::take(&mut self.my_map);
        for (_, list) in &map {
            for s in list {
                let e = s.clone();
                let (v1, v2) = edge_vertices(brep, &e);
                if v1.is_same(&v2) && is_small_closed_edge(brep, &e, &v1) {
                    continue; // remove small closed edges
                }
                if !check_small_param_on_edge(brep, &e) {
                    continue;
                }
                for v in [&v1, &v2] {
                    match mve.iter_mut().find(|(kv, _)| kv.is_same(v)) {
                        Some((_, l)) => l.push(e.clone()),
                        None => mve.push((v.clone(), vec![e.clone()])),
                    }
                }
            }
        }

        //--------------------------------------
        // Creation of parallel wires.
        //--------------------------------------
        let the_wires: Vec<Shape> = Vec::new();
        let mut the_wires = the_wires;

        while !mve.is_empty() {
            let nw = brep.add_twire(Vec::new());

            // find the first vertex with a single edge (L1469-1481)
            let mut i = 0usize;
            while i < mve.len() {
                if mve[i].1.len() == 1 {
                    break;
                }
                i += 1;
            }
            if i >= mve.len() {
                i = 0;
            }

            let mut cv = mve[i].0.clone();
            let vf = cv.clone();
            let mut ce = mve[i].1[0].clone();
            let mut end = false;
            // MVE.ChangeFromKey(CV).RemoveFirst();
            {
                let list = &mut mve[i].1;
                list.remove(0);
            }

            if self.my_is_open_result {
                // if (MVE.FindFromKey(CV).IsEmpty()) MVE.RemoveKey(CV);
                if mve[i].1.is_empty() {
                    mve.remove(i);
                }
            }

            // Modified by Sergey KHROMOV - Thu Mar 14 11:29:59 2002 Begin
            let mut is_closed = false;
            // Modified by Sergey KHROMOV - Thu Mar 14 11:30:00 2002 End

            while !end {
                //-------------------------------
                // Construction of a wire.
                //-------------------------------
                let (v1, v2) = edge_vertices(brep, &ce);
                if !cv.is_same(&v1) {
                    cv = v1.clone();
                } else {
                    cv = v2.clone();
                }

                append_edge_to_wire(brep, &nw, &ce);

                if vf.is_same(&cv) || !mve.iter().any(|(kv, _)| kv.is_same(&cv)) {
                    // Modified by Sergey KHROMOV - Thu Mar 14 11:30:14 2002 Begin
                    is_closed = vf.is_same(&cv);
                    // Modified by Sergey KHROMOV - Thu Mar 14 11:30:15 2002 End
                    end = true;
                    // MVE.RemoveKey(VF);
                    mve.retain(|(kv, _)| !kv.is_same(&vf));
                }

                if !end {
                    // if (MVE.FindFromKey(CV).Extent() > 2) — debug print
                    // remove CE from the CV list
                    let civ = mve.iter().position(|(kv, _)| kv.is_same(&cv)).unwrap();
                    if let Some(pos) = mve[civ].1.iter().position(|x| x.is_same(&ce)) {
                        mve[civ].1.remove(pos);
                    }
                    if !mve[civ].1.is_empty() {
                        ce = mve[civ].1[0].clone();
                        mve[civ].1.remove(0);
                    } else if self.my_is_open_result {
                        // CV was a vertex with one edge
                        end = true;
                    }

                    if mve[civ].1.is_empty() {
                        let cvc = mve[civ].0.clone();
                        mve.retain(|(kv, _)| !kv.is_same(&cvc));
                    }
                }
            }
            // Modified by Sergey KHROMOV - Thu Mar 14 11:29:31 2002 Begin
            //     NW.Closed(true);
            set_wire_closed(brep, &nw, is_closed);
            // Modified by Sergey KHROMOV - Thu Mar 14 11:29:37 2002 End
            the_wires.push(nw);
        }

        // update myShape :
        //      -- if only one wire : myShape is a Wire
        //      -- if several wires : myShape is a Compound.
        if the_wires.len() == 1 {
            self.my_shape = the_wires[0].clone();
        } else {
            let r = brep.add_tcompound(the_wires.clone());
            self.my_shape = r;
        }
        self.my_map = map;
    }

    /// OCCT BRepFill_OffsetWire::FixHoles (L1583-1920) — fix holes between
    /// open wires where it is possible.
    fn fix_holes(&mut self, brep: &mut BRep) {
        let mut closed_wires: Vec<Shape> = Vec::new();
        let mut unclosed_wires: Vec<Shape> = Vec::new();
        let mut isolated_wires: Vec<Shape> = Vec::new();

        let mut max_tol = 0.0f64;

        for v in explored_children(brep, &self.my_spine, ShapeType::Vertex) {
            let tol = brep_tool_tolerance(brep, &v);
            if tol > max_tol {
                max_tol = tol;
            }
        }
        max_tol *= 100.0;

        for a_wire in explored_children(brep, &self.my_shape, ShapeType::Wire) {
            // Remove duplicated edges
            let mut eemap: Vec<(Shape, Vec<Shape>)> = Vec::new();
            for an_edge in wire_edges(brep, &a_wire) {
                match eemap.iter_mut().find(|(k, _)| k.is_same(&an_edge)) {
                    Some((_, l)) => l.push(an_edge.clone()),
                    None => eemap.push((an_edge.clone(), Vec::new())),
                }
            }
            // aWire.Free(true) — the Free flag is not consumed by rcad.
            for (_, le) in &eemap {
                for itl in le {
                    // BB.Remove(aWire, itl.Value())
                    builder_remove_from_wire(brep, &a_wire, itl);
                }
            }
            // Sorting
            if is_closed_wire(brep, &a_wire) {
                closed_wires.push(a_wire);
            } else {
                unclosed_wires.push(a_wire);
            }
        }

        while !unclosed_wires.is_empty() {
            let base = unclosed_wires[0].clone();
            let (vf, vl) = crate::brep_fill::generator::top_exp_wire_vertices(brep, &base);
            if vf.is_null() || vl.is_null() {
                panic!("BRepFill_OffsetWire::FixHoles(): Wrong wire.");
            }
            let pf = vertex_point(brep, &vf);
            let pl = vertex_point(brep, &vl);
            let mut dist_f = REAL_LAST;
            let mut dist_l = REAL_LAST;
            let mut index_f = 0usize;
            let mut index_l = 0usize;
            let mut is_first_f = false;
            let mut is_first_l = false;
            for i in 1..unclosed_wires.len() {
                let a_wire = unclosed_wires[i].clone();
                let (v1, v2) = crate::brep_fill::generator::top_exp_wire_vertices(brep, &a_wire);
                if v1.is_null() || v2.is_null() {
                    panic!("BRepFill_OffsetWire::FixHoles(): Wrong wire.");
                }
                let p1 = vertex_point(brep, &v1);
                let p2 = vertex_point(brep, &v2);
                let mut dist = pf.distance(p1);
                if dist < dist_f {
                    dist_f = dist;
                    index_f = i;
                    is_first_f = true;
                }
                dist = pf.distance(p2);
                if dist < dist_f {
                    dist_f = dist;
                    index_f = i;
                    is_first_f = false;
                }
                dist = pl.distance(p1);
                if dist < dist_l {
                    dist_l = dist;
                    index_l = i;
                    is_first_l = true;
                }
                dist = pl.distance(p2);
                if dist < dist_l {
                    dist_l = dist;
                    index_l = i;
                    is_first_l = false;
                }
            }
            if dist_f > max_tol {
                index_f = 0;
            }
            if dist_l > max_tol {
                index_l = 0;
            }
            let mut try_to_close = true;
            if dist_f <= max_tol && dist_l <= max_tol && index_f == index_l && is_first_f == is_first_l {
                if dist_f < dist_l {
                    dist_l = REAL_LAST;
                    index_l += 1;
                } else {
                    dist_f = REAL_LAST;
                    index_f += 1;
                }
                try_to_close = false;
            }
            if dist_f <= max_tol {
                let the_wire = unclosed_wires[index_f].clone();
                let (v1, v2) = crate::brep_fill::generator::top_exp_wire_vertices(brep, &the_wire);
                let vemap = map_shapes_and_ancestors_ve(brep, &the_wire);
                let edge_of = |v: &Shape| -> Shape {
                    vemap
                        .iter()
                        .find(|(kv, _)| kv.is_same(v))
                        .map(|(_, l)| l[0].clone())
                        .unwrap_or_else(Shape::null)
                };
                let the_edge = if is_first_f { edge_of(&v1) } else { edge_of(&v2) };
                for an_edge in wire_edges(brep, &the_wire) {
                    let mut an_edge = an_edge;
                    if is_first_f {
                        an_edge = shape_reversed(&an_edge);
                    }
                    if !an_edge.is_same(&the_edge) {
                        append_edge_to_wire(brep, &base, &an_edge);
                    }
                }
                let the_vertex = if is_first_f { v1 } else { v2 };
                let common_tol = brep_tool_tolerance(brep, &vf).max(brep_tool_tolerance(brep, &the_vertex));
                if dist_f <= common_tol {
                    // theEdge.Free(true);
                    // Vf.Orientation(theVertex.Orientation());
                    // BB.Remove(theEdge, theVertex); BB.Add(theEdge, Vf);
                    // BB.UpdateVertex(Vf, CommonTol);
                    replace_edge_vertex(brep, &the_edge, &the_vertex, &vf);
                    update_vertex_tolerance(brep, &vf, common_tol);
                    let the_edge2 = if is_first_f {
                        shape_reversed(&the_edge)
                    } else {
                        the_edge.clone()
                    };
                    append_edge_to_wire(brep, &base, &the_edge2);
                } else {
                    let the_edge2 = if is_first_f {
                        shape_reversed(&the_edge)
                    } else {
                        the_edge.clone()
                    };
                    append_edge_to_wire(brep, &base, &the_edge2);
                    // Creating new edge from theVertex to Vf
                    let new_edge = make_edge_vertices(brep, &the_vertex, &vf);
                    append_edge_to_wire(brep, &base, &new_edge);
                }
            }
            if dist_l <= max_tol && index_l != index_f {
                let the_wire = unclosed_wires[index_l].clone();
                let (v1, v2) = crate::brep_fill::generator::top_exp_wire_vertices(brep, &the_wire);
                let vemap = map_shapes_and_ancestors_ve(brep, &the_wire);
                let edge_of = |v: &Shape| -> Shape {
                    vemap
                        .iter()
                        .find(|(kv, _)| kv.is_same(v))
                        .map(|(_, l)| l[0].clone())
                        .unwrap_or_else(Shape::null)
                };
                let the_edge = if is_first_l { edge_of(&v1) } else { edge_of(&v2) };
                for an_edge in wire_edges(brep, &the_wire) {
                    let mut an_edge = an_edge;
                    if !is_first_l {
                        an_edge = shape_reversed(&an_edge);
                    }
                    if !an_edge.is_same(&the_edge) {
                        append_edge_to_wire(brep, &base, &an_edge);
                    }
                }
                let the_vertex = if is_first_l { v1 } else { v2 };
                let common_tol = brep_tool_tolerance(brep, &vl).max(brep_tool_tolerance(brep, &the_vertex));
                if dist_l <= common_tol {
                    replace_edge_vertex(brep, &the_edge, &the_vertex, &vl);
                    update_vertex_tolerance(brep, &vl, common_tol);
                    let the_edge2 = if !is_first_l {
                        shape_reversed(&the_edge)
                    } else {
                        the_edge.clone()
                    };
                    append_edge_to_wire(brep, &base, &the_edge2);
                } else {
                    let the_edge2 = if !is_first_l {
                        shape_reversed(&the_edge)
                    } else {
                        the_edge.clone()
                    };
                    append_edge_to_wire(brep, &base, &the_edge2);
                    // Creating new edge from Vl to theVertex
                    let new_edge = make_edge_vertices(brep, &vl, &the_vertex);
                    append_edge_to_wire(brep, &base, &new_edge);
                }
            }
            // Check if it is possible to close resulting wire
            if try_to_close {
                let (vf2, vl2) = crate::brep_fill::generator::top_exp_wire_vertices(brep, &base);
                let common_tol = brep_tool_tolerance(brep, &vf2).max(brep_tool_tolerance(brep, &vl2));
                let vemap = map_shapes_and_ancestors_ve(brep, &base);
                let edge_of = |v: &Shape| -> Shape {
                    vemap
                        .iter()
                        .find(|(kv, _)| kv.is_same(v))
                        .map(|(_, l)| l[0].clone())
                        .unwrap_or_else(Shape::null)
                };
                let efirst = edge_of(&vf2);
                let elast = edge_of(&vl2);
                let pf2 = vertex_point(brep, &vf2);
                let pl2 = vertex_point(brep, &vl2);
                let dist = pf2.distance(pl2);
                if dist <= common_tol {
                    replace_edge_vertex(brep, &elast, &vl2, &vf2);
                    update_vertex_tolerance(brep, &vf2, common_tol);
                    set_wire_closed(brep, &base, true);
                } else if dist <= max_tol {
                    // Creating new edge from Vl to Vf
                    let new_edge = make_edge_vertices(brep, &vf2, &vl2);
                    append_edge_to_wire(brep, &base, &new_edge);
                    set_wire_closed(brep, &base, true);
                }
                let _ = efirst;
            }
            // Updating sequences ClosedWires and UnclosedWires
            if dist_f <= max_tol {
                unclosed_wires.remove(index_f);
            }
            if dist_l <= max_tol && index_l != index_f {
                if dist_f <= max_tol && index_l > index_f {
                    index_l -= 1;
                }
                unclosed_wires.remove(index_l);
            }
            if is_closed_wire(brep, &base) {
                closed_wires.push(base.clone());
                unclosed_wires.remove(0);
            } else if dist_f > max_tol && dist_l > max_tol {
                isolated_wires.push(base.clone());
                unclosed_wires.remove(0);
            }
        }

        // Updating myShape
        if closed_wires.len() + isolated_wires.len() == 1 {
            if !closed_wires.is_empty() {
                self.my_shape = closed_wires[0].clone();
            } else {
                self.my_shape = isolated_wires[0].clone();
            }
        } else {
            let mut all: Vec<Shape> = Vec::new();
            all.extend(closed_wires.iter().cloned());
            all.extend(isolated_wires.iter().cloned());
            let r = brep.add_tcompound(all);
            self.my_shape = r;
        }
    }
}

/// OCCT TopExp::MapShapesAndAncestors(W, VERTEX, EDGE, Map) for FixHoles.
pub(super) fn map_shapes_and_ancestors_ve(brep: &BRep, w: &Shape) -> Vec<(Shape, Vec<Shape>)> {
    let mut map: Vec<(Shape, Vec<Shape>)> = Vec::new();
    for e in wire_edges(brep, w) {
        let ed = brep.edge(e.clone());
        for v in [&ed.first, &ed.last] {
            match map.iter_mut().find(|(kv, _)| kv.is_same(v)) {
                Some((_, l)) => l.push(e.clone()),
                None => map.push((v.clone(), vec![e.clone()])),
            }
        }
    }
    map
}

/// OCCT BB.Add(Wire, Edge) — append the edge to the wire's child list.
pub(super) fn append_edge_to_wire(brep: &mut BRep, w: &Shape, e: &Shape) {
    let idx = w.index;
    if let TShape::Wire(wd) = Arc::make_mut(&mut brep.tshapes[idx]) {
        wd.edges.push(e.clone());
        wd.my_shapes.push(e.clone());
    }
}

/// OCCT BB.Remove(Edge, Vertex) + BB.Add(Edge, V) — replace one endpoint
/// vertex of an edge.
pub(super) fn replace_edge_vertex(brep: &mut BRep, e: &Shape, old_v: &Shape, new_v: &Shape) {
    let ed = brep.edge_mut_inplace(e.clone());
    if ed.first.is_same(old_v) {
        ed.first = new_v.clone();
    }
    if ed.last.is_same(old_v) {
        ed.last = new_v.clone();
    }
    ed.my_shapes = vec![ed.first.clone(), ed.last.clone()];
}

/// OCCT Geom2d_Curve::FirstParameter for the supported rcad kinds.
fn curve2d_first_parameter(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Line(_) => 0.0,
        Curve2d::Circle(_) => 0.0,
        Curve2d::BSpline(b) => b.knots[b.degree],
        Curve2d::Bezier(_) => 0.0,
        Curve2d::Trimmed(t) => t.t_min,
        _ => 0.0,
    }
}

/// OCCT Geom2d_Curve::LastParameter for the supported rcad kinds.
fn curve2d_last_parameter(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Line(_) => 1.0,
        Curve2d::Circle(_) => 2.0 * std::f64::consts::PI,
        Curve2d::BSpline(b) => b.knots[b.knots.len() - b.degree - 1],
        Curve2d::Bezier(_) => 1.0,
        Curve2d::Trimmed(t) => t.t_max,
        _ => 1.0,
    }
}

/// OCCT Precision::IsInfinite.
pub(super) fn is_infinite_value(v: f64) -> bool {
    v >= INFINITE_VALUE || v <= -INFINITE_VALUE
}

/// OCCT Geom2d_Curve::IsPeriodic for the supported rcad kinds.
pub(super) fn is_curve2d_periodic(c: &Curve2d) -> bool {
    match c {
        Curve2d::Circle(_) => true,
        // rcad BSplineCurve2 has no periodic flag; the bisector trims are
        // non-periodic.
        Curve2d::BSpline(_) => false,
        _ => false,
    }
}

/// Recover a shape handle from its map key (OCCT keeps handles directly).
pub(super) fn find_shape_by_key(brep: &BRep, k: ShapeKey) -> Shape {
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if Arc::as_ptr(ts) as u64 == k.0 {
            return brep.shape_at(i);
        }
    }
    Shape::null()
}

/// OCCT Compute (L643-676) — the zero-offset copy of the spine at altitude.
fn compute(
    brep: &mut BRep,
    spine: &Shape,
    a_shape: &mut Shape,
    map: &mut Vec<(Shape, Vec<Shape>)>,
    alt: f64,
) {
    // B.MakeCompound(aShape)
    *a_shape = brep.add_tcompound(Vec::new());
    let mut a_lt = alt;
    if spine.orientation == Orientation::Reversed {
        a_lt = -alt;
    }
    // T.SetTranslation(gp_Vec(0., 0., ALT)) — the rcad wire copy keeps the
    // 3D curves; the translation is applied to the geometry directly
    // (TopLoc_Location mapping, architecture note).
    let tr = |p: DVec3| DVec3::new(p.x, p.y, p.z + a_lt);

    for cur_w in face_wires(brep, spine) {
        // NewW = CurW.Moved(L) — a translated copy of every edge.
        let mut new_edges: Vec<Shape> = Vec::new();
        let old_edges = wire_edges(brep, &cur_w);
        for e in &old_edges {
            let new_e = brep.empty_copied(e);
            {
                let ed = brep.edge_mut_inplace(new_e.clone());
                if let Some(c) = &ed.curve {
                    let nc = match c {
                        Curve3::Line(l) => {
                            let mut l = *l;
                            l.origin = tr(l.origin);
                            Curve3::Line(l)
                        }
                        Curve3::Circle(c) => {
                            let mut c = *c;
                            c.center = tr(c.center);
                            Curve3::Circle(c)
                        }
                        Curve3::BSpline(b) => {
                            let mut b = b.clone();
                            b.control_points = b.control_points.iter().map(|&p| tr(p)).collect();
                            Curve3::BSpline(b)
                        }
                        other => other.clone(),
                    };
                    ed.curve = Some(nc);
                }
            }
            new_edges.push(new_e);
        }
        let new_w = brep.add_twire(new_edges.clone());
        // B.Add(aShape, NewW)
        {
            let idx = a_shape.index;
            if let TShape::Compound(children) = &mut *Arc::make_mut(&mut brep.tshapes[idx]) {
                children.push(new_w.clone());
            }
        }
        // update Map.
        for (old_e, new_e) in old_edges.iter().zip(new_edges.iter()) {
            map.push((old_e.clone(), vec![new_e.clone()]));
        }
        let _ = &new_edges;
    }
    let _ = GP_RESOLUTION;
}
