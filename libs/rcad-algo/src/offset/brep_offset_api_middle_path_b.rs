//! OCCT BRepOffsetAPI_MiddlePath — Build.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MiddlePath.cxx L265-1010 (Build).
//! The statics + class body live in brep_offset_api_middle_path.rs.
//!
//! Architecture differences:
//! 1. NCollection_Sequence::Iterator Nullify forms -> the Vec entries set
//!    to Shape::null() (the IsNull reads stay the OCCT branch form); the
//!    OCCT `j = 1` counter reset of cxx L712 is the restart flag of the
//!    while-loop.
//! 2. TopExp::Vertices(E, V1, V2, CumOri=true) -> the oriented
//!    top_exp_vertices read; TopExp::MapShapes(F, EDGE, IndexedMap) -> the
//!    ShapeIndexedMap occurrence carrier.
//! 3. The BRepLib_MakeEdge forms -> the local carriers of the class file
//!    (the tool.rs builder style); BRepLib::BuildCurve3d -> the GAP leaf.
//! 4. GC_MakeLine2d(P1, P2) -> the Line2d::new(P, normalize(P2-P1)) form;
//!    gp_Lin/IsParallel/Contains and gp_Vec::AngleWithRef -> the
//!    vector-math re-hosts.
//! 5. BRepLib_MakeWire -> the normal_projection BRepLibMakeWire carrier;
//!    BRepLib_MakeFace(W, OnlyPlane) -> the brep_offset_make_simple_offset
//!    BRepLibMakeFace carrier (IsDone=false keeps the OCCT wire fallback of
//!    the SecFaces construction).

use std::collections::HashMap;

use rcad_kernel::geom::{Circle3, Curve2d, Curve3, Line2d};
use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use glam::DVec2;

use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_curve_on_surface, brep_tool_surface, reversed, sub_shapes,
};
use crate::brep_fill::generator::top_exp_vertices;
use crate::offset::brep_offset_api_middle_path::{
    brep_lib_build_curve3d, geom_lib_inertia, is_valid_edge, tangent_of_edge, type_of_edge,
    BRepGProp, BRepOffsetAPIMiddlePath, GeomAPIInterpolate, GeomAbsCurveType,
};
use rcad_kernel::geom::CurveEval;
use crate::offset::brep_offset_inter2d::BRepToolsWireExplorer;

/// The NCollection_IndexedMap<TopoDS_Shape> occurrence carrier.
struct ShapeIndexedMap {
    items: Vec<Shape>,
}

impl ShapeIndexedMap {
    fn new() -> Self {
        ShapeIndexedMap { items: Vec::new() }
    }

    fn add(&mut self, the_s: &Shape) {
        if !self.items.iter().any(|s| s.is_same(the_s)) {
            self.items.push(the_s.clone());
        }
    }

    fn contains(&self, the_s: &Shape) -> bool {
        self.items.iter().any(|s| s.is_same(the_s))
    }
}

/// OCCT gp_Vec::AngleWithRef(V, Ref) — the signed angle around the
/// reference direction.
fn angle_with_ref(v1: glam::DVec3, v2: glam::DVec3, reference: glam::DVec3) -> f64 {
    let sin = v1.cross(v2).dot(reference);
    let cos = v1.dot(v2);
    sin.atan2(cos)
}

impl BRepOffsetAPIMiddlePath {
    /// OCCT BRepOffsetAPI_MiddlePath::Build(...) (cxx L265-1010).
    #[allow(unused_assignments)]
    pub fn build(&mut self) {
        // OCCT L267-280: the start-section walk.
        let mut start_vertices: Vec<Shape> = Vec::new();
        let mut end_vertices: HashMap<u64, Shape> = HashMap::new();
        let mut end_edges: HashMap<u64, Shape> = HashMap::new();
        let mut sections_edges: Vec<Vec<Shape>> = Vec::new();

        let mut wexp = BRepToolsWireExplorer::new();
        wexp.init(&self.my_start_wire, &Shape::null());
        let mut edge_seq: Vec<Shape> = Vec::new();
        while wexp.more() {
            // OCCT L274: StartVertices.Append(wexp.CurrentVertex()).
            let (_, vl) = top_exp_vertices(&self.my_brep, &wexp.current());
            start_vertices.push(vl);
            edge_seq.push(wexp.current());
            wexp.next();
        }
        if !self.my_closed_section {
            // OCCT L278-279.
            let current = wexp.current();
            let (_, vl) = top_exp_vertices(&self.my_brep, &current);
            start_vertices.push(vl);
        }
        sections_edges.push(edge_seq.clone());

        // OCCT L282-289: the end-section walk.
        wexp.init(&self.my_end_wire, &Shape::null());
        while wexp.more() {
            let current = wexp.current();
            let (_, vl) = top_exp_vertices(&self.my_brep, &current);
            end_vertices.insert(vl.ptr_id(), vl);
            end_edges.insert(current.ptr_id(), current);
            wexp.next();
        }
        if !self.my_closed_section {
            let current = wexp.current();
            let (_, vl) = top_exp_vertices(&self.my_brep, &current);
            end_vertices.insert(vl.ptr_id(), vl);
        }

        // OCCT L291-297: myStartWireEdges / myEndWireEdges.
        for value in sub_shapes(&self.my_start_wire) {
            self.my_start_wire_edges.insert(value.ptr_id(), value);
        }
        for value in sub_shapes(&self.my_end_wire) {
            self.my_end_wire_edges.insert(value.ptr_id(), value);
        }

        // OCCT L299-304: VEmap / EFmap.
        let mut ve_map = indexmap::IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &self.my_initial_shape,
            ShapeType::Vertex,
            ShapeType::Edge,
            &mut ve_map,
        );
        let mut ef_map = indexmap::IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &self.my_initial_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut ef_map,
        );

        // OCCT L306-331: Initialization of myPaths.
        let mut cur_vertices: HashMap<u64, Shape> = HashMap::new();
        {
            let start_vertices_c = start_vertices.clone();
            for a_start_vertex in start_vertices_c.iter() {
                let mut edges: Vec<Shape> = Vec::new();
                let (_, le) = ve_map
                    .get(&(a_start_vertex.ptr_id(), a_start_vertex.location))
                    .cloned()
                    .unwrap_or_else(|| (Shape::null(), Vec::new()));
                for value in le.iter() {
                    let mut an_edge = value.clone();
                    if !self.my_start_wire_edges.contains_key(&an_edge.ptr_id()) {
                        // OCCT L315-324.
                        let (v1, v2) = top_exp_vertices(&self.my_brep, &an_edge);
                        if v1.is_same(a_start_vertex) {
                            cur_vertices.insert(v2.ptr_id(), v2);
                        } else {
                            an_edge = reversed(&an_edge);
                            cur_vertices.insert(v1.ptr_id(), v1);
                        }
                        edges.push(an_edge);
                        break;
                    }
                }
                if !edges.is_empty() {
                    self.my_paths.push(edges);
                } else {
                    // OCCT L330.
                    return;
                }
            }
        }

        // OCCT L334-396: Filling of myPaths.
        let mut next_vertices: Vec<Shape> = Vec::new();
        loop {
            for path in self.my_paths.iter_mut() {
                // OCCT L338-348.
                let the_shape: Shape = path.last().cloned().unwrap_or_else(Shape::null);
                let the_edge: Shape;
                let the_vertex: Shape;
                if the_shape.shape_type() == ShapeType::Edge {
                    the_edge = the_shape.clone();
                    let (_, vl) = top_exp_vertices(&self.my_brep, &the_edge);
                    the_vertex = vl;
                } else {
                    // last segment of path was punctual.
                    let len = path.len();
                    the_edge = path[len - 2].clone();
                    the_vertex = the_shape;
                }

                // OCCT L351-353.
                if end_vertices.contains_key(&the_vertex.ptr_id()) {
                    continue;
                }
                let (_, le) = ve_map
                    .get(&(the_vertex.ptr_id(), the_vertex.location))
                    .cloned()
                    .unwrap_or_else(|| (Shape::null(), Vec::new()));
                // OCCT L356-375.
                let mut next_edge_candidates: HashMap<u64, Shape> = HashMap::new();
                for value in le.iter() {
                    let mut an_edge = value.clone();
                    if an_edge.is_same(&the_edge) {
                        continue;
                    }
                    let (v1, v2) = top_exp_vertices(&self.my_brep, &an_edge);
                    let next_vertex = if v1.is_same(&the_vertex) {
                        v2
                    } else {
                        an_edge = reversed(&an_edge);
                        v1
                    };
                    if !cur_vertices.contains_key(&next_vertex.ptr_id()) {
                        next_edge_candidates.insert(an_edge.ptr_id(), an_edge);
                    }
                }
                // OCCT L376-393.
                if !next_edge_candidates.is_empty() {
                    if next_edge_candidates.len() > 1 {
                        // OCCT L379-380: punctual segment of path.
                        path.push(the_vertex);
                    } else {
                        // OCCT L382-388.
                        let an_edge = next_edge_candidates.values().next().cloned().unwrap();
                        path.push(an_edge.clone());
                        let (_, next_vertex) = top_exp_vertices(&self.my_brep, &an_edge);
                        next_vertices.push(next_vertex);
                    }
                }
            }
            // OCCT L394-395.
            if next_vertices.is_empty() {
                break;
            }
            for value in next_vertices.clone() {
                cur_vertices.insert(value.ptr_id(), value);
            }
            next_vertices.clear();
        }

        // OCCT L436-443: Building of set of sections.
        let nb_e = edge_seq.len();
        let nb_paths = self.my_paths.len();
        let mut nb_ver = self.my_paths.len();
        if self.my_closed_section {
            nb_ver += 1;
        }
        let mut i = 1usize;
        loop {
            // OCCT L440-442: EdgeSeq(j).Nullify().
            for value in edge_seq.iter_mut() {
                *value = Shape::null();
            }

            // OCCT L444.
            let mut to_insert_vertex = false;

            // OCCT L446-762: the main per-vertex loop (the OCCT
            // `j = 1` reset of L712 is the restart flag).
            let mut j = 2usize;
            while j <= nb_ver {
                // OCCT L448-450.
                if !edge_seq[j - 2].is_null() {
                    j += 1;
                    continue;
                }

                // OCCT L454-468: the end-of-initial-shape pads.
                if self.my_paths[j - 2].len() < i {
                    let a_e1 = self.my_paths[j - 2][i - 1].clone();
                    let (_, last_ver) = top_exp_vertices(&self.my_brep, &a_e1);
                    self.my_paths[j - 2].push(last_ver);
                }
                {
                    let idx = if j <= nb_paths { j } else { 1 };
                    if self.my_paths[idx - 1].len() < i {
                        let a_e2 = self.my_paths[idx - 1][i - 1].clone();
                        let (_, last_ver) = top_exp_vertices(&self.my_brep, &a_e2);
                        self.my_paths[idx - 1].push(last_ver);
                    }
                }

                // OCCT L473-488: the vertex insertion.
                if to_insert_vertex {
                    if self.my_paths[j - 2][i - 1].shape_type() == ShapeType::Edge {
                        let a_e1 = self.my_paths[j - 2][i - 1].clone();
                        let (fver, _) = top_exp_vertices(&self.my_brep, &a_e1);
                        self.my_paths[j - 2].insert(i - 1, fver);
                    }
                    {
                        let idx = if j <= nb_paths { j } else { 1 };
                        if self.my_paths[idx - 1][i - 1].shape_type() == ShapeType::Edge {
                            let a_e2 = self.my_paths[idx - 1][i - 1].clone();
                            let (fver, _) = top_exp_vertices(&self.my_brep, &a_e2);
                            self.my_paths[idx - 1].insert(i - 1, fver);
                        }
                    }
                    to_insert_vertex = false;
                }

                // OCCT L490-500: E1 / E2 / E12.
                let mut e1: Shape = Shape::null();
                let mut e2: Shape = Shape::null();
                if self.my_paths[j - 2][i - 1].shape_type() == ShapeType::Edge {
                    e1 = self.my_paths[j - 2][i - 1].clone();
                }
                {
                    let idx = if j <= nb_paths { j } else { 1 };
                    if self.my_paths[idx - 1][i - 1].shape_type() == ShapeType::Edge {
                        e2 = self.my_paths[idx - 1][i - 1].clone();
                    }
                }
                let e12 = sections_edges[i - 1][j - 2].clone();

                // OCCT L502-509: find the face on which (E1 or E2) and E12
                // lie.
                let e1_or_e2 = if e1.is_null() { e2.clone() } else { e1.clone() };
                if e1_or_e2.is_null() {
                    // both E1 and E2 are vertices => proper edge is the edge
                    // of the previous section between them.
                    edge_seq[j - 2] = e12;
                    j += 1;
                    continue;
                }
                let (_, l_f) = ef_map
                    .get(&(e1_or_e2.ptr_id(), e1_or_e2.location))
                    .cloned()
                    .unwrap_or_else(|| (Shape::null(), Vec::new()));
                // OCCT L512-531.
                let mut the_face = Shape::null();
                'outer: for a_face in l_f.iter() {
                    let (_, l_f2) = ef_map
                        .get(&(e12.ptr_id(), e12.location))
                        .cloned()
                        .unwrap_or_else(|| (Shape::null(), Vec::new()));
                    for a_face2 in l_f2.iter() {
                        if a_face.is_same(a_face2) {
                            the_face = a_face.clone();
                            break 'outer;
                        }
                    }
                    if !the_face.is_null() {
                        break;
                    }
                }

                // OCCT L533-536.
                let prev_vertex = if e1.is_null() {
                    self.my_paths[j - 2][i - 1].clone()
                } else {
                    let (_, vl) = top_exp_vertices(&self.my_brep, &e1);
                    vl
                };
                let cur_vertex = if e2.is_null() {
                    let idx = if j <= nb_paths { j } else { 1 };
                    self.my_paths[idx - 1][i - 1].clone()
                } else {
                    let (_, vl) = top_exp_vertices(&self.my_brep, &e2);
                    vl
                };

                // OCCT L538-556: ProperEdge.
                let mut proper_edge = Shape::null();
                {
                    let (_, le) = ve_map
                        .get(&(prev_vertex.ptr_id(), prev_vertex.location))
                        .cloned()
                        .unwrap_or_else(|| (Shape::null(), Vec::new()));
                    let mut edges_of_the_face = ShapeIndexedMap::new();
                    for value in sub_shapes(&the_face) {
                        if value.shape_type() == ShapeType::Edge {
                            edges_of_the_face.add(&value);
                        }
                    }
                    for value in le.iter() {
                        let an_edge = value.clone();
                        let (v1, v2) = top_exp_vertices(&self.my_brep, &an_edge);
                        if ((v1.is_same(&prev_vertex) && v2.is_same(&cur_vertex))
                            || (v1.is_same(&cur_vertex) && v2.is_same(&prev_vertex)))
                            && edges_of_the_face.contains(&an_edge)
                            && !an_edge.is_same(&e1)
                        {
                            proper_edge = an_edge;
                            break;
                        }
                    }
                }

                // OCCT L558-563.
                {
                    let idx2 = if j <= nb_paths { j } else { 1 };
                    if self.my_paths[j - 2][i - 1].shape_type() == ShapeType::Vertex
                        && self.my_paths[idx2 - 1][i - 1].shape_type() == ShapeType::Vertex
                    {
                        edge_seq[j - 2] = proper_edge;
                        j += 1;
                        continue;
                    }
                }

                // OCCT L565-566.
                let prev_prev_ver = if e1.is_null() {
                    prev_vertex.clone()
                } else {
                    let (vf, _) = top_exp_vertices(&self.my_brep, &e1);
                    vf
                };
                let prev_cur_ver = if e2.is_null() {
                    cur_vertex.clone()
                } else {
                    let (vf, _) = top_exp_vertices(&self.my_brep, &e2);
                    vf
                };

                // OCCT L567-757.
                if proper_edge.is_null() {
                    // no connection between these two vertices.
                    // OCCT L570-573.
                    let list_one_face: Vec<Shape> = vec![the_face.clone()];

                    if e1.is_null() || e2.is_null() {
                        // OCCT L575-640.
                        if e1.is_null() {
                            e1 = self.my_paths[j - 2][i - 1].clone();
                        }
                        if e2.is_null() {
                            let idx = if j <= nb_paths { j } else { 1 };
                            e2 = self.my_paths[idx - 1][i - 1].clone();
                        }
                        let Some((pc1, f1p, l1p)) =
                            brep_tool_curve_on_surface(&e1, &the_face)
                        else {
                            j += 1;
                            continue;
                        };
                        let Some((pc2, f2p, l2p)) =
                            brep_tool_curve_on_surface(&e2, &the_face)
                        else {
                            j += 1;
                            continue;
                        };
                        let last_par1 = if e1.orientation == Orientation::Forward {
                            l1p
                        } else {
                            f1p
                        };
                        let last_par2 = if e2.orientation == Orientation::Forward {
                            l2p
                        } else {
                            f2p
                        };
                        use rcad_kernel::geom::Curve2dEval;
                        let first_pnt2d = pc1.point_at(last_par1);
                        let last_pnt2d = pc2.point_at(last_par2);
                        let the_surf = brep_tool_surface(&the_face);
                        // OCCT L600-601: GC_MakeLine2d(FirstPnt2d,
                        // LastPnt2d).
                        let the_line = Line2d::new(
                            first_pnt2d,
                            (last_pnt2d - first_pnt2d).normalize(),
                        );
                        let len_ne = first_pnt2d.distance(last_pnt2d);
                        // OCCT L602-604.
                        let new_edge = Self::brep_lib_make_edge_pcurve(
                            &mut self.my_brep,
                            &Curve2d::Line(the_line),
                            &the_surf.clone().unwrap(),
                            &prev_vertex,
                            &cur_vertex,
                            0.0,
                            len_ne,
                        );
                        brep_lib_build_curve3d(&new_edge);
                        edge_seq[j - 2] = new_edge.clone();
                        // OCCT L606: EFmap.Add(NewEdge, ListOneFace).
                        ef_map_insert(&mut ef_map, &new_edge, &list_one_face);
                    } else {
                        // OCCT L609-757: E1 is edge.
                        let Some((pc1, f1p, l1p)) =
                            brep_tool_curve_on_surface(&e1, &the_face)
                        else {
                            j += 1;
                            continue;
                        };
                        let Some((pc2, f2p, l2p)) =
                            brep_tool_curve_on_surface(&e2, &the_face)
                        else {
                            j += 1;
                            continue;
                        };
                        let (first_par1, last_par1) = if e1.orientation == Orientation::Forward
                        {
                            (f1p, l1p)
                        } else {
                            (l1p, f1p)
                        };
                        let (first_par2, last_par2) = if e2.orientation == Orientation::Forward
                        {
                            (f2p, l2p)
                        } else {
                            (l2p, f2p)
                        };
                        use rcad_kernel::geom::Curve2dEval;
                        let first_pnt2d = pc1.point_at(last_par1);
                        let last_pnt2d = pc2.point_at(last_par2);
                        let the_surf = brep_tool_surface(&the_face);
                        // OCCT L630-637.
                        let the_line = Line2d::new(
                            first_pnt2d,
                            (last_pnt2d - first_pnt2d).normalize(),
                        );
                        let len_ne = first_pnt2d.distance(last_pnt2d);
                        let new_edge = Self::brep_lib_make_edge_pcurve(
                            &mut self.my_brep,
                            &Curve2d::Line(the_line),
                            &the_surf.clone().unwrap(),
                            &prev_vertex,
                            &cur_vertex,
                            0.0,
                            len_ne,
                        );
                        brep_lib_build_curve3d(&new_edge);
                        // OCCT L638-652.
                        let prev_first_pnt2d = pc1.point_at(first_par1);
                        let prev_last_pnt2d = pc2.point_at(first_par2);
                        let line1 = Line2d::new(
                            prev_first_pnt2d,
                            (last_pnt2d - prev_first_pnt2d).normalize(),
                        );
                        let line2 = Line2d::new(
                            first_pnt2d,
                            (prev_last_pnt2d - first_pnt2d).normalize(),
                        );
                        let len_ne1 = prev_first_pnt2d.distance(last_pnt2d);
                        let new_edge1 = Self::brep_lib_make_edge_pcurve(
                            &mut self.my_brep,
                            &Curve2d::Line(line1),
                            &the_surf.clone().unwrap(),
                            &prev_prev_ver,
                            &cur_vertex,
                            0.0,
                            len_ne1,
                        );
                        brep_lib_build_curve3d(&new_edge1);
                        let len_ne2 = first_pnt2d.distance(prev_last_pnt2d);
                        let new_edge2 = Self::brep_lib_make_edge_pcurve(
                            &mut self.my_brep,
                            &Curve2d::Line(line2),
                            &the_surf.clone().unwrap(),
                            &prev_vertex,
                            &prev_cur_ver,
                            0.0,
                            len_ne2,
                        );
                        brep_lib_build_curve3d(&new_edge2);
                        // OCCT L653-658.
                        let good_ne = is_valid_edge(&new_edge, &the_face);
                        let good_ne1 = is_valid_edge(&new_edge1, &the_face);

                        let type_e1 = type_of_edge(&e1);
                        let type_e2 = type_of_edge(&e2);

                        // OCCT L660-694.
                        let mut choose_edge = 0i32;
                        if !good_ne || type_e1 != type_e2 {
                            if type_e1 == type_e2 {
                                // !good_ne
                                if good_ne1 {
                                    choose_edge = 1;
                                } else {
                                    choose_edge = 2;
                                }
                            } else {
                                // types are different
                                if type_e1 == GeomAbsCurveType::Line {
                                    choose_edge = 1;
                                } else if type_e2 == GeomAbsCurveType::Line {
                                    choose_edge = 2;
                                } else {
                                    // to be developed later...
                                }
                            }
                        }

                        // OCCT L696-756.
                        if choose_edge == 0 {
                            edge_seq[j - 2] = new_edge;
                            ef_map_insert(&mut ef_map, &edge_seq[j - 2], &list_one_face);
                        } else if choose_edge == 1 {
                            edge_seq[j - 2] = new_edge1;
                            ef_map_insert(&mut ef_map, &edge_seq[j - 2], &list_one_face);
                            // OCCT L700-708.
                            for k in 1..(j - 1) {
                                edge_seq[k - 1] = Shape::null();
                            }
                            for k in 1..j {
                                let a_last_edge = self.my_paths[k - 1][i - 1].clone();
                                let (vertex_as_edge, _) =
                                    top_exp_vertices(&self.my_brep, &a_last_edge);
                                self.my_paths[k - 1].insert(i - 1, vertex_as_edge);
                            }
                            // OCCT L712: j = 1 — start from beginning (the
                            // while increment restores j = 2).
                            j = 1;
                        } else if choose_edge == 2 {
                            edge_seq[j - 2] = new_edge2;
                            ef_map_insert(&mut ef_map, &edge_seq[j - 2], &list_one_face);
                            // OCCT L721.
                            to_insert_vertex = true;
                        }
                    }
                } else {
                    // OCCT L759-762: connecting edge exists.
                    edge_seq[j - 2] = proper_edge;
                }
                j += 1;
            }

            // OCCT L764: SectionsEdges.Append(EdgeSeq).
            sections_edges.push(edge_seq.clone());

            // OCCT L766-775: check for exit from for(;;).
            let mut nb_end_edges = 0usize;
            for value in edge_seq.iter() {
                if end_edges.contains_key(&value.ptr_id()) {
                    nb_end_edges += 1;
                }
            }
            if nb_end_edges == nb_e {
                break;
            }

            // OCCT L778.
            i += 1;
        }

        // OCCT L781-804: final phase — the section faces.
        let nb_sec_faces = sections_edges.len();
        let mut sec_faces: Vec<Shape> = vec![Shape::null(); nb_sec_faces];
        for i in 1..=nb_sec_faces {
            // OCCT L783-792: BRepLib_MakeWire MW.
            let mut make_wire = crate::brep_algo::normal_projection::BRepLibMakeWire::new();
            for j in 1..=nb_e {
                let an_edge = sections_edges[i - 1][j - 1].clone();
                make_wire.add(&[an_edge]);
            }
            if !self.my_closed_section {
                // OCCT L794-797.
                let wire = make_wire.shape().clone();
                let (v1, v2) = top_exp_vertices(&self.my_brep, &wire);
                let an_edge = Self::brep_lib_make_edge_vertices(&v2, &v1);
                make_wire.add(&[an_edge]);
            }
            // OCCT L798-804: BRepLib_MakeFace MF(aWire, true).
            let a_wire = make_wire.shape().clone();
            let mf =
                super::brep_offset_make_simple_offset::BRepLibMakeFace::from_wire(&a_wire, true);
            if mf.is_done() {
                sec_faces[i - 1] = mf.face();
            } else {
                sec_faces[i - 1] = a_wire;
            }
        }

        // OCCT L806-822: Centers.
        let mut centers: Vec<glam::DVec3> = vec![glam::DVec3::ZERO; nb_sec_faces];
        for i in 1..=nb_sec_faces {
            if sec_faces[i - 1].shape_type() == ShapeType::Face {
                let properties = BRepGProp::surface_properties(&sec_faces[i - 1]);
                centers[i - 1] = properties.centre_of_mass();
            } else {
                // wire.
                let properties = BRepGProp::linear_properties(&sec_faces[i - 1]);
                centers[i - 1] = properties.centre_of_mass();
            }
        }

        // OCCT L824-956: MidEdges.
        let mut mid_edges: Vec<Shape> = vec![Shape::null(); nb_sec_faces.saturating_sub(1)];
        let lin_tol = 1.0e-5;
        let ang_tol = 1.0e-7;
        let mut pnt1 = glam::DVec3::ZERO;
        let mut pnt2 = glam::DVec3::ZERO;
        for i in 1..nb_sec_faces {
            // OCCT L830-852: the mid-edge type.
            let mut type_of_mid_edge = GeomAbsCurveType::OtherCurve;
            for (j, path) in self.my_paths.iter().enumerate() {
                let a_shape = &path[i - 1];
                if a_shape.shape_type() == ShapeType::Vertex {
                    type_of_mid_edge = GeomAbsCurveType::OtherCurve;
                    break;
                }
                let a_type = type_of_edge(a_shape);
                if j == 0 {
                    type_of_mid_edge = a_type;
                } else if a_type != type_of_mid_edge {
                    type_of_mid_edge = GeomAbsCurveType::OtherCurve;
                    break;
                }
            }
            // OCCT L854-856: the line form.
            if type_of_mid_edge == GeomAbsCurveType::Line {
                mid_edges[i - 1] = Self::brep_lib_make_edge_of_points(
                    &mut self.my_brep,
                    &centers[i - 1],
                    &centers[i],
                );
            } else if type_of_mid_edge == GeomAbsCurveType::Circle {
                // OCCT L858-956: the arc consensus.
                let mut the_axis_loc = glam::DVec3::ZERO;
                let mut the_axis_dir = glam::DVec3::Z;
                let mut the_dir1 = glam::DVec3::ZERO;
                let mut the_dir2 = glam::DVec3::ZERO;
                let mut the_angle = 0.0f64;
                let mut the_tangent = glam::DVec3::ZERO;
                let mut similar_arcs = true;
                for (j, path) in self.my_paths.iter().enumerate() {
                    let an_edge = path[i - 1].clone();
                    let Some((a_curve_full, fpar, lpar)) = brep_tool_curve(&an_edge) else {
                        continue;
                    };
                    let a_curve = match a_curve_full {
                        Curve3::Trimmed(tc) => (*tc.curve).clone(),
                        other => other,
                    };
                    pnt1 = a_curve.point_at(fpar);
                    pnt2 = a_curve.point_at(lpar);
                    let Curve3::Circle(a_circ) = &a_curve else {
                        similar_arcs = false;
                        break;
                    };
                    if j == 0 {
                        // OCCT L876-887.
                        the_axis_loc = a_circ.center;
                        the_axis_dir = a_circ.normal;
                        the_dir1 = pnt1 - a_circ.center;
                        the_dir2 = pnt2 - a_circ.center;
                        the_angle = lpar - fpar;
                        let the_param = if an_edge.orientation == Orientation::Forward {
                            fpar
                        } else {
                            lpar
                        };
                        the_tangent = a_curve.tangent_at(the_param);
                        if an_edge.orientation == Orientation::Reversed {
                            the_tangent = -the_tangent;
                        }
                    } else {
                        // OCCT L888-909.
                        let same_axis = (a_circ.center - the_axis_loc).length() <= lin_tol
                            && angle_cmp_parallel(a_circ.normal, the_axis_dir, ang_tol);
                        if !same_axis {
                            similar_arcs = false;
                            break;
                        }
                        let a_dir1 = pnt1 - a_circ.center;
                        let a_dir2 = pnt2 - a_circ.center;
                        let eq11 = angle_cmp_eq(a_dir1, the_dir1, ang_tol);
                        let eq22 = angle_cmp_eq(a_dir2, the_dir2, ang_tol);
                        let eq12 = angle_cmp_eq(a_dir1, the_dir2, ang_tol);
                        let eq21 = angle_cmp_eq(a_dir2, the_dir1, ang_tol);
                        if !((eq11 && eq22) || (eq12 && eq21)) {
                            similar_arcs = false;
                            break;
                        }
                    }
                }
                if similar_arcs {
                    // OCCT L912-920.
                    let axis_loc = the_axis_loc;
                    let axis_dir = the_axis_dir;
                    let parameter = (centers[i - 1] - axis_loc).dot(axis_dir);
                    let the_center_of_circ = axis_loc + parameter * axis_dir;

                    let vec1 = centers[i - 1] - the_center_of_circ;
                    let vec2 = centers[i] - the_center_of_circ;

                    // OCCT L932-940.
                    let mut an_angle = angle_with_ref(vec1, vec2, axis_dir);
                    if an_angle < 0.0 {
                        an_angle += 2.0 * std::f64::consts::PI;
                    }
                    if (an_angle - the_angle).abs() > ang_tol {
                        the_axis_dir = -the_axis_dir;
                    }
                    // OCCT L941-944: GC_MakeCircle(theAx2, Vec1.Magnitude()).
                    let x_ref = vec1.normalize_or_zero();
                    let the_circle = Circle3::new_with_ref_dir(
                        the_center_of_circ,
                        the_axis_dir,
                        vec1.length(),
                        x_ref,
                    );
                    // OCCT L945-953: the tangent orientation fix.
                    let start_tangent = the_circle.x_dir.cross(the_axis_dir);
                    if start_tangent.dot(the_tangent) < 0.0 {
                        the_axis_dir = -the_axis_dir;
                        let the_circle2 = Circle3::new_with_ref_dir(
                            the_center_of_circ,
                            the_axis_dir,
                            vec1.length(),
                            x_ref,
                        );
                        // OCCT L954-955: BRepLib_MakeEdge(aME(theCircle, 0.,
                        // theAngle)).
                        mid_edges[i - 1] = Self::brep_lib_make_edge_curve_range(
                            &mut self.my_brep,
                            &Curve3::Circle(the_circle2),
                            0.0,
                            an_angle,
                        );
                    } else {
                        mid_edges[i - 1] = Self::brep_lib_make_edge_curve_range(
                            &mut self.my_brep,
                            &Curve3::Circle(the_circle),
                            0.0,
                            an_angle,
                        );
                    }
                }
            }
        }

        // OCCT L958-1019: Build missed edges.
        let mut i = 1usize;
        while i < nb_sec_faces {
            if mid_edges[i - 1].is_null() {
                // OCCT L962-969.
                let mut j = i + 1;
                while j < nb_sec_faces {
                    if !mid_edges[j - 1].is_null() {
                        break;
                    }
                    j += 1;
                }
                // OCCT L970-976.
                let the_points: Vec<glam::DVec3> = centers[(i - 1)..j].to_vec();
                let mut the_tangents: Vec<glam::DVec3> = vec![glam::DVec3::ZERO; j - i + 1];
                // OCCT L977-1000.
                for k in i..=j {
                    let mut pnt_seq: Vec<glam::DVec3> = Vec::new();
                    for indp in 0..self.my_paths.len() {
                        let a_tangent: glam::DVec3;
                        if k == i {
                            if self.my_paths[indp][k - 1].shape_type() == ShapeType::Vertex {
                                continue;
                            }
                            // at begin.
                            a_tangent = tangent_of_edge(&self.my_paths[indp][k - 1], true);
                        } else if k == j {
                            if self.my_paths[indp][k - 2].shape_type() == ShapeType::Vertex {
                                continue;
                            }
                            // at end.
                            a_tangent = tangent_of_edge(&self.my_paths[indp][k - 2], false);
                        } else {
                            if self.my_paths[indp][k - 2].shape_type() == ShapeType::Vertex
                                || self.my_paths[indp][k - 1].shape_type() == ShapeType::Vertex
                            {
                                continue;
                            }
                            let tangent1 = tangent_of_edge(&self.my_paths[indp][k - 2], false);
                            let tangent2 = tangent_of_edge(&self.my_paths[indp][k - 1], true);
                            a_tangent = tangent1 + tangent2;
                        }
                        let a_tangent = a_tangent.normalize();
                        pnt_seq.push(a_tangent);
                    }
                    // OCCT L987-997: GeomLib::Inertia.
                    let mut the_bary = glam::DVec3::ZERO;
                    let mut xdir = glam::DVec3::ZERO;
                    let mut ydir = glam::DVec3::ZERO;
                    let mut xgap = 0.0;
                    let mut ygap = 0.0;
                    let mut zgap = 0.0;
                    geom_lib_inertia(
                        &pnt_seq,
                        &mut the_bary,
                        &mut xdir,
                        &mut ydir,
                        &mut xgap,
                        &mut ygap,
                        &mut zgap,
                    );
                    the_tangents[k - i] = the_bary;
                }
                // OCCT L1002-1010.
                let the_flags: Vec<bool> = vec![true; j - i + 1];
                let mut interpol = GeomAPIInterpolate::new(&the_points, false, lin_tol);
                interpol.load(&the_tangents, &the_flags);
                interpol.perform();
                // OCCT L1013-1015: the not-done message path (the rcad
                // carrier keeps the IsDone=false exit).
                let _ = interpol.is_done();
                // OCCT L1017-1018.
                let inter_curve = interpol.curve();
                mid_edges[i - 1] = Self::brep_lib_make_edge_curve_range(
                    &mut self.my_brep,
                    &inter_curve,
                    0.0,
                    1.0,
                );
                // OCCT L1019.
                i = j;
            }
            i += 1;
        }

        // OCCT L1022-1029: MakeFinalWire.
        let mut make_final_wire = crate::brep_algo::normal_projection::BRepLibMakeWire::new();
        for i in 1..nb_sec_faces {
            if !mid_edges[i - 1].is_null() {
                make_final_wire.add(&[mid_edges[i - 1].clone()]);
            }
        }

        // OCCT L1031-1032.
        let _final_wire = make_final_wire.shape().clone();
        self.my_shape = make_final_wire.shape().clone();

        // OCCT L1059: Done().
        self.my_done = true;
    }
}

/// OCCT gp_Dir::IsParallel(TheOther, AngularTolerance).
fn angle_cmp_parallel(a: glam::DVec3, b: glam::DVec3, ang_tol: f64) -> bool {
    let dot = a.normalize_or_zero().dot(b.normalize_or_zero());
    dot.abs() >= 1.0 - ang_tol
}

/// OCCT gp_Dir::IsEqual(TheOther, AngularTolerance).
fn angle_cmp_eq(a: glam::DVec3, b: glam::DVec3, ang_tol: f64) -> bool {
    let a = a.normalize_or_zero();
    let b = b.normalize_or_zero();
    let dot = a.dot(b);
    1.0 - dot <= ang_tol
}

/// OCCT EFmap.Add(S, ListOfShape) — the ancestor append.
fn ef_map_insert(
    ef_map: &mut indexmap::IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_s: &Shape,
    the_list: &[Shape],
) {
    let key = (the_s.ptr_id(), the_s.location);
    match ef_map.entry(key) {
        indexmap::map::Entry::Occupied(mut o) => {
            let (_, l) = o.get_mut();
            for f in the_list {
                if !l.iter().any(|x| x.is_same(f)) {
                    l.push(f.clone());
                }
            }
        }
        indexmap::map::Entry::Vacant(v) => {
            v.insert((the_s.clone(), the_list.to_vec()));
        }
    }
}

// The import anchors of the OCCT forms carried by the module header.
const _: fn(&Shape, &Shape) -> bool = is_valid_edge;
const _: fn(&Shape) -> GeomAbsCurveType = type_of_edge;
const _: Option<DVec2> = None;
