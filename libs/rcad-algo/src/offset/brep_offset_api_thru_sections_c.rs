//! OCCT BRepOffsetAPI_ThruSections — Generated.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_ThruSections.cxx L1385-1740 (Generated).
//! The class body lives in brep_offset_api_thru_sections.rs; the
//! CreateSmoothed batch in brep_offset_api_thru_sections_b.rs.
//!
//! Architecture differences:
//! 1. NCollection_IndexedDataMap<TopoDS_Shape, List, TopTools_ShapeMapHasher>
//!    -> the IndexMap<(u64, u32), (Shape, Vec<Shape>)> of the
//!    loc_ope_glued_shape TopExp::MapShapesAndAncestors re-host (the
//!    FindFromKey reads go through the (ptr_id, location) key);
//!    NCollection_IndexedMap<TopoDS_Shape> -> the dedup occurrence Vec
//!    (Extent/Add/Contains/j-th forms are kept through the local
//!    ShapeIndexedMap carrier below).
//! 2. BRepTools_WireExplorer::CurrentVertex() -> the oriented last vertex
//!    of the current edge (the wexp traversal form).
//! 3. BRepAdaptor_Surface::GetType() == GeomAbs_Plane -> the Surface3
//!    variant discriminant.
//! 4. TopExp::FirstVertex/LastVertex(E) (raw) -> the thru_sections
//!    first_vertex/last_vertex reads.

use std::collections::HashMap;

use rcad_kernel::geom::Surface3;
use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool::{explorer, sub_shapes};
use crate::brep_fill::generator::top_exp_vertices;
use crate::offset::brep_offset_api_thru_sections::{
    brep_tool_degenerated, first_vertex, last_vertex, BRepOffsetAPIThruSections,
};
use crate::offset::brep_offset_inter2d::BRepToolsWireExplorer;

/// The NCollection_IndexedMap<TopoDS_Shape> occurrence carrier (Extent /
/// Add / Contains / j-th — the insertion-ordered dedup).
struct ShapeIndexedMap {
    items: Vec<Shape>,
}

impl ShapeIndexedMap {
    fn new() -> Self {
        ShapeIndexedMap { items: Vec::new() }
    }

    /// OCCT IndexedMap::Add(S) — the 1-based new index on insertion.
    fn add(&mut self, the_s: &Shape) {
        if !self.contains(the_s) {
            self.items.push(the_s.clone());
        }
    }

    /// OCCT IndexedMap::Contains(S).
    fn contains(&self, the_s: &Shape) -> bool {
        self.items.iter().any(|s| s.is_same(the_s))
    }

    /// OCCT IndexedMap::Extent().
    fn extent(&self) -> usize {
        self.items.len()
    }

    /// OCCT IndexedMap::operator(j) — the 1-based key access.
    fn key(&self, the_j: usize) -> Shape {
        self.items[the_j - 1].clone()
    }
}

impl BRepOffsetAPIThruSections {
    /// OCCT BRepOffsetAPI_ThruSections::Generated(S) (cxx L1385-1740).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L1386: myGenerated.Clear().
        self.my_generated.clear();

        // OCCT L1388-1394: AllFaces.
        let mut all_faces: Vec<Shape> = Vec::new();
        for face in explorer(&self.my_shape.clone(), ShapeType::Face, ShapeType::Shape) {
            all_faces.push(face);
        }

        // OCCT L1396-1512: the EDGE case.
        if s.shape_type() == ShapeType::Edge {
            // OCCT L1398-1401.
            let Some(indices) = self.my_edge_new_indices.get(&s.ptr_id()) else {
                return self.my_generated.clone();
            };

            let indices = indices.clone();
            // OCCT L1403-1411: the faces growing from the first section.
            for ind_of_face in indices.iter() {
                let ind_of_face = *ind_of_face;
                if (all_faces.len() as i32) < ind_of_face {
                    continue;
                }
                self.my_generated
                    .push(all_faces[(ind_of_face - 1) as usize].clone());
            }

            // OCCT L1413-1430: the next faces for the ruled case.
            if self.my_is_ruled {
                for i in 2..self.my_wires.len() {
                    for ind_of_face in indices.iter() {
                        let mut ind_of_face = *ind_of_face;
                        ind_of_face += ((i - 1) * self.my_nb_edges_in_section as usize) as i32;
                        if (all_faces.len() as i32) < ind_of_face {
                            continue;
                        }
                        self.my_generated
                            .push(all_faces[(ind_of_face - 1) as usize].clone());
                    }
                }
            }
        } else if s.shape_type() == ShapeType::Vertex {
            // OCCT L1432-1436.
            let Some(eindex0) = self.my_vertex_index.get(&s.ptr_id()) else {
                return self.my_generated.clone();
            };

            // OCCT L1438-1439: VEmap.
            let mut ve_map = indexmap::IndexMap::new();
            crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
                &self.my_shape.clone(),
                ShapeType::Vertex,
                ShapeType::Edge,
                &mut ve_map,
            );

            // OCCT L1441-1466: the degenerated end-section check.
            let mut is_degen = [false, false];
            if self.my_degen1 || self.my_degen2 {
                let end_sections = [
                    self.my_wires[0].clone(),
                    self.my_wires[self.my_wires.len() - 1].clone(),
                ];
                for (i, end_section) in end_sections.iter().enumerate() {
                    if i == 0 && !self.my_degen1 {
                        continue;
                    }
                    if i == 1 && !self.my_degen2 {
                        continue;
                    }

                    let explo = explorer(end_section, ShapeType::Vertex, ShapeType::Shape);
                    if let Some(a_vertex) = explo.first() {
                        if s.is_same(a_vertex) {
                            is_degen[i] = true;
                            break;
                        }
                    }
                }
            }

            // OCCT L1468-1573: the start/end degenerated section.
            if is_degen[0] || is_degen[1] {
                // OCCT L1472-1473.
                let mut ve_map2 = indexmap::IndexMap::new();
                crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
                    &self.my_shape.clone(),
                    ShapeType::Vertex,
                    ShapeType::Edge,
                    &mut ve_map2,
                );
                let mut e_map = ShapeIndexedMap::new();
                // OCCT L1475-1480.
                let mut a_new_shape = s.clone();
                if (self.my_is_ruled || !self.my_mutable_input)
                    && self.my_bf_generator.is_some()
                {
                    a_new_shape = self
                        .my_bf_generator
                        .as_mut()
                        .unwrap()
                        .result_shape(&mut self.my_brep, s);
                }

                // OCCT L1482-1499.
                if let Some((_, an_edge_list)) =
                    ve_map2.get(&(a_new_shape.ptr_id(), a_new_shape.location))
                {
                    for value in an_edge_list.clone() {
                        let an_edge = value;
                        if !brep_tool_degenerated(&an_edge) {
                            let (vv0, vv1) = top_exp_vertices(&self.my_brep, &an_edge);
                            if (is_degen[0] && a_new_shape.is_same(&vv0))
                                || (is_degen[1] && a_new_shape.is_same(&vv1))
                            {
                                e_map.add(&an_edge);
                            }
                        }
                    }
                }

                // OCCT L1500-1536.
                for j in 1..=e_map.extent() {
                    let mut an_edge = e_map.key(j);
                    self.my_generated.push(an_edge.clone());
                    if self.my_is_ruled {
                        let mut i = 2usize;
                        let mut k = self.my_wires.len() as i32 - 1;
                        while i < self.my_wires.len() {
                            let ind_of_sec = if is_degen[0] { i } else { k as usize };
                            let a_vertex = if is_degen[0] {
                                last_vertex(&self.my_brep, &an_edge)
                            } else {
                                first_vertex(&self.my_brep, &an_edge)
                            };
                            let (_, e_elist) = ve_map2
                                .get(&(a_vertex.ptr_id(), a_vertex.location))
                                .cloned()
                                .unwrap_or_else(|| (Shape::null(), Vec::new()));
                            // OCCT L1510-1512: MapShapes(wireSection, EDGE,
                            // EmapOfSection).
                            let mut a_wire_section = self.my_wires[ind_of_sec - 1].clone();
                            if (self.my_is_ruled || !self.my_mutable_input)
                                && self.my_bf_generator.is_some()
                            {
                                a_wire_section = self
                                    .my_bf_generator
                                    .as_mut()
                                    .unwrap()
                                    .result_shape(&mut self.my_brep, &a_wire_section);
                            }
                            let mut e_map_of_section = ShapeIndexedMap::new();
                            for e in sub_shapes(&a_wire_section) {
                                if e.shape_type() == ShapeType::Edge {
                                    e_map_of_section.add(&e);
                                }
                            }
                            // OCCT L1514-1524: find the next edge.
                            let mut next_edge = Shape::null();
                            for value in e_elist.iter() {
                                next_edge = value.clone();
                                if !next_edge.is_same(&an_edge)
                                    && !e_map_of_section.contains(&next_edge)
                                {
                                    break;
                                }
                            }
                            // OCCT L1525-1526.
                            self.my_generated.push(next_edge.clone());
                            an_edge = next_edge;
                            i += 1;
                            k -= 1;
                        }
                    }
                }
                // OCCT L1537.
                return self.my_generated.clone();
            }

            // OCCT L1539-1543: the first longitudinal edge.
            let mut eindex = *eindex0;
            let vindex = if eindex > 0 { 0usize } else { 1usize };
            eindex = eindex.abs();

            let mut first_face = all_faces[(eindex - 1) as usize].clone();
            first_face.orientation = Orientation::Forward;
            let explo = explorer(&first_face, ShapeType::Edge, ShapeType::Shape);
            // OCCT L1546: BRepAdaptor_Surface BAsurf(FirstFace, false).
            let is_plane = matches!(
                crate::brep_algo::tool::brep_tool_surface(&first_face),
                Some(Surface3::Plane(_))
            );

            // OCCT L1548: MapShapesAndAncestors(FirstFace, VERTEX, EDGE,
            // VEmap).
            ve_map.clear();
            crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
                &first_face,
                ShapeType::Vertex,
                ShapeType::Edge,
                &mut ve_map,
            );

            let mut an_edge: Shape;
            if self.my_degen1 && is_plane {
                // OCCT L1550-1559: only 3 edges in the face — take the 1-st
                // or the 3-rd.
                let mut idx = 0usize;
                if vindex == 0 {
                    idx += 2;
                }
                an_edge = explo.get(idx).cloned().unwrap_or_else(Shape::null);
            } else {
                // OCCT L1561-1606.
                let mut first_edge = Shape::null();
                let mut first_vertex_of_first_edge = Shape::null();
                let first_section = self.my_wires[0].clone();
                let mut a_wire_explorer = BRepToolsWireExplorer::new();
                a_wire_explorer.init(&first_section, &Shape::null());
                let mut i = 1usize;
                while a_wire_explorer.more() {
                    first_edge = a_wire_explorer.current();
                    if i == eindex as usize {
                        if (self.my_is_ruled || !self.my_mutable_input)
                            && self.my_bf_generator.is_some()
                        {
                            first_edge = self
                                .my_bf_generator
                                .as_mut()
                                .unwrap()
                                .result_shape(&mut self.my_brep, &first_edge);
                        }
                        // OCCT L1574: CurrentVertex().
                        let (_, vl) = top_exp_vertices(&self.my_brep, &first_edge);
                        first_vertex_of_first_edge = vl;
                        break;
                    }
                    a_wire_explorer.next();
                    i += 1;
                }

                // OCCT L1578-1579.
                let first_edge_in_face = explo.first().cloned().unwrap_or_else(Shape::null);
                let (vv0, vv1) = top_exp_vertices(&self.my_brep, &first_edge);
                // OCCT L1581-1591.
                let first_vertex = if vindex == 0 {
                    if vv0.is_same(&first_vertex_of_first_edge) {
                        vv0
                    } else {
                        vv1
                    }
                } else if vv0.is_same(&first_vertex_of_first_edge) {
                    vv1
                } else {
                    vv0
                };
                // OCCT L1592-1603.
                let (_, e_list) = ve_map
                    .get(&(first_vertex.ptr_id(), first_vertex.location))
                    .cloned()
                    .unwrap_or_else(|| (Shape::null(), Vec::new()));
                let an_edge_or = if vindex == 0 {
                    Orientation::Reversed
                } else {
                    Orientation::Forward
                };
                let mut found = Shape::null();
                for value in e_list.iter() {
                    an_edge = value.clone();
                    if !an_edge.is_same(&first_edge_in_face)
                        && !brep_tool_degenerated(&an_edge)
                        && an_edge.orientation == an_edge_or
                    {
                        found = an_edge;
                        break;
                    }
                }
                an_edge = found;
            }

            // OCCT L1608.
            self.my_generated.push(an_edge.clone());
            // OCCT L1609-1620: the chain of longitudinal edges.
            if self.my_is_ruled {
                for _i in 2..self.my_wires.len() {
                    // OCCT L1611-1614.
                    let first_vertex = last_vertex(&self.my_brep, &an_edge);
                    let (_, e_list1) = ve_map
                        .get(&(first_vertex.ptr_id(), first_vertex.location))
                        .cloned()
                        .unwrap_or_else(|| (Shape::null(), Vec::new()));
                    let first_edge = if !e_list1.is_empty()
                        && an_edge.is_same(&e_list1[0])
                    {
                        e_list1.last().cloned().unwrap_or_else(Shape::null)
                    } else {
                        e_list1.first().cloned().unwrap_or_else(Shape::null)
                    };
                    // OCCT L1615-1617.
                    eindex += self.my_nb_edges_in_section;
                    let mut first_face = Shape::null();
                    if (eindex as usize) <= all_faces.len() {
                        first_face = all_faces[(eindex - 1) as usize].clone();
                    }
                    first_face.orientation = Orientation::Forward;
                    // OCCT L1618-1619.
                    ve_map.clear();
                    crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
                        &first_face,
                        ShapeType::Vertex,
                        ShapeType::Edge,
                        &mut ve_map,
                    );
                    // OCCT L1620-1623.
                    let (_, e_list2) = ve_map
                        .get(&(first_vertex.ptr_id(), first_vertex.location))
                        .cloned()
                        .unwrap_or_else(|| (Shape::null(), Vec::new()));
                    an_edge = if !e_list2.is_empty() && first_edge.is_same(&e_list2[0]) {
                        e_list2.last().cloned().unwrap_or_else(Shape::null)
                    } else {
                        e_list2.first().cloned().unwrap_or_else(Shape::null)
                    };
                    // OCCT L1624.
                    self.my_generated.push(an_edge.clone());
                }
            }
        }

        // OCCT L1740.
        self.my_generated.clone()
    }
}

// The HashMap import anchors the myEdgeNewIndices/myVertexIndex map forms.
const _: Option<HashMap<u64, ()>> = None;
