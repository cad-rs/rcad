//! OCCT BRepBuilderAPI_Sewing.cxx — the shape-analysis group:
//! `FaceAnalysis` (L2597-2876), `FindFreeBoundaries` (L2878-3083) and the
//! vertices-assembling group with its file-scope statics `CreateNewNodes`
//! (L3085-3232), `IsMergedVertices` (L3234-3300), `GlueVertices`
//! (L3302-3547) and `VerticesAssembling` (L3549-3606).

use glam::DVec3;
use rcad_kernel::topo::topods::{
    BRep, BRepBuilder, CurveRepresentation, Orientation, ShapeType,
};
use rcad_kernel::topo_shape::Shape;
use std::sync::Arc;

use crate::brep_algo::tool as bat;

use super::{
    idx_list_add, idx_map_add, make_vertex, set_add,
    set_contains,
};
use super::BRepBuilderAPISewing;
use super::selectors::{CellFilter, VertexInspector};

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::FaceAnalysis(theProgress) (cxx
    /// L2597-2876) — modifies myNbShapes/myOldShapes, constructs
    /// myDegenerated; removes the small edges/faces through the reshape.
    pub(crate) fn face_analysis(&mut self, brep: &mut BRep) {
        // OCCT L2600-2604.
        if !self.my_shape.is_null() && self.my_old_shapes.is_empty() {
            let shape = self.my_shape.clone();
            self.add(brep, &shape);
            self.my_shape = Shape::null();
        }

        // OCCT L2606-2609.
        let mut small_edges: super::ShapeSet = super::ShapeSet::new();
        let mut glued_vertices: super::IdxListMap = super::IdxListMap::new();
        // OCCT L2611-2615: the myOldShapes walk.
        let nb_shapes = self.my_old_shapes.len();
        for i in 1..=nb_shapes {
            // OCCT L2616: myOldShapes(i) is the map VALUE (Add stored
            // Apply(aShape) as the value, cxx L2177-2178).
            let (_, old_shape) = super::idx_shape_get(&self.my_old_shapes, i - 1);
            // OCCT L2616: for (TopExp_Explorer fexp(myOldShapes(i), TopAbs_FACE)...).
            for fexp in bat::explorer(&old_shape, ShapeType::Face, ShapeType::Shape) {
                // Retrieve current face
                // OCCT L2619-2621.
                let face = fexp;
                let mut nb_edges = 0i32;
                let mut nb_small = 0i32;

                // Build replacing face
                // OCCT L2624-2627.
                let mut nface = brep.empty_copied(&face);
                nface.orientation = Orientation::Forward;
                let mut is_face_changed = false;

                // OCCT L2629: TopoDS_Iterator witer(face.Oriented(TopAbs_FORWARD))
                // (the cumOri = true default).
                for witer in bat::sub_shapes(&bat::oriented(&face, Orientation::Forward)) {
                    // Retrieve current wire
                    // OCCT L2632-2637.
                    if witer.shape_type() != ShapeType::Wire {
                        continue;
                    }
                    let wire = witer;

                    // Build replacing wire
                    // OCCT L2640-2643.
                    let mut nwire = brep.empty_copied(&wire);
                    nwire.orientation = Orientation::Forward;
                    let mut is_wire_changed = false;

                    // OCCT L2645: TopoDS_Iterator eiter(wire.Oriented(TopAbs_FORWARD)).
                    for eiter in bat::sub_shapes(&bat::oriented(&wire, Orientation::Forward)) {
                        // Retrieve current edge
                        // OCCT L2648-2651.
                        let edge = eiter;
                        nb_edges += 1;

                        // Process degenerated edge
                        // OCCT L2654-2659.
                        if bat::brep_tool_degenerated(&edge) {
                            let mut b = BRepBuilder::new();
                            b.add_to_wire(brep, nwire.clone(), edge.clone()); // Old edge kept
                            set_add(&mut self.my_degenerated, &edge);
                            nb_small += 1;
                            continue;
                        }

                        // OCCT L2661-2664.
                        let mut is_small = set_contains(&small_edges, &edge);
                        if !is_small {
                            // Check for small edge
                            // OCCT L2666-2669.
                            let c3d = bat::brep_tool_curve(&edge);
                            match c3d {
                                None => {
                                    // OCCT L2671-2673 (debug warning only).
                                }
                                Some((c3d, first, last)) => {
                                    use rcad_kernel::geom::CurveEval;
                                    // Evaluate curve compactness
                                    // OCCT L2680-2692.
                                    let npt = 5usize;
                                    let cp = (c3d.point_at(first) + c3d.point_at(last)) * 0.5;
                                    let mut maxdist = 0.0f64;
                                    let delta = (last - first) / (npt as f64 - 1.0);
                                    for idx in 0..npt {
                                        let dist =
                                            cp.distance(c3d.point_at(first + idx as f64 * delta));
                                        if maxdist < dist {
                                            maxdist = dist;
                                        }
                                    }
                                    // OCCT L2693.
                                    is_small = 2.0 * maxdist <= self.min_tolerance();
                                }
                            }

                            // OCCT L2705-2765.
                            if is_small {
                                // Store small edge in the map
                                set_add(&mut small_edges, &edge);

                                // OCCT L2710-2711.
                                let (v1, v2) = bat::top_exp_vertices_raw(&edge);
                                let v1 = v1.unwrap_or_else(Shape::null);
                                let v2 = v2.unwrap_or_else(Shape::null);
                                let nv1 =
                                    self.my_re_shape.apply(brep, &v1, ShapeType::Shape);
                                let nv2 =
                                    self.my_re_shape.apply(brep, &v2, ShapeType::Shape);

                                // Store glued vertices
                                // OCCT L2714-2764.
                                if !nv1.is_same(&v1) {
                                    // First vertex was already glued
                                    if !nv2.is_same(&v2) {
                                        // Merge lists of glued vertices
                                        if !nv1.is_same(&nv2) {
                                            let vlist2 = glued_vertices
                                                .get(&bat::shape_key(&nv2))
                                                .cloned();
                                            if let Some((_, vlist2)) = vlist2 {
                                                let vlist1 =
                                                    glued_vertices.entry(bat::shape_key(&nv1))
                                                        .or_insert((nv1.clone(), Vec::new()));
                                                for liter in &vlist2 {
                                                    let v = liter.clone();
                                                    let oriented_v =
                                                        bat::oriented(&nv1, v.orientation);
                                                    self.my_re_shape.replace(
                                                        brep, &v, &oriented_v,
                                                    );
                                                    vlist1.1.push(v);
                                                }
                                                glued_vertices
                                                    .shift_remove(&bat::shape_key(&nv2));
                                            }
                                        }
                                    } else {
                                        // Add second vertex to the existing list
                                        let vlist1 = glued_vertices
                                            .entry(bat::shape_key(&nv1))
                                            .or_insert((nv1.clone(), Vec::new()));
                                        vlist1.1.push(v2.clone());
                                        let oriented_v = bat::oriented(&nv1, v2.orientation);
                                        self.my_re_shape.replace(brep, &v2, &oriented_v);
                                    }
                                } else if !nv2.is_same(&v2) {
                                    // Add first vertex to the existing list
                                    let vlist2 = glued_vertices
                                        .entry(bat::shape_key(&nv2))
                                        .or_insert((nv2.clone(), Vec::new()));
                                    vlist2.1.push(v1.clone());
                                    let oriented_v = bat::oriented(&nv2, v1.orientation);
                                    self.my_re_shape.replace(brep, &v1, &oriented_v);
                                } else if !v1.is_same(&v2) {
                                    // Record new glued vertices
                                    let nv = make_vertex(brep);
                                    let vlist = vec![v1.clone(), v2.clone()];
                                    idx_list_add(&mut glued_vertices, &nv, vlist);
                                    let o1 = bat::oriented(&nv, v1.orientation);
                                    self.my_re_shape.replace(brep, &v1, &o1);
                                    let o2 = bat::oriented(&nv, v2.orientation);
                                    self.my_re_shape.replace(brep, &v2, &o2);
                                }
                            }
                        }

                        // Replace small edge
                        // OCCT L2768-2800.
                        if is_small {
                            nb_small += 1;
                            // Create new degenerated edge
                            let fedge = bat::oriented(&edge, Orientation::Forward);
                            let c2d = bat::brep_tool_curve_on_surface(&fedge, &face);
                            if let Some((c2d, pfirst, plast)) = c2d {
                                // OCCT L2782-2796.
                                let nedge = super::same_param::builder_make_edge(brep);
                                let mut a_builder = BRepBuilder::new();
                                a_builder.update_edge_pcurve(
                                    brep,
                                    nedge.clone(),
                                    c2d.clone(),
                                    face.clone(),
                                    rcad_kernel::core::precision::CONFUSION,
                                );
                                a_builder.set_edge_range(brep, nedge.clone(), pfirst, plast);
                                a_builder.set_edge_degenerated(brep, nedge.clone(), true);
                                let (rv1, rv2) = bat::top_exp_vertices_raw(&fedge);
                                let v1 = rv1.unwrap_or_else(Shape::null);
                                let v2 = rv2.unwrap_or_else(Shape::null);
                                let nv1 = self.my_re_shape.apply(brep, &v1, ShapeType::Shape);
                                let nv2 = self.my_re_shape.apply(brep, &v2, ShapeType::Shape);
                                a_builder.add_to_edge(
                                    brep,
                                    nedge.clone(),
                                    bat::oriented(&nv1, v1.orientation),
                                );
                                a_builder.add_to_edge(
                                    brep,
                                    nedge.clone(),
                                    bat::oriented(&nv2, v2.orientation),
                                );
                                let oriented_nedge =
                                    bat::oriented(&nedge, edge.orientation);
                                let mut b = BRepBuilder::new();
                                b.add_to_wire(brep, nwire.clone(), oriented_nedge);
                                set_add(&mut self.my_degenerated, &nedge);
                            }
                            is_wire_changed = true;
                        } else {
                            // OCCT L2797-2800: B.Add(nwire, edge) — old edge kept.
                            let mut b = BRepBuilder::new();
                            b.add_to_wire(brep, nwire.clone(), edge.clone());
                        }
                    }

                    // Record wire in the new face
                    // OCCT L2803-2810.
                    if is_wire_changed {
                        let mut b = BRepBuilder::new();
                        let oriented_wire = bat::oriented(&nwire, wire.orientation);
                        b.add_to_face(brep, nface.clone(), oriented_wire);
                        is_face_changed = true;
                    } else {
                        let mut b = BRepBuilder::new();
                        b.add_to_face(brep, nface.clone(), wire.clone());
                    }
                }

                // Remove small face
                // OCCT L2813-2831.
                if nb_small == nb_edges {
                    set_add(&mut self.my_little_face, &face);
                    self.my_re_shape.remove(brep, &face);
                } else if is_face_changed {
                    let oriented_nface = bat::oriented(&nface, face.orientation);
                    self.my_re_shape.replace(brep, &face, &oriented_nface);
                }
                let _ = &mut nface;
            }
        }

        // Update glued vertices
        // OCCT L2834-2870.
        for idx in 0..glued_vertices.len() {
            let (vnew, vlist) = super::idx_list_get(&glued_vertices, idx);
            let mut coord = DVec3::ZERO;
            let mut nb_points = 0i32;
            for liter1 in vlist.iter() {
                coord += bat::brep_tool_pnt(liter1).unwrap_or(DVec3::ZERO);
                nb_points += 1;
            }
            if nb_points != 0 {
                let vp = coord / (nb_points as f64);
                let mut tol = 0.0f64;
                let mut mtol = 0.0f64;
                for liter2 in vlist.iter() {
                    let mut vtol = bat::brep_tool_tolerance(liter2);
                    if mtol < vtol {
                        mtol = vtol;
                    }
                    vtol = vp.distance(bat::brep_tool_pnt(liter2).unwrap_or(DVec3::ZERO));
                    if tol < vtol {
                        tol = vtol;
                    }
                }
                let mut b = BRepBuilder::new();
                // OCCT L2868: B.UpdateVertex(vnew, vp, tol + mtol).
                b.update_vertex_point(brep, vnew.clone(), vp, tol + mtol);
            }
        }

        // Update input shapes
        // OCCT L2873-2876.
        for i in 0..self.my_old_shapes.len() {
            let (k, old_val) = {
                let e = self.my_old_shapes.get_index(i).unwrap();
                (e.0.clone(), e.1.1.clone())
            };
            let applied = self.my_re_shape.apply(brep, &old_val, ShapeType::Shape);
            if let Some(entry) = self.my_old_shapes.get_mut(&k) {
                entry.1 = applied;
            }
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::FindFreeBoundaries() (cxx L2878-3083) —
    /// constructs myBoundFaces, myVertexNode and myVertexNodeFree.
    pub(crate) fn find_free_boundaries(&mut self, brep: &mut BRep) {
        // Take into account the context shape if needed
        // OCCT L2881-2894.
        let mut new_shapes: super::ShapeSet = super::ShapeSet::new();
        if !self.my_shape.is_null() {
            if self.my_old_shapes.is_empty() {
                let shape = self.my_shape.clone();
                self.add(brep, &shape);
                self.my_shape = Shape::null();
            } else {
                let new_shape = self.my_re_shape.apply(brep, &self.my_shape, ShapeType::Shape);
                if !new_shape.is_null() {
                    set_add(&mut new_shapes, &new_shape);
                }
            }
        }
        // Create map Edge -> Faces
        // OCCT L2896-2915.
        let mut edge_faces: super::IdxListMap = super::IdxListMap::new();
        let nb_shapes = self.my_old_shapes.len();
        for i in 1..=nb_shapes {
            // Retrieve new shape
            // OCCT L2899: const TopoDS_Shape& shape = myOldShapes(i) — the map
            // VALUE (Add stored Apply(aShape) as the value, cxx L2177-2178).
            let (_, shape) = super::idx_shape_get(&self.my_old_shapes, i - 1);
            if shape.is_null() {
                continue;
            }
            set_add(&mut new_shapes, &shape);
            // Explore shape to find all boundaries
            for e_exp in bat::explorer(&shape, ShapeType::Edge, ShapeType::Shape) {
                let edge = e_exp;
                if !edge_faces.contains_key(&bat::shape_key(&edge)) {
                    idx_list_add(&mut edge_faces, &edge, Vec::new());
                }
            }
        }
        // Fill map Edge -> Faces
        // OCCT L2917-2951.
        let nb_shapes = new_shapes.len();
        let mut map_faces: super::ShapeSet = super::ShapeSet::new();
        for i in 1..=nb_shapes {
            // Explore shape to find all faces
            let key_shape = super::set_get(&new_shapes, i - 1);
            for f_exp in bat::explorer(&key_shape, ShapeType::Face, ShapeType::Shape) {
                let face = f_exp;
                if set_contains(&map_faces, &face) {
                    continue;
                } else {
                    set_add(&mut map_faces, &face);
                }
                // Explore face to find all boundaries
                // OCCT L2933-2936: TopoDS_Iterator aIw(face) (cumOri = true).
                for a_iw in bat::sub_shapes(&face) {
                    if a_iw.shape_type() != ShapeType::Wire {
                        continue;
                    }
                    // OCCT L2937-2939: TopoDS_Iterator aIIe(aIw.Value()).
                    for a_iie in bat::sub_shapes(&a_iw) {
                        let edge = a_iie;
                        // OCCT L2943-2947.
                        if edge_faces.contains_key(&bat::shape_key(&edge)) {
                            if let Some(entry) = edge_faces.get_mut(&bat::shape_key(&edge)) {
                                entry.1.push(face.clone());
                            }
                        }
                    }
                }
            }
        }
        // Find free boundaries
        // OCCT L2954-3083.
        for idx in 0..edge_faces.len() {
            let (edge_shape, list_faces) = super::idx_list_get(&edge_faces, idx);
            let nb_faces = list_faces.len();
            let mut edge = edge_shape.clone();
            // OCCT L2959-2962.
            if edge.orientation == Orientation::Internal {
                continue;
            }
            let mut is_seam = false;
            if nb_faces == 1 {
                // OCCT L2965-2967.
                let face0 = list_faces.first().cloned().unwrap_or_else(Shape::null);
                is_seam = bat::brep_tool_is_closed_on_surface(&edge, &face0);
                if is_seam {
                    // OCCT L2972-2999.
                    let mut an_edge = brep.empty_copied(&edge);
                    // OCCT L2975-2980: TopoDS_Iterator aItV(edge) — add the
                    // vertices.
                    let mut a_b = BRepBuilder::new();
                    for a_it_v in bat::sub_shapes(&edge) {
                        a_b.add_to_edge(brep, an_edge.clone(), a_it_v.clone());
                    }

                    // OCCT L2983-2989.
                    let c2dold = bat::brep_tool_curve_on_surface(&edge, &face0);

                    // OCCT L2990: occ::handle<Geom2d_Curve> c2d; — a NULL
                    // pcurve pair.
                    // OCCT L2991 (statement (1)): B.UpdateEdge(anewEdge, c2d,
                    // c2d, face, 0) (BRep_Builder.cxx L702-725 ->
                    // UpdateCurves(TE->ChangeCurves(), C1, C2, S, L)
                    // L251-308): the first on-face representation matching
                    // (S, L) is REMOVED and the null pair appends NOTHING.
                    {
                        let fkeys: Vec<(u64, u32)> = brep
                            .edge_wrapper_locations(&an_edge)
                            .iter()
                            .map(|&el| {
                                (
                                    face0.ptr_id(),
                                    brep.compose_pcurve_location(face0.location, el),
                                )
                            })
                            .collect();
                        let ed = brep.edge_mut_inplace(an_edge.clone());
                        for k in &fkeys {
                            ed.pcurves.shift_remove(k);
                        }
                        ed.representations.retain(|r| match r {
                            CurveRepresentation::CurveOnSurface { face, .. }
                            | CurveRepresentation::CurveOnClosedSurface { face, .. } => {
                                !fkeys.contains(face)
                            }
                            _ => true,
                        });
                    }
                    // OCCT L2992 (statement (2)): B.UpdateEdge(anewEdge,
                    // c2dold, face, 0) (BRep_Builder.cxx L655-671 ->
                    // UpdateCurves(TE->ChangeCurves(), C, S, L) L104-167): the
                    // remaining on-face representations are removed while
                    // scanning, then the single plain CurveOnSurface(c2dold)
                    // is appended.
                    if let Some((c2d, _, _)) = &c2dold {
                        let mut b = BRepBuilder::new();
                        b.update_edge_pcurve(brep, an_edge.clone(), c2d.clone(), face0.clone(), 0.0);
                    }

                    // OCCT L2994-2998.
                    let (a_first, a_last) = bat::brep_tool_range(&edge);
                    a_b.set_edge_range(brep, an_edge.clone(), a_first, a_last);
                    // OCCT L2998: aB.Range(anewEdge, face, first2d, last2d).
                    // (the pcurve-range form; the c2dold range is set.)
                    if let Some((_, first2d, last2d)) = &c2dold {
                        let mut e = an_edge.clone();
                        bat::builder_range_edge_on_face(&mut e, &face0, *first2d, *last2d);
                    }
                    let oriented = bat::oriented(&an_edge, edge.orientation);
                    self.my_re_shape.replace(brep, &edge, &oriented);
                    edge = oriented.clone();
                    an_edge = oriented;

                    // OCCT L3003: isSeam = false;
                    is_seam = false;
                }
            }
            // OCCT L3004-3005.
            let is_bound_float = self.my_floating_edges_mode && nb_faces == 0;
            let is_bound = self.my_face_mode
                && ((self.my_nonmanifold && nb_faces > 0)
                    || (nb_faces == 1 && !is_seam));
            if is_bound || is_bound_float {
                // Ignore degenerated edge
                // OCCT L3008-3011.
                if bat::brep_tool_degenerated(&edge) {
                    continue;
                }
                // Add to BoundFaces
                // OCCT L3016-3018.
                let list_faces_copy = list_faces.clone();
                idx_list_add(&mut self.my_bound_faces, &edge, list_faces_copy);
                // Process edge vertices
                // OCCT L3020-3022.
                let (v_first, v_last) = bat::top_exp_vertices_raw(&edge);
                let v_first = v_first.unwrap_or_else(Shape::null);
                let v_last = v_last.unwrap_or_else(Shape::null);
                if v_first.is_null() || v_last.is_null() {
                    continue;
                }
                // OCCT L3024-3027.
                if v_first.orientation == Orientation::Internal
                    || v_last.orientation == Orientation::Internal
                {
                    continue;
                }
                if is_bound {
                    // Add to VertexNode
                    // OCCT L3030-3037.
                    if !self.my_vertex_node.contains_key(&bat::shape_key(&v_first)) {
                        idx_map_add(&mut self.my_vertex_node, &v_first, v_first.clone());
                    }
                    if !self.my_vertex_node.contains_key(&bat::shape_key(&v_last)) {
                        idx_map_add(&mut self.my_vertex_node, &v_last, v_last.clone());
                    }
                } else {
                    // Add to VertexNodeFree
                    // OCCT L3039-3046.
                    if !self.my_vertex_node_free.contains_key(&bat::shape_key(&v_first)) {
                        idx_map_add(&mut self.my_vertex_node_free, &v_first, v_first.clone());
                    }
                    if !self.my_vertex_node_free.contains_key(&bat::shape_key(&v_last)) {
                        idx_map_add(&mut self.my_vertex_node_free, &v_last, v_last.clone());
                    }
                }
            }
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::VerticesAssembling(theProgress) (cxx
    /// L3549-3606) — glues the node vertices.
    pub(crate) fn vertices_assembling(&mut self, brep: &mut BRep) {
        // OCCT L3550-3551.
        let nb_vert = self.my_vertex_node.len();
        let nb_vert_free = self.my_vertex_node_free.len();
        if nb_vert != 0 || nb_vert_free != 0 {
            // Fill map node -> sections
            // OCCT L3555-3571.
            for i in 1..=self.my_bound_faces.len() {
                let (bound, _) = super::idx_list_get(&self.my_bound_faces, i - 1);
                // OCCT L3558: TopoDS_Iterator itv(bound, false).
                for itv in super::iter_no_cumori(&bound) {
                    let node = itv;
                    match self.my_node_sections.get_mut(&bat::shape_key(&node)) {
                        Some(list) => list.push(bound.clone()),
                        None => {
                            self.my_node_sections
                                .insert(bat::shape_key(&node), vec![bound.clone()]);
                        }
                    }
                }
            }
            // Glue vertices
            if nb_vert != 0 {
                // OCCT L3582-3586: while (GlueVertices(myVertexNode, ...)).
                loop {
                    let again = glue_vertices(
                        brep,
                        &mut self.my_vertex_node,
                        &mut self.my_node_sections,
                        &self.my_bound_faces,
                        self.my_tolerance,
                    );
                    if !again {
                        break;
                    }
                }
            }
            if nb_vert_free != 0 {
                // OCCT L3590-3595.
                loop {
                    let again = glue_vertices(
                        brep,
                        &mut self.my_vertex_node_free,
                        &mut self.my_node_sections,
                        &self.my_bound_faces,
                        self.my_tolerance,
                    );
                    if !again {
                        break;
                    }
                }
            }
        }
    }
}

/// OCCT static CreateNewNodes(NodeNearestNode, NodeVertices, aVertexNode,
/// aNodeEdges) (cxx L3085-3232).
fn create_new_nodes(
    brep: &mut BRep,
    node_nearest_node: &super::IdxShapeMap,
    node_vertices: &super::IdxListMap,
    a_vertex_node: &mut super::IdxShapeMap,
    a_node_edges: &mut super::DataListMap,
) -> bool {
    // Create new nodes
    // OCCT L3090-3092.
    let mut old_node_new_node: super::DataShapeMap = super::DataShapeMap::new();
    let mut new_node_old_nodes: super::IdxListMap = super::IdxListMap::new();
    // OCCT L3093-3096: the NodeNearestNode walk.
    for idx in 0..node_nearest_node.len() {
        // Retrieve a pair of nodes to merge
        let (oldnode1, oldnode2) = super::idx_shape_get(node_nearest_node, idx);
        // Second node should also be in the map
        // OCCT L3100-3103.
        if !node_nearest_node.contains_key(&bat::shape_key(&oldnode2)) {
            continue;
        }
        // Get new node for old node #1
        // OCCT L3106-3155.
        if let Some(newnode1) = old_node_new_node.get(&bat::shape_key(&oldnode1)).cloned() {
            if let Some(newnode2) = old_node_new_node.get(&bat::shape_key(&oldnode2)).cloned() {
                if !newnode1.is_same(&newnode2) {
                    // Change data for new node #2
                    // OCCT L3115-3128.
                    let lnode2 = new_node_old_nodes
                        .get(&bat::shape_key(&newnode2))
                        .cloned();
                    if let Some((_, lnode2)) = lnode2 {
                        let lnode1 = new_node_old_nodes
                            .entry(bat::shape_key(&newnode1))
                            .or_insert((newnode1.clone(), Vec::new()));
                        for itn in &lnode2 {
                            let node2 = itn.clone();
                            lnode1.1.push(node2.clone());
                            old_node_new_node.insert(bat::shape_key(&node2), newnode1.clone());
                        }
                        new_node_old_nodes.shift_remove(&bat::shape_key(&newnode2));
                    }
                }
            } else {
                // Old node #2 is not bound - add to old node #1
                // OCCT L3131-3135.
                old_node_new_node.insert(bat::shape_key(&oldnode2), newnode1.clone());
                let lnode1 = new_node_old_nodes
                    .entry(bat::shape_key(&newnode1))
                    .or_insert((newnode1.clone(), Vec::new()));
                lnode1.1.push(oldnode2.clone());
            }
        } else if let Some(newnode2) = old_node_new_node.get(&bat::shape_key(&oldnode2)).cloned() {
            // Old node #1 is not bound - add to old node #2
            // OCCT L3137-3142.
            old_node_new_node.insert(bat::shape_key(&oldnode1), newnode2.clone());
            let lnode2 = new_node_old_nodes
                .entry(bat::shape_key(&newnode2))
                .or_insert((newnode2.clone(), Vec::new()));
            lnode2.1.push(oldnode1.clone());
        } else {
            // Nodes are not bound - create new node
            // OCCT L3144-3154.
            let newnode = make_vertex(brep);
            old_node_new_node.insert(bat::shape_key(&oldnode1), newnode.clone());
            old_node_new_node.insert(bat::shape_key(&oldnode2), newnode.clone());
            let lnodes = vec![oldnode1.clone(), oldnode2.clone()];
            idx_list_add(&mut new_node_old_nodes, &newnode, lnodes);
        }
    }

    // Stop if no new nodes created
    // OCCT L3158-3161.
    if new_node_old_nodes.is_empty() {
        return false;
    }

    // OCCT L3163-3219.
    for idx in 0..new_node_old_nodes.len() {
        let (newnode, old_nodes) = super::idx_list_get(&new_node_old_nodes, idx);
        // Calculate new node center point
        let mut the_coordinates = DVec3::ZERO;
        let mut lvert: Vec<Shape> = Vec::new(); // Accumulate node vertices
        let mut medge: super::ShapeSet = super::ShapeSet::new();
        let mut ledge: Vec<Shape> = Vec::new(); // Accumulate node edges
        // Iterate on old nodes
        for itn in old_nodes.iter() {
            let oldnode = itn.clone();
            // Iterate on node vertices
            // OCCT L3173-3182.
            let node_verts = node_vertices.get(&bat::shape_key(&oldnode)).cloned();
            if let Some((_, node_verts)) = node_verts {
                for itv in &node_verts {
                    let vertex = itv.clone();
                    // Change node for vertex
                    if let Some(entry) = a_vertex_node.get_mut(&bat::shape_key(&vertex)) {
                        entry.1 = newnode.clone();
                    }
                    // Accumulate coordinates
                    the_coordinates += bat::brep_tool_pnt(&vertex).unwrap_or(DVec3::ZERO);
                    lvert.push(vertex);
                }
            }
            // Iterate on node edges
            // OCCT L3185-3197.
            let edges = a_node_edges.get(&bat::shape_key(&oldnode)).cloned();
            if let Some(edges) = edges {
                for ite in &edges {
                    let edge = ite.clone();
                    if !set_contains(&medge, &edge) {
                        set_add(&mut medge, &edge);
                        ledge.push(edge);
                    }
                }
            }
            // Unbind old node edges
            // OCCT L3198-3199.
            a_node_edges.remove(&bat::shape_key(&oldnode));
        }
        // Bind new node edges
        // OCCT L3201-3202.
        a_node_edges.insert(bat::shape_key(&newnode), ledge);
        // OCCT L3203: gp_Pnt center(theCoordinates / lvert.Extent()).
        let center = the_coordinates / (lvert.len() as f64);
        // Calculate new node tolerance
        // OCCT L3205-3214.
        let mut toler = 0.0f64;
        for itv in &lvert {
            let vertex = itv.clone();
            let t = center.distance(bat::brep_tool_pnt(&vertex).unwrap_or(DVec3::ZERO))
                + bat::brep_tool_tolerance(&vertex);
            if toler < t {
                toler = t;
            }
        }
        // Update new node parameters
        // OCCT L3217: B.UpdateVertex(newnode, center, toler).
        let mut b = BRepBuilder::new();
        b.update_vertex_point(brep, newnode, center, toler);
    }

    // OCCT L3221.
    true
}

/// OCCT static IsMergedVertices(face1, e1, e2, vtx1, vtx2) (cxx L3234-3300).
fn is_merged_vertices(
    brep: &BRep,
    face1: &Shape,
    e1: &Shape,
    e2: &Shape,
    vtx1: &Shape,
    vtx2: &Shape,
) -> i32 {
    // Case of floating edges
    // OCCT L3236-3240.
    if face1.is_null() {
        return if super::closed::is_closed_shape(brep, e1, vtx1, vtx2) { 0 } else { 1 };
    }

    // Find wires containing given edges
    // OCCT L3243-3258.
    let mut wire1 = Shape::null();
    let mut wire2 = Shape::null();
    for itw in bat::explorer(face1, ShapeType::Wire, ShapeType::Shape) {
        if !wire1.is_null() && !wire2.is_null() {
            break;
        }
        // OCCT L3247: TopoDS_Iterator ite(itw.Current(), false).
        for ite in super::iter_no_cumori(&itw) {
            if !wire1.is_null() && !wire2.is_null() {
                break;
            }
            if wire1.is_null() && e1.is_same(&ite) {
                wire1 = itw.clone();
            }
            if wire2.is_null() && e2.is_same(&ite) {
                wire2 = itw.clone();
            }
        }
    }
    // OCCT L3259-3297.
    let mut status = 0i32;
    if !wire1.is_null() && !wire2.is_null() {
        if wire1.is_same(&wire2) {
            // OCCT L3262-3273.
            for a_ite in super::iter_no_cumori(&wire1) {
                let (ve1, ve2) = bat::top_exp_vertices_raw(&a_ite);
                let ve1 = ve1.unwrap_or_else(Shape::null);
                let ve2 = ve2.unwrap_or_else(Shape::null);
                if (ve1.is_same(vtx1) && ve2.is_same(vtx2)) || (ve2.is_same(vtx1) && ve1.is_same(vtx2))
                {
                    return if super::closed::is_closed_shape(brep, &a_ite, vtx1, vtx2) { 0 } else { 1 };
                }
            }
            // OCCT L3274-3288.
            if super::closed::is_closed_shape(brep, &wire1, vtx1, vtx2) {
                let (v1, v2) = bat::top_exp_vertices_wire(&wire1);
                let v1 = v1.unwrap_or_else(Shape::null);
                let v2 = v2.unwrap_or_else(Shape::null);
                let is_end_vertex =
                    (v1.is_same(vtx1) && v2.is_same(vtx2)) || (v2.is_same(vtx1) && v1.is_same(vtx2));
                if !is_end_vertex {
                    status = 1;
                }
            } else {
                status = 1;
            }
        } else {
            // OCCT L3290-3292.
            status = -1;
        }
    }
    // OCCT L3298.
    status
}

/// OCCT static GlueVertices(aVertexNode, aNodeEdges, aBoundFaces, Tolerance,
/// theProgress) (cxx L3302-3547).
fn glue_vertices(
    brep: &mut BRep,
    a_vertex_node: &mut super::IdxShapeMap,
    a_node_edges: &mut super::DataListMap,
    a_bound_faces: &super::IdxListMap,
    tolerance: f64,
) -> bool {
    // Create map of node -> vertices
    // OCCT L3305-3308.
    let mut node_vertices: super::IdxListMap = super::IdxListMap::new();
    // OCCT L3309-3310: the CellFilter + inspector.
    let mut a_filter = CellFilter::new(tolerance);
    let mut an_inspector = VertexInspector::new(tolerance);
    // OCCT L3311-3326.
    for idx in 0..a_vertex_node.len() {
        let (vertex, node) = super::idx_shape_get(a_vertex_node, idx);
        if node_vertices.contains_key(&bat::shape_key(&node)) {
            if let Some(entry) = node_vertices.get_mut(&bat::shape_key(&node)) {
                entry.1.push(vertex);
            }
        } else {
            let vlist = vec![vertex];
            // OCCT L3320-3324: NodeVertices.Add(node, vlist) and the filter
            // registration with the FindIndex rank.
            let a_pnt = bat::brep_tool_pnt(&node).unwrap_or(DVec3::ZERO);
            idx_list_add(&mut node_vertices, &node, vlist);
            let rank = find_index_rank(&node_vertices, &node);
            a_filter.add(rank, a_pnt);
            an_inspector.add(a_pnt);
        }
    }
    let nb_nodes = node_vertices.len();

    // Merge nearest nodes
    // OCCT L3333-3462.
    let mut node_nearest_node: super::IdxShapeMap = super::IdxShapeMap::new();
    for i in 1..=nb_nodes {
        let (node1, _) = super::idx_list_get(&node_vertices, i - 1);
        // Find near nodes
        // OCCT L3337-3343.
        let pt1 = bat::brep_tool_pnt(&node1).unwrap_or(DVec3::ZERO);
        an_inspector.set_current(pt1);
        let a_pnt_min = VertexInspector::shift(pt1, -tolerance);
        let a_pnt_max = VertexInspector::shift(pt1, tolerance);
        a_filter.inspect(a_pnt_min, a_pnt_max, &mut an_inspector);
        if an_inspector.res_ind().is_empty() {
            continue;
        }
        // Retrieve list of edges for the first node
        // OCCT L3346.
        let ledges1 = match a_node_edges.get(&bat::shape_key(&node1)) {
            Some(l) => l.clone(),
            None => Vec::new(),
        };
        // Explore list of near nodes and fill the sequence of glued nodes
        // OCCT L3349-3352.
        let mut seq_nodes: Vec<Shape> = Vec::new();
        let mut list_nodes_same_edge: Vec<Shape> = Vec::new();
        for iter1 in an_inspector.res_ind().clone() {
            // OCCT L3355-3358: node2 = TopoDS::Vertex(NodeVertices.FindKey(iter1)).
            let (node2, _) = super::idx_list_get(&node_vertices, (iter1 - 1) as usize);
            // OCCT L3359-3362.
            if node1.ptr_id() == node2.ptr_id() && node1.location == node2.location {
                continue;
            }
            // Retrieve list of edges for the second node
            // OCCT L3365.
            let ledges2 = match a_node_edges.get(&bat::shape_key(&node2)) {
                Some(l) => l.clone(),
                None => Vec::new(),
            };
            // Check merging condition for the pair of nodes
            // OCCT L3368-3370.
            let mut status = 0i32;
            let mut is_same_edge = false;
            // Explore edges of the first node
            // OCCT L3371-3421.
            for e1 in &ledges1 {
                if status != 0 || is_same_edge {
                    break;
                }
                let e1 = e1.clone();
                // Obtain real vertex from edge
                // OCCT L3375-3390.
                let mut v1 = node1.clone();
                {
                    let (ov1, ov2) = bat::top_exp_vertices_raw(&e1);
                    let ov1 = ov1.unwrap_or_else(Shape::null);
                    let ov2 = ov2.unwrap_or_else(Shape::null);
                    if a_vertex_node.contains_key(&bat::shape_key(&ov1)) {
                        if node1.is_same(&node_of_map(a_vertex_node, &ov1)) {
                            v1 = ov1;
                        }
                    }
                    if a_vertex_node.contains_key(&bat::shape_key(&ov2)) {
                        if node1.is_same(&node_of_map(a_vertex_node, &ov2)) {
                            v1 = ov2;
                        }
                    }
                }
                // Create map of faces for e1
                // OCCT L3393-3403.
                let mut faces1: super::ShapeSet = super::ShapeSet::new();
                let lfac1 = a_bound_faces
                    .get(&bat::shape_key(&e1))
                    .map(|(_, l)| l.clone())
                    .unwrap_or_default();
                if !lfac1.is_empty() {
                    for itf in &lfac1 {
                        if !itf.is_null() {
                            set_add(&mut faces1, itf);
                        }
                    }
                }
                // Explore edges of the second node
                // OCCT L3404-3420.
                for e2 in &ledges2 {
                    if status != 0 || is_same_edge {
                        break;
                    }
                    let e2 = e2.clone();
                    // Obtain real vertex from edge
                    let mut v2 = node2.clone();
                    {
                        let (ov1, ov2) = bat::top_exp_vertices_raw(&e2);
                        let ov1 = ov1.unwrap_or_else(Shape::null);
                        let ov2 = ov2.unwrap_or_else(Shape::null);
                        if a_vertex_node.contains_key(&bat::shape_key(&ov1)) {
                            if node2.is_same(&node_of_map(a_vertex_node, &ov1)) {
                                v2 = ov1;
                            }
                        }
                        if a_vertex_node.contains_key(&bat::shape_key(&ov2)) {
                            if node2.is_same(&node_of_map(a_vertex_node, &ov2)) {
                                v2 = ov2;
                            }
                        }
                    }
                    // Explore faces for e2
                    // OCCT L3415-3434.
                    let lfac2 = a_bound_faces
                        .get(&bat::shape_key(&e2))
                        .map(|(_, l)| l.clone())
                        .unwrap_or_default();
                    if !lfac2.is_empty() {
                        for itf in &lfac2 {
                            if status != 0 || is_same_edge {
                                break;
                            }
                            // Check merging conditions for the same face
                            if set_contains(&faces1, itf) {
                                let stat = is_merged_vertices(brep, itf, &e1, &e2, &v1, &v2);
                                if stat == 1 {
                                    is_same_edge = true;
                                } else {
                                    status = stat;
                                }
                            }
                        }
                    } else if faces1.is_empty() && e1.ptr_id() == e2.ptr_id() && e1.location == e2.location
                    {
                        // OCCT L3435-3443.
                        let stat = is_merged_vertices(brep, &Shape::null(), &e1, &e1, &v1, &v2);
                        if stat == 1 {
                            is_same_edge = true;
                        } else {
                            status = stat;
                        }
                        break;
                    }
                }
            }
            // OCCT L3446-3449.
            if status != 0 {
                continue;
            }
            if is_same_edge {
                list_nodes_same_edge.push(node2.clone());
            }
            // Append near node to the sequence
            // OCCT L3452-3465.
            let pt2 = bat::brep_tool_pnt(&node2).unwrap_or(DVec3::ZERO);
            let dist = pt1.distance(pt2);
            if dist < tolerance {
                let mut is_ins = false;
                for kk in 1..=seq_nodes.len() {
                    if is_ins {
                        break;
                    }
                    let pt = bat::brep_tool_pnt(&seq_nodes[kk - 1]).unwrap_or(DVec3::ZERO);
                    if dist < pt1.distance(pt) {
                        seq_nodes.insert(kk - 1, node2.clone());
                        is_ins = true;
                    }
                }
                if !is_ins {
                    seq_nodes.push(node2.clone());
                }
            }
        }
        // OCCT L3467-3500.
        if !seq_nodes.is_empty() {
            // Remove nodes near to some other from the same edge
            if !list_nodes_same_edge.is_empty() {
                for l_int in &list_nodes_same_edge {
                    let n2 = l_int.clone();
                    let p2 = bat::brep_tool_pnt(&n2).unwrap_or(DVec3::ZERO);
                    let mut k = 1usize;
                    while k <= seq_nodes.len() {
                        let n1 = seq_nodes[k - 1].clone();
                        if !(n1.ptr_id() == n2.ptr_id() && n1.location == n2.location) {
                            let p1 = bat::brep_tool_pnt(&n1).unwrap_or(DVec3::ZERO);
                            if p2.distance(p1) >= pt1.distance(p1) {
                                k += 1;
                                continue;
                            }
                        }
                        seq_nodes.remove(k - 1);
                    }
                }
            }
            // Bind nearest node if at least one exists
            // OCCT L3501-3506.
            if !seq_nodes.is_empty() {
                idx_map_add(&mut node_nearest_node, &node1, seq_nodes[0].clone());
            }
        }
        // OCCT L3508: anInspector.ClearResList().
        an_inspector.clear_res_list();
    }

    // Create new nodes for chained nearest nodes
    // OCCT L3512-3516.
    if node_nearest_node.is_empty() {
        return false;
    }

    // OCCT L3518.
    create_new_nodes(brep, &node_nearest_node, &node_vertices, a_vertex_node, a_node_edges)
}

/// OCCT `NodeVertices.FindFromKey(node)` over the aVertexNode map (the
/// identity when unbound).
fn node_of_map(map: &super::IdxShapeMap, v: &Shape) -> Shape {
    match map.get(&bat::shape_key(v)) {
        Some((_, n)) => n.clone(),
        None => v.clone(),
    }
}

/// OCCT `NCollection_IndexedDataMap::FindIndex(key)` — the 1-based rank of
/// the key.
fn find_index_rank(m: &super::IdxListMap, k: &Shape) -> i32 {
    for (i, (key, _)) in m.iter().enumerate() {
        if *key == bat::shape_key(k) {
            return (i + 1) as i32;
        }
    }
    0
}

// Keep the Arc import referenced (the edge TShape makers live in
// sewing_same_param.rs).
#[allow(unused_imports)]
use Arc as _ArcAnchor;
