//! OCCT ShapeUpgrade_UnifySameDomain.cxx L2677-2948 — static
//! `generateSubSeq` (L2685-2752), `MergeEdges` (L2756-2877), `MergeSeq`
//! (L2884-2906), static `CheckSharedVertices` (L2913-2948).  (The
//! `SubSequenceOfEdges` struct lives in `mod.rs`, cxx L2677-2681.)

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation};

use super::merge_sub_seq::{get_line_edge_points, is_merging_possible};
use super::topexp::vertices;
use super::ShapeUpgradeUnifySameDomain;
use super::{map_add, shape_key, IndexedDataMapOfShapeListOfShape, MapOfShape, SubSequenceOfEdges};
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;

impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L2685-2752: generateSubSeq —
    // generates the sub-sequences of edges that can merge into one edge.
    pub(crate) fn generate_sub_seq(
        brep: &mut BRep,
        an_inp_edge_seq: &[Shape],
        seq_of_sub_seq_of_edges: &mut Vec<SubSequenceOfEdges>,
        is_closed: bool,
        the_ang_tol: f64,
        the_lin_tol: f64,
        avoid_edge_vrt: &MapOfShape,
        the_vfmap: &IndexedDataMapOfShapeListOfShape,
    ) {
        // OCCT L2699-2702.
        let mut a_first_point = glam::DVec3::ZERO;
        let mut a_direction_vec = glam::DVec3::ONE;
        let ref_edge = an_inp_edge_seq[0].clone();
        let mut sub_seq = SubSequenceOfEdges {
            seqs_edges: vec![ref_edge.clone()],
            union_edges: Shape::null(),
        };
        seq_of_sub_seq_of_edges.push(sub_seq);

        let mut is_line_direction_ok =
            get_line_edge_points(brep, &ref_edge, &mut a_first_point, &mut a_direction_vec);

        // OCCT L2708-2732.
        for i in 1..an_inp_edge_seq.len() {
            let edge1 = an_inp_edge_seq[i - 1].clone();
            let edge2 = an_inp_edge_seq[i].clone();
            let is_ok = is_merging_possible(
                brep,
                &edge1,
                &edge2,
                the_ang_tol,
                the_lin_tol,
                avoid_edge_vrt,
                is_line_direction_ok,
                a_first_point,
                a_direction_vec,
                the_vfmap,
            );
            if !is_ok {
                let a_sub_seq = SubSequenceOfEdges {
                    seqs_edges: vec![edge2.clone()],
                    union_edges: Shape::null(),
                };
                seq_of_sub_seq_of_edges.push(a_sub_seq);
                is_line_direction_ok =
                    get_line_edge_points(brep, &edge2, &mut a_first_point, &mut a_direction_vec);
            } else {
                // OCCT L2730: append to the last sub-sequence.
                seq_of_sub_seq_of_edges
                    .last_mut()
                    .unwrap()
                    .seqs_edges
                    .push(edge2);
            }
        }

        // OCCT L2733-2751: the first/last chain segment check.
        if is_closed && seq_of_sub_seq_of_edges.len() > 1 {
            let edge1 = an_inp_edge_seq[an_inp_edge_seq.len() - 1].clone();
            let edge2 = an_inp_edge_seq[0].clone();
            if is_merging_possible(
                brep,
                &edge1,
                &edge2,
                the_ang_tol,
                the_lin_tol,
                avoid_edge_vrt,
                false,
                a_first_point,
                a_direction_vec,
                the_vfmap,
            ) {
                // OCCT L2748-2749: merge the first chain into the last and
                // drop the first.
                let first_edges = seq_of_sub_seq_of_edges[0].seqs_edges.clone();
                let last = seq_of_sub_seq_of_edges.last_mut().unwrap();
                last.seqs_edges.extend(first_edges);
                seq_of_sub_seq_of_edges.remove(0);
            }
        }
    }

    // OCCT ShapeUpgrade_UnifySameDomain.cxx L2756-2877: MergeEdges.
    pub(crate) fn merge_edges(
        &mut self,
        brep: &mut BRep,
        seq_edges: &[Shape],
        the_vfmap: &IndexedDataMapOfShapeListOfShape,
        seq_of_sub_seq_of_edges: &mut Vec<SubSequenceOfEdges>,
        non_merg_vrt: &MapOfShape,
    ) -> bool {
        // OCCT L2764-2767: the map V-E and the vertices to avoid.
        let mut a_map_ve = IndexedDataMapOfShapeListOfShape::new();
        let mut vertices_to_avoid = MapOfShape::new();
        let a_nb_e = seq_edges.len();
        for j in 1..=a_nb_e {
            let an_edge = seq_edges[j - 1].clone();
            // OCCT L2773-2784: the forward-oriented vertex iteration.
            let mut oriented_edge = an_edge.clone();
            oriented_edge.orientation = Orientation::Forward;
            for it in crate::shhealing::shape_build::brep_tool::iter_subshapes(
                brep,
                &oriented_edge,
                true,
                true,
            ) {
                let a_v = it;
                if a_v.orientation == Orientation::Forward
                    || a_v.orientation == Orientation::Reversed
                {
                    let key = shape_key(&a_v);
                    match a_map_ve.get_mut(&key) {
                        Some((_, list)) => list.push(an_edge.clone()),
                        None => {
                            a_map_ve.insert(key, (a_v, vec![an_edge.clone()]));
                        }
                    }
                }
            }
        }
        // OCCT L2786: NCollection_MapAlgo::Unite(VerticesToAvoid, NonMergVrt).
        for (k, v) in non_merg_vrt.iter() {
            vertices_to_avoid.insert(*k, v.clone());
        }

        // OCCT L2788-2789: do loop while there are unused edges.
        let mut a_used_edges = MapOfShape::new();

        for i_e in 1..=a_nb_e {
            let edge = seq_edges[i_e - 1].clone();
            // OCCT L2794-2796.
            if !map_add(&mut a_used_edges, &edge) {
                continue;
            }

            // OCCT L2799-2803: make the chain for the unite.
            let mut a_chain: Vec<Shape> = Vec::new();
            a_chain.push(edge.clone());
            let (v0, v1) = vertices(brep, &edge, true);
            let mut v = [v0, v1];

            // OCCT L2805-2843: connect more edges in both directions.
            for j in 0..2 {
                let mut is_added = true;
                while is_added {
                    is_added = false;
                    if v[j].is_null() {
                        break;
                    }
                    let a_le = a_map_ve
                        .get(&shape_key(&v[j]))
                        .map(|(_, l)| l.clone())
                        .unwrap_or_default();
                    for it_l in &a_le {
                        let edge = it_l.clone();
                        if !a_used_edges.contains_key(&shape_key(&edge)) {
                            let (v20, v21) = vertices(brep, &edge, true);
                            let v2 = [v20, v21];
                            // OCCT L2825: the neighboring edge must have V[j]
                            // reversed and located on the opposite end.
                            if v2[1 - j].is_equal(&reversed_copy(&v[j])) {
                                if j == 0 {
                                    a_chain.insert(0, edge.clone());
                                } else {
                                    a_chain.push(edge.clone());
                                }
                                map_add(&mut a_used_edges, &edge);
                                v[j] = v2[j].clone();
                                is_added = true;
                                break;
                            }
                        }
                    }
                }
            }

            // OCCT L2845-2848.
            if a_chain.len() < 2 {
                continue;
            }

            // OCCT L2850-2854.
            let is_closed = super::topexp::occt_is_same_shape(&v[0], &v[1]);

            // OCCT L2856-2858: split the chain by the non-mergeable vertices.
            let mut a_one_seq: Vec<SubSequenceOfEdges> = Vec::new();
            Self::generate_sub_seq(
                brep,
                &a_chain,
                &mut a_one_seq,
                is_closed,
                self.my_ang_tol,
                self.my_lin_tol,
                &vertices_to_avoid,
                the_vfmap,
            );

            // OCCT L2860-2861: put the sub-chains in the result.
            seq_of_sub_seq_of_edges.extend(a_one_seq);
        }

        // OCCT L2864-2875: the MergeSubSeq pass over the sub-sequences.
        for i in 1..=seq_of_sub_seq_of_edges.len() {
            if seq_of_sub_seq_of_edges[i - 1].seqs_edges.len() < 2 {
                continue;
            }
            let mut ue = Shape::null();
            let seqs = seq_of_sub_seq_of_edges[i - 1].seqs_edges.clone();
            if self.merge_sub_seq(brep, &seqs, the_vfmap, &mut ue) {
                seq_of_sub_seq_of_edges[i - 1].union_edges = ue;
            }
        }
        // OCCT L2876.
        true
    }

    // OCCT ShapeUpgrade_UnifySameDomain.cxx L2884-2906: MergeSeq.
    pub(crate) fn merge_seq(
        &mut self,
        brep: &mut BRep,
        seq_edges: &mut Vec<Shape>,
        the_vfmap: &IndexedDataMapOfShapeListOfShape,
        non_merg_vert: &MapOfShape,
    ) -> bool {
        let mut seq_of_subs_seq_of_edges: Vec<SubSequenceOfEdges> = Vec::new();
        // OCCT L2892.
        if self.merge_edges(
            brep,
            seq_edges,
            the_vfmap,
            &mut seq_of_subs_seq_of_edges,
            non_merg_vert,
        ) {
            for i in 1..=seq_of_subs_seq_of_edges.len() {
                if seq_of_subs_seq_of_edges[i - 1].union_edges.is_null() {
                    continue;
                }
                // OCCT L2901: myContext->Merge(SeqsEdges, UnionEdges).
                self.my_context.merge(
                    brep,
                    &seq_of_subs_seq_of_edges[i - 1].seqs_edges,
                    &seq_of_subs_seq_of_edges[i - 1].union_edges,
                );
            }
            // OCCT L2903.
            return true;
        }
        // OCCT L2905.
        false
    }
}

/// OCCT TopoDS_Shape::Reversed() — the reversed-orientation copy (the
/// `V[j].Reversed()` form of cxx L2825).
fn reversed_copy(s: &Shape) -> Shape {
    let mut c = s.clone();
    c.orientation = super::occt_reverse(c.orientation);
    c
}

/// OCCT static CheckSharedVertices (cxx L2913-2948): checks the sequence of
/// edges on the presence of a shared vertex.
pub fn check_shared_vertices(
    brep: &mut BRep,
    the_seq_edges: &[Shape],
    the_map_edges_vertex: &IndexedDataMapOfShapeListOfShape,
    the_map_keep_shape: &MapOfShape,
    the_share_vert_map: &mut MapOfShape,
) {
    let sae = ShapeAnalysisEdge::new();
    let mut seq_vertexes: Vec<Shape> = Vec::new();
    let mut map_vertexes = MapOfShape::new();
    for k in 1..=the_seq_edges.len() {
        let e = the_seq_edges[k - 1].clone();
        let a_v1 = sae.first_vertex(brep, &e);
        let a_v2 = sae.last_vertex(brep, &e);
        // OCCT L2928-2934: add-or-record duplicates.
        if !map_add(&mut map_vertexes, &a_v1) {
            seq_vertexes.push(a_v1);
        }
        if !map_add(&mut map_vertexes, &a_v2) {
            seq_vertexes.push(a_v2);
        }
    }

    for k in 1..=seq_vertexes.len() {
        let a_v = seq_vertexes[k - 1].clone();
        // OCCT L2940-2945.
        let list_edges_v1 = the_map_edges_vertex
            .get(&shape_key(&a_v))
            .map(|(_, l)| l.len())
            .unwrap_or(0);
        if list_edges_v1 > 2 || the_map_keep_shape.contains_key(&shape_key(&a_v)) {
            map_add(the_share_vert_map, &a_v);
        }
    }
}
