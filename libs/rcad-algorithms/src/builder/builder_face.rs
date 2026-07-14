use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::topods::{self, ShapeRef, Orientation, BRepTool};
use rcad_kernel::geom::*;
use crate::bopds::ds::DS;
use crate::builder::types::*;
use crate::builder::{WireEdgeSource, FaceOrigin};

/// OCCT-aligned: BOPAlgo_BuilderFace (BuilderFace.hxx).
/// Self-contained face splitting class.
/// Two lifetimes: 'a = DS (long-lived), 'b = BRep data (may be shorter).
pub(crate) struct BuilderFace<'a, 'b> {
    ds: &'a DS,
    brep: &'b topods::BRep,
    face_refs: &'b [ShapeRef],
    ic_edge_map: &'b [Option<ShapeRef>],
    my_face_refs: &'a std::cell::RefCell<Vec<ShapeRef>>,
    face_idx: usize,
    is_a: bool,
    /// OCCT-aligned: myShapes (aLE edge list).
    shapes: Option<Vec<ShapeRef>>,
    /// OCCT-aligned: myAreas — resulting split-face TShape refs.
    my_areas: Vec<ShapeRef>,
    /// Internal ResultBuilder for TShape creation during perform().
    my_result: crate::builder::result_builder::ResultBuilder,
}

impl<'a, 'b> BuilderFace<'a, 'b> {
    pub fn new(
        ds: &'a DS,
        brep: &'b topods::BRep,
        face_refs: &'b [ShapeRef],
        ic_edge_map: &'b [Option<ShapeRef>],
        my_face_refs: &'a std::cell::RefCell<Vec<ShapeRef>>,
        face_idx: usize,
        is_a: bool,
    ) -> Self {
        Self {
            ds, brep, face_refs, ic_edge_map, my_face_refs, face_idx, is_a,
            shapes: None,
            my_areas: Vec::new(),
            my_result: crate::builder::result_builder::ResultBuilder::new(),
        }
    }

    pub fn set_shapes(&mut self, shapes: Vec<ShapeRef>) {
        self.shapes = Some(shapes);
    }

    pub fn areas(&self) -> &[ShapeRef] {
        &self.my_areas
    }

    /// ✅ OCCT-aligned: Perform (BuilderFace.cxx L117-148).
    /// Converts myShapes (aLE) to WireSegments, then runs the standard
    /// Avoid → Loops → Areas → InternalShapes pipeline.
    pub(crate) fn perform(&mut self, t: &mut topods::BRep) {
        if self.face_idx >= self.face_refs.len() { return; }
        let face_ref = self.face_refs[self.face_idx];
        if face_ref.index >= self.brep.tshapes.len()
            || !matches!(&*self.brep.tshapes[face_ref.index], topods::TShape::Face(_))
        { return; }

        let shapes = match &self.shapes {
            Some(s) => s,
            None => return,
        };

        let e_base = self.ds.vertices.len();
        let face = &self.ds.faces[self.face_idx];
        let mut segments: Vec<WireSegment> = Vec::with_capacity(shapes.len());

        // OCCT-aligned: convert myShapes (aLE) to WireSegments.
        // Each ShapeRef is a DS edge with orientation from the original wire.
        for &sr in shapes {
            let ei = sr.index.saturating_sub(e_base);
            if ei >= self.ds.edges.len() { continue; }
            let edge = &self.ds.edges[ei];

            // OCCT L371-378: INTERNAL → FORWARD + REVERSED
            // OCCT L379-381: normal → as-is
            // For FORWARD: start→end, for REVERSED: end→start
            let (sv, ev) = match sr.orientation {
                Orientation::Forward | Orientation::Internal => (edge.start_vertex, edge.end_vertex),
                Orientation::Reversed => (edge.end_vertex, edge.start_vertex),
                _ => (edge.start_vertex, edge.end_vertex),
            };
            if sv == ev { continue; }

            let rep = self.ds.edge_on_face(ei, self.face_idx);
            let is_closed = self.ds.edges[ei].start_vertex == self.ds.edges[ei].end_vertex;

            segments.push(WireSegment {
                start_vertex: sv,
                end_vertex: ev,
                source: WireEdgeSource::DsEdge(ei),
                orientation: WireOrientation::Forward,
                is_closed_on_face: is_closed,
                second_pcurve: None,
                first_pcurve: rep.map(|r| r.pcurve.clone()),
                t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
            });

            // OCCT L371-378: INTERNAL edges also need REVERSED copy.
            if sr.orientation == Orientation::Internal {
                segments.push(WireSegment {
                    start_vertex: ev,
                    end_vertex: sv,
                    source: WireEdgeSource::DsEdge(ei),
                    orientation: WireOrientation::Reversed,
                    is_closed_on_face: is_closed,
                    second_pcurve: None,
                    first_pcurve: None,
                    t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
                });
            }
        }

        // OCCT L467-494: IN and SC PaveBlock edges are already in aLE
        // (added by build_split_faces L469-486), so they require no special handling.

        if segments.is_empty() { return; }

        let segments_topo = crate::builder::builder_utils_topo_ds::segments_to_topo_ds(
            &segments, self.ds, self.face_idx, self.face_refs, self.ic_edge_map);
        let tool: &dyn BRepTool = self.brep;

        let (avoided_pids, pid_segs) = crate::builder::wire_splitter::perform_shapes_to_avoid_topo(
            &segments_topo, tool);
        let mut avoided = crate::builder::wire_splitter::expand_avoided_pids(&avoided_pids, &pid_segs);

        let wires = crate::builder::wire_path_topo_ds::build_closed_wires(
            &segments_topo, &avoided, tool);

        let in_loop: std::collections::HashSet<usize> = wires.iter().flatten().copied().collect();
        for si in 0..segments_topo.len() {
            if !in_loop.contains(&si) && !avoided.contains(&si) { avoided.insert(si); }
        }
        let internal_wire_groups = crate::builder::wire_path_topo_ds::build_internal_wires(
            &segments_topo, &avoided);

        let wfs = if !wires.is_empty() {
            crate::builder::wire_path_topo_ds::perform_areas(
                &wires, &internal_wire_groups, &segments_topo, tool, self.face_idx, self.ds)
        } else if !avoided.is_empty() {
            vec![WireFace {
                outer_wire: vec![], inner_wires: vec![],
                internal_wires: segments_topo.iter().enumerate()
                    .filter(|(si, _)| avoided.contains(si)).map(|(si, _)| vec![si]).collect(),
            }]
        } else {
            vec![WireFace {
                outer_wire: (0..segments_topo.len()).collect(),
                inner_wires: vec![], internal_wires: vec![],
            }]
        };
        if wfs.is_empty() { return; }

        let mut wfs = wfs;

        crate::builder::wire_path_topo_ds::perform_internal_shapes(
            &mut wfs, &internal_wire_groups, &segments_topo, tool, self.face_idx, face_ref, self.ds);

        // Store results in myAreas
        let e_base = self.ds.vertices.len();
        let ds_ei_to_sr: HashMap<usize, ShapeRef> = segments.iter()
            .filter_map(|seg| match &seg.source {
                WireEdgeSource::DsEdge(ei) => {
                    let idx = e_base + *ei;
                    let ptr_id = std::sync::Arc::as_ptr(&self.brep.tshapes[idx]) as u64;
                    Some((*ei, ShapeRef { ptr_id, index: idx, orientation: Orientation::Forward, location: 0 }))
                }
                _ => None,
            }).collect();
        let sr_index_to_ds_ei: HashMap<usize, usize> = segments.iter()
            .filter_map(|seg| match &seg.source {
                WireEdgeSource::DsEdge(ei) => Some((e_base + *ei, *ei)),
                _ => None,
            }).collect();

        let origin = if self.is_a {
            FaceOrigin::FromA(self.ds.faces[self.face_idx].source_face_idx)
        } else {
            FaceOrigin::FromB(self.ds.faces[self.face_idx].source_face_idx)
        };
        let ic_curves: HashMap<usize, Curve3> = self.ds.intersection_curves.iter()
            .enumerate().map(|(ci, ic)| (ci, ic.curve.clone())).collect();
        self.my_result.ds_edges = Some(std::sync::Arc::new(self.ds.edges.clone()));
        for wf in &wfs {
            self.my_result.emit_wire_face_topods(self.face_idx, wf, &segments_topo, tool,
                &ic_curves, false, origin, &HashMap::new(),
                self.face_refs[self.face_idx], self.ds.faces[self.face_idx].natural_restriction,
                &ds_ei_to_sr, &sr_index_to_ds_ei, self.ds);
            let fi = self.my_result.faces.len().wrapping_sub(1);
            if fi < self.my_result.faces.len() {
                let mut real_face_refs = Vec::new();
                self.my_result.emit_face_topods(t, &mut real_face_refs);
                if let Some(&real_ref) = real_face_refs.last() {
                    if !real_ref.is_null() {
                        self.my_areas.push(real_ref);
                    }
                } else {
                    let last_idx = t.tshapes.len().wrapping_sub(1);
                    if last_idx < t.tshapes.len() {
                        self.my_areas.push(topods::ShapeRef::synthetic(last_idx));
                    }
                }
            }
        }
    }

}
