use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::topods::{self, ShapeRef, Orientation, BRepTool};
use rcad_kernel::geom::*;
use crate::bopds::ds::DS;
use crate::builder::types::*;
use crate::builder::{collect_face_edge_segments, FaceOrigin, WireEdgeSource};

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
    pub(crate) fn perform(&mut self, t: &mut topods::BRep) {
        if self.face_idx >= self.face_refs.len() { return; }
        let face_ref = self.face_refs[self.face_idx];
        if face_ref.index >= self.brep.tshapes.len()
            || !matches!(&*self.brep.tshapes[face_ref.index], topods::TShape::Face(_))
        { return; }

        let pcurve_lookup = |ci: usize| self.find_pcurve_for_face(ci);
        let segments = collect_face_edge_segments(self.ds, self.face_idx, &pcurve_lookup);
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
                self.my_result.emit_face_topods(t, &mut Vec::new());
                let last_idx = t.tshapes.len().wrapping_sub(1);
                if last_idx < t.tshapes.len() {
                    self.my_areas.push(topods::ShapeRef::synthetic(last_idx));
                }
            }
        }
    }

    fn find_pcurve_for_face(&self, ci: usize) -> Option<Curve2d> {
        if ci >= self.ds.intersection_curves.len() { return None; }
        let ic = &self.ds.intersection_curves[ci];
        let fi = self.face_idx;
        if self.ds.faces.get(fi).map_or(false, |f| f.face_info.curves_sc.contains(&ci)) {
            ic.pcurve_on_a.clone().or_else(|| ic.pcurve_on_b.clone())
        } else {
            ic.pcurve_on_b.clone().or_else(|| ic.pcurve_on_a.clone())
        }
    }
}
