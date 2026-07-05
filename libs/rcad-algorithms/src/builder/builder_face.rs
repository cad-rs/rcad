use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::topods::{self, ShapeRef, Orientation, BRepTool};
use rcad_kernel::geom::*;
use crate::bopds::ds::DS;
use crate::builder::types::*;
use crate::builder::result_builder::ResultBuilder;
use crate::builder::{collect_face_edge_segments, FaceOrigin, WireEdgeSource};

/// ✅ OCCT-aligned: BOPAlgo_BuilderFace — face splitting algorithm.
///   OCCT has a self-contained class; rcad mirrors it as a struct.
pub(crate) struct BuilderFace<'a> {
    ds: &'a DS,
    brep: &'a topods::BRep,
    face_refs: &'a [ShapeRef],
    ic_edge_map: &'a [Option<ShapeRef>],
    my_face_refs: &'a std::cell::RefCell<Vec<ShapeRef>>,
    face_idx: usize,
    is_a: bool,
}

impl<'a> BuilderFace<'a> {
    pub fn new(
        ds: &'a DS,
        brep: &'a topods::BRep,
        face_refs: &'a [ShapeRef],
        ic_edge_map: &'a [Option<ShapeRef>],
        my_face_refs: &'a std::cell::RefCell<Vec<ShapeRef>>,
        face_idx: usize,
        is_a: bool,
    ) -> Self {
        Self { ds, brep, face_refs, ic_edge_map, my_face_refs, face_idx, is_a }
    }

    /// ✅ OCCT-aligned: Perform — main entry (BuilderFace.cxx Perform).
    pub(crate) fn perform(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        // Guard: check face_refs bounds
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
        drop(segments);

        let origin = if self.is_a {
            FaceOrigin::FromA(self.ds.faces[self.face_idx].source_face_idx)
        } else {
            FaceOrigin::FromB(self.ds.faces[self.face_idx].source_face_idx)
        };
        let ic_curves: HashMap<usize, Curve3> = self.ds.intersection_curves.iter()
            .enumerate().map(|(ci, ic)| (ci, ic.curve.clone())).collect();
        result.ds_edges = Some(std::sync::Arc::new(self.ds.edges.clone()));
        for wf in &wfs {
            result.emit_wire_face_topods(self.face_idx, wf, &segments_topo, tool,
                &ic_curves, false, origin, &HashMap::new(),
                self.face_refs[self.face_idx], self.ds.faces[self.face_idx].natural_restriction,
                &ds_ei_to_sr, &sr_index_to_ds_ei, self.ds);
            // Architecture A1: create TShapes immediately.
            result.emit_face_topods(t, &mut *self.my_face_refs.borrow_mut());
        }
    }

    /// Find pcurve on the given intersection curve for this face.
    fn find_pcurve_for_face(&self, ci: usize) -> Option<Curve2d> {
        if ci >= self.ds.intersection_curves.len() { return None; }
        let ic = &self.ds.intersection_curves[ci];
        // Determine which pcurve corresponds to this face via face_info.curves_sc
        let fi = self.face_idx;
        if self.ds.faces.get(fi).map_or(false, |f| f.face_info.curves_sc.contains(&ci)) {
            ic.pcurve_on_a.clone().or_else(|| ic.pcurve_on_b.clone())
        } else {
            ic.pcurve_on_b.clone().or_else(|| ic.pcurve_on_a.clone())
        }
    }
}
