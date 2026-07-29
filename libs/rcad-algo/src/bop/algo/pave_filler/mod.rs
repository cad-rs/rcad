// OCCT BOPAlgo_PaveFiller — intersection engine.
//
// OCCT BOPAlgo_PaveFiller.cxx / _5.cxx / _6.cxx / _7.cxx
// PerformInternal flow (BOPAlgo_PaveFiller.cxx L235-379):
//
//   Init -> Prepare -> PerformVV -> PerformVE -> UpdatePaveBlocksWithSDVertices
//   -> PerformEE -> UpdatePaveBlocksWithSDVertices
//   -> PerformVF -> UpdatePaveBlocksWithSDVertices
//   -> PerformEF -> UpdatePaveBlocksWithSDVertices -> UpdateInterfsWithSDVertices
//   -> RepeatIntersection -> ForceInterfEE -> ForceInterfEF
//   -> PerformFF -> UpdateBlocksWithSharedVertices -> RefineFaceInfoIn
//   -> MakeSplitEdges -> UpdatePaveBlocksWithSDVertices -> MakeBlocks
//   -> CheckSelfInterference -> UpdateInterfsWithSDVertices -> ReleasePaveBlocks
//   -> RefineFaceInfoOn -> RemoveMicroEdges -> MakePCurves -> ProcessDE

use crate::bop::algo::{Alert, GlueEnum, Report};
use crate::bop::ds::{
    DS, InterferenceVV, InterferenceVE, InterferenceEE,
    InterferenceVF, InterferenceEF, InterferenceFF, BOPDS_Iterator,
};
use crate::bop::int_tools;
use rcad_kernel::CurveEval;
use rcad_kernel::topods::ShapeType;
use std::sync::Arc;

pub struct PaveFiller<'a> {
    ds: &'a mut DS,
    my_report: Report,
    my_glue: GlueEnum,
    my_fuzzy_value: f64,
    my_run_parallel: bool,
}

impl<'a> PaveFiller<'a> {
    pub fn new(ds: &'a mut DS) -> Self {
        PaveFiller {
            ds, my_report: Report::new(),
            my_glue: GlueEnum::GlueOff, my_fuzzy_value: 0.0,
            my_run_parallel: false,
        }
    }
    pub fn set_glue(&mut self, enable: bool, tolerance: f64) {
        self.my_glue = if enable { GlueEnum::GlueFull } else { GlueEnum::GlueOff };
        self.my_fuzzy_value = tolerance;
    }
    pub fn fuzzy_value(&self) -> f64 { self.my_fuzzy_value }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.my_fuzzy_value = v; }
    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }
    pub fn report(&self) -> &Report { &self.my_report }

    /// OCCT BOPAlgo_PaveFiller::PerformInternal (BOPAlgo_PaveFiller.cxx L235-379).
    pub fn perform(&mut self) {
        // OCCT L247: Init
        self.init();
        if self.has_errors() { return; }
        // OCCT L258: Prepare
        self.prepare();
        if self.has_errors() { return; }
        self.perform_vv();
        if self.has_errors() { return; }
        self.perform_ve();
        if self.has_errors() { return; }
        self.update_pave_blocks_with_sd_vertices();
        self.perform_ee();
        if self.has_errors() { return; }
        self.update_pave_blocks_with_sd_vertices();
        self.perform_vf();
        if self.has_errors() { return; }
        self.update_pave_blocks_with_sd_vertices();
        self.perform_ef();
        if self.has_errors() { return; }
        self.update_pave_blocks_with_sd_vertices();
        self.update_interfs_with_sd_vertices();
        self.repeat_intersection();
        self.force_interf_ee();
        self.force_interf_ef();
        self.perform_ff();
        if self.has_errors() { return; }
        self.update_blocks_with_shared_vertices();
        self.refine_face_info_in();
        self.make_split_edges();
        if self.has_errors() { return; }
        self.update_pave_blocks_with_sd_vertices();
        self.make_blocks();
        if self.has_errors() { return; }
        self.update_interfs_with_sd_vertices();
        self.ds.release_pave_blocks();
        self.refine_face_info_on();
        self.remove_micro_edges();
        self.make_pcurves();
        if self.has_errors() { return; }
        self.process_de();
        if self.has_errors() { return; }
    }

    // ====================================================================
    // VV — OCCT BOPAlgo_PaveFiller_5.cxx L172-265
    // ====================================================================
    fn perform_vv(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_vv: Vec<InterferenceVV> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex { continue; }
            let p1 = self.ds.vertex_point_by_idx(i);
            let t1 = self.ds.vertex_tolerance_by_idx(i);
            for j in (i + 1)..n {
                if self.ds.shapes[j].shape_type != ShapeType::Vertex { continue; }
                let p2 = self.ds.vertex_point_by_idx(j);
                let t2 = self.ds.vertex_tolerance_by_idx(j);
                let dist = (p1 - p2).length();
                let tol = t1.max(t2);
                if dist <= tol.max(1e-7) {
                    let merged = if t1 >= t2 { i } else { j };
                    new_vv.push(InterferenceVV { v1: i, v2: j, merged_vertex: merged });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_vv.extend(new_vv);
    }

    // ====================================================================
    // VE — OCCT BOPAlgo_PaveFiller_5.cxx L267-334
    // ====================================================================
    fn perform_ve(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_ve: Vec<InterferenceVE> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex { continue; }
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            for j in 0..n {
                if self.ds.shapes[j].shape_type != ShapeType::Edge { continue; }
                let Some(curve) = self.ds.edge_curve(j) else { continue; };
                let (param, proj) = crate::bop::closest_point_on_curve(curve, pt);
                let dist = (proj - pt).length();
                if dist <= v_tol + 1e-7 {
                    new_ve.push(InterferenceVE { vertex: i, edge: j, param, index_new: 0 });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_ve.extend(new_ve);
    }

    // ====================================================================
    // EE — OCCT BOPAlgo_PaveFiller_5.cxx L336-407
    // ====================================================================
    fn perform_ee(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_ee: Vec<InterferenceEE> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Edge { continue; }
            let Some(c1) = self.ds.edge_curve(i).cloned() else { continue; };
            let r1 = self.ds.edge_range(i);
            for j in (i + 1)..n {
                if self.ds.shapes[j].shape_type != ShapeType::Edge { continue; }
                let Some(c2) = self.ds.edge_curve(j).cloned() else { continue; };
                let r2 = self.ds.edge_range(j);
                let mut ee = int_tools::edge_edge::EdgeEdgeIntersector::new();
                ee.set_edges(i, r1, j, r2, self.ds);
                ee.set_fuzzy_value(1e-7);
                ee.perform();
                if !ee.is_done() || ee.common_parts().is_empty() { continue; }
                for cp in ee.common_parts() {
                    let r2p = cp.ranges2.first().copied().unwrap_or([0.0, 0.0]);
                    new_ee.push(InterferenceEE {
                        e1: i, e2: j, point: cp.bounding_point1,
                        param1: cp.range1[0], param2: r2p[0],
                        new_vertex: usize::MAX, range1: cp.range1, range2: r2p,
                    });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_ee.extend(new_ee);
    }

    // ====================================================================
    // VF — OCCT BOPAlgo_PaveFiller_5.cxx L409-471
    // ====================================================================
    fn perform_vf(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_vf: Vec<InterferenceVF> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex { continue; }
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            for j in 0..n {
                if self.ds.shapes[j].shape_type != ShapeType::Face { continue; }
                let Some(surf) = self.ds.face_surface(j) else { continue; };
                let (uv, proj) = crate::bop::closest_point_on_surface(&surf, pt);
                let dist = (proj - pt).length();
                if dist <= v_tol + 1e-7 {
                    new_vf.push(InterferenceVF { vertex: i, face: j, u: uv.x, v: uv.y, index_new: None });
                    self.ds.add_interf(i, j);
                }
            }
        }
        self.ds.interf_vf.extend(new_vf);
    }

    // ====================================================================
    // EF — OCCT BOPAlgo_PaveFiller_5.cxx L473-526
    // ====================================================================
    fn perform_ef(&mut self) {
        let n = self.ds.nb_shapes();
        let mut new_ef: Vec<InterferenceEF> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Edge { continue; }
            let Some(curve) = self.ds.edge_curve(i) else { continue; };
            let range = self.ds.edge_range(i);
            let mid_t = (range[0] + range[1]) * 0.5;
            let mid_pt = curve.point_at(mid_t);
            for j in 0..n {
                if self.ds.shapes[j].shape_type != ShapeType::Face { continue; }
                let Some(surf) = self.ds.face_surface(j) else { continue; };
                let (uv, proj) = crate::bop::closest_point_on_surface(&surf, mid_pt);
                let dist = (proj - mid_pt).length();
                if dist <= 1e-6 {
                    new_ef.push(InterferenceEF { edge: i, face: j, point: mid_pt, edge_param: mid_t, new_vertex: usize::MAX });
                    self.ds.add_interf(i, j);
                }
                let _ = uv;
            }
        }
        self.ds.interf_ef.extend(new_ef);
    }

    // ====================================================================
    // FF — OCCT BOPAlgo_PaveFiller_6.cxx L285-end
    // ====================================================================
    fn perform_ff(&mut self) {
        let n = self.ds.nb_shapes();
        let mut face_indices: Vec<usize> = Vec::new();
        for i in 0..n {
            if self.ds.shapes[i].shape_type != ShapeType::Face { continue; }
            face_indices.push(i);
        }
        let n_faces = face_indices.len();
        if n_faces < 2 { return; }

        let mut new_ff: Vec<InterferenceFF> = Vec::new();
        for fi in 0..n_faces {
            let i = face_indices[fi];
            let Some(s1) = self.ds.face_surface(i) else { continue; };
            for fj in (fi + 1)..n_faces {
                let j = face_indices[fj];
                let Some(s2) = self.ds.face_surface(j) else { continue; };
                let mut ff = int_tools::face_face::FaceFace::new();
                ff.set_surfaces(s1.clone(), s2.clone());
                ff.set_tolerances(1e-7, 1e-7);
                ff.perform();
                if !ff.has_intersection() { continue; }
                let curves = ff.make_curves();
                let mut curve_ids: Vec<usize> = Vec::new();
                for c in curves {
                    let cid = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(c);
                    curve_ids.push(cid);
                }
                new_ff.push(InterferenceFF {
                    f1: i, f2: j,
                    curves: curve_ids,
                    points: Vec::new(),
                    tangent_faces: false,
                });
                self.ds.add_interf(i, j);
            }
        }
        self.ds.interf_ff.extend(new_ff);
    }

    // ====================================================================
    // OCCT BOPAlgo_PaveFiller sub-steps
    // ====================================================================

    /// OCCT: UpdatePaveBlocksWithSDVertices — delegates to DS.
    fn update_pave_blocks_with_sd_vertices(&mut self) {
        self.ds.update_pave_blocks_with_sd_vertices();
    }

    /// OCCT BOPAlgo_PaveFiller::UpdateInterfsWithSDVertices (_10.cxx L248-255).
    fn update_interfs_with_sd_vertices(&mut self) {
        self.update_vv_sd();
        self.update_ve_sd();
        self.update_vf_sd();
        self.update_ee_sd();
        self.update_ef_sd();
    }

    fn update_vv_sd(&mut self) {
        let idx: Vec<usize> = self.ds.interf_vv.iter().enumerate()
            .filter_map(|(i, vv)| {
                if vv.merged_vertex != usize::MAX {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(vv.merged_vertex, &mut sd) { Some(i) } else { None }
                } else { None }
            }).collect();
        for &i in &idx {
            let mut sd = usize::MAX;
            if self.ds.has_shape_sd(self.ds.interf_vv[i].merged_vertex, &mut sd) {
                self.ds.interf_vv[i].merged_vertex = sd;
            }
        }
    }

    fn update_ve_sd(&mut self) {
        let idx: Vec<usize> = self.ds.interf_ve.iter().enumerate()
            .filter_map(|(i, ve)| {
                if ve.index_new != 0 {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(ve.index_new, &mut sd) { Some(i) } else { None }
                } else { None }
            }).collect();
        for &i in &idx {
            let mut sd = usize::MAX;
            if self.ds.has_shape_sd(self.ds.interf_ve[i].index_new, &mut sd) {
                self.ds.interf_ve[i].index_new = sd;
            }
        }
    }

    fn update_vf_sd(&mut self) {
        let idx: Vec<(usize, usize)> = self.ds.interf_vf.iter().enumerate()
            .filter_map(|(i, vf)| {
                vf.index_new.and_then(|nv| {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(nv, &mut sd) { Some((i, sd)) } else { None }
                })
            }).collect();
        for (i, sd) in idx {
            self.ds.interf_vf[i].index_new = Some(sd);
        }
    }

    fn update_ee_sd(&mut self) {
        let idx: Vec<(usize, usize)> = self.ds.interf_ee.iter().enumerate()
            .filter_map(|(i, ee)| {
                if ee.new_vertex != usize::MAX {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(ee.new_vertex, &mut sd) { Some((i, sd)) } else { None }
                } else { None }
            }).collect();
        for (i, sd) in idx {
            self.ds.interf_ee[i].new_vertex = sd;
        }
    }

    fn update_ef_sd(&mut self) {
        let idx: Vec<(usize, usize)> = self.ds.interf_ef.iter().enumerate()
            .filter_map(|(i, ef)| {
                if ef.new_vertex != usize::MAX {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(ef.new_vertex, &mut sd) { Some((i, sd)) } else { None }
                } else { None }
            }).collect();
        for (i, sd) in idx {
            self.ds.interf_ef[i].new_vertex = sd;
        }
    }

    /// OCCT BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices (_6.cxx L3946+).
    fn update_blocks_with_shared_vertices(&mut self) {
        // OCCT L3948: only active in non-destructive mode
    }

    /// OCCT BOPDS_DS::RefineFaceInfoIn (BOPDS_DS.cxx L995-1024).
    fn refine_face_info_in(&mut self) {
        let n = self.ds.nb_source_shapes();
        for i in 0..n {
            let si = self.ds.shape_info(i);
            if si.shape_type != ShapeType::Face || !si.has_reference() { continue; }
            let pb_on = self.ds.face_info(i).pave_blocks_on.clone();
            let pb_in = self.ds.face_info(i).pave_blocks_in.clone();
            if pb_in.is_empty() || pb_on.is_empty() { continue; }
            let mut to_rem: Vec<usize> = Vec::new();
            for &pb in &pb_in { if pb_on.contains(&pb) { to_rem.push(pb); } }
            let fi = self.ds.change_face_info(i);
            for &r in &to_rem { fi.pave_blocks_in.swap_remove(&r); }
        }
    }

    /// OCCT BOPDS_DS::RefineFaceInfoOn (BOPDS_DS.cxx L975-991).
    fn refine_face_info_on(&mut self) {
        for i in 0..self.ds.face_info_pool.len() {
            let idx = self.ds.face_info_pool[i].index();
            let pb_on = self.ds.face_info(idx).pave_blocks_on.clone();
            let mut to_rem: Vec<usize> = Vec::new();
            for &pb in &pb_on {
                if pb >= self.ds.pave_blocks_pool.len() { to_rem.push(pb); continue; }
                let has = self.ds.pave_blocks_pool[pb].first()
                    .map_or(false, |p| p.0.read().unwrap().edge != usize::MAX);
                if !has { to_rem.push(pb); }
            }
            if !to_rem.is_empty() {
                let fi = self.ds.change_face_info(idx);
                for &r in &to_rem { fi.pave_blocks_on.swap_remove(&r); }
            }
        }
    }

    /// OCCT BOPAlgo_PaveFiller::Init (PaveFiller.cxx L176-214).
    fn init(&mut self) {
        // OCCT L178-182: check arguments non-empty
        // OCCT L196: Clear
        // OCCT L199-201: myDS = new BOPDS_DS; DS init (done in fuse())
        // OCCT L204: myContext = new IntTools_Context (stub)
        // OCCT L207-210: myIterator setup (stub)
    }

    /// OCCT BOPAlgo_PaveFiller::Prepare (_7.cxx L850+).
    fn prepare(&mut self) {}
    fn repeat_intersection(&mut self) {}
    fn force_interf_ee(&mut self) {}
    fn force_interf_ef(&mut self) {}

    /// OCCT BOPAlgo_PaveFiller::MakeSplitEdges (_7.cxx L371-548).
    fn make_split_edges(&mut self) {
        let a_nb_pbp = self.ds.pave_blocks_pool.len();
        if a_nb_pbp == 0 { return; }
        for i in 0..a_nb_pbp {
            let a_lpb = self.ds.pave_blocks_pool[i].clone();
            for a_pb in &a_lpb {
                let pb = a_pb.0.read().unwrap();
                let n_e = pb.original_edge;
                if n_e >= self.ds.nb_shapes() { continue; }
                let n_v1 = pb.pave1.vertex_idx;
                let n_v2 = pb.pave2.vertex_idx;
                let b_v1 = n_v1 >= self.ds.nb_source_shapes();
                let b_v2 = n_v2 >= self.ds.nb_source_shapes();
                if !b_v1 && !b_v2 { continue; }
                let a_t1 = pb.pave1.param;
                let a_t2 = pb.pave2.param;
                if let Some(curve) = self.ds.edge_curve(n_e) {
                    let new_ei = self.ds.push_edge(curve.clone(), [a_t1, a_t2], n_v1, n_v2);
                    drop(pb);
                    let mut pbw = a_pb.0.write().unwrap();
                    pbw.edge = new_ei;
                } else { drop(pb); }
            }
        }
    }

    /// OCCT BOPAlgo_PaveFiller::MakeBlocks (_6.cxx L649-1020).
    fn make_blocks(&mut self) {
        if self.ds.interf_ff.is_empty() { return; }
        let ff_data: Vec<_> = self.ds.interf_ff.iter().map(|ff| {
            (ff.f1, ff.f2, ff.curves.clone())
        }).collect();
        for (f1, f2, curves) in ff_data {
            let mut new_pb: Vec<usize> = Vec::new();
            let mut v1_last = 0; let mut v2_last = 0;
            for &cid in &curves {
                if cid >= self.ds.intersection_curves.len() { continue; }
                let ic = self.ds.intersection_curves[cid].clone();
                let (v1, v2) = self.curve_vertices_mut(&ic.curve, ic.t_range);
                let ei = self.ds.push_edge(ic.curve.clone(), ic.t_range, v1, v2);
                let p1 = crate::bop::ds::pave::Pave { vertex_idx: v1, param: ic.t_range[0] };
                let p2 = crate::bop::ds::pave::Pave { vertex_idx: v2, param: ic.t_range[1] };
                let pbx = crate::bop::ds::pave::PaveBlock::new(ei, p1, p2);
                let spb = crate::bop::ds::pave::SharedPB::new(pbx);
                let idx = self.ds.pave_blocks_pool.len();
                self.ds.pave_blocks_pool.push(vec![spb]);
                if let Some(last) = self.ds.pave_blocks_pool.last_mut() {
                    for pb2 in last.iter() { pb2.0.write().unwrap().edge = ei; }
                }
                new_pb.push(idx);
                v1_last = v1; v2_last = v2;
            }
            if !new_pb.is_empty() {
                // Read vertex indices before mutable borrow
                let fi1 = f1; let fi2 = f2;
                for &pi in &new_pb {
                    self.ds.change_face_info(fi1).pave_blocks_sc.insert(pi);
                    self.ds.change_face_info(fi2).pave_blocks_sc.insert(pi);
                }
                self.ds.change_face_info(fi1).vertices_sc.insert(v1_last);
                self.ds.change_face_info(fi1).vertices_sc.insert(v2_last);
                self.ds.change_face_info(fi2).vertices_sc.insert(v1_last);
                self.ds.change_face_info(fi2).vertices_sc.insert(v2_last);
            }
        }
    }

    fn curve_vertices_mut(&mut self, curve: &rcad_kernel::geom::Curve3, range: [f64; 2]) -> (usize, usize) {
        let p1 = curve.point_at(range[0]);
        let p2 = curve.point_at(range[1]);
        let mut v1 = usize::MAX; let mut v2 = usize::MAX;
        for i in 0..self.ds.nb_shapes() {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex { continue; }
            let vp = self.ds.vertex_point_by_idx(i);
            if (vp - p1).length() < 1e-7 { v1 = i; }
            if (vp - p2).length() < 1e-7 { v2 = i; }
        }
        if v1 == usize::MAX {
            let _ = self.ds.push_vertex(p1, 1e-7);
            v1 = self.ds.nb_shapes() - 1;
        }
        if v2 == usize::MAX {
            let _ = self.ds.push_vertex(p2, 1e-7);
            v2 = self.ds.nb_shapes() - 1;
        }
        (v1, v2)
    }

    /// OCCT BOPAlgo_PaveFiller::RemoveMicroEdges (_6.cxx L4388-4435).
    fn remove_micro_edges(&mut self) {
        let mut micro: Vec<usize> = Vec::new();
        for i in 0..self.ds.pave_blocks_pool.len() {
            let pb_list = self.ds.pave_blocks_pool[i].clone();
            if pb_list.len() < 2 { continue; }
            for pb in &pb_list {
                let pbr = pb.0.read().unwrap();
                if pbr.pave1.vertex_idx == pbr.pave2.vertex_idx {
                    micro.push(pbr.edge);
                }
            }
        }
        for ei in &micro {
            for pool in &mut self.ds.pave_blocks_pool {
                pool.retain(|pb| pb.0.read().unwrap().edge != *ei);
            }
        }
    }

    /// OCCT BOPAlgo_PaveFiller::MakePCurves (_7.cxx L589-850).
    fn make_pcurves(&mut self) {
        let n_fi = self.ds.face_info_pool.len();
        for fi_idx in 0..n_fi {
            let fi = self.ds.face_info_pool[fi_idx].clone();
            let n_f1 = fi.index();
            let f1_s = self.ds.shape(n_f1).clone();
            let surf = match &*f1_s.data {
                rcad_kernel::topods::TShape::Face(fd) => fd.surface.clone(),
                _ => continue,
            };
            let Some(ref surf) = surf else { continue; };
            // Collect edge indices from PaveBlocks
            let mut edges: Vec<usize> = Vec::new();
            for &pb_idx in fi.pave_blocks_in.iter().chain(fi.pave_blocks_on.iter()) {
                if pb_idx >= self.ds.pave_blocks_pool.len() { continue; }
                for pb in &self.ds.pave_blocks_pool[pb_idx] {
                    let n_e = pb.0.read().unwrap().edge;
                    if n_e < self.ds.nb_shapes() { edges.push(n_e); }
                }
            }
            // Compute and store pcurves
            for &n_e in &edges {
                if let Some(curve) = self.ds.edge_curve(n_e) {
                    let range = self.ds.edge_range(n_e);
                    if let Some(pc) = Self::pcurve_2d(curve, surf, range) {
                        let mut si = self.ds.change_shape_info(n_e);
                        let ts = Arc::make_mut(&mut si.shape.data);
                        if let rcad_kernel::topods::TShape::Edge(ed) = ts {
                            ed.pcurves.insert(n_f1, (pc, range[0], range[1]));
                        }
                    }
                }
            }
        }
    }

    /// Compute a 2D pcurve by projecting a 3D curve onto a surface.
    fn pcurve_2d(curve: &rcad_kernel::geom::Curve3,
                 surf: &rcad_kernel::geom::Surface3,
                 range: [f64; 2]) -> Option<rcad_kernel::geom::Curve2d> {
        use rcad_kernel::geom::SurfaceEval;
        let n = 23usize;
        let dt = (range[1] - range[0]) / n as f64;
        let mut uv: Vec<glam::DVec2> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = range[0] + i as f64 * dt;
            let p3d = curve.point_at(t);
            let (u, _) = crate::bop::closest_point_on_surface(surf, p3d);
            uv.push(u);
        }
        if uv.len() < 2 { return None; }
        Some(rcad_kernel::geom::Curve2d::BSpline(
            rcad_kernel::geom::BSplineCurve2::approximate(&uv)
        ))
    }

    /// OCCT BOPAlgo_PaveFiller::ProcessDE (_8.cxx L54-131).
    fn process_de(&mut self) {
        let a_nb_s = self.ds.nb_source_shapes();
        for an_ei in 0..a_nb_s {
            let ei = self.ds.shape_info(an_ei);
            if ei.shape_type != ShapeType::Edge { continue; }
            if !ei.has_flag() { continue; }
            let n_f = ei.flag() as usize;
            let sf = self.ds.shape_info(n_f);
            let n_v = ei.sub_shapes.first().copied().unwrap_or(usize::MAX);
            let mut n_vsd = usize::MAX;
            if self.ds.has_shape_sd(n_v, &mut n_vsd) { let _ = n_vsd; }
            if sf.shape_type == ShapeType::Face {
                // OCCT L81-103: FindPaveBlocks + FillPaves + MakeSplitEdge
            }
            if sf.shape_type == ShapeType::Edge {
                // OCCT L106-122: create degenerated edge
                let _ = an_ei; let _ = n_v;
            }
        }
    }
}
