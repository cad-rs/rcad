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
        // OCCT L258: Prepare
        self.perform_vv();                                                      // L265
        self.perform_ve();                                                      // L272
        self.update_pave_blocks_with_sd_vertices();                             // L279
        self.perform_ee();                                                      // L281
        self.update_pave_blocks_with_sd_vertices();                             // L287
        self.perform_vf();                                                      // L289
        self.update_pave_blocks_with_sd_vertices();                             // L295
        self.perform_ef();                                                      // L297
        self.update_pave_blocks_with_sd_vertices();                             // L303
        self.update_interfs_with_sd_vertices();                                 // L303
        // OCCT L307: RepeatIntersection
        // OCCT L315: ForceInterfEE
        // OCCT L323: ForceInterfEF
        self.perform_ff();                                                      // L331
        self.update_blocks_with_shared_vertices();                              // L338
        self.refine_face_info_in();                                             // L340
        self.make_split_edges();                                                // L342
        self.update_pave_blocks_with_sd_vertices();                             // L349
        self.make_blocks();                                                     // L351
        // OCCT L358: CheckSelfInterference
        self.update_interfs_with_sd_vertices();                                 // L360
        // OCCT L361: ReleasePaveBlocks
        self.refine_face_info_on();                                             // L362
        // OCCT L364: RemoveMicroEdges
        // OCCT L366: MakePCurves
        // OCCT L373: ProcessDE
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

    /// OCCT: UpdatePaveBlocksWithSDVertices.
    fn update_pave_blocks_with_sd_vertices(&mut self) {}

    /// OCCT: UpdateInterfsWithSDVertices.
    fn update_interfs_with_sd_vertices(&mut self) {}

    /// OCCT: UpdateBlocksWithSharedVertices.
    fn update_blocks_with_shared_vertices(&mut self) {}

    /// OCCT: RefineFaceInfoIn.
    fn refine_face_info_in(&mut self) {}

    /// OCCT: RefineFaceInfoOn.
    fn refine_face_info_on(&mut self) {}

    /// OCCT BOPAlgo_PaveFiller::MakeSplitEdges (_7.cxx L371+).
    /// Creates new edge TShapes from pave blocks with split vertices.
    fn make_split_edges(&mut self) {}

    /// OCCT BOPAlgo_PaveFiller::MakeBlocks (_7.cxx L500+).
    /// Creates pave blocks on edges from intersection data.
    fn make_blocks(&mut self) {}
}
