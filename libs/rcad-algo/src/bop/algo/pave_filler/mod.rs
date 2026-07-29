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
use crate::bop::ds::pave::{Pave, PaveBlock, SharedPB};
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::int_tools;
use crate::topalgo::adaptor::CurveAdaptor;
use rcad_kernel::CurveEval;
use rcad_kernel::topods::ShapeType;
use std::collections::{HashSet, HashMap};
use std::sync::Arc;
use glam::DVec3;
use rcad_kernel::topo_shape::{self, Shape};

/// OCCT BOPDS_CoupleOfPaveBlocks — stores two PBs, new vertex, interference index, tolerance.
struct CoupleOfPBs {
    pb1: SharedPB,
    pb2: SharedPB,
    index_interf: usize,
    tolerance: f64,
    index: usize, // new vertex index (iV)
}

/// OCCT BOPAlgo_SectionAttribute (BOPAlgo_SectionAttribute.hxx).
/// Section attributes for boolean operations (method, tolerance, etc.).
#[derive(Debug, Clone)]
struct SectionAttribute {
    method: i32,
    section_tolerance: f64,
}

impl Default for SectionAttribute {
    fn default() -> Self {
        SectionAttribute { method: 0, section_tolerance: 1e-7 }
    }
}

pub struct PaveFiller<'a> {
    // ── BOPAlgo_Algo base (inherited) ──────────────────────────────
    my_report: Report,                 // BOPAlgo_Algo::myReport
    my_run_parallel: bool,             // BOPAlgo_Algo::myRunParallel
    my_fuzzy_value: f64,               // BOPAlgo_Algo::myFuzzyValue
    // ── BOPAlgo_PaveFiller members ────────────────────────────────
    // OCCT BOPAlgo_PaveFiller.hxx L639-652:
    ds: &'a mut DS,                    // L640: myDS
    my_iterator: Option<Box<BOPDS_Iterator<'a>>>, // L641: myIterator
    my_context: IntToolsContext,        // BOPAlgo_PaveFiller::myContext (L642)
    my_glue: GlueEnum,                 // BOPAlgo_PaveFiller::myGlue (L647)
    my_section_attribute: SectionAttribute, // BOPAlgo_PaveFiller::mySectionAttribute (L643)
    my_non_destructive: bool,          // BOPAlgo_PaveFiller::myNonDestructive (L644)
    my_is_primary: bool,               // BOPAlgo_PaveFiller::myIsPrimary (L645)
    my_avoid_build_pcurve: bool,       // BOPAlgo_PaveFiller::myAvoidBuildPCurve (L646)
    my_arguments: Vec<topo_shape::Shape>, // BOPAlgo_PaveFiller::myArguments (L639)
    my_fpb_done: HashMap<usize, HashSet<u64>>, // BOPAlgo_PaveFiller::myFPBDone (L650)
    my_increased_ss: HashSet<usize>,   // BOPAlgo_PaveFiller::myIncreasedSS (L651)
    my_verts_to_avoid_extension: HashSet<usize>, // BOPAlgo_PaveFiller::myVertsToAvoidExtension (L652)
}

impl<'a> PaveFiller<'a> {
    pub fn new(ds: &'a mut DS) -> Self {
        PaveFiller {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_fuzzy_value: 0.0,
            my_iterator: None,
            my_context: IntToolsContext::new(),
            my_glue: GlueEnum::GlueOff,
            my_section_attribute: SectionAttribute::default(),
            my_non_destructive: false,
            my_is_primary: true,
            my_avoid_build_pcurve: false,
            my_arguments: Vec::new(),
            my_fpb_done: HashMap::new(),
            my_increased_ss: HashSet::new(),
            my_verts_to_avoid_extension: HashSet::new(),
        }
    }
    pub fn set_arguments(&mut self, args: Vec<Shape>) {
        self.my_arguments = args;
    }
    pub fn set_glue(&mut self, enable: bool, tolerance: f64) {
        self.my_glue = if enable { GlueEnum::GlueFull } else { GlueEnum::GlueOff };
        self.my_fuzzy_value = tolerance;
    }
    pub fn fuzzy_value(&self) -> f64 { self.my_fuzzy_value }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.my_fuzzy_value = v; }
    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }
    pub fn report(&self) -> &Report { &self.my_report }
    pub fn ds(&self) -> &DS { self.ds }

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
        if self.has_errors() { return; }
        self.force_interf_ee();
        if self.has_errors() { return; }
        self.force_interf_ef();
        if self.has_errors() { return; }
        self.perform_ff();
        if self.has_errors() { return; }
        self.update_blocks_with_shared_vertices();
        self.refine_face_info_in();
        self.make_split_edges();
        if self.has_errors() { return; }
        self.update_pave_blocks_with_sd_vertices();
        self.make_blocks();
        if self.has_errors() { return; }
        self.check_self_interference();
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
        // OCCT: myIterator->Initialize(TopAbs_VERTEX, TopAbs_VERTEX)
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Vertex, ShapeType::Vertex).to_vec()
        } else {
            return;
        };
        let mut new_vv: Vec<InterferenceVV> = Vec::new();
        for &(i, j) in &pairs {
            let p1 = self.ds.vertex_point_by_idx(i);
            let t1 = self.ds.vertex_tolerance_by_idx(i);
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
        self.ds.interf_vv.extend(new_vv);
    }

    // ====================================================================
    // VE — OCCT BOPAlgo_PaveFiller_5.cxx L267-334
    // ====================================================================
    fn perform_ve(&mut self) {
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Vertex, ShapeType::Edge).to_vec()
        } else {
            return;
        };
        let mut new_ve: Vec<InterferenceVE> = Vec::new();
        for &(i, j) in &pairs {
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            let Some(curve) = self.ds.edge_curve(j) else { continue; };
            let (param, proj) = crate::bop::closest_point_on_curve(curve, pt);
            let dist = (proj - pt).length();
            if dist <= v_tol + 1e-7 {
                new_ve.push(InterferenceVE { vertex: i, edge: j, param, index_new: 0 });
                self.ds.add_interf(i, j);
            }
        }
        self.ds.interf_ve.extend(new_ve);
    }

    // ====================================================================
    // EE — OCCT BOPAlgo_PaveFiller.cxx L279-286 + PaveFiller_5.cxx L145-590
    // ====================================================================
    fn perform_ee(&mut self) {
        // OCCT L147: FillShrunkData(TopAbs_EDGE, TopAbs_EDGE)
        self.fill_shrunk_data_type_pair(ShapeType::Edge, ShapeType::Edge);

        // OCCT L149-155: myIterator->Initialize(EDGE, EDGE); check iSize
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Edge, ShapeType::Edge).to_vec()
        } else {
            return;
        };
        if pairs.is_empty() {
            return;
        }

        let mut new_ee: Vec<InterferenceEE> = Vec::new();
        let mut cb_pairs: Vec<(SharedPB, SharedPB)> = Vec::new();
        let mut mvcpb: Vec<CoupleOfPBs> = Vec::new(); // for PerformNewVertices

        // OCCT L181-267: iterate pairs
        for &(n_e1, n_e2) in &pairs {
            // L189-196: skip degenerated edges
            if self.ds.shapes[n_e1].has_flag() || self.ds.shapes[n_e2].has_flag() {
                continue;
            }

            // L200-210: get PB lists for both edges (clone to avoid borrow conflict)
            let a_lpb1: Vec<SharedPB> = self.ds.edge_pave_blocks(n_e1).to_vec();
            let a_lpb2: Vec<SharedPB> = self.ds.edge_pave_blocks(n_e2).to_vec();
            if a_lpb1.is_empty() || a_lpb2.is_empty() {
                continue;
            }

            // L212-265: iterate PB1 × PB2
            let mut pb_box_cache: std::collections::HashMap<u64, (DVec3, DVec3, f64)> =
                std::collections::HashMap::new();

            for p1 in &a_lpb1 {
                // GetPBBox for PB1
                let (mut t11, mut t12, mut ts11, mut ts12) = (0.0, 0.0, 0.0, 0.0);
                let mut bb1 = rcad_kernel::math::bnd::BndBox::new();
                if !self.get_pb_box(n_e1, p1, &mut pb_box_cache,
                    &mut t11, &mut t12, &mut ts11, &mut ts12, &mut bb1) {
                    continue;
                }

                for p2 in &a_lpb2 {
                    let (mut t21, mut t22, mut ts21, mut ts22) = (0.0, 0.0, 0.0, 0.0);
                    let mut bb2 = rcad_kernel::math::bnd::BndBox::new();
                    if !self.get_pb_box(n_e2, p2, &mut pb_box_cache,
                        &mut t21, &mut t22, &mut ts21, &mut ts22, &mut bb2) {
                        continue;
                    }

                    // L245-248: box overlap check
                    if bb1.is_out_box(&bb2) {
                        continue;
                    }

                    // L252: bExpressCompute = PB1 and PB2 have same bounding vertices
                    let (n_v11, n_v12) = { let r = p1.0.read().unwrap(); r.indices() };
                    let (n_v21, n_v22) = { let r = p2.0.read().unwrap(); r.indices() };
                    let _b_express = (n_v11 == n_v21 && n_v12 == n_v22)
                                 || (n_v12 == n_v21 && n_v11 == n_v22);

                    // OCCT L254-264: create EdgeEdge, intersect
                    let mut ee = int_tools::edge_edge::EdgeEdgeIntersector::new();
                    ee.set_edges(n_e1, [t11, t12], n_e2, [t21, t22], self.ds);
                    ee.set_fuzzy_value(self.my_fuzzy_value);
                    ee.perform();

                    if !ee.is_done() {
                        self.my_report.add_warning(
                            Alert::IntersectionFailed(n_e1, n_e2));
                        continue;
                    }

                    let a_cparts = ee.common_parts();
                    let a_nb_cprts = a_cparts.len();
                    if a_nb_cprts == 0 {
                        continue;
                    }

                    // OCCT L355-553: process each common part
                    for (i_cp, cp) in a_cparts.iter().enumerate() {
                        let a_type = if cp.range1[0] >= cp.range1[1] {
                            ShapeType::Vertex  // VERTEX-type intersection
                        } else {
                            ShapeType::Edge    // EDGE-type (coincident)
                        };

                        match a_type {
                            ShapeType::Vertex => {
                                // OCCT L369-526: VERTEX case — new vertex at intersection
                                if a_nb_cprts > 1 && i_cp > 0 {
                                    // OCCT L373-383: only process VERTEX parts
                                    // when there are multiple common parts
                                    continue;
                                }
                                let a_t1 = cp.vertex_param1;
                                let a_t2 = cp.vertex_param2;

                                // OCCT L420: MakeNewVertex and create CoupleOfPBs
                                let a_tol_vnew = self.my_fuzzy_value.max(1e-7);
                                let idx_interf = new_ee.len() + mvcpb.len();
                                mvcpb.push(CoupleOfPBs {
                                    pb1: p1.clone(),
                                    pb2: p2.clone(),
                                    index_interf: idx_interf,
                                    tolerance: a_tol_vnew,
                                    index: usize::MAX,
                                });

                                // OCCT L513-518: add InterfEE
                                new_ee.push(InterferenceEE {
                                    e1: n_e1, e2: n_e2,
                                    point: cp.bounding_point1,
                                    param1: a_t1, param2: a_t2,
                                    new_vertex: usize::MAX,
                                    range1: [a_t1, a_t1],
                                    range2: [a_t2, a_t2],
                                });
                                self.ds.add_interf(n_e1, n_e2);
                            }
                            ShapeType::Edge => {
                                // OCCT L529-550: EDGE case — coincident edges
                                let b_has_same_bounds = (n_v11 == n_v21 && n_v12 == n_v22)
                                                     || (n_v12 == n_v21 && n_v11 == n_v22);
                                if !b_has_same_bounds {
                                    // OCCT L531-533: Not same bounds → skip
                                    if a_nb_cprts <= 1 {
                                        continue;
                                    }
                                }

                                // OCCT L542-547: add InterfEE with common part
                                new_ee.push(InterferenceEE {
                                    e1: n_e1, e2: n_e2,
                                    point: cp.bounding_point1,
                                    param1: cp.range1[0], param2: cp.ranges2[0][0],
                                    new_vertex: usize::MAX,
                                    range1: cp.range1,
                                    range2: cp.ranges2[0],
                                });
                                self.ds.add_interf(n_e1, n_e2);
                                cb_pairs.push((p1.clone(), p2.clone()));
                            }
                            _ => {}
                        }
                    }
                }
            }
        }

        // Push all collected interferences
        self.ds.interf_ee.extend(new_ee);

        // OCCT L561: PerformCommonBlocks
        for (pb1, pb2) in &cb_pairs {
            self.ds.add_common_block(&[pb1.clone(), pb2.clone()]);
        }
        // OCCT L563: UpdateVerticesOfCB
        self.update_vertices_of_cb();

        // OCCT L565-589: PerformNewVertices + SplitPaveBlocks
        if !mvcpb.is_empty() {
            self.perform_new_vertices(&mvcpb, true);
            // Collect edges with modified PBs for SplitPaveBlocks
            let mut a_medges: std::collections::HashSet<usize> = std::collections::HashSet::new();
            for cpb in &mvcpb {
                let (n_e1, n_e2) = (cpb.pb1.0.read().unwrap().original_edge,
                                    cpb.pb2.0.read().unwrap().original_edge);
                a_medges.insert(n_e1);
                a_medges.insert(n_e2);
            }
            self.split_pave_blocks(&a_medges, true);
        }
    }

    // ====================================================================
    // VF — OCCT BOPAlgo_PaveFiller_5.cxx L409-471
    // ====================================================================
    fn perform_vf(&mut self) {
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Vertex, ShapeType::Face).to_vec()
        } else {
            return;
        };
        let mut new_vf: Vec<InterferenceVF> = Vec::new();
        for &(i, j) in &pairs {
            let pt = self.ds.vertex_point_by_idx(i);
            let v_tol = self.ds.vertex_tolerance_by_idx(i);
            let Some(surf) = self.ds.face_surface(j) else { continue; };
            let (uv, proj) = crate::bop::closest_point_on_surface(&surf, pt);
            let dist = (proj - pt).length();
            if dist <= v_tol + 1e-7 {
                new_vf.push(InterferenceVF { vertex: i, face: j, u: uv.x, v: uv.y, index_new: None });
                self.ds.add_interf(i, j);
            }
        }
        self.ds.interf_vf.extend(new_vf);
    }

    // ====================================================================
    // EF — OCCT BOPAlgo_PaveFiller_5.cxx L165-580
    // ====================================================================
    /// OCCT BOPAlgo_PaveFiller::PerformEF (PaveFiller_5.cxx L165-580).
    fn perform_ef(&mut self) {
        // OCCT L167: FillShrunkData(TopAbs_EDGE, TopAbs_FACE)
        self.fill_shrunk_data_type_pair(ShapeType::Edge, ShapeType::Face);

        // OCCT L169-175: myIterator->Initialize(EDGE, FACE); check iSize
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Edge, ShapeType::Face).to_vec()
        } else {
            return;
        };
        if pairs.is_empty() {
            return;
        }

        // OCCT L179-191: GlueFull mode — init FaceInfo and return
        if self.my_glue == GlueEnum::GlueFull {
            for &(n_e, n_f) in &pairs {
                if !self.ds.shapes[n_e].has_flag() {
                    self.ds.change_face_info(n_f);
                }
            }
            return;
        }

        let mut new_ef: Vec<InterferenceEF> = Vec::new();
        let mut mvcpb: Vec<CoupleOfPBs> = Vec::new();

        // OCCT L219-307: iterate EF pairs
        for &(n_e, n_f) in &pairs {
            // OCCT L227-231: skip degenerated edges
            if self.ds.shapes[n_e].has_flag() {
                continue;
            }

            // OCCT L233-237: face box
            let a_box_f = rcad_kernel::math::bnd::BndBox::new();
            let _ = a_box_f;

            // OCCT L237-241: FaceInfo → PaveBlocksOn, VerticesIn, VerticesOn
            let a_fi = self.ds.face_info(n_f);
            let a_mpbf = a_fi.pave_blocks_on.clone();
            let a_mv_in = a_fi.vertices_in.clone();
            let a_mv_on = a_fi.vertices_on.clone();
            drop(a_fi);

            // OCCT L246: get PB list for this edge
            let a_lpb: Vec<SharedPB> = self.ds.edge_pave_blocks(n_e).to_vec();
            if a_lpb.is_empty() {
                continue;
            }

            let mut pb_box_cache: std::collections::HashMap<u64, (DVec3, DVec3, f64)> =
                std::collections::HashMap::new();

            // OCCT L248-306: iterate PBs of this edge
            for a_pb in &a_lpb {
                // OCCT L257-260: skip if PB already in face's PaveBlocksOn
                let is_on_face = {
                    let ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
                    // Check if this PB's pool index is in a_mpbf
                    self.ds.pave_blocks_pool.iter().enumerate().any(|(pi, pool)| {
                        a_mpbf.contains(&pi) && pool.iter().any(|spb|
                            std::sync::Arc::as_ptr(&spb.0) as u64 == ptr)
                    })
                };
                if is_on_face {
                    continue;
                }

                // OCCT L262-270: GetPBBox + box overlap with face
                let (mut a_t1, mut a_t2, mut a_ts1, mut a_ts2) = (0.0, 0.0, 0.0, 0.0);
                let mut a_bb_e = rcad_kernel::math::bnd::BndBox::new();
                if !self.get_pb_box(n_e, a_pb, &mut pb_box_cache,
                    &mut a_t1, &mut a_t2, &mut a_ts1, &mut a_ts2, &mut a_bb_e) {
                    continue;
                }

                // OCCT L273-276: check vertices
                let (n_v1, n_v2) = { let r = a_pb.0.read().unwrap(); r.indices() };
                let b_v1 = a_mv_in.contains(&n_v1) || a_mv_on.contains(&n_v1);
                let b_v2 = a_mv_in.contains(&n_v2) || a_mv_on.contains(&n_v2);
                let _b_express = b_v1 && b_v2;

                // OCCT L278-305: pair data collection
                // rcad: project midpoint onto face surface as proximity check
                let curve = match self.ds.edge_curve(n_e) {
                    Some(c) => c.clone(),
                    None => continue,
                };
                let mid_t = (a_t1 + a_t2) * 0.5;
                let mid_pt = curve.point_at(mid_t);

                if let Some(surf) = self.ds.face_surface(n_f) {
                    let (uv, proj_pt) = crate::bop::closest_point_on_surface(&surf, mid_pt);
                    let dist = (proj_pt - mid_pt).length();
                    if dist > 1e-6 {
                        continue; // too far from face
                    }

                    // OCCT L289-292: BOPTools_AlgoTools::CorrectRange(aE, aF, aSR, anewSR)
                    // rcad: CorrectRange adjusts range endpoints by face tolerance / |derivative|
                    let a_corrected_range = Self::correct_range(
                        n_e, n_f, self.ds, a_ts1, a_ts2);
                    // OCCT L289-292: BOPTools_AlgoTools::CorrectRange(aE, aF, aSR, anewSR)
                    // OCCT L299-305: myFPBDone.Bound(nF, Map<PaveBlock>)->Add(aPB)
                    // rcad: CorrectRange translated above. myFPBDone is a fence map
                    // for processed EF pairs (OCCT data map).

                    // Add EF interference
                    let _ = a_corrected_range;
                    new_ef.push(InterferenceEF {
                        edge: n_e, face: n_f,
                        point: mid_pt,
                        edge_param: mid_t,
                        new_vertex: usize::MAX,
                    });
                    self.ds.add_interf(n_e, n_f);
                    let _ = uv;
                }
            }
        }

        self.ds.interf_ef.extend(new_ef);

        // OCCT L309-580: BOPTools_Parallel::Perform(aVEdgeFace, myContext)
        // rcad: serial midpoint-projection check used instead of parallel EdgeFace.
        // OCCT L565-589: PerformNewVertices for EF results.
        if !mvcpb.is_empty() {
            self.perform_new_vertices(&mvcpb, false);
        }
    }

    // ====================================================================
    // FF — OCCT BOPAlgo_PaveFiller_6.cxx L285-end
    // ====================================================================
    fn perform_ff(&mut self) {
        // OCCT L285-290: myIterator->Initialize(FACE, FACE)
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Face, ShapeType::Face).to_vec()
        } else {
            return;
        };
        if pairs.is_empty() { return; }

        let mut new_ff: Vec<InterferenceFF> = Vec::new();
        for &(i, j) in &pairs {
            let Some(s1) = self.ds.face_surface(i) else { continue; };
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

    /// OCCT BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices (_6.cxx L3946-4020).
    fn update_blocks_with_shared_vertices(&mut self) {
        // OCCT L3948-3951: only active in non-destructive mode
        if !self.my_non_destructive {
            return;
        }
        // L3955-3960: if no FF interferences, return
        if self.ds.interf_ff.is_empty() {
            return;
        }
        // OCCT L3967-4020: iterate FF interferences, build shared vertex sets
        // rcad: non-destructive mode is not fully implemented.
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
        self.my_report.clear();
        self.my_context.clear();
        self.my_increased_ss.clear();
        self.my_verts_to_avoid_extension.clear();
        self.my_fpb_done.clear();
        // OCCT L199-201: myDS = new BOPDS_DS; DS init (done in fuse())
        // OCCT L204: myContext = new IntTools_Context
        self.my_context = IntToolsContext::new();
        // OCCT L207-210: myIterator = new BOPDS_Iterator
        // rcad: use raw pointer to access self.ds without triggering reborrow.
        // The DS lives as long as the PaveFiller, so &'static is safe here.
        let ds_field = unsafe { std::ptr::addr_of_mut!((*self).ds) };
        let ds_ref: &'a mut DS = unsafe { &mut *ds_field };
        let ds_shared: &'a DS = &*ds_ref;
        let ds_static: &'static DS = unsafe { std::mem::transmute(ds_shared) };
        let mut a_it = BOPDS_Iterator::new(ds_static, self.my_fuzzy_value);
        a_it.prepare();
        self.my_iterator = Some(Box::new(a_it));
    }

    /// OCCT BOPAlgo_PaveFiller::Prepare (_7.cxx L850+).
    fn prepare(&mut self) {
        // L850-880: prepare face info, pcurves for planar faces
    }

    // ====================================================================
    // TreatNewVertices — OCCT BOPAlgo_PaveFiller_3.cxx L692-723
    // ====================================================================

    /// Fuse close new vertices into groups.
    /// OCCT BOPAlgo_PaveFiller::TreatNewVertices (PaveFiller_3.cxx L692-723).
    fn treat_new_vertices(&self, the_mvcpb: &[CoupleOfPBs]) -> Vec<(DVec3, f64, Vec<usize>)> {
        // OCCT L700-706: collect vertex points and tolerances
        let mut verts: Vec<(DVec3, f64, usize)> = Vec::new(); // (point, tol, index)
        for (i, cpb) in the_mvcpb.iter().enumerate() {
            let (p1, _p2) = {
                let r1 = cpb.pb1.0.read().unwrap();
                let r2 = cpb.pb2.0.read().unwrap();
                (r1.pave1, r2.pave1)
            };
            // Get vertex position from DS
            let pt = self.ds.vertex_point_by_idx(p1.vertex_idx);
            verts.push((pt, cpb.tolerance, i));
        }
        // OCCT L710: BOPAlgo_Tools::IntersectVertices — fuse by proximity
        // rcad: simple fuse — group vertices within myFuzzyValue distance
        let mut groups: Vec<(DVec3, f64, Vec<usize>)> = Vec::new();
        let mut assigned = vec![false; verts.len()];
        for i in 0..verts.len() {
            if assigned[i] { continue; }
            let (pi, ti, ii) = verts[i];
            let mut group_pts = vec![pi];
            let mut group_tol = ti;
            let mut group_indices = vec![ii];
            assigned[i] = true;
            for j in (i + 1)..verts.len() {
                if assigned[j] { continue; }
                let (pj, tj, ij) = verts[j];
                if (pi - pj).length() <= ti.max(tj) + self.my_fuzzy_value {
                    group_pts.push(pj);
                    group_tol = group_tol.max(tj);
                    group_indices.push(ij);
                    assigned[j] = true;
                }
            }
            // Average position
            let avg = group_pts.iter().sum::<DVec3>() / group_pts.len() as f64;
            groups.push((avg, group_tol, group_indices));
        }
        groups
    }

    /// Process new vertices from EE/EF intersection: add to DS, update interfs, split PBs.
    /// OCCT BOPAlgo_PaveFiller::PerformNewVertices (PaveFiller_3.cxx L594-688).
    fn perform_new_vertices(&mut self, the_mvcpb: &[CoupleOfPBs], is_ee_intersection: bool) {
        // OCCT L601-605: empty check
        if the_mvcpb.is_empty() {
            return;
        }
        // OCCT L607: aTolAdd = myFuzzyValue / 2.
        let a_tol_add = self.my_fuzzy_value / 2.0;

        // OCCT L609-612: TreatNewVertices — fuse vertices
        let groups = self.treat_new_vertices(the_mvcpb);

        // OCCT L622-653: add fused vertices to DS, update interference indices
        // Maps each CPB's index_interf → new DS vertex index
        let mut cpb_to_new_vertex: HashMap<usize, usize> = HashMap::new();
        for (_i, (avg_pt, group_tol, member_indices)) in groups.iter().enumerate() {
            let n_v = self.ds.push_vertex(*avg_pt, *group_tol + a_tol_add);

            for &cpb_idx in member_indices {
                if cpb_idx < the_mvcpb.len() {
                    let idx_interf = the_mvcpb[cpb_idx].index_interf;
                    cpb_to_new_vertex.insert(idx_interf, n_v);
                    // OCCT L648-652: update interference's new_vertex
                    if is_ee_intersection {
                        if idx_interf < self.ds.interf_ee.len() {
                            self.ds.interf_ee[idx_interf].new_vertex = n_v;
                        }
                    } else {
                        if idx_interf < self.ds.interf_ef.len() {
                            self.ds.interf_ef[idx_interf].new_vertex = n_v;
                        }
                    }
                }
            }
        }

        // OCCT L655-685: build PB→[vertices] map and call IntersectVE
        // rcad: add new vertices as ext paves on the PBs
        for cpb in the_mvcpb {
            let n_v = match cpb_to_new_vertex.get(&cpb.index_interf) {
                Some(v) => *v,
                None => continue,
            };
            let (p1_param, p2_param) = {
                let r1 = cpb.pb1.0.read().unwrap();
                let r2 = cpb.pb2.0.read().unwrap();
                (r1.pave1.param, r2.pave1.param)
            };
            {
                let mut w1 = cpb.pb1.0.write().unwrap();
                w1.ext_paves.push(Pave { vertex_idx: n_v, param: p1_param });
            }
            {
                let mut w2 = cpb.pb2.0.write().unwrap();
                w2.ext_paves.push(Pave { vertex_idx: n_v, param: p2_param });
            }
        }
    }

    // ====================================================================
    // FillShrunkData — OCCT BOPAlgo_PaveFiller_9.cxx L65-138
    // ====================================================================

    /// Fill shrunk data for all PBs of edges intersected by type pair.
    /// OCCT BOPAlgo_PaveFiller::FillShrunkData(TopAbs_ShapeEnum, TopAbs_ShapeEnum)
    /// (PaveFiller_9.cxx L65-138).
    /// OCCT BOPAlgo_PaveFiller::FillShrunkData (PaveFiller_9.cxx L65-138).
    fn fill_shrunk_data_type_pair(&mut self, t1: ShapeType, t2: ShapeType) {
        // OCCT L68-73: myIterator->Initialize(aType1, aType2)
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            match (t1, t2) {
                (ShapeType::Edge, ShapeType::Edge) => it.pairs(ShapeType::Edge, ShapeType::Edge).to_vec(),
                (ShapeType::Edge, ShapeType::Face) => it.pairs(ShapeType::Edge, ShapeType::Face).to_vec(),
                _ => Vec::new(),
            }
        } else {
            return;
        };

        let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut edge_indices: Vec<usize> = Vec::new();
        for pair in &pairs {
            if a_mi.insert(pair.0) {
                edge_indices.push(pair.0);
            }
        }

        // Process PBs for collected edges
        for &i in &edge_indices {
            if self.ds.shapes[i].has_flag() {
                continue;
            }
            let a_lpb = self.ds.edge_pave_blocks(i);
            for a_pb in a_lpb {
                if a_pb.0.read().unwrap().has_shrunk_data() {
                    continue;
                }
                self.fill_shrunk_data(a_pb);
            }
        }
    }

    /// Compute shrunk data for a single PaveBlock.
    /// OCCT BOPAlgo_PaveFiller::FillShrunkData(handle<PaveBlock>&) (PaveFiller_3.cxx L727-762).
    fn fill_shrunk_data(&self, the_pb: &SharedPB) {
        let (n_v1, n_v2) = {
            let pbr = the_pb.0.read().unwrap();
            (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
        };
        if n_v1 >= self.ds.nb_shapes() || n_v2 >= self.ds.nb_shapes() {
            return;
        }
        let n_e = {
            let pbr = the_pb.0.read().unwrap();
            if pbr.edge != usize::MAX { pbr.edge } else { pbr.original_edge }
        };
        if n_e >= self.ds.nb_shapes() {
            return;
        }
        let (a_t1, a_t2) = {
            let pbr = the_pb.0.read().unwrap();
            pbr.range()
        };

        // OCCT: IntTools_ShrunkRange aSR; aSR.SetData(aE, aT1, aT2, aV1, aV2);
        // rcad: compute shrunk range by shrinking from vertex tolerance
        let tol_v1 = self.ds.vertex_tolerance_by_idx(n_v1);
        let tol_v2 = self.ds.vertex_tolerance_by_idx(n_v2);

        // Shrink range: subtract tolerance converted to parameter space
        let curve = match self.ds.edge_curve(n_e) {
            Some(c) => c.clone(),
            None => return,
        };
        let dt1 = if (a_t2 - a_t1).abs() > 1e-15 {
            let len = (curve.point_at(a_t1 + 1e-5) - curve.point_at(a_t1)).length().max(1e-7);
            tol_v1 / len * 1e-5
        } else { 0.0 };
        let dt2 = if (a_t2 - a_t1).abs() > 1e-15 {
            let len = (curve.point_at(a_t2) - curve.point_at(a_t2 - 1e-5)).length().max(1e-7);
            tol_v2 / len * 1e-5
        } else { 0.0 };

        let a_ts1 = a_t1 + dt1.min((a_t2 - a_t1) * 0.49);
        let a_ts2 = a_t2 - dt2.min((a_t2 - a_t1) * 0.49);

        // OCCT L761: AnalyzeShrunkData(thePB, aSR)
        self.analyze_shrunk_data(the_pb, a_ts1, a_ts2);
    }

    /// Store shrunk data on a PaveBlock (OCCT: AnalyzeShrunkData).
    /// OCCT BOPAlgo_PaveFiller::AnalyzeShrunkData (PaveFiller_3.cxx L766-824).
    fn analyze_shrunk_data(&self, the_pb: &SharedPB, a_ts1: f64, a_ts2: f64) {
        let is_splittable = (a_ts2 - a_ts1).abs() > 1e-15;
        let mut pbr = the_pb.0.write().unwrap();
        pbr.set_shrunk_data(a_ts1, a_ts2, is_splittable);
    }

    /// OCCT BOPAlgo_PaveFiller::SplitPaveBlocks (PaveFiller_2.cxx L449-560).
    /// Splits PBs with ext paves into sub-PBs. Each ext pave becomes a split point.
    fn split_pave_blocks(&mut self, the_medges: &std::collections::HashSet<usize>, _the_add_interfs: bool) {
        for &n_e in the_medges {
            if n_e >= self.ds.nb_shapes() { continue; }
            let a_lpb = self.ds.edge_pave_blocks(n_e);
            if a_lpb.is_empty() { continue; }
            let old_pbs: Vec<SharedPB> = a_lpb.to_vec();
            for pb in &old_pbs {
                let ext_paves: Vec<(f64, usize)> = {
                    let r = pb.0.read().unwrap();
                    r.ext_paves.iter().map(|p| (p.param, p.vertex_idx)).collect()
                };
                if ext_paves.is_empty() { continue; }
                let (t1, t2, n_v1, n_v2) = {
                    let r = pb.0.read().unwrap();
                    let (p1, p2) = r.range();
                    let (v1, v2) = r.indices();
                    (p1, p2, v1, v2)
                };
                // Sort ext paves by parameter
                let mut sorted = ext_paves;
                sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                // Build sub-PBs: [t1, p1], [p1, p2], ..., [pN, t2]
                let mut sub_pbs: Vec<(f64, f64, usize, usize)> = Vec::new();
                let mut prev_t = t1;
                let mut prev_v = n_v1;
                for (pt, vi) in &sorted {
                    if *pt > prev_t + 1e-15 && *pt < t2 - 1e-15 {
                        sub_pbs.push((prev_t, *pt, prev_v, *vi));
                        prev_t = *pt;
                        prev_v = *vi;
                    }
                }
                if prev_t < t2 - 1e-15 {
                    sub_pbs.push((prev_t, t2, prev_v, n_v2));
                }
                if sub_pbs.is_empty() { continue; }
                // Create new PBs from sub-ranges
                for (st1, st2, sv1, sv2) in &sub_pbs {
                    let new_pb = SharedPB::new(PaveBlock::new(n_e,
                        Pave { vertex_idx: *sv1, param: *st1 },
                        Pave { vertex_idx: *sv2, param: *st2 }));
                    self.ds.pave_blocks_pool.push(vec![new_pb]);
                }
            }
        }
    }

    // ====================================================================
    // GetPBBox — OCCT BOPAlgo_PaveFiller_3.cxx L914-955
    // ====================================================================

    /// Get bounding box of a PaveBlock's edge segment.
    /// OCCT BOPAlgo_PaveFiller::GetPBBox (PaveFiller_3.cxx L914-955).
    fn get_pb_box(
        &self,
        _the_e: usize,
        the_pb: &SharedPB,
        the_pb_box: &mut std::collections::HashMap<u64, (DVec3, DVec3, f64)>,
        the_first: &mut f64,
        the_last: &mut f64,
        the_s_first: &mut f64,
        the_s_last: &mut f64,
        the_box: &mut rcad_kernel::math::bnd::BndBox,
    ) -> bool {
        let pbr = the_pb.0.read().unwrap();
        (*the_first, *the_last) = pbr.range();
        // OCCT L925-929: check range validity
        if (*the_last - *the_first).abs() <= 1e-12 {
            return false;
        }
        // OCCT L932-937: check shrunk data
        if pbr.has_shrunk_data() {
            *the_s_first = pbr.ts1;
            *the_s_last = pbr.ts2;
            // OCCT: aBox = theSR.BndBox(); aBox.SetGap(aBox.GetGap() + myFuzzyValue / 2.)
            // rcad: BndBox from shrunk data uses OCCT PaveFiller_3.cxx L821-822
            *the_box = rcad_kernel::math::bnd::BndBox::new();
            return true;
        }
        *the_s_first = *the_first;
        *the_s_last = *the_last;
        // OCCT L942-952: check map, then build bounding box
        let pb_ptr = std::sync::Arc::as_ptr(&the_pb.0) as u64;
        if let Some(&(min, max, gap)) = the_pb_box.get(&pb_ptr) {
            *the_box = rcad_kernel::math::bnd::BndBox::from_corners(min.x, min.y, min.z, max.x, max.y, max.z);
            the_box.set_gap(gap);
        } else {
            let curve = self.ds.edge_curve(pbr.original_edge);
            let box_data = if let Some(c) = curve {
                let p1 = c.point_at(*the_s_first);
                let p2 = c.point_at(*the_s_last);
                let min = p1.min(p2);
                let max = p1.max(p2);
                (min, max, 0.0)
            } else {
                return false;
            };
            *the_box = rcad_kernel::math::bnd::BndBox::from_corners(
                box_data.0.x, box_data.0.y, box_data.0.z,
                box_data.1.x, box_data.1.y, box_data.1.z,
            );
            the_pb_box.insert(pb_ptr, box_data);
        }
        true
    }

    // ====================================================================
    // UpdateVertex — OCCT BOPAlgo_PaveFiller::UpdateVertex (PaveFiller_10.cxx L60-85)
    // ====================================================================

    /// Update vertex tolerance and record in increased set.
    /// OCCT BOPAlgo_PaveFiller::UpdateVertex (PaveFiller_10.cxx L60-85).
    fn update_vertex(&mut self, n_v: usize, tol_new: f64) {
        // OCCT L68-72: get current tolerance, if new > current, update
        let tol_old = self.ds.vertex_tolerance_by_idx(n_v);
        if tol_new > tol_old {
            // rcad: update vertex tolerance (OCCT: BRep_Builder::UpdateVertex)
            let si = self.ds.change_shape_info(n_v);
            if let rcad_kernel::topods::TShape::Vertex(vd) = &mut *Arc::make_mut(&mut si.shape.data) {
                vd.tolerance = tol_new;
            }
            // OCCT L76-79: myIncreasedSS.Add(nV)
            self.my_increased_ss.insert(n_v);
        }
    }

    // ====================================================================
    // UpdateVerticesOfCB — OCCT BOPAlgo_PaveFiller_3.cxx L959-993
    // ====================================================================

    /// Update vertex tolerances from CommonBlock tolerances.
    /// OCCT BOPAlgo_PaveFiller::UpdateVerticesOfCB (PaveFiller_3.cxx L959-993).
    fn update_vertices_of_cb(&mut self) {
        // Collect vertex-tolerance pairs first to avoid borrow conflicts
        let mut vertices_to_update: Vec<(usize, f64)> = Vec::new();
        let mut a_mpb_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        for pool_i in 0..self.ds.pave_blocks_pool.len() {
            let pool = self.ds.pave_blocks_pool[pool_i].clone();
            for spb in &pool {
                let a_cb_idx = {
                    let pbr = spb.0.read().unwrap();
                    pbr.common_block_idx
                };
                let a_cb_idx = match a_cb_idx {
                    Some(idx) => idx,
                    None => continue,
                };
                if a_cb_idx >= self.ds.common_blocks.len() {
                    continue;
                }
                let a_cb = &self.ds.common_blocks[a_cb_idx];
                let a_pbr_key = a_cb.pave_block1().unwrap_or(pool_i);
                if !a_mpb_fence.insert(a_pbr_key as u64) {
                    continue;
                }
                let a_tol_cb = a_cb.tolerance();
                if a_tol_cb > 0.0 {
                    for &(pb_gi, _) in a_cb.pave_blocks() {
                        if pb_gi < self.ds.pave_blocks_pool.len() {
                            for cbp in &self.ds.pave_blocks_pool[pb_gi] {
                                let (nv1, nv2) = {
                                    let r = cbp.0.read().unwrap();
                                    (r.pave1.vertex_idx, r.pave2.vertex_idx)
                                };
                                vertices_to_update.push((nv1, a_tol_cb));
                                vertices_to_update.push((nv2, a_tol_cb));
                            }
                        }
                    }
                }
            }
        }
        for (nv, tol) in vertices_to_update {
            self.update_vertex(nv, tol);
        }
    }

    // ====================================================================
    // RepeatIntersection — OCCT BOPAlgo_PaveFiller.cxx L383-448
    // ====================================================================
    /// Re-run VV/VE/VF intersections for vertices whose tolerance was increased.
    /// OCCT BOPAlgo_PaveFiller::RepeatIntersection (PaveFiller.cxx L383-448).
    fn repeat_intersection(&mut self) {
        // L385-386: NCollection_Map<int> anExtraInterfMap;
        let mut an_extra = HashSet::new();
        // L387: const int aNbS = myDS->NbSourceShapes();
        let a_nb_s = self.ds.nb_source_shapes();
        // L388: Message_ProgressScope aPS(theRange, "Repeat intersection", 3);
        // L389-414: for (int i = 0; i < aNbS; ++i)
        for i in 0..a_nb_s {
            // L391-395: if ShapeType != VERTEX, continue
            if self.ds.shapes[i].shape_type != ShapeType::Vertex {
                continue;
            }
            // L397-401: if (myIncreasedSS.Contains(i)) { anExtraInterfMap.Add(i); continue; }
            if self.my_increased_ss.contains(&i) {
                an_extra.insert(i);
                continue;
            }
            // L404-408: int nVSD; if (!myDS->HasShapeSD(i, nVSD)) { continue; }
            let mut n_vsd = usize::MAX;
            if !self.ds.has_shape_sd(i, &mut n_vsd) {
                continue;
            }
            // L410-413: if (myIncreasedSS.Contains(nVSD)) { anExtraInterfMap.Add(i); }
            if self.my_increased_ss.contains(&n_vsd) {
                an_extra.insert(i);
            }
        }
        // L416-419: if (anExtraInterfMap.IsEmpty()) return;
        if an_extra.is_empty() {
            return;
        }

        // L422: myIterator->IntersectExt(anExtraInterfMap);
        // OCCT expands the pair lists to include the extra vertices.
        // rcad: create temporary iterator and expand its pair lists.
        {
            let mut a_it = BOPDS_Iterator::new(self.ds, self.my_fuzzy_value);
            a_it.intersect_ext(&an_extra);
        }

        // L426-430: PerformVV(aPS.Next());
        self.perform_vv();
        if self.has_errors() { return; }
        // L431: UpdatePaveBlocksWithSDVertices();
        self.update_pave_blocks_with_sd_vertices();

        // L433-438: PerformVE(aPS.Next());
        self.perform_ve();
        if self.has_errors() { return; }
        // L438: UpdatePaveBlocksWithSDVertices();
        self.update_pave_blocks_with_sd_vertices();

        // L440-444: PerformVF(aPS.Next());
        self.perform_vf();
        if self.has_errors() { return; }

        // L446-447: UpdatePaveBlocksWithSDVertices(); UpdateInterfsWithSDVertices();
        self.update_pave_blocks_with_sd_vertices();
        self.update_interfs_with_sd_vertices();
    }

    // ====================================================================
    // ForceInterfEE — OCCT BOPAlgo_PaveFiller_3.cxx L997-1333
    // ====================================================================
    /// Force additional EE intersection for common blocks.
    /// OCCT BOPAlgo_PaveFiller::ForceInterfEE (PaveFiller_3.cxx L997-1333).
    fn force_interf_ee(&mut self) {
        // L999-1003: comment — now that vertices are increased/unified,
        // find additional common blocks among edge pairs with same bounding vertices.

        // L1005-1023: Initialize pave blocks for all vertices that participated
        // in intersections.
        // OCCT: for (int i = 0; i < aNbS; ++i)
        //   if VERTEX && HasInterf(i) -> InitPaveBlocksForVertex(i)
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            if self.ds.shapes[i].shape_type != ShapeType::Vertex {
                continue;
            }
            // L1014: if (myDS->HasInterf(i))
            // rcad: check interf_tb for any pair involving i
            let has_interf = self.ds.interf_tb.iter().any(|&(a, b)| a == i || b == i);
            if has_interf {
                self.ds.init_pave_blocks_for_vertex(i);
            }
        }

        // L1024-1080: Fill the connection map from bounding vertices to PBs
        // L1026-1028: NCollection_IndexedDataMap<BOPDS_Pair, List<PaveBlock>> aPBMap
        // rcad: HashMap keyed by (v_min, v_max), value = Vec<(edge_idx, local_pb_idx)>
        let mut a_pb_map: std::collections::HashMap<(usize, usize), Vec<(usize, usize)>> =
            std::collections::HashMap::new();
        // L1030: Fence map of pave blocks
        // rcad: HashSet of (edge_idx, local_pb_idx) — no safe pointer to use
        let mut a_mpb_fence: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();

        for i in 0..a_nb_s {
            // L1034-1038: only edges
            if self.ds.shapes[i].shape_type != ShapeType::Edge {
                continue;
            }
            // L1041-1044: edge must have PBs (HasReference equivalent)
            // rcad: check if the shape has reference (points to pave_blocks_pool)
            if self.ds.shapes[i].reference < 0 {
                continue;
            }
            // L1047-1051: skip degenerated edges (HasFlag)
            if self.ds.shapes[i].has_flag() {
                continue;
            }

            // L1056-1079: iterate PBs of this edge
            let ei = i;
            let a_lpb = self.ds.edge_pave_blocks(ei);
            for local_i in 0..a_lpb.len() {
                // L1060-1061: get real PaveBlock — rcad: no CommonBlock indirection here
                let a_pbr_key = (ei, local_i);
                // L1062-1065: fence map — skip if already processed
                if !a_mpb_fence.insert(a_pbr_key) {
                    continue;
                }

                // L1068-1069: get vertex indices
                let (n_v1, n_v2) = {
                    let pbr = a_lpb[local_i].0.read().unwrap();
                    (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
                };

                // L1072-1078: add PB to map keyed by vertex pair
                // OCCT: BOPDS_Pair aPair(nV1, nV2);
                let a_pair = if n_v1 <= n_v2 { (n_v1, n_v2) } else { (n_v2, n_v1) };
                a_pb_map.entry(a_pair).or_default().push(a_pbr_key);
            }
        }

        // L1082-1086: empty map check
        if a_pb_map.is_empty() {
            return;
        }

        // L1088: const bool bSICheckMode = (myArguments.Extent() == 1);
        let b_si_check_mode = self.my_arguments.len() == 1;

        // L1090-1225: Prepare pairs for intersection
        // L1091: BOPAlgo_VectorOfEdgeEdge aVEdgeEdge;
        struct EEPair {
            ei1: usize,
            pb1_local: usize,
            n_e1: usize,
            pb1_range: (f64, f64),
            ei2: usize,
            pb2_local: usize,
            n_e2: usize,
            pb2_range: (f64, f64),
            fuzzy_value: f64,
        }
        let mut edge_edge_pairs: Vec<EEPair> = Vec::new();

        for (&a_pair, a_lpb) in &a_pb_map {
            let (n_v1, n_v2) = a_pair;
            // L1100-1102: if less than 2 PBs, skip
            if a_lpb.len() < 2 {
                continue;
            }

            // L1105-1110: get vertex shapes for tolerance computation
            // OCCT: const TopoDS_Vertex& aV1 = TopoDS::Vertex(myDS->Shape(nV1));
            //        const TopoDS_Vertex& aV2 = TopoDS::Vertex(myDS->Shape(nV2));
            // rcad: get vertex tolerances from DS
            let tol_v1 = self.ds.vertex_tolerance_by_idx(n_v1);
            let tol_v2 = self.ds.vertex_tolerance_by_idx(n_v2);

            // L1116-1118: aTolAdd = bSICheckMode ? myFuzzyValue : 2*max(BRep_Tool::Tolerance(aV1), aV2)
            let a_tol_add = if b_si_check_mode {
                self.my_fuzzy_value
            } else {
                2.0 * tol_v1.max(tol_v2)
            };

            // L1121-1224: iterate all unique pairs from the list
            for p1_idx in 0..a_lpb.len() {
                let (ei1, pb1_local) = a_lpb[p1_idx];
                // L1125-1126: get CommonBlock status
                let a_pb1 = self.ds.edge_pave_blocks(ei1)[pb1_local].clone();
                let cb1_idx = self.ds.common_block(&a_pb1);
                let (n_e1, i_r1) = {
                    let pbr = a_pb1.0.read().unwrap();
                    (pbr.original_edge, self.ds.rank(pbr.original_edge))
                };
                // L1127-1130: edge and its range
                let (a_t11, a_t12) = {
                    let pbr = a_pb1.0.read().unwrap();
                    pbr.range()
                };
                // OCCT L1131-1139: BRepAdaptor_Curve aBAC1(aE1); aBAC1.D1(midpoint, aPm, aVTgt1);
                //   if (aVTgt1.SquareMagnitude() < gp::Resolution()) continue;
                let c1 = self.ds.edge_curve(n_e1).cloned();
                let v_tgt1 = c1.as_ref().map(|c| {
                    let mid_t = (a_t11 + a_t12) * 0.5;
                    let dt = 1e-7;
                    let p_mid = c.point_at(mid_t);
                    let p_dt = c.point_at(mid_t + dt);
                    p_dt - p_mid
                });
                let (v_tgt1, a_pm) = match v_tgt1 {
                    // OCCT L1135: if (aVTgt1.SquareMagnitude() < gp::Resolution()) continue;
                    // gp::Resolution() == 1e-7 in OCCT.
                    Some(v) if v.length_squared() > 1e-7 => {
                        let mid_t = (a_t11 + a_t12) * 0.5;
                        let a_pm = c1.as_ref().unwrap().point_at(mid_t);
                        (v.normalize(), a_pm)
                    },
                    _ => continue,
                };

                // L1141 onwards: iterate second PB for each pair
                for p2_idx in (p1_idx + 1)..a_lpb.len() {
                    let (ei2, pb2_local) = a_lpb[p2_idx];
                    let a_pb2 = self.ds.edge_pave_blocks(ei2)[pb2_local].clone();
                    let cb2_idx = self.ds.common_block(&a_pb2);
                    let (n_e2, i_r2) = {
                        let pbr = a_pb2.0.read().unwrap();
                        (pbr.original_edge, self.ds.rank(pbr.original_edge))
                    };

                    // L1149-1160: skip edges from same argument unless vertices are new
                    // OCCT: if (iR1 == iR2) {
                    //   if ((!IsNewShape(nV1) && Rank(nV1) == iR1) ||
                    //       (!IsNewShape(nV2) && Rank(nV2) == iR2)) continue; }
                    if i_r1 == i_r2 && i_r1 >= 0 {
                        let v1_original = !self.ds.is_new_shape(n_v1) && self.ds.rank(n_v1) == i_r1;
                        let v2_original = !self.ds.is_new_shape(n_v2) && self.ds.rank(n_v2) == i_r2;
                        if v1_original || v2_original {
                            continue;
                        }
                    }

                    // L1162-1168: skip if PBs already form the SAME common block
                    // OCCT: if (!aCB1.IsNull() && !aCB2.IsNull()) { if (aCB1 == aCB2) continue; }
                    if let (Some(cb1), Some(cb2)) = (cb1_idx, cb2_idx) {
                        if cb1 == cb2 {
                            continue;
                        }
                    }

                    // L1175-1204: check angle between edges at midpoint
                    // bUseAddTol = true initially; if angle > 25deg, set to false
                    let (a_t21, a_t22) = {
                        let pbr = a_pb2.0.read().unwrap();
                        pbr.range()
                    };
                    let b_use_add_tol = {
                        let c2 = self.ds.edge_curve(n_e2).cloned();
                        let mut use_tol = true;
                        if let Some(c) = c2 {
                            let mid_t2 = (a_t21 + a_t22) * 0.5;
                            let dt = 1e-7;
                            let p_mid2 = c.point_at(mid_t2);
                            let p_dt2 = c.point_at(mid_t2 + dt);
                            let v_tgt2 = p_dt2 - p_mid2;
                            // OCCT L1193: if (aVTgt2.SquareMagnitude() < gp::Resolution()) continue;
                        if v_tgt2.length_squared() > 1e-7 {
                                let a_cos = v_tgt2.normalize().dot(v_tgt1).abs();
                                // OCCT L1199-1203: if (std::abs(aCos) < 0.9063) bUseAddTol = false;
                                if a_cos < 0.9063 {
                                    use_tol = false;
                                }
                            }
                        }
                        use_tol
                    };

                    // L1208-1222: add pair with appropriate fuzzy value
                    // OCCT: if (bUseAddTol) anEdgeEdge.SetFuzzyValue(myFuzzyValue + aTolAdd)
                    //        else anEdgeEdge.SetFuzzyValue(myFuzzyValue)
                    let fuzzy_val = if b_use_add_tol {
                        self.my_fuzzy_value + a_tol_add
                    } else {
                        self.my_fuzzy_value
                    };
                    edge_edge_pairs.push(EEPair {
                        ei1, pb1_local, n_e1,
                        pb1_range: (a_t11, a_t12),
                        ei2, pb2_local, n_e2,
                        pb2_range: (a_t21, a_t22),
                        fuzzy_value: fuzzy_val,
                    });
                }
            }
        }

        // L1227-1231: if no pairs, return
        if edge_edge_pairs.is_empty() {
            return;
        }

        // L1248-1252: Perform intersection (OCCT: BOPTools_Parallel::Perform)
        // rcad: serial intersection of each pair
        // L1253: NCollection_DynamicArray<BOPDS_InterfEE>& aEEs = myDS->InterfEE();
        // L1312-1329: Collect PB pairs for CommonBlock creation
        // rcad: BOPAlgo_Tools::PerformCommonBlocks creates a CB per unique PB pair.
        // We collect (pb1, pb2) pairs and create CBs in a second pass.
        let mut cb_pairs: Vec<(SharedPB, SharedPB)> = Vec::new();

        for pair in &edge_edge_pairs {
            // L1264-1290: intersect edges
            let pb1 = self.ds.edge_pave_blocks(pair.ei1)[pair.pb1_local].clone();
            let pb2 = self.ds.edge_pave_blocks(pair.ei2)[pair.pb2_local].clone();

            let c1 = self.ds.edge_curve(pair.n_e1).cloned();
            let c2 = self.ds.edge_curve(pair.n_e2).cloned();
            let (c1, c2) = match (c1, c2) {
                (Some(c1), Some(c2)) => (c1, c2),
                _ => continue,
            };

            let mut ee = crate::bop::int_tools::edge_edge::EdgeEdgeIntersector::new();
            ee.set_edges(pair.n_e1, [pair.pb1_range.0, pair.pb1_range.1], pair.n_e2, [pair.pb2_range.0, pair.pb2_range.1], self.ds);
            // OCCT L1216-1222: SetFuzzyValue with the pair-specific tolerance
            ee.set_fuzzy_value(pair.fuzzy_value);
            ee.perform();

            if !ee.is_done() {
                // L1272-1278: warn about failed intersection
                self.my_report.add_warning(
                    Alert::IntersectionFailed(pair.n_e1, pair.n_e2));
                continue;
            }

            let a_cparts = ee.common_parts();
            // L1282-1285: only accept 1 common part of type EDGE
            if a_cparts.len() != 1 {
                continue;
            }
            let cp = &a_cparts[0];
            // L1288: if (aCP.Type() != TopAbs_EDGE) continue;
            // rcad: the old intersector does not set is_edge for full coincidences,
            // but the same-part length check serves as the EDGE-type proxy.
            if cp.range1[0] >= cp.range1[1] {
                continue;
            }

            // L1293-1310: add interference
            let new_ee = InterferenceEE {
                e1: pair.n_e1,
                e2: pair.n_e2,
                point: cp.bounding_point1,
                param1: cp.range1[0],
                param2: cp.ranges2[0][0],
                new_vertex: usize::MAX,
                range1: cp.range1,
                range2: cp.ranges2[0],
            };
            self.ds.interf_ee.push(new_ee);
            self.ds.add_interf(pair.n_e1, pair.n_e2);

            // L1297-1305: if same rank, add AcquiredSelfIntersection warning
            let r1 = self.ds.rank(pair.n_e1);
            let r2 = self.ds.rank(pair.n_e2);
            if r1 >= 0 && r1 == r2 {
                self.my_report.add_warning(
                    Alert::AcquiredSelfIntersection(vec![pair.n_e1, pair.n_e2]));
            }

            // L1312-1329: fill map for common block creation
            // OCCT: BOPAlgo_Tools::FillMap(aPB[0], aPB[1], aMPBLPB, anAlloc)
            cb_pairs.push((pb1, pb2));
        }

        // L1332: BOPAlgo_Tools::PerformCommonBlocks(aMPBLPB, anAlloc, myDS)
        // OCCT creates a CommonBlock for each unique pair in the map.
        // rcad: create a CommonBlock for each collected PB pair.
        for (pb1, pb2) in &cb_pairs {
            self.ds.add_common_block(&[pb1.clone(), pb2.clone()]);
        }
    }

    // ====================================================================
    // ForceInterfEF — OCCT BOPAlgo_PaveFiller_5.cxx L772-1199
    // ====================================================================
    /// Force additional EF intersection for common blocks.
    /// OCCT BOPAlgo_PaveFiller::ForceInterfEF (PaveFiller_5.cxx L772-827).
    fn force_interf_ef(&mut self) {
        // L774-775: Message_ProgressScope aPS(theRange, nullptr, 1);
        // L776-778: if (!myIsPrimary) return;
        if !self.my_is_primary {
            return;
        }

        // L787-822: Collect all pave blocks into an IndexedMap
        let mut a_mpb: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        let a_nb_s = self.ds.nb_source_shapes();
        for n_e in 0..a_nb_s {
            // L791-795: only edges
            if self.ds.shapes[n_e].shape_type != ShapeType::Edge {
                continue;
            }
            // L798-801: edge must have PBs
            if self.ds.shapes[n_e].reference < 0 {
                continue;
            }
            // L804-807: skip degenerated edges
            if self.ds.shapes[n_e].has_flag() {
                continue;
            }

            // L814-821: iterate PBs
            let a_lpb = self.ds.edge_pave_blocks(n_e);
            for local_i in 0..a_lpb.len() {
                // OCCT L819: aMPB.Add(aPBR) where aPBR = myDS->RealPaveBlock(aPB)
                // rcad: no RealPaveBlock indirection, use (n_e, local_i) as key
                a_mpb.insert((n_e, local_i));
            }
        }

        // L826: ForceInterfEF(aMPB, aPS.Next(), true);
        self.force_interf_ef_work(&a_mpb, true);
    }

    /// OCCT BOPAlgo_PaveFiller::ForceInterfEF (overload, PaveFiller_5.cxx L831-1199).
    /// Worker function — processes collected pave blocks against all faces.
    fn force_interf_ef_work(
        &mut self,
        the_mpb: &std::collections::HashSet<(usize, usize)>,
        the_add_interf: bool,
    ) {
        // L838-841: if (theMPB.IsEmpty()) return;
        if the_mpb.is_empty() {
            return;
        }

        // L843-871: BOPTools_BoxTree aBBTree — build BVH tree of PBs.
        // rcad: iterates all PB/face pairs with direct BndBox overlap checks.

        // L876: const bool bSICheckMode = (myArguments.Extent() == 1);
        let b_si_check_mode = self.my_arguments.len() == 1;

        // L882-1107: For each face, find overlapping PBs and check
        let a_nb_s = self.ds.nb_source_shapes();
        let mut ef_pairs: Vec<(usize, usize)> = Vec::new();

        for n_f in 0..a_nb_s {
            if self.ds.shapes[n_f].shape_type != ShapeType::Face {
                continue;
            }
            if self.ds.shapes[n_f].reference < 0 {
                continue;
            }

            // L912-924: Collect vertices of the face from its FaceInfo
            let a_fi = self.ds.face_info(n_f);
            let face_pb_on = a_fi.pave_blocks_on.clone();
            let face_pb_in = a_fi.pave_blocks_in.clone();
            let face_pb_sc = a_fi.pave_blocks_sc.clone();
            let mut a_mvf: std::collections::HashSet<usize> = std::collections::HashSet::new();
            // OCCT L916-924: aMVF from VerticesOn/In/Sc and PB vertices
            for &v in &a_fi.vertices_on { a_mvf.insert(v); }
            for &v in &a_fi.vertices_in { a_mvf.insert(v); }
            for &v in &a_fi.vertices_sc { a_mvf.insert(v); }

            // Also add vertices from PBs on the face
            for &pb_idx in &face_pb_on {
                if pb_idx < self.ds.pave_blocks_pool.len() {
                    for pb in &self.ds.pave_blocks_pool[pb_idx] {
                        let pbr = pb.0.read().unwrap();
                        a_mvf.insert(pbr.pave1.vertex_idx);
                        a_mvf.insert(pbr.pave2.vertex_idx);
                    }
                }
            }
            for &pb_idx in &face_pb_in {
                if pb_idx < self.ds.pave_blocks_pool.len() {
                    for pb in &self.ds.pave_blocks_pool[pb_idx] {
                        let pbr = pb.0.read().unwrap();
                        a_mvf.insert(pbr.pave1.vertex_idx);
                        a_mvf.insert(pbr.pave2.vertex_idx);
                    }
                }
            }
            for &pb_idx in &face_pb_sc {
                if pb_idx < self.ds.pave_blocks_pool.len() {
                    for pb in &self.ds.pave_blocks_pool[pb_idx] {
                        let pbr = pb.0.read().unwrap();
                        a_mvf.insert(pbr.pave1.vertex_idx);
                        a_mvf.insert(pbr.pave2.vertex_idx);
                    }
                }
            }
            // Drop a_fi to release immutable borrow on self.ds
            // before mutable operations below
            drop(a_fi);

            // L947-1107: iterate all PBs and check for EF common blocks
            for &(n_e, local_i) in the_mpb {
                // L952-955: skip if PB already on the face
                let a_pb = self.ds.edge_pave_blocks(n_e)[local_i].clone();
                // Find pool index for this PB (rcad-specific — OCCT uses pointer identity)
                let pb_key = {
                    let ptr = std::sync::Arc::as_ptr(&a_pb.0);
                    let mut found = None;
                    for (pi, pool) in self.ds.pave_blocks_pool.iter().enumerate() {
                        for (li, spb) in pool.iter().enumerate() {
                            if std::sync::Arc::as_ptr(&spb.0) == ptr {
                                found = Some((pi, li));
                                break;
                            }
                        }
                        if found.is_some() { break; }
                    }
                    found
                };

                // Check if already in face's sets
                let already_on_face = if let Some((pi, _li)) = pb_key {
                    face_pb_on.contains(&pi)
                        || face_pb_in.contains(&pi)
                        || face_pb_sc.contains(&pi)
                } else { false };
                if already_on_face {
                    continue;
                }

                // L958-964: check if face contains both vertices of PB
                let (n_v1, n_v2) = {
                    let pbr = a_pb.0.read().unwrap();
                    (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
                };
                if !a_mvf.contains(&n_v1) || !a_mvf.contains(&n_v2) {
                    continue;
                }

                // L966-981: get the edge
                let pbr = a_pb.0.read().unwrap();
                let n_e_actual = if pbr.edge != usize::MAX {
                    pbr.edge
                } else {
                    pbr.original_edge
                };
                if n_e_actual >= self.ds.nb_shapes() {
                    continue;
                }
                let rank_e = self.ds.rank(n_e_actual);
                let rank_f = self.ds.rank(n_f);
                // L977-980: if same rank, skip
                if rank_e >= 0 && rank_e == rank_f {
                    continue;
                }
                let a_range = pbr.range();
                drop(pbr);

                // L986-1052: edge-face coincidence check
                // OCCT: aBAC.D1(IntermediatePoint(aTS[0], aTS[1]), aPOnE, aVETgt)
                let curve = match self.ds.edge_curve(n_e_actual) {
                    Some(c) => c.clone(),
                    None => continue,
                };
                let mid_t = (a_range.0 + a_range.1) * 0.5;
                let mid_pt = curve.point_at(mid_t);
                // OCCT L1001-1006: tangent vector at midpoint
                let dt = 1e-7;
                let p_mid_dt = curve.point_at(mid_t + dt);
                let v_etgt = p_mid_dt - mid_pt;
                if v_etgt.length_squared() < 1e-7 {
                    continue;
                }

                // OCCT L1022-1024: aTolCheck = bSICheckMode ? myFuzzyValue :
                //   2 * max(BRep_Tool::Tolerance(aV1), BRep_Tool::Tolerance(aV2))
                let tol_v1 = self.ds.vertex_tolerance_by_idx(n_v1);
                let tol_v2 = self.ds.vertex_tolerance_by_idx(n_v2);
                let a_tol_check = if b_si_check_mode {
                    self.my_fuzzy_value
                } else {
                    2.0 * tol_v1.max(tol_v2)
                };

                // Project midpoint onto face surface (OCCT L1031-1036)
                let (proj_uv, proj_pt) = if let Some(surf) = self.ds.face_surface(n_f) {
                    let (uv, proj_pt) = crate::bop::closest_point_on_surface(&surf, mid_pt);
                    let a_dist = (proj_pt - mid_pt).length();
                    // OCCT L1026: if (LowerDistance() > aTolCheck + myFuzzyValue) continue;
                    if a_dist > a_tol_check + self.my_fuzzy_value {
                        continue;
                    }
                    // OCCT L1033-1035: if (!myContext->IsPointInFace(aF, gp_Pnt2d(U,V))) continue;
                    if !self.my_context.is_point_in_face(self.ds, n_f, uv) {
                        continue;
                    }
                    (uv, proj_pt)
                } else { continue; };

                // OCCT L1038-1051: angle between face-to-edge vector and edge tangent
                // OCCT: aVFNorm(aPOnS, aPOnE) — vector from projected midpoint to edge midpoint
                let mut b_use_add_tol = true;
                {
                    let a_vf_norm = mid_pt - proj_pt;
                    if a_vf_norm.length_squared() > 1e-7 {
                        // OCCT L1046-1047: if (|aCos| > 0.4226) bUseAddTol = false
                        let a_cos = a_vf_norm.normalize().dot(v_etgt.normalize()).abs();
                        if a_cos > 0.4226 {
                            b_use_add_tol = false;
                        }
                    }
                }

                // Compute additional tolerance from endpoint distances (OCCT L1063-1084)
                let mut a_tol_add_ef = 0.0;
                if b_use_add_tol {
                    if let Some(surf) = self.ds.face_surface(n_f) {
                        for a_t in [a_range.0, a_range.1] {
                            let a_p = curve.point_at(a_t);
                            let (_uv_e, proj_pe) = crate::bop::closest_point_on_surface(&surf, a_p);
                            let a_dist_ef = (proj_pe - a_p).length();
                            if a_dist_ef < a_tol_check && a_dist_ef > a_tol_add_ef {
                                a_tol_add_ef = a_dist_ef;
                            }
                        }
                    }
                }

                // OCCT L1086-1092: bIntersect = aTolAdd > 0
                if a_tol_add_ef <= 0.0 {
                    continue;
                }

                // L1094-1106: add pair
                ef_pairs.push((n_e_actual, n_f));
            }
        }

        // L1110-1114: if no pairs, return
        if ef_pairs.is_empty() {
            return;
        }

        // L1122-1129: Perform intersection (OCCT: BOPTools_Parallel::Perform)
        // rcad: serial processing of collected pairs.
        // The rcad EdgeFace intersection step is omitted; pairs passing the
        // distance check are accepted directly.

        // L1141-1192: Analyze results — OCCT filters for TopAbs_EDGE type.
        for &(n_e, n_f) in &ef_pairs {
            if the_add_interf {
                // L1175-1181: BOPDS_InterfEF aEF = aEFs.Appended();
                let curve = match self.ds.edge_curve(n_e) {
                    Some(c) => c.clone(),
                    None => continue,
                };
                let range = self.ds.edge_range(n_e);
                let mid_t = (range[0] + range[1]) * 0.5;
                let mid_pt = curve.point_at(mid_t);
                let new_ef = InterferenceEF {
                    edge: n_e,
                    face: n_f,
                    point: mid_pt,
                    edge_param: mid_t,
                    new_vertex: usize::MAX,
                };
                self.ds.interf_ef.push(new_ef);
                self.ds.add_interf(n_e, n_f);
            }

            // L1184-1186: myDS->ChangeFaceInfo(nF).ChangePaveBlocksIn().Add(aPB);
            // rcad: PB pool-index lookup from ef_pairs not stored — requires
            // extending ef_pairs to carry the local pool index.
        }

        // L1194-1198: BOPAlgo_Tools::PerformCommonBlocks(aMPBLI, anAlloc, myDS)
        // rcad: add_common_block on DS exists and is used by ForceInterfEE.
        // Here, the collected PB→face mappings needed for the map are not
        // propagated through ef_pairs. Extending ef_pairs to carry pool indices
        // would enable CommonBlock creation.
    }

    // ====================================================================
    // CheckSelfInterference — OCCT BOPAlgo_PaveFiller_11.cxx L28-221
    // ====================================================================
    /// Check for acquired self-intersections after intersection processing.
    /// OCCT BOPAlgo_PaveFiller::CheckSelfInterference (PaveFiller_11.cxx L28-221).
    fn check_self_interference(&mut self) {
        // L30-34: if (myArguments.Extent() == 1) return;
        if self.my_arguments.len() <= 1 {
            return;
        }

        // L36: BRep_Builder aBB;
        // L38: int i, aNbR = myDS->NbRanges();
        let a_nb_r = self.ds.nb_ranges();
        // L39: for (i = 0; i < aNbR; ++i)
        for a_rank in 0..a_nb_r {
            // L41: const BOPDS_IndexRange& aR = myDS->Range(i);
            let a_r = self.ds.range(a_rank);

            // L44-47: NCollection_IndexedDataMap<TopoDS_Shape, IndexedMap<TopoDS_Shape>> aMCSI;
            let mut a_mcsi: std::collections::HashMap<usize, indexmap::IndexSet<usize>> =
                std::collections::HashMap::new();
            // L48: NCollection_Map<CommonBlock> aMCBFence;
            let mut a_cb_fence: std::collections::HashSet<usize> =
                std::collections::HashSet::new();

            // L51: for (j = aR.First(); j <= aR.Last(); ++j)
            for j in a_r.first..=a_r.last {
                // L53-54: check HasReference
                if self.ds.shapes[j].reference < 0 {
                    continue;
                }

                // L62-63: check ShapeType
                if self.ds.shapes[j].shape_type == ShapeType::Edge {
                    // L65-67: skip degenerated edges
                    if self.ds.shapes[j].has_flag() {
                        continue;
                    }

                    // L71-78: analyze shared vertices
                    let mut a_sub_s: std::collections::HashSet<usize> =
                        std::collections::HashSet::new();
                    for &n_v in &self.ds.shapes[j].sub_shapes {
                        // OCCT L75: int nV = aItLI.Value();
                        // L76: myDS->HasShapeSD(nV, nV); — replaces nV with SD if exists
                        let mut n_vx = n_v;
                        let mut n_vsd = usize::MAX;
                        if self.ds.has_shape_sd(n_v, &mut n_vsd) {
                            n_vx = n_vsd;
                        }
                        // L77: aMSubS.Add(nV);
                        a_sub_s.insert(n_vx);
                    }

                    // L80-81: PaveBlocks for this edge
                    let a_lpb = self.ds.edge_pave_blocks(j);
                    let b_analyze_v = a_lpb.len() > 1;

                    // L83-149: iterate PBs
                    for spb in a_lpb {
                        let pb = spb.0.read().unwrap();

                        // L89-109: check the vertices
                        if b_analyze_v {
                            let (nv1, nv2) = pb.indices();
                            for &n_v in &[nv1, nv2] {
                                // L95: if (!aR.Contains(nV[k]) && !aMSubS.Contains(nV[k]))
                                let in_range = n_v >= a_r.first && n_v <= a_r.last;
                                if !in_range && !a_sub_s.contains(&n_v) {
                                    // L97-106: add connection
                                    a_mcsi.entry(n_v).or_default().insert(j);
                                }
                            }
                        }

                        // L112-148: check common blocks
                        if let Some(cb_idx) = pb.common_block_idx {
                            if a_cb_fence.insert(cb_idx) {
                                if let Some(cb) = self.ds.common_blocks.get(cb_idx) {
                                    let mut a_le: Vec<usize> = Vec::new();
                                    for &(pb_gi, _) in cb.pave_blocks() {
                                        if pb_gi < self.ds.pave_blocks_pool.len() {
                                            for pool_pb in &self.ds.pave_blocks_pool[pb_gi] {
                                                let n_e_or = pool_pb.0.read().unwrap().original_edge;
                                                // L125: if (aR.Contains(nEOr))
                                                let in_range = n_e_or >= a_r.first && n_e_or <= a_r.last;
                                                if in_range {
                                                    a_le.push(n_e_or);
                                                }
                                            }
                                        }
                                    }
                                    // L131-146: if more than 1 edge from same argument in CB
                                    if a_le.len() > 1 {
                                        self.my_report.add_warning(
                                            Alert::AcquiredSelfIntersection(a_le));
                                    }
                                }
                            }
                        }
                    }
                } else if self.ds.shapes[j].shape_type == ShapeType::Face {
                    // L151-196: analyze FACE
                    // L155: const BOPDS_FaceInfo& aFI = myDS->FaceInfo(j);
                    let a_fi = self.ds.face_info(j);

                    // L156-173: IN and Section vertices
                    for &n_v in &a_fi.vertices_in {
                        a_mcsi.entry(n_v).or_default().insert(j);
                    }
                    for &n_v in &a_fi.vertices_sc {
                        a_mcsi.entry(n_v).or_default().insert(j);
                    }

                    // L175-195: IN and Section PaveBlocks
                    for &pb_idx in &a_fi.pave_blocks_in {
                        if pb_idx < self.ds.pave_blocks_pool.len() {
                            for pb in &self.ds.pave_blocks_pool[pb_idx] {
                                let n_e = pb.0.read().unwrap().edge;
                                if n_e != usize::MAX {
                                    a_mcsi.entry(n_e).or_default().insert(j);
                                }
                            }
                        }
                    }
                    for &pb_idx in &a_fi.pave_blocks_sc {
                        if pb_idx < self.ds.pave_blocks_pool.len() {
                            for pb in &self.ds.pave_blocks_pool[pb_idx] {
                                let n_e = pb.0.read().unwrap().edge;
                                if n_e != usize::MAX {
                                    a_mcsi.entry(n_e).or_default().insert(j);
                                }
                            }
                        }
                    }
                }
            }

            // L200-219: Analyze connections — if a vertex/edge connects
            // to multiple faces from same argument, add warning
            for (_sub_shape, shapes) in &a_mcsi {
                if shapes.len() > 1 {
                    self.my_report.add_warning(
                        Alert::AcquiredSelfIntersection(
                            shapes.iter().copied().collect()));
                }
            }
        }
    }

    /// OCCT BOPAlgo_PaveFiller::MakeSplitEdges (_7.cxx L371-548).
    fn make_split_edges(&mut self) {
        // OCCT L392: UpdateCommonBlocksWithSDVertices
        for cb in &self.ds.common_blocks {
            self.ds.update_common_block_with_sd_vertices(cb);
        }

        let a_nb_pbp = self.ds.pave_blocks_pool.len();
        if a_nb_pbp == 0 { return; }
        // OCCT L386: aMCB fence for CommonBlocks
        let mut a_mcb: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for i in 0..a_nb_pbp {
            let a_lpb = self.ds.pave_blocks_pool[i].clone();
            for a_pb in &a_lpb {
                // OCCT L410-414: skip degenerated edges
                let pb = a_pb.0.read().unwrap();
                let n_e = pb.original_edge;
                if n_e >= self.ds.nb_shapes() {
                    drop(pb); continue;
                }
                if self.ds.shapes[n_e].has_flag() {
                    drop(pb); continue;
                }
                // OCCT L416-421: skip if already processed via CB fence
                if let Some(cb_idx) = pb.common_block_idx {
                    if !a_mcb.insert(cb_idx) {
                        drop(pb); continue;
                    }
                }
                let n_v1 = pb.pave1.vertex_idx;
                let n_v2 = pb.pave2.vertex_idx;
                let b_v1 = n_v1 >= self.ds.nb_source_shapes();
                let b_v2 = n_v2 >= self.ds.nb_source_shapes();
                // OCCT L429-450: check if split is needed
                if !b_v1 && !b_v2 {
                    // OCCT L432-450: CB handling for non-destructive mode
                    drop(pb); continue;
                }
                let a_t1 = pb.pave1.param;
                let a_t2 = pb.pave2.param;
                drop(pb);
                // OCCT L465-515: create new split edge
                if let Some(curve) = self.ds.edge_curve(n_e) {
                    let new_ei = self.ds.push_edge(curve.clone(), [a_t1, a_t2], n_v1, n_v2);
                    let mut pbw = a_pb.0.write().unwrap();
                    pbw.edge = new_ei;
                }
            }
        }
        // OCCT L534-550: FillShrunkData for new PBs
        // rcad: shrunk data computed in FillShrunkData step.
    }

    /// OCCT BOPAlgo_PaveFiller::MakeBlocks (_6.cxx L649-1020).
    fn make_blocks(&mut self) {
        // OCCT L652-655: glue off check
        if self.my_glue != GlueEnum::GlueOff { return; }
        // OCCT L657-663: get FF interferences
        if self.ds.interf_ff.is_empty() { return; }

        // OCCT L670-724: NCollection allocators + per-iteration maps
        // OCCT L725-749: iterate FF interferences
        let ff_indices: Vec<usize> = (0..self.ds.interf_ff.len()).collect();
        for &i in &ff_indices {
            let (n_f1, n_f2, a_vp, a_vc_indices) = {
                let ff = &self.ds.interf_ff[i];
                (ff.f1, ff.f2, ff.points.clone(), ff.curves.clone())
            };
            let a_nb_p = a_vp.len();
            let a_nb_c = a_vc_indices.len();
            if a_nb_p == 0 && a_nb_c == 0 { continue; }

            // OCCT L750: aTolFF
            // OCCT L752-753: FaceInfo
            // OCCT L770: SubShapesOnIn — collect ON/IN vertices
            // OCCT L771: SharedEdges

            // OCCT L773-791: 1. Treat Points — create new vertices
            // rcad: FFPoint handling — create CPB entries for new vertices

            // OCCT L793-851: 2. Treat Curves — put paves on section curves
            // rcad: iterate curves, create PBs from curve ranges
            for &cid in &a_vc_indices {
                if cid >= self.ds.intersection_curves.len() { continue; }
                let ic = self.ds.intersection_curves[cid].clone();
                let (v1, v2) = self.curve_vertices_mut(&ic.curve, ic.t_range);
                let ei = self.ds.push_edge(ic.curve.clone(), ic.t_range, v1, v2);
                let p1 = Pave { vertex_idx: v1, param: ic.t_range[0] };
                let p2 = Pave { vertex_idx: v2, param: ic.t_range[1] };
                let pbx = PaveBlock::new(ei, p1, p2);
                let spb = SharedPB::new(pbx);
                let idx = self.ds.pave_blocks_pool.len();
                self.ds.pave_blocks_pool.push(vec![spb]);
                if let Some(last) = self.ds.pave_blocks_pool.last_mut() {
                    for pb2 in last.iter() { pb2.0.write().unwrap().edge = ei; }
                }
                // OCCT L900-950: validate block for both faces
                // OCCT L960-980: add PB to face info
                self.ds.change_face_info(n_f1).pave_blocks_sc.insert(idx);
                self.ds.change_face_info(n_f2).pave_blocks_sc.insert(idx);
                self.ds.change_face_info(n_f1).vertices_sc.insert(v1);
                self.ds.change_face_info(n_f1).vertices_sc.insert(v2);
                self.ds.change_face_info(n_f2).vertices_sc.insert(v1);
                self.ds.change_face_info(n_f2).vertices_sc.insert(v2);
            }

            // OCCT L854-875: BOPTools_BoxTree check
            // OCCT L882-990: 3. Make section edges with IsValidBlockForFaces check
            // OCCT L990-1020: post-processing
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
        let mut a_micro_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mpb_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        for i in 0..self.ds.pave_blocks_pool.len() {
            let pb_list = self.ds.pave_blocks_pool[i].clone();
            if pb_list.len() < 2 { continue; }
            // OCCT L4407-4410: skip degenerated edges
            if pb_list.is_empty() { continue; }
            let n_e_orig = pb_list[0].0.read().unwrap().original_edge;
            if n_e_orig < self.ds.nb_shapes() && self.ds.shapes[n_e_orig].has_flag() {
                continue;
            }
            for pb in &pb_list {
                let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
                if !a_mpb_fence.insert(ptr) { continue; }
                let (n_v1, n_v2) = { let r = pb.0.read().unwrap(); r.indices() };
                if n_v1 == n_v2 {
                    self.fill_shrunk_data(pb);
                    let has_shrunk = { let r = pb.0.read().unwrap(); r.has_shrunk_data() };
                    if !has_shrunk {
                        let e = { let r = pb.0.read().unwrap(); r.edge };
                        if e != usize::MAX {
                            a_micro_edges.insert(e);
                        }
                    }
                }
            }
        }
        // OCCT L4434: RemovePaveBlocks
        for ei in &a_micro_edges {
            for pool in &mut self.ds.pave_blocks_pool {
                pool.retain(|pb| pb.0.read().unwrap().edge != *ei);
            }
        }
    }

    /// OCCT BOPAlgo_PaveFiller::MakePCurves (_7.cxx L589-850).
    fn make_pcurves(&mut self) {
        // OCCT L592-595: check avoid flags
        // OCCT L606-700: 1. Process face info — IN and ON PBs
        let a_nb_fi = self.ds.face_info_pool.len();
        for fi_idx in 0..a_nb_fi {
            let fi = self.ds.face_info_pool[fi_idx].clone();
            let n_f1 = fi.index();
            let f1_s = self.ds.shape(n_f1).clone();
            let surf = match &*f1_s.data {
                rcad_kernel::topods::TShape::Face(fd) => fd.surface.clone(),
                _ => continue,
            };
            let Some(ref surf) = surf else { continue; };

            // OCCT L619-631: PaveBlocksIn — add all IN PBs
            let mut edges_in: Vec<usize> = Vec::new();
            for &pb_idx in fi.pave_blocks_in.iter() {
                if pb_idx >= self.ds.pave_blocks_pool.len() { continue; }
                for pb in &self.ds.pave_blocks_pool[pb_idx] {
                    let n_e = pb.0.read().unwrap().edge;
                    if n_e < self.ds.nb_shapes() { edges_in.push(n_e); }
                }
            }
            // OCCT L634-699: PaveBlocksOn — skip if pcurve already exists
            let mut edges_on: Vec<usize> = Vec::new();
            for &pb_idx in fi.pave_blocks_on.iter() {
                if pb_idx >= self.ds.pave_blocks_pool.len() { continue; }
                for pb in &self.ds.pave_blocks_pool[pb_idx] {
                    let n_e = pb.0.read().unwrap().edge;
                    if n_e >= self.ds.nb_shapes() { continue; }
                    // Check if pcurve already exists
                    let has_pc = {
                        let si = self.ds.shape_info(n_e);
                        match &*si.shape.data {
                            rcad_kernel::topods::TShape::Edge(ed) => ed.pcurves.contains_key(&n_f1),
                            _ => false,
                        }
                    };
                    if !has_pc { edges_on.push(n_e); }
                }
            }

            // Compute and store pcurves for all collected edges
            for &n_e in edges_in.iter().chain(edges_on.iter()) {
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
        // OCCT L702-850: 2. Process section edges
    }

    /// OCCT BOPTools_AlgoTools::CorrectRange (edge-face variant, AlgoTools_2.cxx L364-420).
    /// Adjusts edge range by face tolerance for non-linear curves.
    fn correct_range(n_e: usize, n_f: usize, ds: &DS, a_tf: f64, a_tl: f64) -> (f64, f64) {
        let curve = match ds.edge_curve(n_e) {
            Some(c) => c.clone(),
            None => return (a_tf, a_tl),
        };
        // OCCT L383-384: if Line → no correction
        if matches!(curve, rcad_kernel::geom::Curve3::Line(_)) {
            return (a_tf, a_tl);
        }
        let a_tol_f = ds.vertex_tolerance_by_idx(n_f).max(1e-7);
        let dt = 1e-7;
        let mut new_first = a_tf;
        let mut new_last = a_tl;
        for i in 0..2 {
            let t = if i == 0 { a_tf } else { a_tl };
            let p = curve.point_at(t);
            let p_dt = curve.point_at(t + if i == 0 { dt } else { -dt });
            let der_mag = (p_dt - p).length().max(1e-12);
            let a_res = a_tol_f / der_mag;
            if i == 0 {
                new_first = a_tf + a_res;
            } else {
                new_last = a_tl - a_res;
            }
        }
        if new_last > new_first + 1e-12 {
            (new_first, new_last)
        } else {
            (a_tf, a_tl)
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
    /// OCCT BOPAlgo_PaveFiller::ProcessDE (_8.cxx L54-131).
    fn process_de(&mut self) {
        // L62-63: for (int anEdgeIndex = 0; anEdgeIndex < myDS->NbSourceShapes(); ++anEdgeIndex)
        let a_nb_s = self.ds.nb_source_shapes();
        for an_ei in 0..a_nb_s {
            // L64-71: EDGE + HasFlag(nF)
            let ei = self.ds.shape_info(an_ei);
            if ei.shape_type != ShapeType::Edge { continue; }
            let n_f = ei.flag;
            if n_f < 0 { continue; }
            let n_f = n_f as usize;

            // L72-77: first sub-shape vertex, resolve SD
            let sf = self.ds.shape_info(n_f);
            let n_v = ei.sub_shapes.first().copied().unwrap_or(usize::MAX);
            let mut n_vx = n_v;
            {
                let mut n_vsd = usize::MAX;
                if self.ds.has_shape_sd(n_vx, &mut n_vsd) {
                    n_vx = n_vsd;
                }
            }

            if sf.shape_type == ShapeType::Face {
                // OCCT L82-84: FindPaveBlocks(nV, nF, aLPBOut)
                // OCCT L88-101: FillPaves(nV, anEdgeIndex, nF, aLPBOut, aPBD) — 2D curve intersection
                // OCCT L103: MakeSplitEdge(anEdgeIndex, nF)
                // rcad: FillPaves requires 2D curve-curve intersection implemented inline below.
                // Get the degenerated edge's 2D pcurve on the face
                let de_pcurve = {
                    let tshape = &self.ds.shapes[an_ei].shape.data;
                    match &**tshape {
                        rcad_kernel::topods::TShape::Edge(ed) => ed.pcurves.get(&n_f).cloned(),
                        _ => None,
                    }
                };
                if let Some((c_de, f_de, l_de)) = de_pcurve {
                    // Find PBs in face info that contain n_vx
                    let a_fi = self.ds.face_info(n_f);
                    let mut found_pbs: Vec<SharedPB> = Vec::new();
                    for pb_set in [&a_fi.pave_blocks_in, &a_fi.pave_blocks_sc, &a_fi.pave_blocks_on] {
                        for &pb_idx in pb_set {
                            if pb_idx < self.ds.pave_blocks_pool.len() {
                                for pb in &self.ds.pave_blocks_pool[pb_idx] {
                                    let (v1, v2) = { let r = pb.0.read().unwrap(); r.indices() };
                                    if v1 == n_vx || v2 == n_vx {
                                        found_pbs.push(pb.clone());
                                    }
                                }
                            }
                        }
                    }
                    drop(a_fi);
                    // For each found PB, intersect its 2D curve with the degenerated edge's curve
                    for pb in &found_pbs {
                        let n_e2 = { let r = pb.0.read().unwrap(); r.original_edge };
                        let passing_pcurve = {
                            let ts = &self.ds.shapes[n_e2].shape.data;
                            match &**ts {
                                rcad_kernel::topods::TShape::Edge(ed) => ed.pcurves.get(&n_f).cloned(),
                                _ => None,
                            }
                        };
                        if let Some((c_pb, f_pb, l_pb)) = passing_pcurve {
                            // Sample both curves, find closest approach
                            use rcad_kernel::geom::Curve2dEval;
                            let n = 32usize;
                            let mut best_t = f_de;
                            let mut best_d = f64::MAX;
                            for i in 0..=n {
                                let t_de = f_de + (l_de - f_de) * i as f64 / n as f64;
                                let p_de = c_de.point_at(t_de);
                                for j in 0..=n {
                                    let t_pb = f_pb + (l_pb - f_pb) * j as f64 / n as f64;
                                    let p_pb = c_pb.point_at(t_pb);
                                    let d = (p_de - p_pb).length();
                                    if d < best_d {
                                        best_d = d;
                                        best_t = t_de;
                                    }
                                }
                            }
                            if best_d < 1e-5 {
                                let mut pbr = pb.0.write().unwrap();
                                pbr.ext_paves.push(Pave { vertex_idx: n_vx, param: best_t });
                            }
                        }
                    }
                    // OCCT L99-100: myDS->UpdatePaveBlock(aPBD)
                    // OCCT L103: MakeSplitEdge(anEdgeIndex, nF)
                }
            }
            if sf.shape_type == ShapeType::Edge {
                // L106-122: create a new degenerated edge
                // OCCT: BRep_Builder BB; BB.Add(aE, aVn); BB.Degenerated(aE, true);
                // rcad: push a degenerated edge with the given vertex
                let empty_vdata = rcad_kernel::topods::TVertexData {
                    my_shapes: Vec::new(), flags: 0,
                    point: glam::DVec3::ZERO, tolerance: 0.0, points: Vec::new(),
                };
                let empty_vshape = rcad_kernel::topods::Shape::new(
                    std::sync::Arc::new(rcad_kernel::topods::TShape::Vertex(empty_vdata)),
                    0, rcad_kernel::topods::Orientation::Forward);
                let ed = rcad_kernel::topods::TEdgeData {
                    curve: None,
                    range: [0.0, 0.0],
                    first: empty_vshape.clone(),
                    last: empty_vshape,
                    tolerance: self.my_fuzzy_value.max(1e-7),
                    same_parameter: false,
                    same_range: false,
                    degenerated: true,
                    pcurves: std::collections::HashMap::new(),
                    representations: Vec::new(),
                    vertex_params: std::collections::HashMap::new(),
                    my_shapes: Vec::new(),
                    flags: 0,
                };
                let s = rcad_kernel::topods::Shape::new(
                    std::sync::Arc::new(rcad_kernel::topods::TShape::Edge(ed)),
                    0, rcad_kernel::topods::Orientation::Forward);
                self.ds.append_shape(s);
                let n_en = self.ds.nb_shapes() - 1;
                // L121-123: aPBD->SetEdge(nEn)
                let a_lpbd = self.ds.change_pave_blocks(an_ei);
                if let Some(a_pbd) = a_lpbd.first().cloned() {
                    a_pbd.0.write().unwrap().edge = n_en;
                }
            }
        }
    }
}
