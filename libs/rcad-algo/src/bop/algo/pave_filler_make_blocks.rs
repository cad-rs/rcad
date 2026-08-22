// OCCT BOPAlgo_PaveFiller::MakeBlocks (_6.cxx L649-1137) 1:1 translation.
//
// MakeBlocks builds the section edges from the Face-Face intersection
// curves/points produced by PerformFF:
//   - Treat Points:  promote FF tangent-touch points to DS vertices.
//   - Treat Curves:  put existing ON/IN/stick/EF/bound vertices on the section
//                    curves as paves, splitting each curve's initial pave block.
//   - Make section edges: for every sub-block of a curve, check validity for
//                    both faces, reuse an existing edge if possible, otherwise
//                    create a new section edge + pave block.
//   - Post treatment (PostTreatFF): fuse the created section vertices/edges
//                    with the existing shapes, resolve SD vertices.
//   - UpdateFaceInfo / UpdatePaveBlocks / PutSEInOtherFaces: refresh DS.
//
// OCCT refs:
//   BOPAlgo_PaveFiller_6.cxx  L649-1137 (MakeBlocks)
//   BOPAlgo_PaveFiller_6.cxx  L1141-1161 (MakeSDVerticesFF)
//   BOPAlgo_PaveFiller_6.cxx  L1165-1669 (PostTreatFF)
//   BOPAlgo_PaveFiller_6.cxx  L1673-1946 (UpdateFaceInfo)
//   BOPAlgo_PaveFiller_6.cxx  L1950-2251 (IsExistingVertex, IsExistingPaveBlock x2)
//   BOPAlgo_PaveFiller_6.cxx  L2255-2368 (getBoundPaves, PutBoundPaveOnCurve)
//   BOPAlgo_PaveFiller_6.cxx  L2372-2538 (PutPavesOnCurve, FilterPavesOnCurves)
//   BOPAlgo_PaveFiller_6.cxx  L2542-3068 (ExtendedTolerance, GetEFPnts,
//                                  PutEFPavesOnCurve, PutStickPavesOnCurve,
//                                  GetStickVertices, GetFullShapeMap,
//                                  RemoveUsedVertices, PutPaveOnCurve)
//   BOPAlgo_PaveFiller_6.cxx  L3072-3496 (ProcessExistingPaveBlocks x2,
//                                  UpdateExistingPaveBlocks)
//   BOPAlgo_PaveFiller_6.cxx  L3500-3635 (PutClosingPaveOnCurve, PreparePostTreatFF)
//   BOPAlgo_PaveFiller_6.cxx  L3679-3915 (UpdatePaveBlocks, RemovePaveBlocks)
//   BOPAlgo_PaveFiller_6.cxx  L4072-4304 (CorrectToleranceOfSE, PutSEInOtherFaces,
//                                  RemoveMicroSectionEdges)
//
// Architecture notes (rcad vs OCCT):
//   - BOPDS_Curve::myPaveBlocks  -> face_face::IntersectionCurve::pave_blocks
//     (Vec<SharedPB>, the curve's sub-blocks).  BOPDS_Curve::myPaveBlock1 is
//     simply pave_blocks[0] (InitPaveBlock1 appends an empty PB when empty).
//   - PB-handle-keyed maps (NCollection_Map/IndexedMap<handle<PB>>) -> rcad
//     uses the Arc pointer (u64) as the handle key, mirroring map_pb_cb.
//   - BOPTools_BoxTree (BVH) -> linear scan over the same candidate set
//     (correctness-equal; the BVH is only a performance accelerator).

use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::ds::common_block::CommonBlock;
use crate::bop::ds::pave::{Pave, PaveBlock, SharedPB};
use crate::bop::int_tools::face_face::IntersectionCurve;
use glam::DVec3;
use indexmap::{IndexMap, IndexSet};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::geom::{Curve2dEval, Curve3, SurfaceEval};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, ShapeType, TShape};
use rcad_kernel::CurveEval;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// PB handle (OCCT NCollection handle) = Arc pointer for identity.
fn pb_ptr(pb: &SharedPB) -> u64 {
    Arc::as_ptr(&pb.0) as u64
}

/// BOPTools_AlgoTools::PointOnEdge(aE, aT, aP) — point on the edge's curve.
fn point_on_edge(ds: &crate::bop::ds::DS, n_e: usize, t: f64) -> DVec3 {
    match ds.edge_curve(n_e) {
        Some(c) => c.point_at(t),
        None => DVec3::ZERO,
    }
}

/// Direct children of a Shape (TopoDS_Iterator equivalent).
fn shape_sub_shapes(s: &Shape) -> Vec<Shape> {
    let cp = |sr: &Shape| Shape::from_parts(sr.data.clone(), sr.index, sr.location, sr.orientation);
    match &*s.data {
        TShape::Vertex(_) => vec![],
        TShape::Edge(ed) => vec![cp(&ed.first), cp(&ed.last)],
        TShape::Wire(wd) => wd.edges.iter().map(cp).collect(),
        TShape::Face(fd) => {
            let mut v = vec![cp(&fd.outer_wire)];
            v.extend(fd.inner_wires.iter().map(cp));
            v
        }
        TShape::Shell(sd) => sd.faces.iter().map(cp).collect(),
        TShape::Solid(sd) => sd.shells.iter().map(cp).collect(),
        TShape::CompSolid(cd) => cd.iter().map(cp).collect(),
        TShape::Compound(cd) => cd.iter().map(cp).collect(),
    }
}

/// BOPTools_AlgoTools::MakeVertex / BRepLib::BoundingVertex — the bounding
/// vertex of a list of vertex shapes (average point, max distance tolerance).
fn bounding_vertex(a_lv: &[Shape]) -> (DVec3, f64) {
    if a_lv.is_empty() {
        return (DVec3::ZERO, 0.0);
    }
    if a_lv.len() == 1 {
        let pt = a_lv[0].as_vertex().map(|v| v.point).unwrap_or(DVec3::ZERO);
        // OCCT BRepLib::BoundingVertex uses BRep_Tool::Tolerance (clamped to
        // Precision::Confusion minimum, BRep_Tool.cxx L1314-1333).
        let tol = a_lv[0].as_vertex().map(|v| v.tolerance).unwrap_or(0.0)
            .max(rcad_kernel::precision::CONFUSION);
        return (pt, tol);
    }
    let mut acc = DVec3::ZERO;
    let mut n = 0usize;
    let mut max_tol = 0.0f64;
    for s in a_lv {
        if let Some(v) = s.as_vertex() {
            acc += v.point;
            n += 1;
            max_tol = max_tol.max(v.tolerance.max(rcad_kernel::precision::CONFUSION));
        }
    }
    if n == 0 {
        return (DVec3::ZERO, 0.0);
    }
    let center = acc / n as f64;
    let mut a_max = 0.0f64;
    for s in a_lv {
        if let Some(v) = s.as_vertex() {
            let d = center.distance(v.point);
            if d > a_max { a_max = d; }
        }
    }
    (center, a_max.max(max_tol))
}

/// OCCT BOPDS_CoupleOfPaveBlocks (PaveFiller_6.cxx L705-706): per section
/// shape (vertex or edge) the originating FF interference + curve index and
/// the pave block created from it.
#[derive(Debug, Clone)]
struct CoupleOfPBs {
    index_interf: usize,
    index: usize,
    pb: Option<SharedPB>,
}

impl CoupleOfPBs {
    fn new(index_interf: usize, index: usize) -> Self {
        CoupleOfPBs { index_interf, index, pb: None }
    }
    fn set_pb(&mut self, pb: SharedPB) { self.pb = Some(pb); }
    fn pb1(&self) -> Option<&SharedPB> { self.pb.as_ref() }
}

impl PaveFiller {
    /// OCCT BOPAlgo_PaveFiller::MakeBlocks (_6.cxx L649-1137).
    pub(crate) fn make_blocks(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L652-655: glue off check
        if self.my_glue != crate::bop::algo::GlueEnum::GlueOff { return; }
        // OCCT L657-663: get FF interferences
        let a_ffs = self.ds.interf_ff.clone();
        let mut a_nb_ff = a_ffs.len();
        if a_nb_ff == 0 { return; }

        // OCCT L665-667: locals
        // OCCT L679-681: iterators
        // OCCT L687-717: per-iteration / cross-iteration collections.
        //   aMVTol (UnBind) and aLPB (Remove) need default allocator; the rest
        //   are reused per iteration.  rcad uses plain collections.
        let mut a_lse: Vec<i64> = Vec::new();      // OCCT aLSE (list of int, may hold -1)
        let mut a_lbv: Vec<usize> = Vec::new();    // OCCT aLBV
        let mut a_mv_on_in: HashSet<usize> = HashSet::new();
        let mut a_mv_common: HashSet<usize> = HashSet::new();
        let mut a_mv_stick: HashSet<usize> = HashSet::new();
        let mut a_mv_ef: HashSet<usize> = HashSet::new();
        let mut a_mv_bounds: HashSet<usize> = HashSet::new();
        let mut a_mi: HashSet<usize> = HashSet::new();
        let mut a_mpb_on_in: Vec<SharedPB> = Vec::new(); // IndexedMap<PB>
        let mut a_mpb_common: HashSet<u64> = HashSet::new();
        let mut a_dm_bv: crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>> =
            crate::bop::algo::occt_map::OcctDataMapInt::new(); // j -> LBV (OCCT L705, DataMap bucket order)
        let mut a_mv_tol: crate::bop::algo::occt_map::OcctDataMapInt<usize, f64> =
            crate::bop::algo::occt_map::OcctDataMapInt::new(); // OCCT L707, DataMap bucket order
        // Cross-iteration collections.
        let mut a_mpb_add: HashSet<u64> = HashSet::new();
        let mut a_lpb: Vec<SharedPB> = Vec::new();
        // aMSCPB: IndexedDataMap<TopoDS_Shape, CoupleOfPaveBlocks>
        let mut a_ms_cpb: Vec<(Shape, CoupleOfPBs)> = Vec::new();
        // aMVI: DataMap<TopoDS_Shape, int> — shape identity -> DS index.
        let mut a_mvi: HashMap<(u64, u32), usize> = HashMap::new();
        let mut a_dm_ex_edges: crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<u64>> =
            crate::bop::algo::occt_map::OcctDataMapInt::new();
        let mut a_dm_new_sd: crate::bop::algo::occt_map::OcctDataMapInt<usize, usize> =
            crate::bop::algo::occt_map::OcctDataMapInt::new();
        let mut a_dm_vlv: crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>> =
            crate::bop::algo::occt_map::OcctDataMapInt::new();
        let mut a_micro_pb: Vec<SharedPB> = Vec::new(); // IndexedMap<PB>
        let mut a_verts_on_rejected_pb: Vec<Shape> = Vec::new();
        // OCCT L725: aPBFacesMap is BOPAlgo_DataMapOfPaveBlockListOfInteger —
        // NCollection_DataMap, bucket iteration order.
        let mut a_pb_faces_map: crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<usize>> =
            crate::bop::algo::occt_map::OcctDataMapInt::new();

        // OCCT L720-724: aFFToRecheck — potentially problematic FF pairs to reprocess.
        let mut a_ff_to_recheck: Vec<usize> = Vec::new();
        let a_nb_ff_prev = a_nb_ff;

        let mut i = 0;
        while i < a_nb_ff {
            if the_range.user_break() { return; }
            // OCCT L733: aCurInd = i < aNbFFPrev ? i : aFFToRecheck[i - aNbFFPrev]
            let a_cur_ind = if i < a_nb_ff_prev { i } else { a_ff_to_recheck[i - a_nb_ff_prev] };

            let (n_f1, n_f2, a_vp, a_vc_indices) = {
                let ff = &a_ffs[a_cur_ind];
                (ff.f1, ff.f2, ff.points.clone(), ff.curves.clone())
            };
            let a_nb_p = a_vp.len();
            let a_nb_c = a_vc_indices.len();
            if a_nb_p == 0 && a_nb_c == 0 {
                i += 1;
                continue;
            }
            // OCCT L750: aTolFF = max(Tolerance(aF1), Tolerance(aF2))
            let a_tol_ff = self.ds.face_tolerance(n_f1).max(self.ds.face_tolerance(n_f2));

            // OCCT L755-768: clear per-iteration collections.
            a_mv_on_in.clear();
            a_mv_common.clear();
            a_mpb_on_in.clear();
            a_mpb_common.clear();
            a_dm_bv.clear();
            a_mv_tol.clear();
            a_lse.clear();
            a_lbv.clear();
            a_mv_stick.clear();
            a_mv_ef.clear();
            a_mv_bounds.clear();
            a_mi.clear();

            // OCCT L770-771: SubShapesOnIn + SharedEdges
            self.sub_shapes_on_in(n_f1, n_f2, &mut a_mv_on_in, &mut a_mv_common,
                                  &mut a_mpb_on_in, &mut a_mpb_common);
            self.shared_edges(n_f1, n_f2, &mut a_lse);

            // 1. Treat Points (OCCT L773-791)
            for (j, np) in a_vp.iter().enumerate() {
                let a_p = np.pnt;
                let b_exist = self.is_existing_vertex(a_p, a_tol_ff, &a_mv_on_in);
                if !b_exist {
                    // BOPTools_AlgoTools::MakeNewVertex(aP, aTolFF, aV)
                    let n_v = self.append_vertex(a_p, a_tol_ff);
                    let mut a_cpb = CoupleOfPBs::new(a_cur_ind, j);
                    a_cpb.index = j;
                    let v_shape = self.ds.shape(n_v).clone();
                    a_ms_cpb.push((v_shape, a_cpb));
                }
            }

            // 2. Treat Curves (OCCT L793-851)
            self.get_stick_vertices(n_f1, n_f2, &mut a_mv_stick, &mut a_mv_ef, &mut a_mi);

            for &cid in &a_vc_indices {
                if cid >= self.ds.intersection_curves.len() { continue; }
                // OCCT L800: aNC.InitPaveBlock1() — append an empty PB when empty.
                self.init_pave_block1(cid);
                self.put_paves_on_curve(&a_mv_on_in, &a_mv_common, cid, &a_mi,
                                        &a_mv_ef, &mut a_mv_tol, &mut a_dm_vlv);
            }
            // OCCT L814: FilterPavesOnCurves
            self.filter_paves_on_curves(&a_vc_indices, &mut a_mv_tol);

            for (j, &cid) in a_vc_indices.iter().enumerate() {
                if cid >= self.ds.intersection_curves.len() { continue; }
                // OCCT L821-826: PutStickPavesOnCurve + (aNbC == 1) PutEFPavesOnCurve
                self.put_stick_paves_on_curve(n_f1, n_f2, cid, &a_mi, &a_vc_indices,
                                              &a_mv_stick, &mut a_mv_tol, &mut a_dm_vlv);
                if a_nb_c == 1 {
                    self.put_ef_paves_on_curve(&a_vc_indices, j, &a_mi, &a_mv_ef,
                                               &mut a_mv_tol, &mut a_dm_vlv);
                }
                // OCCT L828-843: PutBoundPaveOnCurve when the curve has bounds.
                let has_bounds = {
                    let ic = &self.ds.intersection_curves[cid];
                    ic.t_range[0].is_finite() && ic.t_range[1].is_finite()
                };
                if has_bounds {
                    a_lbv.clear();
                    self.put_bound_pave_on_curve(n_f1, n_f2, cid, &mut a_lbv);
                    if !a_lbv.is_empty() {
                        a_dm_bv.insert(j, a_lbv.clone());
                        for &v in &a_lbv {
                            a_mv_bounds.insert(v);
                        }
                    }
                }
            }

            // OCCT L847-851: PutClosingPaveOnCurve for each curve.
            for &cid in &a_vc_indices {
                if cid >= self.ds.intersection_curves.len() { continue; }
                self.put_closing_pave_on_curve(cid);
            }

            // OCCT L854-875: build the pave-block tree (linear scan candidate set).
            let mut a_pb_candidates: Vec<usize> = Vec::new(); // indices into a_mpb_on_in
            for (i_pb, a_pb) in a_mpb_on_in.iter().enumerate() {
                let pb = a_pb.0.read().unwrap();
                if !pb.has_edge() { continue; }
                let n_e_orig = pb.original_edge;
                if n_e_orig < self.ds.nb_shapes() && self.ds.shapes[n_e_orig].has_flag() {
                    continue;
                }
                a_pb_candidates.push(i_pb);
            }

            // OCCT L879: isToRecheck
            let mut is_to_recheck = a_nb_c > 0 && i < a_nb_ff_prev;

            // 3. Make section edges (OCCT L882-1066)
            for (j, &cid) in a_vc_indices.iter().enumerate() {
                if cid >= self.ds.intersection_curves.len() { continue; }
                let a_tol_r3d = {
                    let ic = &self.ds.intersection_curves[cid];
                    ic.tolerance.max(ic.tang_tolerance)
                };
                // OCCT L891: aLPB.Clear()
                a_lpb.clear();
                // OCCT L892: aPB1->Update(aLPB, false) — split PaveBlock1 by ext paves.
                {
                    let pb1 = self.ds.intersection_curves[cid].pave_blocks.first().cloned();
                    if let Some(pb1) = pb1 {
                        pb1.0.write().unwrap().update(&mut a_lpb, false);
                    }
                }
                if !a_lpb.is_empty() {
                    is_to_recheck = false;
                }

                let n_lpb = a_lpb.len();
                for k in 0..n_lpb {
                    let a_pb = a_lpb[k].clone();
                    let (n_v1, n_v2) = { let r = a_pb.0.read().unwrap(); r.indices() };
                    let (a_t1, a_t2) = { let r = a_pb.0.read().unwrap(); r.range() };
                    // OCCT L906-909: skip degenerate (zero-range) blocks.
                    if (a_t1 - a_t2).abs() < rcad_kernel::PCONFUSION {
                        continue;
                    }
                    // OCCT L914-918: IsValidBlockForFaces
                    let b_valid_2d = {
                        let ic = self.ds.intersection_curves[cid].clone();
                        self.my_context.is_valid_block_for_faces(
                            a_t1, a_t2, &ic, n_f1, n_f2, &self.ds, a_tol_r3d)
                    };
                    if !b_valid_2d {
                        continue;
                    }
                    // OCCT L920-930: IsExistingPaveBlock (with shared edges).
                    let mut n_e_out = usize::MAX;
                    let mut a_tol_new = -1.0;
                    let b_exist = self.is_existing_pave_block_lse(
                        &a_pb, cid, &a_lse, &mut n_e_out, &mut a_tol_new);
                    if b_exist {
                        // OCCT L925-929: update edge + saved tolerances.
                        self.update_edge_tolerance(n_e_out, a_tol_new);
                        self.update_saved_tolerance(n_e_out, a_tol_new, &mut a_mv_tol);
                        continue;
                    }
                    // OCCT L932-960: FindValidRange for the block.
                    let a_v1_p = self.ds.vertex_point_by_idx(n_v1);
                    let a_v2_p = self.ds.vertex_point_by_idx(n_v2);
                    let a_tol_v1 = self.ds.vertex_tolerance_by_idx(n_v1);
                    let a_tol_v2 = self.ds.vertex_tolerance_by_idx(n_v2);
                    let (a_first, a_last) = {
                        let ic = &self.ds.intersection_curves[cid];
                        let mut sr = crate::bop::algo::pave_filler::ShrunkRange::new(
                            &a_pb, n_v1, n_v2, a_t1, a_t2);
                        if sr.find_valid_range(&ic.curve, a_tol_r3d, a_v1_p,
                                               a_tol_r3d.max(a_tol_v1), a_v2_p,
                                               a_tol_r3d.max(a_tol_v2)) {
                            (sr.shrunk_range().0, sr.shrunk_range().1)
                        } else {
                            (f64::NAN, f64::NAN)
                        }
                    };
                    if a_first.is_nan() {
                        // OCCT L952-958: micro block — keep for post treatment.
                        if !a_mv_bounds.contains(&n_v1) && !a_mv_bounds.contains(&n_v2) {
                            a_micro_pb.push(a_pb.clone());
                            a_mvi.insert((self.ds.shape(n_v1).ptr_id(), self.ds.shape(n_v1).location), n_v1);
                            a_mvi.insert((self.ds.shape(n_v2).ptr_id(), self.ds.shape(n_v2).location), n_v2);
                        }
                        continue;
                    }
                    // OCCT L962-1021: IsExistingPaveBlock (with tree/common).
                    let mut a_pb_out: Option<SharedPB> = None;
                    let mut a_tol_new2 = -1.0;
                    let b_exist2 = self.is_existing_pave_block(
                        &a_pb, cid, a_tol_r3d, &a_mpb_on_in, &a_pb_candidates,
                        &a_mpb_common, &mut a_pb_out, &mut a_tol_new2);
                    if b_exist2 {
                        let a_pb_out = a_pb_out.unwrap();
                        let pb_out_key = pb_ptr(&a_pb_out);
                        // OCCT L966-969
                        let b_in_f1 = self.pb_in_face(n_f1, &a_pb_out);
                        let b_in_f2 = self.pb_in_face(n_f2, &a_pb_out);
                        if !b_in_f1 || !b_in_f2 {
                            // OCCT L972-985: update edge tolerance to touch both faces.
                            let n_e = { let r = a_pb_out.0.read().unwrap(); r.edge };
                            let a_tol_e = self.ds.edge_tolerance(n_e);
                            if a_tol_new2 < a_tol_r3d {
                                a_tol_new2 = a_tol_r3d;
                            }
                            if a_tol_new2 > a_tol_e {
                                self.update_edge_tolerance(n_e, a_tol_new2);
                                self.update_saved_tolerance(n_e, a_tol_new2, &mut a_mv_tol);
                            }
                            // OCCT L988-998: face without pave block.
                            let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                            let p_faces = a_pb_faces_map.bound(pb_out_key);
                            if !p_faces.contains(&n_f) {
                                p_faces.push(n_f);
                            }
                            // OCCT L1000-1012: vertices on rejected PB.
                            let (n_v_out1, n_v_out2) = {
                                let r = a_pb_out.0.read().unwrap();
                                (r.pave1.vertex_idx, r.pave2.vertex_idx)
                            };
                            if n_v1 != n_v_out1 && n_v1 != n_v_out2 && !a_mv_bounds.contains(&n_v1) {
                                a_verts_on_rejected_pb.push(self.ds.shape(n_v1).clone());
                            }
                            if n_v2 != n_v_out1 && n_v2 != n_v_out2 && !a_mv_bounds.contains(&n_v2) {
                                a_verts_on_rejected_pb.push(self.ds.shape(n_v2).clone());
                            }
                            // OCCT L1014-1018: PreparePostTreatFF if newly added.
                            if a_mpb_add.insert(pb_out_key) {
                                self.prepare_post_treat_ff(
                                    a_cur_ind, j, &a_pb_out, &mut a_ms_cpb, &mut a_mvi,
                                    cid);
                            }
                        }
                        continue;
                    }
                    // OCCT L1023-1024: Make Edge
                    let ic_curve = self.ds.intersection_curves[cid].curve.clone();
                    let n_e = {
                        self.append_edge(&ic_curve, a_t1, a_t2, n_v1, n_v2, a_tol_r3d)
                    };
                    // OCCT L1026-1032: BOPTools_AlgoTools::MakePCurve(aES, aF1,
                    // aF2, aIC, PCurveOnS1(), PCurveOnS2(), ctx) — attach the exact
                    // FF 2D intersection curves (aIC.FirstCurve2d/SecondCurve2d)
                    // to the section edge, trimmed to the edge's range [aT1, aT2]
                    // (BOPTools_AlgoTools.cxx L1657-1725).  When the pcurve is null
                    // OCCT falls back to BuildPCurveForEdgeOnFace (projection);
                    // rcad leaves the edge without a pcurve in that case and
                    // make_pcurves part 2 (projection) supplies it later.
                    let a_c2d1 = self.ds.intersection_curves[cid].pcurve1.clone();
                    if let Some(a_c2d) = a_c2d1 {
                        let fk1 = self.ds.face_key(n_f1);
                        self.ds.mutate_shape_data(n_e, |ts| {
                            if let topods::TShape::Edge(ed) = ts {
                                if let Some(k) = fk1 {
                                    ed.pcurves.insert(k, (a_c2d, a_t1, a_t2));
                                }
                            }
                        });
                        self.ds.remap_shape_idx(n_e);
                    }
                    let a_c2d2 = self.ds.intersection_curves[cid].pcurve2.clone();
                    if let Some(a_c2d) = a_c2d2 {
                        let fk2 = self.ds.face_key(n_f2);
                        self.ds.mutate_shape_data(n_e, |ts| {
                            if let topods::TShape::Edge(ed) = ts {
                                if let Some(k) = fk2 {
                                    ed.pcurves.insert(k, (a_c2d, a_t1, a_t2));
                                }
                            }
                        });
                        self.ds.remap_shape_idx(n_e);
                    }
                    // OCCT L1035: aLPBC.Append(aPB)
                    self.ds.intersection_curves[cid].pave_blocks.push(a_pb.clone());
                    // OCCT L1074-1075: the PB is appended to the curve's PB list.
                    // Its edge is NOT set here — the section edge TShape created by
                    // MakeEdge is appended to the DS only in PostTreatFF (L1536-1540,
                    // bOld==false branch: myDS->Append(aSx); aPB1->SetEdge(iE)).
                    // Setting the edge early would make PostTreatFF take the
                    // bOld==true branch and add the PB to aDMExEdges, which then
                    // excludes it from FaceInfo::PaveBlocksSc (UpdateFaceInfo L1757).
                    // rcad keeps an orphan pool entry so FaceInfo::PaveBlocksSc can
                    // reference the PB by pool index; the edge is left unset until
                    // PostTreatFF sets it (mirroring OCCT, where the section edge's
                    // ShapeInfo has no PaveBlocks reference).
                    self.ds.pave_blocks_pool.entry(usize::MAX).or_default().push(a_pb.clone());
                    // OCCT L1038-1045: keep info for post treatment.
                    let e_shape = self.ds.shape(n_e).clone();
                    let mut a_cpb = CoupleOfPBs::new(a_cur_ind, j);
                    a_cpb.set_pb(a_pb.clone());
                    a_ms_cpb.push((e_shape, a_cpb));
                    a_mvi.insert((self.ds.shape(n_v1).ptr_id(), self.ds.shape(n_v1).location), n_v1);
                    a_mvi.insert((self.ds.shape(n_v2).ptr_id(), self.ds.shape(n_v2).location), n_v2);
                    a_mv_tol.remove(n_v1);
                    a_mv_tol.remove(n_v2);
                    // OCCT L1050-1062: ProcessExistingPaveBlocks.
                    self.process_existing_pave_blocks_1(
                        a_cur_ind, j, n_f1, n_f2, n_e, &a_mpb_on_in, &a_pb_candidates,
                        &mut a_ms_cpb, &mut a_mvi, cid, &mut a_pb_faces_map, &mut a_mpb_add);
                }
                // OCCT L1065: aLPBC.RemoveFirst() — remove the initial (split) PB.
                if let Some(pbs) = self.ds.intersection_curves.get_mut(cid) {
                    if !pbs.pave_blocks.is_empty() {
                        pbs.pave_blocks.remove(0);
                    }
                }
            }
            // OCCT L1067-1071: recheck.
            if is_to_recheck {
                a_ff_to_recheck.push(a_cur_ind);
                a_nb_ff += 1;
            }
            // OCCT L1073-1095: restore tolerance of unused vertices — aMVTol is
            // NCollection_DataMap iterated in bucket order (L1120 aItMV).
            for (n_v1, &a_tol) in a_mv_tol.iter() {
                self.set_vertex_tolerance(n_v1, a_tol);
                // reset bounding box
                let pt = self.ds.vertex_point_by_idx(n_v1);
                let vt = self.ds.vertex_tolerance_by_idx(n_v1);
                let mut a_box = BndBox::from_point(pt);
                a_box.set_gap(vt + rcad_kernel::CONFUSION);
                self.ds.change_shape_info(n_v1).bbox = a_box;
                if a_dm_vlv.contains(n_v1) {
                    a_dm_vlv.remove(n_v1);
                }
            }
            // OCCT L1097-1106: ProcessExistingPaveBlocks (with bound vertices).
            self.process_existing_pave_blocks_2(
                a_cur_ind, n_f1, n_f2, &a_mpb_on_in, &a_pb_candidates, &a_dm_bv,
                &mut a_ms_cpb, &mut a_mvi, &mut a_pb_faces_map, &mut a_mpb_add);
            i += 1;
        }

        // OCCT L1110: RemoveMicroSectionEdges
        self.remove_micro_section_edges(&mut a_ms_cpb, &mut a_micro_pb);
        // OCCT L1113: MakeSDVerticesFF
        self.make_sd_vertices_ff(&a_dm_vlv, &mut a_dm_new_sd);
        // OCCT L1114-1120: PostTreatFF
        self.post_treat_ff(&mut a_ms_cpb, &mut a_dm_ex_edges, &mut a_dm_new_sd,
                           &a_micro_pb, &a_verts_on_rejected_pb);
        if self.has_errors() {
            return;
        }
        // OCCT L1126: CorrectToleranceOfSE
        self.correct_tolerance_of_se();
        // OCCT L1129: UpdateFaceInfo
        self.update_face_info(&a_dm_ex_edges, &a_dm_new_sd, &a_pb_faces_map);
        // OCCT L1131: UpdatePaveBlocks
        self.update_pave_blocks(&a_dm_new_sd);
        // OCCT L1136: PutSEInOtherFaces
        self.put_se_in_other_faces();
    }

    // ====================================================================
    // Helper: append a vertex with a box, = BOPTools_AlgoTools::MakeNewVertex
    // followed by BOPDS_ShapeInfo::SetShape + BRepBndLib::Add + SetGap.
    // ====================================================================
    pub(crate) fn append_vertex(&mut self, p: DVec3, tol: f64) -> usize {
        let idx = self.ds.push_vertex(p, tol);
        let mut a_box = BndBox::from_point(p);
        a_box.set_gap(rcad_kernel::CONFUSION);
        self.ds.change_shape_info(idx).bbox = a_box;
        idx
    }

    /// BOPTools_AlgoTools::MakeEdge — create a section edge with the vertex
    /// tolerance raised to theTolR3D + DTolerance and the edge tolerance
    /// raised to theTolR3D.
    fn append_edge(&mut self, curve: &rcad_kernel::geom::Curve3, t1: f64, t2: f64,
                   v1: usize, v2: usize, tol_r3d: f64) -> usize {
        // OCCT L1738: aNeedTol = theTolR3D + DTolerance()
        let a_need_tol = tol_r3d + crate::bop::tools::algo_tools::d_tolerance();
        self.set_vertex_tolerance(v1, self.ds.vertex_tolerance_by_idx(v1).max(a_need_tol));
        self.set_vertex_tolerance(v2, self.ds.vertex_tolerance_by_idx(v2).max(a_need_tol));
        // OCCT L1743: MakeSectEdge — edge with curve + range + vertices.
        let n_e = self.ds.push_edge(curve.clone(), [t1, t2], v1, v2);
        // OCCT L1745: aBB.UpdateEdge(theE, theTolR3D)
        let si = self.ds.change_shape_info(n_e);
        let ts = Arc::make_mut(&mut si.shape.data);
        if let topods::TShape::Edge(ref mut ed) = *ts {
            ed.tolerance = tol_r3d;
        }
        self.ds.remap_shape_idx(n_e);
        n_e
    }

    fn set_vertex_tolerance(&mut self, idx: usize, tol: f64) {
        if idx < self.ds.shapes.len() {
            // In-place (OCCT BRep_Builder UpdateVertex): every reference to the
            // vertex observes the new tolerance without cloning the TShape.
            // Arc::make_mut would clone a shared vertex and split its identity
            // (DS entry vs face-wire references), disconnecting the WireSplitter
            // vertex map. Safe because the DS owns a private input copy.
            self.ds.mutate_shape_data(idx, |ts| {
                if let topods::TShape::Vertex(vd) = ts {
                    vd.tolerance = tol;
                }
            });
            self.ds.remap_shape_idx(idx);
        }
    }

    /// pb_in_face — is the PB's pool entry in the face's On/In sets (OCCT
    /// FaceInfo::PaveBlocksOn/In().Contains(aPB)).
    fn pb_in_face(&self, n_f: usize, pb: &SharedPB) -> bool {
        let key = pb_ptr(pb);
        let fi = self.ds.face_info(n_f);
        fi.pave_blocks_on.contains(&key) || fi.pave_blocks_in.contains(&key)
    }

    // ====================================================================
    // SubShapesOnIn — OCCT BOPDS_DS::SubShapesOnIn (BOPDS_DS.cxx L1066-1143)
    // ====================================================================
    #[allow(clippy::too_many_arguments)]
    fn sub_shapes_on_in(
        &self,
        n_f1: usize,
        n_f2: usize,
        the_mv_on_in: &mut HashSet<usize>,
        the_mv_common: &mut HashSet<usize>,
        the_pb_on_in: &mut Vec<SharedPB>,
        the_common_pave_blocks: &mut HashSet<u64>,
    ) {
        // OCCT L1074-1082: ON/IN PB sets of both faces.
        let a_pb_on1 = self.face_pbs(n_f1, true);
        let a_pb_in1 = self.face_pbs(n_f1, false);
        let a_pb_on2 = self.face_pbs(n_f2, true);
        let a_pb_in2 = self.face_pbs(n_f2, false);
        // OCCT L1084-1102: processMap for all four maps.
        let mut process_map = |the_map: &[SharedPB],
                               the_pb_on_in: &mut Vec<SharedPB>,
                               the_mv_on_in: &mut HashSet<usize>| {
            for a_pb in the_map {
                let key = pb_ptr(a_pb);
                if !the_pb_on_in.iter().any(|p| pb_ptr(p) == key) {
                    the_pb_on_in.push(a_pb.clone());
                }
                let (a_v1, a_v2) = { let r = a_pb.0.read().unwrap(); r.indices() };
                the_mv_on_in.insert(a_v1);
                the_mv_on_in.insert(a_v2);
            }
        };
        process_map(&a_pb_on1, the_pb_on_in, the_mv_on_in);
        process_map(&a_pb_in1, the_pb_on_in, the_mv_on_in);
        process_map(&a_pb_on2, the_pb_on_in, the_mv_on_in);
        process_map(&a_pb_in2, the_pb_on_in, the_mv_on_in);
        // OCCT L1104-1122: find common pave blocks (in Face1 that are also in Face2).
        let mut find_common = |the_map: &[SharedPB],
                               the_common: &mut HashSet<u64>,
                               the_mv_common: &mut HashSet<usize>| {
            for a_pb in the_map {
                let key = pb_ptr(a_pb);
                if a_pb_on2.iter().any(|p| pb_ptr(p) == key)
                    || a_pb_in2.iter().any(|p| pb_ptr(p) == key)
                {
                    the_common.insert(key);
                    let (a_v1, a_v2) = { let r = a_pb.0.read().unwrap(); r.indices() };
                    the_mv_common.insert(a_v1);
                    the_mv_common.insert(a_v2);
                }
            }
        };
        find_common(&a_pb_on1, the_common_pave_blocks, the_mv_common);
        find_common(&a_pb_in1, the_common_pave_blocks, the_mv_common);
        // OCCT L1124-1142: vertices of Face1 that are also in Face2.
        let a_mv_on1: HashSet<usize> = self.ds.face_info(n_f1).vertices_on.iter().copied().collect();
        let a_mv_in1: HashSet<usize> = self.ds.face_info(n_f1).vertices_in.iter().copied().collect();
        let a_mv_on2: HashSet<usize> = self.ds.face_info(n_f2).vertices_on.iter().copied().collect();
        let a_mv_in2: HashSet<usize> = self.ds.face_info(n_f2).vertices_in.iter().copied().collect();
        for &a_v in a_mv_on1.iter().chain(a_mv_in1.iter()) {
            if a_mv_on2.contains(&a_v) || a_mv_in2.contains(&a_v) {
                the_mv_on_in.insert(a_v);
                the_mv_common.insert(a_v);
            }
        }
    }

    /// PBs of a face from one of the On/In sets (PB pointer ids -> SharedPB).
    fn face_pbs(&self, n_f: usize, on: bool) -> Vec<SharedPB> {
        let fi = self.ds.face_info(n_f);
        let set: Vec<u64> = if on {
            fi.pave_blocks_on.iter().copied().collect()
        } else {
            fi.pave_blocks_in.iter().copied().collect()
        };
        let mut out = Vec::new();
        for &pb_ptr in &set {
            if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                out.push(pb);
            }
        }
        out
    }

    // ====================================================================
    // SharedEdges — OCCT BOPDS_DS::SharedEdges (BOPDS_DS.cxx L1147-1208)
    // ====================================================================
    fn shared_edges(&self, n_f1: usize, n_f2: usize, the_edge_list: &mut Vec<i64>) {
        let mut a_first_face_edges: HashSet<usize> = HashSet::new();
        // OCCT L1155-1176: collect edges of the first face.
        for &a_sub in &self.ds.shape_info(n_f1).sub_shapes {
            let a_si = &self.ds.shapes[a_sub];
            if a_si.shape_type != ShapeType::Edge { continue; }
            let a_sub_pbs = self.ds.pave_blocks(a_sub);
            if a_sub_pbs.is_empty() {
                a_first_face_edges.insert(a_sub);
            } else {
                for a_pb in a_sub_pbs {
                    let rpb = self.ds.real_pave_block(a_pb);
                    a_first_face_edges.insert(rpb.0.read().unwrap().edge);
                }
            }
        }
        // OCCT L1179-1207: add edges of the second face contained in the first.
        for &a_sub in &self.ds.shape_info(n_f2).sub_shapes {
            let a_si = &self.ds.shapes[a_sub];
            if a_si.shape_type != ShapeType::Edge { continue; }
            let a_sub_pbs = self.ds.pave_blocks(a_sub);
            if a_sub_pbs.is_empty() {
                if a_first_face_edges.contains(&a_sub) {
                    the_edge_list.push(a_sub as i64);
                }
            } else {
                for a_pb in a_sub_pbs {
                    let rpb = self.ds.real_pave_block(a_pb);
                    let an_edge = rpb.0.read().unwrap().edge;
                    if a_first_face_edges.contains(&an_edge) {
                        the_edge_list.push(an_edge as i64);
                    }
                }
            }
        }
    }

    // ====================================================================
    // IsExistingVertex — OCCT BOPAlgo_PaveFiller::IsExistingVertex
    // (PaveFiller_6.cxx L1950-1984)
    // ====================================================================
    fn is_existing_vertex(&self, a_p: DVec3, the_tol_r3d: f64, a_mv_on_in: &HashSet<usize>) -> bool {
        // OCCT L1960-1964: aTolCheck + aBoxP.
        let a_tol_check = the_tol_r3d + self.my_fuzzy_value;
        let mut a_box_p = BndBox::from_point(a_p);
        a_box_p.enlarge(the_tol_r3d);
        for &n_v in a_mv_on_in {
            if n_v >= self.ds.nb_shapes() { continue; }
            let a_si_v = self.ds.shape_info(n_v);
            let a_box_v = a_si_v.bbox.clone();
            if !a_box_p.is_out_box(&a_box_v) {
                let a_v_pt = self.ds.vertex_point_by_idx(n_v);
                let a_v_tol = self.ds.vertex_tolerance_by_idx(n_v);
                let i_flag = crate::bop::tools::algo_tools::compute_vv_vertex_point(
                    a_v_tol, a_v_pt, a_p, a_tol_check);
                if i_flag == 0 {
                    return true;
                }
            }
        }
        false
    }

    // ====================================================================
    // GetFullShapeMap — OCCT BOPAlgo_PaveFiller::GetFullShapeMap
    // (PaveFiller_6.cxx L2909-2924)
    // ====================================================================
    fn get_full_shape_map(&self, n_f: usize, a_mi: &mut HashSet<usize>) {
        a_mi.insert(n_f);
        for &n_s in &self.ds.shape_info(n_f).sub_shapes {
            a_mi.insert(n_s);
        }
    }

    // ====================================================================
    // GetStickVertices — OCCT BOPAlgo_PaveFiller::GetStickVertices
    // (PaveFiller_6.cxx L2847-2905)
    // ====================================================================
    fn get_stick_vertices(&self, n_f1: usize, n_f2: usize,
                          a_mv_stick: &mut HashSet<usize>,
                          a_mv_ef: &mut HashSet<usize>,
                          a_mi: &mut HashSet<usize>) {
        // OCCT L2861-2865: collect all shapes of the two faces into aMI.
        a_mi.clear();
        self.get_full_shape_map(n_f1, a_mi);
        self.get_full_shape_map(n_f2, a_mi);
        // OCCT L2868-2888: VV, VE, EE, VF interferences.
        // VV
        for a_int in &self.ds.interf_vv {
            if a_int.merged_vertex != usize::MAX {
                let (n_s1, n_s2) = (a_int.v1, a_int.v2);
                if a_mi.contains(&n_s1) && a_mi.contains(&n_s2) {
                    let mut n_v_new = a_int.merged_vertex;
                    self.ds.has_shape_sd(n_v_new, &mut n_v_new);
                    a_mv_stick.insert(n_v_new);
                }
            }
        }
        // VE
        for a_int in &self.ds.interf_ve {
            if a_int.index_new != 0 {
                let (n_s1, n_s2) = (a_int.vertex, a_int.edge);
                if a_mi.contains(&n_s1) && a_mi.contains(&n_s2) {
                    let mut n_v_new = a_int.index_new;
                    self.ds.has_shape_sd(n_v_new, &mut n_v_new);
                    a_mv_stick.insert(n_v_new);
                }
            }
        }
        // EE
        for a_int in &self.ds.interf_ee {
            if a_int.new_vertex != usize::MAX {
                let (n_s1, n_s2) = (a_int.e1, a_int.e2);
                if a_mi.contains(&n_s1) && a_mi.contains(&n_s2) {
                    let mut n_v_new = a_int.new_vertex;
                    self.ds.has_shape_sd(n_v_new, &mut n_v_new);
                    a_mv_stick.insert(n_v_new);
                }
            }
        }
        // VF
        for a_int in &self.ds.interf_vf {
            if let Some(idx_new) = a_int.index_new {
                let (n_s1, n_s2) = (a_int.vertex, a_int.face);
                if a_mi.contains(&n_s1) && a_mi.contains(&n_s2) {
                    let mut n_v_new = idx_new;
                    self.ds.has_shape_sd(n_v_new, &mut n_v_new);
                    a_mv_stick.insert(n_v_new);
                }
            }
        }
        // OCCT L2890-2904: EF interferences.
        for a_int in &self.ds.interf_ef {
            if a_int.new_vertex != usize::MAX {
                let (n_s1, n_s2) = (a_int.edge, a_int.face);
                if a_mi.contains(&n_s1) && a_mi.contains(&n_s2) {
                    let mut n_v_new = a_int.new_vertex;
                    self.ds.has_shape_sd(n_v_new, &mut n_v_new);
                    a_mv_stick.insert(n_v_new);
                    a_mv_ef.insert(n_v_new);
                }
            }
        }
    }

    // ====================================================================
    // RemoveUsedVertices — OCCT BOPAlgo_PaveFiller::RemoveUsedVertices
    // (PaveFiller_6.cxx L2928-2955)
    // ====================================================================
    fn remove_used_vertices(&self, a_vc: &[usize], a_mv: &mut HashSet<usize>) {
        if a_mv.is_empty() { return; }
        for &cid in a_vc {
            if cid >= self.ds.intersection_curves.len() { continue; }
            let a_lpbc = self.ds.intersection_curves[cid].pave_blocks.clone();
            for a_pb in &a_lpbc {
                let (pave1_v, pave2_v) = { let r = a_pb.0.read().unwrap(); r.indices() };
                let ext: Vec<usize> = {
                    let r = a_pb.0.read().unwrap();
                    r.ext_paves.iter().map(|p| p.vertex_idx).collect()
                };
                for v in ext { a_mv.remove(&v); }
                a_mv.remove(&pave1_v);
                a_mv.remove(&pave2_v);
            }
        }
    }

    // ====================================================================
    // InitPaveBlock1 — OCCT BOPDS_Curve::InitPaveBlock1 (BOPDS_Curve.lxx L85-92)
    // ====================================================================
    fn init_pave_block1(&mut self, cid: usize) {
        if cid >= self.ds.intersection_curves.len() { return; }
        if self.ds.intersection_curves[cid].pave_blocks.is_empty() {
            self.ds.intersection_curves[cid].pave_blocks.push(SharedPB::new(PaveBlock::new_curve_block()));
        }
    }

    /// The curve's bounding box (OCCT BOPDS_Curve::Box()).
    /// OCCT PaveFiller_6.cxx L599-606: CheckCurve computes aBox, then
    /// aBox.Enlarge(aBoxExpandValue) and SetBox(aBox) — PutPavesOnCurve reads
    /// this ENLARGED box (L2383: aBoxC = theNC.Box()), not a freshly computed
    /// one.  rcad stores the enlarged box in IntersectionCurve::bbox.
    fn curve_bbox(&self, cid: usize) -> BndBox {
        let ic = &self.ds.intersection_curves[cid];
        if let Some((mn, mx)) = ic.bbox {
            BndBox::from_corners(mn.x, mn.y, mn.z, mx.x, mx.y, mx.z)
        } else {
            let tol = ic.tolerance.max(ic.tang_tolerance);
            match rcad_kernel::curve_bounding_box_range(&ic.curve, ic.t_range[0], ic.t_range[1], tol) {
                Some([mn, mx]) => BndBox::from_corners(mn.x, mn.y, mn.z, mx.x, mx.y, mx.z),
                None => BndBox::new(),
            }
        }
    }

    // ====================================================================
    // PutPavesOnCurve — OCCT BOPAlgo_PaveFiller::PutPavesOnCurve
    // (PaveFiller_6.cxx L2372-2421)
    // ====================================================================
    #[allow(clippy::too_many_arguments)]
    fn put_paves_on_curve(
        &mut self,
        the_mv_on_in: &HashSet<usize>,
        the_mv_common: &HashSet<usize>,
        cid: usize,
        the_mi: &HashSet<usize>,
        the_mv_ef: &HashSet<usize>,
        the_mv_tol: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, f64>,
        the_dm_vlv: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>>,
    ) {
        if cid >= self.ds.intersection_curves.len() { return; }
        let a_box_c = self.curve_bbox(cid);
        let a_tol_r3d = {
            let ic = &self.ds.intersection_curves[cid];
            ic.tolerance.max(ic.tang_tolerance)
        };
        // OCCT L2386-2392: Put EF vertices first.
        let mv_ef_list: Vec<usize> = the_mv_ef.iter().copied().collect();
        for n_v in mv_ef_list {
            self.put_pave_on_curve(n_v, a_tol_r3d, cid, the_mi, the_mv_tol, the_dm_vlv, 2);
        }
        // OCCT L2394-2420: Put all other vertices.
        let mv_on_in_list: Vec<usize> = the_mv_on_in.iter().copied().collect();
        for n_v in mv_on_in_list {
            if the_mv_ef.contains(&n_v) { continue; }
            if !the_mv_common.contains(&n_v) {
                if n_v >= self.ds.nb_shapes() { continue; }
                let a_si_v = self.ds.shape_info(n_v);
                let a_box_v = a_si_v.bbox.clone();
                if a_box_c.is_out_box(&a_box_v) { continue; }
                if !self.ds.is_new_shape(n_v) { continue; }
            }
            self.put_pave_on_curve(n_v, a_tol_r3d, cid, the_mi, the_mv_tol, the_dm_vlv, 1);
        }
    }

    // ====================================================================
    // PutPaveOnCurve — OCCT BOPAlgo_PaveFiller::PutPaveOnCurve
    // (PaveFiller_6.cxx L2959-3068)
    // ====================================================================
    #[allow(clippy::too_many_arguments)]
    fn put_pave_on_curve(
        &mut self,
        n_v: usize,
        a_tol_r3d: f64,
        cid: usize,
        a_mi: &HashSet<usize>,
        a_mv_tol: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, f64>,
        a_dm_vlv: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>>,
        i_check_extend: i32,
    ) {
        if cid >= self.ds.intersection_curves.len() { return; }
        let a_ic = self.ds.intersection_curves[cid].clone();
        let a_pb = a_ic.pave_blocks.first().cloned();
        let Some(a_pb) = a_pb else { return; };
        let a_v_pt = self.ds.vertex_point_by_idx(n_v);
        let a_v_tol = self.ds.vertex_tolerance_by_idx(n_v);
        let mut a_tol_v = a_mv_tol.get(n_v).copied().unwrap_or(a_v_tol);
        // OCCT L2976: IsVertexOnLine
        let mut a_t = 0.0;
        let mut b_is_vertex_on_line = self.my_context.is_vertex_on_line(
            a_v_pt, a_tol_v, &a_ic.curve, a_tol_r3d + self.my_fuzzy_value, &mut a_t, a_ic.t_range);
        // OCCT L2977-2990: extended tolerance check.
        if !b_is_vertex_on_line && i_check_extend != 0 && !self.my_verts_to_avoid_extension.contains(&n_v) {
            let mut an_extra_tol = a_tol_v;
            if self.extended_tolerance(n_v, a_mi, &mut an_extra_tol, i_check_extend) {
                b_is_vertex_on_line = self.my_context.is_vertex_on_line(
                    a_v_pt, an_extra_tol, &a_ic.curve, a_tol_r3d + self.my_fuzzy_value, &mut a_t, a_ic.t_range);
                if b_is_vertex_on_line {
                    let a_p_on_c = a_ic.curve.point_at(a_t);
                    a_tol_v = a_p_on_c.distance(a_v_pt);
                }
            }
        }
        if b_is_vertex_on_line {
            // OCCT L2994-3003: aDTol + aPTol = Resolution(max(aTolR3D, aTolV)).
            let a_dtol = crate::bop::tools::algo_tools::d_tolerance();
            let a_ptol = crate::bop::algo::pave_filler::shrunk_range_resolution(
                &a_ic.curve, a_ic.t_range[0], a_ic.t_range[1], a_tol_r3d.max(a_tol_v));
            // OCCT L3004-3039: ContainsParameter.
            let mut n_v_used = usize::MAX;
            let b_exist = { a_pb.0.read().unwrap().contains_parameter(a_t, a_ptol, &mut n_v_used) };
            if b_exist {
                // OCCT L3008-3032: use existing pave.
                let list = a_dm_vlv.bound(n_v_used);
                if list.is_empty() {
                    list.push(n_v_used);
                    if !a_mv_tol.contains(n_v_used) {
                        let a_tol_used = self.ds.vertex_tolerance_by_idx(n_v_used);
                        a_tol_v = a_tol_used;
                        a_mv_tol.insert(n_v_used, a_tol_v);
                    }
                }
                if !list.contains(&n_v) {
                    list.push(n_v);
                }
                if !a_mv_tol.contains(n_v) {
                    a_tol_v = self.ds.vertex_tolerance_by_idx(n_v);
                    a_mv_tol.insert(n_v, a_tol_v);
                }
            } else {
                // OCCT L3042-3066: add new pave.
                let mut a_pave = Pave::new(n_v, a_t);
                a_pb.0.write().unwrap().append_ext_pave(a_pave);
                let a_p1 = a_ic.curve.point_at(a_t);
                a_tol_v = self.ds.vertex_tolerance_by_idx(n_v);
                let a_p2 = a_v_pt;
                let a_dist = a_p1.distance(a_p2);
                if a_tol_v < a_dist + a_dtol {
                    // OCCT L3054-3055: BRep_Builder().UpdateVertex(aV, aDist + aDTol)
                    self.set_vertex_tolerance(n_v, a_dist + a_dtol);
                    if !a_mv_tol.contains(n_v) {
                        a_mv_tol.insert(n_v, a_tol_v);
                    }
                    // OCCT L3061-3064: rebuild the vertex box.
                    let mut a_box = BndBox::from_point(a_v_pt);
                    a_box.set_gap(rcad_kernel::CONFUSION);
                    self.ds.change_shape_info(n_v).bbox = a_box;
                }
            }
        }
    }

    // ====================================================================
    // ExtendedTolerance — OCCT BOPAlgo_PaveFiller::ExtendedTolerance
    // (PaveFiller_6.cxx L2542-2604)
    // ====================================================================
    fn extended_tolerance(&self, n_v: usize, a_mi: &HashSet<usize>,
                          a_tol_v_ext: &mut f64, a_type: i32) -> bool {
        // OCCT L2548-2551: only new shapes.
        if !self.ds.is_new_shape(n_v) {
            return false;
        }
        let mut k = 0;
        let mut a_nb_int = 2;
        if a_type == 1 {
            a_nb_int = 1;
        } else if a_type == 2 {
            k = 1;
        }
        let a_pv = self.ds.vertex_point_by_idx(n_v);
        while k < a_nb_int {
            if k == 0 {
                // EE interferences
                for a_int in &self.ds.interf_ee {
                    if a_int.new_vertex == n_v {
                        if a_mi.contains(&a_int.e1) && a_mi.contains(&a_int.e2) {
                            let (a_t11, a_t12) = (a_int.range1[0], a_int.range1[1]);
                            let a_p11 = point_on_edge(&self.ds, a_int.e1, a_t11);
                            let a_p12 = point_on_edge(&self.ds, a_int.e1, a_t12);
                            let a_d1 = a_pv.distance(a_p11);
                            let a_d2 = a_pv.distance(a_p12);
                            let a_d = a_d1.max(a_d2);
                            if a_d > *a_tol_v_ext {
                                *a_tol_v_ext = a_d;
                            }
                            return true;
                        }
                    }
                }
            } else {
                // EF interferences
                for a_int in &self.ds.interf_ef {
                    if a_int.new_vertex == n_v {
                        if a_mi.contains(&a_int.edge) && a_mi.contains(&a_int.face) {
                            // OCCT uses the common part Range1; rcad EF stores the
                            // single intersection point, so use its distance.
                            let a_d = a_pv.distance(a_int.point);
                            if a_d > *a_tol_v_ext {
                                *a_tol_v_ext = a_d;
                            }
                            return true;
                        }
                    }
                }
            }
            k += 1;
        }
        false
    }

    // ====================================================================
    // FilterPavesOnCurves — OCCT BOPAlgo_PaveFiller::FilterPavesOnCurves
    // (PaveFiller_6.cxx L2437-2538)
    // ====================================================================
    fn filter_paves_on_curves(&mut self, the_vnc: &[usize], the_mv_tol: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, f64>) {
        // OCCT L2427-2435: struct PaveBlockDist.
        struct PaveBlockDist {
            pb: SharedPB,
            square_dist: f64,
            sin_angle: f64,
            tolerance: f64,
        }
        let an_eps = 1e-18; // gp::Resolution()
        // OCCT L2442: aIDMVertPBs — IndexedDataMap<int, List<PaveBlockDist>>.
        let mut a_idm_vert_pbs: Vec<(usize, Vec<PaveBlockDist>)> = Vec::new();
        for &cid in the_vnc {
            if cid >= self.ds.intersection_curves.len() { continue; }
            let a_ic = self.ds.intersection_curves[cid].clone();
            let a_tol_r3d = a_ic.tolerance.max(a_ic.tang_tolerance);
            let a_pb = match a_ic.pave_blocks.first() {
                Some(pb) => pb.clone(),
                None => continue,
            };
            let a_paves: Vec<Pave> = { let r = a_pb.0.read().unwrap(); r.ext_paves.clone() };
            for a_pave in &a_paves {
                let n_v = a_pave.vertex_idx;
                let a_pv = self.ds.vertex_point_by_idx(n_v);
                let a_par = a_pave.param;
                let a_p_on_c = a_ic.curve.point_at(a_par);
                let a_d1 = a_ic.curve.derivative_at(a_par);
                let a_proj_vec = a_p_on_c - a_pv;
                let a_sq_dist = a_proj_vec.length_squared();
                let a_sq_d1_mod = a_d1.length_squared();
                let mut a_sin = a_proj_vec.cross(a_d1).length_squared();
                if a_sq_dist > an_eps && a_sq_d1_mod > an_eps {
                    a_sin = (a_sin / a_sq_dist / a_sq_d1_mod).sqrt();
                }
                match a_idm_vert_pbs.iter_mut().find(|(v, _)| *v == n_v) {
                    Some(e) => e.1.push(PaveBlockDist {
                        pb: a_pb.clone(), square_dist: a_sq_dist, sin_angle: a_sin,
                        tolerance: a_tol_r3d,
                    }),
                    None => a_idm_vert_pbs.push((n_v, vec![PaveBlockDist {
                        pb: a_pb.clone(), square_dist: a_sq_dist, sin_angle: a_sin,
                        tolerance: a_tol_r3d,
                    }])),
                }
            }
        }
        // OCCT L2485: aSinAngleMin = 0.5.
        const A_SIN_ANGLE_MIN: f64 = 0.5;
        // OCCT L2486-2536: process each vertex.
        let a_nb = a_idm_vert_pbs.len();
        for idx in 0..a_nb {
            let (n_v, a_list) = &a_idm_vert_pbs[idx];
            let n_v = *n_v;
            // OCCT L2491-2502: find the pave with minimal distance.
            let mut a_min_dist = f64::INFINITY;
            for a_pbd in a_list {
                if a_pbd.square_dist < a_min_dist {
                    a_min_dist = a_pbd.square_dist;
                }
            }
            // OCCT L2510-2525: reduce tolerance / remove paves.
            let mut a_max_dist_kept = -1.0;
            let mut is_removed = false;
            for a_pbd in a_list {
                let a_check_dist = 100.0 * (a_pbd.tolerance * a_pbd.tolerance).max(a_min_dist);
                if a_pbd.square_dist > a_check_dist && a_pbd.sin_angle < A_SIN_ANGLE_MIN {
                    a_pbd.pb.0.write().unwrap().remove_ext_pave(n_v);
                    is_removed = true;
                } else if a_pbd.square_dist > a_max_dist_kept {
                    a_max_dist_kept = a_pbd.square_dist;
                }
            }
            // OCCT L2527-2536.
            if is_removed && a_max_dist_kept > 0.0 {
                if let Some(&p_tol) = the_mv_tol.get(n_v) {
                    let a_real_tol = p_tol.max(a_max_dist_kept.sqrt() + rcad_kernel::CONFUSION);
                    self.set_vertex_tolerance(n_v, a_real_tol);
                }
            }
        }
    }

    // ====================================================================
    // PutStickPavesOnCurve — OCCT BOPAlgo_PaveFiller::PutStickPavesOnCurve
    // (PaveFiller_6.cxx L2748-2843)
    // ====================================================================
    #[allow(clippy::too_many_arguments)]
    fn put_stick_paves_on_curve(
        &mut self,
        n_f1: usize,
        n_f2: usize,
        cid: usize,
        a_mi: &HashSet<usize>,
        the_vc: &[usize],
        a_mv_stick: &HashSet<usize>,
        a_mv_tol: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, f64>,
        a_dm_vlv: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>>,
    ) {
        if cid >= self.ds.intersection_curves.len() { return; }
        let a_bnd_nv = self.get_bound_paves(cid);
        // OCCT L2762-2766: both curve ends already have vertices.
        if a_bnd_nv[0] >= 0 && a_bnd_nv[1] >= 0 { return; }
        let mut a_mv: HashSet<usize> = a_mv_stick.clone();
        self.remove_used_vertices(the_vc, &mut a_mv);
        if a_mv.is_empty() { return; }
        let (a_s1, a_s2) = match (self.ds.face_surface(n_f1), self.ds.face_surface(n_f2)) {
            (Some(s1), Some(s2)) => (s1, s2),
            _ => return,
        };
        let a_ic = self.ds.intersection_curves[cid].clone();
        let a_c2d = [a_ic.pcurve1.clone(), a_ic.pcurve2.clone()];
        if let (Some(pc0), Some(pc1)) = (a_c2d[0].clone(), a_c2d[1].clone()) {
            // OCCT L2793-2795: the rich / creasing criteria.
            let a_dt2 = 2e-7;
            let a_dsc_pr = 5e-9;
            let a_tc = a_ic.t_range;
            let a_pc = [a_ic.curve.point_at(a_tc[0]), a_ic.curve.point_at(a_tc[1])];
            let mv_list: Vec<usize> = a_mv.iter().copied().collect();
            for n_v in mv_list {
                let a_pv = self.ds.vertex_point_by_idx(n_v);
                for m in 0..2 {
                    if a_bnd_nv[m] >= 0 { continue; }
                    let a_d2 = a_pc[m].distance_squared(a_pv);
                    if a_d2 > a_dt2 { continue; }
                    // OCCT L2816-2822: surface normals at the curve end point.
                    let mut a_dn = [DVec3::ZERO; 2];
                    for n in 0..2 {
                        let a_s = if n == 0 { &a_s1 } else { &a_s2 };
                        let a_p2d = if n == 0 { pc0.point_at(a_tc[m]) } else { pc1.point_at(a_tc[m]) };
                        a_dn[n] = a_s.normal_at(a_p2d.x, a_p2d.y);
                    }
                    // OCCT L2824-2834: creasing check.
                    let mut a_sc_pr = a_dn[0].dot(a_dn[1]);
                    if a_sc_pr < 0.0 { a_sc_pr = -a_sc_pr; }
                    a_sc_pr = 1.0 - a_sc_pr;
                    if a_sc_pr > a_dsc_pr { continue; }
                    let a_d = a_d2.sqrt();
                    self.put_pave_on_curve(n_v, a_d, cid, a_mi, a_mv_tol, a_dm_vlv, 0);
                }
            }
        }
    }

    // ====================================================================
    // PutEFPavesOnCurve — OCCT BOPAlgo_PaveFiller::PutEFPavesOnCurve
    // (PaveFiller_6.cxx L2692-2744)
    // ====================================================================
    fn put_ef_paves_on_curve(
        &mut self,
        the_vc: &[usize],
        the_index: usize,
        a_mi: &HashSet<usize>,
        a_mv_ef: &HashSet<usize>,
        a_mv_tol: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, f64>,
        a_dm_vlv: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>>,
    ) {
        if a_mv_ef.is_empty() { return; }
        if the_index >= the_vc.len() { return; }
        let cid = the_vc[the_index];
        if cid >= self.ds.intersection_curves.len() { return; }
        let a_ic = self.ds.intersection_curves[cid].clone();
        // OCCT L2707-2711: only Bezier/BSpline curves.
        if !matches!(a_ic.curve, Curve3::BSpline(_) | Curve3::Bezier(_)) { return; }
        let mut a_mv: HashSet<usize> = a_mv_ef.clone();
        self.remove_used_vertices(the_vc, &mut a_mv);
        if a_mv.is_empty() { return; }
        let mv_list: Vec<usize> = a_mv.iter().copied().collect();
        for n_v in mv_list {
            let a_pv = self.ds.vertex_point_by_idx(n_v);
            // OCCT L2726-2741: GeomAPI_ProjectPointOnCurve (ProjPT).
            let a_proj = closest_point_on_curve_range(
                &a_ic.curve, a_pv, a_ic.t_range[0], a_ic.t_range[1], 64);
            let a_dist = a_proj.distance;
            self.put_pave_on_curve(n_v, a_dist, cid, a_mi, a_mv_tol, a_dm_vlv, 0);
        }
    }

    // ====================================================================
    // getBoundPaves — OCCT static getBoundPaves (PaveFiller_6.cxx L2255-2304)
    // ====================================================================
    fn get_bound_paves(&self, cid: usize) -> [i64; 2] {
        let mut the_nv = [-1i64; 2];
        if cid >= self.ds.intersection_curves.len() { return the_nv; }
        let a_ic = self.ds.intersection_curves[cid].clone();
        let a_pb = match a_ic.pave_blocks.first() {
            Some(pb) => pb.clone(),
            None => return the_nv,
        };
        let a_lp: Vec<Pave> = { let r = a_pb.0.read().unwrap(); r.ext_paves.clone() };
        if a_lp.is_empty() { return the_nv; }
        // OCCT L2267-2285: extreme paves.
        let mut a_t_min = f64::INFINITY;
        let mut a_t_max = -f64::INFINITY;
        for a_pv in &a_lp {
            let n_v = a_pv.vertex_idx;
            let a_tv = a_pv.param;
            if a_tv < a_t_min {
                the_nv[0] = n_v as i64;
                a_t_min = a_tv;
            }
            if a_tv > a_t_max {
                the_nv[1] = n_v as i64;
                a_t_max = a_tv;
            }
        }
        // OCCT L2288-2303: compare extreme vertices with the curve ends.
        let a_t = a_ic.t_range;
        let a_p = [a_ic.curve.point_at(a_t[0]), a_ic.curve.point_at(a_t[1])];
        let a_tol = a_ic.tolerance.max(a_ic.tang_tolerance) + rcad_kernel::CONFUSION;
        for j in 0..2 {
            if the_nv[j] < 0 { continue; }
            let n_v = the_nv[j] as usize;
            let a_v_pt = self.ds.vertex_point_by_idx(n_v);
            let a_v_tol = self.ds.vertex_tolerance_by_idx(n_v);
            let i_flag = crate::bop::tools::algo_tools::compute_vv_vertex_point(
                a_v_tol, a_v_pt, a_p[j], a_tol);
            if i_flag != 0 { the_nv[j] = -1; }
        }
        the_nv
    }

    // ====================================================================
    // PutBoundPaveOnCurve — OCCT BOPAlgo_PaveFiller::PutBoundPaveOnCurve
    // (PaveFiller_6.cxx L2308-2368)
    // ====================================================================
    fn put_bound_pave_on_curve(&mut self, n_f1: usize, n_f2: usize, cid: usize,
                               a_lvb: &mut Vec<usize>) {
        if cid >= self.ds.intersection_curves.len() { return; }
        let a_ic = self.ds.intersection_curves[cid].clone();
        let a_t = a_ic.t_range;
        let a_p = [a_ic.curve.point_at(a_t[0]), a_ic.curve.point_at(a_t[1])];
        let a_tol_r3d = a_ic.tolerance.max(a_ic.tang_tolerance);
        let a_pb = match a_ic.pave_blocks.first() {
            Some(pb) => pb.clone(),
            None => return,
        };
        let a_bnd_nv = self.get_bound_paves(cid);
        // OCCT L2323-2328: closed curve with bound paves — nothing to do.
        let a_tol_v_new = rcad_kernel::CONFUSION;
        let is_closed = a_p[1].distance(a_p[0]) <= a_tol_v_new;
        if is_closed && (a_bnd_nv[0] > 0 || a_bnd_nv[1] > 0) { return; }
        for j in 0..2 {
            if a_bnd_nv[j] < 0 {
                // OCCT L2335-2339: closed curve — only process one bound.
                if j == 1 && is_closed { continue; }
                let b_vf = self.my_context.is_valid_point_for_faces(a_p[j], n_f1, n_f2, &self.ds, a_tol_r3d);
                if !b_vf { continue; }
                // OCCT L2345-2347: MakeNewVertex + UpdateVertex (move to curve point).
                let n_vn = self.append_vertex(a_p[j], a_tol_r3d);
                // OCCT L2348: aTolVnew = Tolerance(aVn) — used by the closure test.
                let _ = a_tol_v_new;
                // OCCT L2350-2363: append a new extreme pave.
                let a_pn = Pave::new(n_vn, a_t[j]);
                a_pb.0.write().unwrap().append_ext_pave(a_pn);
                a_lvb.push(n_vn);
            }
        }
    }

    // ====================================================================
    // PutClosingPaveOnCurve — OCCT BOPAlgo_PaveFiller::PutClosingPaveOnCurve
    // (PaveFiller_6.cxx L3500-3605)
    // ====================================================================
    fn put_closing_pave_on_curve(&mut self, cid: usize) {
        if cid >= self.ds.intersection_curves.len() { return; }
        let a_ic = self.ds.intersection_curves[cid].clone();
        // OCCT L3503-3514: check 3d curve and bounds.
        if !(a_ic.t_range[0].is_finite() && a_ic.t_range[1].is_finite()) { return; }
        let a_t = a_ic.t_range;
        let a_p = [a_ic.curve.point_at(a_t[0]), a_ic.curve.point_at(a_t[1])];
        let a_pb = match a_ic.pave_blocks.first() {
            Some(pb) => pb.clone(),
            None => return,
        };
        // OCCT L3521-3547: find the pave put at one of the ends.
        let mut n_v = usize::MAX;
        let mut a_t_op = 0.0;
        let mut a_p_op = DVec3::ZERO;
        let a_lp: Vec<Pave> = { let r = a_pb.0.read().unwrap(); r.ext_paves.clone() };
        for a_pv in &a_lp {
            if n_v != usize::MAX { break; }
            let a_tc = a_pv.param;
            for j in 0..2 {
                if (a_tc - a_t[j]).abs() < rcad_kernel::PCONFUSION {
                    n_v = a_pv.vertex_idx;
                    a_t_op = if j == 0 { a_t[1] } else { a_t[0] };
                    a_p_op = if j == 0 { a_p[1] } else { a_p[0] };
                    break;
                }
            }
        }
        if n_v == usize::MAX { return; }
        // OCCT L3557-3569: check if the curve is closed.
        let a_v_pt = self.ds.vertex_point_by_idx(n_v);
        let mut a_tol_v = self.ds.vertex_tolerance_by_idx(n_v);
        let a_tol_p = a_ic.tolerance.max(a_ic.tang_tolerance) + rcad_kernel::CONFUSION;
        let a_dist_vp = a_v_pt.distance(a_p_op);
        if a_dist_vp > a_tol_v + a_tol_p { return; }
        // OCCT L3572-3587: valid range check.
        let a_new_tol_v = a_tol_v.max(a_dist_vp + crate::bop::tools::algo_tools::d_tolerance());
        let mut sr = crate::bop::algo::pave_filler::ShrunkRange::new(
            &a_pb, n_v, n_v, a_t[0], a_t[1]);
        if !sr.find_valid_range(&a_ic.curve, a_ic.tolerance, a_p[0], a_new_tol_v, a_p[1], a_new_tol_v) {
            return;
        }
        // OCCT L3589-3598: UpdateVertex — may replace the vertex with a new
        // one (old vertex path, PaveFiller_10.cxx L121-150); the closing pave
        // must reference the (possibly new) vertex so its index differs from
        // the t=0 endpoint pave and passes the ext-pave fence.
        if a_new_tol_v > a_tol_v {
            let n_vn = self.update_vertex(n_v, a_new_tol_v);
            if n_vn != n_v {
                n_v = n_vn;
            }
            a_tol_v = self.ds.vertex_tolerance_by_idx(n_v);
        }
        // OCCT L3601-3604: add the closing pave.  OCCT appends it directly
        // (aLP.Append(aNewPave)) — NOT through AppendExtPave, so the
        // vertex-index fence (myMFence) does not reject a closing pave that
        // reuses the bound vertex index (closed curve: same vertex at both
        // ends).  rcad append_ext_pave1 is the fence-less equivalent.
        let a_new_pave = Pave::new(n_v, a_t_op);
        a_pb.0.write().unwrap().append_ext_pave1(a_new_pave);
    }

    // ====================================================================
    // UpdateEdgeTolerance — OCCT BOPAlgo_PaveFiller::UpdateEdgeTolerance
    // (PaveFiller_10.cxx L63-101)
    // ====================================================================
    pub(crate) fn update_edge_tolerance(&mut self, n_e: usize, the_tol: f64) {
        // OCCT L69-85: safe input mode — avoid modifying the input shapes.
        if self.my_non_destructive && !self.ds.is_new_shape(n_e) {
            return;
        }
        // OCCT L88-92: update edge + rebuild its box. In-place edit of the
        // shared TShape (OCCT BRep_Builder::UpdateEdge semantics) — an edge
        // referenced with different Locations (prism sweep copies) must keep
        // its TShape identity, or the result gains duplicate edges.
        self.ds.mutate_shape_data(n_e, |ts| {
            if let topods::TShape::Edge(ref mut ed) = *ts {
                if the_tol > ed.tolerance {
                    ed.tolerance = the_tol;
                }
            }
        });
        self.ds.remap_shape_idx(n_e);
        self.rebuild_edge_box(n_e);
        // OCCT L94-100: update sub-vertices.
        let sub_shapes: Vec<usize> = self.ds.shape_info(n_e).sub_shapes.clone();
        for n_v in sub_shapes {
            self.update_vertex(n_v, the_tol);
        }
    }

    /// Rebuild an edge's bounding box (OCCT BRepBndLib::Add(aE, aBoxE)).
    pub(crate) fn rebuild_edge_box(&mut self, n_e: usize) {
        let curve = match self.ds.edge_curve(n_e) { Some(c) => c.clone(), None => return };
        let range = self.ds.edge_range(n_e);
        let tol = self.ds.edge_tolerance(n_e);
        let mut b = match rcad_kernel::curve_bounding_box_range(&curve, range[0], range[1], tol) {
            Some([mn, mx]) => BndBox::from_corners(mn.x, mn.y, mn.z, mx.x, mx.y, mx.z),
            None => BndBox::new(),
        };
        b.set_gap(b.get_gap() + rcad_kernel::CONFUSION);
        self.ds.change_shape_info(n_e).bbox = b;
    }

    /// Edge bounding box for a section edge.
    fn edge_bbox(&self, n_e: usize) -> BndBox {
        self.ds.shape_info(n_e).bbox.clone()
    }

    // ====================================================================
    // UpdateSavedTolerance — OCCT static UpdateSavedTolerance (PaveFiller_6.cxx L629-645)
    // ====================================================================
    fn update_saved_tolerance(&self, n_e: usize, the_tol_new: f64,
                              the_mv_tol: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, f64>) {
        let sub_shapes: Vec<usize> = self.ds.shape_info(n_e).sub_shapes.clone();
        for n_v in sub_shapes {
            if let Some(p_tol_saved) = the_mv_tol.get_mut(n_v) {
                if *p_tol_saved < the_tol_new {
                    *p_tol_saved = the_tol_new;
                }
            }
        }
    }

    // ====================================================================
    // IsExistingPaveBlock (with SharedEdges) — OCCT BOPAlgo_PaveFiller::
    // IsExistingPaveBlock (PaveFiller_6.cxx L1988-2043)
    // ====================================================================
    fn is_existing_pave_block_lse(&mut self, the_pb: &SharedPB, cid: usize,
                                  the_lse: &[i64], the_n_e_out: &mut usize,
                                  the_tol_new: &mut f64) -> bool {
        if the_lse.is_empty() {
            return false;
        }
        let (a_t1, a_t2) = { let r = the_pb.0.read().unwrap(); r.range() };
        let (n_v1, n_v2) = { let r = the_pb.0.read().unwrap(); r.indices() };
        let a_tol_v1 = self.ds.vertex_tolerance_by_idx(n_v1);
        let a_tol_v2 = self.ds.vertex_tolerance_by_idx(n_v2);
        let a_tol = a_tol_v1.max(a_tol_v2);
        let a_tm = crate::bop::int_tools::face_make_curve::intermediate_point(a_t1, a_t2);
        let a_pm = {
            let ic = &self.ds.intersection_curves[cid];
            ic.curve.point_at(a_tm)
        };
        let mut a_box_pm = BndBox::from_point(a_pm);
        a_box_pm.enlarge(a_tol);
        for &n_e in the_lse {
            if n_e < 0 { continue; }
            let n_e = n_e as usize;
            let a_si_e = self.ds.shape_info(n_e);
            let a_box_e = a_si_e.bbox.clone();
            if !a_box_e.is_out_box(&a_box_pm) {
                // OCCT L2030-2039: ComputePE.
                let a_tol_e = self.ds.edge_tolerance(n_e);
                let a_tol_check = a_tol_e.max(a_tol) + self.my_fuzzy_value;
                let mut a_tx = 0.0;
                let mut a_dist = 0.0;
                let i_flag = self.my_context.compute_pe(a_pm, a_tol_check, n_e, &self.ds, &mut a_tx, &mut a_dist);
                if i_flag == 0 {
                    *the_n_e_out = n_e;
                    *the_tol_new = a_dist;
                    return true;
                }
            }
        }
        false
    }

    // ====================================================================
    // IsExistingPaveBlock (BoxTree variant) — OCCT BOPAlgo_PaveFiller::
    // IsExistingPaveBlock (PaveFiller_6.cxx L2047-2251)
    // ====================================================================
    fn is_existing_pave_block(&mut self, the_pb: &SharedPB, cid: usize, the_tol_r3d: f64,
                              the_mpb_on_in: &[SharedPB], pb_candidates: &[usize],
                              the_mpb_common: &HashSet<u64>,
                              a_pb_out: &mut Option<SharedPB>, the_tol_new: &mut f64) -> bool {
        let a_ic = self.ds.intersection_curves[cid].clone();
        let (a_t1, a_t2) = { let r = the_pb.0.read().unwrap(); r.range() };
        let (n_v11, n_v12) = { let r = the_pb.0.read().unwrap(); r.indices() };
        // OCCT L2066-2072: first point box.
        let a_p1 = a_ic.curve.point_at(a_t1);
        let mut a_box_p1 = BndBox::from_point(a_p1);
        let a_tol_v11 = self.ds.vertex_tolerance_by_idx(n_v11);
        a_box_p1.enlarge(a_tol_v11);
        // OCCT L2074-2080: BoxTree selector for the first point.
        let mut a_candidates: Vec<usize> = Vec::new();
        for &i_pb in pb_candidates {
            let a_pb2 = &the_mpb_on_in[i_pb];
            let n_e2 = { let r = a_pb2.0.read().unwrap(); r.edge };
            let a_box_sp = self.ds.shape_info(n_e2).bbox.clone();
            if !a_box_sp.is_out_box(&a_box_p1) {
                a_candidates.push(i_pb);
            }
        }
        if a_candidates.is_empty() {
            return false;
        }
        // OCCT L2082-2094: intermediate point.
        let a_tm = crate::bop::int_tools::face_make_curve::intermediate_point(a_t1, a_t2);
        let a_pm = a_ic.curve.point_at(a_tm);
        let a_vtgt1 = a_ic.curve.derivative_at(a_tm);
        let mut a_box_pm = BndBox::from_point(a_pm);
        let _ = &mut a_box_pm;
        let is_vtgt1_valid = a_vtgt1.length_squared() > 1e-18;
        let a_vtgt1 = if is_vtgt1_valid { a_vtgt1.normalize() } else { a_vtgt1 };
        // OCCT L2097-2102: last point box.
        let a_p2 = a_ic.curve.point_at(a_t2);
        let mut a_box_p2 = BndBox::from_point(a_p2);
        let a_tol_v12 = self.ds.vertex_tolerance_by_idx(n_v12);
        a_box_p2.enlarge(a_tol_v12);
        let a_tol_v1 = a_tol_v11.max(a_tol_v12) + self.my_fuzzy_value;
        let a_tol_check = the_tol_r3d + self.my_fuzzy_value;
        // OCCT L2110-2112: thin-face limit values.
        let mut a_max_tol_add = 0.001f64;
        let a_coeff_tol_add = 10.0;
        a_max_tol_add = a_max_tol_add.min(a_coeff_tol_add * a_tol_check);
        let mut b_found = false;
        *the_tol_new = f64::MAX;
        for &i_pb in &a_candidates {
            let a_pb = &the_mpb_on_in[i_pb];
            let (n_v21, n_v22) = { let r = a_pb.0.read().unwrap(); r.indices() };
            let a_tol_v21 = self.ds.vertex_tolerance_by_idx(n_v21);
            let a_tol_v22 = self.ds.vertex_tolerance_by_idx(n_v22);
            let a_tol_v2 = a_tol_v21.max(a_tol_v22) + self.my_fuzzy_value;
            let n_e_sp = { let r = a_pb.0.read().unwrap(); r.edge };
            let a_si_sp = self.ds.shape_info(n_e_sp);
            let a_box_sp = a_si_sp.bbox.clone();
            let i_flag1 = if n_v11 == n_v21 || n_v11 == n_v22 { 2 } else { 1 };
            let i_flag2 = if n_v12 == n_v21 || n_v12 == n_v22 {
                2
            } else if !a_box_sp.is_out_box(&a_box_p2) {
                1
            } else {
                0
            };
            if i_flag2 == 0 { continue; }
            let mut a_dist = 0.0;
            let mut a_coeff = 1.0;
            let mut a_dist_m1m2 = 0.0;
            let mut a_pe_status = 1;
            let mut a_real_tol = a_tol_check;
            if self.ds.is_common_block(a_pb) {
                a_real_tol = a_real_tol.max(a_tol_v1.max(a_tol_v2));
                if the_mpb_common.contains(&pb_ptr(a_pb)) {
                    a_real_tol *= 2.0;
                }
            } else if i_flag1 == 2 && i_flag2 == 2 {
                // OCCT L2163-2164: skip if one edge is closed but the other is not.
                let b_skip_processing =
                    (n_v11 == n_v12 && n_v21 != n_v22) || (n_v11 != n_v12 && n_v21 == n_v22);
                if !b_skip_processing && is_vtgt1_valid {
                    // OCCT L2171-2173: only if not both lines.
                    let edge_curve = match self.ds.edge_curve(n_e_sp) { Some(c) => c.clone(), None => continue };
                    if !matches!(a_ic.curve, Curve3::Line(_)) || !matches!(edge_curve, Curve3::Line(_)) {
                        let a_tol_add = 2.0 * a_max_tol_add.min(a_real_tol.max(a_tol_v1.max(a_tol_v2)));
                        let mut a_tldp = 0.0;
                        let mut a_dist_pe = 0.0;
                        a_pe_status = self.my_context.compute_pe(
                            a_pm, a_tol_add, n_e_sp, &self.ds, &mut a_tldp, &mut a_dist_pe);
                        if a_pe_status == 0 {
                            let a_vtgt2 = edge_curve.derivative_at(a_tldp);
                            if a_vtgt2.length_squared() > 1e-18 {
                                let a_cos = a_vtgt1.dot(a_vtgt2.normalize());
                                if a_cos.abs() >= 0.9063 {
                                    a_real_tol = a_tol_add;
                                    a_coeff = 2.0;
                                }
                            }
                        }
                    }
                }
            }
            let mut a_box_tmp = a_box_pm.clone();
            a_box_tmp.enlarge(a_real_tol);
            let mut a_dist_to_sp = 0.0;
            if a_box_sp.is_out_box(&a_box_tmp) || a_pe_status < 0 {
                continue;
            } else if a_pe_status == 0 {
                a_dist_to_sp = a_dist_m1m2;
            } else if a_pe_status == 1 {
                let mut a_tx = 0.0;
                a_pe_status = self.my_context.compute_pe(
                    a_pm, a_real_tol, n_e_sp, &self.ds, &mut a_tx, &mut a_dist_to_sp);
                if a_pe_status < 0 {
                    continue;
                }
            }
            let mut i_flag1 = i_flag1;
            if i_flag1 == 1 {
                let mut a_tx = 0.0;
                let f = self.my_context.compute_pe(
                    a_p1, a_real_tol, n_e_sp, &self.ds, &mut a_tx, &mut a_dist);
                i_flag1 = if f == 0 { 1 } else { 0 };
                if i_flag1 != 0 && a_dist_to_sp < a_dist {
                    a_dist_to_sp = a_dist;
                }
            }
            let mut i_flag2 = i_flag2;
            if i_flag2 == 1 {
                let mut a_tx = 0.0;
                let f = self.my_context.compute_pe(
                    a_p2, a_real_tol, n_e_sp, &self.ds, &mut a_tx, &mut a_dist);
                i_flag2 = if f == 0 { 1 } else { 0 };
                if i_flag2 != 0 && a_dist_to_sp < a_dist {
                    a_dist_to_sp = a_dist;
                }
            }
            if i_flag1 != 0 && i_flag2 != 0 {
                if a_dist_to_sp < *the_tol_new {
                    *a_pb_out = Some(a_pb.clone());
                    *the_tol_new = a_coeff * a_dist_to_sp;
                    b_found = true;
                }
            }
        }
        b_found
    }

    // ====================================================================
    // PreparePostTreatFF — OCCT BOPAlgo_PaveFiller::PreparePostTreatFF
    // (PaveFiller_6.cxx L3609-3635)
    // ====================================================================
    #[allow(clippy::too_many_arguments)]
    fn prepare_post_treat_ff(&mut self, a_int: usize, a_cur: usize, a_pb: &SharedPB,
                             a_ms_cpb: &mut Vec<(Shape, CoupleOfPBs)>,
                             a_mvi: &mut HashMap<(u64, u32), usize>, cid: usize) {
        // OCCT L3620: aLPBC.Append(aPB).
        if cid < self.ds.intersection_curves.len() {
            self.ds.intersection_curves[cid].pave_blocks.push(a_pb.clone());
        }
        let (n_v1, n_v2) = { let r = a_pb.0.read().unwrap(); r.indices() };
        let n_e = { let r = a_pb.0.read().unwrap(); r.edge };
        // OCCT L3626-3632: keep info for post treatment.
        let e_shape = self.ds.shape(n_e).clone();
        let mut a_cpb = CoupleOfPBs::new(a_int, a_cur);
        a_cpb.set_pb(a_pb.clone());
        a_ms_cpb.push((e_shape, a_cpb));
        a_mvi.insert((self.ds.shape(n_v1).ptr_id(), self.ds.shape(n_v1).location), n_v1);
        a_mvi.insert((self.ds.shape(n_v2).ptr_id(), self.ds.shape(n_v2).location), n_v2);
    }

    // ====================================================================
    // ProcessExistingPaveBlocks (after MakeEdge) — OCCT BOPAlgo_PaveFiller::
    // ProcessExistingPaveBlocks (PaveFiller_6.cxx L3072-3167)
    // ====================================================================
    #[allow(clippy::too_many_arguments)]
    fn process_existing_pave_blocks_1(
        &mut self,
        the_int: usize,
        the_cur: usize,
        n_f1: usize,
        n_f2: usize,
        n_e_es: usize,
        the_mpb_on_in: &[SharedPB],
        pb_candidates: &[usize],
        a_ms_cpb: &mut Vec<(Shape, CoupleOfPBs)>,
        a_mvi: &mut HashMap<(u64, u32), usize>,
        cid: usize,
        a_pb_faces_map: &mut crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<usize>>,
        a_mpb_add: &mut HashSet<u64>,
    ) {
        let a_box_es = self.edge_bbox(n_e_es);
        // OCCT L3090-3096: BoxTree selector for the section edge box.
        let mut sel_indices: Vec<usize> = Vec::new();
        for &i_pb in pb_candidates {
            let a_pbf = &the_mpb_on_in[i_pb];
            let n_e = { let r = a_pbf.0.read().unwrap(); r.edge };
            let a_si_e = self.ds.shape_info(n_e);
            if !a_si_e.bbox.is_out_box(&a_box_es) {
                sel_indices.push(i_pb);
            }
        }
        if sel_indices.is_empty() { return; }
        let a_tol_es = self.ds.edge_tolerance(n_e_es);
        for &i_pb in &sel_indices {
            let a_pbf = &the_mpb_on_in[i_pb];
            if a_mpb_add.contains(&pb_ptr(a_pbf)) { continue; }
            let b_in_f1 = self.pb_in_face(n_f1, a_pbf);
            let b_in_f2 = self.pb_in_face(n_f2, a_pbf);
            if b_in_f1 && b_in_f2 {
                // OCCT L3115-3118: add all common edges for post treatment.
                a_mpb_add.insert(pb_ptr(a_pbf));
                self.prepare_post_treat_ff(the_int, the_cur, a_pbf, a_ms_cpb, a_mvi, cid);
                continue;
            }
            // OCCT L3121-3124: myDistances lookup by (OriginalEdge, nF).
            let n_f = if b_in_f1 { n_f2 } else { n_f1 };
            let a_pbf_orig_edge = { let r = a_pbf.0.read().unwrap(); r.original_edge };
            let (a_t1, a_t2) = { let r = a_pbf.0.read().unwrap(); r.range() };
            let p_list = match self.my_distances.get(&(a_pbf_orig_edge, n_f)) {
                Some(l) => l,
                None => continue,
            };
            let mut a_dist = f64::INFINITY;
            for a_range_dist in p_list {
                if (a_t1 <= a_range_dist.first && a_range_dist.first <= a_t2)
                    || (a_t1 <= a_range_dist.last && a_range_dist.last <= a_t2)
                    || (a_range_dist.first <= a_t1 && a_t1 <= a_range_dist.last)
                    || (a_range_dist.first <= a_t2 && a_t2 <= a_range_dist.last)
                {
                    a_dist = a_range_dist.distance;
                    break;
                }
            }
            if a_dist < f64::INFINITY {
                let a_ef_edge = { let r = a_pbf.0.read().unwrap(); r.edge };
                let a_tol_sum = a_tol_es + self.ds.edge_tolerance(a_ef_edge);
                if a_dist <= a_tol_sum {
                    a_mpb_add.insert(pb_ptr(a_pbf));
                    self.prepare_post_treat_ff(the_int, the_cur, a_pbf, a_ms_cpb, a_mvi, cid);
                    let p_faces = a_pb_faces_map.bound(pb_ptr(a_pbf));
                    if !p_faces.contains(&n_f) {
                        p_faces.push(n_f);
                    }
                }
            }
        }
    }

    // ====================================================================
    // ProcessExistingPaveBlocks (with bound vertices) — OCCT BOPAlgo_PaveFiller::
    // ProcessExistingPaveBlocks (PaveFiller_6.cxx L3171-3274)
    // ====================================================================
    #[allow(clippy::too_many_arguments)]
    fn process_existing_pave_blocks_2(
        &mut self,
        the_int: usize,
        n_f1: usize,
        n_f2: usize,
        a_mpb_on_in: &[SharedPB],
        pb_candidates: &[usize],
        a_dm_bv: &crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>>,
        a_ms_cpb: &mut Vec<(Shape, CoupleOfPBs)>,
        a_mvi: &mut HashMap<(u64, u32), usize>,
        a_pb_faces_map: &mut crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<usize>>,
        a_mpb_add: &mut HashSet<u64>,
    ) {
        if a_dm_bv.is_empty() { return; }
        // OCCT L3194-3196: the FF's curves.
        let a_vc: Vec<usize> = if the_int < self.ds.interf_ff.len() {
            self.ds.interf_ff[the_int].curves.clone()
        } else {
            return;
        };
        let a_keys: Vec<usize> = a_dm_bv.iter_keys().collect();
        for i_c in a_keys {
            let a_lbv = match a_dm_bv.get(i_c) { Some(l) => l.clone(), None => continue };
            if i_c >= a_vc.len() { continue; }
            let cid = a_vc[i_c];
            for n_v in &a_lbv {
                let n_v = *n_v;
                let a_si_v = self.ds.shape_info(n_v);
                let a_box_v = a_si_v.bbox.clone();
                let a_v_shape = a_si_v.shape.clone();
                if !a_mvi.contains_key(&(a_v_shape.ptr_id(), a_v_shape.location)) {
                    continue;
                }
                // OCCT L3222-3228: BoxTree selector for the vertex box.
                let mut sel: Vec<usize> = Vec::new();
                for &i_pb in pb_candidates {
                    let a_pb = &a_mpb_on_in[i_pb];
                    let n_e = { let r = a_pb.0.read().unwrap(); r.edge };
                    let a_si_e = self.ds.shape_info(n_e);
                    if !a_si_e.bbox.is_out_box(&a_box_v) {
                        sel.push(i_pb);
                    }
                }
                for i_pb in sel {
                    let a_pb = &a_mpb_on_in[i_pb];
                    let (pv1, pv2) = { let r = a_pb.0.read().unwrap(); (r.pave1.vertex_idx, r.pave2.vertex_idx) };
                    if pv1 == n_v || pv2 == n_v { continue; }
                    if a_mpb_add.contains(&pb_ptr(a_pb)) { continue; }
                    let n_e = { let r = a_pb.0.read().unwrap(); r.edge };
                    // OCCT L3246: ComputeVE.
                    let (i_flag, _a_t, _a_tol) = self.my_context.compute_ve(
                        n_v, n_e, &self.ds, self.my_fuzzy_value);
                    if i_flag == 0 {
                        a_mpb_add.insert(pb_ptr(a_pb));
                        self.prepare_post_treat_ff(the_int, i_c, a_pb, a_ms_cpb, a_mvi, cid);
                        // OCCT L3252-3268: add faces to the PB.
                        let b_in_f1 = self.pb_in_face(n_f1, a_pb);
                        let b_in_f2 = self.pb_in_face(n_f2, a_pb);
                        if !b_in_f1 || !b_in_f2 {
                            let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                            let p_faces = a_pb_faces_map.bound(pb_ptr(a_pb));
                            if !p_faces.contains(&n_f) {
                                p_faces.push(n_f);
                            }
                        }
                    }
                }
            }
        }
    }

    // ====================================================================
    // Append a Shape (already containing its TShape data) into the DS with a box.
    // OCCT BOPDS_ShapeInfo::SetShapeType + SetShape + myDS->Append.
    // ====================================================================
    fn append_shape_with_box(&mut self, s: &Shape) -> usize {
        let idx = self.ds.append_shape(s.clone());
        match s.shape_type() {
            ShapeType::Vertex => {
                let pt = s.as_vertex().map(|v| v.point).unwrap_or(DVec3::ZERO);
                // OCCT BRepBndLib::Add(vertex) uses BRep_Tool::Tolerance
                // (clamped to Precision::Confusion minimum).
                let tol = s.as_vertex().map(|v| v.tolerance).unwrap_or(0.0)
                    .max(rcad_kernel::precision::CONFUSION);
                let mut b = BndBox::from_point(pt);
                b.set_gap(tol + rcad_kernel::CONFUSION);
                self.ds.change_shape_info(idx).bbox = b;
            }
            ShapeType::Edge => {
                self.rebuild_edge_box(idx);
            }
            _ => {}
        }
        idx
    }

    // ====================================================================
    // RemoveMicroSectionEdges — OCCT BOPAlgo_PaveFiller::RemoveMicroSectionEdges
    // (PaveFiller_6.cxx L4308-4384)
    // ====================================================================
    fn remove_micro_section_edges(&mut self, the_ms_cpb: &mut Vec<(Shape, CoupleOfPBs)>,
                                  the_micro_pb: &mut Vec<SharedPB>) {
        if the_ms_cpb.is_empty() { return; }
        // OCCT L4326-4376: build the new map of section edges avoiding micro edges.
        let mut a_se_pb_map: Vec<(Shape, CoupleOfPBs)> = Vec::new();
        let mut any_removed = false;
        let a_nb_cpb = the_ms_cpb.len();
        let mut removed: Vec<bool> = vec![false; a_nb_cpb];
        for (i, (a_si, a_cpb)) in the_ms_cpb.iter().enumerate() {
            if a_si.shape_type() != ShapeType::Edge {
                continue;
            }
            let a_pb = match a_cpb.pb1() {
                Some(pb) => pb.clone(),
                None => continue,
            };
            if a_pb.0.read().unwrap().has_edge() {
                continue;
            }
            // OCCT L4348: BOPTools_AlgoTools::IsMicroEdge(aSI, ctx, false) —
            // checked on the section edge shape, since the section PB carries no
            // edge reference yet (SetEdge is deferred to PostTreatFF).
            let n_e = { let r = a_pb.0.read().unwrap(); r.original_edge };
            let is_micro = self.is_micro_section_edge(a_si);
            if !is_micro {
                continue;
            }
            // Micro edge: remove from the FF intersection info.
            let a_ff_idx = a_cpb.index_interf;
            let a_cur_idx = a_cpb.index;
            let cid = if a_ff_idx < self.ds.interf_ff.len()
                && a_cur_idx < self.ds.interf_ff[a_ff_idx].curves.len()
            {
                self.ds.interf_ff[a_ff_idx].curves[a_cur_idx]
            } else {
                usize::MAX
            };
            if cid != usize::MAX && cid < self.ds.intersection_curves.len() {
                let a_lpbc = &mut self.ds.intersection_curves[cid].pave_blocks;
                if let Some(pos) = a_lpbc.iter().position(|p| pb_ptr(p) == pb_ptr(&a_pb)) {
                    a_lpbc.remove(pos);
                }
            }
            // OCCT L4376: add the micro PB for vertex unification.
            the_micro_pb.push(a_pb.clone());
            removed[i] = true;
            any_removed = true;
        }
        // OCCT L4380-4383: overwrite the old map if necessary.
        if any_removed {
            for (i, item) in the_ms_cpb.drain(..).enumerate() {
                if i < removed.len() && removed[i] {
                    continue;
                }
                a_se_pb_map.push(item);
            }
            *the_ms_cpb = a_se_pb_map;
        }
    }

    /// OCCT BOPTools_AlgoTools::IsMicroEdge — the edge has no valid shrunk range.
    fn is_micro_edge(&mut self, a_pb: &SharedPB) -> bool {
        self.fill_shrunk_data_pb(a_pb);
        let has_shrunk = { a_pb.0.read().unwrap().has_shrunk_data() };
        !has_shrunk
    }

    // ====================================================================
    // MakeSDVertices — OCCT BOPAlgo_PaveFiller::MakeSDVertices (PaveFiller_1.cxx L136-233)
    // ====================================================================
    pub(crate) fn make_sd_vertices(&mut self, the_vert_indices: &[usize], the_add_interfs: bool) -> usize {        // OCCT L143-161: collect vertices (resolving SD).
        let mut n_sd = usize::MAX;
        let mut a_vsd = None;
        let mut a_lv: Vec<Shape> = Vec::new();
        for &n_x in the_vert_indices {
            let mut n_sd1 = usize::MAX;
            if self.ds.has_shape_sd(n_x, &mut n_sd1) {
                let a_vsd1 = self.ds.shape(n_sd1).clone();
                if n_sd == usize::MAX {
                    a_vsd = Some(a_vsd1.clone());
                    n_sd = n_sd1;
                } else {
                    a_lv.push(a_vsd1);
                }
            }
            a_lv.push(self.ds.shape(n_x).clone());
        }
        // BOPTools_AlgoTools::MakeVertex — bounding vertex of the list.
        let (a_vn, a_vn_tol) = bounding_vertex(&a_lv);
        let (n_v, a_vn) = if n_sd != usize::MAX {
            // update old SD vertex with the new value.
            let a_vsd = a_vsd.unwrap();
            self.set_vertex_tolerance(n_sd, a_vn_tol);
            self.set_vertex_point(n_sd, a_vn);
            (n_sd, a_vsd)
        } else {
            let n_v = self.append_vertex(a_vn, a_vn_tol);
            (n_v, self.ds.shape(n_v).clone())
        };
        // OCCT L181-184: vertex box.
        let pt = a_vn.as_vertex().map(|v| v.point).unwrap_or(DVec3::ZERO);
        // OCCT BRepBndLib::Add(vertex) uses BRep_Tool::Tolerance (clamped).
        let tol = a_vn.as_vertex().map(|v| v.tolerance).unwrap_or(0.0)
            .max(rcad_kernel::precision::CONFUSION);
        let mut b = BndBox::from_point(pt);
        b.set_gap(tol + rcad_kernel::CONFUSION);
        self.ds.change_shape_info(n_v).bbox = b;
        // OCCT L186-231: fill ShapesSD + optional VV interferences.
        for (ii, &n1) in the_vert_indices.iter().enumerate() {
            self.ds.add_shape_sd(n1, n_v);
            if the_add_interfs {
                let i_r1 = self.ds.rank(n1);
                for &n2 in &the_vert_indices[ii + 1..] {
                    if i_r1 >= 0 && i_r1 == self.ds.rank(n2) {
                        self.my_report.add_warning(
                            crate::bop::algo::Alert::SelfInterferingShape(vec![n1, n2]));
                    }
                    if self.ds.add_interf(n1, n2) {
                        self.ds.interf_vv.push(crate::bop::ds::InterferenceVV {
                            v1: n1, v2: n2, merged_vertex: n_v,
                        });
                    }
                }
            }
        }
        n_v
    }

    /// Set a vertex's point (OCCT BRep_TVertex::Pnt).
    fn set_vertex_point(&mut self, idx: usize, pt: DVec3) {
        if idx < self.ds.shapes.len() {
            let si = self.ds.change_shape_info(idx);
            let ts = Arc::make_mut(&mut si.shape.data);
            if let topods::TShape::Vertex(ref mut vd) = *ts {
                vd.point = pt;
            }
            self.ds.remap_shape_idx(idx);
        }
    }

    // ====================================================================
    // MakeSDVerticesFF — OCCT BOPAlgo_PaveFiller::MakeSDVerticesFF
    // (PaveFiller_6.cxx L1141-1161)
    // ====================================================================
    fn make_sd_vertices_ff(&mut self, the_dm_vlv: &crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<usize>>,
                           the_dm_new_sd: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, usize>) {
        // OCCT L1190-1209: aItG(theDMVLV) — NCollection_DataMap<int, List<int>>
        // iterated in bucket order; each group makes one SD vertex.
        for (_n_v, a_list) in the_dm_vlv.iter() {
            let pts: Vec<(usize, DVec3)> = a_list.iter().map(|&i| (i, self.ds.vertex_point_by_idx(i))).collect();
            // OCCT L1152: MakeSDVertices(aList, false).
            let n_sd = self.make_sd_vertices(a_list, false);
            for &n_vx in a_list {
                the_dm_new_sd.insert(n_vx, n_sd);
            }
        }
    }

    // ====================================================================
    // PostTreatFF — OCCT BOPAlgo_PaveFiller::PostTreatFF (PaveFiller_6.cxx L1165-1669)
    // ====================================================================
    fn post_treat_ff(&mut self, the_ms_cpb: &mut Vec<(Shape, CoupleOfPBs)>,
                     a_dm_ex_edges: &mut crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<u64>>,
                     a_dm_new_sd: &mut crate::bop::algo::occt_map::OcctDataMapInt<usize, usize>,
                     the_micro_pb: &Vec<SharedPB>,
                     the_verts_on_rejected_pb: &[Shape]) {
        let a_nb_s = the_ms_cpb.len();
        if a_nb_s == 0 { return; }
        // OCCT L1203-1231: find unused vertices.
        let mut verts_unused: Vec<Shape> = Vec::new();
        let mut ind_map: HashSet<usize> = HashSet::new();
        let a_ffs = self.ds.interf_ff.clone();
        for ff in &a_ffs {
            let (n_f1, n_f2) = (ff.f1, ff.f2);
            let mut a_mv = HashSet::new();
            let mut a_mv_ef = HashSet::new();
            let mut a_mi = HashSet::new();
            self.get_stick_vertices(n_f1, n_f2, &mut a_mv, &mut a_mv_ef, &mut a_mi);
            let a_vc = ff.curves.clone();
            self.remove_used_vertices(&a_vc, &mut a_mv);
            let mv_list: Vec<usize> = a_mv.iter().copied().collect();
            for ind_v in mv_list {
                let a_vertex = self.ds.shape(ind_v).clone();
                if ind_map.insert(ind_v) {
                    verts_unused.push(a_vertex);
                } else {
                    verts_unused.retain(|s| {
                        !(s.ptr_id() == a_vertex.ptr_id() && s.location == a_vertex.location)
                    });
                }
            }
        }
        // OCCT L1234-1276: fast path (single shape, no micro/verts).
        let a_nb_me = the_micro_pb.len();
        let a_nb_v_on_rpb = the_verts_on_rejected_pb.len();
        if a_nb_s == 1 && a_nb_me == 0 && a_nb_v_on_rpb == 0 && verts_unused.is_empty() {
            let (a_s, a_cpb) = the_ms_cpb[0].clone();
            let a_type = a_s.shape_type();
            if a_type == ShapeType::Vertex {
                let i_v = self.append_shape_with_box(&a_s);
                let i_x = a_cpb.index_interf;
                let i_p = a_cpb.index;
                if i_x < self.ds.interf_ff.len() {
                    if let Some(np) = self.ds.interf_ff[i_x].points.get_mut(i_p) {
                        np.vertex_index = i_v;
                    }
                }
            } else if a_type == ShapeType::Edge {
                let a_pb1 = a_cpb.pb1().cloned();
                if let Some(a_pb1) = a_pb1 {
                    if a_pb1.0.read().unwrap().has_edge() {
                        a_dm_ex_edges.insert(pb_ptr(&a_pb1), vec![pb_ptr(&a_pb1)]);
                    } else {
                        // rcad: the section edge was already appended by append_edge
                        // in MakeBlocks; reuse that entry (OCCT appends aSx here —
                        // PostTreatFF L1536-1540).
                        let mut i_e = self.ds.index(&a_s);
                        if i_e < 0 {
                            i_e = self.append_shape_with_box(&a_s) as isize;
                        } else {
                            self.rebuild_edge_box(i_e as usize);
                        }
                        a_pb1.0.write().unwrap().edge = i_e as usize;
                    }
                }
            }
            return;
        }
        // OCCT L1279-1316: 1. prepare arguments for the fuse operation.
        let mut a_ls: Vec<Shape> = Vec::new();
        let mut an_added_sd: HashSet<(u64, u32)> = HashSet::new();
        let mut existing_edges: Vec<Shape> = Vec::new();
        for k in (0..a_nb_s).rev() {
            let (a_s, a_cpb) = &the_ms_cpb[k];
            let a_pb = a_cpb.pb1().cloned();
            if let Some(a_pb) = &a_pb {
                if a_pb.0.read().unwrap().has_edge() {
                    existing_edges.push(a_s.clone());
                } else {
                    a_ls.push(a_s.clone());
                }
                // OCCT L1297-1311: add SD-vertex candidates.
                let sub_shapes: Vec<Shape> = shape_sub_shapes(a_s);
                for a_ver in &sub_shapes {
                    let i_ver = self.ds.index(a_ver);
                    if i_ver < 0 { continue; }
                    let i_ver = i_ver as usize;
                    if let Some(&p_sd) = a_dm_new_sd.get(i_ver) {
                        let a_vsd = self.ds.shape(p_sd).clone();
                        if an_added_sd.insert((a_vsd.ptr_id(), a_vsd.location)) {
                            a_ls.push(a_vsd);
                        }
                    }
                }
            } else {
                a_ls.push(a_s.clone());
            }
        }
        if !existing_edges.is_empty() {
            // OCCT L1313-1316: add the existing edges as a compound.
            a_ls.push(Shape::new(Arc::new(TShape::Compound(existing_edges)), 0, topods::Orientation::Forward));
        }
        // OCCT L1324-1359: micro edges — add their vertices.
        for a_pb in the_micro_pb {
            let (n_v0, n_v1) = { let r = a_pb.0.read().unwrap(); r.indices() };
            let verts = [n_v0, n_v1];
            let mut a_verts = [self.ds.shape(verts[0]).clone(), self.ds.shape(verts[1]).clone()];
            for i in 0..2 {
                if let Some(&p_sd) = a_dm_new_sd.get(verts[i]) {
                    a_verts[i] = self.ds.shape(p_sd).clone();
                }
                if an_added_sd.insert((a_verts[i].ptr_id(), a_verts[i].location)) {
                    a_ls.push(a_verts[i].clone());
                }
            }
            // OCCT L1339-1358: ensure these vertices will be united.
            if a_verts[0].ptr_id() == a_verts[1].ptr_id() && a_verts[0].location == a_verts[1].location {
                continue;
            }
            let a_p1 = a_verts[0].as_vertex().map(|v| v.point).unwrap_or(DVec3::ZERO);
            let a_p2 = a_verts[1].as_vertex().map(|v| v.point).unwrap_or(DVec3::ZERO);
            // OCCT L1344-1347: BRep_Tool::Tolerance(aV1/2) — clamped to
            // Precision::Confusion minimum.
            let a_tol_v1 = a_verts[0].as_vertex().map(|v| v.tolerance).unwrap_or(0.0)
                .max(rcad_kernel::precision::CONFUSION);
            let a_tol_v2 = a_verts[1].as_vertex().map(|v| v.tolerance).unwrap_or(0.0)
                .max(rcad_kernel::precision::CONFUSION);
            let mut a_dist = a_p1.distance(a_p2);
            a_dist -= (a_tol_v1 + a_tol_v2);
            if a_dist > 0.0 {
                a_dist /= 2.0;
                self.set_vertex_tolerance_by_shape(&a_verts[0], a_tol_v1 + a_dist);
                self.set_vertex_tolerance_by_shape(&a_verts[1], a_tol_v2 + a_dist);
            }
        }
        // OCCT L1361-1385: add vertices put on the rejected section curves and unused vertices.
        for ver_map in [the_verts_on_rejected_pb, &verts_unused] {
            for a_ver in ver_map {
                let mut a_ver = a_ver.clone();
                let i_ver = self.ds.index(&a_ver);
                if i_ver >= 0 {
                    if let Some(&p_sd) = a_dm_new_sd.get(i_ver as usize) {
                        a_ver = self.ds.shape(p_sd).clone();
                    }
                }
                if an_added_sd.insert((a_ver.ptr_id(), a_ver.location)) {
                    a_ls.push(a_ver);
                }
            }
        }
        // OCCT L1389-1397: 2. Fuse shapes — run a nested PaveFiller.
        if a_ls.is_empty() {
            return;
        }
        let mut a_pf = PaveFiller::new();
        a_pf.my_run_parallel = self.my_run_parallel;
        a_pf.my_non_destructive = self.my_non_destructive;
        // OCCT PostTreatFF runs a non-primary BOPAlgo_PaveFiller to fuse the
        // section edges; it keeps the arguments by reference (no deep clone) so
        // aPDS->Index(aSx) resolves the section edges back to the main DS.
        a_pf.my_is_primary = false;
        a_pf.set_arguments(a_ls.clone());
        let prog = NoopProgress;
        let a_ps = ProgressScope::new(&prog, "Intersection of section edges", 100);
        a_pf.perform(&a_ps);
        if a_pf.has_errors() {
            self.my_report.add_error(crate::bop::algo::Alert::PostTreatFF);
            return;
        }
        // OCCT L1398-1405: tolerance cache for common blocks.
        // aMCBTol — compute-on-demand; aMEPB — one PB per intersection edge,
        // shared across all section shapes (OCCT L1450, outside the loop).
        let mut a_me_pb: HashMap<usize, SharedPB> = HashMap::new();
        // OCCT L1407-1656: map the fused shapes back into myDS.
        let a_ls2: Vec<Shape> = a_ls.clone();
        for a_sx in &a_ls2 {
            if a_sx.shape_type() == ShapeType::Compound {
                if let TShape::Compound(children) = &*a_sx.data {
                    for child in children {
                        a_ls.push(child.clone());
                    }
                }
                continue;
            }
            let n_sx = a_pf.ds().index(a_sx);
            if n_sx < 0 { continue; }
            let n_sx = n_sx as usize;
            let a_si_x = a_pf.ds().shape_info(n_sx).clone();
            let a_type = a_si_x.shape_type;
            if a_type == ShapeType::Vertex {
                let b_intersection_point = the_ms_cpb.iter().any(|(s, _)| {
                    s.ptr_id() == a_sx.ptr_id() && s.location == a_sx.location
                });
                let mut a_sd = usize::MAX;
                let a_v = if a_pf.ds().has_shape_sd(n_sx, &mut a_sd) {
                    a_pf.ds().shape(a_sd).clone()
                } else {
                    a_sx.clone()
                };
                let i_v = self.ds.index(&a_v);
                let i_v = if i_v < 0 {
                    self.append_shape_with_box(&a_v)
                } else {
                    i_v as usize
                };
                if !b_intersection_point {
                    let n_sx2 = self.ds.index(a_sx);
                    if n_sx2 >= 0 && (n_sx2 as usize) != i_v {
                        let n_sx2 = n_sx2 as usize;
                        a_dm_new_sd.insert(n_sx2, i_v);
                        self.ds.add_shape_sd(n_sx2, i_v);
                    }
                } else {
                    // update the FF interference point index.
                    if let Some((_, a_cpb)) = the_ms_cpb.iter().find(|(s, _)| {
                        s.ptr_id() == a_sx.ptr_id() && s.location == a_sx.location
                    }) {
                        let i_x = a_cpb.index_interf;
                        let i_p = a_cpb.index;
                        if i_x < self.ds.interf_ff.len() {
                            if let Some(np) = self.ds.interf_ff[i_x].points.get_mut(i_p) {
                                np.vertex_index = i_v;
                            }
                        }
                    }
                }
            } else if a_type == ShapeType::Edge {
                let b_has_pave_blocks = a_pf.ds().has_pave_blocks(n_sx);
                let a_cpb = match the_ms_cpb.iter().find(|(s, _)| {
                    s.ptr_id() == a_sx.ptr_id() && s.location == a_sx.location
                }) {
                    Some((_, cpb)) => cpb.clone(),
                    None => continue,
                };
                let i_x = a_cpb.index_interf;
                let i_c = a_cpb.index;
                let a_pb1 = a_cpb.pb1().cloned();
                let Some(a_pb1) = a_pb1 else { continue };
                let b_old = a_pb1.0.read().unwrap().has_edge();
                if b_old {
                    a_dm_ex_edges.bound(pb_ptr(&a_pb1));
                }
                if !b_has_pave_blocks {
                    if b_old {
                        a_dm_ex_edges.get_mut(pb_ptr(&a_pb1)).unwrap().push(pb_ptr(&a_pb1));
                    } else {
                        // rcad: the section edge was already appended to the DS by
                        // append_edge in MakeBlocks (a_sx is a clone of that entry),
                        // whereas OCCT appends aSx here for the first time
                        // (PostTreatFF L1536-1540). Reuse the existing entry to avoid
                        // a duplicate DS edge; fall back to appending if absent.
                        let mut i_e = self.ds.index(a_sx);
                        if i_e < 0 {
                            i_e = self.append_shape_with_box(a_sx) as isize;
                        } else {
                            self.rebuild_edge_box(i_e as usize);
                        }
                        a_pb1.0.write().unwrap().edge = i_e as usize;
                    }
                } else {
                    let a_lpbx: Vec<SharedPB> = a_pf.ds().pave_blocks(n_sx).to_vec();
                    let a_nb_lpbx = a_lpbx.len();
                    // micro edge check.
                    let is_micro = a_nb_lpbx == 0
                        || (a_nb_lpbx == 1 && !a_lpbx[0].0.read().unwrap().has_shrunk_data());
                    if is_micro {
                        // remove aPB1 from the curve's PB list.
                        let cid = if i_x < self.ds.interf_ff.len()
                            && i_c < self.ds.interf_ff[i_x].curves.len()
                        {
                            self.ds.interf_ff[i_x].curves[i_c]
                        } else {
                            usize::MAX
                        };
                        if cid != usize::MAX && cid < self.ds.intersection_curves.len() {
                            let a_lpbc = &mut self.ds.intersection_curves[cid].pave_blocks;
                            if let Some(pos) = a_lpbc.iter().position(|p| pb_ptr(p) == pb_ptr(&a_pb1)) {
                                a_lpbc.remove(pos);
                            }
                        }
                        // append the edge's vertices for SD.
                        if let Some(si) = self.ds.shapes.get(
                            self.ds.index(a_sx) as usize
                        ).cloned() {
                            for sub in &si.sub_shapes {
                                a_ls.push(self.ds.shape(*sub).clone());
                            }
                        }
                        continue;
                    }
                    if b_old && a_nb_lpbx == 0 {
                        a_dm_ex_edges.get_mut(pb_ptr(&a_pb1)).unwrap().push(pb_ptr(&a_pb1));
                        continue;
                    }
                    if !b_old {
                        let cid = if i_x < self.ds.interf_ff.len()
                            && i_c < self.ds.interf_ff[i_x].curves.len()
                        {
                            self.ds.interf_ff[i_x].curves[i_c]
                        } else {
                            usize::MAX
                        };
                        if cid != usize::MAX && cid < self.ds.intersection_curves.len() {
                            let a_lpbc = &mut self.ds.intersection_curves[cid].pave_blocks;
                            if let Some(pos) = a_lpbc.iter().position(|p| pb_ptr(p) == pb_ptr(&a_pb1)) {
                                a_lpbc.remove(pos);
                            }
                        }
                    }
                    if a_nb_lpbx > 0 {
                        for a_pbx in &a_lpbx {
                            let a_pbrx = a_pf.ds().real_pave_block(a_pbx);
                            let (a_pave0, a_pave1) = {
                                let r = a_pbrx.0.read().unwrap();
                                (r.pave1.clone(), r.pave2.clone())
                            };
                            let mut a_pave = [a_pave0, a_pave1];
                            for j in 0..2 {
                                let n_v = a_pave[j].vertex_idx;
                                let a_v = a_pf.ds().shape(n_v).clone();
                                let i_v = self.ds.index(&a_v);
                                let i_v = if i_v < 0 {
                                    self.append_shape_with_box(&a_v)
                                } else {
                                    i_v as usize
                                };
                                let a_p1 = if j == 0 {
                                    a_pb1.0.read().unwrap().pave1.clone()
                                } else {
                                    a_pb1.0.read().unwrap().pave2.clone()
                                };
                                if a_p1.vertex_idx != i_v {
                                    if (a_p1.param - a_pave[j].param).abs() < 1e-12 {
                                        a_dm_new_sd.insert(a_p1.vertex_idx, i_v);
                                        self.ds.add_shape_sd(a_p1.vertex_idx, i_v);
                                    } else {
                                        // check aPDS for the SD connection.
                                        let a_v_pave = self.ds.shape(a_p1.vertex_idx).clone();
                                        let n_v_new = a_pf.ds().index(&a_v_pave);
                                        if n_v_new >= 0 {
                                            let mut n_v_new_sd = usize::MAX;
                                            if a_pf.ds().has_shape_sd(n_v_new as usize, &mut n_v_new_sd)
                                                && n_v_new_sd == n_v
                                            {
                                                a_dm_new_sd.insert(a_p1.vertex_idx, i_v);
                                                self.ds.add_shape_sd(a_p1.vertex_idx, i_v);
                                            }
                                        }
                                    }
                                }
                                a_pave[j].set_index(i_v);
                            }
                            // add edge.
                            let n_e_pbrx = { a_pbrx.0.read().unwrap().edge };
                            let a_e = a_pf.ds().shape(n_e_pbrx).clone();
                            let i_e = self.ds.index(&a_e);
                            let i_e = if i_e < 0 {
                                self.append_shape_with_box(&a_e)
                            } else {
                                i_e as usize
                            };
                            // update the curve tolerance from the common block if any.
                            if a_pf.ds().is_common_block(&a_pbrx) {
                                if let Some(cb_idx) = a_pf.ds().common_block(&a_pbrx) {
                                    let a_cb = a_pf.ds().common_blocks[cb_idx].clone();
                                    let a_tol = a_cb.tolerance();
                                    if a_tol > 0.0 {
                                        let cid2 = if i_x < self.ds.interf_ff.len()
                                            && i_c < self.ds.interf_ff[i_x].curves.len()
                                        {
                                            self.ds.interf_ff[i_x].curves[i_c]
                                        } else {
                                            usize::MAX
                                        };
                                        if cid2 != usize::MAX && cid2 < self.ds.intersection_curves.len() {
                                            let a_nc = &mut self.ds.intersection_curves[cid2];
                                            if a_nc.tolerance < a_tol {
                                                a_nc.tolerance = a_tol;
                                            }
                                        }
                                    }
                                }
                            }
                            // OCCT L1628-1651: append new PaveBlock to aLPBC / aDMExEdges.
                            let p_pbc = a_me_pb.entry(i_e).or_insert_with(|| {
                                let a_pave_r1 = Pave::new(a_pave[0].vertex_idx, a_pave[0].param);
                                let a_pave_r2 = Pave::new(a_pave[1].vertex_idx, a_pave[1].param);
                                let mut pb = PaveBlock::new(i_e, a_pave_r1, a_pave_r2);
                                pb.original_edge = i_e;
                                let spb = SharedPB::new(pb);
                                // OCCT L1672-1696: the new PB replaces the removed aPB1
                                // in the curve's PB list / aDMExEdges; UpdateFaceInfo
                                // adds it by handle to PaveBlocksSc (L1778). rcad stores
                                // pool indices, so allocate a pool entry for it.
                                self.ds.pave_blocks_pool.entry(i_e).or_default().push(spb.clone());
                                spb
                            });
                            if b_old {
                                p_pbc.0.write().unwrap().set_original_edge(
                                    a_pb1.0.read().unwrap().original_edge);
                                a_dm_ex_edges.get_mut(pb_ptr(&a_pb1)).unwrap().push(pb_ptr(p_pbc));
                            } else {
                                let cid3 = if i_x < self.ds.interf_ff.len()
                                    && i_c < self.ds.interf_ff[i_x].curves.len()
                                {
                                    self.ds.interf_ff[i_x].curves[i_c]
                                } else {
                                    usize::MAX
                                };
                                if cid3 != usize::MAX && cid3 < self.ds.intersection_curves.len() {
                                    self.ds.intersection_curves[cid3].pave_blocks.push(p_pbc.clone());
                                }
                            }
                        }
                    }
                }
            }
        }
        // OCCT L1658-1668: update SD for vertices that did not participate.
        let keys: Vec<usize> = a_dm_new_sd.iter_keys().collect();
        for k in keys {
            if let Some(&v) = a_dm_new_sd.get(k) {
                if let Some(&v2) = a_dm_new_sd.get(v) {
                    a_dm_new_sd.insert(k, v2);
                    self.ds.add_shape_sd(k, v2);
                }
            }
        }
    }

    /// Set vertex tolerance by shape (used to raise tolerances before the fuse).
    fn set_vertex_tolerance_by_shape(&mut self, s: &Shape, tol: f64) {
        let idx = self.ds.index(s);
        if idx >= 0 {
            self.set_vertex_tolerance(idx as usize, tol);
        }
    }

    /// Pool key of a SharedPB (OCCT: PaveBlock handle → pool position).
    fn pb_pool_index(&self, pb: &SharedPB) -> Option<(usize, usize)> {
        let ptr = Arc::as_ptr(&pb.0);
        for (&key, pool) in self.ds.pave_blocks_pool.iter() {
            for (li, spb) in pool.iter().enumerate() {
                if Arc::as_ptr(&spb.0) == ptr {
                    return Some((key, li));
                }
            }
        }
        None
    }

    // ====================================================================
    // UpdateFaceInfo — OCCT BOPAlgo_PaveFiller::UpdateFaceInfo (PaveFiller_6.cxx L1673-1946)
    // ====================================================================
    fn update_face_info(&mut self, the_dm_e: &crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<u64>>,
                        the_dm_v: &crate::bop::algo::occt_map::OcctDataMapInt<usize, usize>,
                        the_pb_faces_map: &crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<usize>>) {
        // OCCT L1729: anEdgeLPB is NCollection_DataMap<int, List<PB>> —
        // bucket iteration order; key is the edge index.
        let mut an_edge_lpb: crate::bop::algo::occt_map::OcctDataMapInt<usize, Vec<u64>> =
            crate::bop::algo::occt_map::OcctDataMapInt::new();
        let a_ffs = self.ds.interf_ff.clone();
        let a_nb_ff = a_ffs.len();
        // OCCT L1726: aMF is NCollection_Map<int> — bucket iteration order
        // (used at L1919 to update the face info of the affected faces).
        let mut a_mf: crate::bop::algo::occt_map::OcctMapInt =
            crate::bop::algo::occt_map::OcctMapInt::new();
        // 1. Sections (curves, points).
        for i in 0..a_nb_ff {
            let (n_f1, n_f2) = (a_ffs[i].f1, a_ffs[i].f2);
            // 1.1. Section edges.
            let a_vnc = a_ffs[i].curves.clone();
            for cid in &a_vnc {
                if *cid >= self.ds.intersection_curves.len() { continue; }
                let old_pbs = self.ds.intersection_curves[*cid].pave_blocks.clone();
                let mut new_pbs: Vec<SharedPB> = Vec::new();
                for a_pb in &old_pbs {
                    let key = pb_ptr(a_pb);
                    // OCCT L1712-1731: treat existing pave blocks.
                    if let Some(a_lpb) = the_dm_e.get(key) {
                        // OCCT: UpdateExistingPaveBlocks(aPB, aLPB, thePBFacesMap).
                        let a_lpb_pbs: Vec<SharedPB> = a_lpb.iter().filter_map(|k| {
                            self.find_pb_by_key(*k)
                        }).collect();
                        self.update_existing_pave_blocks(a_pb, &a_lpb_pbs, the_pb_faces_map);
                        for pbe in &a_lpb_pbs {
                            let n_e = pbe.0.read().unwrap().edge;
                            an_edge_lpb.bound(n_e).push(pb_ptr(pbe));
                        }
                        continue; // removed from aLPBC
                    }
                    // OCCT L1733-1734: add section PB to both faces.
                    let n_e = a_pb.0.read().unwrap().edge;
                    // OCCT: ChangePaveBlocksSc().Add(aPB) — keyed by PB handle.
                    self.ds.change_face_info(n_f1).pave_blocks_sc.insert(key);
                    self.ds.change_face_info(n_f2).pave_blocks_sc.insert(key);
                    an_edge_lpb.bound(n_e).push(key);
                    new_pbs.push(a_pb.clone());
                }
                self.ds.intersection_curves[*cid].pave_blocks = new_pbs;
            }
            // 1.2. Section vertices.
            let points = a_ffs[i].points.clone();
            for np in &points {
                if np.vertex_index != usize::MAX {
                    let n_v1 = np.vertex_index;
                    self.ds.change_face_info(n_f1).vertices_sc.insert(n_v1);
                    self.ds.change_face_info(n_f2).vertices_sc.insert(n_v1);
                }
            }
            a_mf.add(n_f1);
            a_mf.add(n_f2);
        }
        // OCCT L1767-1858: create new common blocks from unified edge PBs.
        // OCCT anEdgeLPB (L1729) is NCollection_DataMap<int, List<PB>> —
        // iterated in bucket order (L1817 MakeCommonBlocks call).
        for (n_e, a_lpb_keys) in an_edge_lpb.iter() {
            if a_lpb_keys.len() == 1 { continue; }
            let mut a_cb_idx: Option<usize> = None;
            // OCCT L1831: aMFaces is NCollection_Map<int> — bucket iteration
            // order feeds SetFaces (L1896-1899).
            let mut a_m_faces: crate::bop::algo::occt_map::OcctMapInt =
                crate::bop::algo::occt_map::OcctMapInt::new();
            let mut a_mpave_blocks: Vec<SharedPB> = Vec::new();
            for &key in a_lpb_keys {
                let Some(a_pb) = self.find_pb_by_key(key) else { continue };
                if !a_mpave_blocks.iter().any(|p| pb_ptr(p) == key) {
                    a_mpave_blocks.push(a_pb.clone());
                }
                if self.ds.is_common_block(&a_pb) {
                    let pbcb_idx = self.ds.common_block(&a_pb).unwrap();
                    let a_pbcb = self.ds.common_blocks[pbcb_idx].clone();
                    for (p, _) in a_pbcb.pave_blocks() {
                        if !a_mpave_blocks.iter().any(|q| pb_ptr(q) == pb_ptr(p)) {
                            a_mpave_blocks.push(p.clone());
                        }
                    }
                    for &f in a_pbcb.faces() {
                        a_m_faces.add(f);
                    }
                    if a_cb_idx.is_none() {
                        a_cb_idx = Some(pbcb_idx);
                    }
                }
            }
            if a_cb_idx.is_none() {
                // OCCT L1821-1833: none of the PBs is a common block — create a new one.
                let a_pbs_ref: Vec<SharedPB> = a_mpave_blocks.clone();
                let cb_idx = self.ds.add_common_block(&a_pbs_ref);
                for a_pb in &a_pbs_ref {
                    self.ds.set_common_block(a_pb, cb_idx);
                }
            } else {
                // OCCT L1834-1856: update the existing common block.
                let cb_idx = a_cb_idx.unwrap();
                for a_pb in &a_mpave_blocks {
                    self.ds.set_common_block(a_pb, cb_idx);
                }
                let a_lpb_new: Vec<(SharedPB, usize)> =
                    a_mpave_blocks.iter().map(|p| (p.clone(), 0)).collect();
                self.ds.common_blocks[cb_idx].set_pave_blocks(a_lpb_new);
                let a_l_faces: Vec<usize> = a_m_faces.iter_keys().collect();
                self.ds.common_blocks[cb_idx].set_faces(a_l_faces);
            }
        }
        // OCCT L1860-1945: update face info with new vertices and PBs.
        let b_verts = !the_dm_v.is_empty();
        let b_edges = !the_dm_e.is_empty() || {
            let mut any = false;
            for cb in &self.ds.common_blocks {
                if cb.pave_blocks().len() > 1 { any = true; break; }
            }
            any
        };
        if !b_verts && !b_edges {
            return;
        }
        // OCCT L1919: aItMF.Initialize(aMF) — NCollection_Map bucket order.
        for n_f1 in a_mf.iter_keys() {
            // 2.1. update vertices.
            if b_verts {
                let mv_on = self.ds.change_face_info(n_f1).vertices_on.clone();
                let mv_in = self.ds.change_face_info(n_f1).vertices_in.clone();
                // OCCT L1698-1705: aDMNewSD is NCollection_DataMap<int,int>
                // — bucket iteration order.
                for (n_v1, &n_v2) in the_dm_v.iter() {
                    if mv_on.contains(&n_v1) {
                        self.ds.change_face_info(n_f1).vertices_on.remove(&n_v1);
                        self.ds.change_face_info(n_f1).vertices_on.insert(n_v2);
                    }
                    if mv_in.contains(&n_v1) {
                        self.ds.change_face_info(n_f1).vertices_in.remove(&n_v1);
                        self.ds.change_face_info(n_f1).vertices_in.insert(n_v2);
                    }
                }
            }
            // 2.2. update pave blocks.
            if b_edges {
                // OCCT L1906-1944: rebuild each PB set replacing PBs with their
                // RealPaveBlock (dedup via aMPBFence).
                let fi = self.ds.face_info(n_f1);
                let sets_copy = [
                    fi.pave_blocks_on.clone(),
                    fi.pave_blocks_in.clone(),
                    fi.pave_blocks_sc.clone(),
                ];
                drop(fi);
                let mut new_sets: Vec<IndexSet<u64>> = Vec::new();
                for copy in &sets_copy {
                    let mut a_mpb_fence: HashSet<u64> = HashSet::new();
                    let mut new_set: IndexSet<u64> = IndexSet::new();
                    for &pb_key in copy {
                        if let Some(a_pb) = self.ds.pb_from_ptr(pb_key) {
                            let rpb = self.ds.real_pave_block(&a_pb);
                            let rkey = pb_ptr(&rpb);
                            if a_mpb_fence.insert(rkey) {
                                // OCCT: Add(RealPaveBlock(aPB)) — PB handle key.
                                new_set.insert(rkey);
                            }
                        }
                    }
                    new_sets.push(new_set);
                }
                let fi = self.ds.change_face_info(n_f1);
                fi.pave_blocks_on = new_sets[0].clone();
                fi.pave_blocks_in = new_sets[1].clone();
                fi.pave_blocks_sc = new_sets[2].clone();
            }
        }
    }

    /// Resolve a PB key (u64 pointer) to a SharedPB.
    fn find_pb_by_key(&self, key: u64) -> Option<SharedPB> {
        for pool in self.ds.pave_blocks_pool.values() {
            for p in pool {
                if pb_ptr(p) == key {
                    return Some(p.clone());
                }
            }
        }
        None
    }

    // ====================================================================
    // UpdateExistingPaveBlocks — OCCT BOPAlgo_PaveFiller::UpdateExistingPaveBlocks
    // (PaveFiller_6.cxx L3278-3496)
    // ====================================================================
    fn update_existing_pave_blocks(&mut self, a_pbf: &SharedPB, a_lpb: &[SharedPB],
                                   the_pb_faces_map: &crate::bop::algo::occt_map::OcctDataMapInt<u64, Vec<usize>>) {
        if a_lpb.is_empty() { return; }
        // OCCT L3295-3324: 1. remove old pave blocks.
        let a_cb1 = self.ds.common_block(a_pbf);
        let b_cb = a_cb1.is_some();
        let mut a_lpb1: Vec<SharedPB> = Vec::new();
        if let Some(cb_idx) = a_cb1 {
            let cb = self.ds.common_blocks[cb_idx].clone();
            for (p, _) in cb.pave_blocks() {
                a_lpb1.push(p.clone());
            }
        } else {
            a_lpb1.push(a_pbf.clone());
        }
        // remove old PBs from the pool (by original edge).
        for a_pb1 in &a_lpb1 {
            let n_e = a_pb1.0.read().unwrap().original_edge;
            if let Some(pool) = self.ds.pave_blocks_pool.get_mut(&n_e) {
                if let Some(pos) = pool.iter().position(|p| pb_ptr(p) == pb_ptr(a_pb1)) {
                    pool.remove(pos);
                }
            }
        }
        // OCCT L3327-3446: 2. update pave blocks (create new common blocks).
        if b_cb {
            let cb1_idx = a_cb1.unwrap();
            let a_faces: Vec<usize> = self.ds.common_blocks[cb1_idx].faces().to_vec();
            let mut a_lpb_new: Vec<SharedPB> = Vec::new();
            for a_pb_value in a_lpb {
                let (vp0, vp1) = { let r = a_pb_value.0.read().unwrap(); (r.pave1.clone(), r.pave2.clone()) };
                let a_pb_value_paves = [vp0, vp1];
                for a_pb2 in &a_lpb1 {
                    let n_e = a_pb2.0.read().unwrap().original_edge;
                    let mut a_pb2n = PaveBlock::new(usize::MAX,
                        Pave::new(0, 0.0), Pave::new(0, 0.0));
                    let a_pb_value_oe = a_pb_value.0.read().unwrap().original_edge;
                    if a_pb_value_oe == n_e {
                        a_pb2n.pave1 = a_pb_value_paves[0];
                        a_pb2n.pave2 = a_pb_value_paves[1];
                    } else {
                        // compute paves for the different original edge.
                        let mut a_pave = [Pave::new(0, 0.0), Pave::new(0, 0.0)];
                        let (pb2v1, pb2v2) = { let r = a_pb2.0.read().unwrap(); (r.pave1.clone(), r.pave2.clone()) };
                        if a_pb_value_paves[0].vertex_idx == a_pb_value_paves[1].vertex_idx
                            && pb2v1.vertex_idx == pb2v2.vertex_idx
                        {
                            a_pave[0] = Pave::new(a_pb_value_paves[0].vertex_idx, pb2v1.param);
                            a_pave[1] = Pave::new(a_pb_value_paves[1].vertex_idx, pb2v2.param);
                        } else {
                            for i in 0..2 {
                                let n_v = a_pb_value_paves[i].vertex_idx;
                                a_pave[i] = Pave::new(n_v, 0.0);
                                if n_v == pb2v1.vertex_idx {
                                    a_pave[i].param = pb2v1.param;
                                } else if n_v == pb2v2.vertex_idx {
                                    a_pave[i].param = pb2v2.param;
                                } else {
                                    // project the vertex onto the original edge.
                                    let mut a_t_out = 0.0;
                                    let mut a_dist = 0.0;
                                    let (i_err, t, _) = self.my_context.compute_ve(
                                        n_v, n_e, &self.ds, self.my_fuzzy_value);
                                    a_t_out = t;
                                    a_dist = 0.0;
                                    if i_err == 0 {
                                        a_pave[i].param = a_t_out;
                                    } else {
                                        // closest boundary parameter.
                                        let p1 = self.ds.vertex_point_by_idx(pb2v1.vertex_idx);
                                        let p2 = self.ds.vertex_point_by_idx(pb2v2.vertex_idx);
                                        let pv = self.ds.vertex_point_by_idx(n_v);
                                        let d1 = pv.distance_squared(p1);
                                        let d2 = pv.distance_squared(p2);
                                        a_pave[i].param = if d1 < d2 { pb2v1.param } else { pb2v2.param };
                                    }
                                    let _ = a_dist;
                                }
                            }
                            if a_pave[1].param < a_pave[0].param {
                                a_pave.swap(0, 1);
                            }
                        }
                        a_pb2n.pave1 = a_pave[0];
                        a_pb2n.pave2 = a_pave[1];
                    }
                    a_pb2n.edge = a_pb_value.0.read().unwrap().edge;
                    a_pb2n.original_edge = n_e;
                    let spb = SharedPB::new(a_pb2n);
                    let cb_idx = self.ds.add_common_block(&[spb.clone()]);
                    self.ds.set_common_block(&spb, cb_idx);
                    self.ds.common_blocks[cb_idx].set_faces(a_faces.clone());
                    // myDS->ChangePaveBlocks(nE).Append(aPB2n) — IndexedDataMap
                    // grows on demand; the key may be usize::MAX ("no original edge").
                    self.ds.pave_blocks_pool.entry(n_e).or_default().push(spb.clone());
                }
                // aLPBNew.Append(aCB->PaveBlock1())
                let first = a_lpb1.first().cloned();
                if let Some(f) = first {
                    let key = pb_ptr(&f);
                    if let Some(found) = self.find_pb_by_key(key) {
                        a_lpb_new.push(found);
                    }
                }
            }
            let _ = a_lpb_new;
        } else {
            let n_e = a_pbf.0.read().unwrap().original_edge;
            for a_pb in a_lpb {
                self.ds.pave_blocks_pool.entry(n_e).or_default().push(a_pb.clone());
            }
        }
        // OCCT L3448-3496: project the edge on the faces.
        if let Some(p_l_faces) = the_pb_faces_map.get(pb_ptr(a_pbf)) {
            for &n_f in p_l_faces {
                for a_pb in a_lpb {
                    if self.pb_in_face(n_f, a_pb) {
                        continue;
                    }
                    // OCCT: IntTools_EdgeFace coincidence check → CommonBlock.
                    // rcad: approximate with the EF intersection.
                    let n_e = a_pb.0.read().unwrap().edge;
                    let (i_flag, _t, _tol) = self.my_context.compute_ef(
                        n_e, n_f, 0.0, 0.0, false, &self.ds, self.my_fuzzy_value);
                    if i_flag == 0 {
                        let cb_idx = if let Some(cb) = self.ds.common_block(a_pb) {
                            cb
                        } else {
                            let idx = self.ds.add_common_block(&[a_pb.clone()]);
                            self.ds.set_common_block(a_pb, idx);
                            idx
                        };
                        self.ds.common_blocks[cb_idx].add_face(n_f);
                        // OCCT: ChangePaveBlocksIn().Add(aPB) — PB handle key.
                        self.ds.change_face_info(n_f).pave_blocks_in
                            .insert(std::sync::Arc::as_ptr(&a_pb.0) as u64);
                    }
                }
            }
        }
    }

    // ====================================================================
    // UpdatePaveBlocks — OCCT BOPAlgo_PaveFiller::UpdatePaveBlocks (PaveFiller_6.cxx L3679-3811)
    // ====================================================================
    fn update_pave_blocks(&mut self, a_dm_new_sd: &crate::bop::algo::occt_map::OcctDataMapInt<usize, usize>) {
        if a_dm_new_sd.is_empty() { return; }
        let mut a_mpb: HashSet<u64> = HashSet::new();
        let mut a_micro_edges: HashSet<usize> = HashSet::new();
        // Collect all PBs: from section curves + the pool.
        let mut an_all_pbs: Vec<SharedPB> = Vec::new();
        let a_ffs = self.ds.interf_ff.clone();
        for ff in &a_ffs {
            for &cid in &ff.curves {
                if cid >= self.ds.intersection_curves.len() { continue; }
                an_all_pbs.extend(self.ds.intersection_curves[cid].pave_blocks.clone());
            }
        }
        // OCCT L3764-3766: myDS->ChangePaveBlocksPool() is a DynamicArray —
        // iterate the pool in ascending edge-key order (a HashMap values()
        // order would randomize the an_all_pbs sequence and hence the
        // split-edge indices created by SplitEdge below).
        let mut pb_keys: Vec<usize> = self.ds.pave_blocks_pool.keys().copied().collect();
        pb_keys.sort_unstable();
        for k in pb_keys {
            let a_lpb = match self.ds.pave_blocks_pool.get(&k) {
                Some(v) => v.clone(),
                None => continue,
            };
            an_all_pbs.extend(a_lpb);
        }
        for a_pb in &an_all_pbs {
            let mut a_pb = a_pb.clone();
            let a_cb = self.ds.common_block(&a_pb);
            let b_cb = a_cb.is_some();
            if let Some(cb_idx) = a_cb {
                if let Some(fpb) = self.ds.common_blocks[cb_idx].pave_block1() {
                    a_pb = fpb;
                }
            }
            if !a_mpb.insert(pb_ptr(&a_pb)) { continue; }
            let mut b_rebuild = false;
            let (mut n_v, mut a_t) = {
                let r = a_pb.0.read().unwrap();
                ([r.pave1.vertex_idx, r.pave2.vertex_idx], [r.pave1.param, r.pave2.param])
            };
            let was_regular_edge = n_v[0] != n_v[1];
            for j in 0..2 {
                if let Some(&sd) = a_dm_new_sd.get(n_v[j]) {
                    n_v[j] = sd;
                    b_rebuild = true;
                    let mut pbw = a_pb.0.write().unwrap();
                    let p = Pave::new(n_v[j], a_t[j]);
                    if j == 0 { pbw.pave1 = p; } else { pbw.pave2 = p; }
                }
            }
            if b_rebuild {
                let mut n_e = a_pb.0.read().unwrap().edge;
                if n_e == usize::MAX {
                    n_e = a_pb.0.read().unwrap().original_edge;
                }
                let is_deg_edge = n_e < self.ds.nb_shapes() && self.ds.shapes[n_e].has_flag();
                if was_regular_edge && !is_deg_edge && n_v[0] == n_v[1] {
                    self.fill_shrunk_data_pb(&a_pb);
                    if !a_pb.0.read().unwrap().has_shrunk_data() {
                        a_micro_edges.insert(n_e);
                        continue;
                    }
                }
                let n_sp = self.split_edge(n_e, n_v[0], a_t[0], n_v[1], a_t[1]);
                if n_sp != usize::MAX {
                    if let Some(cb_idx) = a_cb {
                        self.ds.common_blocks[cb_idx].set_edge(n_sp);
                    } else {
                        a_pb.0.write().unwrap().edge = n_sp;
                    }
                }
            }
        }
        if !a_micro_edges.is_empty() {
            self.remove_pave_blocks(&a_micro_edges);
        }
    }

    /// OCCT BOPAlgo_PaveFiller::SplitEdge (PaveFiller_7.cxx L553-585).
    fn split_edge(&mut self, n_e: usize, n_v1: usize, a_t1: f64, n_v2: usize, a_t2: f64) -> usize {
        let curve = match self.ds.edge_curve(n_e) {
            Some(c) => c.clone(),
            None => return usize::MAX,
        };
        let n_sp = self.ds.push_edge_inherit(curve, [a_t1, a_t2], n_v1, n_v2, Some(n_e));
        let a_tol = self.ds.edge_tolerance(n_e);
        self.ds.mutate_shape_data(n_sp, |ts| {
            if let topods::TShape::Edge(ed) = ts {
                ed.tolerance = a_tol;
            }
        });
        self.ds.remap_shape_idx(n_sp);
        self.rebuild_edge_box(n_sp);
        n_sp
    }

    // ====================================================================
    // RemovePaveBlocks — OCCT BOPAlgo_PaveFiller::RemovePaveBlocks (PaveFiller_6.cxx L3815-3915)
    // ====================================================================
    pub(crate) fn remove_pave_blocks(&mut self, the_edges: &HashSet<usize>) {
        // 1. from the Pave Blocks Pool.
        for pool in self.ds.pave_blocks_pool.values_mut() {
            pool.retain(|pb| !the_edges.contains(&pb.0.read().unwrap().edge));
        }
        // 2. from section curves.
        let a_ffs = self.ds.interf_ff.clone();
        for ff in &a_ffs {
            for &cid in &ff.curves {
                if cid >= self.ds.intersection_curves.len() { continue; }
                let old = self.ds.intersection_curves[cid].pave_blocks.clone();
                self.ds.intersection_curves[cid].pave_blocks = old.into_iter()
                    .filter(|pb| !the_edges.contains(&pb.0.read().unwrap().edge))
                    .collect();
            }
        }
        // 3. from Face Info.
        // OCCT removes the PB handles whose Edge() is in theEdges (L3951-3968);
        // rcad face info stores pool indices. Step 1 already emptied the pool
        // entries of the removed PBs, so a now-empty pool entry means every PB
        // it held referenced a removed edge — drop that index from all three sets.
        let a_nb_src = self.ds.nb_source_shapes();
        for i in 0..a_nb_src {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != ShapeType::Face || !a_si.has_reference() { continue; }
            // Collect the PB ids whose pave block no longer resolves.
            let to_remove: Vec<u64> = {
                let fi = self.ds.face_info(i);
                let mut tr = Vec::new();
                for pb_key in fi.pave_blocks_in.iter()
                    .chain(fi.pave_blocks_on.iter())
                    .chain(fi.pave_blocks_sc.iter())
                {
                    let pb_key = *pb_key;
                    // OCCT: PBs are handles — the pool entry cannot go empty
                    // while the PB exists; drop ids that no longer resolve.
                    if self.ds.pb_from_ptr(pb_key).is_none() {
                        tr.push(pb_key);
                    }
                }
                tr
            };
            let fi = self.ds.change_face_info(i);
            for pb_key in to_remove {
                fi.pave_blocks_in.remove(&pb_key);
                fi.pave_blocks_on.remove(&pb_key);
                fi.pave_blocks_sc.remove(&pb_key);
            }
        }
    }

    // ====================================================================
    // CorrectToleranceOfSE — OCCT BOPAlgo_PaveFiller::CorrectToleranceOfSE
    // (PaveFiller_6.cxx L4072-4273)
    // ====================================================================
    fn correct_tolerance_of_se(&mut self) {
        let a_ffs = self.ds.interf_ff.clone();
        let mut a_mpb: HashSet<u64> = HashSet::new();
        let mut a_mvi_to_reduce: HashSet<usize> = HashSet::new();
        // 1. iterate on all sections F-F.
        for ff in &a_ffs {
            for &cid in &ff.curves {
                if cid >= self.ds.intersection_curves.len() { continue; }
                let a_lpb = self.ds.intersection_curves[cid].pave_blocks.clone();
                for a_pb in &a_lpb {
                    let n_e = { let r = a_pb.0.read().unwrap(); r.edge };
                    if n_e == usize::MAX { continue; }
                    if !a_mpb.insert(pb_ptr(a_pb)) { continue; }
                    // OCCT L4110-4132: reduce the section edge tolerance.
                    let a_tol_c = self.ds.intersection_curves[cid].tolerance;
                    let a_tol_tang = self.ds.intersection_curves[cid].tang_tolerance;
                    let mut b_is_reduced = false;
                    if a_tol_c < a_tol_tang {
                        let a_tol_e = self.ds.edge_tolerance(n_e);
                        if a_tol_c < a_tol_e {
                            // In-place edit of the shared TShape (OCCT
                            // BRep_Builder::UpdateEdge) — same rationale as
                            // update_edge_tolerance: do not split an edge's
                            // TShape identity across its Location copies.
                            self.ds.mutate_shape_data(n_e, |ts| {
                                if let topods::TShape::Edge(ref mut ed) = *ts {
                                    ed.tolerance = a_tol_c;
                                }
                            });
                            self.ds.remap_shape_idx(n_e);
                            b_is_reduced = true;
                        }
                    }
                    // fill vertex -> PB map.
                    let (v0, v1) = { let r = a_pb.0.read().unwrap(); (r.pave1.vertex_idx, r.pave2.vertex_idx) };
                    for n_v in [v0, v1] {
                        let mut n_v = n_v;
                        self.ds.has_shape_sd(n_v, &mut n_v);
                        if b_is_reduced {
                            a_mvi_to_reduce.insert(n_v);
                        }
                    }
                }
            }
        }
        if a_mvi_to_reduce.is_empty() { return; }
        // 2. find the max tolerance of edges containing the vertices.
        let mut a_mvi_tol: HashMap<usize, f64> = HashMap::new();
        let mut a_mvi_pbs: HashMap<usize, Vec<SharedPB>> = HashMap::new();
        for a_lpb in self.ds.pave_blocks_pool.values().cloned() {
            for a_pb in &a_lpb {
                let n_e = { let r = a_pb.0.read().unwrap(); r.edge };
                if n_e == usize::MAX { continue; }
                let a_tol_e = self.ds.edge_tolerance(n_e);
                let (v0, v1) = { let r = a_pb.0.read().unwrap(); (r.pave1.vertex_idx, r.pave2.vertex_idx) };
                for n_v in [v0, v1] {
                    if a_mvi_to_reduce.contains(&n_v) {
                        let max = a_mvi_tol.entry(n_v).or_insert(a_tol_e);
                        if a_tol_e > *max { *max = a_tol_e; }
                        a_mvi_pbs.entry(n_v).or_default().push(a_pb.clone());
                    }
                }
            }
        }
        // 2.2 reduce tolerances if possible.
        for &n_v in &a_mvi_to_reduce {
            let a_v = self.ds.vertex_point_by_idx(n_v);
            let a_tol_v = self.ds.vertex_tolerance_by_idx(n_v);
            let a_max_tol = a_mvi_tol.get(&n_v).copied().unwrap_or(0.0);
            if a_tol_v - a_max_tol < 0.001 * a_tol_v { continue; }
            let mut a_max_tol = a_max_tol;
            let a_pbs = a_mvi_pbs.get(&n_v).cloned().unwrap_or_default();
            let mut fence: HashSet<u64> = HashSet::new();
            for a_pb in &a_pbs {
                if !fence.insert(pb_ptr(a_pb)) { continue; }
                let n_e = a_pb.0.read().unwrap().edge;
                let curve = match self.ds.edge_curve(n_e) { Some(c) => c.clone(), None => continue };
                for iPave in 0..2 {
                    let a_pave = if iPave == 0 {
                        a_pb.0.read().unwrap().pave1.clone()
                    } else {
                        a_pb.0.read().unwrap().pave2.clone()
                    };
                    let mut n_vsd = a_pave.vertex_idx;
                    self.ds.has_shape_sd(n_vsd, &mut n_vsd);
                    if n_vsd != n_v { continue; }
                    let a_p_on_e = curve.point_at(a_pave.param);
                    let a_dist = a_v.distance(a_p_on_e) + self.ds.edge_tolerance(n_e);
                    if a_dist > a_max_tol {
                        a_max_tol = a_dist;
                    }
                }
            }
            if a_max_tol < a_tol_v {
                self.set_vertex_tolerance(n_v, a_max_tol);
            }
        }
    }

    // ====================================================================
    // PutSEInOtherFaces — OCCT BOPAlgo_PaveFiller::PutSEInOtherFaces
    // (PaveFiller_6.cxx L4277-4304)
    // ====================================================================
    fn put_se_in_other_faces(&mut self) {
        // OCCT L4338: NCollection_IndexedMap<handle<BOPDS_PaveBlock>> aMPBScAll
        // — insertion order feeds ForceInterfEF → InterfEF array order.
        let mut a_mpb_sc_all: indexmap::IndexSet<(usize, usize)> = indexmap::IndexSet::new();
        let a_ffs = self.ds.interf_ff.clone();
        for ff in &a_ffs {
            for &cid in &ff.curves {
                if cid >= self.ds.intersection_curves.len() { continue; }
                let a_lpbc = self.ds.intersection_curves[cid].pave_blocks.clone();
                for a_pb in &a_lpbc {
                    let n_e = { let r = a_pb.0.read().unwrap(); r.edge };
                    if n_e != usize::MAX {
                        a_mpb_sc_all.insert((n_e, 0));
                    }
                }
            }
        }
        // OCCT L4303: ForceInterfEF(aMPBScAll, aPS.Next(), false).
        self.force_interf_ef_work(&a_mpb_sc_all, false);
    }
}
