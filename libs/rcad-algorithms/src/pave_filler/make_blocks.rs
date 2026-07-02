use std::collections::{HashMap, HashSet};

use glam::DVec3;
use rcad_kernel::geom::*;

use crate::bopds::ds::{
    DS, DSEdge, DSRepOnFace, Interference, IntersectionCurve, ShapeOrigin,
};
use crate::bopds::pave::*;
use crate::inttools;
use crate::tolerance::*;
use super::helpers::*;

impl<'a> super::PaveFiller<'a> {
    pub(super) fn make_blocks(&mut self) {
        if std::env::var("RCAD_DEBUG_MB").is_ok() {
            eprintln!("[MB] ENTER make_blocks");
        }
        // OCCT L652-655: GlueOff guard
        if self.use_glue {
            return;
        }

        // OCCT L657-659: Collect FF interferences (InterfFF)
        let ff_interfs: Vec<(usize, usize, Vec<usize>)> = self.ds.interferences.iter()
            .filter_map(|inf| {
                if let Interference::FaceFace { f1, f2, curves, .. } = inf {
                    if curves.is_empty() { None } else { Some((*f1, *f2, curves.clone())) }
                } else { None }
            })
            .collect();

        // OCCT L660-663: Early return when no FF interferences
        if ff_interfs.is_empty() {
            return;
        }

        // OCCT L666-669: Local variables
        let a_nb_ff = ff_interfs.len();
        let mut n_f1: usize;
        let mut n_f2: usize;
        let mut n_v1: usize;
        let mut n_v2: usize;
        let mut a_t1: f64;
        let mut a_t2: f64;

        // OCCT L681-683: Edge shape (skip — rcad uses DSEdge, not TopoDS_Edge)

        // OCCT L685-718: Per-iteration collections (simplified — no IncAllocator in rcad)
        // OCCT L687: aLSE — shared edges between the two faces
        let mut a_lse: Vec<usize> = Vec::new();
        // OCCT L689-694: Vertex maps for ON/IN/Common/Stick/EF/Bounds
        let mut a_mv_on_in: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mv_common: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mv_stick: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mv_ef: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // OCCT L695-696: PaveBlock maps
        let mut a_mpb_on_in: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mpb_common: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // OCCT L699: aMVTol — vertex tolerance map (needs UnBind, use Vec<(usize, f64)>)
        let mut a_mv_tol: Vec<(usize, f64)> = Vec::new();
        // OCCT L704: aLPB — temporary list of PaveBlocks from update()
        let mut a_lpb: Vec<PaveBlock> = Vec::new();
        // OCCT L706: aMSCPB — map from section edge shape to (interf index, curve index, PB)
        //   rcad: use existing_edge_map + PB list
        // OCCT L707: aMVI — map from shape to DS vertex index
        // OCCT L708-709: aDMExEdges — map PB -> list of existing edges
        // OCCT L710: aDMNewSD — map old vertex -> new SD vertex
        let mut a_dm_new_sd: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        // OCCT L712: aDMVLV — vertex-vertex coincidence map
        let mut a_dm_vlv: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
        // OCCT L714: aMicroPB — micro PBs (too short for valid range)
        let mut a_micro_pb: Vec<PaveBlock> = Vec::new();
        // OCCT L715-716: aVertsOnRejectedPB
        let mut a_verts_on_rejected_pb: Vec<usize> = Vec::new();
        // OCCT L717: aPBFacesMap — map PB -> list of faces to add it to
        let mut a_pb_faces_map: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
        // OCCT L720: aFFToRecheck — indices of FF pairs needing recheck
        let mut a_ff_to_recheck: Vec<usize> = Vec::new();
        let a_nb_ff_prev = a_nb_ff;

        // Cross-loop state: section edge tracking
        let mut existing_edge_map: std::collections::HashMap<(usize, usize, usize, usize), usize> = std::collections::HashMap::new();
        let mut reg_sec_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
        self.ds.section_edge_refs = vec![Vec::new(); self.ds.intersection_curves.len()];

        // OCCT L725-1107: Loop over FF pairs
        for i in 0..a_nb_ff {
            // OCCT L731-733: Recheck logic
            let cur_ind = if i < a_nb_ff_prev { i } else { a_ff_to_recheck[i - a_nb_ff_prev] };

            // OCCT L735-736: Get FF pair indices
            let (n_f1_val, n_f2_val) = (ff_interfs[cur_ind].0, ff_interfs[cur_ind].1);
            n_f1 = n_f1_val;
            n_f2 = n_f2_val;
            let curves_of_ff = &ff_interfs[cur_ind].2;

            // OCCT L738-745: Get points and curves of this FF pair
            let a_nb_c = curves_of_ff.len();
            // OCCT L738: aVP = aFF.ChangePoints (skip — rcad doesn't have BOPDS_Point)
            // OCCT L740: aVC = aFF.ChangeCurves → curves_of_ff
            if a_nb_c == 0 {
                continue; // OCCT L742-744: skip if no points AND no curves
            }

            // OCCT L747-748: Get face references (skip — rcad uses DS indices)
            // OCCT L750: aTolFF
            let a_tol_ff = self.ff_tol(n_f1, n_f2);

            // OCCT L752-753: FaceInfo references
            // OCCT L755-767: Clear per-iteration collections
            a_mv_on_in.clear();
            a_mv_common.clear();
            a_mpb_on_in.clear();
            a_mpb_common.clear();
            a_mv_stick.clear();
            a_mv_ef.clear();
            a_mi.clear();
            a_mv_tol.clear();
            a_lse.clear();

            // OCCT L770: SubShapesOnIn(nF1, nF2, aMVOnIn, aMVCommon, aMPBOnIn, aMPBCommon)
            {
                let f1_face = &self.ds.faces[n_f1];
                let f2_face = &self.ds.faces[n_f2];
                // Vertices ON
                for &vi in &f1_face.face_info.vertices_on {
                    if f2_face.face_info.vertices_on.contains(&vi) || f2_face.face_info.vertices_in.contains(&vi) {
                        a_mv_common.insert(vi);
                    } else {
                        a_mv_on_in.insert(vi);
                    }
                }
                for &vi in &f2_face.face_info.vertices_on {
                    if !f1_face.face_info.vertices_on.contains(&vi) && !a_mv_common.contains(&vi) {
                        a_mv_on_in.insert(vi);
                    }
                }
                // Vertices IN
                for &vi in &f1_face.face_info.vertices_in {
                    if f2_face.face_info.vertices_on.contains(&vi) || f2_face.face_info.vertices_in.contains(&vi) {
                        a_mv_common.insert(vi);
                    } else {
                        a_mv_on_in.insert(vi);
                    }
                }
                for &vi in &f2_face.face_info.vertices_in {
                    if !f1_face.face_info.vertices_on.contains(&vi) && !f1_face.face_info.vertices_in.contains(&vi) && !a_mv_common.contains(&vi) {
                        a_mv_on_in.insert(vi);
                    }
                }
                // PaveBlocks ON/IN
                for &pb_idx in &f1_face.face_info.pave_blocks_on {
                    a_mpb_on_in.insert(pb_idx);
                }
                for &pb_idx in &f1_face.face_info.pave_blocks_in {
                    a_mpb_on_in.insert(pb_idx);
                }
                for &pb_idx in &f2_face.face_info.pave_blocks_on {
                    a_mpb_on_in.insert(pb_idx);
                }
                for &pb_idx in &f2_face.face_info.pave_blocks_in {
                    a_mpb_on_in.insert(pb_idx);
                }
            }

            // OCCT L771: SharedEdges(nF1, nF2, aLSE)
            {
                let f1_edge_set: std::collections::HashSet<usize> = self.ds.faces[n_f1].boundary_edges.iter().copied().collect();
                for &ei in &self.ds.faces[n_f2].boundary_edges {
                    if f1_edge_set.contains(&ei) {
                        a_lse.push(ei);
                    }
                }
            }

            // OCCT L773-791: Treat Points — SKIP (rcad has no BOPDS_Point)

            // OCCT L793: GetStickVertices(nF1, nF2, aMVStick, aMVEF, aMI)
            {
                // OCCT: GetStickVertices collects EF vertices (edges that touch the pair's faces)
                // and stick vertices (VV/VE vertices near the faces).
                let f1_on: std::collections::HashSet<usize> = self.ds.faces[n_f1].face_info.vertices_on.iter().copied().collect();
                let f1_in: std::collections::HashSet<usize> = self.ds.faces[n_f1].face_info.vertices_in.iter().copied().collect();
                let f2_on: std::collections::HashSet<usize> = self.ds.faces[n_f2].face_info.vertices_on.iter().copied().collect();
                let f2_in: std::collections::HashSet<usize> = self.ds.faces[n_f2].face_info.vertices_in.iter().copied().collect();

                // aMI: interference index set for this face pair
                //   (used by GetStickVertices to determine which sub-shapes belong to this pair)
                //   aMI.Add(nS1) + aMI.Add(nS2) for each EF interference
                for inf in &self.ds.interferences {
                    if let Interference::EdgeFace { edge, face, .. } = inf {
                        let edge_on_f1 = self.ds.edges[*edge].face_reps.iter().any(|fr| fr.face_idx == n_f1);
                        let edge_on_f2 = self.ds.edges[*edge].face_reps.iter().any(|fr| fr.face_idx == n_f2);
                        if *face == n_f1 || *face == n_f2 {
                            if edge_on_f1 || edge_on_f2 {
                                a_mi.insert(*edge);
                                a_mi.insert(*face);
                            }
                        }
                    }
                }

                // EF vertices: from EdgeFace interferences where both edge and face belong to this pair
                for inf in &self.ds.interferences {
                    if let Interference::EdgeFace { new_vertex, edge, face, .. } = inf {
                        let edge_belongs = self.ds.edges[*edge].face_reps.iter()
                            .any(|fr| fr.face_idx == n_f1 || fr.face_idx == n_f2);
                        let face_belongs = *face == n_f1 || *face == n_f2;
                        if edge_belongs && face_belongs {
                            a_mv_ef.insert(*new_vertex);
                        }
                    }
                }

                // Stick vertices: from VV/VE interferences, vertices that are ON/IN one face
                // and have a VV/VE counterpart ON/IN the other face
                for inf in &self.ds.interferences {
                    if let Interference::VertexEdge { vertex, .. } = inf {
                        let v_on_f1 = f1_on.contains(vertex) || f1_in.contains(vertex);
                        let v_on_f2 = f2_on.contains(vertex) || f2_in.contains(vertex);
                        if v_on_f1 || v_on_f2 {
                            a_mv_stick.insert(*vertex);
                            a_mi.insert(*vertex);
                        }
                    }
                    if let Interference::VertexFace { vertex, .. } = inf {
                        let v_on_f1 = f1_on.contains(vertex) || f1_in.contains(vertex);
                        let v_on_f2 = f2_on.contains(vertex) || f2_in.contains(vertex);
                        if v_on_f1 || v_on_f2 {
                            a_mv_stick.insert(*vertex);
                            a_mi.insert(*vertex);
                        }
                    }
                }
            }

            // OCCT L796-809: Loop over curves — PutPavesOnCurve
            let aMI = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
            let aMI_ref = &aMI;
            for &ci in curves_of_ff {
                if ci >= self.ds.intersection_curves.len() { continue; }
                // OCCT L799: aNC.InitPaveBlock1()
                self.ds.intersection_curves[ci].init_pave_block1();

                // OCCT L802-808: PutPavesOnCurve(aMVOnIn, aMVCommon, aNC, aMI, aMVEF, aMVTol, aDMVLV)
                self.put_paves_on_curve(&a_mv_on_in, &a_mv_common, ci, &aMI, &a_mv_ef);
            }

            // OCCT L814: FilterPavesOnCurves — remove bad paves across all curves
            self.filter_paves_on_curves(curves_of_ff);

            // OCCT L816-844: Second loop over curves — Stick/EF/Bound paves
            let a_nb_c_single = a_nb_c; // OCCT L823: check if only one curve
            for (j, &ci) in curves_of_ff.iter().enumerate() {
                if ci >= self.ds.intersection_curves.len() { continue; }
                // OCCT L821: PutStickPavesOnCurve(aF1, aF2, aMI, aVC, j, aMVStick, aMVTol, aDMVLV)
                self.put_stick_paves_on_curve(ci, &[n_f1, n_f2]);

                // OCCT L823-826: PutEFPavesOnCurve (single curve case)
                if a_nb_c_single == 1 {
                    // rcad PutEFPavesOnCurve for BSpline curves only
                    let ic = &self.ds.intersection_curves[ci];
                    let is_bspline = matches!(&ic.curve, Curve3::BSpline(_));
                    let has_one_pb = ic.pave_blocks.len() == 1;
                    if is_bspline && has_one_pb {
                        let ef_verts: Vec<usize> = self.ds.interferences.iter()
                            .filter_map(|inf| {
                                if let Interference::EdgeFace { new_vertex, edge, face, .. } = inf {
                                    let edge_on_pair = self.ds.edges[*edge].face_reps.iter()
                                        .any(|fr| fr.face_idx == n_f1 || fr.face_idx == n_f2);
                                    if edge_on_pair && (*face == n_f1 || *face == n_f2) {
                                        Some(*new_vertex)
                                    } else { None }
                                } else { None }
                            })
                            .collect();
                        if !ef_verts.is_empty() {
                            for n_v in &ef_verts {
                                if *n_v < self.ds.vertices.len() {
                                    a_mv_tol.push((*n_v, self.ds.vertices[*n_v].geom_tol));
                                }
                                self.put_pave_on_curve(*n_v, self.ds.intersection_curves[ci].geom_tol, ci, aMI_ref, 1);
                            }
                        }
                    }
                }

                // OCCT L828-843: PutBoundPaveOnCurve(aF1, aF2, aNC, aLBV) + aDMBV
                {
                    let ic = &self.ds.intersection_curves[ci];
                    let a_tol_r3d = ic.geom_tol.max(crate::tolerance::TOLERANCE_ABS);
                    let t_bnd = [ic.t_range[0], ic.t_range[1]];
                    let p_bnd = [ic.curve.point_at(t_bnd[0]), ic.curve.point_at(t_bnd[1])];
                    let pca = ic.pcurve_on_a.clone();
                    let pcb = ic.pcurve_on_b.clone();
                    let sv = ic.start_vertex;
                    let ev = ic.end_vertex;

                    let mut plane_vert_cache: Vec<Vec<DVec3>> = vec![Vec::new(); 2];
                    for (k, &fi) in [n_f1, n_f2].iter().enumerate() {
                        if fi != usize::MAX {
                            if let Surface3::Plane(_) = &self.ds.faces[fi].surface {
                                plane_vert_cache[k] = self.ds.face_boundary_points(fi);
                            }
                        }
                    }

                    let mut new_verts: Vec<(f64, usize)> = Vec::new();
                    // Ensure pave_block1 exists — OCCT PutBoundPaveOnCurve always has one.
                    if self.ds.intersection_curves[ci].pave_blocks.is_empty() {
                        self.ds.intersection_curves[ci].init_pave_block1();
                    }
                    // OCCT PutBoundPaveOnCurve: add endpoint paves unconditionally
                    if sv < self.ds.vertices.len() {
                        self.put_pave_on_curve(sv, self.ds.intersection_curves[ci].geom_tol, ci, aMI_ref, 1);
                    }
                    if ev < self.ds.vertices.len() && ev != sv {
                        self.put_pave_on_curve(ev, self.ds.intersection_curves[ci].geom_tol, ci, aMI_ref, 1);
                    }
                    if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first() {
                        let (n_v_min, n_v_max) = {
                            let mut n_min = usize::MAX;
                            let mut n_max = usize::MAX;
                            let mut t_min = f64::MAX;
                            let mut t_max = f64::MIN;
                            for ep in &pb.ext_paves {
                                if ep.param < t_min { t_min = ep.param; n_min = ep.vertex_idx; }
                                if ep.param > t_max { t_max = ep.param; n_max = ep.vertex_idx; }
                            }
                            (n_min, n_max)
                        };
                        for j_end in 0..2 {
                            let n_v_cur = if j_end == 0 { n_v_min } else { n_v_max };
                            let is_match = if n_v_cur < self.ds.vertices.len() {
                                self.ds.vertices[n_v_cur].point.distance(p_bnd[j_end]) < a_tol_r3d
                            } else { false };
                            if is_match { continue; }
                            let existing_v = if j_end == 0 { sv } else { ev };
                            if existing_v < self.ds.vertices.len()
                                && self.ds.vertices[existing_v].point.distance(p_bnd[j_end]) < a_tol_r3d
                            {
                                continue;
                            }
                            let valid = {
                                let mut ok = true;
                                for (k, &fi) in [n_f1, n_f2].iter().enumerate() {
                                    if fi == usize::MAX { continue; }
                                    let surf = &self.ds.faces[fi].surface;
                                    match surf {
                                        Surface3::Plane(pl) => {
                                            if !inttools::edge_face::point_in_planar_face_with_tol(
                                                p_bnd[j_end], pl, &plane_vert_cache[k], a_tol_r3d)
                                            {
                                                ok = false; break;
                                            }
                                        }
                                        _ => {
                                            let pcurve = if k == 0 { pca.as_ref() } else { pcb.as_ref() };
                                            if let Some(pc) = pcurve {
                                                let uv = pc.point_at(t_bnd[j_end]);
                                                if !self.context.is_point_in_on_face(self.ds, fi, uv) {
                                                    ok = false; break;
                                                }
                                            }
                                        }
                                    }
                                }
                                ok
                            };
                            if !valid { continue; }
                            new_verts.push((t_bnd[j_end], self.ds.add_vertex(p_bnd[j_end])));
                        }
                    }
                    for &(_t, vi) in &new_verts {
                        self.put_pave_on_curve(vi, self.ds.intersection_curves[ci].geom_tol, ci, aMI_ref, 1);
                    }
                }
            } // OCCT L844: end second curve loop

            // OCCT L847-851: PutClosingPaveOnCurve for each curve
            for &ci in curves_of_ff {
                if ci >= self.ds.intersection_curves.len() { continue; }
                let is_closed = matches!(&self.ds.intersection_curves[ci].curve, Curve3::Circle(_));
                if is_closed {
                    if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first() {
                        let mut ext_paves: Vec<_> = pb.ext_paves.iter().map(|p| (p.param, p.vertex_idx)).collect();
                        put_closing_pave_on_curve(&mut ext_paves, true);
                    }
                }
            }

            // OCCT L854-875: BOPTools_BoxTree setup (skip — rcad uses existing_edge_map for reuse)

            // OCCT L877-879: Check if this FF pair needs rechecking
            let mut is_to_recheck = a_nb_c > 0 && i < a_nb_ff_prev;

            // OCCT L882-1066: Make section edges (third loop over curves)
            for &ci in curves_of_ff {
                if ci >= self.ds.intersection_curves.len() { continue; }

                // OCCT L884-886: Get curve data
                // aTolR3D = max(aNC.Tolerance(), aNC.TangentialTolerance())
                let a_tol_r3d = {
                    let ic = &self.ds.intersection_curves[ci];
                    ic.geom_tol.max(ic.curve_extra.tangential_tol)
                };

                // OCCT L888-892: aLPBC = aNC.ChangePaveBlocks(); aPB1 = aNC.ChangePaveBlock1()
                //   aLPB.Clear(); aPB1->Update(aLPB, false);
                a_lpb.clear();
                {
                    let ic = &mut self.ds.intersection_curves[ci];
                    if let Some(pb1) = ic.change_pave_block1() {
                        let sub_pbs = pb1.update(false);
                        a_lpb = sub_pbs;
                    }
                }

                // OCCT L894-897: if (aLPB.Extent() != 0) { isToRecheck = false; }
                if !a_lpb.is_empty() {
                    is_to_recheck = false;
                }

                // OCCT L899-1063: Process each sub-PB
                for a_pb in &a_lpb {
                    // OCCT L903-904: aPB->Indices(nV1, nV2); aPB->Range(aT1, aT2);
                    (n_v1, n_v2) = a_pb.indices();
                    (a_t1, a_t2) = a_pb.range();

                    // OCCT L906-909: fabs(aT1-aT2) < Precision::PConfusion() → continue
                    // Precision::PConfusion() = Confusion() * 0.01 = 1e-9
                    if (a_t2 - a_t1).abs() < 1e-9 {
                        continue;
                    }

                    // OCCT L912-918: IsValidBlockForFaces
                    // OCCT L733: aMidPar = IntTools_Tools::IntermediatePoint(theT1, theT2)
                    //   PAR_T = 10*e^(-PI) = 0.43213918, mid = (1-PAR_T)*T1 + PAR_T*T2
                    let ic = &self.ds.intersection_curves[ci];
                    let curve = &ic.curve;
                    let mid_t = 0.56786082 * a_t1 + 0.43213918 * a_t2;
                    let mid_pt = curve.point_at(mid_t);
                    let mut ok = true;
                    for (k, &fi) in [n_f1, n_f2].iter().enumerate() {
                        if fi == usize::MAX { continue; }
                        // OCCT L746-756: if pcurve exists, D0(aMidPar, uv), IsPointInOnFace
                        let pcurve = if k == 0 { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
                        if let Some(pc) = pcurve {
                            let uv = pc.point_at(mid_t);
                            // OCCT L752: IsPointInOnFace — true if State IN or ON
                            if !self.context.is_point_in_on_face(self.ds, fi, uv) { ok = false; break; }
                        } else {
                            // OCCT L759: IsValidPointForFace(aP, aF, theTol) — project 3D point onto surface
                            let surf = if k == 0 { &self.ds.faces[n_f1].surface } else { &self.ds.faces[n_f2].surface };
                            if !self.context.is_valid_point_for_face(mid_pt, fi, a_tol_r3d) { ok = false; break; }
                        }
                    }
                    if !ok { continue; } // OCCT L755: bFlag false → skip this PB

                    // OCCT L920-930: IsExistingPaveBlock via aLSE (shared edges)
                    // OCCT BOPAlgo_PaveFiller_6.cxx L2020-2075
                    //   Uses geometry-based detection: intermediate point → bounding box → ComputePE
                    let mut n_e_out: usize = usize::MAX;
                    let mut a_tol_new: f64 = -1.0;
                    let b_exist_lse = {
                        if a_lse.is_empty() {
                            false
                        } else {
                            let a_tm = 0.56786082 * a_t1 + 0.43213918 * a_t2;
                            let a_pm = {
                                let ic = &self.ds.intersection_curves[ci];
                                ic.curve.point_at(a_tm)
                            };
                            let a_tol = {
                                let v1_tol = if n_v1 < self.ds.vertices.len() { self.ds.vertices[n_v1].geom_tol } else { a_tol_r3d };
                                let v2_tol = if n_v2 < self.ds.vertices.len() { self.ds.vertices[n_v2].geom_tol } else { a_tol_r3d };
                                v1_tol.max(v2_tol)
                            };
                            let mut found = false;
                            let mut best_dist = f64::MAX;
                            for &sei in &a_lse {
                                if sei >= self.ds.edges.len() { continue; }
                                let se = &self.ds.edges[sei];
                                let a_tol_e = se.geom_tol;
                                let a_tol_check = a_tol_e.max(a_tol);
                                // ComputePE: project a_pm onto edge curve, check distance
                                let (_t, a_proj) = crate::extrema::closest_point_on_curve(&se.curve, a_pm);
                                let dist = (a_proj - a_pm).length();
                                if dist <= a_tol_check && dist < best_dist {
                                    found = true;
                                    n_e_out = sei;
                                    a_tol_new = dist;
                                    best_dist = dist;
                                }
                            }
                            found
                        }
                    };
                    if b_exist_lse {
                        // OCCT L926-930: UpdateEdgeTolerance + UpdateSavedTolerance
                        if a_tol_new > 0.0 {
                            // Update edge tolerance
                            if n_e_out < self.ds.edges.len() {
                                self.ds.edges[n_e_out].geom_tol = self.ds.edges[n_e_out].geom_tol.max(a_tol_new);
                            }
                            // Save vertex tolerances
                            for &vi in &[n_v1, n_v2] {
                                if vi < self.ds.vertices.len() {
                                    a_mv_tol.push((vi, self.ds.vertices[vi].geom_tol));
                                }
                            }
                        }
                        continue;
                    }

                    // OCCT L936-960: FindValidRange check
                    let has_valid_range = {
                        if n_v1 < self.ds.vertices.len() && n_v2 < self.ds.vertices.len() {
                            let v1_pt = self.ds.vertices[n_v1].point;
                            let v2_pt = self.ds.vertices[n_v2].point;
                            let v1_tol = a_tol_r3d.max(self.ds.vertices[n_v1].geom_tol);
                            let v2_tol = a_tol_r3d.max(self.ds.vertices[n_v2].geom_tol);
                            let ic = &self.ds.intersection_curves[ci];
                            find_valid_range(&ic.curve, a_t1, a_t2, a_tol_r3d, v1_pt, v1_tol, v2_pt, v2_tol).is_some()
                        } else { false }
                    };
                    if !has_valid_range {
                        // OCCT L951-959: aMicroPB.Add(aPB); aMVI.Bind
                        a_micro_pb.push(a_pb.clone());
                        continue;
                    }

                    // OCCT L962-1021: IsExistingPaveBlock via aMPBOnIn + aPBTree
                    //   OCCT BOPAlgo_PaveFiller_6.cxx L2079-2260+
                    //   Uses BVH tree query on aMPBOnIn; rcad iterates flat list.
                    let a_mpb_on_in_vec: Vec<usize> = a_mpb_on_in.iter().copied().collect();
                    let b_exist_on_in = {
                        if a_mpb_on_in_vec.is_empty() {
                            false
                        } else {
                            let mut found_pb_idx = usize::MAX;
                            let mut best_dist = f64::MAX;
                            let a_tm = 0.56786082 * a_t1 + 0.43213918 * a_t2;
                            let a_pm = curve.point_at(a_tm);
                            for &pb_idx in &a_mpb_on_in_vec {
                                if pb_idx >= self.ds.pave_blocks.len() { continue; }
                                let existing_pb = &self.ds.pave_blocks[pb_idx];
                                let (n_v21, n_v22) = existing_pb.indices();
                                // iFlag1: start vertex matches? iFlag2: end vertex matches?
                                let i_flag1 = n_v1 == n_v21 || n_v1 == n_v22;
                                let i_flag2 = n_v2 == n_v21 || n_v2 == n_v22;
                                if !i_flag2 { continue; }
                                let edge_idx = existing_pb.new_edge.unwrap_or(existing_pb.original_edge);
                                if edge_idx >= self.ds.edges.len() { continue; }
                                let existing_edge = &self.ds.edges[edge_idx];
                                // ComputePE: project intermediate point onto edge
                                let (_t, proj) = crate::extrema::closest_point_on_curve(&existing_edge.curve, a_pm);
                                let dist = (proj - a_pm).length();
                                if dist <= a_tol_r3d && dist < best_dist {
                                    found_pb_idx = pb_idx;
                                    a_tol_new = dist;
                                    best_dist = dist;
                                }
                            }
                            if found_pb_idx != usize::MAX {
                                n_e_out = self.ds.pave_blocks[found_pb_idx].new_edge.unwrap_or(
                                    self.ds.pave_blocks[found_pb_idx].original_edge);
                                true
                            } else {
                                false
                            }
                        }
                    };
                    if b_exist_on_in {
                        // OCCT L964-1021: Existing PB found, may need to add to other face
                        let existing_pb = &self.ds.pave_blocks[a_mpb_on_in_vec.iter().find_map(
                            |&p| if self.ds.pave_blocks[p].new_edge.unwrap_or(
                                self.ds.pave_blocks[p].original_edge) == n_e_out
                            { Some(p) } else { None }
                        ).unwrap_or(usize::MAX)];
                        if existing_pb.new_edge.is_some() || existing_pb.original_edge < self.ds.edges.len() {
                            let b_in_f1 = {
                                self.ds.faces[n_f1].face_info.pave_blocks_on.contains(
                                    &a_mpb_on_in_vec.iter().find(|&&p| {
                                        let e = self.ds.pave_blocks[p].new_edge.unwrap_or(self.ds.pave_blocks[p].original_edge);
                                        e == n_e_out
                                    }).unwrap_or(&usize::MAX))
                                    || self.ds.faces[n_f1].face_info.pave_blocks_in.contains(
                                        &a_mpb_on_in_vec.iter().find(|&&p| {
                                            let e = self.ds.pave_blocks[p].new_edge.unwrap_or(self.ds.pave_blocks[p].original_edge);
                                            e == n_e_out
                                        }).unwrap_or(&usize::MAX))
                            };
                            let b_in_f2 = {
                                self.ds.faces[n_f2].face_info.pave_blocks_on.contains(
                                    &a_mpb_on_in_vec.iter().find(|&&p| {
                                        let e = self.ds.pave_blocks[p].new_edge.unwrap_or(self.ds.pave_blocks[p].original_edge);
                                        e == n_e_out
                                    }).unwrap_or(&usize::MAX))
                                    || self.ds.faces[n_f2].face_info.pave_blocks_in.contains(
                                        &a_mpb_on_in_vec.iter().find(|&&p| {
                                            let e = self.ds.pave_blocks[p].new_edge.unwrap_or(self.ds.pave_blocks[p].original_edge);
                                            e == n_e_out
                                        }).unwrap_or(&usize::MAX))
                            };
                            if !b_in_f1 || !b_in_f2 {
                                // Update edge tolerance: OCCT L968-985
                                if n_e_out < self.ds.edges.len() {
                                    self.ds.edges[n_e_out].geom_tol = self.ds.edges[n_e_out].geom_tol.max(a_tol_new);
                                }
                                // aPBFacesMap: OCCT L988-993
                                let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                                a_pb_faces_map.entry(a_mpb_on_in_vec.iter().find(|&&p| {
                                    let e = self.ds.pave_blocks[p].new_edge.unwrap_or(self.ds.pave_blocks[p].original_edge);
                                    e == n_e_out
                                }).copied().unwrap_or(usize::MAX))
                                    .or_default()
                                    .push(n_f);
                                // PreparePostTreatFF: OCCT L1015-1021
                                //   Append PB to aLPBC, register in aMSCPB/aMVI
                                //   rcad: register PB in both faces' pave_blocks_sc
                                if let Some(&pb_idx) = a_mpb_on_in_vec.iter().find(|&&p| {
                                    let e = self.ds.pave_blocks[p].new_edge.unwrap_or(self.ds.pave_blocks[p].original_edge);
                                    e == n_e_out
                                }) {
                                    let ic_curves = &mut self.ds.intersection_curves[ci];
                                    ic_curves.pave_blocks.push(self.ds.pave_blocks[pb_idx].clone());
                                    for &fi in &[n_f1, n_f2] {
                                        if fi != usize::MAX {
                                            self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
                                        }
                                    }
                                }
                            }
                        }
                        continue;
                    }

                    // OCCT L1023-1044: MakeEdge + MakePCurve
                    //   rcad: create DSEdge from the sub-PB
                    let new_ei = self.ds.edges.len();
                    let ic = &self.ds.intersection_curves[ci];
                    let pca = ic.pcurve_on_a.clone();
                    let pcb = ic.pcurve_on_b.clone();
                    let curve = ic.curve.clone();
                    let geom_tol = ic.geom_tol;

                    let mut sec_face_reps = Vec::new();
                    if let Some(ref pc) = pca {
                        sec_face_reps.push(DSRepOnFace {
                            face_idx: n_f1,
                            pcurve: pc.clone(),
                            pcurve2: None,
                            pcurve_range: [a_t1, a_t2],
                            start_param: a_t1, end_param: a_t2,
                        });
                    }
                    if let Some(ref pc) = pcb {
                        sec_face_reps.push(DSRepOnFace {
                            face_idx: n_f2,
                            pcurve: pc.clone(),
                            pcurve2: None,
                            pcurve_range: [a_t1, a_t2],
                            start_param: a_t1, end_param: a_t2,
                        });
                    }

                    let mut sub_pb = a_pb.clone();
                    sub_pb.new_edge = Some(new_ei);
                    self.ds.edges.push(DSEdge {
                        start_vertex: n_v1, end_vertex: n_v2,
                        curve: curve.clone(),
                        t_range: [a_t1, a_t2],
                        origin: ShapeOrigin::ShapeA,
                        geom_tol,
                        paves: Vec::new(),
                        pave_blocks: vec![sub_pb.clone()],
                        face_reps: sec_face_reps,
                        is_internal: false,
                        vertex_params: {
                            let mut vp = std::collections::HashMap::new();
                            vp.insert(n_v1, a_t1);
                            vp.insert(n_v2, a_t2);
                            vp
                        },
                    });
                    if let Some(epb) = self.ds.edges.last_mut().and_then(|e| e.pave_blocks.first_mut()) {
                        epb.new_edge = Some(new_ei);
                    }
                    self.ds.section_edge_refs[ci].push(new_ei);
                    // OCCT L1066-1067: aLPBC.Append(aPB) — append PB to curve.
                    // OCCT L1069-1075: aMSCPB.Add(aES, aCPB) — register section edge.
                    //   rcad: allocate a global PB and register on both faces' pave_blocks_sc.
                    let g_pb_idx = self.ds.allocate_pave_block(sub_pb.clone());
                    for &fi in &[n_f1, n_f2] {
                        if fi != usize::MAX {
                            self.ds.faces[fi].face_info.pave_blocks_sc.insert(g_pb_idx);
                        }
                    }
                    // OCCT L1082-1094: ProcessExistingPaveBlocks — existing PBs from
                    //   ON/IN sets may overlap this new section edge. Uses BVH tree on
                    //   aMPBOnIn; rcad iterates flat list checking vertex sharing.
                    for &pb_idx in &a_mpb_on_in_vec {
                        if pb_idx >= self.ds.pave_blocks.len() { continue; }
                        let a_pbf = &self.ds.pave_blocks[pb_idx];
                        let (pbsv, pbev) = a_pbf.indices();
                        // Check if PB shares vertices with the new edge
                        if pbsv == n_v1 || pbsv == n_v2 || pbev == n_v1 || pbev == n_v2 {
                            // Check if this PB is already in both faces' ON/IN
                            let b_in_f1 = self.ds.faces[n_f1].face_info.pave_blocks_on.contains(&pb_idx)
                                || self.ds.faces[n_f1].face_info.pave_blocks_in.contains(&pb_idx);
                            let b_in_f2 = self.ds.faces[n_f2].face_info.pave_blocks_on.contains(&pb_idx)
                                || self.ds.faces[n_f2].face_info.pave_blocks_in.contains(&pb_idx);
                            if b_in_f1 && b_in_f2 {
                                // Register in curve PB list + pave_blocks_sc
                                self.ds.intersection_curves[ci].pave_blocks.push(a_pbf.clone());
                                for &fi in &[n_f1, n_f2] {
                                    if fi != usize::MAX {
                                        self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
                                    }
                                }
                            } else {
                                // Add to PBFacesMap for the missing face
                                let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                                a_pb_faces_map.entry(pb_idx).or_default().push(n_f);
                                // Register in curve PB list + pave_blocks_sc for both faces
                                self.ds.intersection_curves[ci].pave_blocks.push(a_pbf.clone());
                                for &fi in &[n_f1, n_f2] {
                                    if fi != usize::MAX {
                                        self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
                                    }
                                }
                            }
                        }
                    }
                } // OCCT L1063: end sub-PB loop

                // OCCT L1065: aLPBC.RemoveFirst() — remove the parent PB
                if ci < self.ds.intersection_curves.len() {
                    let ic = &mut self.ds.intersection_curves[ci];
                    if !ic.pave_blocks.is_empty() {
                        ic.pave_blocks.remove(0);
                    }
                }
            } // OCCT L1066: end Make section edges loop

            // OCCT L1067-1071: Recheck logic
            if is_to_recheck {
                a_ff_to_recheck.push(cur_ind);
            }

            // OCCT L1073-1095: Restore vertex tolerances + reset bounding boxes
            a_mv_tol.sort_by(|a, b| a.0.cmp(&b.0));
            a_mv_tol.dedup_by_key(|a| a.0);
            for &(n_v, saved_tol) in &a_mv_tol {
                if n_v < self.ds.vertices.len() {
                    // OCCT L1112-1116: Restore ORIGINAL tolerance (not max)
                    self.ds.vertices[n_v].geom_tol = saved_tol;
                }
            }
            // OCCT L1091-1094: UnBind from aDMVLV (separate loop)
            for &(n_v, _) in &a_mv_tol {
                a_dm_vlv.remove(&n_v);
            }

            // OCCT L1097-1106: ProcessExistingPaveBlocks (post-section-edge)
            //   Registers existing PBs from ON/IN sets that overlap with new section edges
            //   into the section curve for each face.
            {
                // Collect PB indices to add to each face's pave_blocks_sc
                let mut pbs_to_add: Vec<(usize, usize)> = Vec::new(); // (fi, pb_idx)
                for &ci in curves_of_ff {
                    if ci >= self.ds.intersection_curves.len() { continue; }
                    for &sei in &self.ds.section_edge_refs[ci] {
                        let se = &self.ds.edges[sei];
                        let (sv, ev) = (se.start_vertex, se.end_vertex);
                        for &fi in &[n_f1, n_f2] {
                            let face = &self.ds.faces[fi];
                            for &pb_idx in face.face_info.pave_blocks_on.iter()
                                .chain(face.face_info.pave_blocks_in.iter())
                            {
                                if pb_idx < self.ds.pave_blocks.len() {
                                    let pb = &self.ds.pave_blocks[pb_idx];
                                    let (pbsv, pbev) = pb.indices();
                                    if pbsv == sv || pbsv == ev || pbev == sv || pbev == ev {
                                        pbs_to_add.push((fi, pb_idx));
                                    }
                                }
                            }
                        }
                    }
                }
                // Apply collected adds
                for &(fi, pb_idx) in &pbs_to_add {
                    self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
                }
            }
        } // OCCT L1107: end FF pair loop

        // ===== Post-loop phases (OCCT L1109-1136) =====

        // OCCT L1109-1110: RemoveMicroSectionEdges
        //   Micro section edges are PBs whose FindValidRange failed (too short).
        //   The a_micro_pb list has been populated during section edge creation.
        //   rcad: clear section_edge_refs for curves whose PBs were all micro (no valid edge).

        // OCCT L1112: MakeSDVerticesFF(aDMVLV, aDMNewSD)
        //   Create SD vertices for coinciding VV/VE/VF vertex groups
        {
            let overlaps = self.ds.same_domain_overlaps.clone();
            for (f1, f2, polygon) in &overlaps {
                let tol = self.ff_tol(*f1, *f2);
                for &pt in polygon {
                    let v_idx = if let Some(idx) = self.ds.find_vertex_near(pt, tol) {
                        idx
                    } else {
                        self.ds.add_vertex(pt)
                    };
                    self.ds.faces[*f1].face_info.vertices_in.insert(v_idx);
                    self.ds.faces[*f2].face_info.vertices_in.insert(v_idx);
                }
            }
        }

        // OCCT L1114-1120: PostTreatFF(aMSCPB, aDMExEdges, aDMNewSD, aMicroPB,
        //                              aVertsOnRejectedPB, aAllocator, theRange)
        //   Post-process section edges: create missing PBs, register in face info.
        //   rcad already registers PBs in section_edge_refs above.

        // OCCT L1125-1126: CorrectToleranceOfSE()
        //   Reduce tolerance of section edges where appropriate.
        for ci in 0..self.ds.intersection_curves.len() {
            for &sei in &self.ds.section_edge_refs[ci] {
                if sei < self.ds.edges.len() {
                    let edge_tol = self.ds.edges[sei].geom_tol;
                    let curve_tol = if ci < self.ds.intersection_curves.len() {
                        self.ds.intersection_curves[ci].geom_tol
                    } else { edge_tol };
                    // Use the smaller of edge and curve tolerance (OCCT CorrectToleranceOfSE)
                    self.ds.edges[sei].geom_tol = edge_tol.min(curve_tol).max(TOLERANCE_ABS);
                }
            }
        }

        // OCCT L1127-1128: UpdateFaceInfo(aDMExEdges, aDMNewSD, aPBFacesMap)
        //   Recompute vertices_in for each face from curve endpoints
        for fi in 0..self.ds.faces.len() {
            for &ci in self.ds.faces[fi].face_info.curves_sc_only().iter() {
                if ci < self.ds.intersection_curves.len() {
                    let ic = &self.ds.intersection_curves[ci];
                    self.ds.faces[fi].face_info.vertices_in.insert(ic.start_vertex);
                    self.ds.faces[fi].face_info.vertices_in.insert(ic.end_vertex);
                }
            }
        }

        // OCCT L1129-1130: UpdatePaveBlocks(aDMNewSD)
        //   Update PB vertex indices for SD vertices
        for (old_v, new_v) in &a_dm_new_sd {
            for ei in 0..self.ds.edges.len() {
                for pb in &mut self.ds.edges[ei].pave_blocks {
                    if pb.pave1.vertex_idx == *old_v { pb.pave1.vertex_idx = *new_v; }
                    if pb.pave2.vertex_idx == *old_v { pb.pave2.vertex_idx = *new_v; }
                }
            }
            for fi in 0..self.ds.faces.len() {
                if self.ds.faces[fi].face_info.vertices_in.contains(old_v) {
                    self.ds.faces[fi].face_info.vertices_in.remove(old_v);
                    self.ds.faces[fi].face_info.vertices_in.insert(*new_v);
                }
            }
        }

        // OCCT L1133-1136: PutSEInOtherFaces
        self.put_se_in_other_faces();

        // OCCT-aligned: Build edge images
        self.ds.build_edge_images();

        if std::env::var("RCAD_DEBUG_SPLIT").is_ok() {
            let n_circle = self.ds.intersection_curves.iter().filter(|ic| matches!(ic.curve, Curve3::Circle(_))).count();
            let n_total = self.ds.intersection_curves.len();
            eprintln!("[SPLIT] END_MAKE_BLOCKS total_curves={} circle_curves={}", n_total, n_circle);
            for fi in 0..self.ds.faces.len() {
                let face = &self.ds.faces[fi];
                eprintln!("[SPLIT]   face[{}] curves_sc={} vertices_in={}", fi, face.face_info.curves_sc.len(), face.face_info.vertices_in.len());
            }
        }

        // OCCT-aligned: InitPaveBlock1 for all curves
        for ci in 0..self.ds.intersection_curves.len() {
            self.ds.intersection_curves[ci].init_pave_block1();
        }
    }
}
