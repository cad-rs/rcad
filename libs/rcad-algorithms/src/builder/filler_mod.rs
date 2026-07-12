use super::{BooleanBuilder, SourceSide};
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use indexmap::IndexMap;
use glam::{DVec2, DVec3};
use rayon::prelude::*;
use rcad_kernel::topods;
use rcad_kernel::PCurve;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval, *};
use rcad_kernel::topology::*;
use crate::bvh::{Aabb, DsBvh};
use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::bopalgo::{GlueEnum, Alert, Report};
use crate::builder::types::*;
use crate::builder::wire_splitter::{EdgeInfo, build_closed_wires, world_to_uv};
use super::ResultBuilder;
use crate::history::{BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin};
use crate::inttools::context::Context;
use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;

impl<'a> BooleanBuilder<'a> {
    // ====================================================================
    // --OCCT-aligned: dimension-by-dimension pipeline (PerformInternal1)
    //   BOPAlgo_Builder.cxx L310-440
    // ====================================================================

    /// --OCCT-aligned: FillImagesVertices (BOPAlgo_Builder_1.cxx L40-67).
    ///   Iterates ShapesSD --builds myImages(VERTEX) + myShapesSD + myOrigins.
    /// ✅ OCCT-aligned: FillImagesVertices (BOPAlgo_Builder_1.cxx L42-65).
    /// Maps each SD vertex pair as myImages[source]→[target], myShapesSD, myOrigins.
    pub(super) fn fill_images_vertices(&self) {
        for &(va, vb) in self.ds.shape_sd.sd_vertices_iter() {
            if va >= vb { continue; }
            let src = va;
            let sd  = vb;
            self.my_images.borrow_mut().entry(self.brep_sr(src)).or_default().push(self.brep_sr(sd));
            self.my_shapes_sd.borrow_mut().insert(self.brep_sr(src), self.brep_sr(sd));
            self.my_origins.borrow_mut().entry(self.brep_sr(sd)).or_default().push(self.brep_sr(src));
        }
    }

    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-125).
    /// Maps source edges → split images via pave-block new_edge.
    /// Also handles CommonBlocks via myShapesSD.
    pub(super) fn fill_images_edges(&self) {
        let e_base = self.ds.vertices.len();
        for (ei, edge) in self.ds.edges.iter().enumerate() {
            if edge.pave_blocks.is_empty() { continue; }
            let aE = self.brep_sr(e_base + ei);
            for pb in &edge.pave_blocks {
                let nSpR = self.ds.real_pave_block_edge(ei, pb)
                    .or(pb.0.read().unwrap().new_edge)
                    .unwrap_or(ei);
                let aSpR = self.brep_sr(e_base + nSpR);
                self.my_images.borrow_mut().entry(aE).or_default().push(aSpR);
                self.my_origins.borrow_mut().entry(aSpR).or_default().push(aE);
                if pb.0.read().unwrap().common_block_idx.is_some() {
                    if let Some(nSp) = pb.0.read().unwrap().new_edge {
                        let aSp = self.brep_sr(e_base + nSp);
                        self.my_shapes_sd.borrow_mut().insert(aSp, aSpR);
                    }
                }
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(WIRE) (Builder_1.cxx L221-276).
    /// Rebuilds wire edge lists with split sub-edges when any edge has images.
    pub(super) fn fill_images_container_wire(&self, _result: &ResultBuilder) {
        let e_base = self.ds.vertices.len();
        let mut pending: Vec<(rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>)> = Vec::new();
        {
            let my_images = self.my_images.borrow();
            for (wi, wire) in self.ds.wires.iter().enumerate() {
                let w_ref = self.brep_sr(
                    e_base + self.ds.edges.len() + wi);
                let mut a_c_im: Vec<rcad_kernel::topods::ShapeRef> = Vec::new();
                let mut has_images = false;
                for &ei in &wire.edges {
                    let e_ref = self.brep_sr(e_base + ei);
                    if let Some(imgs) = my_images.get(&e_ref) {
                        has_images = true;
                        for &img_sr in imgs {
                            if !a_c_im.contains(&img_sr) {
                                a_c_im.push(img_sr);
                            }
                        }
                    } else {
                        if !a_c_im.contains(&e_ref) {
                            a_c_im.push(e_ref);
                        }
                    }
                }
                if has_images {
                    pending.push((w_ref, a_c_im));
                }
            }
        }
        for (w_ref, a_c_im) in pending {
            self.my_images.borrow_mut().entry(w_ref).or_default().extend(a_c_im);
        }
    }
    /// ✅ OCCT-aligned: FillImagesFaces (Builder_2.cxx L215-229).
    /// 3-step dispatcher: BuildSplitFaces -> FillSameDomainFaces -> FillInternalVertices.
    pub(super) fn fill_images_faces(&self) {
        let mut result = crate::builder::result_builder::ResultBuilder::new();
        let mut t = self.my_shape.borrow_mut();
        self.build_split_faces(&mut result, &mut *t);
        if self.has_errors { return; }
        self.fill_same_domain_faces(&mut result);
        if self.has_errors { return; }
        self.fill_internal_vertices(&mut result);
    }

    /// ✅ OCCT-aligned: PostTreat (BOPAlgo_Builder.cxx L456-481).
    /// Corrects tolerances of the result shape after building.
    pub(super) fn post_treat(&mut self) {
        
        let _a_ma: std::collections::HashSet<usize> = if self.my_non_destructive {
            (0..self.ds.nb_source_shapes)
                .filter(|&i| {
                    if i >= self.ds.shape_info.len() { return false; }
                    let si = &self.ds.shape_info[i];
                    si.shape_type == rcad_kernel::topods::ShapeType::Vertex
                        || si.shape_type == rcad_kernel::topods::ShapeType::Edge
                        || si.shape_type == rcad_kernel::topods::ShapeType::Face
                })
                .collect()
        } else {
            std::collections::HashSet::new()
        };
        
        let e_base = self.ds.vertices.len();
        let mut edge_updates: Vec<(usize, f64)> = Vec::new();
        let mut vert_updates: Vec<(usize, f64)> = Vec::new();
        {
            let t = self.my_shape.borrow();
            for (ti, ts) in t.tshapes.iter().enumerate() {
                if let topods::TShape::Edge(_ed) = &**ts {
                    let ei = ti.saturating_sub(e_base);
                    if ei < self.ds.edges.len() {
                        edge_updates.push((ti, self.ds.edges[ei].geom_tol.max(0.05)));
                    }
                } else if let topods::TShape::Vertex(_vd) = &**ts {
                    if ti < self.ds.vertices.len() {
                        vert_updates.push((ti, self.ds.vertices[ti].geom_tol.max(0.05)));
                    }
                }
            }
        }
        {
            let mut t = self.my_shape.borrow_mut();
            for (ti, tol) in &edge_updates {
                if let topods::TShape::Edge(ed) = &*t.tshapes[*ti].clone() {
                    t.tshapes[*ti] = std::sync::Arc::new(topods::TShape::Edge(topods::TEdgeData {
                        tolerance: *tol,
                        ..ed.clone()
                    }));
                }
            }
            for (ti, tol) in &vert_updates {
                if let topods::TShape::Vertex(vd) = &*t.tshapes[*ti].clone() {
                    t.tshapes[*ti] = std::sync::Arc::new(topods::TShape::Vertex(topods::TVertexData {
                        tolerance: *tol,
                        ..vd.clone()
                    }));
                }
            }
        }

        
        // rcad: CorrectShapeTolerances is a BRep-level operation not yet translated.
    }

    /// --OCCT-aligned: BuildSplitFaces (Builder_2.cxx L233-374).
    ///   Iterates source faces — splits each along intersection curves.
    ///   For faces with IN/SC PBs: full BuilderFace::Perform (builder_face_perform).
    ///   For ON-only faces: BuildDraftFace.
    ///   Faces with no interferences — skipped (no images).
    pub(super) fn build_split_faces(
        &self,
        result: &mut ResultBuilder,
        t: &mut topods::BRep,
    ) {
        let a_nb_s = self.ds.nb_source_shapes;
        let mut face_counter = 0usize;
        let brep_snapshot = self.brep.borrow().clone().unwrap_or_default();
        let (brep_owned, face_refs_owned, ic_edge_map_owned) = brep_snapshot;
        let mut a_vbf: Vec<crate::builder::BuilderFace> = Vec::new();
        let mut a_vbf_face_srs: Vec<topods::ShapeRef> = Vec::new();
        // OCCT aFacesIm: draft face results keyed by source face ref.
        let mut a_faces_im_draft: std::collections::HashMap<topods::ShapeRef, Vec<topods::ShapeRef>> =
            std::collections::HashMap::new();

        for i in 0..a_nb_s {
            if i >= self.ds.shape_info.len() { continue; }
            let si = &self.ds.shape_info[i];
            if si.shape_type != rcad_kernel::topods::ShapeType::Face { continue; }
            let fi = face_counter;
            face_counter += 1;
            if fi >= self.ds.faces.len() { continue; }
            let is_a = self.ds.faces[fi].origin == ShapeOrigin::ShapeA;

            let has_info = self.ds.faces[fi].face_info.has_any_interference();
            let has_pb_in = !self.ds.faces[fi].face_info.pave_blocks_in.is_empty();
            let has_pb_sc = !self.ds.faces[fi].face_info.pave_blocks_sc.is_empty();
            let has_pb_on = !self.ds.faces[fi].face_info.pave_blocks_on.is_empty();
            let a_nb_av = self.ds.faces[fi].face_info.vertices_on.len();

            // OCCT L275-279 + L293-296 combined skip: no face info or no PBs and no AV.
            if !has_pb_in && !has_pb_sc && !has_pb_on && a_nb_av == 0 && !has_info {
                continue;
            }

            let sf_idx = self.ds.faces[fi].source_face_idx;
            let f_base = self.ds.vertices.len() + self.ds.edges.len();
            let side_offset = if is_a { 0usize } else { self.ds.a_face_count };
            let f_sr = self.brep_sr(f_base + side_offset + sf_idx);

            // OCCT L298-351: No IN/SC PBs branch.
            if !has_pb_in && !has_pb_sc {
                // OCCT L309: hasInternals, initially false.
                let mut has_internals = false;
                // OCCT L310-334: check internals and modified wires when no alone vertices.
                if a_nb_av == 0 {
                    // OCCT L321: internal edges check (first edge of wire with INTERNAL orientation).
                    for &bei in &self.ds.faces[fi].boundary_edges {
                        if bei < self.ds.edges.len() && self.ds.edges[bei].is_internal {
                            has_internals = true;
                            break;
                        }
                    }
                    // OCCT L327: modified wires check (myImages.IsBound(wire)).
                    let mut has_modified = false;
                    if !has_internals {
                        let wi_base = self.ds.vertices.len() + self.ds.edges.len();
                        for &wi in std::iter::once(&self.ds.faces[fi].outer_wire_idx)
                            .flatten()
                            .chain(self.ds.faces[fi].inner_wire_idxs.iter())
                        {
                            let w_ref = self.brep_sr(wi_base + wi);
                            if self.my_images.borrow().contains_key(&w_ref) {
                                has_modified = true;
                                break;
                            }
                        }
                    }
                    // OCCT L330-333: no internals and no modified → skip this face.
                    if !has_internals && !has_modified {
                        continue;
                    }
                }
                // OCCT L336-350: if no internals, attempt BuildDraftFace.
                if !has_internals {
                    if let Some(draft) = self.build_draft_face(fi) {
                        let (segments, wfs, vertex_positions) = draft;
                        for wf in &wfs {
                            let origin = if is_a {
                                FaceOrigin::FromA(self.ds.faces[fi].source_face_idx)
                            } else {
                                FaceOrigin::FromB(self.ds.faces[fi].source_face_idx)
                            };
                            // OCCT L344-347: store draft face in aFacesIm and continue.
                            result.emit_wire_face(fi, wf, &segments, self.ds, false, origin,
                                &vertex_positions);
                            // Create TShape and capture ref for my_images (OCCT L344-347).
                            let mut tmp_fr = Vec::new();
                            result.emit_face_topods(t, &mut tmp_fr);
                            if let Some(&draft_ref) = tmp_fr.last() {
                                if !draft_ref.is_null() {
                                    a_faces_im_draft.entry(f_sr).or_default().push(draft_ref);
                                }
                            }
                        }
                        continue; // OCCT L347: draft face stored → skip BuilderFace.
                    }
                    // OCCT L349: BuildDraftFace returned null → fall through to BuilderFace.
                }
                // OCCT L351 (implicit): has_internals → fall through to BuilderFace.
            }

            let face_sr = self.my_face_refs.borrow().get(fi).copied().unwrap_or(topods::ShapeRef::NULL);
            let mut a_le: Vec<topods::ShapeRef> = Vec::new();
            // OCCT L353: aMFence — fence for SEAM edge dedup.
            let mut a_m_fence_local: std::collections::HashSet<u64> = std::collections::HashSet::new();
            // OCCT L387-393: surface closed state (computed once per face).
            let (is_u_closed, is_v_closed) = match &self.ds.faces[fi].surface {
                s if s.is_u_closed() && s.is_v_closed() => (true, true),
                s if s.is_u_closed() => (true, false),
                s if s.is_v_closed() => (false, true),
                _ => (false, false),
            };
            let e_base = self.ds.vertices.len();
            {
                let t_shape: &topods::BRep = &*t;
                if face_sr.index < t_shape.tshapes.len() {
                    if let topods::TShape::Face(fd) = &*t_shape.tshapes[face_sr.index] {
                        // OCCT L362-363: aExp.Init(aFF, TopAbs_EDGE).
                        // rcad: iterate edges via wire topology.
                        for &wi in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                            if wi.index >= t_shape.tshapes.len() { continue; }
                            if let topods::TShape::Wire(wd) = &*t_shape.tshapes[wi.index] {
                                for &e_sr in &wd.edges {
                                    // OCCT L367: anOriE = edge orientation in this wire.
                                    let an_ori_e = e_sr.orientation;
                                    let my_images = self.my_images.borrow();

                                    // OCCT L369-385: edge NOT in myImages.
                                    if !my_images.contains_key(&e_sr) {
                                        // OCCT L371-378: INTERNAL → add FORWARD + REVERSED.
                                        if an_ori_e == topods::Orientation::Internal {
                                            let mut fwd = e_sr;
                                            fwd.orientation = topods::Orientation::Forward;
                                            a_le.push(fwd);
                                            let mut rev = e_sr;
                                            rev.orientation = topods::Orientation::Reversed;
                                            a_le.push(rev);
                                        } else {
                                            // OCCT L379-381: normal → add as-is.
                                            a_le.push(e_sr);
                                        }
                                        continue; // OCCT L383-384.
                                    }

                                    // OCCT L386+: edge has split images.
                                    let ei = e_sr.index.saturating_sub(e_base);
                                    let b_is_degenerated = ei < self.ds.edges.len()
                                        && self.ds.is_edge_degenerated(ei);
                                    // OCCT L395-404: bIsClosed — periodic-surface closed edge.
                                    let b_is_closed = if !b_is_degenerated && (is_u_closed || is_v_closed)
                                        && ei < self.ds.edges.len()
                                    {
                                        let ed = &self.ds.edges[ei];
                                        // Check if edge is closed on this face (OCCT BRep_Tool::IsClosed).
                                        let is_edge_closed_on_face = ed.start_vertex == ed.end_vertex
                                            || (is_u_closed || is_v_closed);
                                        if is_edge_closed_on_face {
                                            // OCCT L400-403: IsEdgeIsoline.
                                            let (is_ui, is_vi) = if let Some(rep) = ed.face_reps.iter().find(|r| r.face_idx == fi) {
                                                crate::builder::wire_splitter::is_edge_isoline(&rep.pcurve, ed.t_range)
                                            } else {
                                                (false, false)
                                            };
                                            (is_u_closed && is_ui) || (is_v_closed && is_vi)
                                        } else {
                                            false
                                        }
                                    } else {
                                        false
                                    };

                                    // OCCT L408-464: iterate split image edges.
                                    if let Some(imgs) = my_images.get(&e_sr) {
                                        for &sp_sr in imgs {
                                            let mut a_sp = sp_sr;

                                            // OCCT L413-417: degenerated → set orientation, push.
                                            if b_is_degenerated {
                                                a_sp.orientation = an_ori_e;
                                                a_le.push(a_sp);
                                                continue;
                                            }

                                            // OCCT L420-426: INTERNAL → push FORWARD + REVERSED.
                                            if an_ori_e == topods::Orientation::Internal {
                                                a_sp.orientation = topods::Orientation::Forward;
                                                a_le.push(a_sp);
                                                let mut rev = sp_sr;
                                                rev.orientation = topods::Orientation::Reversed;
                                                a_le.push(rev);
                                                continue;
                                            }

                                            // OCCT L429-454: SEAM / closed edge handling.
                                            if b_is_closed {
                                                if a_m_fence_local.insert(a_sp.ptr_id) {
                                                    // OCCT L433-447: DoSplitSEAMOnFace if not closed.
                                                    let sp_ei = a_sp.index.saturating_sub(e_base);
                                                    if sp_ei < self.ds.edges.len() {
                                                        crate::boptools::do_split_seam_on_face(sp_ei, fi, self.ds);
                                                    }
                                                    // OCCT L449-452: push FORWARD + REVERSED.
                                                    a_sp.orientation = topods::Orientation::Forward;
                                                    a_le.push(a_sp);
                                                    let mut rev = sp_sr;
                                                    rev.orientation = topods::Orientation::Reversed;
                                                    a_le.push(rev);
                                                }
                                                continue;
                                            }

                                            // OCCT L457-462: normal image edge → IsSplitToReverse.
                                            a_sp.orientation = an_ori_e;
                                            let sp_ei = a_sp.index.saturating_sub(e_base);
                                            if ei < self.ds.edges.len() && sp_ei < self.ds.edges.len() {
                                                let needs_rev = crate::builder::edge_builders::is_split_to_reverse(
                                                    self.ds, sp_ei, ei);
                                                if needs_rev {
                                                    a_sp.orientation = match a_sp.orientation {
                                                        topods::Orientation::Forward => topods::Orientation::Reversed,
                                                        topods::Orientation::Reversed => topods::Orientation::Forward,
                                                        other => other,
                                                    };
                                                }
                                            }
                                            a_le.push(a_sp);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for &pb_idx in &self.ds.faces[fi].face_info.pave_blocks_in {
                if pb_idx < self.ds.pave_blocks.len() {
                    let pb_ei = self.ds.pave_blocks[pb_idx].0.read().unwrap()
                        .new_edge.unwrap_or(self.ds.pave_blocks[pb_idx].0.read().unwrap().original_edge);
                    let e_sr = self.brep_sr(self.ds.vertices.len() + pb_ei);
                    a_le.push(e_sr);
                    a_le.push(topods::ShapeRef { index: e_sr.index, orientation: topods::Orientation::Reversed, ..e_sr });
                }
            }
            for &pb_idx in &self.ds.faces[fi].face_info.pave_blocks_sc {
                if pb_idx < self.ds.pave_blocks.len() {
                    let pb_ei = self.ds.pave_blocks[pb_idx].0.read().unwrap()
                        .new_edge.unwrap_or(self.ds.pave_blocks[pb_idx].0.read().unwrap().original_edge);
                    let e_sr = self.brep_sr(self.ds.vertices.len() + pb_ei);
                    a_le.push(e_sr);
                    a_le.push(topods::ShapeRef { index: e_sr.index, orientation: topods::Orientation::Reversed, ..e_sr });
                }
            }
            // OCCT L496-500: BuildPCurveForEdgesOnPlane — speed up for planar faces.
            if !self.my_non_destructive {
                if matches!(self.ds.faces[fi].surface, rcad_kernel::geom::Surface3::Plane(_)) {
                    for &e_sr in &a_le {
                        // Skip synthetic refs with no real TShape in the BRep.
                        if e_sr.index >= t.tshapes.len() { continue; }
                        let ei = e_sr.index.saturating_sub(e_base);
                        if ei >= self.ds.edges.len() { continue; }
                        // Skip edges that already have a pcurve for this face.
                        let has_pc = match &*t.tshapes[e_sr.index] {
                            topods::TShape::Edge(ed) => ed.pcurves.contains_key(&face_sr.index),
                            _ => continue,
                        };
                        if has_pc { continue; }
                        // Project 3D curve to 2D pcurve on the plane surface.
                        if let Some(pc) = crate::geom2d_api::project_curve_to_plane(
                            &self.ds.edges[ei].curve, &self.ds.faces[fi].surface)
                        {
                            // OCCT pattern: clone TEdgeData → insert → replace Arc.
                            if let topods::TShape::Edge(ed) = &*t.tshapes[e_sr.index].clone() {
                                let mut new_ed = ed.clone();
                                new_ed.pcurves.insert(face_sr.index, (pc, ed.range[0], ed.range[1]));
                                t.tshapes[e_sr.index] = std::sync::Arc::new(topods::TShape::Edge(new_ed));
                            }
                        }
                    }
                }
            }
            let mut bf = crate::builder::BuilderFace::new(
                self.ds,
                &brep_owned,
                &face_refs_owned,
                &ic_edge_map_owned,
                &self.my_face_refs,
                fi,
                is_a,
            );
            bf.set_shapes(a_le);
            a_vbf.push(bf);
            a_vbf_face_srs.push(f_sr);
        }

        let a_nb_bf = a_vbf.len();
        for k in 0..a_nb_bf {
            a_vbf[k].perform(t);
        }

        if self.has_errors { return; }

        let mut a_faces_im: std::collections::HashMap<topods::ShapeRef, Vec<topods::ShapeRef>> =
            std::collections::HashMap::new();
        for (bf, f_sr) in a_vbf.iter().zip(a_vbf_face_srs.iter()) {
            if bf.areas().is_empty() { continue; }
            let entry = a_faces_im.entry(*f_sr).or_default();
            for &sr in bf.areas() {
                entry.push(sr);
            }
        }

        // OCCT L534-552: merge draft face results into aFacesIm before storing to myImages.
        for (src_sr, draft_refs) in &a_faces_im_draft {
            let entry = a_faces_im.entry(*src_sr).or_default();
            entry.extend(draft_refs.iter().copied());
        }

        // OCCT L534-552: store aFacesIm entries into myImages with source-face orientation.
        for (src_face_sr, a_lfr) in &a_faces_im {
            // OCCT L537-538: get source face's original orientation from the DS/BRep.
            let an_ori_f = face_refs_owned.get(src_face_sr.index)
                .map(|sr| sr.orientation)
                .unwrap_or(topods::Orientation::Forward);
            let mut my_images = self.my_images.borrow_mut();
            let p_lf_im = my_images.entry(*src_face_sr).or_insert_with(Vec::new);
            for &a_fr in a_lfr {
                let mut out_sr = a_fr;
                if an_ori_f == topods::Orientation::Reversed {
                    out_sr.orientation = topods::Orientation::Reversed;
                }
                p_lf_im.push(out_sr);
            }
        }
    }

    /// ✅ OCCT-aligned: FillInternalVertices (BOPAlgo_Builder_2.cxx L929-1008).
    ///   L937-980: For each source FACE with split images:
    ///     a) Get alone vertices (myDS->AloneVertices = VerticesIn + VerticesSc
    ///        minus endpoints of PaveBlocksIn + PaveBlocksSc).
    ///     b) For each alone vertex, create (vertex, split_face) pairs for
    ///        classification via FClass2d.
    ///   L997-1007: For pairs classified as INTERNAL — add vertex to face.
    /// rcad: alone vertices computed from FaceInfo; classification via
    ///   FClass2d::perform; results stored in face_internal_vtx.
    pub(super) fn fill_internal_vertices(&self, result: &mut ResultBuilder) {
        // OCCT L937-980: iterate source faces with split images.
        for (ds_fi, ds_face) in self.ds.faces.iter().enumerate() {
            // OCCT L952-956: skip if no split images for this source face.
            let image_rfis: Vec<usize> = result.face_origins.iter().enumerate()
                .filter(|(_, origin)| match origin {
                    FaceOrigin::FromA(sfi) =>
                        ds_face.origin == ShapeOrigin::ShapeA && ds_face.source_face_idx == *sfi,
                    FaceOrigin::FromB(sfi) =>
                        ds_face.origin == ShapeOrigin::ShapeB && ds_face.source_face_idx == *sfi,
                    _ => false,
                })
                .map(|(rfi, _)| rfi)
                .collect();
            if image_rfis.is_empty() { continue; }

            // OCCT L958-960: Get alone vertices = (VerticesIn + VerticesSc)
            // minus endpoints of (PaveBlocksIn + PaveBlocksSc).
            let fi = &ds_face.face_info;
            let mut pb_endpoints: HashSet<usize> = HashSet::new();
            for &pb_idx in fi.pave_blocks_in.iter().chain(fi.pave_blocks_sc.iter()) {
                if pb_idx < self.ds.pave_blocks.len() {
                    let (nV1, nV2) = self.ds.pave_blocks[pb_idx].0.read().unwrap().indices();
                    pb_endpoints.insert(nV1);
                    pb_endpoints.insert(nV2);
                }
            }
            let alone: Vec<usize> = fi.vertices_in.iter()
                .chain(fi.vertices_sc.iter())
                .copied()
                .filter(|vi| !pb_endpoints.contains(vi))
                .collect();
            if alone.is_empty() { continue; }

            // OCCT L963-979: classify each alone vertex against each split face.
            for &vi in &alone {
                if vi >= self.ds.vertices.len() { continue; }
                let v_pt = self.ds.vertices[vi].point;

                for &rfi in &image_rfis {
                    if rfi >= result.faces.len() { continue; }

                    // Find DS face index for classification tolerance lookup.
                    let ds_fi_for_classify = match &result.face_origins[rfi] {
                        FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                        FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                        _ => None,
                    };
                    let Some(cfi) = ds_fi_for_classify else { continue };
                    if cfi >= self.ds.faces.len() { continue; }

                    // OCCT L974-977: tolerance = MAX(tolV, tolF) + fuzzyValue.
                    let tol_v = self.ds.vertices.get(vi).map_or(crate::tolerance::TOLERANCE_ABS, |v| v.geom_tol);
                    let tol_f = self.ds.faces[cfi].geom_tol;
                    let class_tol = tol_v.max(tol_f) + self.ds.fuzzy_tol;

                    let fs = &self.ds.faces[cfi].surface;
                    if let Some(uv) = world_to_uv(fs, v_pt) {
                        let fclass = crate::inttools::fclass2d::FClass2d::new(
                            self.ds, cfi, class_tol);
                        if fclass.perform(uv, true) == crate::inttools::fclass2d::State::In {
                            if rfi < result.face_internal_vtx.len() {
                                result.face_internal_vtx[rfi].push(vi);
                            }
                        }
                    }
                }
            }
        }
    }

    /// --OCCT-aligned: FillSameDomainFaces (BOPAlgo_Builder_2.cxx L580-925).
    ///   OCCT structure:
    ///   1. L584-589: Check FF interferences --return if none.
    ///   2. L597-648: Build aFaceToParent map (source solid --face) + propagate
    ///      to split images.  Prevents merging faces from the same operand solid.
    ///   3. L659-684: Collect FF-interfering face indices into aFIVec.
    ///   4. L690-739: Build edge-set map (BOPTools_Set) + planar-face set.
    ///   5. L740+: Group by edge set, check AreFacesSameDomain, remove duplicates.
    pub(super) fn fill_same_domain_faces(&self, result: &mut ResultBuilder) {
        let nf = result.faces.len();
        if nf < 2 { return; }

        
        let has_ff = !self.ds.interf_ff.is_empty();
        if !has_ff { return; }

        
        //   solid are NOT SD merged (prevents zero-thickness interior).
        //   OCCT: iterate NbSourceShapes --filter TopAbs_SOLID --TopExp_Explorer
        //   collect sub-faces --aFaceToParent.Bind(aF, aSolid) --propagate to images.
        //   rcad: use DSFace.source_solid_idx as parent-solid identity.  Result faces
        //   with the same (operand, source_solid_idx) share a parent and are NOT merged.
        let face_parent = |fi: usize| -> Option<(bool, usize)> {
            let origin = match &result.face_origins[fi] {
                FaceOrigin::FromA(_) => ShapeOrigin::ShapeA,
                FaceOrigin::FromB(_) => ShapeOrigin::ShapeB,
                _ => return None,
            };
            let ds_fi = self.ds.faces.iter().position(|f| {
                f.origin == origin && f.source_face_idx == match &result.face_origins[fi] {
                    FaceOrigin::FromA(sfi) => *sfi,
                    FaceOrigin::FromB(sfi) => *sfi,
                    _ => unreachable!(),
                }
            })?;
            let solid_idx = self.ds.faces.get(ds_fi)?.source_solid_idx?;
            Some((origin == ShapeOrigin::ShapeA, solid_idx))
        };

        
        // rcad: build (origin, source_face_idx) set from FF interferences,
        // then filter result faces to only those matching the FF set.
        let mut ff_source_set: std::collections::HashSet<(bool, usize)> = std::collections::HashSet::new();
        for ff in &self.ds.interf_ff {
            if true {
                for &dfi in &[ff.f1, ff.f2] {
                    if let Some(df) = self.ds.faces.get(dfi) {
                        ff_source_set.insert((df.origin == ShapeOrigin::ShapeA, df.source_face_idx));
                    }
                }
            }
        }
        // OCCT aFence: skip repeated checks.  Also skip result faces whose
        // source DS face has no FF interference (not in aFIVec).
        let face_origin_pair = |fi: usize| -> (bool, usize) {
            match &result.face_origins[fi] {
                FaceOrigin::FromA(sfi) => (true, *sfi),
                FaceOrigin::FromB(sfi) => (false, *sfi),
                _ => (false, usize::MAX),
            }
        };
        let mut result_fi_filtered: Vec<usize> = (0..nf)
            .filter(|fi| ff_source_set.contains(&face_origin_pair(*fi)))
            .collect();
        if result_fi_filtered.len() < 2 { return; }

        // ── Edge-set signature per face (OCCT BOPTools_Set ──
        
        // rcad: use edge index ei directly (add_edge already deduplicates
        // by vertex pair, making ei a stable identity).  Exclude degenerate
        // edges (matching OCCT's BRep_Tool::Degenerated skip).
        let face_edge_set: std::collections::HashMap<usize, Vec<usize>> =
            result_fi_filtered.iter().map(|&fi| {
                let entry = &result.faces[fi];
                let collect_ids = |edges: &[(usize, bool)]| -> Vec<usize> {
                    edges.iter()
                        .filter(|(ei, _)| !result.deg_edge_indices.contains(ei))
                        .map(|&(ei, _)| ei)
                        .collect()
                };
                let mut ids: Vec<usize> = collect_ids(&entry.0);
                for iw_es in &entry.1 {
                    ids.extend(collect_ids(iw_es));
                }
                for iw_es in &entry.9 {
                    ids.extend(collect_ids(iw_es));
                }
                ids.sort_unstable();
                ids.dedup();
                (fi, ids)
            }).collect();

        
        let mut planars: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &fi in &result_fi_filtered {
            if matches!(result.faces[fi].4, Surface3::Plane(_)) {
                // Check boundedness: non-natural-restriction faces are bounded
                let is_bounded = result.faces[fi].5.map_or(true, |uv| {
                    uv[0].is_finite() && uv[1].is_finite()
                });
                if is_bounded {
                    planars.insert(fi);
                }
            }
        }

        // ── Group by edge-set signature ──
        let mut groups: std::collections::BTreeMap<Vec<usize>, Vec<usize>> =
            std::collections::BTreeMap::new();
        for &fi in &result_fi_filtered {
            if let Some(sig) = face_edge_set.get(&fi) {
                if sig.is_empty() { continue; }
                groups.entry(sig.clone()).or_default().push(fi);
            }
        }

        // ── AreFacesSameDomain: projection-based (OCCT BOPTools_AlgoTools.cxx L1131-1197) ──
        // OCCT: PointInFace(F1) --IsValidPointForFace(point, F2, aTol)
        //   where aTol = aTolF1 + aTolF2 + max(theFuzz, Precision::Confusion())
        //   and aTolF = max(face_tolerance, max_edge_tolerance_on_face)
        // rcad: use sample_pt from result face + projection + FClass2d.
        // Map result face index --DS face index for tolerance lookup
        let ds_face_idx = |rfi: usize| -> Option<usize> {
            match &result.face_origins[rfi] {
                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                _ => None,
            }
        };
        // OCCT aTolF = max(face_tol, max_edge_tol_on_face) per face
        let face_tol_with_edges = |dsfi: usize| -> f64 {
            let mut a_tol = self.ds.faces[dsfi].geom_tol;
            for &ei in &self.ds.faces[dsfi].boundary_edges {
                if ei < self.ds.edges.len() {
                    let e_tol = self.ds.edges[ei].geom_tol;
                    if e_tol > a_tol { a_tol = e_tol; }
                }
            }
            a_tol
        };
        let mut to_remove = vec![false; nf];
        for (_edge_set, members) in groups.iter() {
            if members.len() < 2 { continue; }
            let survivors: Vec<usize> = members.iter().filter(|&&fi| !to_remove[fi]).copied().collect();
            for i in 0..survivors.len() {
                for j in (i + 1)..survivors.len() {
                    let fi = survivors[i];
                    let fj = survivors[j];
                    if face_parent(fi) == face_parent(fj) {
                        continue;
                    }
                    
                    if planars.contains(&fi) && planars.contains(&fj) {
                        to_remove[fj] = true;
                        continue;
                    }
                    // Get interior point from result face fi
                    let pt_i = result.faces[fi].8; // sample_pt
                    let pt_j = result.faces[fj].8; // sample_pt
                    let surf_j = &result.faces[fj].4;
                    let surf_i = &result.faces[fi].4;
                    // Compute tolerance: aTolF1 + aTolF2 + fuzzy
                    let ds_i = ds_face_idx(fi);
                    let ds_j = ds_face_idx(fj);
                    let a_tol = match (ds_i, ds_j) {
                        (Some(di), Some(dj)) => {
                            face_tol_with_edges(di) + face_tol_with_edges(dj) + self.ds.fuzzy_tol
                        }
                        _ => continue,
                    };
                    // OCCT: project point from fi onto fj's surface, check distance + classification
                    let (uv_j, proj_j) = crate::extrema::closest_point_on_surface(surf_j, pt_i);
                    let dist_j = proj_j.distance(pt_i);
                    let valid_j = if dist_j <= a_tol {
                        if let Some(dj) = ds_j {
                            self.context.borrow_mut().is_point_in_on_face(self.ds, dj, uv_j)
                        } else { false }
                    } else { false };
                    // Reverse: project point from fj onto fi's surface
                    let (uv_i, proj_i) = crate::extrema::closest_point_on_surface(surf_i, pt_j);
                    let dist_i = proj_i.distance(pt_j);
                    let valid_i = if dist_i <= a_tol {
                        if let Some(di) = ds_i {
                            self.context.borrow_mut().is_point_in_on_face(self.ds, di, uv_i)
                        } else { false }
                    } else { false };
                    if valid_j && valid_i {
                        // OCCT: face with smaller DS index survives.
                        // rcad: higher-index result face is removed.
                        to_remove[fj] = true;
                    }
                }
            }
        }

        // ── Apply removals ──
        let removed = to_remove.iter().filter(|&&r| r).count();
        if removed == 0 { return; }

        for fi in 0..nf {
            if to_remove[fi] {
                result.co_face_origins.push((fi, result.face_origins[fi]));
            }
        }
        let old_faces = std::mem::take(&mut result.faces);
        let old_origins = std::mem::take(&mut result.face_origins);
        for (fi, face) in old_faces.into_iter().enumerate() {
            if !to_remove[fi] {
                result.faces.push(face);
                result.face_origins.push(old_origins[fi]);
            }
        }
    }

    /// --OCCT-aligned: FillImagesContainers (Builder.cxx L363-422).
    ///   Unified dispatch matching OCCT's FillImagesContainers(TopAbs_ShapeEnum).
    ///
    /// OCCT: single function called with WIRE, SHELL, or COMPSOLID type.
    ///   Iterates source shapes, filters by type, calls FillImagesContainer.
    ///   rcad: dispatches to type-specific implementations.
    /// --OCCT-aligned: FillImagesContainers (Builder_1.cxx L172-193).
    ///   OCCT: iterates source shapes --filters by TopAbs_ShapeEnum --    ///   FillImagesContainer for each.  rcad: dispatches to type-specific handlers.
    pub(super) fn fill_images_container(&self, shape_type: ShapeType, result: &mut ResultBuilder) {
        let mut t = self.my_shape.borrow_mut();
        match shape_type {
            ShapeType::Wire => self.fill_images_container_wire(result),
            ShapeType::Shell => self.fill_images_container_shell(result, &mut *t),
            ShapeType::CompSolid => self.fill_images_container_compsolid(result, &mut *t),
            _ => {}
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(SHELL) (Builder_1.cxx L221-276).
    /// Builds TShape::Shell from result faces per DS shell.
    pub(super) fn fill_images_container_shell(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        if self.ds.shells.is_empty() { return; }
        for ds_shell in &self.ds.shells {
            let mut shell_faces: Vec<topods::ShapeRef> = Vec::new();
            for &dsfi in &ds_shell.faces {
                if dsfi >= self.ds.faces.len() { continue; }
                let origin = self.ds.faces[dsfi].origin;
                let source_fi = self.ds.faces[dsfi].source_face_idx;
                for (rfi, fo) in result.face_origins.iter().enumerate() {
                    let matches = match (fo, origin) {
                        (FaceOrigin::FromA(s), ShapeOrigin::ShapeA) => *s == source_fi,
                        (FaceOrigin::FromB(s), ShapeOrigin::ShapeB) => *s == source_fi,
                        _ => false,
                    };
                    if matches {
                        if let Some(&sr) = self.my_face_refs.borrow().get(rfi) {
                            if !shell_faces.contains(&sr) {
                                shell_faces.push(sr);
                            }
                        }
                    }
                }
            }
            if shell_faces.is_empty() { continue; }

            // Check closure via edge valence (each edge appears exactly twice)
            let is_closed = {
                let mut ecount: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
                for &fsr in &shell_faces {
                    if let Some(tf) = t.tshapes.get(fsr.index) {
                        if let topods::TShape::Face(fd) = &**tf {
                            let count_wire_edges = |wsr: topods::ShapeRef| {
                                t.tshapes.get(wsr.index).and_then(|tw| {
                                    if let topods::TShape::Wire(wd) = &**tw { Some(&wd.edges) } else { None }
                                })
                            };
                            if let Some(edges) = count_wire_edges(fd.outer_wire) {
                                for esr in edges { *ecount.entry(esr.index).or_default() += 1; }
                            }
                            for iwsr in &fd.inner_wires {
                                if let Some(edges) = count_wire_edges(*iwsr) {
                                    for esr in edges { *ecount.entry(esr.index).or_default() += 1; }
                                }
                            }
                        }
                    }
                }
                !ecount.is_empty() && ecount.values().all(|&c| c == 2)
            };

            let shell_ref = t.add_tshell(shell_faces);
            if is_closed { t.shell_mut(shell_ref).flags |= rcad_kernel::topods::tshape_flags::CLOSED; }
            let skey = topods::ShapeRef::synthetic(usize::MAX - self.my_shells.borrow().len());
            self.my_images.borrow_mut().entry(skey).or_default().push(shell_ref);
            self.my_shells.borrow_mut().push(shell_ref);
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(COMPSOLID) (Builder_1.cxx L221-276).
    /// Builds TShape::CompSolid from result solids per DS compsolid.
    pub(super) fn fill_images_container_compsolid(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        if self.ds.comp_solids.is_empty() { return; }
        for (csi, cs) in self.ds.comp_solids.iter().enumerate() {
            let mut solid_refs: Vec<topods::ShapeRef> = Vec::new();
            for &soi in &cs.solids {
                if soi >= self.ds.solids.len() { continue; }
                for &shi in &self.ds.solids[soi].shells {
                    if shi >= self.ds.shells.len() { continue; }
                    for &dsfi in &self.ds.shells[shi].faces {
                        if dsfi >= self.ds.faces.len() { continue; }
                        let origin = self.ds.faces[dsfi].origin;
                        let sfi = self.ds.faces[dsfi].source_face_idx;
                        for (rfi, fo) in result.face_origins.iter().enumerate() {
                            let matches = match (fo, origin) {
                                (FaceOrigin::FromA(s), ShapeOrigin::ShapeA) => *s == sfi,
                                (FaceOrigin::FromB(s), ShapeOrigin::ShapeB) => *s == sfi,
                                _ => false,
                            };
                            if matches {
                                if let Some(&fsr) = self.my_face_refs.borrow().get(rfi) {
                                    for &ssr in &*self.my_solids.borrow() {
                                        if !solid_refs.contains(&ssr) {
                                            if let Some(ts) = t.tshapes.get(ssr.index) {
                                                if let topods::TShape::Solid(sd) = &**ts {
                                                    if sd.shells.iter().any(|sh_sr| {
                                                        t.tshapes.get(sh_sr.index).map_or(false, |tsh| {
                                                            if let topods::TShape::Shell(shd) = &**tsh {
                                                                shd.faces.contains(&fsr)
                                                            } else { false }
                                                        })
                                                    }) {
                                                        solid_refs.push(ssr);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if !solid_refs.is_empty() {
                let cs_ref = t.add_tcompsolid(solid_refs.clone());
                self.my_compsolid_groups.borrow_mut().push(cs_ref);
                let cskey = topods::ShapeRef::synthetic(usize::MAX - self.my_compsolid_groups.borrow_mut().len());
                self.my_images.borrow_mut().entry(cskey).or_default().push(cs_ref);
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesSolids (BOPAlgo_Builder_3.cxx L60-93).
    /// 3-step: FillIn3DParts → BuildSplitSolids → FillInternalShapes.
    pub(super) fn fill_images_solids(&self, result: &mut ResultBuilder) {
        let mut t = self.my_shape.borrow_mut();
        let has_solid = self.ds.faces.iter().any(|f| f.source_solid_idx.is_some());
        if !has_solid { return; }
        let shell_assignments = self.fill_in_3d_parts(result);
        self.build_split_solids(result, &shell_assignments, &mut *t);
        self.fill_internal_shapes(result);
    }

    /// ✅ OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-232).
    /// Classifies result faces against draft solids, returns shell assignments.
    ///
    /// ✅ OCCT-aligned: BuildDraftSolid (Builder_3.cxx L267-368).
    /// Builds a draft solid from a source solid, replacing split faces with their images.
    pub(super) fn build_draft_solid(&self, result: &ResultBuilder, side: usize)
        -> (Vec<Vec<usize>>, Vec<usize>)
    {
        let origin_side = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
        let mut draft_shells: Vec<Vec<usize>> = Vec::new();
        let mut the_lif: Vec<usize> = Vec::new();

        // Iterate sub-shapes (shells) of the solid.
        for ds_shell in &self.ds.shells {
            let belongs = ds_shell.faces.iter().any(|&dsfi|
                self.ds.faces.get(dsfi).map_or(false, |f| f.origin == origin_side));
            if !belongs { continue; }

            let mut a_sh_d: Vec<usize> = Vec::new();
            let mut i_flag = false;

            // Iterate sub-shapes (faces) of the shell.
            for &dsfi in &ds_shell.faces {
                let dsf = &self.ds.faces[dsfi];
                if dsf.origin != origin_side { continue; }

                // Check if face has split images.
                let result_faces: Vec<usize> = result.face_origins.iter().enumerate()
                    .filter(|(_, fo)| match fo {
                        FaceOrigin::FromA(sfi) => dsf.origin == ShapeOrigin::ShapeA && dsf.source_face_idx == *sfi,
                        FaceOrigin::FromB(sfi) => dsf.origin == ShapeOrigin::ShapeB && dsf.source_face_idx == *sfi,
                        _ => false,
                    })
                    .map(|(i, _)| i)
                    .collect();
                if result_faces.is_empty() { continue; }

                if result_faces.len() > 1 {
                    
                    for &a_fx in &result_faces {
                        let a_fx_dfi_opt = match &result.face_origins[a_fx] {
                            FaceOrigin::FromA(s) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *s),
                            FaceOrigin::FromB(s) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *s),
                            _ => None,
                        };
                        let is_sd = a_fx_dfi_opt.map_or(false, |fx_dfi|
                            self.ds.shape_sd.has_sd_face(dsfi, fx_dfi)
                                || self.ds.shape_sd.has_sd_face(fx_dfi, dsfi));

                        if is_sd {
                            // Same-domain image face --check reverse
                            let b_to_reverse = a_fx_dfi_opt.map_or(false, |fx_dfi|
                                crate::boptools::is_split_to_reverse(
                                    self.ds.faces[dsfi].normal, self.ds.faces[fx_dfi].normal));
                            if !b_to_reverse {
                                i_flag = true;
                                if !a_sh_d.contains(&a_fx) { a_sh_d.push(a_fx); }
                            }
                        } else {
                            // Not same-domain --use original orientation
                            i_flag = true;
                            if !a_sh_d.contains(&a_fx) { a_sh_d.push(a_fx); }
                        }
                    }
                } else {
                    // No images --add original face directly
                    let fi = result_faces[0];
                    i_flag = true;
                    if !a_sh_d.contains(&fi) { a_sh_d.push(fi); }
                }
            }

            if i_flag && !a_sh_d.is_empty() {
                // Map result face indices -> DS face indices for classify_point
                let mut ds_sh = Vec::with_capacity(a_sh_d.len());
                for &rfi in &a_sh_d {
                    let dfi_opt = match result.face_origins.get(rfi) {
                        Some(FaceOrigin::FromA(sfi)) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                        Some(FaceOrigin::FromB(sfi)) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                        _ => None,
                    };
                    if let Some(dfi) = dfi_opt {
                        if !ds_sh.contains(&dfi) {
                            ds_sh.push(dfi);
                        }
                    }
                }
                if !ds_sh.is_empty() {
                    draft_shells.push(ds_sh);
                }
            }
        }

        (draft_shells, the_lif)
    }

    /// ✅ OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-263).
    ///   Phase 1: Collect all result faces (aLFaces).
    ///   Phase 2: Build draft solids from each source solid (BuildDraftSolid).
    ///   Phase 3: ClassifyFaces against draft solids.
    ///   Phase 4: Analyze results — store in myInParts + return shell assignments.
    pub(super) fn fill_in_3d_parts(&self, result: &mut ResultBuilder) -> Vec<(usize, usize, &'static str)> {
        // === Phase 1: Collect all faces ===
        let mut a_l_faces: Vec<usize> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<usize> =
            std::collections::HashSet::new();

        for (fi, fo) in result.face_origins.iter().enumerate() {
            let is_face = match fo {
                FaceOrigin::FromA(_) | FaceOrigin::FromB(_) => true,
                _ => false,
            };
            if !is_face { continue; }
            if a_m_fence.insert(fi) {
                a_l_faces.push(fi);
            }
        }

        // === Phase 2: Build draft solids ===
        let mut a_l_solids: Vec<Vec<Vec<usize>>> = Vec::new();
        let mut a_solids_if: Vec<Vec<usize>> = Vec::new();
        let mut draft_solid_origin: Vec<(usize, usize)> = Vec::new();

        for side in 0..2 {
            let (draft_shells, the_lif) = self.build_draft_solid(result, side);
            if draft_shells.is_empty() { continue; }
            a_l_solids.push(draft_shells);
            a_solids_if.push(the_lif);
            // Find the DS shell(s) matching this side (OCCT: iterate DS ShapeInfo SOLID).
            let origin_side = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
            for (si, ds_shell) in self.ds.shells.iter().enumerate() {
                if ds_shell.faces.iter().any(|&dfi|
                    self.ds.faces.get(dfi).map_or(false, |f| f.origin == origin_side))
                {
                    draft_solid_origin.push((si, side));
                    break;
                }
            }
        }

        // === Phase 3: ClassifyFaces ===
        let face_samples: Vec<DVec3> = a_l_faces.iter()
            .map(|&fi| if fi < result.faces.len() { result.faces[fi].8 } else { DVec3::ZERO })
            .collect();
        let aabb_of_face: Vec<Aabb> = a_l_faces.iter().map(|&fi| {
            // Build minimal AABB from face boundary vertices via DS
            if fi < result.face_origins.len() {
                let dfi_opt = match &result.face_origins[fi] {
                    FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                        f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                    FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                        f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                    _ => None,
                };
                if let Some(dfi) = dfi_opt {
                    let mut aabb = Aabb::empty();
                    for &vi in &self.ds.faces[dfi].boundary_verts {
                        if vi < self.ds.vertices.len() {
                            aabb.expand_point(self.ds.vertices[vi].point);
                        }
                    }
                    aabb
                } else { Aabb::empty() }
            } else { Aabb::empty() }
        }).collect();
        let aabb_of_solid: Vec<Aabb> = a_l_solids.iter().map(|shells| {
            let mut aabb = Aabb::empty();
            for sh in shells {
                for &dfi in sh {
                    if dfi < self.ds.faces.len() {
                        for &vi in &self.ds.faces[dfi].boundary_verts {
                            if vi < self.ds.vertices.len() {
                                aabb.expand_point(self.ds.vertices[vi].point);
                            }
                        }
                    }
                }
            }
            aabb
        }).collect();
        let an_in_parts_list = crate::bopalgo::classify_faces(
            &a_l_faces, &face_samples, &a_l_solids, self.ds,
            &aabb_of_face, &aabb_of_solid,
        );
        let mut an_in_parts: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (dsi, in_faces) in an_in_parts_list.into_iter().enumerate() {
            if !in_faces.is_empty() {
                an_in_parts.insert(dsi, in_faces);
            }
        }

        // === Phase 4: Analyze classification results ===
        let mut assignments: Vec<(usize, usize, &'static str)> = Vec::new();

        for (dsi, &(si, side)) in draft_solid_origin.iter().enumerate() {
            let in_faces: Vec<usize> = an_in_parts.get(&dsi).cloned().unwrap_or_default();
            let n_in = in_faces.len();

            if n_in == 0 {
                let mut has_image = false;
                if let Some(ds_shell) = self.ds.shells.get(si) {
                    let v_base = self.ds.vertices.len();
                    for &dsfi in &ds_shell.faces {
                        if let Some(dsf) = self.ds.faces.get(dsfi) {
                            for &ei in &dsf.boundary_edges {
                                if self.my_images.borrow().contains_key(
                                    &self.brep_sr(v_base + ei))
                                {
                                    has_image = true; break;
                                }
                            }
                            if has_image { break; }
                        }
                    }
                }
                if !has_image { continue; }
            }

            let state: &'static str = if n_in > 0 { "IN" } else { "OUT" };
            assignments.push((si, side, state));

            let mut my_in_parts = self.my_in_parts.borrow_mut();
            let a_nb_int = a_solids_if.get(dsi).map_or(0, |v| v.len());
            if a_nb_int > 0 || n_in > 0 {
                let p_lin = my_in_parts.entry(side).or_default();
                for &fi in &in_faces {
                    if !p_lin.contains(&fi) {
                        p_lin.push(fi);
                    }
                }
                if let Some(lif) = a_solids_if.get(dsi) {
                    for &lif_fi in lif {
                        if !p_lin.contains(&lif_fi) {
                            p_lin.push(lif_fi);
                        }
                    }
                }
            }
        }
        assignments
    }

    /// ✅ OCCT-aligned: BuildSplitSolids (Builder_3.cxx L413-618).
    /// Build result solids from draft solids and IN faces.
    pub(super) fn build_split_solids(&self, result: &mut ResultBuilder,
                          assignments: &[(usize, usize, &'static str)],
                          t: &mut topods::BRep) {
        let my_in_parts = self.my_in_parts.borrow();
        let has_in_faces = !my_in_parts.is_empty();

        let mut a_mst: Vec<std::collections::BTreeSet<usize>> = Vec::new();
        let mut result_solids: Vec<Vec<usize>> = Vec::new();

        // Helper: result face index → DS face index
        let result_to_ds = |rfi: usize, expected_origin: ShapeOrigin| -> Option<usize> {
            let fo = result.face_origins.get(rfi)?;
            let sfi = match (expected_origin, fo) {
                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => sfi,
                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => sfi,
                _ => return None,
            };
            self.ds.faces.iter().position(|f| f.origin == expected_origin && f.source_face_idx == *sfi)
        };
        // Inverse: DS face index --result face index
        let ds_to_result = |dfi: usize| -> Option<usize> {
            let dsf = self.ds.faces.get(dfi)?;
            result.face_origins.iter().position(|fo| match (dsf.origin, fo) {
                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => dsf.source_face_idx == *sfi,
                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => dsf.source_face_idx == *sfi,
                _ => false,
            })
        };

        // === Phase 0: Non-interfered solids --aMST (OCCT L431-461) ===
        //   OCCT: iterate DS ShapeInfo for TopAbs_SOLID NOT in theDraftSolids --        //         build BOPTools_Set of faces, add to aMST.
        //   rcad: shells WITHOUT IN faces are "non-interfered" --a_mst + stored as solids.
        //   --OCCT iterates DS shape info for TopAbs_SOLID entries; rcad uses assignments.
        for &(si, side, _state) in assignments {
            
            
            
            let in_faces_this: Vec<usize> = my_in_parts.get(&side).cloned().unwrap_or_default();
            if has_in_faces && !in_faces_this.is_empty() {
                continue;
            }

            
            if let Some(ds_shell) = self.ds.shells.get(si) {
                let ds_set: std::collections::BTreeSet<usize> = ds_shell.faces.iter().copied().collect();
                if ds_set.is_empty() { continue; }
                a_mst.push(ds_set);

                
                let result_faces: Vec<usize> = ds_shell.faces.iter()
                    .flat_map(|&dsfi| {
                        let dsf = &self.ds.faces[dsfi];
                        result.face_origins.iter().enumerate()
                            .filter(|(_, fo)| match (dsf.origin, fo) {
                                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => dsf.source_face_idx == *sfi,
                                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => dsf.source_face_idx == *sfi,
                                _ => false,
                            })
                            .map(|(fi, _)| fi)
                    })
                    .collect();
                if result_faces.is_empty() { continue; }
                // OCCT-aligned: create TShape::Solid in my_images.
                // Build real face refs for the shell. Use existing real refs from
                // my_face_refs (split faces from BuilderFace). For non-split faces
                // whose ref in my_face_refs is a flat-DS-index synthetic, create a
                // real tface from DS data so that solids() can traverse the topology.
                {
                    let mut sf: Vec<topods::ShapeRef> = Vec::new();
                    for &rfi in &result_faces {
                        let sr = self.my_face_refs.borrow().get(rfi).copied().unwrap_or(topods::ShapeRef::NULL);
                        if sr.ptr_id != 0 {
                            sf.push(sr);
                        } else {
                            // Synthetic ref — build real tface from DS data.
                            let origin = &result.face_origins[rfi];
                            let dsfi = match origin {
                                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                                _ => None,
                            };
                            if let Some(dfi) = dsfi {
                                let face = &self.ds.faces[dfi];
                                let e_base = self.ds.vertices.len();
                                let mut outer_edges: Vec<topods::ShapeRef> = Vec::new();
                                for &ei in &face.boundary_edges {
                                    if ei >= self.ds.edges.len() { continue; }
                                    let e = &self.ds.edges[ei];
                                    let sv_sr = t.add_tvertex(self.ds.vertices[e.start_vertex].point);
                                    let ev_sr = t.add_tvertex(self.ds.vertices[e.end_vertex].point);
                                    let e_sr = t.add_tedge(Some(e.curve.clone()), sv_sr, ev_sr, e.t_range);
                                    outer_edges.push(e_sr);
                                }
                                if outer_edges.len() >= 3 {
                                    let ow = t.add_twire(outer_edges);
                                    let real_sr = t.add_tface(Some(face.surface.clone()), ow, vec![],
                                        Some(self.ds.vertices[face.boundary_verts[0]].point),
                                        None, vec![], face.natural_restriction);
                                    if let Some(slot) = self.my_face_refs.borrow_mut().get_mut(rfi) {
                                        *slot = real_sr;
                                    }
                                    sf.push(real_sr);
                                }
                            }
                        }
                    }
                    if !sf.is_empty() {
                        let shell_ref = t.add_tshell(sf);
                        let solid_ref = t.add_tsolid(vec![shell_ref]);
                        let so_key = topods::ShapeRef::synthetic(usize::MAX - 1 - self.my_solids.borrow().len());
                        self.my_images.borrow_mut().entry(so_key).or_default().push(solid_ref);
                        self.my_solids.borrow_mut().push(solid_ref);
                    }
                }
                let csi = result.tmp_shells.len();
                result.tmp_shells.push(result_faces);
                result_solids.push(vec![csi]);
                result.solid_side_origin.push(side);
            }
        }

        // === Phase 1a: Collect interfered solid tasks (OCCT L467-518: build aVBS) ===
        struct SolidTask {
            side: usize,
            builder_solid: crate::bopds::builder_solid::BuilderSolid,
        }
        let mut tasks: Vec<SolidTask> = Vec::new();
        for &(si, side, _state) in assignments {
            let in_faces_this: Vec<usize> = my_in_parts.get(&side).cloned().unwrap_or_default();
            if !has_in_faces || in_faces_this.is_empty() { continue; }
            let origin = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
            let other_origin = if side == 0 { ShapeOrigin::ShapeB } else { ShapeOrigin::ShapeA };

            
            let mut ds_face_set: Vec<usize> = Vec::new();
            if let Some(ds_shell) = self.ds.shells.get(si) {
                for &dsfi in &ds_shell.faces {
                    if self.ds.faces.get(dsfi).map_or(false, |f| f.origin == origin) {
                        ds_face_set.push(dsfi);
                    }
                }
            }
            
            for &rfi in &in_faces_this {
                if let Some(dfi) = result_to_ds(rfi, other_origin) {
                    ds_face_set.push(dfi);
                    ds_face_set.push(dfi);
                }
            }
            ds_face_set.sort_unstable();
            ds_face_set.dedup();
            if ds_face_set.is_empty() { continue; }

            let mut bs = crate::bopds::builder_solid::BuilderSolid::new();
            bs.set_shapes(&ds_face_set);
            tasks.push(SolidTask { side, builder_solid: bs });
        }

        // === Phase 1b: BOPTools_Parallel::Perform (OCCT L531) ===
        // OCCT runs aVBS in parallel when myRunParallel==true.
        // rcad uses rayon::par_iter_mut for the same effect.
        use rayon::prelude::*;
        let ds_ref: &crate::bopds::ds::DS = self.ds;
        tasks.par_iter_mut().for_each(|task| {
            task.builder_solid.perform(ds_ref);
        });

        // === Phase 2: Collect areas --aSolidsIm + myImages (OCCT L539-617) ===
        for task in &tasks {
            for area_ds in task.builder_solid.areas() {
                
                let ds_set: std::collections::BTreeSet<usize> = area_ds.iter().copied().collect();
                if a_mst.iter().any(|s| s == &ds_set) { continue; }
                a_mst.push(ds_set);

                
                let mut result_faces: Vec<usize> = Vec::new();
                let mut mapped: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
                for &dfi in area_ds {
                    if let Some(rfi) = ds_to_result(dfi) {
                        if mapped.insert(rfi) { result_faces.push(rfi); }
                    }
                }
                if result_faces.is_empty() { continue; }

                
                {
                    let sf: Vec<topods::ShapeRef> = result_faces.iter()
                        .filter_map(|&rfi| self.my_face_refs.borrow().get(rfi).copied())
                        .collect();
                    if !sf.is_empty() {
                        let shell_ref = t.add_tshell(sf);
                        let solid_ref = t.add_tsolid(vec![shell_ref]);
                        let so_key = topods::ShapeRef::synthetic(usize::MAX - 1 - self.my_solids.borrow().len());
                        self.my_images.borrow_mut().entry(so_key).or_default().push(solid_ref);
                        self.my_solids.borrow_mut().push(solid_ref);
                    }
                }
                let csi = result.tmp_shells.len();
                result.tmp_shells.push(result_faces);
                result_solids.push(vec![csi]);
                result.solid_side_origin.push(task.side);
            }
        }

        
        result.tmp_solids = result_solids;

        // OCCT BuilderSolid::PerformAreas (BuilderSolid.cxx L397-576): void detection.
        //   --rcad: separate post-step because BuilderSolid does not perform
        //     internal void detection during bs.perform().
        self.detect_internal_voids(result, assignments);
    }
}

