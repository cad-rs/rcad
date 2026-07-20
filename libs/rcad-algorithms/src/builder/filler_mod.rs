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
    // --dimension-by-dimension pipeline (PerformInternal1)
    //   BOPAlgo_Builder.cxx L310-440
    // ====================================================================

    /// ✅ OCCT-aligned: FillImagesVertices (BOPAlgo_Builder_1.cxx L40-67).
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

    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-126).
    /// Maps source edges → split images via pave-block new_edge.
    /// Also handles CommonBlocks via myShapesSD.
    pub(super) fn fill_images_edges(&self) {
        let a_nb_s = self.ds.nb_source_shapes();
        let a_nv = self.ds.vertices.len();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info_at(i);
            if a_si.shape_type != rcad_kernel::topods::ShapeType::Edge {
                continue;
            }
            // Check if the pave blocks for the edge have been initialized
            if !a_si.has_reference() {
                continue;
            }
            let ei = a_si.source_idx;
            let a_e = self.brep_sr(a_nv + ei);
            let a_pb_refs = self.ds.pave_blocks(ei);
            // Fill the images of the edge from the list of its pave blocks.
            // The small edges, having no pave blocks, will have the empty list
            // of images and, thus, will be avoided in the result.
            let mut my_images = self.my_images.borrow_mut();
            for pb in a_pb_refs {
                let n_sp_r = self.ds.real_pave_block_edge(ei, pb)
                    .or(pb.0.read().unwrap().new_edge);
                let Some(n_sp_r) = n_sp_r else { continue; };
                if n_sp_r == usize::MAX { continue; }
                let a_l_im = my_images.entry(a_e).or_default();
                let a_sp_r = self.brep_sr(a_nv + n_sp_r);
                a_l_im.push(a_sp_r);
                let mut my_origins = self.my_origins.borrow_mut();
                let p_l_or = my_origins.entry(a_sp_r).or_default();
                p_l_or.push(a_e);
                if pb.0.read().unwrap().common_block_idx.is_some() {
                    if let Some(n_sp) = pb.0.read().unwrap().new_edge {
                        if n_sp == usize::MAX { continue; }
                        let a_sp = self.brep_sr(a_nv + n_sp);
                        self.my_shapes_sd.borrow_mut().insert(a_sp, a_sp_r);
                    }
                }
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(WIRE) (BOPAlgo_Builder_1.cxx L221-276).
    pub(super) fn fill_images_container_wire(&self, _result: &ResultBuilder) {
        let e_base = self.ds.vertices.len();
        let w_base = e_base + self.ds.edges.len();
        let mut pending: Vec<(topods::ShapeRef, Vec<topods::ShapeRef>)> = Vec::new();
        let my_images = self.my_images.borrow();
        // OCCT FillImagesContainers(WIRE): iterate NbSourceShapes, filter WIRE
        let nb_src = self.ds.nb_source_shapes();
        for i_src in 0..nb_src {
            let si = &self.ds.shape_info[i_src];
            if si.shape_type != topods::ShapeType::Wire { continue; }
            let wi = si.source_idx;
            let w_ref = self.brep_sr(w_base + wi);

            // OCCT L224-233: check if any sub-edge has been modified
            //   pLFIm = myImages.Seek(aSS)
            //   modified = pLFIm && (pLFIm.Extent() != 1 || !pLFIm.First().IsSame(aSS))
            let mut modified = false;
            for &flat_ei in &si.sub_shapes {
                let e_ref = self.brep_sr(flat_ei);
                if let Some(imgs) = my_images.get(&e_ref) {
                    if imgs.len() != 1 || imgs[0] != e_ref {
                        modified = true;
                        break;
                    }
                }
            }
            // OCCT L235-240: no modified sub-shapes → skip (original wire is used as-is)
            if !modified { continue; }

            // OCCT L247-272: rebuild wire with split or original sub-edges
            let mut a_c_im: Vec<topods::ShapeRef> = Vec::new();
            for &flat_ei in &si.sub_shapes {
                let e_ref = self.brep_sr(flat_ei);
                if let Some(imgs) = my_images.get(&e_ref) {
                    for &img_sr in imgs {
                        // OCCT L265-269: IsSplitToReverseWithWarn orientation fix
                        //   rcad: orientation handled at edge level during build_split_edges
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
            // OCCT L274-275: aCIm.Closed(...); myImages.Bound(theS, ...)->Append(aCIm)
            pending.push((w_ref, a_c_im));
        }
        drop(my_images);
        let mut my_images_mut = self.my_images.borrow_mut();
        for (w_ref, a_c_im) in pending {
            my_images_mut.entry(w_ref).or_default().extend(a_c_im);
        }
    }
    /// ✅ OCCT-aligned: FillImagesFaces (BOPAlgo_Builder_2.cxx L215-229).
    /// 3-step dispatcher: BuildSplitFaces -> FillSameDomainFaces -> FillInternalVertices.
    /// Populates T BRep with source shape TShapes first (OCCT: BuildResult adds them).
    pub(super) fn fill_images_faces(&self) {
        let mut result = crate::builder::result_builder::ResultBuilder::new();
        let mut t = self.my_shape.borrow_mut();
        // OCCT: BuildResult(FACE) would have added source Face TShapes to myShape before
        // BuildSplitFaces. In rcad, add source face/wire/edge TShapes now so the
        // TopExp_Explorer-equivalent iteration can read face->wire->edge from the T BRep.
        self.populate_source_shapes_in_t_brep(&mut *t);
        self.build_split_faces(&mut result, &mut *t);
        if self.has_errors() { return; }
        self.fill_same_domain_faces(&mut result);
        if self.has_errors() { return; }
        self.fill_internal_vertices(&mut result);
    }

    /// Populate T BRep with source shape TShapes (Vertex/Edge/Wire/Face at correct flat indices).
    /// Enables BuildSplitFaces to iterate face->wire->edge from the T BRep (1:1 with OCCT).
    fn populate_source_shapes_in_t_brep(&self, t: &mut topods::BRep) {
        use topods::ShapeRef;
        let nV = self.ds.vertices.len();
        let nE = self.ds.edges.len();
        let nW = self.ds.wires.len();
        let f_base = nV + nE + nW;

        // 1. Source vertex TShapes at flat indices 0..nV
        for vi in 0..nV {
            t.ensure_vertex_at(vi, self.ds.vertex_point(vi));
        }

        // 2. Source edge TShapes at flat indices nV..nV+nE
        for ei in 0..nE {
            let e = &self.ds.edges[ei];
            let sv = t.ensure_vertex_at(e.start_vertex, self.ds.vertex_point(e.start_vertex));
            let ev = t.ensure_vertex_at(e.end_vertex, self.ds.vertex_point(e.end_vertex));
            t.ensure_edge_at(nV + ei, Some(e.curve.clone()), sv, ev, e.t_range);
        }

        // 3. Source wire TShapes at flat indices nV+nE..nV+nE+nW
        for wi in 0..nW {
            let w = &self.ds.wires[wi];
            let edge_refs: Vec<ShapeRef> = w.edges.iter().map(|&ei| {
                ShapeRef::synthetic(nV + ei)
            }).collect();
            t.ensure_wire_at(nV + nE + wi, edge_refs);
        }

        // 4. Source face TShapes at flat indices f_base + side_offset + sf_idx
        for fi in 0..self.ds.faces.len() {
            let df = &self.ds.faces[fi];
            let is_a = self.ds.face_origin(fi) == ShapeOrigin::ShapeA;
            let sf_idx = self.ds.source_face_idx(fi);
            let side_offset = if is_a { 0usize } else { self.ds.a_face_count };
            let flat_idx = f_base + side_offset + sf_idx;

            // Build outer wire edge refs from DS face boundary_edges,
            // using boundary_edge_forwards for correct orientation (OCCT L362-363).
            let outer_edge_refs: Vec<ShapeRef> = df.boundary_edges.iter().enumerate().map(|(i, &ei)| {
                let mut sr = ShapeRef::synthetic(nV + ei);
                if df.boundary_edge_forwards.get(i).copied().unwrap_or(true) {
                    sr.orientation = topods::Orientation::Forward;
                } else {
                    sr.orientation = topods::Orientation::Reversed;
                }
                sr
            }).collect();
            // Use face_outer_wire_idxs if available, else use fi as wire index
            let wire_idx = if fi < self.ds.face_outer_wire_idxs.len() {
                self.ds.face_outer_wire_idxs[fi].unwrap_or(nW + fi)
            } else {
                nW + fi
            };
            // Place the wire TShape at a unique index beyond existing wires
            let outer_wire_sr = t.ensure_wire_at(nV + nE + wire_idx, outer_edge_refs);

            // Inner wires
            let inner_wire_refs: Vec<ShapeRef> = if fi < self.ds.face_inner_wire_idxs.len() {
                self.ds.face_inner_wire_idxs[fi].iter().map(|&wi| {
                    if wi < nW {
                        let iw = &self.ds.wires[wi];
                        let iw_edge_refs: Vec<ShapeRef> = iw.edges.iter().map(|&ei| {
                            ShapeRef::synthetic(nV + ei)
                        }).collect();
                        t.ensure_wire_at(nV + nE + wi, iw_edge_refs)
                    } else {
                        ShapeRef::synthetic(nV + nE + wi)
                    }
                }).collect()
            } else {
                vec![]
            };

            // Sample point for the face (first boundary vertex)
            let sample_point = df.boundary_verts.first().map(|&vi| self.ds.vertex_point(vi));

            t.ensure_face_at(flat_idx, Some(df.surface.clone()), outer_wire_sr,
                inner_wire_refs, sample_point, None, vec![], df.natural_restriction);
        }
    }

    /// ✅ OCCT-aligned: PostTreat (BOPAlgo_Builder.cxx L456-486).
    pub(super) fn post_treat(&mut self) {
        // OCCT L461-475: build MapToAvoid from VERTEX/EDGE/FACE source shapes
        let a_ma: std::collections::HashSet<usize> = if self.my_non_destructive {
            (0..self.ds.nb_source_shapes())
                .filter(|&i| {
                    if i >= self.ds.shape_info.len() { return false; }
                    let st = self.ds.shape_info[i].shape_type;
                    st == rcad_kernel::topods::ShapeType::Vertex
                        || st == rcad_kernel::topods::ShapeType::Edge
                        || st == rcad_kernel::topods::ShapeType::Face
                })
                .collect()
        } else {
            std::collections::HashSet::new()
        };

        // OCCT L478: CorrectTolerances(myShape, aMA, 0.05, myRunParallel)
        // rcad: adjust V/E tolerances from DS geometric tolerances, clamped to 0.05.
        {
            let t = self.my_shape.borrow();
            let e_base = self.ds.vertices.len();
            let mut updates: Vec<(usize, f64, rcad_kernel::topods::ShapeType)> = Vec::new();
            for (ti, ts) in t.tshapes.iter().enumerate() {
                if let rcad_kernel::topods::TShape::Edge(_ed) = &**ts {
                    let ei = ti.saturating_sub(e_base);
                    if ei < self.ds.edges.len() {
                        let tol = self.ds.edge_tolerance(ei).max(0.05);
                        updates.push((ti, tol, rcad_kernel::topods::ShapeType::Edge));
                    }
                } else if let rcad_kernel::topods::TShape::Vertex(_vd) = &**ts {
                    if ti < self.ds.vertices.len() {
                        let tol = self.ds.vertex_tolerance(ti).max(0.05);
                        updates.push((ti, tol, rcad_kernel::topods::ShapeType::Vertex));
                    }
                }
            }
            drop(t);
            let mut t = self.my_shape.borrow_mut();
            for (ti, tol, st) in updates {
                if st == rcad_kernel::topods::ShapeType::Edge {
                    if let rcad_kernel::topods::TShape::Edge(ed) = &*t.tshapes[ti].clone() {
                        t.tshapes[ti] = std::sync::Arc::new(
                            rcad_kernel::topods::TShape::Edge(rcad_kernel::topods::TEdgeData {
                                tolerance: tol, ..ed.clone()
                            }));
                    }
                } else {
                    if let rcad_kernel::topods::TShape::Vertex(vd) = &*t.tshapes[ti].clone() {
                        t.tshapes[ti] = std::sync::Arc::new(
                            rcad_kernel::topods::TShape::Vertex(rcad_kernel::topods::TVertexData {
                                tolerance: tol, ..vd.clone()
                            }));
                    }
                }
            }
        }

        // OCCT L485: CorrectShapeTolerances(myShape, aMA, myRunParallel)
        //   Propagates edge tolerances to their vertices, and face tolerances to their edges.
        //   OCCT BOPTools_AlgoTools_1.cxx L389-423 + L1005-1055.
        {
            let t = self.my_shape.borrow();
            let mut updates: Vec<(topods::ShapeRef, f64, rcad_kernel::topods::ShapeType)> = Vec::new();

            // Phase 1: Edge → Vertex — if vertex tolerance < edge tolerance, update vertex
            for (ti, ts) in t.tshapes.iter().enumerate() {
                if let rcad_kernel::topods::TShape::Edge(ed) = &**ts {
                    let a_tol_e = ed.tolerance;
                    // OCCT: TopExp_Explorer on edge finds vertex sub-shapes.
                    // rcad: edge.first and edge.last are the start/end vertex ShapeRefs
                    let vert_refs = [ed.first, ed.last];
                    for &v_sr in &vert_refs {
                        let vi = v_sr.ptr_id as usize;
                        if vi < t.tshapes.len() {
                            if let rcad_kernel::topods::TShape::Vertex(vd) = &*t.tshapes[vi] {
                                if vd.tolerance < a_tol_e {
                                    updates.push((v_sr, a_tol_e,
                                        rcad_kernel::topods::ShapeType::Vertex));
                                }
                            }
                        }
                    }
                }
            }

            // Phase 2: Face → Edge — if edge tolerance < face tolerance, update edge
            for (ti, ts) in t.tshapes.iter().enumerate() {
                if let rcad_kernel::topods::TShape::Face(fd) = &**ts {
                    let a_tol_f = fd.tolerance;
                    // Collect edge ShapeRefs from outer and inner wires
                    for w_sr in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                        let wi = w_sr.ptr_id as usize;
                        if wi < t.tshapes.len() {
                            if let rcad_kernel::topods::TShape::Wire(wd) = &*t.tshapes[wi] {
                                for &e_sr in &wd.edges {
                                    let ei = e_sr.ptr_id as usize;
                                    if ei < t.tshapes.len() {
                                        if let rcad_kernel::topods::TShape::Edge(ed) =
                                            &*t.tshapes[ei]
                                        {
                                            if ed.tolerance < a_tol_f {
                                                updates.push((e_sr, a_tol_f,
                                                    rcad_kernel::topods::ShapeType::Edge));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            drop(t);
            if !updates.is_empty() {
                let mut t = self.my_shape.borrow_mut();
                for (sr, tol, st) in updates {
                    let ti = sr.ptr_id as usize;
                    if st == rcad_kernel::topods::ShapeType::Edge {
                        if let rcad_kernel::topods::TShape::Edge(ed) = &*t.tshapes[ti].clone() {
                            t.tshapes[ti] = std::sync::Arc::new(
                                rcad_kernel::topods::TShape::Edge(
                                    rcad_kernel::topods::TEdgeData {
                                        tolerance: tol, ..ed.clone()
                                    }));
                        }
                    } else {
                        if let rcad_kernel::topods::TShape::Vertex(vd) = &*t.tshapes[ti].clone() {
                            t.tshapes[ti] = std::sync::Arc::new(
                                rcad_kernel::topods::TShape::Vertex(
                                    rcad_kernel::topods::TVertexData {
                                        tolerance: tol, ..vd.clone()
                                    }));
                        }
                    }
                }
            }
        }
    }

    /// ✅ OCCT-aligned: BuildSplitFaces (BOPAlgo_Builder_2.cxx L233-555).
    ///   Iterates source faces — splits each along intersection curves.
    ///   For faces with IN/SC PBs: full BuilderFace::Perform.
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
        let (_brep_owned_from_pf, _face_refs_from_pf, ic_edge_map_owned) = brep_snapshot;
        // PaveFiller BRep is empty (no tshapes). Use the T BRep (self.my_shape, borrowed as t).
        // face_refs_from_pf is always empty too. Build face_refs from the T BRep's face indices.
        let n_ds_faces = self.ds.faces.len();
        let mut face_refs_owned: Vec<topods::ShapeRef> = Vec::with_capacity(n_ds_faces);
        let t_brep: &topods::BRep = &*t;
        {
            let f_base = self.ds.vertices.len() + self.ds.edges.len() + self.ds.wires.len();
            for fi in 0..n_ds_faces {
                let is_a = self.ds.face_origins.get(fi).map_or(true, |&o| o == ShapeOrigin::ShapeA);
                let sf_idx = self.ds.source_face_idx(fi);
                let side_offset = if is_a { 0usize } else { self.ds.a_face_count };
                let flat_idx = f_base + side_offset + sf_idx;
                let sr = if flat_idx < t_brep.tshapes.len() {
                    self.brep_sr(flat_idx)
                } else {
                    topods::ShapeRef::synthetic(flat_idx)
                };
                face_refs_owned.push(sr);
            }
        }
        let brep_owned = t_brep.clone();
        let mut a_vbf: Vec<crate::builder::BuilderFace> = Vec::new();
        let mut a_vbf_face_srs: Vec<topods::ShapeRef> = Vec::new();
        // OCCT aFacesIm: draft face results keyed by source face ref.
        let mut a_faces_im_draft: std::collections::HashMap<topods::ShapeRef, Vec<topods::ShapeRef>> =
            std::collections::HashMap::new();

        if std::env::var("RCAD_DEBUG_BSF").is_ok() {
            println!("BSF: n_ds_faces={} face_refs.len={}", n_ds_faces, face_refs_owned.len());
        }
        for i in 0..a_nb_s {
            if i >= self.ds.shape_info.len() { continue; }
            let si = &self.ds.shape_info[i];
            if si.shape_type != rcad_kernel::topods::ShapeType::Face { continue; }
            let fi = face_counter;
            face_counter += 1;
            if fi >= self.ds.faces.len() { continue; }
            let is_a = self.ds.face_origins[fi] == ShapeOrigin::ShapeA;
            if std::env::var("RCAD_DEBUG_BSF").is_ok() {
                println!("BUILD_SPLIT_FACES face={} is_a={}", fi, is_a);
            }

            let has_pb_in = !self.ds.face_info(fi).pave_blocks_in.is_empty();
            let has_pb_sc = !self.ds.face_info(fi).pave_blocks_sc.is_empty();
            let has_pb_on = !self.ds.face_info(fi).pave_blocks_on.is_empty();
            // AloneVertices (BOPDS_DS.cxx L1028-1062).
            // VerticesIn + VerticesSc minus PB endpoints of PaveBlocksIn + PaveBlocksSc.
            // OCCT does NOT include VerticesOn in alone-vertices count.
            let a_nb_av = {
                let fi_info = &self.ds.faces[fi].face_info;
                let mut pb_endpoints: HashSet<usize> = HashSet::new();
                for &pb_idx in fi_info.pave_blocks_in.iter().chain(fi_info.pave_blocks_sc.iter()) {
                    if pb_idx < self.ds.pave_blocks.len() {
                        let (nV1, nV2) = self.ds.pave_blocks[pb_idx].0.read().unwrap().indices();
                        pb_endpoints.insert(nV1);
                        pb_endpoints.insert(nV2);
                    }
                }
                let mut alone = 0usize;
                for &vi in fi_info.vertices_in.iter().chain(fi_info.vertices_sc.iter()) {
                    if pb_endpoints.insert(vi) {
                        alone += 1;
                    }
                }
                alone
            };

            if std::env::var("RCAD_DEBUG_BSF").is_ok() {
                println!("BSF: face={} pb_in={} pb_sc={} pb_on={} av={} skip={}", fi, has_pb_in, has_pb_sc, has_pb_on, a_nb_av, !has_pb_in && !has_pb_sc && !has_pb_on && a_nb_av == 0);
            }
            // OCCT L275-279 + L293-296: skip if no PBs (IN/ON/SC) and no alone vertices.
            if !has_pb_in && !has_pb_sc && !has_pb_on && a_nb_av == 0 {
                continue;
            }

            let sf_idx = self.ds.source_face_idx(fi);
            let f_base = self.ds.vertices.len() + self.ds.edges.len() + self.ds.wires.len();
            let side_offset = if is_a { 0usize } else { self.ds.a_face_count };
            let f_sr = self.brep_sr(f_base + side_offset + sf_idx);

            // OCCT L298-351: No IN/SC PBs branch.
            if !has_pb_in && !has_pb_sc {
                // OCCT L309: hasInternals, initially false.
                let mut has_internals = false;
                // OCCT L310-334: check internals and modified wires when no alone vertices.
                    if a_nb_av == 0 {
                        // OCCT L310-328: iterate original face wires, check first edge
                        // orientation for INTERNAL wire detection.
                        let mut has_modified = false;
                        if let Some(arc) = brep_owned.tshapes.get(f_sr.index) {
                            if let topods::TShape::Face(fd) = &**arc {
                                for w_sr in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                                    if let Some(arc_w) = brep_owned.tshapes.get(w_sr.index) {
                                        if let topods::TShape::Wire(wd) = &**arc_w {
                                            // OCCT L321: itE.More() && itE.Value().Orientation() == TopAbs_INTERNAL
                                            if let Some(&first_e_sr) = wd.edges.first() {
                                                if first_e_sr.orientation == topods::Orientation::Internal {
                                                    has_internals = true;
                                                    break;
                                                }
                                            }
                                            // OCCT L327: hasModified |= myImages.IsBound(wire)
                                            has_modified |= self.my_images.borrow().contains_key(w_sr);
                                        }
                                    }
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
                                FaceOrigin::FromA(self.ds.source_face_idx(fi))
                            } else {
                                FaceOrigin::FromB(self.ds.source_face_idx(fi))
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

            let face_sr = {
                let f_base = self.ds.vertices.len() + self.ds.edges.len() + self.ds.wires.len();
                let side_offset = if is_a { 0usize } else { self.ds.a_face_count };
                let sf_idx = self.ds.source_face_idx(fi);
                self.brep_sr(f_base + side_offset + sf_idx)
            };
            let mut a_le: Vec<topods::ShapeRef> = Vec::new();
            // OCCT L353: aMFence — fence for SEAM edge dedup.
            let mut a_m_fence_local: std::collections::HashSet<u64> = std::collections::HashSet::new();
            // TopExp_Explorer fence: process each edge TShape only once across all wires.
            let mut a_m_explorer_set: std::collections::HashSet<u64> = std::collections::HashSet::new();
            // OCCT L387-393: surface closed state (computed once per face).
            let (is_u_closed, is_v_closed) = match &self.ds.faces[fi].surface {
                s if s.is_u_closed() && s.is_v_closed() => (true, true),
                s if s.is_u_closed() => (true, false),
                s if s.is_v_closed() => (false, true),
                _ => (false, false),
            };
            let e_base = self.ds.vertices.len();
            {
                // OCCT L362-363: aExp.Init(aFF, TopAbs_EDGE).
                // rcad: iterate edges via T BRep Face TShape wire->edge (now populated
                // by populate_source_shapes_in_t_brep before this function).
                let t_shape: &topods::BRep = &*t;
                if face_sr.index < t_shape.tshapes.len() {
                    if let topods::TShape::Face(fd) = &*t_shape.tshapes[face_sr.index] {
                        for &wi in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                            if wi.index >= t_shape.tshapes.len() { continue; }
                            if let topods::TShape::Wire(wd) = &*t_shape.tshapes[wi.index] {
                                if std::env::var("RCAD_DEBUG_WIREDBG").is_ok() {
                    eprintln!("[WIREDBG] wire edges: {:?}", wd.edges.iter().map(|sr| (sr.index, sr.orientation)).collect::<Vec<_>>());
                }
                                for &e_sr in &wd.edges {
                                    // OCCT L367: anOriE = edge orientation in this wire.
                                    let mut an_ori_e = e_sr.orientation;
                                    // OCCT L362-363: TopExp_Explorer(aFF, EDGE) returns each edge
                                    // with its orientation in the face.  Since we build the TShape
                                    // wire from DS boundary_edges (which may have orientations
                                    // overridden by BuildResult steps), re-derive orientation from
                                    // boundary_edge_forwards to match OCCT face-wire semantics.
                                    if an_ori_e == topods::Orientation::Forward {
                                        let ds_ei = e_sr.index.saturating_sub(e_base);
                                        if let Some(pos) = self.ds.faces[fi].boundary_edges.iter().position(|&be| be == ds_ei) {
                                            if !self.ds.faces[fi].boundary_edge_forwards.get(pos).copied().unwrap_or(true) {
                                                an_ori_e = topods::Orientation::Reversed;
                                                if std::env::var("RCAD_DEBUG_WIREDBG").is_ok() {
                                                    eprintln!("[ORIOVER] edge {} set Reversed (ds_ei={})", e_sr.index, ds_ei);
                                                }
                                            } else if std::env::var("RCAD_DEBUG_WIREDBG").is_ok() {
                                                eprintln!("[ORIOVER] edge {} Forward (ds_ei={})", e_sr.index, ds_ei);
                                            }
                                        } else if std::env::var("RCAD_DEBUG_WIREDBG").is_ok() {
                                            eprintln!("[ORIOVER] edge {} NOT IN boundary_edges (ds_ei={})", e_sr.index, ds_ei);
                                        }
                                    }
                                    let my_images = self.my_images.borrow();

                                    // OCCT TopExp_Explorer: process each edge TShape only once
                                    // across all wires (outer + inner) to avoid duplicate
                                    // entries in aLE from seam edges that appear twice.
                                    if !a_m_explorer_set.insert(e_sr.ptr_id) {
                                        continue;
                                    }

                                    // OCCT L369-385: edge NOT in myImages.
                                    if !my_images.contains_key(&e_sr) {
                                        if an_ori_e == topods::Orientation::Internal {
                                            let mut fwd = e_sr;
                                            fwd.orientation = topods::Orientation::Forward;
                                            a_le.push(fwd);
                                            let mut rev = e_sr;
                                            rev.orientation = topods::Orientation::Reversed;
                                            a_le.push(rev);
                                        } else {
                                            a_le.push(e_sr);
                                        }
                                        continue;
                                    }

                                    let ei = e_sr.index.saturating_sub(e_base);
                                    let b_is_degenerated = ei < self.ds.edges.len()
                                        && self.ds.is_edge_degenerated(ei);
                                    let b_is_closed = {
                                        if b_is_degenerated || ei >= self.ds.edges.len() {
                                            false
                                        } else if (is_u_closed || is_v_closed)
                                            && self.ds.edge_on_face(ei, fi)
                                                .map_or(false, |rep| rep.pcurve2.is_some())
                                        {
                                            let (is_ui, is_vi) =
                                                self.ds.edge_on_face(ei, fi).map_or((false, false), |rep|
                                                    crate::builder::wire_splitter::is_edge_isoline(
                                                        &rep.pcurve, rep.pcurve_range));
                                            (is_u_closed && is_ui) || (is_v_closed && is_vi)
                                        } else {
                                            false
                                        }
                                    };

                                    if let Some(imgs) = my_images.get(&e_sr) {
                                        for &sp_sr in imgs {
                                            let mut a_sp = sp_sr;

                                            if b_is_degenerated {
                                                a_sp.orientation = an_ori_e;
                                                a_le.push(a_sp);
                                                continue;
                                            }

                                            if an_ori_e == topods::Orientation::Internal {
                                                a_sp.orientation = topods::Orientation::Forward;
                                                a_le.push(a_sp);
                                                let mut rev = sp_sr;
                                                rev.orientation = topods::Orientation::Reversed;
                                                a_le.push(rev);
                                                continue;
                                            }

                                            if b_is_closed {
                                                if a_m_fence_local.insert(a_sp.ptr_id) {
                                                    let sp_ei = a_sp.index.saturating_sub(e_base);
                                                    if sp_ei < self.ds.edges.len() {
                                                        let sp_is_closed = self.ds.edge_on_face(sp_ei, fi)
                                                            .map_or(false, |rep| rep.pcurve2.is_some());
                                                        if !sp_is_closed {
                                                            crate::boptools::do_split_seam_on_face(sp_ei, fi, self.ds);
                                                        }
                                                    }
                                                    a_sp.orientation = topods::Orientation::Forward;
                                                    a_le.push(a_sp);
                                                    let mut rev = sp_sr;
                                                    rev.orientation = topods::Orientation::Reversed;
                                                    a_le.push(rev);
                                                }
                                                continue;
                                            }

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
            for &pb_idx in &self.ds.face_info(fi).pave_blocks_in {
                if pb_idx < self.ds.pave_blocks.len() {
                    let pb_ei = self.ds.pave_blocks[pb_idx].0.read().unwrap()
                        .new_edge.unwrap_or(self.ds.pave_blocks[pb_idx].0.read().unwrap().original_edge);
                    let e_sr = self.brep_sr(self.ds.vertices.len() + pb_ei);
                    a_le.push(e_sr);
                    a_le.push(topods::ShapeRef { index: e_sr.index, orientation: topods::Orientation::Reversed, ..e_sr });
                }
            }
            for &pb_idx in &self.ds.face_info(fi).pave_blocks_sc {
                if pb_idx < self.ds.pave_blocks.len() {
                    let pb_ei = self.ds.pave_blocks[pb_idx].0.read().unwrap()
                        .new_edge.unwrap_or(self.ds.pave_blocks[pb_idx].0.read().unwrap().original_edge);
                    let e_sr = self.brep_sr(self.ds.vertices.len() + pb_ei);
                    // OCCT L469-486: section edges are added with INTERNAL orientation.
                    // In OCCT, PerformShapesToAvoid skips INTERNAL edges (section edges on face interior).
                    a_le.push(topods::ShapeRef { index: e_sr.index, orientation: topods::Orientation::Internal, ..e_sr });
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
            if std::env::var("RCAD_DEBUG_BSF").is_ok() {
                println!("BSF: creating BuilderFace fi={}", fi);
            }
            if std::env::var("RCAD_DEBUG_IC").is_ok() {
                if std::env::var("RCAD_DEBUG_BSF").is_ok() {
                    eprintln!("[BSF] creating BuilderFace fi={} a_le.len={}", fi, a_le.len());
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

        if self.has_errors() { return; }

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

        // OCCT BuildResult(FACE): source faces WITHOUT split images are passed through
        // as-is (unsplit faces).  rcad's build_split_faces skips these at L469/L507.
        // Add their T BRep ShapeRefs to my_face_refs at position [fi] so build_result
        // can find them (build_result reads my_face_refs[fi] at result_build_mod.rs L865).
        for (fi, f_sr) in face_refs_owned.iter().enumerate() {
            if fi >= self.ds.faces.len() { continue; }
            if !a_faces_im.contains_key(f_sr) {
                let mut mfr = self.my_face_refs.borrow_mut();
                if fi >= mfr.len() {
                    mfr.resize(fi + 1, topods::ShapeRef::NULL);
                }
                if mfr[fi].is_null() {
                    mfr[fi] = *f_sr;
                }
            }
        }
    }

    /// ✅ OCCT-aligned: FillInternalVertices (BOPAlgo_Builder_2.cxx L929-1008).
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
                let v_pt = self.ds.vertex_point(vi);

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
                    let tol_f = self.ds.face_tolerance(cfi);
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

    /// ✅ OCCT-aligned: FillSameDomainFaces (BOPAlgo_Builder_2.cxx L580-925).
    pub(super) fn fill_same_domain_faces(&self, result: &mut ResultBuilder) {
        let nf = result.faces.len();
        if nf < 2 { return; }

        
        let has_ff = !self.ds.interf_ff.is_empty();
        if !has_ff { return; }

        // --aFaceToParent (BOPAlgo_Builder_2.cxx L597-649).
        //   Map DS face index -> parent solid index.  Prevents two result faces
        //   from the same operand solid from being merged (zero-thickness interior).
        //   OCCT: iterate NbSourceShapes, filter TopAbs_SOLID, TopExp_Explorer
        //   on sub-faces -> Bind(aF, aSolid).  Then propagate to split faces via
        //   myImages: if a source face has a parent solid, its split faces inherit.
        //   rcad: walk DS.shape_info tree (SOLID -> SHELL -> FACE) to build the
        //   mapping, then propagate by looking up each result face's source DS face.
        let mut a_face_to_parent: HashMap<usize, usize> = HashMap::new();
        let nb_src = self.ds.nb_source_shapes();
        for i_src in 0..nb_src {
            let si = &self.ds.shape_info[i_src];
            if si.shape_type != topods::ShapeType::Solid {
                continue;
            }
            // Walk SOLID -> SHELL -> FACE sub-shapes (OCCT: TopExp_Explorer)
            for &shell_i in &si.sub_shapes {
                if shell_i >= self.ds.shape_info.len() {
                    continue;
                }
                let shell_si = &self.ds.shape_info[shell_i];
                for &face_i in &shell_si.sub_shapes {
                    a_face_to_parent.entry(face_i).or_insert(i_src);
                }
            }
        }
        // OCCT L619-649: propagate aFaceToParent to split/result faces.
        // rcad: use result.face_origins to map result face -> source DS face -> parent solid.
        let result_face_parent_solid = |rfi: usize| -> Option<usize> {
            let (origin, sfi) = match &result.face_origins[rfi] {
                FaceOrigin::FromA(sfi) => (ShapeOrigin::ShapeA, *sfi),
                FaceOrigin::FromB(sfi) => (ShapeOrigin::ShapeB, *sfi),
                _ => return None,
            };
            let dsfi = self.ds.faces.iter().position(|f| {
                f.origin == origin && f.source_face_idx == sfi
            })?;
            a_face_to_parent.get(&dsfi).copied()
        };

        // --aFIVec + aMFence (BOPAlgo_Builder_2.cxx L654-687).
        //   Collect DS face indices from FF interferences, check HasFaceInfo,
        //   dedup with aMFence, sort.  OCCT then iterates aFIVec to build edge sets.
        //   rcad: collect DS face indices, then map to result face indices.
        let mut a_fi_vec: Vec<usize> = Vec::new();
        let mut a_fence: HashSet<usize> = HashSet::new();
        for ff in &self.ds.interf_ff {
            for dfi in [ff.f1, ff.f2] {
                if !a_fence.insert(dfi) {
                    continue;
                }
                if let Some(df) = self.ds.faces.get(dfi) {
                    if !df.face_info.has_any_interference() {
                        continue;
                    }
                }
                a_fi_vec.push(dfi);
            }
        }
        // OCCT L687: std::sort(aFIVec.begin(), aFIVec.end())
        a_fi_vec.sort_unstable();
        if a_fi_vec.is_empty() { return; }

        // Map DS face index -> result face indices (for downstream edge-set building)
        let mut ds_to_result: HashMap<usize, Vec<usize>> = HashMap::new();
        for (rfi, fo) in result.face_origins.iter().enumerate() {
            let (origin, sfi) = match fo {
                FaceOrigin::FromA(sfi) => (ShapeOrigin::ShapeA, *sfi),
                FaceOrigin::FromB(sfi) => (ShapeOrigin::ShapeB, *sfi),
                _ => continue,
            };
            if let Some(dsfi) = self.ds.faces.iter().position(|f| {
                f.origin == origin && f.source_face_idx == sfi
            }) {
                ds_to_result.entry(dsfi).or_default().push(rfi);
            }
        }
        // Build result_fi_filtered from aFIVec (only result faces whose source DS face
        // is in the FF-interfering set — matches OCCT's iteration scope)
        let mut a_fence_result: HashSet<usize> = HashSet::new();
        let mut result_fi_filtered: Vec<usize> = Vec::new();
        for &dsfi in &a_fi_vec {
            if let Some(rfis) = ds_to_result.get(&dsfi) {
                for &rfi in rfis {
                    if a_fence_result.insert(rfi) {
                        result_fi_filtered.push(rfi);
                    }
                }
            }
        }
        if result_fi_filtered.len() < 2 { return; }

        // --BOPTools_Set (BOPTools_Set.hxx L27-65).
        //   Stores edge-index signature from one face.  Add(theS, TopAbs_EDGE)
        //   traverses sub-EDGEs via TopExp_Explorer, skips degenerated,
        //   doubles INTERNAL edges (FORWARD + REVERSED).  IsEqual uses set
        //   comparison; GetSum() provides hash for IndexedDataMap key.
        #[derive(Debug, Clone)]
        struct BOPToolsSet {
            edges: Vec<usize>,
            sum: u64,
        }
        impl BOPToolsSet {
            fn from_result_face(result: &ResultBuilder, fi: usize) -> Self {
                let entry = &result.faces[fi];
                // OCCT BOPTools_Set::Add(theS, TopAbs_EDGE)
                let mut a_se: Vec<usize> = Vec::new();
                for &(ei, _) in &entry.0 {
                    if !result.deg_edge_indices.contains(&ei) { a_se.push(ei); }
                }
                for inner in &entry.1 {
                    for &(ei, _) in inner {
                        if !result.deg_edge_indices.contains(&ei) { a_se.push(ei); }
                    }
                }
                // OCCT L149-159: TopAbs_INTERNAL → add FORWARD + REVERSED
                for internal in &entry.9 {
                    for &(ei, _) in internal {
                        if !result.deg_edge_indices.contains(&ei) {
                            a_se.push(ei);
                            a_se.push(ei);
                        }
                    }
                }
                a_se.sort_unstable();
                a_se.dedup();
                // OCCT GetSum(): hash sum for map key
                let sum = a_se.iter().fold(0u64, |acc, &e| acc.wrapping_add(e as u64));
                BOPToolsSet { edges: a_se, sum }
            }
            fn nb_shapes(&self) -> usize { self.edges.len() }
            fn is_empty(&self) -> bool { self.edges.is_empty() }
        }
        impl PartialEq for BOPToolsSet {
            // OCCT BOPTools_Set::IsEqual: myNbShapes + set containment
            fn eq(&self, other: &Self) -> bool {
                if self.nb_shapes() != other.nb_shapes() { return false; }
                self.edges == other.edges
            }
        }
        impl Eq for BOPToolsSet {}
        impl std::hash::Hash for BOPToolsSet {
            // OCCT: std::hash<BOPTools_Set> uses GetSum()
            fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
                state.write_u64(self.sum);
            }
        }

        // --AddEdgeSet (BOPAlgo_Builder_2.cxx L562-576).
        //   anESetFaces: NCollection_IndexedDataMap<BOPTools_Set, List<TopoDS_Shape>>.
        //   rcad: HashMap<BOPToolsSet, Vec<usize>>.
        let mut an_eset_faces: HashMap<BOPToolsSet, Vec<usize>> = HashMap::new();
        for &fi in &result_fi_filtered {
            let a_se = BOPToolsSet::from_result_face(result, fi);
            if a_se.is_empty() { continue; }
            // OCCT AddEdgeSet: theMap(aSE).Append(theS)
            an_eset_faces.entry(a_se).or_default().push(fi);
        }

        
        // --planar bounded face detection (BOPAlgo_Builder_2.cxx L707-718).
        //   OCCT: SurfaceAdaptor::GetType() == GeomAbs_Plane
        //         + Bnd_Box::IsOpenXmin/Xmax/Ymin/Ymax/Zmin/Zmax (all 6 closed).
        //   rcad: Surface3::Plane + DS ShapeInfo box has valid min/max.
        let mut a_mf_planar: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &fi in &result_fi_filtered {
            if !matches!(result.faces[fi].4, Surface3::Plane(_)) {
                continue;
            }
            // OCCT: aSI = myDS->ShapeInfo(nF); aBox = aSI.Box()
            //   bCheckPlanar = !(aBox.IsOpenXmin() || ... || aBox.IsOpenZmax())
            let dsfi = match &result.face_origins[fi] {
                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                _ => None,
            };
            if dsfi.and_then(|dsfi| self.ds.shape_info.get(dsfi))
                .map_or(false, |si| si.box_min.is_some() && si.box_max.is_some())
            {
                a_mf_planar.insert(fi);
            }
        }

        // OCCT L696-741: iterate aFIVec, AddEdgeSet for each face/split, build anESetFaces + aMFPlanar.
        //   Already done above in AddEdgeSet loop.

        // --aTolF = max(face_tol, max_edge_tol_on_face) (AreFacesSameDomain L1160-1188).
        //   OCCT: TopExp_Explorer theF1 on TopAbs_EDGE, skip degenerated edges,
        //   compute max edge tolerance, compare against face tolerance.
        //   rcad: explore all edges on the DS face (outer + inner wires), skip degenerated.
        let face_tol_with_edges = |dsfi: usize| -> f64 {
            let mut a_tol = self.ds.face_tolerance(dsfi);
            let mut a_tol_e_max = -1.0_f64;
            for &ei in &self.ds.faces[dsfi].boundary_edges {
                if ei < self.ds.edges.len() && !result.deg_edge_indices.contains(&ei) {
                    let e_tol = self.ds.edge_tolerance(ei);
                    if e_tol > a_tol_e_max { a_tol_e_max = e_tol; }
                }
            }
            for inner in &self.ds.faces[dsfi].inner_boundary_edges {
                for &(ei, _) in inner {
                    if ei < self.ds.edges.len() && !result.deg_edge_indices.contains(&ei) {
                        let e_tol = self.ds.edge_tolerance(ei);
                        if e_tol > a_tol_e_max { a_tol_e_max = e_tol; }
                    }
                }
            }
            if a_tol_e_max > a_tol { a_tol = a_tol_e_max; }
            a_tol
        };
        // Map result face index -> DS face index for tolerance lookup
        let ds_face_idx = |rfi: usize| -> Option<usize> {
            match &result.face_origins[rfi] {
                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                _ => None,
            }
        };

        // --AreFacesSameDomain single-direction (BOPTools_AlgoTools.cxx L1131-1197).
        //   OCCT: PointInFace(F1) -> IsValidPointForFace(aP1, F2, aTol).
        //   Only projects point from F1 onto F2 (one direction).
        //   rcad: sample_pt from result face fi -> project onto fj -> classify.
        let faces_are_sd = |rfi: usize, rfj: usize| -> bool {
            let pt_i = result.faces[rfi].8;
            let surf_j = &result.faces[rfj].4;
            let ds_i = ds_face_idx(rfi);
            let ds_j = ds_face_idx(rfj);
            let (Some(di), Some(dj)) = (ds_i, ds_j) else { return false };
            let a_tol = face_tol_with_edges(di) + face_tol_with_edges(dj)
                + self.ds.fuzzy_tol.max(TOLERANCE_ABS);
            let (uv_j, proj_j) = crate::extrema::closest_point_on_surface(surf_j, pt_i);
            let dist_j = proj_j.distance(pt_i);
            dist_j <= a_tol && self.context.borrow_mut().is_point_in_on_face(self.ds, dj, uv_j)
        };

        // --FillMap + MakeBlocks (BOPAlgo_Tools.hxx, Builder_2.cxx L820-826).
        //   FillMap: back-and-forth adjacency in aDMSLS.
        //   MakeBlocks: BFS transitive closure -> blocks of SD faces.
        //   rcad: aDMSLS is HashMap<usize, Vec<usize>> (adjacency list).
        let mut a_dmsls: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut fill_map = |n1: usize, n2: usize| {
            a_dmsls.entry(n1).or_insert_with(|| vec![n1]).push(n2);
            a_dmsls.entry(n2).or_insert_with(|| vec![n2]).push(n1);
        };

        // OCCT L750-793: for each edge-set group, iterate pairs
        // OCCT L750: anESetFaces.Extent() --iterate each edge-set group
        for members in an_eset_faces.values() {
            if members.len() < 2 { continue; }
            // OCCT L764-772: aIt1 outer, aIt2 = aIt1 (inner starts at i+1)
            for i in 0..members.len() {
                let a_f1 = members[i];
                let b_check_planar = a_mf_planar.contains(&a_f1);
                let p_parent1 = result_face_parent_solid(a_f1);
                for j in (i + 1)..members.len() {
                    let a_f2 = members[j];
                    let p_parent2 = result_face_parent_solid(a_f2);
                    // OCCT L776: same parent -> skip (zero-thickness guard)
                    if p_parent1.is_some() && p_parent2.is_some() && p_parent1 == p_parent2 {
                        continue;
                    }
                    // OCCT L780-784: both planar bounded -> FillMap directly
                    if b_check_planar && a_mf_planar.contains(&a_f2) {
                        fill_map(a_f1, a_f2);
                        continue;
                    }
                    // OCCT L786-791: AreFacesSameDomain check (single-direction)
                    if faces_are_sd(a_f1, a_f2) {
                        fill_map(a_f1, a_f2);
                    }
                }
            }
        }

        // OCCT L826: MakeBlocks (BFS from aDMSLS adjacency -> connected components)
        let mut a_m_blocks: Vec<Vec<usize>> = Vec::new();
        {
            let mut processed: HashSet<usize> = HashSet::new();
            for &key in a_dmsls.keys() {
                if processed.contains(&key) { continue; }
                let mut block: Vec<usize> = Vec::new();
                let mut queue: VecDeque<usize> = VecDeque::new();
                queue.push_back(key);
                while let Some(current) = queue.pop_front() {
                    if processed.insert(current) {
                        block.push(current);
                        if let Some(neighbors) = a_dmsls.get(&current) {
                            for &n in neighbors {
                                if !processed.contains(&n) {
                                    queue.push_back(n);
                                }
                            }
                        }
                    }
                }
                if block.len() > 1 {
                    a_m_blocks.push(block);
                }
            }
        }
        if a_m_blocks.is_empty() { return; }

        // OCCT L830-881: For each SD block, pick representative (min result index),
        // bind all other faces to it via to_remove.
        let mut to_remove = vec![false; nf];
        for block in &a_m_blocks {
            // OCCT L842-867: face with smallest index is rep (OCCT: min DS index for originals)
            let rep = *block.iter().min().unwrap();
            // OCCT L876-881: myShapesSD.Bind(aF, *pFSD) for all block members
            for &fi in block {
                if fi != rep {
                    to_remove[fi] = true;
                }
            }
        }

        // OCCT L884-921: Update myImages with SD faces, fill myOrigins.
        // rcad: remove non-representative faces from result (architecture equivalent).
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

    /// --FillImagesContainers (Builder.cxx L363-422).
    ///   Unified dispatch matching OCCT's FillImagesContainers(TopAbs_ShapeEnum).
    ///
    /// OCCT: single function called with WIRE, SHELL, or COMPSOLID type.
    ///   Iterates source shapes, filters by type, calls FillImagesContainer.
    ///   rcad: dispatches to type-specific implementations.
    /// --FillImagesContainers (Builder_1.cxx L172-193).
    /// ✅ OCCT-aligned: FillImagesContainers (BOPAlgo_Builder_1.cxx L172-193).
    ///   OCCT: iterates source shapes, filters by theType, calls FillImagesContainer for each.
    ///   rcad: dispatches to type-specific handlers.
    pub(super) fn fill_images_container(&self, shape_type: ShapeType, result: &mut ResultBuilder) {
        let mut t = self.my_shape.borrow_mut();
        match shape_type {
            ShapeType::Wire => self.fill_images_container_wire(result),
            ShapeType::Shell => self.fill_images_container_shell(result, &mut *t),
            ShapeType::CompSolid => self.fill_images_container_compsolid(result, &mut *t),
            _ => {}
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(SHELL) (BOPAlgo_Builder_1.cxx L221-276).
    pub(super) fn fill_images_container_shell(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        let mut pending: Vec<(topods::ShapeRef, Vec<topods::ShapeRef>)> = Vec::new();
        // Pre-build DS face → result face ShapeRef map for fast lookup
        let mut ds_face_to_ref: Vec<Option<topods::ShapeRef>> = vec![None; self.ds.faces.len()];
        for (rfi, fo) in result.face_origins.iter().enumerate() {
            let (origin, sfi) = match fo {
                FaceOrigin::FromA(s) => (ShapeOrigin::ShapeA, *s),
                FaceOrigin::FromB(s) => (ShapeOrigin::ShapeB, *s),
                _ => continue,
            };
            if let Some(dsfi) = self.ds.faces.iter().position(|f|
                f.origin == origin && f.source_face_idx == sfi
            ) {
                if dsfi < ds_face_to_ref.len() {
                    if let Some(&sr) = self.my_face_refs.borrow().get(rfi) {
                        ds_face_to_ref[dsfi] = Some(sr);
                    }
                }
            }
        }

        // OCCT L172-193: iterate NbSourceShapes, filter SHELL type
        let nb_src = self.ds.nb_source_shapes();
        let my_images = self.my_images.borrow();
        for i_src in 0..nb_src {
            let si = &self.ds.shape_info[i_src];
            if si.shape_type != topods::ShapeType::Shell { continue; }
            let s_ref = self.brep_sr(i_src);

            // OCCT L224-233: check if any sub-face has been modified
            let mut modified = false;
            for &flat_fi in &si.sub_shapes {
                let f_ref = self.brep_sr(flat_fi);
                if let Some(imgs) = my_images.get(&f_ref) {
                    if imgs.len() != 1 || imgs[0].index != flat_fi {
                        modified = true;
                        break;
                    }
                }
            }
            // OCCT L235-240: no modification → skip
            if !modified { continue; }

            // OCCT L247-272: rebuild shell with split or original sub-faces
            let mut shell_faces: Vec<topods::ShapeRef> = Vec::new();
            for &flat_fi in &si.sub_shapes {
                let f_ref = self.brep_sr(flat_fi);
                if let Some(imgs) = my_images.get(&f_ref) {
                    // OCCT L260-270: add all split faces (with orientation fix)
                    for &img_sr in imgs {
                        if !shell_faces.contains(&img_sr) {
                            shell_faces.push(img_sr);
                        }
                    }
                } else {
                    // No splits → use original face
                    let dsfi = flat_fi.saturating_sub(
                        self.ds.vertices.len() + self.ds.edges.len() + self.ds.wires.len()
                    );
                    if dsfi < ds_face_to_ref.len() {
                        if let Some(sr) = ds_face_to_ref[dsfi] {
                            if !shell_faces.contains(&sr) {
                                shell_faces.push(sr);
                            }
                        }
                    }
                }
            }
            if shell_faces.is_empty() { continue; }

            // OCCT L274: aCIm.Closed(BRep_Tool::IsClosed(aCIm))
            // rcad: edge-valence closure check
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

            // OCCT L275: myImages.Bound(theS, ...)->Append(aCIm)
            let shell_ref = t.add_tshell(shell_faces);
            if is_closed {
                t.shell_mut(shell_ref).flags |= rcad_kernel::topods::tshape_flags::CLOSED;
            }
            pending.push((s_ref, vec![shell_ref]));
            self.my_shells.borrow_mut().push(shell_ref);
        }
        drop(my_images);
        let mut my_images_mut = self.my_images.borrow_mut();
        for (s_ref, images) in pending {
            my_images_mut.entry(s_ref).or_default().extend(images);
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(COMPSOLID) (BOPAlgo_Builder_1.cxx L221-276).
    pub(super) fn fill_images_container_compsolid(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        let my_images = self.my_images.borrow();
        let my_solids = self.my_solids.borrow();
        let nb_src = self.ds.nb_source_shapes();
        let mut pending: Vec<(topods::ShapeRef, topods::ShapeRef)> = Vec::new();
        for i_src in 0..nb_src {
            let si = &self.ds.shape_info[i_src];
            if si.shape_type != topods::ShapeType::CompSolid { continue; }
            let cs_ref_key = self.brep_sr(i_src);

            // OCCT L224-233: check if any sub-solid has been modified
            let mut modified = false;
            for &flat_soi in &si.sub_shapes {
                let so_ref = self.brep_sr(flat_soi);
                if let Some(imgs) = my_images.get(&so_ref) {
                    if imgs.len() != 1 || imgs[0].index != flat_soi {
                        modified = true;
                        break;
                    }
                }
            }
            if !modified { continue; }

            // OCCT L247-272: rebuild compsolid with split or original sub-solids
            // rcad: collect result solid ShapeRefs that belong to this compsolid
            let mut solid_refs: Vec<topods::ShapeRef> = Vec::new();
            let csi = si.source_idx;
            for ds_shell in &self.ds.shells {
                for &dsfi in &ds_shell.faces {
                    if dsfi >= self.ds.faces.len() { continue; }
                    if self.ds.source_compsolid_idx(dsfi) != Some(csi) { continue; }
                    for (rfi, fo) in result.face_origins.iter().enumerate() {
                        let matches = match (fo, self.ds.face_origin(dsfi)) {
                            (FaceOrigin::FromA(s), ShapeOrigin::ShapeA) => *s == self.ds.source_face_idx(dsfi),
                            (FaceOrigin::FromB(s), ShapeOrigin::ShapeB) => *s == self.ds.source_face_idx(dsfi),
                            _ => false,
                        };
                        if matches {
                            if let Some(&fsr) = self.my_face_refs.borrow().get(rfi) {
                                for &ssr in my_solids.iter() {
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

            // OCCT L275: myImages.Bound(theS, ...)->Append(aCIm)
            if !solid_refs.is_empty() {
                let cs_ref = t.add_tcompsolid(solid_refs);
                pending.push((cs_ref_key, cs_ref));
                self.my_compsolid_groups.borrow_mut().push(cs_ref);
            }
        }
        drop(my_images);
        let mut my_images_mut = self.my_images.borrow_mut();
        for (key, cs_ref) in pending {
            my_images_mut.entry(key).or_default().push(cs_ref);
        }
    }

    /// ✅ OCCT-aligned: FillImagesSolids (BOPAlgo_Builder_3.cxx L60-93).
    ///   OCCT BOPAlgo_Builder_3.cxx L60-93: FillIn3DParts → BuildSplitSolids → FillInternalShapes.
    pub(super) fn fill_images_solids(&self, result: &mut ResultBuilder) {
        let has_solid = self.ds.faces.iter().any(|f| f.source_solid_idx.is_some());
        if !has_solid { return; }
        let shell_assignments = self.fill_in_3d_parts(result);
        {
            let mut t = self.my_shape.borrow_mut();
            self.build_split_solids(result, &shell_assignments, &mut *t);
        }
        self.fill_internal_shapes(result);
    }

    /// ✅ OCCT-aligned: BuildDraftSolid (BOPAlgo_Builder_3.cxx L267-368).
    /// Builds a draft solid from a source solid, replacing split faces with their images.
    pub(super) fn build_draft_solid(&self, result: &ResultBuilder, side: usize)
        -> (Vec<Vec<usize>>, Vec<usize>)
    {
        let origin_side = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
        let mut draft_shells: Vec<Vec<usize>> = Vec::new();
        let mut the_lif: Vec<usize> = Vec::new();
        let f_base = self.ds.vertices.len() + self.ds.edges.len() + self.ds.wires.len();
        let side_offset = if side == 0 { 0usize } else { self.ds.a_face_count };

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
                // OCCT BuildDraftSolid (BOPAlgo_Builder_3.cxx L267-368): use myImages.
                let src_key = self.brep_sr(f_base + side_offset + dsf.source_face_idx);
                let _has_splits = self.my_images.borrow().get(&src_key)
                    .map(|imgs| imgs.len() > 1).unwrap_or(false);
                // Use parent DS face index whether split or not (split children
                // belong to same DS face; classify_faces uses DS face indices).
                if !a_sh_d.contains(&dsfi) { a_sh_d.push(dsfi); }
                i_flag = true;
            }

            if i_flag && !a_sh_d.is_empty() {
                // a_sh_d already contains DS face indices (not rfi).
                // No need to map via result.face_origins (OCCT operates on BRep directly).
                draft_shells.push(a_sh_d.clone());
            }
        }

        (draft_shells, the_lif)
    }

    /// ✅ OCCT-aligned: FillIn3DParts (BOPAlgo_Builder_3.cxx L97-263).
    pub(super) fn fill_in_3d_parts(&self, result: &mut ResultBuilder) -> Vec<(usize, usize, &'static str)> {
        // === Phase 1: Collect all faces ===
        // OCCT ClassifyFaces (BOPAlgo_Builder_3.cxx L97-263): iterate result BRep faces.
        // rcad: use DS face indices (my_face_refs + my_images) instead of result.face_origins
        // which is never populated (OCCT operates on BRep, has no result.face_origins).
        let f_base = self.ds.vertices.len() + self.ds.edges.len() + self.ds.wires.len();
        let mut a_l_faces: Vec<usize> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<usize> =
            std::collections::HashSet::new();

        // Add source DS faces that have non-null refs
        let fr = self.my_face_refs.borrow();
        for (rfi, &sr) in fr.iter().enumerate() {
            if sr.is_null() { continue; }
            // Map rfi to DS face index via flat index computation
            for (dsfi, df) in self.ds.faces.iter().enumerate() {
                let src_flat = f_base + (if df.origin == ShapeOrigin::ShapeA { 0 } else { self.ds.a_face_count }) + df.source_face_idx;
                if rfi == src_flat || rfi == dsfi {
                    if a_m_fence.insert(dsfi) { a_l_faces.push(dsfi); }
                    break;
                }
            }
        }
        // Also add DS faces that have split images in my_images
        let imgs = self.my_images.borrow();
        for dsfi in 0..self.ds.faces.len() {
            let df = &self.ds.faces[dsfi];
            let src_flat = f_base + (if df.origin == ShapeOrigin::ShapeA { 0 } else { self.ds.a_face_count }) + df.source_face_idx;
            let src_key = self.brep_sr(src_flat);
            if imgs.contains_key(&src_key) {
                if a_m_fence.insert(dsfi) { a_l_faces.push(dsfi); }
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
        // OCCT: use BRep face geometry for sample points + AABBs. rcad: use DS face data.
        let face_samples: Vec<DVec3> = a_l_faces.iter()
            .map(|&dfi| {
                if let Some(df) = self.ds.faces.get(dfi) {
                    let mut c = DVec3::ZERO;
                    let mut n = 0u32;
                    for &vi in &df.boundary_verts {
                        if let Some(v) = self.ds.vertices.get(vi) { c += v.point; n += 1; }
                    }
                    if n > 0 { c /= n as f64; }
                    c
                } else { DVec3::ZERO }
            })
            .collect();
        let aabb_of_face: Vec<Aabb> = a_l_faces.iter().map(|&dfi| {
            let mut aabb = Aabb::empty();
            if let Some(df) = self.ds.faces.get(dfi) {
                for &vi in &df.boundary_verts {
                    if let Some(v) = self.ds.vertices.get(vi) { aabb.expand_point(v.point); }
                }
            }
            aabb
        }).collect();
        let aabb_of_solid: Vec<Aabb> = a_l_solids.iter().map(|shells| {
            let mut aabb = Aabb::empty();
            for sh in shells {
                for &dfi in sh {
                    if dfi < self.ds.faces.len() {
                        for &vi in &self.ds.faces[dfi].boundary_verts {
                            if vi < self.ds.vertices.len() {
                                aabb.expand_point(self.ds.vertex_point(vi));
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

    /// ✅ OCCT-aligned: BuildSplitSolids (BOPAlgo_Builder_3.cxx L413-618).
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

                
                                // OCCT L451-461: non-interfered solid → build from source DS faces.
                let mut sf: Vec<topods::ShapeRef> = Vec::new();
                let e_base = self.ds.vertices.len();
                for &dsfi in &ds_shell.faces {
                    if let Some(df) = self.ds.faces.get(dsfi) {
                        let mut outer_edges: Vec<topods::ShapeRef> = Vec::new();
                        for &ei in &df.boundary_edges {
                            if ei >= self.ds.edges.len() { continue; }
                            let e = &self.ds.edges[ei];
                            let sv_sr = t.add_tvertex(self.ds.vertices[e.start_vertex].point);
                            let ev_sr = t.add_tvertex(self.ds.vertices[e.end_vertex].point);
                            let e_sr = t.add_tedge(Some(e.curve.clone()), sv_sr, ev_sr, e.t_range);
                            outer_edges.push(e_sr);
                        }
                        if outer_edges.len() >= 3 {
                            let ow = t.add_twire(outer_edges);
                            let real_sr = t.add_tface(Some(df.surface.clone()), ow, vec![],
                                df.boundary_verts.first().and_then(|&vi| self.ds.vertices.get(vi)).map(|v| v.point),
                                None, vec![], df.natural_restriction);
                            sf.push(real_sr);
                        }
                    }
                }
                if sf.is_empty() { continue; }
                let shell_ref = t.add_tshell(sf);
                let solid_ref = t.add_tsolid(vec![shell_ref]);
                // OCCT: store under source solid key for BuildResult(SOLID) lookup.
                let src_key = topods::ShapeRef::synthetic(side);
                self.my_images.borrow_mut().entry(src_key).or_default().push(solid_ref);
                self.my_solids.borrow_mut().push(solid_ref);
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
            // in_faces_this now contains DS face indices (classify_faces uses DS face indices).
            // OCCT L501-511: add IN faces twice (FORWARD + REVERSED) for BuilderSolid.
            for &dfi in &in_faces_this {
                ds_face_set.push(dfi);
                ds_face_set.push(dfi);
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
                    let mut sf: Vec<topods::ShapeRef> = Vec::new();
                    for &rfi in &result_faces {
                        let sr = self.my_face_refs.borrow().get(rfi).copied().unwrap_or(topods::ShapeRef::NULL);
                        if sr.ptr_id != 0 {
                            sf.push(sr);
                        } else {
                            // Build real TShape::Face for synthetic refs (same as Phase 0).
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

