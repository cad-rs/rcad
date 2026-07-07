impl<'a> BooleanBuilder<'a> {
    // ====================================================================
    // ✅ OCCT-aligned: dimension-by-dimension pipeline (PerformInternal1)
    //   BOPAlgo_Builder.cxx L310-440
    // ====================================================================

    /// ✅ OCCT-aligned: FillImagesVertices (BOPAlgo_Builder_1.cxx L40-67).
    ///   Iterates ShapesSD → builds myImages(VERTEX) + myShapesSD + myOrigins.
    ///   OCCT L42: NCollection_DataMap<int,int>::Iterator aIt(myDS->ShapesSD())
    ///   rcad: symmetric HashSet<(usize,usize)> → process once per pair (a<b).
    fn fill_images_vertices(&self) {
        // OCCT L43-48: for (; aIt.More(); aIt.Next())
        for &(va, vb) in self.ds.shape_sd.sd_vertices_iter() {
            // rcad stores symmetric pairs; process each pair once (a < b).
            if va >= vb { continue; }
            let src = va;   // OCCT: nV = aIt.Key()
            let sd  = vb;   // OCCT: nVSD = aIt.Value()

            // OCCT L56: myImages.Bound(aV, ...)->Append(aVSD)
            self.my_images.borrow_mut().entry(self.brep_sr(src)).or_default().push(self.brep_sr(sd));
            // OCCT L58: myShapesSD.Bind(aV, aVSD)
            self.my_shapes_sd.borrow_mut().insert(self.brep_sr(src), self.brep_sr(sd));
            // OCCT L60-65: myOrigins.ChangeSeek(aVSD).Append(aV)
            self.my_origins.borrow_mut().entry(self.brep_sr(sd)).or_default().push(self.brep_sr(src));
        }
    }

    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-126).
    ///   Iterates source edges → populates myImages(EDGE) + myOrigins(EDGE).
    ///   OCCT L73: aNbS = myDS->NbSourceShapes()
    ///   OCCT L78-80: filter TopAbs_EDGE
    ///   OCCT L84-86: filter HasReference (has pave blocks)
    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-126).
    ///   Reads split edges created by MakeSplitEdges (build_split_edges in PaveFiller)
    ///   via pb.0.read().unwrap().new_edge, matching OCCT's aPBR->Edge() pattern.
    ///   Creates myImages(EDGE) and myOrigins(EDGE) mappings.
    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-125).
    ///   L75-81: iterate source shapes → filter TopAbs_EDGE
    ///   L83-87: HasReference (pave blocks exist) → skip if none
    ///   L89-90: aE = aSI.Shape(); aLPB = myDS->PaveBlocks(i)
    ///   L95:    myImages.Bound(aE, ...)
    ///   L97-119: for each pave block:
    ///     L101:   aPBR = myDS->RealPaveBlock(aPB)
    ///     L103:   nSpR = aPBR->Edge()
    ///     L104:   aSpR = myDS->Shape(nSpR)
    ///     L105:   pLS->Append(aSpR)  → myImages[source].Append(split)
    ///     L107-112: myOrigins[split].Append(source)
    ///     L114-118: IsCommonBlockOnEdge → myShapesSD.Bind(source, split)
    fn fill_images_edges(&self) {
        // OCCT L73: aNbS = myDS->NbSourceShapes()
        // OCCT L75-81: iterate source shapes, filter TopAbs_EDGE
        // OCCT L83-87: filter HasReference (has pave blocks)
        // OCCT L89-90: aE = aSI.Shape(); aLPB = myDS->PaveBlocks(i)
        // OCCT L95:    myImages.Bound(aE, ...)
        // OCCT L97-119: for each pave block:
        //   L101:   aPBR = myDS->RealPaveBlock(aPB)
        //   L103:   nSpR = aPBR->Edge()
        //   L104-105: aSpR = myDS->Shape(nSpR); pLS->Append(aSpR)
        //   L107-112: myOrigins[split].Append(source)
        //   L114-118: IsCommonBlockOnEdge → myShapesSD.Bind(aPB->Edge(), aSpR)
        //              where aPB->Edge() = split edge, aSpR = real edge
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
                // OCCT L114-118: if IsCommonBlockOnEdge → myShapesSD.Bind(aPB->Edge(), aSpR)
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
    /// For each DSWire, check if any edge has split images (my_images.Seek).
    /// If so, rebuild the wire edge list with split sub-edges and store in
    /// my_images[wire_ref]. BuildResult(WIRE) reads these to create TShape::Wire.
    fn fill_images_container_wire(&self, _result: &ResultBuilder) {
        let e_base = self.ds.vertices.len();
        // Collect wire image data from immutable borrow, then apply mutations.
        let mut pending: Vec<(rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>)> = Vec::new();
        {
            let my_images = self.my_images.borrow();
            for (wi, wire) in self.ds.wires.iter().enumerate() {
                let w_ref = self.brep_sr(
                    e_base + self.ds.edges.len() + wi);
                // OCCT L224-233: check if any sub-edge has been modified
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
                // OCCT L235-240: if no sub-edge modified — skip (no wire image needed).
                // OCCT L274-275: store new wire image in myImages.
                if has_images {
                    pending.push((w_ref, a_c_im));
                }
            }
        }
        for (w_ref, a_c_im) in pending {
            self.my_images.borrow_mut().entry(w_ref).or_default().extend(a_c_im);
        }
    }
    /// ✅ OCCT-aligned: FillImagesFaces (BOPAlgo_Builder_1.cxx L376-386).
    ///   Phase 3: splits each face via WireSplitter → classifies → emits
    ///   via emit_wire_face.  rcad equivalent: for each face with IC data,
    ///   call builder_face_perform (TopoDS-based BuilderFace::Perform), then
    ///   classify_against_solid_for_boolean + classification_keep_policy.
    /// ✅ OCCT-aligned: FillImagesFaces (BOPAlgo_Builder_2.cxx L215-229).
    ///   Equivalent to BuildSplitFaces + FillSameDomainFaces + FillInternalVertices.
    ///   OCCT L258: aNbS = myDS->NbSourceShapes()
    ///   OCCT L260-266: iterates all source shapes, filters TopAbs_FACE.
    ///   OCCT L275-279: HasFaceInfo check.
    ///   OCCT L283-287: PaveBlocksIn/On/Sc + AloneVertices.
    ///   OCCT L293-296: if no PBs and no AV → skip.
    /// ✅ OCCT-aligned: FillImagesFaces (Builder_2.cxx L215-229).
    ///   Calls BuildSplitFaces → FillSameDomainFaces → FillInternalVertices.
    /// OCCT-aligned: FillImagesFaces (Builder_2.cxx L215-229).
    ///   3-step dispatcher: BuildSplitFaces → FillSameDomainFaces → FillInternalVertices.
    fn fill_images_faces(
        &self,
        result: &mut ResultBuilder,
        a_faces: &[usize],
        b_faces: &[usize],
    ) {
        let mut t = self.my_shape.borrow_mut();
        self.build_split_faces(result, a_faces, b_faces, &mut *t);
        if self.has_errors { return; }
        self.fill_same_domain_faces(result);
        if self.has_errors { return; }
        self.fill_internal_vertices(result);
    }

    /// ✅ OCCT-aligned: BuildSplitFaces (Builder_2.cxx L233-374).
    ///   Iterates source faces → splits each along intersection curves.
    ///   For faces with IN/SC PBs: full BuilderFace::Perform (builder_face_perform).
    ///   For ON-only faces: BuildDraftFace.
    ///   Faces with no interferences → skipped (no images).
    fn build_split_faces(
        &self,
        result: &mut ResultBuilder,
        a_faces: &[usize],
        b_faces: &[usize],
        t: &mut topods::BRep,
    ) {
        // OCCT L258-266: iterate all source shapes → filter TopAbs_FACE.
        for fi in 0..self.ds.faces.len() {
            let is_a = a_faces.contains(&fi);
            if !is_a && !b_faces.contains(&fi) { continue; }

            // OCCT L275: bHasFaceInfo = myDS->HasFaceInfo(i)
            let has_info = self.ds.faces[fi].face_info.has_any_interference();

            // OCCT L283-287: PBsIn → curves_sc, PBsOn → curves_on.
            let has_pb_in = !self.ds.faces[fi].face_info.pave_blocks_in.is_empty()
                || !self.ds.faces[fi].face_info.pave_blocks_sc.is_empty();
            let has_pb_sc = !self.ds.faces[fi].face_info.curves_sc.is_empty();
            let has_pb_on = !self.ds.faces[fi].face_info.pave_blocks_on.is_empty();

            // OCCT L293-296: if (!aNbPBIn && !aNbPBOn && !aNbPBSc && !aNbAV) continue.
            if !has_pb_in && !has_pb_sc && !has_pb_on && !has_info {
                continue;
            }

            // OCCT L298-332: no IN/SC PBs → BuildDraftFace for ON PBs / alone vertices.
            // OCCT L332+:    has IN/SC PBs → full BuilderFace::Perform.
            if !has_pb_in && !has_pb_sc {
                let has_internals = self.ds.faces[fi].boundary_edges.iter().any(|&ei| {
                    self.ds.edges.get(ei).map_or(false, |e| e.is_internal)
                });
                let has_modified = self.ds.faces[fi].boundary_edges.iter().any(|&ei| {
                    let e_ref = self.brep_sr(self.ds.vertices.len() + ei);
                    self.my_images.borrow().get(&e_ref).map_or(false, |imgs| {
                        imgs.len() != 1 || imgs[0].index != self.ds.vertices.len() + ei
                    })
                });
                if !has_internals && !has_modified && !has_pb_on {
                    continue;
                }
                // OCCT L336-350: if no internals → BuildDraftFace.
                if !has_internals && has_info {
                    if let Some(draft) = self.build_draft_face(fi) {
                        let (_segments, wfs, _vp) = draft;
                        for wf in &wfs {
                            let origin = if is_a {
                                FaceOrigin::FromA(self.ds.faces[fi].source_face_idx)
                            } else {
                                FaceOrigin::FromB(self.ds.faces[fi].source_face_idx)
                            };
                            result.emit_wire_face(fi, wf, &[], self.ds, false, origin,
                                &std::collections::HashMap::new());
                        }
                    }
                }
                continue;
            }

            // Has IN or SC pave blocks → full BuilderFace::Perform (TopoDS path).
            // OCCT-aligned: record source face split in myImages.
            let sf_idx = self.ds.faces[fi].source_face_idx;
            let f_base = self.ds.vertices.len() + self.ds.edges.len();
            let side_offset = if is_a { 0usize } else { self.ds.a_face_count };
            self.my_images.borrow_mut()
                .entry(self.brep_sr(f_base + side_offset + sf_idx))
                .or_insert_with(Vec::new);
            // Architecture A1: pass t so split faces create TShapes incrementally.
            // OCCT-style: BuilderFace::Perform with self-contained struct.
            if let Some(ref brep_data) = *self.brep.borrow() {
                let bf = crate::builder::BuilderFace::new(
                    self.ds,
                    &brep_data.0,
                    &brep_data.1,
                    &brep_data.2,
                    &self.my_face_refs,
                    fi,
                    is_a,
                );
                bf.perform(result, t);
            }
        }
    }

    /// ✅ OCCT-aligned: FillInternalVertices (Builder_2.cxx L929-1008).
    ///   Settle alone vertices into split faces as INTERNAL sub-shapes.
    ///
    /// OCCT flow:
    ///   L937-980: For each source FACE with split images:
    ///     a) Get alone vertices (myDS->AloneVertices → vertices ON face, not on any edge)
    ///     b) For each alone vertex, create (vertex, split_face) pairs for classification
    ///   L982-991: Classify each pair via BOPAlgo_VFI (IntTools_FClass2d)
    ///   L997-1007: For pairs classified as INTERNAL → BRep_Builder.Add(aF, aV)
    ///
    /// rcad: alone vertices = FaceInfo.vertices_on.  For each result face,
    ///   classify alone vertices from its source DS face.  If the vertex
    ///   falls inside the result face's UV boundary → add to face_internal_vtx.
    fn fill_internal_vertices(&self, result: &mut ResultBuilder) {
        // OCCT L935: BOPAlgo_VectorOfVFI aVVFI — build vertex-face pairs.
        // OCCT L937-944: iterate source shapes, filter TopAbs_FACE.
        for (ds_fi, ds_face) in self.ds.faces.iter().enumerate() {
            // OCCT L941-944: skip non-face shapes (DS only has faces here).

            // OCCT L951-956: find images (split result faces) for this source face.
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

            // OCCT L959-960: AloneVertices(i, aLIAV).
            //   Alone vertices = (VerticesIn + VerticesSc) minus endpoints of
            //   (PaveBlocksIn + PaveBlocksSc), matching BOPDS_DS.cxx L1028-1062.
            let fi = &ds_face.face_info;
            let mut pb_endpoints: HashSet<usize> = HashSet::new();
            for &pb_idx in fi.pave_blocks_in.iter().chain(fi.pave_blocks_sc.iter()) {
                if pb_idx < self.ds.pave_blocks.len() {
                    let (nV1, nV2) = self.ds.pave_blocks[pb_idx].indices();
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

            // OCCT L964-978: for each alone vertex × each image face → classify.
            for &vi in &alone {
                if vi >= self.ds.vertices.len() { continue; }
                let v_pt = self.ds.vertices[vi].point;

                for &rfi in &image_rfis {
                    if rfi >= result.faces.len() { continue; }

                    // OCCT L972: classify against split face aFIm.
                    let ds_fi_for_classify = match &result.face_origins[rfi] {
                        FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                        FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                        _ => None,
                    };
                    let Some(cfi) = ds_fi_for_classify else { continue };
                    if cfi >= self.ds.faces.len() { continue; }

                    let fs = &self.ds.faces[cfi].surface;
                    if let Some(uv) = world_to_uv(fs, v_pt) {
                        let fclass = crate::inttools::fclass2d::FClass2d::new(
                            self.ds, cfi, crate::tolerance::TOLERANCE_ABS * 100.0);
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
    ///   OCCT structure:
    ///   1. L584-589: Check FF interferences → return if none.
    ///   2. L597-648: Build aFaceToParent map (source solid → face) + propagate
    ///      to split images.  Prevents merging faces from the same operand solid.
    ///   3. L659-684: Collect FF-interfering face indices into aFIVec.
    ///   4. L690-739: Build edge-set map (BOPTools_Set) + planar-face set.
    ///   5. L740+: Group by edge set, check AreFacesSameDomain, remove duplicates.
    fn fill_same_domain_faces(&self, result: &mut ResultBuilder) {
        let nf = result.faces.len();
        if nf < 2 { return; }

        // OCCT L584-589: Check FF interferences — if none, nothing to merge.
        let has_ff = !self.ds.interf_ff.is_empty();
        if !has_ff { return; }

        // OCCT L597-648: Build aFaceToParent map — faces from the same parent
        //   solid are NOT SD merged (prevents zero-thickness interior).
        //   OCCT: iterate NbSourceShapes → filter TopAbs_SOLID → TopExp_Explorer
        //   collect sub-faces → aFaceToParent.Bind(aF, aSolid) → propagate to images.
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

        // OCCT L659-684: Collect FF-interfering DS face indices into aFIVec.
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
        // OCCT L689-741: BOPTools_Set uses TopoDS_Edge identity.
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

        // OCCT L694: aMFPlanar — track bounded planar faces for fast-path SD
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
        // OCCT: PointInFace(F1) → IsValidPointForFace(point, F2, aTol)
        //   where aTol = aTolF1 + aTolF2 + max(theFuzz, Precision::Confusion())
        //   and aTolF = max(face_tolerance, max_edge_tolerance_on_face)
        // rcad: use sample_pt from result face + projection + FClass2d.
        // Map result face index → DS face index for tolerance lookup
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
                    // OCCT L780-784: bounded planar faces with same edge set → SD fast path
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

    /// ✅ OCCT-aligned: FillImagesContainers (Builder.cxx L363-422).
    ///   Unified dispatch matching OCCT's FillImagesContainers(TopAbs_ShapeEnum).
    ///
    /// OCCT: single function called with WIRE, SHELL, or COMPSOLID type.
    ///   Iterates source shapes, filters by type, calls FillImagesContainer.
    ///   rcad: dispatches to type-specific implementations.
    /// ✅ OCCT-aligned: FillImagesContainers (Builder_1.cxx L172-193).
    ///   OCCT: iterates source shapes → filters by TopAbs_ShapeEnum →
    ///   FillImagesContainer for each.  rcad: dispatches to type-specific handlers.
    fn fill_images_container(&self, shape_type: ShapeType, result: &mut ResultBuilder) {
        let mut t = self.my_shape.borrow_mut();
        match shape_type {
            ShapeType::Wire => self.fill_images_container_wire(result),
            ShapeType::Shell => self.fill_images_container_shell(result, &mut *t),
            ShapeType::CompSolid => self.fill_images_container_compsolid(result, &mut *t),
            _ => {}
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(SHELL) (Builder_1.cxx L221-276).
    ///   L224-240: check if any sub-shape (FACE) has been modified via myImages.Seek.
    ///   L242-275: build new container (TShape::Shell) from sub-shape images.
    ///   L274: aCIm.Closed(BRep_Tool::IsClosed(aCIm)).
    ///   L275: myImages.Bound(theS).Append(aCIm) — store shell in myImages.
    fn fill_images_container_shell(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        // No DS shells → no SHELL container to process
        if self.ds.shells.is_empty() { return; }
        // OCCT L247-249: aIt.Initialize(theS) — re-iterate sub-shapes to build container.
        //   rcad: iterate DS shells. For each, find result face ShapeRefs.
        for ds_shell in &self.ds.shells {
            // Collect TShape::Face refs from self.my_face_refs.borrow() for this shell's DS faces
            let mut shell_faces: Vec<topods::ShapeRef> = Vec::new();
            for &dsfi in &ds_shell.faces {
                if dsfi >= self.ds.faces.len() { continue; }
                let origin = self.ds.faces[dsfi].origin;
                let source_fi = self.ds.faces[dsfi].source_face_idx;
                for (rfi, fo) in result.face_origins.iter().enumerate() {
                    let matches = match (fo, origin) {
                        (FaceOrigin::FromA(s), crate::bopds::ds::ShapeOrigin::ShapeA) => *s == source_fi,
                        (FaceOrigin::FromB(s), crate::bopds::ds::ShapeOrigin::ShapeB) => *s == source_fi,
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

            // OCCT L274: aCIm.Closed(BRep_Tool::IsClosed(aCIm))
            //   Check closure via edge valence (each edge appears exactly twice).
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

            // OCCT L274-275: create TShape::Shell → store in myImages
            let shell_ref = t.add_tshell(shell_faces);
            if is_closed { t.shell_mut(shell_ref).flags |= rcad_kernel::topods::tshape_flags::CLOSED; }
            let skey = topods::ShapeRef::synthetic(usize::MAX - self.my_shells.borrow().len());
            self.my_images.borrow_mut().entry(skey).or_default().push(shell_ref);
            self.my_shells.borrow_mut().push(shell_ref);
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(COMPSOLID) (Builder_1.cxx L221-276).
    ///   L224-233: iterate sub-shapes (SOLIDs), check if any has been modified.
    ///   L235-240: if none modified → early return.
    ///   L242-275: build new container from sub-shape images.
    ///
    /// rcad: iterate DS faces → find those from CompSolids → for each unique
    ///   source compsolid, check if any sub-solid has split images.  If modified,
    ///   group result solids by their source compsolid → result.tmp_compsolid_groups.
    /// ✅ OCCT-aligned: FillImagesContainer(COMPSOLID) (Builder_1.cxx L221-276).
    ///   L224-233: iterate sub-shapes (SOLIDs), check if any has been modified.
    ///   L235-240: if none modified → early return.
    ///   L242-275: build new container from sub-shape images.
    fn fill_images_container_compsolid(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        if self.ds.comp_solids.is_empty() { return; }

        // OCCT L224-233: iterate sub-shapes (SOLIDs), check myImages.IsBound.
        //   rcad: a sub-solid is modified if it has >1 result solid.
        for (csi, cs) in self.ds.comp_solids.iter().enumerate() {
            let mut solid_refs: Vec<topods::ShapeRef> = Vec::new();
            for &soi in &cs.solids {
                if soi >= self.ds.solids.len() { continue; }
                // Collect TShape::Solid refs for this DS solid's result faces
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
                                    // Find the result solid that contains this face
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
            // OCCT L242-275: create TShape::CompSolid from sub-solid images
            if !solid_refs.is_empty() {
                let cs_ref = t.add_tcompsolid(solid_refs.clone());
                self.my_compsolid_groups.borrow_mut().push(cs_ref);
                // OCCT L275: myImages.Bound(theS).Append(aCIm)
                let cskey = topods::ShapeRef::synthetic(usize::MAX - self.my_compsolid_groups.borrow_mut().len());
                self.my_images.borrow_mut().entry(cskey).or_default().push(cs_ref);
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesSolids (BOPAlgo_Builder_3.cxx L60-93).
    ///   Phase 6: group shells into solids.
    ///
    /// OCCT flow:
    ///   L60-73: check if any source shape is TopAbs_SOLID → skip if none.
    ///   L77-83: FillIn3DParts — build draft solids from each source SOLID,
    ///           classify all result faces IN/OUT of each draft solid.
    ///   L86:   BuildSplitSolids — group (draft_solid, IN/OUT) into result solids.
    ///   L92:   FillInternalShapes — add internal sub-shapes.
    ///
    /// rcad: reads source face indices from DS internally (OCCT does not pass
    ///   A/B lists as parameters — FillIn3DParts iterates myDS->ShapeInfo()).
    ///   OCCT L60-73 check: rcad's CheckData (L320-325) has already ensured
    ///   both operands have faces, so the source-solid skip never triggers.
    fn fill_images_solids(&self, result: &mut ResultBuilder) {
        let mut t = self.my_shape.borrow_mut();
        let has_solid = self.ds.faces.iter().any(|f| f.source_solid_idx.is_some());
        if !has_solid { return; }
        let shell_assignments = self.fill_in_3d_parts(result);
        self.build_split_solids(result, &shell_assignments, &mut *t);
        self.fill_internal_shapes(result);
    }

    /// ✅ OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-232).
    ///   Classify each result face against the other source solid,
    ///   store IN faces in myInParts per source solid.
    ///
    /// OCCT L107-150: collect all result faces (images + originals)
    /// OCCT L164-195: for each source SOLID, build draft solid
    /// OCCT L201-204: ClassifyFaces against all draft solids → anInParts
    /// OCCT L215-232: for each source solid with IN faces,
    ///                store in myInParts[solid] = IN_faces + INTERNAL_faces
    ///
    /// ✅ OCCT-aligned: BuildDraftSolid (Builder_3.cxx L267-368).
    ///   Build a draft solid from a source solid, replacing split faces with their
    ///   image sub-faces and collecting INTERNAL faces into theLIF.
    ///   OCCT L283-367: iterate source solid sub-shapes (shells→faces), myImages.Seek
    ///   for each face → replace with images if bound, add INTERNAL faces to theLIF.
    ///   rcad: iterates DS shells filtered by source side, finds matching result faces.
    fn build_draft_solid(&self, result: &ResultBuilder, side: usize)
        -> (Vec<Vec<usize>>, Vec<usize>)
    {
        // OCCT L280-281: aOrSd = theSolid.Orientation(); theDraftSolid.Orientation(aOrSd).
        //   rcad: solid orientation tracked per-face via FaceOrigin.
        let origin_side = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
        let mut draft_shells: Vec<Vec<usize>> = Vec::new();
        let mut the_lif: Vec<usize> = Vec::new();

        // OCCT L283-367: iterate sub-shapes (shells) of the solid.
        for ds_shell in &self.ds.shells {
            let belongs = ds_shell.faces.iter().any(|&dsfi|
                self.ds.faces.get(dsfi).map_or(false, |f| f.origin == origin_side));
            if !belongs { continue; }

            // OCCT L292-295: MakeShell(aShD); aShD.Orientation(aOrSh); iFlag = 0.
            let mut a_sh_d: Vec<usize> = Vec::new();
            let mut i_flag = false;

            // OCCT L297-360: iterate sub-shapes (faces) of the shell.
            for &dsfi in &ds_shell.faces {
                let dsf = &self.ds.faces[dsfi];
                if dsf.origin != origin_side { continue; }

                // OCCT L301: aOrF = aF.Orientation() — rcad: all DS faces are FORWARD.
                //   INTERNAL orientation is not tracked in rcad's DS, so all faces
                //   are treated as non-INTERNAL (the common case).

                // OCCT L303: if (myImages.IsBound(aF)) — check if face has split images.
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
                    // OCCT L305-346: face has images → iterate image faces
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
                            // OCCT L311-330: same-domain image face — IsSplitToReverse check
                            //   rcad: approximate with normal comparison
                            let b_to_reverse = a_fx_dfi_opt.map_or(false, |fx_dfi|
                                crate::boptools::is_split_to_reverse(
                                    self.ds.faces[dsfi].normal, self.ds.faces[fx_dfi].normal));
                            if !b_to_reverse {
                                i_flag = true;
                                if !a_sh_d.contains(&a_fx) { a_sh_d.push(a_fx); }
                            }
                            // OCCT L321-326: if bToReverse → aFx.Reverse(); then add to shell
                            //   rcad: reversed normal means the face goes to shell either way
                        } else {
                            // OCCT L333-344: not same-domain → use original orientation
                            i_flag = true;
                            if !a_sh_d.contains(&a_fx) { a_sh_d.push(a_fx); }
                        }
                    }
                } else {
                    // OCCT L348-359: no images → add original face directly
                    let fi = result_faces[0];
                    i_flag = true;
                    if !a_sh_d.contains(&fi) { a_sh_d.push(fi); }
                }
            }

            // OCCT L362-366: if (iFlag) { aShD.Closed(...); aBB.Add(theDraftSolid, aShD); }
            if i_flag && !a_sh_d.is_empty() {
                draft_shells.push(a_sh_d);
            }
        }

        (draft_shells, the_lif)
    }

    /// ✅ OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-263).
    ///   Phase 1: collect all result faces (aLFaces).
    ///   Phase 2: build draft solids from each source solid (BuildDraftSolid).
    ///   Phase 3: classify faces against each draft solid (per-face classify_point
    ///            approximates OCCT's BVH-based BOPAlgo_Tools::ClassifyFaces).
    ///   Phase 4: analyze results → store in myInParts + return assignments.
    fn fill_in_3d_parts(&self, result: &mut ResultBuilder) -> Vec<(usize, usize, &'static str)> {
        // OCCT L101: Message_ProgressScope — rcad: skipped.
        // OCCT L103: NCollection_IncAllocator — rcad: Rust allocator.

        // === Phase 1: Collect all faces (OCCT L107-150) ===
        // OCCT L107-108: aShapeBoxMap — bounding boxes for shape acceleration.
        // OCCT L111: aMFence — fence map to prevent duplicate face entries.
        // OCCT L114: aLFaces — list of all faces to classify.
        let mut a_l_faces: Vec<usize> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<usize> =
            std::collections::HashSet::new();

        // OCCT L116-150: Iterate all source FACE shapes via DS ShapeInfo.
        //   rcad: iterate result.face_origins (all result faces already resolved).
        for (fi, fo) in result.face_origins.iter().enumerate() {
            let is_face = match fo {
                FaceOrigin::FromA(_) | FaceOrigin::FromB(_) => true,
                _ => false,
            };
            if !is_face { continue; }
            // OCCT L131-149: if myImages bound → add images (with fence); else add original.
            if a_m_fence.insert(fi) {
                a_l_faces.push(fi);
            }
        }

        // === Phase 2: Build draft solids (OCCT L152-195) ===
        // OCCT L152: BRep_Builder aBB;
        // OCCT L155: aLSolids — list of draft solids for classification.
        // OCCT L157-158: aSolidsIF — internal faces per draft solid.
        // OCCT L160-162: aDraftSolid — map: source solid → draft solid.
        //   rcad: each draft solid = Vec of shell groups of DS face indices.
        let mut a_l_solids: Vec<Vec<Vec<usize>>> = Vec::new();
        let mut a_solids_if: Vec<Vec<usize>> = Vec::new();
        // (shell_idx, side) for each draft solid (replaces OCCT's source→draft map).
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

        // === Phase 3: ClassifyFaces (OCCT L197-208) ===
        // OCCT L197-199: LOCAL anInParts — classification result map: draft solid → IN faces.
        // OCCT L201-208: BOPAlgo_Tools::ClassifyFaces(aLFaces, aLSolids,...) batch BVH.
        //   rcad: using bopalgo::classify_faces with per-face classify_point.
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

        // === Phase 4: Analyze classification results (OCCT L210-262) ===
        let mut assignments: Vec<(usize, usize, &'static str)> = Vec::new();

        // OCCT L211: aNbSol = aDraftSolid.Extent()
        for (dsi, &(si, side)) in draft_solid_origin.iter().enumerate() {
            // OCCT L220: aLInFaces = IN faces for this draft solid (from anInParts).
            let in_faces: Vec<usize> = an_in_parts.get(&dsi).cloned().unwrap_or_default();
            let n_in = in_faces.len();

            // OCCT L225-238: if no IN faces, check if shell has images → skip if none.
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

            // OCCT L241: theDraftSolids.Bind(aSolid, aSDraft)
            let state: &'static str = if n_in > 0 { "IN" } else { "OUT" };
            assignments.push((si, side, state));

            // OCCT L243-261: myInParts[source] = IN_faces + INTERNAL_faces
            let mut my_in_parts = self.my_in_parts.borrow_mut();
            let a_nb_int = a_solids_if.get(dsi).map_or(0, |v| v.len());
            if a_nb_int > 0 || n_in > 0 {
                let p_lin = my_in_parts.entry(side).or_default();
                // OCCT L250-254: append IN faces
                for &fi in &in_faces {
                    if !p_lin.contains(&fi) {
                        p_lin.push(fi);
                    }
                }
                // OCCT L256-260: append INTERNAL faces (aLInternal)
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
    ///
    ///   Build result solids from draft solids and IN faces.
    ///
    ///   Phase 0 (L431-461):  Non-interfered solids → aMST (face-set dedup).
    ///   Phase 1 (L467-518):  Interfered solids → BOPAlgo_SplitSolid → collect areas.
    ///   Phase 2 (L531-537):  Parallel execution (rcad: sequential).
    ///   Phase 3 (L539-577):  Collect results + merge alerts.
    ///   Phase 4 (L580-617):  Dedup via aMST, store in myImages / myOrigins / myShapesSD.
    ///
    ///   rcad: results stored in result.tmp_solids (BuildRC applies boolean filtering).
    ///   ⏳ myImages / myOrigins / myShapesSD storage deferred to BuildRC / build_topods.
    fn build_split_solids(&self, result: &mut ResultBuilder,
                          assignments: &[(usize, usize, &'static str)],
                          t: &mut topods::BRep) {
        // OCCT L413-415: void BuildSplitSolids(theDraftSolids, theRange)
        //   rcad: assignments + saved_shells + my_in_parts replace theDraftSolids + myInParts.
        // OCCT L417-428: local variables (aAlr0, aSFS, aLSEmpty, aMFence, aMST, aVBS)
        let my_in_parts = self.my_in_parts.borrow();
        let has_in_faces = !my_in_parts.is_empty();

        // OCCT L425: aSFS — list of all faces for building new solid
        // OCCT L426: aMFence — fence to avoid processing same solid twice
        //   rcad: implicit in assignments iteration + in_faces_this filter.
        // OCCT L427: aMST — BOPTools_Set for same-domain detection (dedup).
        //   rcad: BTreeSet<usize> of DS face indices per registered set.
        let mut a_mst: Vec<std::collections::BTreeSet<usize>> = Vec::new();

        // OCCT L463-466: aSolidsIm — indexed map: source solid → list of result solids.
        //   rcad: result_solids accumulates all shells → tmp_solids at end.
        let mut result_solids: Vec<Vec<usize>> = Vec::new();

        // Helper: result face index → DS face index
        let result_to_ds = |rfi: usize, expected_origin: ShapeOrigin| -> Option<usize> {
            let fo = result.face_origins.get(rfi)?;
            let sfi = match (expected_origin, fo) {
                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => *sfi,
                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => *sfi,
                _ => return None,
            };
            self.ds.faces.iter().position(|f| f.origin == expected_origin && f.source_face_idx == sfi)
        };
        // Inverse: DS face index → result face index
        let ds_to_result = |dfi: usize| -> Option<usize> {
            let dsf = self.ds.faces.get(dfi)?;
            result.face_origins.iter().position(|fo| match (dsf.origin, fo) {
                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => dsf.source_face_idx == *sfi,
                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => dsf.source_face_idx == *sfi,
                _ => false,
            })
        };

        // === Phase 0: Non-interfered solids → aMST (OCCT L431-461) ===
        //   OCCT: iterate DS ShapeInfo for TopAbs_SOLID NOT in theDraftSolids →
        //         build BOPTools_Set of faces, add to aMST.
        //   rcad: shells WITHOUT IN faces are "non-interfered" → a_mst + stored as solids.
        //   ⏳ OCCT iterates DS shape info for TopAbs_SOLID entries; rcad uses assignments.
        for &(si, side, _state) in assignments {
            // OCCT L437-440: if (aSI.ShapeType() != TopAbs_SOLID) continue;
            // OCCT L447: if (!aMFence.Add(aS)) continue; — fence dedup.
            // OCCT L451-454: if (theDraftSolids.IsBound(aS)) continue; — skip interfered.
            let in_faces_this: Vec<usize> = my_in_parts.get(&side).cloned().unwrap_or_default();
            if has_in_faces && !in_faces_this.is_empty() {
                continue;
            }

            // OCCT L456-459: BOPTools_Set aST; aST.Add(aS, TopAbs_FACE); aMST.Add(aST);
            if let Some(ds_shell) = self.ds.shells.get(si) {
                let ds_set: std::collections::BTreeSet<usize> = ds_shell.faces.iter().copied().collect();
                if ds_set.is_empty() { continue; }
                a_mst.push(ds_set);

                // OCCT L487-488: aSolidsIm.Add(aS).Append(aSD) — store non-interfered draft solid.
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
                // OCCT-aligned: create TShape::Solid in my_images
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

            // OCCT L491-499: 1.1 Fill Shell Faces Set
            let mut ds_face_set: Vec<usize> = Vec::new();
            if let Some(ds_shell) = self.ds.shells.get(si) {
                for &dsfi in &ds_shell.faces {
                    if self.ds.faces.get(dsfi).map_or(false, |f| f.origin == origin) {
                        ds_face_set.push(dsfi);
                    }
                }
            }
            // OCCT L501-511: 1.2 Fill internal faces (FWD + REV orientations)
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

        // === Phase 2: Collect areas → aSolidsIm + myImages (OCCT L539-617) ===
        for task in &tasks {
            for area_ds in task.builder_solid.areas() {
                // OCCT L590-602: BOPTools_Set dedup via aMST.Contains / aMST.Added.
                let ds_set: std::collections::BTreeSet<usize> = area_ds.iter().copied().collect();
                if a_mst.iter().any(|s| s == &ds_set) { continue; }
                a_mst.push(ds_set);

                // OCCT L590-602: map DS faces to result faces.
                let mut result_faces: Vec<usize> = Vec::new();
                let mut mapped: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
                for &dfi in area_ds {
                    if let Some(rfi) = ds_to_result(dfi) {
                        if mapped.insert(rfi) { result_faces.push(rfi); }
                    }
                }
                if result_faces.is_empty() { continue; }

                // OCCT L603-614: store in myImages + myOrigins + myShapesSD.
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

        // OCCT L580-617: aMST-based dedup already applied per-area above.
        result.tmp_solids = result_solids;

        // OCCT BuilderSolid::PerformAreas (BuilderSolid.cxx L397-576): void detection.
        //   ⏳ rcad: separate post-step because BuilderSolid does not perform
        //     internal void detection during bs.perform().
        self.detect_internal_voids(result, assignments);
    }
}
