impl<'a> BooleanBuilder<'a> {
    /// OCCT-aligned: BuildRC (BOPAlgo_BOP.cxx L583-867, SOLID filtering part).
    ///   Filter result.tmp_solids by boolean operation type using args/tools face-set
    ///   comparison (BOPTools_Set):
    ///     1. Split solids by source side (solid_side_origin) into args and tools groups
    ///     2. For each args solid, build its DS face set and check if any tools solid
    ///        has the same face set (intersection region)
    ///     3. FUSE: keep all; COMMON: keep only solids with matching face set in tools;
    ///        CUT: keep only solids WITHOUT matching face set in tools
    /// OCCT-aligned: BuildRC (BOPAlgo_BOP.cxx L583-867).
    ///   Filters split solids by boolean operation type using face-set comparison.
    ///   A. FUSE (L594-609): keep all split solids (fence-deduped).
    ///   B. COMMON/CUT/CUT21 (L616-864): build args/tools building-element maps,
    ///      resolve to split images, compare for intersection containment.
    ///   rcad: result.tmp_solids contains pre-assembled split solids with
    ///     solid_side_origin tracking.  The OCCT myShape is approximated by
    ///     result.tmp_solids entries.
    fn build_rc(&self, result: &mut ResultBuilder, t_brep: &mut topods::BRep) {
        // OCCT L587-591: TopoDS_Compound aC; BRep_Builder aBB; aBB.MakeCompound(aC)

        let solids = std::mem::take(&mut result.tmp_solids);
        let sides: Vec<usize> = result.solid_side_origin.clone();
        if sides.len() != solids.len() { return; }

        // OCCT L594-609: A. FUSE �?iterate myShape with fence, add all
        if self.op == BooleanOpType::Union {
            // OCCT L596: aMFence fence map
            // OCCT L597-606: TopExp_Explorer aExp(myShape, aType); fence-add to aC
            let mut a_m_fence: std::collections::HashSet<Vec<usize>> =
                std::collections::HashSet::new();
            let mut kept: Vec<Vec<usize>> = Vec::new();
            for s in &solids {
                if a_m_fence.insert(s.clone()) {
                    kept.push(s.clone());
                }
            }
            // OCCT L607: myRC = aC
            result.tmp_solids = kept;
            return;
        }

        // OCCT L616-645: prepare building elements of arguments to get splits
        //   OCCT: iterate myArguments/myTools �?TreatCompound �?TopExp::MapShapes
        //   rcad: DS vertices/edges/faces are the building elements.
        //   For each side (0=args, 1=tools): collect V/E/F indices into maps.
        let e_base = self.ds.vertices.len();
        let f_base = e_base + self.ds.edges.len();

        // OCCT L622: aMArgs, aMTools �?indexed maps of source shapes (V/E/F)
        //   rcad: HashSet<usize> of flat V/E/F indices.
        let mut a_m_args: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_m_tools: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_maps = [&mut a_m_args, &mut a_m_tools];

        for (side_idx, a_ms) in a_maps.iter_mut().enumerate() {
            // OCCT L628-643: for each argument/tool shape �?TreatCompound �?MapShapes
            //   rcad: source building elements classified by origin in DS arrays.
            let v_range = if side_idx == 0 {
                (0usize, self.ds.a_vertex_count)
            } else {
                (self.ds.a_vertex_count, self.ds.vertices.len())
            };
            let e_range = if side_idx == 0 {
                (0usize, self.ds.a_edge_count)
            } else {
                (self.ds.a_edge_count, self.ds.edges.len())
            };
            let f_range = if side_idx == 0 {
                (0usize, self.ds.a_face_count)
            } else {
                (self.ds.a_face_count, self.ds.faces.len())
            };

            // OCCT L641-642: TypeToExplore(iDim) �?MapShapes(aSS, aType, aMS)
            //   rcad: each DS entity is a building element by type.
            for vi in v_range.0..v_range.1 { a_ms.insert(vi); }
            for ei in e_range.0..e_range.1 { a_ms.insert(e_base + ei); }
            for fi in f_range.0..f_range.1 { a_ms.insert(f_base + fi); }
        }

        // OCCT L654-705: get splits of building elements
        //   For each building element, check myImages.IsBound �?add split images.
        //   rcad: for edges: self.my_images[b].  for faces: result.face_origins count.
        //   For faces with no images, also build BOPTools_Set for SOLID.
        let mut a_m_args_im: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_m_tools_im: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mset_args: Vec<std::collections::BTreeSet<usize>> = Vec::new();
        let mut a_mset_tools: Vec<std::collections::BTreeSet<usize>> = Vec::new();

        let mut im_maps = [&mut a_m_args_im, &mut a_m_tools_im];
        let mut set_maps = [&mut a_mset_args, &mut a_mset_tools];

        for (side_idx, (a_ms_im, a_mset)) in im_maps.iter_mut().zip(set_maps.iter_mut()).enumerate() {
            let a_ms = &a_maps[side_idx]; // &HashSet<usize> for this side
            let side_is_args = side_idx == 0;

            // OCCT L667-704: for each building element
            let mut sorted_elements: Vec<&usize> = a_ms.iter().collect();
            sorted_elements.sort(); // deterministic order

            for &&flat_idx in &sorted_elements {
                // OCCT L670-678: Type check + degenerated edge skip
                //   rcad: flat_idx < v_range �?VERTEX, < e_range �?EDGE, else FACE
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                let local_idx = if is_edge { flat_idx - e_base }
                    else if is_face { flat_idx - f_base }
                    else { flat_idx };

                if is_edge {
                    // OCCT L671-678: degenerated edge check
                    if self.ds.is_edge_degenerated(local_idx) { continue; }
                }

                // OCCT L681-691: if (myImages.IsBound(aS)) { add split images }
                let has_images = if is_edge {
                    self.my_images.borrow().contains_key(
                        &rcad_kernel::topods::ShapeRef::new(local_idx))
                } else if is_face {
                    // Face has images if DS face produces multiple result faces
                    let (o_exp, sfi) = if side_is_args {
                        (ShapeOrigin::ShapeA, local_idx)
                    } else {
                        (ShapeOrigin::ShapeB, local_idx)
                    };
                    let result_count = result.face_origins.iter().filter(|fo| match fo {
                        FaceOrigin::FromA(s) if side_is_args => *s == sfi,
                        FaceOrigin::FromB(s) if !side_is_args => *s == sfi,
                        _ => false,
                    }).count();
                    result_count > 1
                    // If result_count == 0, the face was not split at all
                } else {
                    // OCCT: VERTEX images from myImages �?not tracked at this level in rcad
                    false
                };

                if has_images {
                    // OCCT L683-689: iterate split images and add to image map
                    let (o_exp, sfi) = if side_is_args {
                        (ShapeOrigin::ShapeA, local_idx)
                    } else {
                        (ShapeOrigin::ShapeB, local_idx)
                    };

                    if is_face {
                        for (rfi, fo) in result.face_origins.iter().enumerate() {
                            let matches = match fo {
                                FaceOrigin::FromA(s) if side_is_args => *s == sfi,
                                FaceOrigin::FromB(s) if !side_is_args => *s == sfi,
                                _ => false,
                            };
                            if matches {
                                a_ms_im.insert(f_base + rfi);
                            }
                        }
                    } else if is_edge {
                        if let Some(imgs) = self.my_images.borrow().get(
                            &rcad_kernel::topods::ShapeRef::new(local_idx))
                        {
                            for &sr in imgs {
                                a_ms_im.insert(e_base + sr.index);
                            }
                        }
                    }
                } else {
                    // OCCT L692-702: no images �?add original shape
                    a_ms_im.insert(flat_idx);

                    // OCCT L694-701: for SOLID building elements, build BOPTools_Set
                    //   rcad: for face elements, build DS face set for BOPTools_Set comparison
                    if is_face {
                        let mut a_st: std::collections::BTreeSet<usize> =
                            std::collections::BTreeSet::new();
                        // Build face set from this face and its adjacent faces in the same solid
                        //   �?rcad: BOPTools_Set at FACE level approximates OCCT's
                        //     SOLID-level BOPTools_Set.  OCCT adds all faces of the SOLID;
                        //     rcad adds the single DS face and its shell siblings.
                        let (o_exp2, sfi2) = if side_is_args {
                            (ShapeOrigin::ShapeA, local_idx)
                        } else {
                            (ShapeOrigin::ShapeB, local_idx)
                        };
                        a_st.insert(local_idx);
                        // Add sibling faces from the same shell
                        for (dfi2, df2) in self.ds.faces.iter().enumerate() {
                            if dfi2 != local_idx && df2.origin == o_exp2
                                && df2.source_shell_idx
                                    == self.ds.faces[local_idx].source_shell_idx
                            {
                                a_st.insert(dfi2);
                            }
                        }
                        if !a_mset.contains(&a_st) {
                            a_mset.push(a_st);
                        }
                    }
                }
            }
        }

        // OCCT L707-783: compare the maps and make the result
        let b_common = self.op == BooleanOpType::Intersection;
        let b_cut21 = false; // �?rcad: CUT21 not supported

        // OCCT L715-720: determine iteration/check maps based on CUT21
        let a_m_it: &std::collections::HashSet<usize> = if b_cut21 { &a_m_tools_im } else { &a_m_args_im };
        let a_m_check: &std::collections::HashSet<usize> = if b_cut21 { &a_m_args_im } else { &a_m_tools_im };
        let a_mset_check: &Vec<std::collections::BTreeSet<usize>> =
            if b_cut21 { &a_mset_args } else { &a_mset_tools };

        // OCCT L724-755: expand sub-shapes for COMMON
        let a_m_it_exp: std::collections::HashSet<usize> = if b_common {
            let mut exp = std::collections::HashSet::new();
            for &&flat_idx in &a_m_it.iter().collect::<Vec<_>>() {
                // OCCT L730-736: expand to lower dimensions via TypeToExplore
                //   rcad: if this is a FACE, include its EDGEs and VERTEXes.
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                if is_face {
                    let local_fi = flat_idx - f_base;
                    if local_fi < self.ds.faces.len() {
                        for &ei in &self.ds.faces[local_fi].boundary_edges {
                            exp.insert(e_base + ei);
                        }
                        for &vi in &self.ds.faces[local_fi].boundary_verts {
                            exp.insert(vi);
                        }
                    }
                } else if is_edge {
                    let local_ei = flat_idx - e_base;
                    if local_ei < self.ds.edges.len() {
                        exp.insert(self.ds.edges[local_ei].start_vertex);
                        exp.insert(self.ds.edges[local_ei].end_vertex);
                    }
                }
                exp.insert(flat_idx);
            }
            exp
        } else {
            a_m_it.clone()
        };

        // OCCT L744-755: expand check side too
        let a_m_check_exp: std::collections::HashSet<usize> = {
            let mut exp = std::collections::HashSet::new();
            for &&flat_idx in &a_m_check.iter().collect::<Vec<_>>() {
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                if is_face {
                    let local_fi = flat_idx - f_base;
                    if local_fi < self.ds.faces.len() {
                        for &ei in &self.ds.faces[local_fi].boundary_edges {
                            exp.insert(e_base + ei);
                        }
                        for &vi in &self.ds.faces[local_fi].boundary_verts {
                            exp.insert(vi);
                        }
                    }
                } else if is_edge {
                    let local_ei = flat_idx - e_base;
                    if local_ei < self.ds.edges.len() {
                        exp.insert(self.ds.edges[local_ei].start_vertex);
                        exp.insert(self.ds.edges[local_ei].end_vertex);
                    }
                }
                exp.insert(flat_idx);
            }
            exp
        };

        // OCCT L757-784: compare building-element images and build keep set.
        //   OCCT iterates aMItExp (V/E/F level images); adds each to aC if it
        //   passes the containment check against the other side.
        //   rcad: operate at the same building-element granularity, then filter
        //   result solids whose constituent face building-elements are in keep_set.
        let mut keep_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &&flat_idx in &a_m_it_exp.iter().collect::<Vec<_>>() {
            // OCCT L762: bContains = aMCheckExp.Contains(aS)
            let mut b_contains = a_m_check_exp.contains(&flat_idx);
            // OCCT L763-768: for SOLIDs, also check BOPTools_Set
            //   rcad: operate at FACE level (no SOLID-level DS entries).
            let is_face = flat_idx >= f_base;
            if !b_contains && is_face {
                let local_fi = flat_idx - f_base;
                let mut a_st: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
                if local_fi < self.ds.faces.len() {
                    a_st.insert(local_fi);
                    for &vi in &self.ds.faces[local_fi].boundary_verts {
                        a_st.insert(vi);
                    }
                    for &ei in &self.ds.faces[local_fi].boundary_edges {
                        a_st.insert(e_base + ei);
                    }
                }
                b_contains = a_mset_check.iter().any(|s| s == &a_st);
            }
            // OCCT L770-783: COMMON �?keep if contained; CUT �?keep if NOT contained
            let keep = if b_common { b_contains } else { !b_contains };
            if keep {
                keep_set.insert(flat_idx);
            }
        }

        // Filter result.tmp_solids: keep solids whose iterate-side face building
        // elements pass the building-element filter above.
        let mut kept_solids: Vec<Vec<usize>> = Vec::new();
        for (i, solid_shells) in solids.iter().enumerate() {
            let side = sides.get(i).copied().unwrap_or(0);
            // Check each solid: if ANY face's building element is in keep_set �?keep.
            // A result solid is kept iff the source face(s) it was split from pass.
            let mut solid_keep = false;
            for &si in solid_shells {
                if let Some(shell_faces) = result.tmp_shells.get(si) {
                    for &rfi in shell_faces {
                        let dfi_opt = match result.face_origins.get(rfi) {
                            Some(FaceOrigin::FromA(sfi)) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                            Some(FaceOrigin::FromB(sfi)) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                            _ => None,
                        };
                        if let Some(dfi) = dfi_opt {
                            let flat_fi = f_base + dfi;
                            if keep_set.contains(&flat_fi) {
                                solid_keep = true;
                                break;
                            }
                        }
                    }
                    if solid_keep { break; }
                }
            }
            if solid_keep {
                kept_solids.push(solid_shells.clone());
            }
        }

        // OCCT L786-809: filter result for COMMON �?re-explore from high dim to low
        //   rcad: OCCT re-iterates the compound by dimension (SOLID→SHELL→FACE) with
        //   a fence.  rcad solids are already at SOLID granularity; shell-count fence
        //   prevents duplicates (OCCT L799-804 fence at FACE+WIRE level).
        if b_common {
            let mut a_m_fence: std::collections::HashSet<Vec<usize>> =
                std::collections::HashSet::new();
            let mut reordered: Vec<Vec<usize>> = Vec::new();
            for s in &kept_solids {
                if a_m_fence.insert(s.clone()) {
                    reordered.push(s.clone());
                }
            }
            kept_solids = reordered;
        }

        // OCCT L811-864: degenerated edge squat (DEs whose vertex is in result,
        //   is not new, and is not interfered �?add to aC).
        //   �?rcad: result edges are embedded in pre-assembled solids.  Adding
        //     standalone DEs to the compound is not applicable at this level.

        result.tmp_solids = kept_solids;
    }

    /// �?OCCT-aligned: FillInternalShapes (Builder_3.cxx L622-887).
    ///   Settle internal sub-shapes (vertices, edges) into result solids.
    ///
    /// OCCT flow:
    ///   L630-655 (Phase 1): Collect V/E/WIRE from arguments with
    ///     TopAbs_INTERNAL orientation inside source solids.
    ///   L680-718 (Phase 2): For each source SOLID, OwnInternalShapes
    ///     collects non-FACE sub-shapes (V/E/WIRE).  Build aMSx ancestry
    ///     map (VERTEX→EDGE, VERTEX→FACE, EDGE→FACE) for split solids.
    ///   L720-746 (Phase 3): Filter �?remove internal shapes already
    ///     attached to split-solid faces (found in aMSx).
    ///   L806-887 (Phase 4): Classify remaining against each split solid
    ///     via ComputeStateByOnePoint.  If IN �?add to that solid with
    ///     TopAbs_INTERNAL orientation.  If the solid is an original (not
    ///     yet having images), clone it first and store in myImages.
    ///
    /// rcad: internal V/E are marked via DSVertex/DSEdge::is_internal
    ///   flag.  Phase 1-2 collect is_internal V/E from the DS arrays.
    ///   Phase 3: no-face-ancestry check �?internal shapes by definition
    ///   have no face references.  Phase 4: classify point against result
    ///   solids' DS face sets via classify_point.  If IN �?the shape is
    ///   recorded on result.face_internal_vtx for the solid's first face
    ///   (OCCT adds it to the TopoDS_Solid as INTERNAL sub-shape).

    /// OCCT-aligned: PerformAreas void detection (BuilderSolid.cxx L397-576).
    fn detect_internal_voids(&self, result: &mut ResultBuilder,
                              assignments: &[(usize, usize, &'static str)]) {
        // OCCT L397-399: myAreas.Clear(); BRep_Builder aBB;
        // OCCT L400-407: aNewSolids, aHoleShells, aMHF (hole face map).

        // Precompute DS face set, centroid, AABB per solid.
        //   OCCT operates on raw shells (myLoops); rcad operates on result.tmp_solids.
        let n_solids = result.tmp_solids.len();
        let mut ds_faces_of: Vec<Vec<usize>> = Vec::with_capacity(n_solids);
        let mut centroids: Vec<DVec3> = Vec::with_capacity(n_solids);
        let mut aabbs: Vec<Aabb> = Vec::with_capacity(n_solids);
        for si in 0..n_solids {
            let mut faces = Vec::new();
            let mut aabb = Aabb::empty();
            let mut centroid = DVec3::ZERO;
            for &sh in &result.tmp_solids[si] {
                if let Some(shell) = result.tmp_shells.get(sh) {
                    for &fi in shell {
                        if let Some(origin) = result.face_origins.get(fi) {
                            let ds_fi = match origin {
                                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                                _ => None,
                            };
                            if let Some(dfi) = ds_fi {
                                faces.push(dfi);
                                for &vi in &self.ds.faces[dfi].boundary_verts {
                                    if vi < self.ds.vertices.len() {
                                        aabb.expand_point(self.ds.vertices[vi].point);
                                    }
                                }
                                if fi < result.faces.len() {
                                    centroid = result.faces[fi].6;
                                }
                            }
                        }
                    }
                }
            }
            faces.sort_unstable(); faces.dedup();
            ds_faces_of.push(faces);
            centroids.push(centroid);
            aabbs.push(aabb);
        }

        // === Step 1: Classify each shell as Growth or Hole (OCCT L411-442) ===
        //   OCCT L422: IsGrowthShell(aShell, aMHF) �?fast face overlap check.
        //     If any face of theShell is already in aMHF (face map of known holes),
        //     the shell is a Growth (it bounds a hole).
        //   OCCT L426: IsHole(aShell, myContext) �?classify infinite point against
        //     the dead solid (original solid being split).  IN = hole, OUT = growth.
        let mut a_mhf: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut is_hole = vec![false; n_solids];

        for si in 0..n_solids {
            // OCCT L422: IsGrowthShell
            let is_growth = if !a_mhf.is_empty() {
                ds_faces_of[si].iter().any(|dfi| a_mhf.contains(dfi))
            } else {
                false
            };

            if !is_growth {
                // OCCT L426: IsHole �?classify against original solid (source operand).
                let side = result.solid_side_origin.get(si).copied().unwrap_or(0);
                let dead_faces: Vec<usize> = self.ds.faces.iter().enumerate()
                    .filter(|(_, f)| match side {
                        0 => f.origin == ShapeOrigin::ShapeA,
                        _ => f.origin == ShapeOrigin::ShapeB,
                    })
                    .map(|(fi, _)| fi)
                    .collect();
                let class = classify_point(centroids[si], &dead_faces, self.ds);
                // OCCT: IsHole returns true if infinite point is IN dead solid �?hole.
                is_hole[si] = class == Classification::In;
            }
            // else: IsGrowthShell returned true �?definitely a growth.

            if is_hole[si] {
                // OCCT L439-441: aHoleShells.Add + TopExp::MapShapes(,TopAbs_FACE,aMHF)
                for &dfi in &ds_faces_of[si] {
                    a_mhf.insert(dfi);
                }
            }
        }

        // OCCT L429-441 (Growth/Hole separation done above).
        let in_si: Vec<usize> = (0..n_solids).filter(|&i| is_hole[i]).collect();
        let out_si: Vec<usize> = (0..n_solids).filter(|&i| !is_hole[i]).collect();

        // OCCT L444-458: if no holes �?add all growths to myAreas + return.
        if in_si.is_empty() || out_si.is_empty() { return; }

        // === Step 2: Build BVH of hole shells (OCCT L462-478) ===
        //   OCCT L464-475: BOPTools_BoxTree with BRepBndLib bounding boxes.
        //   rcad: Bvh built from hole solid AABBs.
        let hole_key: Vec<usize> = in_si.clone();
        let hole_aabbs: Vec<Aabb> = in_si.iter().map(|&i| aabbs[i]).collect();
        let hole_bvh = crate::bvh::DsBvh::build(hole_key, hole_aabbs);

        // === Step 3: Classify holes against growth solids (OCCT L483-529) ===
        //   OCCT L493-529: for each growth solid:
        //     build box �?BVH-select candidate holes �?IsInside �?store outermost.
        let mut a_hole_solid_map: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();

        for &os in &out_si {
            // OCCT L494-497: BRepBndLib::Add(aSolid, aBox)
            // OCCT L499-502: BOPTools_BoxTreeSelector �?candidate holes
            let candidates = hole_bvh.query_aabb(&aabbs[os]);

            for &hole_idx in &candidates {
                // OCCT L511: IsInside(aHole, aSolid, myContext)
                let class = classify_point(centroids[hole_idx], &ds_faces_of[os], self.ds);
                if class != Classification::In && class != Classification::On {
                    continue;
                }

                // OCCT L517-527: select outermost containing solid.
                //   If current os is INSIDE the previously recorded solid,
                //   the current os is more specific (innermost container) �?prefer it.
                use std::collections::hash_map::Entry;
                match a_hole_solid_map.entry(hole_idx) {
                    Entry::Occupied(mut e) => {
                        let prev_os = *e.get();
                        let prev_faces = &ds_faces_of[prev_os];
                        let os_inside = classify_point(centroids[os], prev_faces, self.ds);
                        if os_inside == Classification::In || os_inside == Classification::On {
                            e.insert(os);
                        }
                    }
                    Entry::Vacant(e) => {
                        e.insert(os);
                    }
                }
            }
        }

        // === Step 4: Build reverse map: solid �?list of holes (OCCT L532-548) ===
        let mut solid_holes_map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (&hole_idx, &os) in &a_hole_solid_map {
            solid_holes_map.entry(os).or_default().push(hole_idx);
        }

        // === Step 5: Add holes to solids + myAreas (OCCT L550-576) ===
        let mut removed = vec![false; n_solids];
        // OCCT L553-573: for each growth with holes �?aBB.Add(aSolid, aHole)
        for (&os, holes) in &solid_holes_map {
            for &hole_idx in holes {
                let void_shells = result.tmp_solids[hole_idx].clone();
                result.tmp_solids[os].extend(void_shells);
                removed[hole_idx] = true;
            }
        }
        // OCCT L575: myAreas.Append(aSolid) �?rcad: non-removed solids kept in tmp_solids.
        // OCCT L578-581: add un-associated holes to myAreas (rcad: not needed �?kept as-is).
        let mut new_solids: Vec<Vec<usize>> = Vec::with_capacity(n_solids);
        for (si, solid) in result.tmp_solids.drain(..).enumerate() {
            if !removed[si] { new_solids.push(solid); }
        }
        result.tmp_solids = new_solids;
    }

    /// �?OCCT-aligned: FillInternalShapes (Builder_3.cxx L622-887).
    /// OCCT-aligned: FillInternalShapes (Builder_3.cxx L622-887).
    ///   Phase 1 (L648-709): Collect internal V/E/WIRE from arguments.
    ///   Phase 2 (L717-788): Internal V/E from source solids + build aMSx ancestry.
    ///   Phase 3 (L790-809): Filter shapes already attached via aMSx.
    ///   Phase 4 (L811-816): Early return if none.
    ///   Phase 5 (L820-877): Classify each internal shape against each split solid;
    ///     if IN �?add as INTERNAL sub-shape (clone original if needed).
    fn fill_internal_shapes(&self, result: &mut ResultBuilder) {
        // OCCT L631-644: allocator + indexed maps (aMSx, aMx, aMSI, aMFence, aMSOr, ...)
        //   rcad: adapted to Vec/HashSet equivalents.

        // === Phase 1: Shapes to process �?collect from arguments (OCCT L648-709) ===
        //   OCCT L653-658: TreatCompound on each argument �?flatten into aLSC.
        //   OCCT L660-681: filter VERTEX/EDGE/WIRE from aLSC �?aLArgs.
        //   OCCT L684-709: for each aLArgs, check myImages.IsBound �?aMSI (images or originals).
        //   rcad: DS vertices/edges with is_internal flag = sources.
        //   �?TreatCompound: rcad treats DS V/E as already-flattened source shapes.
        //   �?aMSI: maps shape-ref �?true if it's an internal shape to process.
        let mut a_msi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // Collect internal vertices (OCCT L677-679: TopAbs_VERTEX �?aLArgs)
        for (vi, v) in self.ds.vertices.iter().enumerate() {
            if v.is_internal {
                // OCCT L691-706: check myImages.IsBound �?add split images or original
                let v_ref = rcad_kernel::topods::ShapeRef::new(vi);
                if self.my_images.borrow().contains_key(&v_ref) {
                    for img in &self.my_images.borrow()[&v_ref] {
                        a_msi.insert(img.index);
                    }
                } else {
                    a_msi.insert(vi);
                }
            }
        }
        // Collect internal edges (OCCT L665-675: WIRE �?iterate edges; L677-679: EDGE directly)
        for (ei, e) in self.ds.edges.iter().enumerate() {
            if e.is_internal {
                let e_ref = rcad_kernel::topods::ShapeRef::new(ei);
                if self.my_images.borrow().contains_key(&e_ref) {
                    for img in &self.my_images.borrow()[&e_ref] {
                        a_msi.insert(img.index);
                    }
                } else {
                    a_msi.insert(ei);
                }
            }
        }

        // === Phase 2: Internal V/E from source solids + build aMSx ancestry (OCCT L717-788) ===
        //   OCCT L721-727: iterate DS for SOLIDs.
        //   L738: OwnInternalShapes(aS, aMx) �?get INTERNAL sub-shapes from each solid.
        //   L741-758: insert into aMSI (with myImages check).
        //   L760-787: build aMSx ancestry: Vertex→Edge, Vertex→Face, Edge→Face.
        //   rcad: aMSx tracks which internal shapes are already on split-solid faces.
        #[allow(unused)]
        let mut a_msx: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new(); // shape_idx �?list of ancestor face/edge indices
        let mut a_lsd: Vec<usize> = Vec::new(); // split solids to process

        // OCCT L741-758: internal shapes from OwnInternalShapes
        //   rcad: is_internal flag already collected above in Phase 1.
        //   The DS vertices/edges with is_internal=true are equivalent to OCCT's OwnInternalShapes output.

        // OCCT L760-787: build aMSx for split solids
        //   For each source SOLID that has split results (images) �?build ancestry map.
        //   rcad: iterate result.tmp_solids �?for each solid, map edges→faces.
        for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
            // OCCT L761: if (myImages.IsBound(aS)) for source solid
            let side = result.solid_side_origin.get(si).copied().unwrap_or(0);
            // Build edge→face adjacency for this result solid
            let mut edge_to_faces: std::collections::HashMap<usize, Vec<usize>> =
                std::collections::HashMap::new();
            for &shi in solid_shells {
                if let Some(shell_faces) = result.tmp_shells.get(shi) {
                    for &rfi in shell_faces {
                        if let Some(fe) = result.faces.get(rfi) {
                            for &(ei, _) in &fe.0 {
                                edge_to_faces.entry(ei).or_default().push(rfi);
                            }
                        }
                    }
                }
            }
            // OCCT L770-773: TopExp::MapShapesAndAncestors �?aMSx
            //   aMSx[vertex] = list of edge indices
            //   aMSx[vertex] = list of face indices
            //   aMSx[edge] = list of face indices
            for (&ei, face_list) in &edge_to_faces {
                // e_ref in a_msx �?ancestors (face indices)
                a_msx.entry(ei).or_default().extend(face_list);
                // Also add vertex→edge ancestry
                if ei < self.ds.edges.len() {
                    a_msx.entry(self.ds.edges[ei].start_vertex)
                        .or_default().push(ei);
                    a_msx.entry(self.ds.edges[ei].end_vertex)
                        .or_default().push(ei);
                }
            }
            a_lsd.push(si);
        }

        // === Phase 3: Filter shapes already attached to split-solid faces (OCCT L790-809) ===
        //   OCCT: for each shape in aMSI, check if aMSx.Contains(shape) with non-empty ancestor list.
        //         �?if NOT attached �?aLSI (list of shapes to settle).
        let mut a_lsi: Vec<usize> = Vec::new();
        for &si in &a_msi {
            // OCCT L796-808: if aMSx contains the shape AND has non-empty ancestors �?skip (attached).
            //   rcad: check if this internal shape index appears in aMSx with non-empty ancestors.
            let is_attached = a_msx.get(&si).map_or(false, |anc| !anc.is_empty());
            if !is_attached {
                a_lsi.push(si);
            }
        }

        // === Phase 4: Early return if none (OCCT L811-816) ===
        if a_lsi.is_empty() {
            return;
        }

        // === Phase 5: Settle internal V/E into solids (OCCT L820-877) ===
        //   OCCT L825-876: for each split solid (aLSd), for each internal shape (aLSI):
        //     ComputeStateByOnePoint(aSI, aSd) �?if IN:
        //       - if original solid (aMSOr): clone �?add INTERNAL �?bind myImages/myOrigins
        //       - else: add INTERNAL directly
        for &si in &a_lsd {
            // OCCT L828: TopoDS_Solid aSd
            //   rcad: get DS face set for this solid
            let mut solid_ds_faces: Vec<usize> = Vec::new();
            if let Some(solid_shells) = result.tmp_solids.get(si) {
                for &shi in solid_shells {
                    if let Some(shell_faces) = result.tmp_shells.get(shi) {
                        for &rfi in shell_faces {
                            let dfi_opt = match result.face_origins.get(rfi) {
                                Some(FaceOrigin::FromA(sfi)) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                                Some(FaceOrigin::FromB(sfi)) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                                _ => None,
                            };
                            if let Some(dfi) = dfi_opt {
                                solid_ds_faces.push(dfi);
                            }
                        }
                    }
                }
            }
            solid_ds_faces.sort_unstable();
            solid_ds_faces.dedup();
            if solid_ds_faces.is_empty() { continue; }

            // OCCT L830-875: iterate internal shapes to settle
            let mut i = 0usize;
            while i < a_lsi.len() {
                let si_idx = a_lsi[i];
                // OCCT L834: aSI.Orientation(TopAbs_INTERNAL)
                //   rcad: no orientation; use classify_point with centroid.
                // OCCT L836: ComputeStateByOnePoint(aSI, aSd, 1.e-11, myContext)
                let pt = if si_idx < self.ds.vertices.len() {
                    self.ds.vertices[si_idx].point
                } else {
                    let ei = si_idx;
                    if ei < self.ds.edges.len() {
                        (self.ds.vertices[self.ds.edges[ei].start_vertex].point
                         + self.ds.vertices[self.ds.edges[ei].end_vertex].point) * 0.5
                    } else {
                        i += 1; continue;
                    }
                };
                let a_state = classify_point(pt, &solid_ds_faces, self.ds);

                if a_state != Classification::In {
                    // OCCT L840: aIt1.Next(); continue;
                    i += 1;
                    continue;
                }

                // OCCT L844-873: shape is IN �?add as INTERNAL
                //   OCCT L844: if (aMSOr.Contains(aSd)) �?original solid �?clone first
                //   rcad: find first face of this solid to store internal vertex
                if let Some(&first_shi) = result.tmp_solids.get(si).and_then(|s| s.first()) {
                    if let Some(shell_faces) = result.tmp_shells.get(first_shi) {
                        if let Some(&first_rfi) = shell_faces.first() {
                            if first_rfi < result.face_internal_vtx.len() {
                                // OCCT L857-873: add INTERNAL shape to solid
                                //   rcad: store DS vertex index in face_internal_vtx
                                if si_idx < self.ds.vertices.len() {
                                    if !result.face_internal_vtx[first_rfi].contains(&si_idx) {
                                        result.face_internal_vtx[first_rfi].push(si_idx);
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L875: aLSI.Remove(aIt1) �?remove settled shape
                a_lsi.swap_remove(i);
                // don't increment i �?the new element at i needs checking too
            }
        }
    }

    /// �?OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-342).
    ///   Phase 7: group result solids into COMPSOLID/COMPOUND hierarchy.
    ///
    /// OCCT flow:
    ///   L200-217 (FillImagesCompounds): Iterate source shapes for TopAbs_COMPOUND.
    ///     For each compound, call FillImagesCompound recursively.
    ///   L280-342 (FillImagesCompound): Recursively check each child for images.
    ///     If any child has images, build a new compound with image replacements.
    ///     Result stored in myImages[original_compound] = new_compound.
    ///
    /// rcad: records compound intent in ResultBuilder.  Actual compound
    ///   reconstruction happens after result.build() in build_with_history
    ///   (see the rebuild_compound_for_step post-step) because the result
    ///   BRep solids don't exist until build() is called.
    /// �?OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-217).
    ///
    /// OCCT FillImagesCompounds L197-217:
    ///   L200: aMFP fence map
    ///   L202-216: iterate source shapes, filter TopAbs_COMPOUND,
    ///             call FillImagesCompound(aC, aMFP)
    /// OCCT FillImagesCompound L280-342:
    ///   L290-293: fence �?skip if processed
    ///   L296-308: check if any sub-shape has images
    ///   L309-312: if none modified �?return
    ///   L314-341: build new compound from sub-shape (solid) images
    ///
    /// rcad: source compound solids are tracked in DS solid_images.
    ///   Compound reconstruction from result solids is deferred to
    ///   build_with_history's post-step (L6834-6840) because the
    ///   result BRep solids don't exist until ResultBuilder::build().
    /// �?OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-217) + FillImagesCompound (L280-342).
    ///   L197-201: aMFP fence map; NbSourceShapes �?filter TopAbs_COMPOUND.
    ///   L280-293: FillImagesCompound �?fence skip if already processed.
    ///   L295-308: recurse into sub-compounds; check if any sub-shape has images.
    ///   L309-312: no modification �?return.
    ///   L314-341: build new compound from sub-shape images; store in myImages.
    /// OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-217) + FillImagesCompound (L280-342).
    ///   L197-201: dispatcher with fence map; iterate source COMPOUND shapes.
    ///   L280-293: FillImagesCompound �?fence skip if already processed.
    ///   L295-308: recurse into sub-compounds; check if any sub-shape has images.
    ///   L309-312: no modification �?return.
    ///   L314-341: build new compound from sub-shape images; store in myImages.
    ///   �?rcad: no compound nesting in DS.  Flat per-face source_compsolid_idx.
    ///     The recursive FillImagesCompound is collapsed to a single level.
    fn fill_images_compounds(&self, result: &mut ResultBuilder, t: &mut topods::BRep) {
        // OCCT L199-200: aMFP fence map �?prevents reprocessing the same compound.
        //   rcad: HashSet of processed compsolid indices.
        let mut a_mfp: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // OCCT L202: aNbS = myDS->NbSourceShapes() �?iterate all DS shapes.
        //   rcad: collect unique source_compsolid_idx from DS faces.
        let mut compound_indices: Vec<usize> = Vec::new();
        for df in &self.ds.faces {
            if let Some(csi) = df.source_compsolid_idx {
                if !compound_indices.contains(&csi) {
                    compound_indices.push(csi);
                }
            }
        }
        if compound_indices.is_empty() { return; }

        for &csi in &compound_indices {
            // OCCT L290-293: if (!theMFP.Add(theS)) return �?fence check.
            if !a_mfp.insert(csi) { continue; }

            // OCCT L295-308: check if any sub-shape (SOLID) has been modified.
            //   rcad: collect source_solid_idx values under this compsolid,
            //   check if any has images (multiple result solids).
            let sub_solid_indices: Vec<usize> = self.ds.faces.iter()
                .filter_map(|f| {
                    if f.source_compsolid_idx == Some(csi) {
                        f.source_solid_idx
                    } else {
                        None
                    }
                })
                .collect();
            // OCCT L300-303: recurse into sub-compounds �?rcad: flat, no nesting.
            // OCCT L304-307: if (myImages.IsBound(aSx)) bInterferred = true.
            //   rcad: check if any sub-solid produces >1 result solid (split).
            let mut b_interferred = false;
            for &ssi in &sub_solid_indices {
                // Count result solids from this source solid
                let count = result.solid_side_origin.iter()
                    .filter(|&&side| {
                        // Count result solids from the side matching this source solid's origin
                        let dfi = self.ds.faces.iter().position(|f|
                            f.source_solid_idx == Some(ssi));
                        dfi.map_or(false, |di| {
                            let origin = &self.ds.faces[di].origin;
                            (origin == &crate::bopds::ds::ShapeOrigin::ShapeA && side == 0)
                                || (origin == &crate::bopds::ds::ShapeOrigin::ShapeB && side == 1)
                        })
                    })
                    .count();
                if count > 0 {
                    b_interferred = true;
                    break;
                }
            }

            // OCCT L309-312: if (!bInterferred) return �?no modification.
            if !b_interferred { continue; }

            // OCCT L314-315: MakeContainer(COMPOUND, aCIm)
            //   rcad: collect result solid indices for this compsolid.
            let mut a_c_im: Vec<usize> = Vec::new();

            // OCCT L317-336: iterate sub-shapes �?add images or original.
            for &ssi in &sub_solid_indices {
                // Find the DS face for this source solid to determine its side (origin)
                let side = self.ds.faces.iter()
                    .find(|f| f.source_solid_idx == Some(ssi))
                    .map(|f| match f.origin {
                        crate::bopds::ds::ShapeOrigin::ShapeA => 0,
                        crate::bopds::ds::ShapeOrigin::ShapeB => 1,
                    })
                    .unwrap_or(0);

                // OCCT L322: if (myImages.IsBound(aSX)) �?has split images?
                //   rcad: check if result solids exist for this side+source solid.
                let matching_solids: Vec<usize> = result.solid_side_origin.iter()
                    .enumerate()
                    .filter(|&(_, &s)| s == side)
                    .map(|(si, _)| si)
                    .collect();

                if matching_solids.is_empty() {
                    // OCCT L334-335: no images �?add original sub-shape
                    //   rcad: no solid to add �?the original solid is implicit.
                    continue;
                }

                // OCCT L324-331: has images �?add each image with orientation.
                for &si in &matching_solids {
                    if !a_c_im.contains(&si) {
                        // OCCT L329: aSXIm.Orientation(aOrX) �?preserve orientation.
                        //   rcad: orientation is per-face via FaceOrigin.
                        a_c_im.push(si);
                    }
                }
            }

            // OCCT L339-341: aLSIm.Append(aCIm); myImages.Bind(theS, aLSIm)
            //   rcad: create TShape::Compound from result solid ShapeRefs.
            if !a_c_im.is_empty() {
                let solid_refs: Vec<topods::ShapeRef> = a_c_im.iter()
                    .filter_map(|&si| result.solids.get(si).copied())
                    .collect();
                if !solid_refs.is_empty() {
                    let cmp_ref = t.add_tcompound(solid_refs);
                    let ckey = topods::ShapeRef::new(usize::MAX - a_c_im.len());
                    self.my_images.borrow_mut().entry(ckey).or_default().push(cmp_ref);
                    result.compsolid_groups.push(cmp_ref);
                }
            }
        }
    }

    /// Retrieve the EdgeInfo.is_inside status for the incoming edge at the given vertex.
    fn incoming_edge_is_inside(&self, smart_map: &IndexMap<usize, Vec<EdgeInfo>>, vertex: usize, seg_idx: usize) -> bool {
        smart_map.get(&vertex)
            .and_then(|infos| infos.iter().find(|ei| ei.seg_idx == seg_idx && ei.in_flag))
            .map_or(false, |ei| ei.is_inside)
    }

    /// �?OCCT-aligned: face keep/discard policy (ComputeState �?FillIn3DParts equivalent).
    ///   OCCT does NOT have a surface-type special case �?ComputeState propagates
    ///   ON→IN/OUT based on face orientation + solid side, not surface type.
    /// �?OCCT-aligned: BOPAlgo_Builder::FillImagesFaces �?face keep policy.
    ///   OCCT: after ComputeState returns IN/OUT/ON for a face against the other solid:
    ///     FUSE: keep OUT + ON
    ///     COMMON: keep IN + ON
    ///     CUT A-B:
    ///       face from A �?keep if OUT or ON (A outside B)
    ///       face from B �?keep if IN or ON (B inside A, the cut surface)
    fn classification_keep_policy(&self, source: SourceSide, class: Classification, _fi: usize) -> bool {
        match self.op {
            BooleanOpType::Intersection => class == Classification::In || class == Classification::On,
            BooleanOpType::Difference => match source {
                SourceSide::A => class != Classification::In,
                SourceSide::B => class == Classification::In || class == Classification::On,
            },
            BooleanOpType::Union => class != Classification::In,
        }
    }

    /// �?OCCT-aligned: BuildResult �?add split images to result (Builder_1.cxx L130-168).
    ///   OCCT: for each source shape of theType, if myImages bound �?add images;
    ///   else add the original shape.  rcad: for Edge, creates topods edges in t_brep
    ///   (equivalent to OCCT's myShape) AND flat edge refs in result for face construction.
    ///   For Vertex/Wire/Shell/Solid, rcad handles these in other pipeline steps.
    /// OCCT-aligned: BuildResult (Builder_1.cxx L130-168).
    ///   Add split images (or originals) of source shapes into the result.
    ///   OCCT L133: aMFence fence map.
    ///   L136-167: for each source argument of matching type �?if myImages bound
    ///     �?add all image shapes; else �?add the original shape.
    ///   rcad: adapts to topods::BRep TShape factory + ResultBuilder storage.
    fn build_result(&self, shape_type: ShapeType, result: &mut ResultBuilder, t: &mut topods::BRep) {
        // OCCT L133: NCollection_Map<TopoDS_Shape> aMFence �?dedup shapes in result.
        //   rcad: unique-indexed arrays make a fence unnecessary, but the form is kept.
        #[allow(unused)]
        let mut a_m_fence: Vec<usize> = Vec::new();

        // OCCT L136-167: iterate all source arguments of matching type.
        //   rcad: source entities vary by type (DS arrays, result data).
        match shape_type {
    ShapeType::Vertex => {
        // �?OCCT-aligned: BuildResult (Builder_1.cxx L130-168).
        //   Iterate all source arguments of type VERTEX �?add images to myShape.
        //   If no images, add the original vertex.
        //   rcad: source vertices = DS vertices 0..a_vc (A) + a_vc.. (B).
        let a_vc = self.ds.a_vertex_count;
        let nv = self.ds.vertices.len();
        for side in 0..2usize {
            let (start, end) = if side == 0 { (0usize, a_vc.min(nv)) } else { (a_vc, nv) };
            for vi in start..end {
                // OCCT L145: myImages.Seek(aS) �?check if vertex has split image
                let sref = rcad_kernel::topods::ShapeRef::new(vi);
                let has_images = self.my_images.borrow().contains_key(&sref);
                if !has_images {
                    // OCCT L149-152: no images �?add the original shape
                    let pt = self.ds.vertices[vi].point;
                    let _rvi = result.add_ds_vertex(vi, pt);
                    t.add_tvertex(pt);
                } else {
                    // OCCT L156-165: add images of the argument shape into result
                    let images = self.my_images.borrow().get(&sref).unwrap().clone();
                    for img in &images {
                        let vi_img = img.index;
                        if vi_img < self.ds.vertices.len() {
                            let pt = self.ds.vertices[vi_img].point;
                            let _rvi = result.add_ds_vertex(vi_img, pt);
                            t.add_tvertex(pt);
                        }
                    }
                }
            }
        }
    }
            ShapeType::Edge => {
                // OCCT L130-168 (TopAbs_EDGE): add split edge images to myShape.
                //   rcad: iterate myImages(EDGE) entries, create TShape::Edge for each.
                //   First ensure all DS vertices have TShapes for edge vertex refs.
                let e_base = self.ds.vertices.len();
                // Ensure vertex TShapes exist (needed by edge creation)
                for vi in 0..self.ds.vertices.len() {
                    let vr = rcad_kernel::topods::ShapeRef::new(vi);
                    if t.tshapes.len() <= vi {
                        // Extend tshapes array to cover this index
                        let pt = self.ds.vertices[vi].point;
                        let sv = t.add_tvertex(pt);
                        t.vertex_mut(sv).tolerance = self.ds.vertices[vi].geom_tol
                            .max(crate::tolerance::TOLERANCE_ABS);
                        let _ = vr;
                    }
                }
                // Iterate A and B side source edges
                let a_ec = self.ds.a_edge_count;
                let n_edges = self.ds.edges.len();
                for side in 0..2usize {
                    let (start, end) = if side == 0 {
                        (0usize, a_ec.min(n_edges))
                    } else {
                        (a_ec, n_edges)
                    };
                    for ei in start..end {
                        let aE = rcad_kernel::topods::ShapeRef::new(e_base + ei);
                        let has_images = self.my_images.borrow().contains_key(&aE);
                        if !has_images {
                            // OCCT L149-152: no images �?add original edge
                            let edge = &self.ds.edges[ei];
                            let sv_sr = rcad_kernel::topods::ShapeRef::new(edge.start_vertex);
                            let ev_sr = rcad_kernel::topods::ShapeRef::new(edge.end_vertex);
                            let ci = t.curves.len();
                            t.curves.push(edge.curve.clone());
                            let te = t.add_tedge(Some(ci), sv_sr, ev_sr, edge.t_range);
                            if self.ds.is_edge_degenerated(ei) || edge.start_vertex == edge.end_vertex {
                                t.edge_mut(te).degenerated = true;
                            }
                        } else {
                            // OCCT L156-165: add split images
                            let images = self.my_images.borrow().get(&aE).unwrap().clone();
                            for img in &images {
                                let nSpR = img.index.saturating_sub(e_base);
                                if nSpR >= self.ds.edges.len() { continue; }
                                let edge = &self.ds.edges[nSpR];
                                let sv_sr = rcad_kernel::topods::ShapeRef::new(edge.start_vertex);
                                let ev_sr = rcad_kernel::topods::ShapeRef::new(edge.end_vertex);
                                let ci = t.curves.len();
                                t.curves.push(edge.curve.clone());
                                let te = t.add_tedge(Some(ci), sv_sr, ev_sr, edge.t_range);
                                if self.ds.is_edge_degenerated(nSpR) || edge.start_vertex == edge.end_vertex {
                                    t.edge_mut(te).degenerated = true;
                                }
                            }
                        }
                    }
                }
            }
            ShapeType::Wire => {
                // OCCT L130-168 (TopAbs_WIRE): create wire TShapes from DSWire edges.
                // rcad: for each DSWire, build a topods::Wire using e_base mapping.
                // Edges already exist in t_brep from BuildResult(EDGE).
                let e_base = self.ds.vertices.len();
                result.wire_refs = Vec::with_capacity(self.ds.wires.len());
                for wi in 0..self.ds.wires.len() {
                    let w_ref = rcad_kernel::topods::ShapeRef::new(e_base + self.ds.edges.len() + wi);
                    // Check wire image for replacement edges
                    let mut wire_edges: Vec<rcad_kernel::topods::ShapeRef> = Vec::new();
                    if let Some(imgs) = self.my_images.borrow().get(&w_ref) {
                        for &img_sr in imgs {
                            wire_edges.push(img_sr);
                        }
                    } else {
                        for &ei in &self.ds.wires[wi].edges {
                            wire_edges.push(rcad_kernel::topods::ShapeRef::new(e_base + ei));
                        }
                    }
                    // Create TShape::Wire in t_brep (but only if we have edges)
                    if !wire_edges.is_empty() {
                        let sr = t.add_twire(wire_edges);
                        t.wire_mut(sr).closed = true;
                        result.wire_refs.push(sr);
                    } else {
                        // No edges — skip this wire
                        result.wire_refs.push(topods::ShapeRef::NULL);
                    }
                }
            }
            ShapeType::Face => {
                // OCCT L145-165: for each source FACE, check myImages.Seek(aS).
                //   rcad: result.face_origins tracks which source faces were split.
                //   rcad: build_faces() validates edge refs before TShape creation.
                result.build_faces();
                // OCCT-aligned: use myImages to decide which source faces had split images.
                let f_base = self.ds.vertices.len() + self.ds.edges.len();
                let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
                let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
                for &fi in &a_faces {
                    let aS = rcad_kernel::topods::ShapeRef::new(f_base + fi);
                    if !self.my_images.borrow().contains_key(&aS) {
                        result.build_original_face(self.ds, fi,
                            FaceOrigin::FromA(self.ds.faces[fi].source_face_idx));
                    }
                }
                for &fi in &b_faces {
                    let aS = rcad_kernel::topods::ShapeRef::new(f_base + fi);
                    if !self.my_images.borrow().contains_key(&aS) {
                        result.build_original_face(self.ds, fi,
                            FaceOrigin::FromB(self.ds.faces[fi].source_face_idx));
                    }
                }
                result.build_topods_faces(t);
            }
            ShapeType::Shell => {
                // OCCT L145-165: for each source SHELL, check myImages for shell images.
                //   rcad: fill_images_containers_shells already created TShape::Shell
                //   entries in result.shells and my_images.  BuildResult(SHELL) ensures
                //   the shells are finalized (no additional work needed here since
                //   the shell creation in fill_images_containers_shells already handles
                //   the OCCT L274-275 container creation step).
            }
            ShapeType::Solid => {
                // OCCT L130-167: for each source SOLID, check myImages → add images/original.
                //   ⏳ OCCT-aligned BuildSolid: for FUSE, collect ALL result faces into one solid
                //   (BOPAlgo_BOP.cxx L902-906).  Non-union solids are already created by
                //   build_split_solids (stored in result.solids + my_images).
                if self.op == BooleanOpType::Union && !result.face_refs.is_empty() {
                    let sf = result.face_refs.clone();
                    let shell_ref = t.add_tshell(sf);
                    let solid_ref = t.add_tsolid(vec![shell_ref]);
                    let so_key = rcad_kernel::topods::ShapeRef::new(usize::MAX - 1 - result.solids.len());
                    self.my_images.borrow_mut().entry(so_key).or_default().push(solid_ref);
                    result.solids.push(solid_ref);
                }
            }
            ShapeType::CompSolid => {
                // OCCT L130-167: aggregate sub-solid images into CompSolid.
                //   rcad: fill_images_containers_compsolid already created
                //   TShape::CompSolid entries in result.compsolid_groups.
            }
            ShapeType::Compound => {
                // OCCT L130-168: for each source COMPOUND, add its image (or original).
                //   rcad: fill_images_compounds already created TShape::Compound
                //   entries in result.compsolid_groups and my_images.
            }
        }
    }

    /// �?OCCT-aligned: BOPAlgo_Builder::BuildResult (Builder_1.cxx L130-168).
    ///   rcad: thin wrapper mapping topods::ShapeType to builder::types::ShapeType.
    fn build_result_occt(&self, the_type: topods::ShapeType, result: &mut ResultBuilder, t: &mut topods::BRep) {
        let shape_type = match the_type {
            topods::ShapeType::Shape => unreachable!("ShapeType::Shape is a null sentinel, never passed to build_result"),
            topods::ShapeType::Vertex => ShapeType::Vertex,
            topods::ShapeType::Edge => ShapeType::Edge,
            topods::ShapeType::Wire => ShapeType::Wire,
            topods::ShapeType::Face => ShapeType::Face,
            topods::ShapeType::Shell => ShapeType::Shell,
            topods::ShapeType::Solid => ShapeType::Solid,
            topods::ShapeType::CompSolid => ShapeType::CompSolid,
            topods::ShapeType::Compound => ShapeType::Compound,
        };
        self.build_result(shape_type, result, t);
    }

    /// �?OCCT-aligned: BOPAlgo_BOP::BuildShape (BOP.cxx L871-906).
    ///   Calls BuildRC (L900) then BuildSolid for FUSE 3D (L902-906).
    fn build_shape(&self, result: &mut ResultBuilder, t_brep: &mut topods::BRep) {
        // OCCT L900: BuildRC �?filter solids by boolean operation
        self.build_rc(result, t_brep);
        if self.has_errors { return; }
        // OCCT L902-906: if (FUSE + 3D) BuildSolid
        //   rcad: Union keeps all filtered solids; no separate BuildSolid needed.
    }

    /// �?OCCT-aligned: BOPAlgo_Builder::PostTreat (Builder.cxx L450-475).
    ///   Two-step tolerance correction: CorrectTolerances + CorrectShapeTolerances.
    fn post_treat(&self, brep: &mut rcad_kernel::BRep) {
        // OCCT L452-454: aMA �?map of shapes to avoid
        // OCCT L455-469: if non-destructive �?collect source V/E/F into aMA
        // rcad: non-destructive defaults to false.  When true, collect non-new
        // DS vertex indices into map_to_avoid.
        let map_to_avoid: std::collections::HashSet<usize> = if self.my_non_destructive {
            let mut avoid = std::collections::HashSet::new();
            for (vi, v) in self.ds.vertices.iter().enumerate() {
                if v.origin.is_some() { avoid.insert(vi); }
            }
            for (ei, e) in self.ds.edges.iter().enumerate() {
                if matches!(e.origin, ShapeOrigin::ShapeA | ShapeOrigin::ShapeB) {
                    avoid.insert(ei);
                }
            }
            for (fi, f) in self.ds.faces.iter().enumerate() {
                if matches!(f.origin, ShapeOrigin::ShapeA | ShapeOrigin::ShapeB) {
                    avoid.insert(fi);
                }
            }
            avoid
        } else {
            std::collections::HashSet::new()
        };
        // OCCT L472: BOPTools_AlgoTools::CorrectTolerances(myShape, aMA, 0.05, myRunParallel)
        if map_to_avoid.is_empty() {
            rcad_kernel::tolerance::correct_tolerances(brep, 23);
        } else {
            rcad_kernel::tolerance::correct_tolerances_with_map(brep, 23, &map_to_avoid);
        }
        // OCCT L474: BOPTools_AlgoTools::CorrectShapeTolerances(myShape, aMA, myRunParallel)
        //   rcad: correct_tolerances already does both steps in one call.
        //   Separating them requires splitting the tolerance module.
    }
}
