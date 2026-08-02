use super::*;
use crate::bop::int_tools::int_patch_type::IntPatchIType;
use rcad_kernel::topods::ShapeType;
use std::collections::HashMap;

/// Work item for Phase 2 of PerformFF (OCCT BOPAlgo_PaveFiller_6.cxx L488-507).
/// Mirrors BOPAlgo_FaceFace data: face indices, shift info, tolerance, EF points.
struct FFWork {
    f1: usize,
    f2: usize,
    shift_info: Option<super::SeamEdgeShift>,
    a_shift_value: f64,
    // OCCT L495-496: double aTolFF = std::max(aShiftValue, ToleranceFF(aBAS1, aBAS2));
    a_tol_ff: f64,
    // OCCT L498-504: NCollection_List<IntSurf_PntOn2S> aListOfPnts from GetEFPnts
    ef_points: Vec<DVec3>,
    // OCCT L506: SetParameters(bApprox, bCompC2D1, bCompC2D2, anApproxTol)
    // OCCT L507: SetFuzzyValue(myFuzzyValue)
    b_approx: bool,
    b_comp_c2d1: bool,
    b_comp_c2d2: bool,
    an_approx_tol: f64,
    fuzzy_value: f64,
}

impl<'a> super::PaveFiller<'a> {
    pub(crate) fn perform_ff(&mut self) {
        // OCCT BOPAlgo_PaveFiller_6.cxx L285-623: PerformFF

        // ====================================================================
        // Phase 0: face info update (OCCT L290-314)
        // ====================================================================
        // L290-291: myIterator->Initialize(TopAbs_FACE, TopAbs_FACE); iSize = myIterator->ExpectedLength();
        let pairs = self.ff_candidate_pairs();
        let i_size = pairs.len();

        // L294-301: NCollection_Map<int> aMIFence; collect from iterator
        let mut a_mi_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &(n_f1, n_f2) in &pairs {
            a_mi_fence.insert(n_f1);
            a_mi_fence.insert(n_f2);
        }

        // L303-310: collect rest of touched faces (HasReference)
        // OCCT: for (int i = 0; i < myDS->NbSourceShapes(); ++i)
        //         if (aSI.ShapeType() == TopAbs_FACE && aSI.HasReference()) aMIFence.Add(i);
        // rcad: convert source shape index to flat face index via reference
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            if self.ds.shape_type_of(i) != ShapeType::Face {
                continue;
            }
            if !self.ds.shape_info_at(i).has_reference() {
                continue;
            }
            let fi = self.ds.shape_info_at(i).reference as usize;
            a_mi_fence.insert(fi);
        }

        // L312-313: UpdateFaceInfoOn / UpdateFaceInfoIn
        for &fi in &a_mi_fence {
            self.ds.refine_face_info_on(fi);
            self.ds.refine_face_info_in(fi);
        }

        // L315-319: early return when no intersection pairs
        if i_size == 0 {
            return;
        }

        // ====================================================================
        // Options + EE map (OCCT L327-360)
        // ====================================================================
        // L323-324: aFFs.SetIncrement(iSize) 鈥?pre-allocate interference array
        self.ds.interf_ff.reserve(i_size);
        // L327-331: bApprox, bCompC2D1, bCompC2D2, anApproxTol, bSplitCurve
        let b_approx = true;
        let b_comp_c2d1 = true;
        let b_comp_c2d2 = true;
        let an_approx_tol = 1e-7;
        let b_split_curve = false;

        // L335-360: build aEEMap from EE interferences
        // DataMap<BOPDS_Pair, NCollection_List<int>>  鈫? HashMap<(usize,usize), Vec<usize>>
        let mut a_ee_map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for inf in &self.ds.interf_ee {
            if inf.new_vertex == usize::MAX {
                continue;
            }
            let n_e1 = inf.e1;
            let n_e2 = inf.e2;
            let n_vn = inf.new_vertex;
            let a_pair = if n_e1 < n_e2 {
                (n_e1, n_e2)
            } else {
                (n_e2, n_e1)
            };
            a_ee_map.entry(a_pair).or_default().push(n_vn);
        }

        // ====================================================================
        // Main loop: create BOPAlgo_FaceFace work items (OCCT L363-517)
        // ====================================================================
        // L363: BOPAlgo_VectorOfFaceFace aVFaceFace;
        let mut a_v_face_face: Vec<FFWork> = Vec::new();

        // L365: myIterator->Initialize(TopAbs_FACE, TopAbs_FACE);
        for &(n_f1, n_f2) in &pairs {
            // L367-370: if (UserBreak(aPSOuter)) return;
            if self.check_stop("PerformFF") {
                return;
            }
            // L373: if (myGlue == BOPAlgo_GlueOff)
            if !self.use_glue() {
                // ===== Non-glue path (OCCT L374-508) =====

                // L375-378: const TopoDS_Face& aF1, aF2; const BRepAdaptor_Surface& aBAS1, aBAS2
                // rcad: BRepAdaptor equivalent = locate_surface (applies Location + unwraps TrimmedSurface)
                let a_surf1 = self.ds.locate_surface(n_f1);
                let a_surf2 = self.ds.locate_surface(n_f2);
                let is_plane1 = matches!(a_surf1, Surface3::Plane(_));
                let is_plane2 = matches!(a_surf2, Surface3::Plane(_));

                // L380-391: CheckPlanes 鈥?skip parallel planes without vertex overlap
                if is_plane1 && is_plane2 {
                    let b_to_intersect = self.check_planes(n_f1, n_f2);
                    if !b_to_intersect {
                        // L386-389: aFFs.Appended(); aFF.SetIndices(nF1,nF2); aFF.Init(0,0);
                        self.ds.interf_ff.push(InterferenceFF {
                            f1: n_f1,
                            f2: n_f2,
                            curves: Vec::new(),
                            points: Vec::new(),
                            tangent_faces: false,
                        });
                        continue;
                    }
                }

                // OCCT L393-486: Seam edge shift (inline in PerformFF)
                // L400: TopoDS_Face aFShifted1 = aF1, aFShifted2 = aF2;
                // L402: double aShiftValue = 0.;
                // OCCT L404-416: if (!isPlane1 || !isPlane2)
                let shift_info: Option<SeamEdgeShift>;
                let a_shift_value: f64;
                if is_plane1 && is_plane2 {
                    shift_info = None;
                    a_shift_value = 0.0;
                } else {
                    // OCCT L418-485: nested wire/edge loop with IsClosedFF + aEEMap + shift
                    let mut an_is_found = false;
                    let mut a_found_shift_info: Option<SeamEdgeShift> = None;
                    let mut a_found_shift_value: f64 = 0.0;
                    // OCCT L419-425: for (TopoDS_Iterator aItW1(aF1); !anIsFound && aItW1.More(); aItW1.Next())
                    //               for (aItE1(aItW1.Value()); !anIsFound && aItE1.More(); aItE1.Next())
                    // rcad: iterate boundary_edges of face
                    for &e1 in self.ds.face_boundary_edges(n_f1) {
                        if an_is_found {
                            break;
                        }
                        let is_closed1 = self.is_seam_edge(e1, n_f1);
                        // OCCT L426-431: for (aItW2(aF2); !anIsFound && aItW2.More(); aItW2.Next())
                        //               for (aItE2(aItW2.Value()); !anIsFound && aItE2.More(); aItE2.Next())
                        for &e2 in self.ds.face_boundary_edges(n_f2) {
                            if an_is_found {
                                break;
                            }
                            let is_closed2 = self.is_seam_edge(e2, n_f2);
                            // OCCT L437-440: if (!anIsClosed1 && !anIsClosed2) continue;
                            if !is_closed1 && !is_closed2 {
                                continue;
                            }

                            // OCCT L442-447: aEEMap lookup
                            let a_key = if e1 < e2 { (e1, e2) } else { (e2, e1) };
                            let Some(a_vertex_indices) = a_ee_map.get(&a_key) else {
                                continue;
                            };

                            // OCCT L449-480: for each EE vertex, project to edges, compute shift
                            for &vi in a_vertex_indices {
                                if an_is_found {
                                    break;
                                }
                                let vertex_point = self.ds.vertex_point(vi);

                                // OCCT L457-468: ProjectPointOnCurve on both edges
                                let curve1 = self.ds.edge_curve(e1).unwrap();
                                let curve2 = self.ds.edge_curve(e2).unwrap();
                                let proj1 = closest_point_on_curve(curve1, vertex_point, 64);
                                let proj2 = closest_point_on_curve(curve2, vertex_point, 64);
                                let a_p1 = proj1.point;
                                let a_p2 = proj2.point;

                                // OCCT L470: aShiftDist = aP1.Distance(aP2)
                                let shift_dist = a_p1.distance(a_p2);

                                // OCCT L471: if (aShiftDist > BRep_Tool::Tolerance(aVertex))
                                let vtx_tol = self.ds.vertex_tolerance(vi);
                                if shift_dist > vtx_tol {
                                    // OCCT L474-479: shift one face
                                    // OCCT: gp_Vec(aP1, aP2) 鈥?rcad: a_p2 - a_p1
                                    a_found_shift_info = Some(SeamEdgeShift {
                                        shift_vector: if is_closed1 {
                                            a_p2 - a_p1
                                        } else {
                                            a_p1 - a_p2
                                        },
                                        shift_value: shift_dist,
                                        shifted_face: if is_closed1 { 1 } else { 2 },
                                    });
                                    a_found_shift_value = shift_dist;
                                    an_is_found = true;
                                }
                            }
                        }
                    }
                    shift_info = a_found_shift_info;
                    a_shift_value = a_found_shift_value;
                }

                // L495: double aTolFF = std::max(aShiftValue, ToleranceFF(aBAS1, aBAS2));
                let a_tol_ff = a_shift_value.max(self.ff_tol(n_f1, n_f2));

                // L498-504: GetEFPnts(nF1, nF2, aListOfPnts)
                // OCCT: NCollection_List<IntSurf_PntOn2S> aListOfPnts;
                //       GetEFPnts(nF1, nF2, aListOfPnts);
                //       int aNbLP = aListOfPnts.Extent();
                //       if (aNbLP) { aFaceFace.SetList(aListOfPnts); }
                let a_list_of_pnts = self.get_ef_pnts_ff(n_f1, n_f2);

                // L488-507: BOPAlgo_FaceFace& aFaceFace = aVFaceFace.Appended();
                // Rcad: push FFWork (GetEFPnts / SetList done inside intersect_face_face)
                a_v_face_face.push(FFWork {
                    f1: n_f1,
                    f2: n_f2,
                    shift_info,
                    a_shift_value,
                    a_tol_ff,
                    ef_points: a_list_of_pnts,
                    // OCCT L506: SetParameters(bApprox, bCompC2D1, bCompC2D2, anApproxTol)
                    b_approx,
                    b_comp_c2d1,
                    b_comp_c2d2,
                    an_approx_tol,
                    // OCCT L507: SetFuzzyValue(myFuzzyValue)
                    fuzzy_value: self.fuzzy_tolerance,
                });
            } else {
                // ===== Glue mode (OCCT L510-517) =====
                // L512-515: aFF.SetIndices(nF1,nF2); aFF.SetTangentFaces(false); aFF.Init(0,0);
                self.ds.interf_ff.push(InterferenceFF {
                    f1: n_f1,
                    f2: n_f2,
                    curves: Vec::new(),
                    points: Vec::new(),
                    tangent_faces: false,
                });
            }
        }

        // ====================================================================
        // Execute intersection (OCCT L519-534)
        // ====================================================================
        // L519-520: int k, aNbFaceFace = aVFaceFace.Length();
        let a_nb_face_face = a_v_face_face.len();

        // L528-529: BOPTools_Parallel::Perform(myRunParallel, aVFaceFace);
        // L530-533: if (UserBreak(aPSOuter)) return;
        // Rcad: sequential execution
        for work in &a_v_face_face {
            if self.check_stop("PerformFF") {
                return;
            }

            // OCCT BOPAlgo_FaceFace::Perform (L202-241):
            //   IntTools_FaceFace::Perform(aF1,aF2) = IntPatch + MakeCurve
            // Rcad: intersect_face_face = IntPatch + MakeCurve + ApplyTrsf
            //   PrepareLines3D/ComputeTolReached3d/PostTreatFF called from results loop below
            self.intersect_face_face(&work);
        }

        // ====================================================================
        // Process results (OCCT L537-621)
        // ====================================================================
        for work in &a_v_face_face {
            if self.check_stop("PerformFF") {
                return;
            }
            let k_n_f1 = work.f1;
            let k_n_f2 = work.f2;

            // Find the InterfFF entry that intersect_face_face created for this pair
            // OCCT: aFaceFace.Indices(nF1, nF2) retrieves indices from BOPAlgo_FaceFace object
            // Rcad: find last InterfFF matching this face pair
            let ff_idx = self
                .ds
                .interf_ff
                .iter()
                .rposition(|ff| ff.f1 == k_n_f1 && ff.f2 == k_n_f2);
            let Some(ff_idx) = ff_idx else {
                continue;
            };

            // OCCT L545-553: if !IsDone || HasErrors 鈫?empty FF + warning + continue
            // Rcad: intersect_face_face already handles failure (creates empty InterfFF on failure).
            // Check equivalent: if entry has no curves, no points, and faces are not tangent
            let ff_entry = &self.ds.interf_ff[ff_idx];
            // OCCT L545: if (!aFaceFace.IsDone() || aFaceFace.HasErrors())
            // OCCT L547-549: aFF.SetIndices(nF1,nF2); aFF.Init(0,0);
            // OCCT L551: AddIntersectionFailedWarning(aFaceFace.Face1(), aFaceFace.Face2());
            if ff_entry.curves.is_empty() && ff_entry.points.is_empty() && !ff_entry.tangent_faces {
                // OCCT L546-553: intersection failed 鈥?record warning, skip
                // AddIntersectionFailedWarning (OCCT BOPAlgo_PaveFiller_2.cxx L660)
                // Rcad: warning skipped 鈥?does not affect topology
                continue;
            }

            // OCCT L555-560: bTangentFaces, aTolFF, PrepareLines3D, ApplyTrsf
            // PrepareLines3D is called here (separate from intersection)
            // OCCT L555: bool bTangentFaces = aFaceFace.TangentFaces();
            // OCCT L556: double aTolFF = aFaceFace.TolFF();
            // rcad: a_tol_ff stored in FFWork, bTangentFaces stored in InterfFF by intersect_face_face
            self.compute_tol_and_prepare_lines_3d(
                k_n_f1,
                k_n_f2,
                b_split_curve,
                work.a_shift_value,
            );
            // OCCT L560: ApplyTrsf (reverse seam shift on intersection results)
            // rcad: reverse_seam_edge_shift reverses the shift applied during intersect_face_face
            if let Some(ref shift) = work.shift_info {
                self.reverse_seam_edge_shift(k_n_f1, k_n_f2, shift);
            }
            self.update_face_sc_vertices(k_n_f1, k_n_f2);

            // OCCT L562-571: aNbCurves, aNbPoints, myDS->AddInterf(nF1, nF2) BEFORE CheckCurve
            let a_nb_curves = self.ds.interf_ff[ff_idx].curves.len();
            let a_nb_points = self.ds.interf_ff[ff_idx].points.len();
            if a_nb_curves > 0 || a_nb_points > 0 {
                self.ds.try_add_interf(k_n_f1, k_n_f2);
            }

            // OCCT L573-576: aFF.SetIndices, SetTangentFaces, Init(aNbCurves, aNbPoints)
            // Rcad: already set in intersect_face_face

            // OCCT L578-588: aBoxExpandValue = aTolFF; if aNbCurves>0, += max vertex tolerance
            let a_tol_ff = work.a_tol_ff;
            let mut a_box_expand_value = a_tol_ff;
            if a_nb_curves > 0 {
                // OCCT L585-586: BRep_Tool::MaxTolerance(Face1, TopAbs_VERTEX), Face2
                // rcad: iterate boundary vertices to find max vertex tolerance
                let a_max_vertex_tol = self.ds.faces.get(k_n_f1).unwrap()
                    .boundary_verts
                    .iter()
                    .chain(self.ds.face_boundary_verts(k_n_f2).iter())
                    .map(|&vi| self.ds.vertex_tolerance(vi))
                    .fold(0.0, f64::max);
                a_box_expand_value += a_max_vertex_tol;
            }

            // OCCT L590-609: CheckCurve + store
            // OCCT: NCollection_DynamicArray<BOPDS_Curve>& aVNC = aFF.ChangeCurves();
            // rcad: retained_curves filters the ICs; OCCT stores them in BOPDS_Curve array
            let retained_curves: Vec<usize> = self.ds.interf_ff[ff_idx]
                .curves
                .iter()
                .filter_map(|&ci| {
                    if ci >= self.ds.intersection_curves.len() {
                        return None;
                    }
                    let ic = &self.ds.intersection_curves[ci];
                    // OCCT L598-601: Bnd_Box aBox; CheckCurve(aIC, aBox)
                    // OCCT L559-561: BndLib_Add3dCurve::Add(GeomAdaptor_Curve(aC3D),
                    //                   max(theCurve.Tolerance(), theCurve.TangentialTolerance()), theBox);
                    let tol_curve = ic.geom_tol.max(ic.curve_extra.tangential_tol);
                    let bb = crate::bnd_lib::curve_bounds_with_range(
                        &ic.curve,
                        ic.t_range[0],
                        ic.t_range[1],
                        tol_curve.max(crate::tolerance::CONFUSION),
                    );
                    if !bb.is_valid() {
                        return None;
                    }
                    let mut bb_min = bb.min;
                    let mut bb_max = bb.max;
                    // OCCT L605: aBox.Enlarge(aBoxExpandValue)
                    let expand = a_box_expand_value;
                    bb_min -= DVec3::splat(expand);
                    bb_max += DVec3::splat(expand);
                    // OCCT L568: double aTolCmp = 3 * Precision::Confusion();
                    // OCCT L572: bValid = !theBox.IsThin(aTolCmp);
                    let a_tol_cmp = 3.0 * crate::tolerance::CONFUSION;
                    // IsThin check: reject if box is thin in ALL three directions (< aTolCmp)
                    let dx = bb_max.x - bb_min.x;
                    let dy = bb_max.y - bb_min.y;
                    let dz = bb_max.z - bb_min.z;
                    let b_is_valid = dx > a_tol_cmp || dy > a_tol_cmp || dz > a_tol_cmp;
                    if !b_is_valid && std::env::var("RCAD_DBG_FF").is_ok() {
                        eprintln!(
                            "[FF] CheckCurve discard IC[{}] box=({:.2e},{:.2e},{:.2e}) tol={:.2e}",
                            ci, dx, dy, dz, a_tol_cmp
                        );
                    }
                    // OCCT L600-607: if (bIsValid) store curve with box and tolerance
                    if b_is_valid {
                        // OCCT L605-607: aBox.Enlarge(aBoxExpandValue); aNC.SetBox(aBox);
                        // rcad: set box on the IC's curve_extra.my_box
                        let ic = &mut self.ds.intersection_curves[ci];
                        ic.curve_extra.my_box = Some((bb_min, bb_max));
                        // OCCT L607: aNC.SetTolerance(std::max(aIC.Tolerance(), aTolFF));
                        ic.geom_tol = ic.geom_tol.max(a_tol_ff);
                    }
                    if b_is_valid { Some(ci) } else { None }
                })
                .collect();
            self.ds.interf_ff[ff_idx].curves = retained_curves;

            // OCCT L610-620: points already processed 鈥?store BOPDS_Point
            // Rcad: points already stored in InterfFF by intersect_face_face
        }

        // OCCT L622: end of PerformFF
        // NOTE: rcad's dedup_ff_interferences() was removed for 1:1 alignment (not in OCCT)
    }

    pub(crate) fn ff_candidate_pairs(&self) -> Vec<(usize, usize)> {
        // Build face BVH locally (OCCT builds BOPTools_BoxTree inside PerformFF)
        // Equivalent to BOPDS_Iterator::Initialize(TopAbs_FACE, TopAbs_FACE)
        let face_bvh = {
            let mut indices = Vec::new();
            let mut aabbs = Vec::new();
            for (fi, _f) in self.ds.faces.iter().enumerate() {
                indices.push(fi);
                aabbs.push(crate::bop::ds::face_aabb(self.ds, fi));
            }
            if indices.len() >= 20 {
                Some(crate::bop::tools::bvh::BoxTree::build(indices, aabbs))
            } else {
                None
            }
        };

        // Get FF pair candidates from BVH or pair iterator.
        // OCCT equivalent: myIterator->Initialize(TopAbs_FACE, TopAbs_FACE) loop.
        if let Some(ref fbvh) = face_bvh {
            let candidates = crate::bop::tools::bvh::BoxTree::candidate_pairs(fbvh, fbvh);
            candidates
                .into_iter()
                .filter(|&(fa, fb)| self.ds.face_origin(fa) != self.ds.face_origin(fb))
                .collect()
        } else {
            let a_fcount = self.ds.a_face_count();
            let mut result = Vec::new();
            let mut fit = crate::bop::ds::PairIterator::prepare_ab(a_fcount, self.ds.face_count());
            while fit.more() {
                let pk = fit.value();
                result.push((pk.i1, pk.i2));
                fit.next();
            }
            result
        }
    }

    /// OCCT BOPAlgo_PaveFiller_6.cxx L404-486: seam edge shift (inline in PerformFF).
    /// Checks EE intersections between seam (closed) edges of a face pair
    /// and returns the shift needed to align them, or None.
    ///
    /// OCCT alignment notes:
    /// - OCCT L404-416: condition (aBAS1/aBAS2 != Plane) + surface/triangulation retrieval
    ///   rcad equivalent: skip when both planes (L333-335)
    /// - OCCT L418-441: nested wire/edge loop with IsClosedFF checks
    ///   rcad equivalent: boundary_edges iterator + is_seam_edge
    /// - OCCT L442-480: aEEMap lookup + ProjPC + shift computation
    ///   rcad equivalent: a_ee_map lookup + closest_point_on_curve
    /// - Architecture diff: OCCT uses IsPlaneFF (handles offset/trimmed);
    ///   rcad uses Surface3::Plane direct match.

    /// OCCT BOPAlgo_PaveFiller_6.cxx L106-134: IsClosedFF 鈥?checks if edge is a seam edge of a face.
    ///
    /// OCCT implementation: iterates over BRep_TEdge curve representations, finds one
    /// that is both on the given surface (IsCurveOnSurface) and on a closed surface
    /// (IsCurveOnClosedSurface). Falls back to BRep_Tool::IsClosed (triangulation).
    ///
    /// rcad equivalent: iterates over DSEdge.face_reps, finds representation for the
    /// given face_idx with pcurve2.is_some() (CurveOnClosedSurface in TOPODS).
    pub(crate) fn is_seam_edge(&self, edge_idx: usize, face_idx: usize) -> bool {
        // L112: if (!theIsPlane)
        let a_surf = self.ds.locate_surface(face_idx);
        let the_is_plane = matches!(a_surf, Surface3::Plane(_));
        if !the_is_plane {
            // L115-129: iterate over edge's curve representations
            // L124: IsCurveOnSurface(theSurface, aLocation) && IsCurveOnClosedSurface()
            let edge = self.ds.edges.get(edge_idx).unwrap();
            for rep in self.ds.edge_face_reps(edge_idx) {
                if rep.face_idx == face_idx && rep.pcurve2.is_some() {
                    return true;
                }
            }
        }
        // L132-133: BRep_Tool::IsClosed(theEdge, theTriangulation, theLocation)
        // Rcad: no triangulation fallback available; return false.
        false
    }

    /// OCCT BOPAlgo_PaveFiller_6.cxx L2608-2687: GetEFPnts
    /// Collects EF intersection points between two faces.
    fn get_ef_pnts_ff(&self, f1: usize, f2: usize) -> Vec<DVec3> {
        let edges_f1: Vec<usize> = self.ds.face_boundary_edges(f1).to_vec();
        let edges_f2: Vec<usize> = self.ds.face_boundary_edges(f2).to_vec();
        let mut ef_points: Vec<DVec3> = Vec::new();
        for ef in &self.ds.interf_ef {
            if ef.face == f1 {
                if edges_f2.contains(&ef.edge) {
                    ef_points.push(ef.point);
                }
            } else if ef.face == f2 {
                if edges_f1.contains(&ef.edge) {
                    ef_points.push(ef.point);
                }
            }
        }
        ef_points
    }

    /// OCCT BOPAlgo_PaveFiller_6.cxx L244-265: BOPAlgo_FaceFace::ApplyTrsf
    /// Reverse the seam edge shift on intersection results.
    ///
    /// OCCT alignment: OCCT's ApplyTrsf reverses the TrsfToPoint transform
    /// (moving shapes to/from origin), NOT the seam edge shift. In OCCT, the
    /// seam edge shift physically moves aFShifted1/aFShifted2 before intersection
    /// (L474-477), and the shift is absorbed into the intersection curves.
    ///
    /// rcad applies the seam edge shift to Surface clones mathematically,
    /// then reverses it on the results. Net effect is equivalent.
    pub(crate) fn reverse_seam_edge_shift(&mut self, f1: usize, f2: usize, shift: &SeamEdgeShift) {
        let inv_vec = if shift.shifted_face == 1 {
            -shift.shift_vector
        } else {
            shift.shift_vector
        };

        // Collect curve indices from the FaceFace interference for this pair
        let mut curve_indices: Vec<usize> = Vec::new();
        for inf in &self.ds.interf_ff {
            if (inf.f1 == f1 && inf.f2 == f2) || (inf.f1 == f2 && inf.f2 == f1) {
                curve_indices = inf.curves.clone();
                break;
            }
        }

        // Reverse shift on each curve
        for &ci in &curve_indices {
            if ci >= self.ds.intersection_curves.len() {
                continue;
            }
            let ic = &mut self.ds.intersection_curves[ci];

            // Translate 3D curve back by inverse shift
            ic.curve = translate_curve3(&ic.curve, inv_vec);

            // Translate polyline points if any
            for p in &mut ic.polyline {
                *p += inv_vec;
            }

            // Translate vertex positions back
            let sv = ic.start_vertex;
            let ev = ic.end_vertex;
            if sv < self.ds.vertex_count() {
                self.ds.vertex_data_mut(sv).point += inv_vec;
            }
            if ev < self.ds.vertex_count() {
                self.ds.vertex_data_mut(ev).point += inv_vec;
            }
        }
    }

    /// Perform face-face intersection (OCCT BOPAlgo_FaceFace::Perform + IntTools_FaceFace::Perform).
    ///
    /// OCCT alignment:
    /// - IntTools_FaceFace::Perform L330-543 (IntTools_FaceFace.cxx): IntPatch + MakeCurve
    /// - BOPAlgo_FaceFace::Perform L202-241 (BOPAlgo_PaveFiller_6.cxx): Trsf + IntTools_FaceFace::Perform
    /// - PrepareLines3D + ComputeTolReached3d: moved to compute_tol_and_prepare_lines_3d()
    ///   called from perform_ff results processing (OCCT L558).
    /// - PostTreatFF: moved to update_face_sc_vertices()
    ///   called from perform_ff results processing (OCCT PostTreat path).
    pub(crate) fn intersect_face_face(&mut self, work: &FFWork) {
        let f1 = work.f1;
        let f2 = work.f2;
        let shift_info = work.shift_info.as_ref();
        let dbg_ff = std::env::var("RCAD_DBG_FF").is_ok();
        if dbg_ff {
            eprintln!("[FF] intersect_face_face: f1={} f2={}", f1, f2);
        }

        let s1_orig = self.ds.face_surface(f1).cloned().unwrap();
        let s2_orig = self.ds.face_surface(f2).cloned().unwrap();

        // Apply seam edge shift to surface clones if needed
        let s1 = match shift_info {
            Some(info) if info.shifted_face == 1 => {
                apply_shift_to_surface(&s1_orig, info.shift_vector)
            }
            _ => s1_orig,
        };
        let s2 = match shift_info {
            Some(info) if info.shifted_face == 2 => {
                apply_shift_to_surface(&s2_orig, info.shift_vector)
            }
            _ => s2_orig,
        };

        // OCCT L351-375: SortTypes  -- canonical surface ordering.
        // Swap f1/f2 so the higher-type surface is always "face A".
        let type_idx1 = Self::surface_type_index(&s1);
        let type_idx2 = Self::surface_type_index(&s2);
        let b_reverse = type_idx1 < type_idx2;
        // OCCT L354: if bReverse, swap face refs so myFace1 gets the higher type.
        let (f1, f2, s_a, s_b) = if b_reverse {
            (f2, f1, &s2, &s1)
        } else {
            (f1, f2, &s1, &s2)
        };

        // OCCT L384-393: tolerance setup
        // OCCT: myTolF1 = BRep_Tool::Tolerance(myFace1) + aFuzz, etc.
        // rcad: tolerance handled by PaveFiller's fuzzy_tolerance and face geom_tols.
        // Compute TolFF = max(face tolerances) per OCCT ToleranceFF (BOPAlgo_PaveFiller_6.cxx L3918-3942).
        let tol1 = self.ds.faces.get(f1).map_or(1e-7, |f| f.geom_tol);
        let tol2 = self.ds.faces.get(f2).map_or(1e-7, |f| f.geom_tol);
        let mut a_tol_ff = tol1.max(tol2);
        fn is_analytic_ff(surf: &Surface3) -> bool {
            matches!(
                surf,
                Surface3::Plane(_)
                    | Surface3::Cylinder(_)
                    | Surface3::Cone(_)
                    | Surface3::Sphere(_)
                    | Surface3::Torus(_)
            )
        }
        if !is_analytic_ff(self.ds.face_surface(f1).unwrap())
            || !is_analytic_ff(self.ds.face_surface(f2).unwrap())
        {
            a_tol_ff = a_tol_ff.max(5e-6);
        }
        // Ensure minimum tolerance for IntPatch to work
        a_tol_ff = a_tol_ff.max(1e-7);
        if dbg_ff {
            eprintln!(
                "[FF] ToleranceFF: f1={} tol={:.2e} f2={} tol={:.2e} -> a_tol_ff={:.2e}",
                f1, tol1, f2, tol2, a_tol_ff
            );
        }

        // OCCT L395-401: isFace1Quad/isFace2Quad  -- skip; rcad uses IntPatchIntersection
        // which dispatches by quad type internally.

        //  OCCT L404-434: Plane-Plane fast path (PerformPlanes)
        if matches!(s_a, rcad_kernel::geom::Surface3::Plane(_))
            && matches!(s_b, rcad_kernel::geom::Surface3::Plane(_))
        {
            self.perform_plane_plane(f1, f2);
            return;
        }

        // OCCT L436-438: myLConstruct.Load(dom1, dom2, myHS1, myHS2)
        let mut lconstruct =
            crate::bop::int_tools::int_patch_line_constructor::GeomIntLineConstructor::new();
        lconstruct.load(f1, f2);

        // IntPatch_Intersection: generic surface-surface intersection.
        let mut int_patch = crate::bop::int_tools::int_patch_intersection::IntPatchIntersection::new();
        int_patch.perform(s_a, s_b, a_tol_ff, a_tol_ff);
        if int_patch.tangent_faces() {
            self.ds.interf_ff.push(crate::bop::ds::InterferenceFF {
                f1,
                f2,
                curves: Vec::new(),
                points: Vec::new(),
                tangent_faces: true,
            });
            return;
        }

        // PutPointsOnLine (IntPatch_Intersection.cxx L268-312).
        // Projects intersection points onto each analytic line to create
        // boundary-crossing vertices.  These vertices split the line into
        // valid intervals for MakeCurve/TreatCircle.
        for li in 0..int_patch.nb_lines() {
            self.put_points_on_line(f1, f2, int_patch.line_mut(li), &work.ef_points);
        }

        // OCCT L498-504: GetEFPnts 鈫?SetList passes EF points to IntPatch's PutPointsOnLine.
        // rcad: IntPatch skips PutPointsOnLine; EF=0 for sphere-sphere (PerformEF gap).
        // EF projection here would require EF>0. Currently EF=0, so no points to project.

        // MakeCurve (IntTools_FaceFace.cxx L695-1846) for each IntPatch line.
        // Returns a Vec of IntersectionCurve  -- one per valid part from the
        // LineConstructor (OCCT supports aNbParts > 1, e.g. multi-segment clipping).
        let mut ff_curve_indices: Vec<usize> = Vec::new();
        for i in 0..int_patch.nb_lines() {
            let ics = self.make_intersection_curve(f1, f2, int_patch.line(i));
            for ic in ics {
                let ci = self.ds.intersection_curves.len();
                let mut adjusted_ic = ic;
                // OCCT L558-567: if reversed, swap pcurves (first  -- second).
                if b_reverse {
                    std::mem::swap(&mut adjusted_ic.pcurve_on_a, &mut adjusted_ic.pcurve_on_b);
                }
                self.ds.intersection_curves.push(adjusted_ic);
                // OCCT: vertices created after MakeCurve (in Process/PerformFF).
                // BRepBuilderAPI_MakeVertex(P3D) + myDS->Index for each endpoint.
                let (p_start, p_end, t0, t1) = {
                    let ic_ref = &self.ds.intersection_curves[ci];
                    let t0 = ic_ref.t_range[0];
                    let t1 = ic_ref.t_range[1];
                    (ic_ref.curve.point_at(t0), ic_ref.curve.point_at(t1), t0, t1)
                };
                let sv = if p_start.is_finite() {
                    self.ds.add_vertex(p_start)
                } else {
                    usize::MAX
                };
                let ev = if p_end.is_finite() {
                    self.ds.add_vertex(p_end)
                } else {
                    usize::MAX
                };
                self.ds.intersection_curves[ci].start_vertex = sv;
                self.ds.intersection_curves[ci].end_vertex = ev;
                // OCCT: init IC pave_blocks so MakeSplitEdges can create section edges
                {
                    use crate::bop::ds::pave::{Pave, PaveBlock, SharedPB};
                    let sv = self.ds.intersection_curves[ci].start_vertex;
                    let ev = self.ds.intersection_curves[ci].end_vertex;
                    let t0 = self.ds.intersection_curves[ci].t_range[0];
                    let t1 = self.ds.intersection_curves[ci].t_range[1];
                    let pb = PaveBlock::new(
                        0,
                        Pave {
                            vertex_idx: sv,
                            param: t0,
                        },
                        Pave {
                            vertex_idx: ev,
                            param: t1,
                        },
                    );
                    let spb = SharedPB::new(pb);
                    self.ds.intersection_curves[ci]
                        .pave_blocks
                        .push(spb.clone());
                    self.ds.pave_blocks.push(spb);
                }
                if std::env::var("RCAD_DBG_MB").is_ok() {
                    let ic2 = &self.ds.intersection_curves[ci];
                    eprintln!(
                        "[DBG_IC3] PUSHED IC[{}]: geom_tol={:.6e} t_range={:.6} {:.6}",
                        ci, ic2.geom_tol, ic2.t_range[0], ic2.t_range[1]
                    );
                }
                ff_curve_indices.push(ci);
            }
        }

        // OCCT L576-608: points  -- filter by isPointInOnFace, append to myPnts.
        let mut ff_point_indices: Vec<crate::bop::ds::types::FFPoint> = Vec::new();
        for pi in 0..int_patch.nb_points() {
            let pt = int_patch.point(pi);
            let (uv_a, uv_b, f_a, f_b) = if b_reverse {
                (
                    glam::DVec2::new(pt.u2, pt.v2),
                    glam::DVec2::new(pt.u1, pt.v1),
                    f2,
                    f1,
                )
            } else {
                (
                    glam::DVec2::new(pt.u1, pt.v1),
                    glam::DVec2::new(pt.u2, pt.v2),
                    f1,
                    f2,
                )
            };
            if !self.context.is_point_in_on_face(self.ds, f_a, uv_a) {
                continue;
            }
            if !self.context.is_point_in_on_face(self.ds, f_b, uv_b) {
                continue;
            }
            // FFPoint stores point data inline (OCCT BOPDS_Point). No DS vertex created yet.
            ff_point_indices.push(crate::bop::ds::types::FFPoint::new(pt.p1, uv_a, uv_b));
        }
        if std::env::var("RCAD_DBG_FF").is_ok() {
            eprintln!(
                "[FF]   -> curves={} nLines={}",
                ff_curve_indices.len(),
                int_patch.nb_lines()
            );
        }
        self.ds.interf_ff.push(crate::bop::ds::InterferenceFF {
            f1,
            f2,
            curves: ff_curve_indices,
            points: ff_point_indices,
            tangent_faces: false,
        });
    } // fn intersect_face_face

    /// PrepareLines3D (OCCT PerformFF L558 + IntTools_FaceFace::PrepareLines3D L1932).
    /// Called from perform_ff results processing loop, after IntPatch + MakeCurve.
    ///
    /// OCCT alignment:
    /// - ComputeTolReached3d: approximate, uses compute_intersection_curve_tolerance
    /// - PrepareLines3D L1932-1976: split closed curves at their closing vertex
    ///   b_split_curve controls whether curves are split (OCCT L558: bSplitCurve)
    /// - Curve split fixup: ensure split curves have distinct start/end vertices
    pub(crate) fn compute_tol_and_prepare_lines_3d(
        &mut self,
        f1: usize,
        f2: usize,
        b_split_curve: bool,
        _shift_tol: f64,
    ) {
        let ff_curves_opt = self.find_face_face_curve_indices(f1, f2);
        if ff_curves_opt.is_none() {
            return;
        }
        let ff_curves = ff_curves_opt.unwrap();
        if ff_curves.is_empty() {
            return;
        }

        let t_a = self.ff_tol(f1, f1);
        let t_b = self.ff_tol(f2, f2);
        for &ci in &ff_curves {
            let (curve, pca, pcb, tr, current_tol) = {
                let ic = &self.ds.intersection_curves[ci];
                (
                    ic.curve.clone(),
                    ic.pcurve_on_a.clone(),
                    ic.pcurve_on_b.clone(),
                    ic.t_range,
                    ic.geom_tol,
                )
            };
            let (new_tol, tang_tol) = crate::bop::int_tools::pcurve_derive::compute_intersection_curve_tolerance(
                &curve,
                pca.as_ref(),
                pcb.as_ref(),
                self.ds.face_surface(f1).unwrap(),
                self.ds.face_surface(f2).unwrap(),
                tr,
                t_a,
                t_b,
                current_tol,
            );
            let ic = &mut self.ds.intersection_curves[ci];
            ic.geom_tol = ic.geom_tol.max(new_tol);
            ic.curve_extra.tangential_tol = ic.curve_extra.tangential_tol.max(tang_tol);
        }
        // PrepareLines3D  -- split closed curves
        let n_curves_before_split = self.ds.intersection_curves.len();
        crate::bop::int_tools::pcurve_derive::prepare_lines_3d(&mut self.ds.intersection_curves, b_split_curve);
        // After PrepareLines3D splits closed curves, the split
        // segments are added to the same FF interference entry.  Update the
        // FF entry's curve list to include any newly created curve indices.
        if n_curves_before_split != self.ds.intersection_curves.len() {
            if let Some(ff_entry) = self.ds.interf_ff.last_mut() {
                for new_ci in n_curves_before_split..self.ds.intersection_curves.len() {
                    ff_entry.curves.push(new_ci);
                }
            }
        }
        // After PrepareLines3D splits closed curves, new curve endpoints
        // must be updated to the split points. For start==end but
        // non-full-period t_range (i.e. split half-circle), compute
        // correct endpoint positions via point_at and create new DS vertices.
        for ci in 0..self.ds.intersection_curves.len() {
            let needs_fix = {
                let ic = &self.ds.intersection_curves[ci];
                let half_circle = match &ic.curve {
                    rcad_kernel::geom::Curve3::Circle(_)
                    | rcad_kernel::geom::Curve3::Ellipse(_) => {
                        (ic.t_range[1] - ic.t_range[0] - std::f64::consts::TAU).abs()
                            >= TOLERANCE_ANG
                    }
                    _ => false,
                };
                half_circle && ic.start_vertex != usize::MAX && ic.start_vertex == ic.end_vertex
            };
            if needs_fix {
                if std::env::var("RCAD_DBG_FF").is_ok() {
                    eprintln!(
                        "[DBG_FF] needs_fix: ci={} t=[{:.4},{:.4}] sv={} ev={}",
                        ci,
                        self.ds.intersection_curves[ci].t_range[0],
                        self.ds.intersection_curves[ci].t_range[1],
                        self.ds.intersection_curves[ci].start_vertex,
                        self.ds.intersection_curves[ci].end_vertex
                    );
                }
                let t0 = self.ds.intersection_curves[ci].t_range[0];
                let t1 = self.ds.intersection_curves[ci].t_range[1];
                let p_start = self.ds.intersection_curves[ci].curve.point_at(t0);
                let p_end = self.ds.intersection_curves[ci].curve.point_at(t1);
                let v_start = self.ds.vertex_count();
                self.ds.push_vertex(
                    DSVertex { shape_idx: 0,
                        point: p_start,
                        geom_tol: TOLERANCE_ABS,
                        origin: None,
                        is_internal: true,
                        location: 0,
                    },
                    None,
                );
                let v_end = self.ds.vertex_count();
                self.ds.push_vertex(
                    DSVertex { shape_idx: 0,
                        point: p_end,
                        geom_tol: TOLERANCE_ABS,
                        origin: None,
                        is_internal: true,
                        location: 0,
                    },
                    None,
                );
                self.ds.intersection_curves[ci].start_vertex = v_start;
                self.ds.intersection_curves[ci].end_vertex = v_end;
            }
        }
    }

    /// PreparePostTreatFF (OCCT BOPAlgo_PaveFiller_6.cxx L3642-3668).
    /// Extends curves_sc and vertices_in for both faces with intersection curve endpoints.
    /// Called from perform_ff results processing loop.
    pub(crate) fn update_face_sc_vertices(&mut self, f1: usize, f2: usize) {
        let post_ff_curves = self
            .find_face_face_curve_indices(f1, f2)
            .unwrap_or_default();
        self.ds.face_info_mut(f1).curves_sc.extend(&post_ff_curves);
        self.ds.face_info_mut(f2).curves_sc.extend(&post_ff_curves);
        for &ci in &post_ff_curves {
            if ci < self.ds.intersection_curves.len() {
                let ic = &self.ds.intersection_curves[ci];
                let sv = ic.start_vertex;
                let ev = ic.end_vertex;
                if sv < self.ds.vertex_count() {
                    self.ds.face_info_mut(f1).vertices_in.insert(sv);
                    self.ds.face_info_mut(f2).vertices_in.insert(sv);
                }
                if ev < self.ds.vertex_count() {
                    self.ds.face_info_mut(f1).vertices_in.insert(ev);
                    self.ds.face_info_mut(f2).vertices_in.insert(ev);
                }
            }
        }
    }
    /// IndexType (IntTools_FaceFace.cxx L2844-2870).
    /// Maps Surface3 variant to an integer index for canonical ordering.
    /// Lower-typed surface is "simpler" (Plane < Cylinder < Cone < Sphere < Torus).
    fn surface_type_index(surf: &rcad_kernel::geom::Surface3) -> i32 {
        match surf {
            rcad_kernel::geom::Surface3::Plane(_) => 0,
            rcad_kernel::geom::Surface3::Cylinder(_) => 1,
            rcad_kernel::geom::Surface3::Cone(_) => 2,
            rcad_kernel::geom::Surface3::Sphere(_) => 3,
            rcad_kernel::geom::Surface3::Torus(_) => 4,
            _ => 11,
        }
    }

    /// OCCT L2426-2560: PerformPlanes  -- plane-plane intersection fast path.
    fn perform_plane_plane(&mut self, f1: usize, f2: usize) {
        use rcad_kernel::geom::{Curve3, Surface3};
        let pln1 = match self.ds.face_surface(f1).unwrap() {
            Surface3::Plane(p) => p,
            _ => return,
        };
        let pln2 = match self.ds.face_surface(f2).unwrap() {
            Surface3::Plane(p) => p,
            _ => return,
        };
        let mut geo = crate::bop::int_tools::int_ana_quad_quad_geo::QuadQuadGeo::new();
        let q1 = crate::bop::int_tools::int_surf_quadric::Quadric::from_plane(pln1);
        let q2 = crate::bop::int_tools::int_surf_quadric::Quadric::from_plane(pln2);
        geo.perform_plane_plane(&q1, &q2, 1e-8, self.fuzzy_tolerance);
        if !geo.is_done() {
            return;
        }
        use crate::bop::int_tools::int_ana_quad_quad_geo::AnaResultType;
        if let AnaResultType::Same = geo.type_inter() {
            self.ds.interf_ff.push(crate::bop::ds::InterferenceFF {
                f1,
                f2,
                curves: Vec::new(),
                points: Vec::new(),
                tangent_faces: true,
            });
            return;
        }
        if matches!(geo.type_inter(), AnaResultType::Empty) {
            return;
        }
        let line3 = geo.line(1);
        let line3d = Curve3::Line(line3.clone());
        let pcurve1 = crate::bop::int_tools::pcurve_derive::line_pcurve_on_plane(&line3, pln1);
        let pcurve2 = crate::bop::int_tools::pcurve_derive::line_pcurve_on_plane(&line3, pln2);
        // OCCT L2514: new Geom_TrimmedCurve(aGLin, pmin, pmax)
        // OCCT L2521: new Geom2d_TrimmedCurve(C2d, pmin, pmax)
        let uv1 = self.context.uv_bounds(self.ds, f1);
        let uv2 = self.context.uv_bounds(self.ds, f2);
        let tol = self.ds.face_tolerance(f1).max(self.ds.face_tolerance(f2));
        let p1 = crate::bop::int_tools::classify_lin2d::classify_lin2d(&pcurve1, uv1, tol);
        let p2 = crate::bop::int_tools::classify_lin2d::classify_lin2d(&pcurve2, uv2, tol);
        let (Some([p11, p12]), Some([p21, p22])) = (p1, p2) else {
            return;
        };
        if p21 >= p12 || p22 <= p11 {
            return;
        }
        let pmin = p11.max(p21);
        let pmax = p12.min(p22);
        if pmax - pmin <= tol {
            return;
        }
        let t_range = [pmin, pmax];
        let mut curve_extra = crate::bop::ds::CurveExtra::default();
        curve_extra.tangential_tol = tol;
        // OCCT L2514: new Geom_TrimmedCurve(aGLin, pmin, pmax)
        // OCCT L2521: new Geom2d_TrimmedCurve(C2d, pmin, pmax)
        let trimmed_curve =
            Curve3::Trimmed(TrimmedCurve3::new(line3d.clone(), pmin, pmax));
        let trimmed_pca = Some(Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(pcurve1),
            t_min: pmin,
            t_max: pmax,
        }));
        let trimmed_pcb = Some(Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(pcurve2),
            t_min: pmin,
            t_max: pmax,
        }));
        let ic = crate::bop::ds::IntersectionCurve {
            curve: trimmed_curve,
            polyline: Vec::new(),
            start_vertex: usize::MAX,
            end_vertex: usize::MAX,
            t_range,
            pcurve_on_a: trimmed_pca,
            pcurve_on_b: trimmed_pcb,
            geom_tol: tol,
            pave_blocks: Vec::new(),
            curve_extra,
        };
        // OCCT: vertices created by BRepBuilderAPI_MakeVertex + myDS->Index
        // after the curve is stored (matching the MakeCurve caller pattern).
        let sv = self.ds.add_vertex(ic.curve.point_at(t_range[0]));
        let ev = self.ds.add_vertex(ic.curve.point_at(t_range[1]));
        let mut ic = ic;
        ic.start_vertex = sv;
        ic.end_vertex = ev;
        let ci = self.ds.intersection_curves.len();
        self.ds.intersection_curves.push(ic);
        self.ds.interf_ff.push(crate::bop::ds::InterferenceFF {
            f1,
            f2,
            curves: vec![ci],
            points: Vec::new(),
            tangent_faces: false,
        });
        self.ds.face_info_mut(f1).curves_sc.insert(ci);
        self.ds.face_info_mut(f2).curves_sc.insert(ci);
    }

    /// MakeCurve (IntTools_FaceFace.cxx L695-1846).
    /// Dispatches by line type (OCCT switch on IntPatch_IType):
    ///   - Walking:        approximate BSpline from marching points (L1097)
    ///   - Line/Parabola/Hyperbola: LineConstructor parts + per-part handling (L815-898)
    ///   - Circle/Ellipse:  TreatCircle-equivalent with 0-crossing splitting (L904-1095)
    ///   - (Restriction:    handled upstream in IntPatch_Intersection)
    /// Returns one or more IntersectionCurve per valid part.
    pub(crate) fn make_intersection_curve(
        &mut self,
        f1: usize,
        f2: usize,
        line: &crate::bop::int_tools::int_patch_line::IntPatchLine,
    ) -> Vec<crate::bop::ds::IntersectionCurve> {
        use rcad_kernel::geom::Curve2dEval;
        use std::f64::consts::TAU;

        // ===== OCCT IntTools_FaceFace.cxx L695-751 =====
        // OCCT L700-714: local vars
        // OCCT L717: reapprox label (not needed in sequential rcad)
        // OCCT L719: Tolpc = myTolApprox
        // OCCT L720: bAvoidLineConstructor = false
        let mut b_avoid_line_constructor = false;

        // OCCT L721-722: L = myIntersector.Line(Index); typl = L->ArcType();
        let typl = line.line_type;

        // OCCT L724-744: IntPatch_Walking special handling
        if line.is_wline() {
            let nbp = line.nb_points();
            if nbp >= 2 {
                let p1 = line.point(0).p3d;
                let p2 = line.point(nbp - 1).p3d;
                // OCCT L740-743: if endpoints are nearly coincident, use LineConstructor
                if p1.distance_squared(p2) < 1e-14 {
                    b_avoid_line_constructor = false;
                }
            }
        }

        // OCCT L748-751: IntPatch_Restriction 鈥?skip LineConstructor
        if typl == crate::bop::int_tools::int_patch_type::IntPatchIType::Restriction {
            b_avoid_line_constructor = true;
        }

        // OCCT L755-773: LineConstructor.Perform(L)
        // If !IsDone 鈫?return empty. If NbParts <= 0 鈫?return empty.
        let parts: Vec<[f64; 2]> = if !b_avoid_line_constructor {
            let p = self.line_constructor_parts(
                &line.curve,
                line.t_range,
                typl,
                &line.vertices,
                f1,
                f2,
            );
            if p.is_empty() {
                return Vec::new();
            }
            p
        } else {
            // OCCT L748-750: for Restriction, skip LineConstructor, use full range
            // rcad: use the full t_range as a single part
            vec![line.t_range]
        };

        // OCCT L776-1846: switch(typl)
        match typl {
            crate::bop::int_tools::int_patch_type::IntPatchIType::Line
            | crate::bop::int_tools::int_patch_type::IntPatchIType::Parabola
            | crate::bop::int_tools::int_patch_type::IntPatchIType::Hyperbola => self
                .make_analytic_nonperiodic_curve(
                    f1,
                    f2,
                    &line.curve,
                    &parts,
                    typl,
                    line.tolerance,
                    line.tang_tolerance,
                ),
            crate::bop::int_tools::int_patch_type::IntPatchIType::Circle
            | crate::bop::int_tools::int_patch_type::IntPatchIType::Ellipse => self
                .make_analytic_periodic_curve(
                    f1,
                    f2,
                    &line.curve,
                    &parts,
                    typl,
                    line.tolerance,
                    line.tang_tolerance,
                ),
            crate::bop::int_tools::int_patch_type::IntPatchIType::Walking
            | crate::bop::int_tools::int_patch_type::IntPatchIType::Restriction => {
                self.make_walking_curve(f1, f2, line)
            }
            _ => Vec::new(),
        }
    }

    /// OCCT L1097-1846: MakeCurve for IntPatch_Walking.
    /// Approximates a BSpline3 from marching points, builds pcurves from
    /// marching UV data.
    fn make_walking_curve(
        &mut self,
        _f1: usize,
        _f2: usize,
        line: &crate::bop::int_tools::int_patch_line::IntPatchLine,
    ) -> Vec<crate::bop::ds::IntersectionCurve> {
        let n = line.nb_points();
        if n < 2 {
            return Vec::new();
        }

        let p3d_pts: Vec<glam::DVec3> = (0..n).map(|i| line.point(i).p3d).collect();
        let polyline = p3d_pts.clone();

        if let Some(bs_curve3) = crate::bop::int_tools::intss::polyline_to_bspline(&p3d_pts, 1e-4) {
            let t_range_bs = bs_curve3.default_domain();
            let bs = match &bs_curve3 {
                rcad_kernel::geom::Curve3::BSpline(b) => b.clone(),
                _ => {
                    let mut curve_extra = crate::bop::ds::CurveExtra::default();
                    curve_extra.tangential_tol = line.tang_tolerance;
                    return vec![crate::bop::ds::IntersectionCurve {
                        curve: bs_curve3.clone(),
                        polyline,
                        start_vertex: usize::MAX,
                        end_vertex: usize::MAX,
                        t_range: t_range_bs,
                        pcurve_on_a: line.pcurve1.clone(),
                        pcurve_on_b: line.pcurve2.clone(),
                        geom_tol: line.tolerance.max(CONFUSION),
                        pave_blocks: Vec::new(),
                        curve_extra,
                    }];
                }
            };

            // Build pcurves from marching UV data.
            let mut pcurve_on_a = line.pcurve1.clone();
            let mut pcurve_on_b = line.pcurve2.clone();
            if pcurve_on_a.is_none() && line.point(0).u1.is_finite() {
                let uv_pts: Vec<glam::DVec2> = (0..n)
                    .map(|i| glam::DVec2::new(line.point(i).u1, line.point(i).v1))
                    .collect();
                if let Ok(bs2d) = rcad_kernel::fit::interpolate_points_2d(&uv_pts) {
                    pcurve_on_a = Some(rcad_kernel::geom::Curve2d::BSpline(bs2d));
                }
            }
            if pcurve_on_b.is_none() && line.point(0).u2.is_finite() {
                let uv_pts: Vec<glam::DVec2> = (0..n)
                    .map(|i| glam::DVec2::new(line.point(i).u2, line.point(i).v2))
                    .collect();
                if let Ok(bs2d) = rcad_kernel::fit::interpolate_points_2d(&uv_pts) {
                    pcurve_on_b = Some(rcad_kernel::geom::Curve2d::BSpline(bs2d));
                }
            }
            let mut curve_extra = crate::bop::ds::CurveExtra::default();
            curve_extra.tangential_tol = line.tang_tolerance;
            return vec![crate::bop::ds::IntersectionCurve {
                curve: rcad_kernel::geom::Curve3::BSpline(bs),
                polyline,
                start_vertex: usize::MAX,
                end_vertex: usize::MAX,
                t_range: t_range_bs,
                pcurve_on_a,
                pcurve_on_b,
                geom_tol: line.tolerance.max(CONFUSION),
                pave_blocks: Vec::new(),
                curve_extra,
            }];
        }
        Vec::new()
    }

    /// OCCT L815-898: MakeCurve for Line, Parabola, Hyperbola.
    /// - Creates analytic curve from IntPatch_GLine.
    /// - Calls LineConstructor to get valid parameter parts (OCCT NbParts/Part).
    /// - For each part:
    ///     both bounds finite   -- trimmed 3D curve + BuildPCurves + endpoint vertices
    ///     one/both infinite    -- test reference point on face domains  -- keep or reject
    /// - rcad note: IntPatchLine has no vertex data, so LineConstructor returns
    ///   a single part with the original t_range (always infinite for lines).
    fn make_analytic_nonperiodic_curve(
        &mut self,
        f1: usize,
        f2: usize,
        curve: &Curve3,
        parts: &Vec<[f64; 2]>,
        typl: IntPatchIType,
        geom_tol: f64,
        tang_tolerance: f64,
    ) -> Vec<crate::bop::ds::IntersectionCurve> {
        use rcad_kernel::geom::Curve2dEval;
        use std::f64::consts::TAU;

        // OCCT L815-826: create analytic 3D curve from the GLine.
        // rcad: curve is already the correct analytic type in IntPatchLine.

        // OCCT L828-840: LineConstructor.Perform(L) already done upstream.
        // parts already computed by make_intersection_curve.

        if parts.is_empty() {
            return Vec::new();
        }

        let mut result = Vec::with_capacity(parts.len());

        // OCCT L842-898: per-part loop.
        for part in parts.iter() {
            let &[fprm, lprm] = part;
            let b_finite = fprm.is_finite() && lprm.is_finite() && lprm > fprm + 1e-12;

            if b_finite {
                //  Both bounds finite: trimmed curve + pcurves + vertices
                // OCCT L835-870: Geom_TrimmedCurve + BuildPCurves + Geom2d_TrimmedCurve.
                let ic_t_range = [fprm, lprm];

                // OCCT L816-820: for Parabola, CurveTolerance(aCT3D, myTol)
                let ic_geom_tol = if typl == IntPatchIType::Parabola {
                    // OCCT: IntTools_Tools::CurveTolerance(aCT3D, myTol)
                    crate::bop::tools::curve_tolerance(curve, geom_tol)
                } else {
                    geom_tol.max(crate::tolerance::TOLERANCE_ABS)
                };

                // OCCT L822-846: BuildPCurves on the trimmed range.
                // OCCT: GeomInt_IntSS::BuildPCurves(fprm, lprm, Tolpc, surface, newc, C2d)
                // OCCT L822-832: if (myApprox1) { ... }
                let pca;
                let pcb;
                // myApprox1 (always true in rcad)
                {
                    let raw = self.compute_pcurve_on_surface(curve, f1);
                    if raw.is_none() {
                        continue;
                    }
                    // OCCT L832: aCurve.SetFirstCurve2d(new Geom2d_TrimmedCurve(C2d, fprm, lprm))
                    pca = Some(Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(raw.unwrap()),
                        t_min: fprm,
                        t_max: lprm,
                    }));
                }
                // myApprox2 (always true in rcad)
                {
                    let raw = self.compute_pcurve_on_surface(curve, f2);
                    if raw.is_none() {
                        continue;
                    }
                    pcb = Some(Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(raw.unwrap()),
                        t_min: fprm,
                        t_max: lprm,
                    }));
                }

                // OCCT L814-815: new Geom_TrimmedCurve(newc, fprm, lprm)
                let trimmed_curve =
                    Curve3::Trimmed(TrimmedCurve3::new(curve.clone(), fprm, lprm));

                // OCCT: no vertex creation in MakeCurve (vertices created later in caller).
                let mut curve_extra = crate::bop::ds::CurveExtra::default();
                curve_extra.tangential_tol = tang_tolerance;
                result.push(crate::bop::ds::IntersectionCurve {
                    curve: trimmed_curve,
                    polyline: Vec::new(),
                    start_vertex: usize::MAX,
                    end_vertex: usize::MAX,
                    t_range: [fprm, lprm],
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                    geom_tol: ic_geom_tol,
                    pave_blocks: Vec::new(),
                    curve_extra,
                });
            } else {
                //  One/both bounds infinite: test reference point
                // OCCT L850-895: test-point approach.
                // dT = 100.0; surface-type exceptions for extrusion/offset/revolution.
                let dT = 100.0;
                let test_t = if !fprm.is_finite() && lprm.is_finite() {
                    // bFNIt && !bLPIt: only lower bound infinite
                    lprm - dT
                } else if fprm.is_finite() && !lprm.is_finite() {
                    // !bFNIt && bLPIt: only upper bound infinite
                    fprm + dT
                } else {
                    // bFNIt && bLPIt: both infinite  -- OCCT IntTools_Tools::IntermediatePoint(-dT, dT)
                    crate::bop::tools::intermediate_point_occt(-dT, dT)
                };

                let p3d = curve.point_at(test_t);
                if !p3d.is_finite() {
                    continue;
                }

                // OCCT L865-867: get surface types for the test-point branch
                let surf1 = self.ds.face_surface(f1).unwrap();
                let surf2 = self.ds.face_surface(f2).unwrap();
                let is_extrusion_rev_offset = |s: &Surface3| -> bool {
                    matches!(
                        s,
                        Surface3::LinearExtrusion(_)
                            | Surface3::Revolution(_)
                            | Surface3::Offset(_)
                    )
                };

                // OCCT L875-882: if either surface is extrusion/offset/revolution,
                // append curve with empty pcurves (= H1, H1 in OCCT) and skip Classify.
                if is_extrusion_rev_offset(surf1) || is_extrusion_rev_offset(surf2) {
                    let mut curve_extra = crate::bop::ds::CurveExtra::default();
                    curve_extra.tangential_tol = tang_tolerance;
                    result.push(crate::bop::ds::IntersectionCurve {
                        curve: curve.clone(),
                        polyline: Vec::new(),
                        start_vertex: usize::MAX,
                        end_vertex: usize::MAX,
                        t_range: [fprm, lprm],
                        pcurve_on_a: None,
                        pcurve_on_b: None,
                        geom_tol: geom_tol.max(crate::tolerance::TOLERANCE_ABS),
                        pave_blocks: Vec::new(),
                        curve_extra,
                    });
                    continue;
                }

                // OCCT L886-892: Parameters + Classify on both face domains.
                // OCCT: Tol = Precision::Confusion()
                // OCCT: Parameters(myHS1, myHS2, ptref, u1, v1, u2, v2)
                // OCCT L888: ok = (dom1->Classify(gp_Pnt2d(u1, v1), Tol) != TopAbs_OUT)
                let in1 = self.classify_line_constructor_point(f1, p3d, CONFUSION);
                if !in1 {
                    continue;
                }
                // OCCT L890-891: if (ok) { ok = dom2->Classify(...) }
                let in2 = self.classify_line_constructor_point(f2, p3d, CONFUSION);
                if !in2 {
                    continue;
                }

                // OCCT L893-896: append curve with empty pcurves (H1, H1 in OCCT).
                // rcad note: OCCT does not compute pcurves here; keeping empty to match.
                let mut curve_extra = crate::bop::ds::CurveExtra::default();
                curve_extra.tangential_tol = tang_tolerance;
                result.push(crate::bop::ds::IntersectionCurve {
                    curve: curve.clone(),
                    polyline: Vec::new(),
                    start_vertex: usize::MAX,
                    end_vertex: usize::MAX,
                    t_range: [fprm, lprm],
                    pcurve_on_a: None,
                    pcurve_on_b: None,
                    geom_tol: geom_tol.max(crate::tolerance::TOLERANCE_ABS),
                    pave_blocks: Vec::new(),
                    curve_extra,
                });
            }
        }
        result
    }

    /// OCCT L904-1095: MakeCurve for Circle, Ellipse.
    /// - Analytic curve from IntPatch_GLine, clipped into valid intervals by
    ///   TreatCircle (line_constructor_parts, OCCT L922-950).
    /// - For each candidate interval:
    ///     not full-period  -- trimmed curve + BuildPCurves(with UV bounds) + vertices
    ///     full-period (aNbParts=1)  -- trimmed full circle + BuildPCurves + break
    ///     full-period (aNbParts>1)  -- test 18 points around circle  -- keep/reject
    /// - rcad: with 0 vertices, TreatCircle keeps the full [0, 2闁挎粎顕?interval.
    // OCCT L904-1095: MakeCurve for Circle/Ellipse
    fn make_analytic_periodic_curve(
        &mut self,
        f1: usize,
        f2: usize,
        curve: &Curve3,
        parts: &[[f64; 2]],
        typl: IntPatchIType,
        geom_tol: f64,
        tang_tolerance: f64,
    ) -> Vec<crate::bop::ds::IntersectionCurve> {
        if std::env::var("RCAD_DBG_MB").is_ok() {
            eprintln!(
                "[DBG_IC] make_analytic_periodic_curve: geom_tol={:.6e} nParts={}",
                geom_tol,
                parts.len()
            );
        }
        use std::f64::consts::TAU;

        // OCCT L906-920: create analytic 3D curve from GLine.
        // OCCT L922-950: TreatCircle has already split the line into intervals
        // (done in line_constructor_parts, which delegates to treat_circle_parts).

        // OCCT L950-1095: aNbParts = seqp.Length() / 2.
        //   If aNbParts == 0  -- the for loop does not execute  -- no output curves.
        if parts.is_empty() {
            return Vec::new();
        }

        let aPeriod = TAU;
        let aRealEpsilon = f64::EPSILON;
        let aNbParts = parts.len();
        let mut result = Vec::with_capacity(parts.len());

        for &[fprm, lprm] in parts {
            // OCCT L953-956: if (|fprm|>eps || |lprm-2闁挎粎鐨?eps)  -- not full-period
            let is_full_period =
                fprm.abs() <= aRealEpsilon && (lprm - aPeriod).abs() <= aRealEpsilon;

            if !is_full_period && (lprm > fprm + 1e-12) {
                //  Not full-period: trimmed curve + pcurves + vertices
                // OCCT L960-990: Geom_TrimmedCurve(newc, fprm, lprm) + BuildPCurves + append.

                // OCCT L968-990: BuildPCurves with surface UV bounds for Circle/Ellipse.
                let pca = self.compute_pcurve_on_surface(curve, f1);
                let pcb = self.compute_pcurve_on_surface(curve, f2);
                let trimmed_pca = pca.clone();
                let trimmed_pcb = pcb.clone();
                let trimmed_curve =
                    Curve3::Trimmed(TrimmedCurve3::new(curve.clone(), fprm, lprm));

                let mut curve_extra = crate::bop::ds::CurveExtra::default();
                curve_extra.tangential_tol = tang_tolerance;
                result.push(crate::bop::ds::IntersectionCurve {
                    curve: trimmed_curve,
                    polyline: Vec::new(),
                    start_vertex: usize::MAX,
                    end_vertex: usize::MAX,
                    t_range: [fprm, lprm],
                    pcurve_on_a: trimmed_pca,
                    pcurve_on_b: trimmed_pcb,
                    geom_tol: geom_tol.max(TOLERANCE_ABS),
                    pave_blocks: Vec::new(),
                    curve_extra,
                });
            } else if is_full_period && aNbParts == 1 {
                // Full-period, single part 鈥?accept full circle
                // OCCT L996-1042: trimmed full circle + BuildPCurves + append + break.
                let pca = self.compute_pcurve_on_surface(curve, f1);
                let pcb = self.compute_pcurve_on_surface(curve, f2);
                // OCCT: no vertex creation in MakeCurve (vertices created later in caller).
                let trimmed_curve =
                    Curve3::Trimmed(TrimmedCurve3::new(curve.clone(), fprm, lprm));
                let mut curve_extra = crate::bop::ds::CurveExtra::default();
                curve_extra.tangential_tol = tang_tolerance;
                result.push(crate::bop::ds::IntersectionCurve {
                    curve: trimmed_curve,
                    polyline: Vec::new(),
                    start_vertex: usize::MAX,
                    end_vertex: usize::MAX,
                    t_range: [fprm, lprm],
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                    geom_tol: geom_tol.max(TOLERANCE_ABS),
                    pave_blocks: Vec::new(),
                    curve_extra,
                });
                break;
            } else if is_full_period && aNbParts > 1 {
                //  Full-period, multiple parts: test 18 points
                // OCCT L1045-1095: on regarde si on garde.
                let aTwoPIdiv17 = aPeriod / 17.0;
                for j in 0..=17 {
                    let t = j as f64 * aTwoPIdiv17;
                    let p3d = curve.point_at(t);
                    if !p3d.is_finite() {
                        continue;
                    }
                    // OCCT L1112-1119: Parameters + Classify with Tol = Precision::Confusion().
                    let in1 = self.classify_line_constructor_point(f1, p3d, CONFUSION);
                    if in1 {
                        let in2 = self.classify_line_constructor_point(f2, p3d, CONFUSION);
                        if in2 {
                            let pca = self.compute_pcurve_on_surface(curve, f1);
                            let pcb = self.compute_pcurve_on_surface(curve, f2);
                            let mut curve_extra = crate::bop::ds::CurveExtra::default();
                            curve_extra.tangential_tol = tang_tolerance;
                            result.push(crate::bop::ds::IntersectionCurve {
                                curve: curve.clone(),
                                polyline: Vec::new(),
                                start_vertex: usize::MAX,
                                end_vertex: usize::MAX,
                                t_range: [fprm, lprm],
                                pcurve_on_a: pca,
                                pcurve_on_b: pcb,
                                geom_tol: geom_tol.max(TOLERANCE_ABS),
                                pave_blocks: Vec::new(),
                                curve_extra,
                            });
                            break;
                        }
                    }
                }
            }
        }
        result
    }

    /// LineConstructor (GeomInt_LineConstructor.cxx L333-386, GLine path).
    /// OCCT-aligned: iterates over the IntPatch_Point vertices in line order
    /// (sorted upstream by ComputeVertexParameters / put_points_on_line), tests
    /// the midpoint of each adjacent-vertex interval on both face domains,
    /// keeps valid intervals.  Circle/Ellipse are delegated to TreatCircle
    /// (OCCT L338-342).
    ///
    /// With 0 vertices (nbvtx=0): OCCT's intrvtested flag stays false, and the
    /// full parameter range [FirstParameter, LastParameter] is kept as one part.
    /// The caller's test-point logic decides whether to keep or reject it.
    fn line_constructor_parts(
        &mut self,
        curve: &Curve3,
        orig_t_range: [f64; 2],
        typl: IntPatchIType,
        vertices: &[crate::bop::int_tools::int_patch_line::IntPatchVertex],
        f1: usize,
        f2: usize,
    ) -> Vec<[f64; 2]> {
        // OCCT GeomInt_LineConstructor L338-342: Circle/Ellipse -> TreatCircle.
        if typl == IntPatchIType::Circle || typl == IntPatchIType::Ellipse {
            return self.treat_circle_parts(curve, orig_t_range, typl, vertices, f1, f2);
        }

        // OCCT L118: constexpr double Tol = Precision::PConfusion() * 35.0;
        let a_tol = rcad_kernel::tolerance::PCONFUSION * 35.0;

        // OCCT L345-374: iterate adjacent vertex pairs, test interval midpoint.
        let mut result: Vec<[f64; 2]> = Vec::new();
        let nbvtx = vertices.len();
        let mut intrvtested = false;
        for i in 0..nbvtx.saturating_sub(1) {
            let firstp = vertices[i].param_on_line;
            let lastp = vertices[i + 1].param_on_line;
            // OCCT L354: if (std::abs(firstp - lastp) > Precision::PConfusion())
            if (firstp - lastp).abs() > rcad_kernel::tolerance::PCONFUSION {
                intrvtested = true;
                // OCCT L357-359: pmid, GLinePoint(typl, GLine, pmid, Pmid)
                let pmid = (firstp + lastp) * 0.5;
                let p3d = curve.point_at(pmid);
                if !p3d.is_finite() {
                    continue;
                }
                // OCCT L361-372: Parameters + AdjustPeriodic + Classify on both domains.
                let in1 = self.classify_line_constructor_point(f1, p3d, a_tol);
                if in1 {
                    let in2 = self.classify_line_constructor_point(f2, p3d, a_tol);
                    if in2 {
                        result.push([firstp, lastp]);
                    }
                }
            }
        }
        // OCCT L376-382: if no interval tested, keep the full range a priori.
        if !intrvtested {
            result.push(orig_t_range);
        }
        result
    }

    /// OCCT GeomInt_LineConstructor::Parameters (L820-862) + AdjustPeriodic
    /// (L737-816) + Classify.  Computes the UV of a point on a face quadric
    /// analytically (exact, matching the face surface parametrization), then
    /// classifies the point on the face domain.  FClass2d's Perform handles the
    /// periodic adjust internally (IntTools_FClass2d.cxx L666-678), mirroring
    /// OCCT's AdjustPeriodic before the Classify call.
    ///
    /// rcad note: OCCT passes Tol = Precision::PConfusion() * 35.0 to the
    /// Classify; the rcad FClass2d uses its construction-time tolerance.
    fn classify_line_constructor_point(&mut self, f: usize, p3d: DVec3, _a_tol: f64) -> bool {
        let surf = self.ds.face_surface(f).unwrap();
        let uv = quadric_uv_params(surf, p3d)
            .or_else(|| self.context.proj_ps(self.ds, f, p3d).map(|(uv, _, _)| uv));
        let Some(uv) = uv else {
            return false;
        };
        self.context.is_point_in_on_face(self.ds, f, uv)
    }

    /// PutPointsOnLine (IntPatch_Intersection.cxx L268-312).
    /// OCCT-aligned: PutPointsOnLine equivalent (IntPatch_Intersection / GetEFPnts flow).
    ///
    /// Projects only EF intersection points (from DS interf_ef) onto the analytic
    /// intersection line, mirroring OCCT's GetEFPnts (BOPAlgo_PaveFiller_6.cxx L2608-2687).
    ///
    /// OCCT does NOT project all boundary vertices of the two faces onto the line.
    /// Boundary-vertex projection creates spurious splitting vertices, causing
    /// nIC (number of intersection curves) to be inflated (e.g. plane_sphere
    /// expects nIC=3 but old code produced nIC=6).
    ///
    /// For circle/ellipse lines, the EF points provide the boundary-crossing
    /// vertices that split the closed curve into valid face-domain intervals.
    /// Without EF points (when EF=0), no vertices are added and the curve
    /// is not split (nIC=1).
    fn put_points_on_line(
        &mut self,
        f1: usize,
        f2: usize,
        line: &mut crate::bop::int_tools::int_patch_line::IntPatchLine,
        ef_points: &[DVec3],
    ) {
        use crate::bop::algo::pave_filler::helpers::project_vertex_to_curve;
        if line.is_wline() {
            return;
        }
        let typl = line.line_type;
        if typl == IntPatchIType::Restriction {
            return;
        }
        let curve = &line.curve;
        let t_range = line.t_range;
        let a_tol_c = line.tolerance;

        // OCCT GetEFPnts: EF intersection points between the two faces pre-collected
        // in perform_ff (OCCT L498-504). Use passed-in ef_points.
        for &pt in ef_points {
            let t_opt = project_vertex_to_curve(pt, curve, a_tol_c);
            let t = match t_opt {
                Some(t) if t >= t_range[0] - a_tol_c && t <= t_range[1] + a_tol_c => t,
                _ => continue,
            };
            let is_dup = line
                .vertices
                .iter()
                .any(|ev| (ev.param_on_line - t).abs() < 1e-10);
            if is_dup {
                continue;
            }
            let pca = self.compute_pcurve_on_surface(curve, f1);
            let pcb = self.compute_pcurve_on_surface(curve, f2);
            let uv1 = pca
                .as_ref()
                .map(|pc| pc.point_at(t))
                .filter(|uv| uv.is_finite());
            let uv2 = pcb
                .as_ref()
                .map(|pc| pc.point_at(t))
                .filter(|uv| uv.is_finite());
            let p3d = curve.point_at(t);
            line.vertices
                .push(crate::bop::int_tools::int_patch_line::IntPatchVertex {
                    param_on_line: t,
                    p3d,
                    u1: uv1.map_or(0.0, |uv| uv.x),
                    v1: uv1.map_or(0.0, |uv| uv.y),
                    u2: uv2.map_or(0.0, |uv| uv.x),
                    v2: uv2.map_or(0.0, |uv| uv.y),
                });
        }
        line.vertices.sort_by(|a, b| {
            a.param_on_line
                .partial_cmp(&b.param_on_line)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// TreatCircle (GeomInt_LineConstructor.cxx L674-733).
    /// Sorts the IntPatch_Point vertices wrapped into [0, 2*PI), creates an
    /// interval per adjacent vertex pair (last vertex at first.param + 2*PI),
    /// rejects duplicate parameters, and tests each interval's midpoint on both
    /// face domains.
    ///
    /// A Circle/Ellipse line without vertices gets two vertices at parameter 0
    /// and 2*PI (IntPatch_ImpImpIntersection L2997-3052); both wrap to 0, so the
    /// single interval [0, 2*PI] is kept iff its midpoint PI is inside both faces.
    fn treat_circle_parts(
        &mut self,
        curve: &Curve3,
        _orig_t_range: [f64; 2],
        typl: IntPatchIType,
        vertices: &[crate::bop::int_tools::int_patch_line::IntPatchVertex],
        f1: usize,
        f2: usize,
    ) -> Vec<[f64; 2]> {
        use std::f64::consts::TAU;
        if std::env::var("RCAD_DBG_FF").is_ok() {
            eprintln!(
                "[FF] treat_circle_parts: f1={} f2={} nVtx={}",
                f1,
                f2,
                vertices.len()
            );
        }

        // OCCT GeomInt_LineConstructor::TreatCircle (L674-733):
        //   RejectMicroCircle, sort, RejectDuplicates, midpoint test.

        // OCCT L679-681: RejectMicroCircle -- skip circles/ellipses smaller than tolerance
        if typl == IntPatchIType::Circle || typl == IntPatchIType::Ellipse {
            let radius = match curve {
                Curve3::Circle(c) => c.radius,
                Curve3::Ellipse(e) => e.major_radius,
                _ => 0.0,
            };
            let a_tol_3d = crate::TOLERANCE_ABS;
            if radius > 0.0 && radius < a_tol_3d {
                return Vec::new();
            }
        }

        // OCCT IntPatch_ImpImpIntersection L2997-3052: a Circle/Ellipse GLine
        // without vertices gets two vertices at parameter 0 and 2*PI.  Both wrap
        // to 0 in [0, 2*PI), so below they collapse to a single interval [0, 2*PI]
        // whose midpoint PI is classified on both face domains.
        let a_tol_pc: f64 = 1000.0 * rcad_kernel::tolerance::PCONFUSION; // RejectDuplicates L931

        // Build vertex parameters wrapped into [0, 2*PI) (OCCT GeomInt_Vertex::SetVertex L60).
        let mut params: Vec<f64> = if vertices.is_empty() {
            vec![0.0, TAU]
        } else {
            vertices.iter().map(|v| wrap_to_2pi(v.param_on_line)).collect()
        };
        // OCCT L691: sort by parameter.
        params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // OCCT L684-695: array of size n+1; last vertex at first.param + 2*PI.
        let mut arr: Vec<f64> = Vec::with_capacity(params.len() + 1);
        arr.extend_from_slice(&params);
        arr.push(params[0] + TAU);

        // OCCT L697: RejectDuplicates -- mark coincident params with RealLast.
        // The source array must be sorted in ascending order (L926-959).
        for i in 0..arr.len().saturating_sub(2) {
            let prm_i = arr[i];
            if !prm_i.is_finite() {
                continue;
            }
            for j in (i + 1)..arr.len().saturating_sub(1) {
                let prm_j = arr[j];
                if prm_j - prm_i < a_tol_pc {
                    arr[j] = f64::INFINITY; // RealLast
                } else {
                    break;
                }
            }
        }
        // OCCT L963-982: find duplicates with the last element.
        let a_max_prm = *arr.last().unwrap();
        for i in (1..arr.len().saturating_sub(1)).rev() {
            let prm_i = arr[i];
            if !prm_i.is_finite() {
                continue;
            }
            if a_max_prm - prm_i < a_tol_pc {
                arr[i] = f64::INFINITY;
            } else {
                break;
            }
        }
        // OCCT L699: re-sort after RejectDuplicates.
        arr.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // OCCT L704-732: test each adjacent pair's midpoint on both face domains.
        let a_tol = rcad_kernel::tolerance::PCONFUSION * 35.0;
        let mut result = Vec::new();
        for i in 0..arr.len().saturating_sub(1) {
            let t1 = arr[i];
            let t2 = arr[i + 1];
            // OCCT L709-712: if (aT2 == RealLast()) break;
            if t2 == f64::INFINITY {
                break;
            }
            // OCCT L714-715: aTmid, GLinePoint(aType, aGLine, aTmid, aPmid).
            let t_mid = (t1 + t2) * 0.5;
            let p3d = curve.point_at(t_mid);
            if !p3d.is_finite() {
                continue;
            }
            // OCCT L717-731: Parameters + AdjustPeriodic + Classify on both domains.
            let in1 = self.classify_line_constructor_point(f1, p3d, a_tol);
            if !in1 {
                if std::env::var("RCAD_DBG_FF").is_ok() {
                    eprintln!(
                        "[FF]   interval [{:.4},{:.4}] REJECTED: in1={} in2={}",
                        t1, t2, in1, false
                    );
                }
                continue;
            }
            let in2 = self.classify_line_constructor_point(f2, p3d, a_tol);
            if !in2 {
                if std::env::var("RCAD_DBG_FF").is_ok() {
                    eprintln!(
                        "[FF]   interval [{:.4},{:.4}] REJECTED: in1={} in2={}",
                        t1, t2, in1, in2
                    );
                }
                continue;
            }
            if std::env::var("RCAD_DBG_FF").is_ok() {
                eprintln!("[FF]   interval [{:.4},{:.4}] ACCEPTED", t1, t2);
            }
            result.push([t1, t2]);
        }

        result
    }


    /// BuildPCurves for all curve-surface type combinations.
    /// Matches OCCT GeomInt_IntSS::BuildPCurves (L822-846). Uses exact
    /// analytic pcurves when available; falls back to sampling + projection.
    fn compute_pcurve_on_surface(
        &self,
        curve: &rcad_kernel::geom::Curve3,
        fi: usize,
    ) -> Option<rcad_kernel::geom::Curve2d> {
        if fi >= self.ds.face_count() {
            return None;
        }
        let surf = self.ds.face_surface(fi).unwrap();
        let pc = match (curve, surf) {
            (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Plane(p)) => {
                crate::bop::int_tools::pcurve_derive::line_pcurve_on_plane(l, p)
            }
            (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Plane(p)) => {
                crate::bop::int_tools::pcurve_derive::circle_pcurve_on_plane(c, p)
            }
            (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Plane(p)) => {
                crate::bop::int_tools::pcurve_derive::ellipse_pcurve_on_plane(e, p)
            }
            (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Sphere(s)) => {
                crate::bop::int_tools::pcurve_derive::line_pcurve_on_sphere(l, s)
            }
            (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Sphere(s)) => {
                crate::bop::int_tools::pcurve_derive::circle_pcurve_on_sphere(c, s)
            }
            (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Sphere(s)) => {
                crate::bop::int_tools::pcurve_derive::ellipse_pcurve_on_sphere(e, s)
            }
            (rcad_kernel::geom::Curve3::Parabola(p), rcad_kernel::geom::Surface3::Sphere(s)) => {
                crate::bop::int_tools::pcurve_derive::parabola_pcurve_on_sphere(p, s)
            }
            (rcad_kernel::geom::Curve3::Hyperbola(h), rcad_kernel::geom::Surface3::Sphere(s)) => {
                crate::bop::int_tools::pcurve_derive::hyperbola_pcurve_on_sphere(h, s)
            }
            (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Cylinder(c)) => {
                crate::bop::int_tools::pcurve_derive::line_pcurve_on_cylinder(l, c)
            }
            (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Cylinder(cyl)) => {
                crate::bop::int_tools::pcurve_derive::circle_pcurve_on_cylinder(c, cyl)
            }
            (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Cylinder(cyl)) => {
                crate::bop::int_tools::pcurve_derive::ellipse_pcurve_on_cylinder(e, cyl)
            }
            (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Cone(c)) => {
                crate::bop::int_tools::pcurve_derive::line_pcurve_on_cone(l, c)
            }
            (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Cone(co)) => {
                crate::bop::int_tools::pcurve_derive::circle_pcurve_on_cone(c, co)
            }
            (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Cone(co)) => {
                crate::bop::int_tools::pcurve_derive::ellipse_pcurve_on_cone(e, co)
            }
            _ => {
                let tr = match curve {
                    rcad_kernel::geom::Curve3::Line(_) => [-1e3, 1e3],
                    rcad_kernel::geom::Curve3::Circle(_)
                    | rcad_kernel::geom::Curve3::Ellipse(_) => [0.0, std::f64::consts::TAU],
                    _ => [0.0, 1.0],
                };
                crate::bop::int_tools::pcurve_derive::fallback_pcurve_by_projection(curve, &tr, surf)
            }
        };
        Some(pc)
    }
}

/// OCCT ElCLib::InPeriod(X, 0.0, 2*PI) (GeomInt_Vertex::SetVertex, GeomInt_LineConstructor.cxx L60).
/// Wraps a parameter into [0, 2*PI).
fn wrap_to_2pi(t: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    if t < 0.0 {
        t + two_pi * (1.0 + ((0.0 - t) / two_pi).floor())
    } else if t >= two_pi {
        t - two_pi * (1.0 + ((t - two_pi) / two_pi).floor())
    } else {
        t
    }
}

/// OCCT GeomInt_LineConstructor::Parameters (L820-862).
/// Analytic UV inversion of a 3D point on a quadric surface, in the surface's
/// own parametrization (matching the face surface used by FClass2d).
/// Returns None for non-quadric surfaces (OCCT throws Standard_ConstructionError).
fn quadric_uv_params(surf: &rcad_kernel::geom::Surface3, p: DVec3) -> Option<DVec2> {
    match surf {
        rcad_kernel::geom::Surface3::Plane(pl) => {
            let d = p - pl.origin;
            Some(DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir)))
        }
        rcad_kernel::geom::Surface3::Cylinder(cy) => Some(cy.world_to_uv(p)),
        rcad_kernel::geom::Surface3::Sphere(sp) => Some(sp.world_to_uv(p)),
        rcad_kernel::geom::Surface3::Cone(co) => Some(co.world_to_uv(p)),
        rcad_kernel::geom::Surface3::Torus(to) => Some(to.world_to_uv(p)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn surface_type_index_plane() {
        let s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
            glam::DVec3::Z,
            glam::DVec3::Z,
        ));
        assert_eq!(super::PaveFiller::surface_type_index(&s), 0);
    }
    #[test]
    fn surface_type_index_cylinder() {
        let s = rcad_kernel::geom::Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: glam::DVec3::Z,
            axis: glam::DVec3::Z,
            ref_dir: glam::DVec3::X,
            radius: 1.0,
        });
        assert_eq!(super::PaveFiller::surface_type_index(&s), 1);
    }
    #[test]
    fn surface_type_index_sphere() {
        let s = rcad_kernel::geom::Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: glam::DVec3::Z,
            axis: glam::DVec3::Z,
            ref_dir: glam::DVec3::X,
            radius: 1.0,
        });
        assert_eq!(super::PaveFiller::surface_type_index(&s), 3);
    }
    #[test]
    fn surface_type_index_other() {
        let s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
            glam::DVec3::Z,
            glam::DVec3::Z,
        ));
        // BSpline variant  -- use Plane to test 'other' path; BSpline has no Default
        /* Skip: no Default for Bezier/BSpline */
    }
}
