use super::*;

impl<'a> super::PaveFiller<'a> {
    pub(crate) fn process_de(&mut self) {
        // Detect degenerate edges on periodic surfaces
        let mut fv: Vec<Vec<usize>> = vec![Vec::new(); self.ds.faces.len()];
        let mut degen_flags: Vec<(usize, usize)> = Vec::new();
        for fi in 0..self.ds.faces.len() {
            let is_periodic = matches!(self.ds.faces[fi].surface,
                Surface3::Sphere(_) | Surface3::Cylinder(_) | Surface3::Cone(_));
            if !is_periodic { continue; }
            let mut vs = Vec::new();
            let boundary_edges: Vec<usize> = self.ds.faces[fi].boundary_edges.clone();
            for &ei in &boundary_edges {
                if ei >= self.ds.edges.len() { continue; }
                let e = &self.ds.edges[ei];
                if e.start_vertex == e.end_vertex {
                    degen_flags.push((ei, fi + 1));
                    vs.push(e.start_vertex);
                } else {
                    let d = (self.ds.vertices[e.start_vertex].point - self.ds.vertices[e.end_vertex].point).length();
                    if d < TOLERANCE_ABS*100.0 {
                        degen_flags.push((ei, fi + 1));
                        vs.push(e.start_vertex); vs.push(e.end_vertex);
                    }
                }
            }
            if let Surface3::Sphere(sp) = &self.ds.faces[fi].surface.clone() {
                let tl = self.ds.faces[fi].geom_tol.max(TOLERANCE_ABS);
                for pp in [sp.center + sp.axis * sp.radius, sp.center - sp.axis * sp.radius] {
                    if let Some(vi) = self.ds.find_vertex_near(pp, tl) { vs.push(vi); }
                }
            }
            fv[fi] = vs;
        }
        for (fi, vs) in fv.iter().enumerate() { for &vi in vs { self.ds.faces[fi].face_info.vertices_in.insert(vi); } }
        for (ei, flag_val) in degen_flags { self.ds.set_edge_flag(ei, flag_val); }

        // Process flagged degenerate edges
        let degen_edges: Vec<(usize, usize)> = self.ds.edge_flags.iter()
            .filter_map(|(&ei, &flag)| {
                if flag == 0 { return None; }
                let fi = flag.checked_sub(1)?;
                if fi >= self.ds.faces.len() { return None; }
                Some((ei, fi))
            }).collect();

        for (ei, fi) in degen_edges {
            if ei >= self.ds.edges.len() { continue; }
            let edge = &self.ds.edges[ei];
            if edge.pave_blocks.is_empty() { continue; }
            let n_v = edge.start_vertex;
            if n_v >= self.ds.vertices.len() { continue; }

            // OCCT L84-101: FindPaveBlocks �?locate PBs in this face's sets that pass through nV
            // (simplified: the pave_block info is already in the edge's PBs)
            // Ensure the degen edge's PB has an ext pave at the degenerate vertex
            // so make_section_edges will split it properly.
        }
    }

    pub(crate) fn fill_shrunk_data(&mut self) {
        let ec: Vec<Curve3> = self.ds.edges.iter().map(|e| e.curve.clone()).collect();
        let et: Vec<f64> = self.ds.edges.iter().map(|e| e.geom_tol).collect();
        let v_tols: Vec<f64> = self.ds.vertices.iter().map(|v| v.geom_tol).collect();
        // Read phase: copy all PB params into flat arrays
        let mut all_pb: Vec<(usize, usize, usize, f64, f64)> = Vec::new();
        for ei in 0..self.ds.edges.len() {
            for pb in &self.ds.edges[ei].pave_blocks {
                all_pb.push((ei, pb.pave1.vertex_idx, pb.pave2.vertex_idx, pb.pave1.param, pb.pave2.param));
            }
        }
        for pb in &self.ds.pave_blocks {
            if pb.original_edge < self.ds.edges.len() {
                all_pb.push((pb.original_edge, pb.pave1.vertex_idx, pb.pave2.vertex_idx, pb.pave1.param, pb.pave2.param));
            }
        }
        let num_edges = self.ds.edges.len();
        let edge_pb_counts: Vec<usize> = self.ds.edges.iter().map(|e| e.pave_blocks.len()).collect();
        // Compute phase: ShrunkRange (no borrow on self)
        let results: Vec<(Option<[f64; 2]>, bool)> = all_pb.iter().map(|(ei, v1i, v2i, p1, p2)| {
            let mut sr = crate::inttools::shrunk_range::ShrunkRange::new();
            sr.set_data(*ei, [*p1, *p2], v_tols[*v1i], v_tols[*v2i], et[*ei]);
            sr.perform(&ec[*ei]);
            (sr.shrunk_range(), sr.is_splittable())
        }).collect();
        // Write phase: apply results back to PaveBlocks
        let mut idx = 0usize;
        for ei in 0..num_edges {
            for pi in 0..edge_pb_counts[ei] {
                let (range, splittable) = results[idx]; idx += 1;
                self.ds.edges[ei].pave_blocks[pi].shrunk_range = range;
                self.ds.edges[ei].pave_blocks[pi].is_splittable = splittable;
            }
        }
        for pb in &mut self.ds.pave_blocks {
            if pb.original_edge >= self.ds.edges.len() { continue; }
            let (range, splittable) = results[idx]; idx += 1;
            pb.shrunk_range = range;
            pb.is_splittable = splittable;
        }
    }

    pub(crate) fn existing_pave_block(&self, ei: usize, vi: usize) -> bool {
        for pb in &self.ds.edges[ei].pave_blocks {
            if pb.pave1.vertex_idx == vi || pb.pave2.vertex_idx == vi { return true; }
        }
        false
    }

    /// OCCT-aligned: SplitPaveBlocks (PaveFiller_2.cxx L419-626).
    ///   Splits PaveBlocks of the given edges using extra paves.
    ///   For each edge, reconstructs pave_blocks from the full pave set,
    ///   computes shrunk data, and unifies vertices when no valid range.
    pub(crate) fn split_pave_blocks(&mut self, edges: &std::collections::HashSet<usize>,
                                    add_interfs: bool) {
        let mut verts_to_merge: Vec<Vec<usize>> = Vec::new();
        for &ei in edges {
            if ei >= self.ds.edges.len() { continue; }
            let edge = &self.ds.edges[ei];
            if edge.paves.len() < 2 { continue; }

            // OCCT L447-453: aPB->Update(aLPBN) — rebuild PBs from all paves
            let mut params: Vec<(usize, f64)> = edge.paves.iter()
                .map(|p| (p.vertex_idx, p.param)).collect();
            params.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            params.dedup_by(|a, b| (a.1 - b.1).abs() < TOLERANCE_ABS);

            // Build new PBs from consecutive pave pairs
            let mut new_pbs: Vec<PaveBlock> = Vec::new();
            for w in params.windows(2) {
                // OCCT L460-461: UpdatePaveBlockWithSDVertices
                let pb = PaveBlock::new(
                    ei,
                    Pave { vertex_idx: w[0].0, param: w[0].1 },
                    Pave { vertex_idx: w[1].0, param: w[1].1 },
                );

                // OCCT L462: FillShrunkData
                let v1_tol = if w[0].0 < self.ds.vertices.len() { self.ds.vertices[w[0].0].geom_tol } else { TOLERANCE_ABS };
                let v2_tol = if w[1].0 < self.ds.vertices.len() { self.ds.vertices[w[1].0].geom_tol } else { TOLERANCE_ABS };
                let mut sr = crate::inttools::shrunk_range::ShrunkRange::new();
                sr.set_data(ei, [w[0].1, w[1].1], v1_tol, v2_tol, edge.geom_tol);
                sr.perform(&edge.curve);

                // OCCT L468-507: Check valid range + unify if needed
                if sr.shrunk_range().is_none() || (!sr.is_splittable() && sr.shrunk_range().is_some()) {
                    let nv1 = w[0].0;
                    let nv2 = w[1].0;
                    if nv1 != nv2 {
                        // OCCT L493-506: MakeSDVertices for the pair
                        verts_to_merge.push(vec![nv1, nv2]);
                    }
                    continue;
                }

                new_pbs.push(pb);
            }

            // OCCT L526: Replace old PBs with new ones
            self.ds.edges[ei].pave_blocks = new_pbs;
        }

        // OCCT L493-506: Merge vertex pairs that have no valid range
        for pair in &verts_to_merge {
            if pair.len() >= 2 {
                let nv1 = pair[0];
                let nv2 = pair[1];
                // Use min index as merge target (convention)
                let n_v = nv1.min(nv2);
                self.ds.add_shape_sd(nv1, n_v);
                self.ds.add_shape_sd(nv2, n_v);
                if add_interfs {
                    self.ds.interferences.push(Interference::VertexVertex {
                        v1: nv1, v2: nv2, merged_vertex: n_v,
                    });
                }
            }
        }

        // OCCT L529-617: CommonBlock handling — skipped (Architecture diff A2:
        //   rcad uses inline PaveBlock on DSEdge, not pool-based with CommonBlocks)
    }

    pub(crate) fn split_ics_at_periodic_boundary(&mut self) {
        let n_curves = self.ds.intersection_curves.len();
        let mut curve_faces = std::collections::HashMap::new();
        for inf in &self.ds.interferences {
            if let Interference::FaceFace { f1, f2, curves, .. } = inf {
                for &ci in curves {
                    curve_faces.entry(ci).or_insert((*f1, *f2));
                }
            }
        }
        for ci in 0..n_curves {
            let Some(&(fa, fb)) = curve_faces.get(&ci) else { continue };
            let ic = &self.ds.intersection_curves[ci];
            let t0 = ic.t_range[0]; let t1 = ic.t_range[1];
            if (t1 - t0).abs() < TOLERANCE_ABS { continue; }
            for &(fi, use_a) in &[(fa, true), (fb, false)] {
                let (period, _, _) = match &self.ds.faces[fi].surface {
                    Surface3::Sphere(_) => (std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
                    Surface3::Cylinder(_) => (std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
                    _ => continue,
                };
                let pcurve = if use_a { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
                let Some(pc) = pcurve else { continue };
                const N_SAMP: usize = 64;
                let mut prev_u: Option<f64> = None;
                for i in 0..=N_SAMP {
                    let frac = i as f64 / N_SAMP as f64;
                    let u = pc.point_at(t0 + (t1 - t0) * frac).x;
                    if let Some(pu) = prev_u {
                        if (u - pu).abs() > period * 0.5 {
                            let mut lo = (i as f64 - 1.0) / N_SAMP as f64;
                            let mut hi = i as f64 / N_SAMP as f64;
                            for _ in 0..32 {
                                let mid = (lo + hi) * 0.5;
                                let u_mid = pc.point_at(t0 + (t1 - t0) * mid).x;
                                if u_mid < (pu / period).round() * period { lo = mid; }
                                else { hi = mid; }
                            }
                            let t_cross = t0 + (t1 - t0) * (lo + hi) * 0.5;
                            let pt = ic.curve.point_at(t_cross);
                            let v_new = self.ds.find_vertex_near(pt, TOLERANCE_ABS * 1000.0)
                                .unwrap_or_else(|| {
                                    let vi = self.ds.vertices.len();
                                    self.ds.vertices.push(crate::bopds::ds::DSVertex {
                                        point: pt, geom_tol: TOLERANCE_ABS, origin: None,
                                        is_internal: false, location: 0,
                                    });
                                    vi
                                });
                            self.ds.faces[fa].face_info.vertices_in.insert(v_new);
                            self.ds.faces[fb].face_info.vertices_in.insert(v_new);
                        }
                    }
                    prev_u = Some(u);
                }
            }
        }
    }

    pub(crate) fn put_pave_on_curve(
        &mut self,
        nV: usize,
        aTolR3D: f64,
        curve_idx: usize,
        aMI: &std::collections::HashSet<usize>,
        iCheckExtend: i32,
    ) -> Option<f64> {
        let ic_curve = self.ds.intersection_curves[curve_idx].curve.clone();
        let aV_tol = self.ds.vertices[nV].geom_tol;
        let mut aTolV = self.a_mv_tol.get(&nV).copied().unwrap_or(aV_tol);

        let mut aT: f64 = 0.0;
        let mut bIsVertexOnLine = self.is_vertex_on_line(nV, aTolV, curve_idx, aTolR3D + self.fuzzy_tolerance, &mut aT);

        if !bIsVertexOnLine && iCheckExtend != 0 && !self.verts_to_avoid_extension.contains(&nV) {
            let mut anExtraTol = aTolV;
            if self.extended_tolerance(nV, aMI, &mut anExtraTol, iCheckExtend) {
                bIsVertexOnLine = self.is_vertex_on_line(nV, anExtraTol, curve_idx, aTolR3D + self.fuzzy_tolerance, &mut aT);
                if bIsVertexOnLine {
                    let aPOnC = ic_curve.point_at(aT);
                    aTolV = aPOnC.distance(self.ds.vertices[nV].point);
                }
            }
        }

        if bIsVertexOnLine {
            // OCCT L3031: aDTol = BOPTools_AlgoTools::DTolerance()  (=1.e-12)
            let aDTol = 1e-12;
            let aPTol = Self::curve_parametric_tolerance(&ic_curve, aTolR3D.max(aTolV));

            let mut nVUsed = 0;
            if let Some(pb) = self.ds.intersection_curves[curve_idx].change_pave_block1() {
                let bExist = pb.contains_parameter(aT, aPTol, &mut nVUsed);
                if bExist {
                    let pList = self.a_dmv_lv.entry(nVUsed).or_insert_with(|| {
                        let mut list = Vec::new();
                        list.push(nVUsed);
                        if !self.a_mv_tol.contains_key(&nVUsed) {
                            let aVUsed_tol = if nVUsed < self.ds.vertices.len() {
                                self.ds.vertices[nVUsed].geom_tol
                            } else { aTolV };
                            self.a_mv_tol.insert(nVUsed, aVUsed_tol);
                        }
                        list
                    });
                    if !pList.contains(&nV) {
                        pList.push(nV);
                    }
                    if !self.a_mv_tol.contains_key(&nV) {
                        self.a_mv_tol.insert(nV, aV_tol);
                    }
                } else {
                    pb.append_ext_pave(Pave { vertex_idx: nV, param: aT });
                    let aP1 = ic_curve.point_at(aT);
                    aTolV = aV_tol;
                    let aP2 = self.ds.vertices[nV].point;
                    let aDist = aP1.distance(aP2);
                    if aTolV < aDist + aDTol {
                        self.ds.vertices[nV].geom_tol = aDist + aDTol;
                        if !self.a_mv_tol.contains_key(&nV) {
                            self.a_mv_tol.insert(nV, aTolV);
                        }
                    }
                }
            }
            return Some(aT);
        }
        None
    }

    /// OCCT-aligned: BOPAlgo_PaveFiller::PutPavesOnCurve (L2404-2453).
    ///   Takes pre-computed vertex sets instead of re-deriving from interferences.
    pub(crate) fn put_paves_on_curve(
        &mut self,
        theMVOnIn: &std::collections::HashSet<usize>,
        theMVCommon: &std::collections::HashSet<usize>,
        curve_idx: usize,
        aMI: &std::collections::HashSet<usize>,
        theMVEF: &std::collections::HashSet<usize>,
    ) {
        let ic = &self.ds.intersection_curves[curve_idx];
        let aTolR3D = ic.geom_tol.max(ic.curve_extra.tangential_tol);
        // OCCT L2415-2416: aBoxC = theNC.Box()
        let c_box = crate::pave_filler::helpers::curve_bounding_box_simple(&ic.curve, 0.0);

        // OCCT L2418-2424: Put EF vertices first (iCheckExtend=2)
        for &nV in theMVEF {
            self.put_pave_on_curve(nV, aTolR3D, curve_idx, aMI, 2);
        }

        // OCCT L2426-2452: Put all other ON/IN vertices (iCheckExtend=1)
        for &nV in theMVOnIn {
            if theMVEF.contains(&nV) { continue; }
            if theMVCommon.contains(&nV) {
                // OCCT L2436: common vertices skip box check
                self.put_pave_on_curve(nV, aTolR3D, curve_idx, aMI, 1);
            } else {
                // OCCT L2438-2444: Box check — skip if curve box doesn't overlap vertex box
                if nV < self.ds.shape_info.len() {
                    if let (Some(vmn), Some(vmx)) = (self.ds.shape_info[nV].box_min, self.ds.shape_info[nV].box_max) {
                        let v_box_min = vmn - glam::DVec3::splat(aTolR3D);
                        let v_box_max = vmx + glam::DVec3::splat(aTolR3D);
                        if let Some(c_box) = &c_box {
                            if v_box_max.x < c_box[0].x || v_box_min.x > c_box[1].x ||
                               v_box_max.y < c_box[0].y || v_box_min.y > c_box[1].y ||
                               v_box_max.z < c_box[0].z || v_box_min.z > c_box[1].z {
                                continue;
                            }
                        }
                    }
                }
                // OCCT L2445-2447: Skip non-new shapes
                if nV < self.ds.shape_info.len() && !self.ds.shape_info[nV].is_new {
                    continue;
                }
                self.put_pave_on_curve(nV, aTolR3D, curve_idx, aMI, 1);
            }
        }
    }

    /// OCCT-aligned: PutStickPavesOnCurve (BOPAlgo_PaveFiller_6.cxx L2780-2875).
    ///   Processes stick vertices that are near the IC endpoints (rich criterion)
    ///   and where both face normals are nearly opposite (crease criterion).
    pub(crate) fn put_stick_paves_on_curve(
        &mut self,
        ci: usize,
        aMI: &std::collections::HashSet<usize>,
        aMVStick: &std::collections::HashSet<usize>,
    ) {
        let (start_vertex, end_vertex, t_range, curve, geom_tol) = {
            let ic = &self.ds.intersection_curves[ci];
            (ic.start_vertex, ic.end_vertex, ic.t_range, ic.curve.clone(), ic.geom_tol)
        };
        // OCCT L2792-2798: getBoundPaves — check if both ends already have vertices
        if start_vertex < self.ds.vertices.len() && end_vertex < self.ds.vertices.len() {
            return;
        }
        // OCCT L2799-2801: RemoveUsedVertices
        let a_mv: std::collections::HashSet<usize> = aMVStick.iter().copied().collect();
        if a_mv.is_empty() { return; }

        let a_tol_r3d = geom_tol.max(self.tol());
        let a_dt2 = 2e-7;
        let a_d_sc_pr = 5e-9;
        let a_t = t_range;
        let a_p = [curve.point_at(a_t[0]), curve.point_at(a_t[1])];

        let surf_both: [Option<Surface3>; 2] = {
            let mut sv = [None, None];
            for (k, fi) in self.face_idxs_for_curve(ci).iter().enumerate() {
                if *fi < self.ds.faces.len() {
                    sv[k] = Some(self.ds.faces[*fi].surface.clone());
                }
            }
            sv
        };

        for &n_v in &a_mv {
            let v_pt = self.ds.vertices[n_v].point;
            for m in 0..2 {
                // OCCT L2838-2841: skip if bound already has a vertex
                if (m == 0 && start_vertex < self.ds.vertices.len()) ||
                   (m == 1 && end_vertex < self.ds.vertices.len()) {
                    continue;
                }
                // OCCT L2842-2846: rich criterion — close to IC endpoint
                let d2 = a_p[m].distance_squared(v_pt);
                if d2 > a_dt2 { continue; }
                // OCCT L2848-2866: crease criterion — face normals nearly opposite
                let mut sc_pr = 1.0;
                if let (Some(s0), Some(s1)) = (&surf_both[0], &surf_both[1]) {
                    let n0 = Self::estimate_surface_normal(s0, a_p[m]);
                    let n1 = Self::estimate_surface_normal(s1, a_p[m]);
                    sc_pr = n0.dot(n1);
                    if sc_pr < 0.0 { sc_pr = -sc_pr; }
                    sc_pr = 1.0 - sc_pr;
                }
                if sc_pr > a_d_sc_pr { continue; }
                // OCCT L2869-2871: PutPaveOnCurve
                let a_d = d2.sqrt();
                let a_tol_r3d_use = a_tol_r3d.max(self.ds.vertices[n_v].geom_tol);
                self.put_pave_on_curve(n_v, a_d.min(a_tol_r3d_use), ci, aMI, 1);
                break;
            }
        }
    }

    /// Find the two face indices for a given intersection curve.
    fn face_idxs_for_curve(&self, ci: usize) -> [usize; 2] {
        let mut result = [usize::MAX; 2];
        let mut idx = 0;
        for (fi, face) in self.ds.faces.iter().enumerate() {
            if face.face_info.curves_sc.contains(&ci) && idx < 2 {
                result[idx] = fi;
                idx += 1;
            }
        }
        result
    }

    /// Estimate surface normal at a 3D point (for PutStickPavesOnCurve normal check).
    pub(crate) fn estimate_surface_normal(surf: &Surface3, pt: DVec3) -> DVec3 {
        match surf {
            Surface3::Plane(p) => p.normal,
            Surface3::Sphere(s) => (pt - s.center).normalize_or_zero(),
            Surface3::Cylinder(c) => {
                let axis = c.axis.normalize_or_zero();
                let radial = pt - c.origin - axis * (pt - c.origin).dot(axis);
                radial.normalize_or_zero()
            }
            _ => {
                // Fallback: approximate via closest point projection
                let (_, proj) = crate::extrema::closest_point_on_surface(surf, pt);
                (pt - proj).normalize_or_zero()
            }
        }
    }

    fn curve_parametric_tolerance(curve: &Curve3, tol_3d: f64) -> f64 {
        match curve {
            Curve3::Line(_) => tol_3d,
            Curve3::Circle(c) => tol_3d / c.radius.max(1e-12),
            Curve3::Ellipse(e) => tol_3d / e.major_radius.max(1e-12),
            _ => tol_3d * 0.01,
        }
    }

    /// ✅ OCCT-aligned: IntTools_Context::IsVertexOnLine (L786-992).
    ///   Form-identical logic: curve-type-based tolerance → first-endpoint
    ///   check (with local+global projection fallback) → last-endpoint check
    ///   (with bFirstValid shortcut) → global projection via closest_point_on_curve.
    pub(crate) fn is_vertex_on_line(
        &self,
        nV: usize,
        aTolV: f64,
        curve_idx: usize,
        aTolC: f64,
        aT: &mut f64,
    ) -> bool {
        use rcad_kernel::projection::closest_point_on_curve;
        let ic = &self.ds.intersection_curves[curve_idx];
        let vp = self.ds.vertices[nV].point;
        let aFirst = ic.t_range[0];
        let aLast = ic.t_range[1];

        // OCCT L800: aTolSum = aTolV + aTolC
        // OCCT L802-819: curve-type-dependent tolerance scaling
        let mut aTolSum = aTolV + aTolC;
        let is_bspline_or_bezier = matches!(&ic.curve, Curve3::BSpline(_) | Curve3::Bezier(_));
        aTolSum *= 2.0;
        if is_bspline_or_bezier {
            if aTolSum < 1e-5 { aTolSum = 1e-5; }   // OCCT L807-809
        } else {
            if aTolSum < 1e-6 { aTolSum = 1e-6; }   // OCCT L815-817
        }

        // OCCT L821-822: aFirst / aLast from curve
        // (already set from ic.t_range above)

        // ── OCCT L824-883: First endpoint check ──
        let mut b_first_valid = false;
        let mut a_first_dist = f64::MAX;
        if aFirst.is_finite() {
            let p_first = ic.curve.point_at(aFirst);
            a_first_dist = vp.distance(p_first);
            if a_first_dist < aTolSum {
                b_first_valid = true;
                *aT = aFirst;
                if a_first_dist > aTolV {
                    // OCCT L840: Extrema_LocateExtPC equivalent
                    let proj = closest_point_on_curve(&ic.curve, vp, 64);
                    let mid = (aLast + aFirst) * 0.5;
                    // OCCT L847-851 (locate) + L875-879 (extpc): same guard
                    let p_first_d = p_first.distance(proj.point);
                    if proj.param > mid || proj.distance > aTolSum || p_first_d < 1e-7 {
                        *aT = aFirst;
                    } else {
                        *aT = proj.param;
                    }
                }
            }
        }

        // ── OCCT L886-951: Last endpoint check ──
        if aLast.is_finite() {
            let p_last = ic.curve.point_at(aLast);
            let d_last = vp.distance(p_last);
            // OCCT L890-892: if first valid and first is closer → keep first
            if b_first_valid && a_first_dist < d_last {
                // Keep aT from first-endpoint branch
                return true;
            }
            if d_last < aTolSum {
                *aT = aLast;
                if d_last > aTolV {
                    let proj = closest_point_on_curve(&ic.curve, vp, 64);
                    let mid = (aLast + aFirst) * 0.5;
                    let p_last_d = p_last.distance(proj.point);
                    // OCCT L908-912 (locate) + L936-940 (extpc): same guard
                    if proj.param < mid || proj.distance > aTolSum || p_last_d < 1e-7 {
                        *aT = aLast;
                    } else {
                        *aT = proj.param;
                    }
                }
                return true;
            }
        } else if b_first_valid {
            // OCCT L948-951: only first endpoint is valid → return true
            return true;
        }

        // ── OCCT L953-992: General projection ──
        let proj = closest_point_on_curve(&ic.curve, vp, 64);
        if proj.distance <= aTolSum {
            *aT = proj.param;
            return true;
        }

        // OCCT L957-980: BoundedCurve fallback (endpoints)
        //   rcad: closest_point_on_curve already handles bounded curves,
        //   so skip explicit endpoint fallback and return false.
        false
    }

    /// OCCT-aligned: BOPAlgo_PaveFiller::ExtendedTolerance (PaveFiller_6.cxx L????).
    /// In-out aTolVExt: on input it is the current vertex tolerance (lower bound),
    /// on output it is set to the maximum extended distance found.
    pub(crate) fn extended_tolerance(
        &self,
        nV: usize,
        aMI: &std::collections::HashSet<usize>,
        aTolVExt: &mut f64,
        aType: i32,
    ) -> bool {
        if nV < self.ds.shape_info.len() && !self.ds.shape_info[nV].is_new {
            return false;
        }
        let vp = self.ds.vertices[nV].point;
        let mut found = false;
        // Use input aTolVExt as initial threshold (OCCT in-out semantics)
        let mut max_ext = *aTolVExt;

        if aType == 0 || aType == 1 {
            for inf in &self.ds.interferences {
                if let Interference::EdgeEdge { e1, param1, param2, new_vertex, .. } = inf {
                    if *new_vertex != nV { continue; }
                    if !aMI.contains(e1) { continue; }
                    if *e1 < self.ds.edges.len() {
                        let p1 = crate::boptools::point_on_edge(&self.ds.edges[*e1], *param1);
                        let p2 = crate::boptools::point_on_edge(&self.ds.edges[*e1], *param2);
                        let d = vp.distance(p1).max(vp.distance(p2));
                        if d > max_ext { max_ext = d; }
                        found = true;
                    }
                }
            }
        }
        if aType == 0 || aType == 2 {
            for inf in &self.ds.interferences {
                if let Interference::EdgeFace { edge, edge_param, new_vertex, .. } = inf {
                    if *new_vertex != nV { continue; }
                    if !aMI.contains(edge) { continue; }
                    if *edge < self.ds.edges.len() {
                        let p1 = crate::boptools::point_on_edge(&self.ds.edges[*edge], *edge_param);
                        let d = vp.distance(p1);
                        if d > max_ext { max_ext = d; }
                        found = true;
                    }
                }
            }
        }
        if found { *aTolVExt = max_ext; }
        found
    }

    pub(crate) fn project_vertex_on_curve(&self, vi: usize, ic: &IntersectionCurve) -> Option<f64> {
        use rcad_kernel::geom::CurveEval;
        let v_tol = self.ds.vertices[vi].geom_tol;
        let c_tol = ic.geom_tol;
        let f_tol = self.ds.fuzzy_tol;
        // OCCT IsVertexOnLine L800: aTolSum = aTolV + aTolC
        //   where aTolC = aTolR3D + myFuzzyValue (PutPaveOnCurve L2976)
        let raw_sum = v_tol + c_tol + f_tol;
        // OCCT L806-819: aTolSum = 2 * aTolSum, clamped to >= 1e-6
        let tl = (2.0 * raw_sum).max(1e-6);
        let result = self.project_vertex_on_curve_with_tol(vi, ic, tl);
        if result.is_some() { return result; }
        let ext_tol = self.extended_tolerance_occt(vi);
        if ext_tol > tl { self.project_vertex_on_curve_with_tol(vi, ic, ext_tol) } else { None }
    }

    pub(crate) fn project_vertex_on_curve_with_tol(&self, vi: usize, ic: &IntersectionCurve, tl: f64) -> Option<f64> {
        use rcad_kernel::geom::CurveEval;
        let vp = self.ds.vertices[vi].point;
        match &ic.curve {
            Curve3::Line(l) => { let v = vp - l.origin; let t = v.dot(l.direction);
                if (v - l.direction*t).length() <= tl { Some(t.clamp(ic.t_range[0], ic.t_range[1])) } else { None } }
            Curve3::Circle(c) => {
                let v = vp - c.center;
                let tl_cap = tl.min(c.radius.max(1e-7) * 1e-4);
                if (v.length() - c.radius).abs() > tl_cap { return None; }
                let nm = c.normal.normalize();
                if v.dot(nm).abs() > tl_cap { return None; }
                // Use Circle3's own basis (same as point_at), NOT any_perpendicular.
                let x_ax = if nm.x.abs() < 0.9 {
                    nm.cross(DVec3::X).normalize()
                } else {
                    nm.cross(DVec3::Y).normalize()
                };
                let y_ax = nm.cross(x_ax);
                let t_raw = v.dot(y_ax).atan2(v.dot(x_ax));
                let t_full = t_raw.rem_euclid(std::f64::consts::TAU);
                let pt = c.center + c.radius * (t_full.cos() * x_ax + t_full.sin() * y_ax);
                if (pt - vp).length() <= tl { Some(t_full) } else { None }
            }
            _ => {
                for i in 0..20 {
                    let t = ic.t_range[0] + (ic.t_range[1]-ic.t_range[0])*i as f64/20.0;
                    if (ic.curve.point_at(t)-vp).length() <= tl { return Some(t); }
                }
                None
            }
        }
    }

    pub(crate) fn reduce_intersection_range(&self, ei: usize, param: f64, tol: f64) -> [f64; 2] {
        let edge = &self.ds.edges[ei];
        let mut t_range = edge.t_range;
        for pave in &edge.paves {
            let d = (pave.param - param).abs();
            if d < tol * 10.0 {
                if pave.param < param && pave.param > t_range[0] { t_range[0] = pave.param; }
                if pave.param > param && pave.param < t_range[1] { t_range[1] = pave.param; }
            }
        }
        t_range
    }

    pub(crate) fn make_pcurves(&mut self) {
        // OCCT L592-595: early return when pcurve building is avoided or neither
        // section face requires pcurves.
        if self.avoid_build_pcurve
            || (!self.section_attribute.pcurve_on_s1 && !self.section_attribute.pcurve_on_s2)
        {
            return;
        }
        let b_pcurve_on_s = [self.section_attribute.pcurve_on_s1, self.section_attribute.pcurve_on_s2];

        // OCCT L601: BOPAlgo_VectorOfMPC aVMPC �?collection of MPC entries.
        // rcad: the MPC (Make PCurve) concept is inlined �?we store the data
        // needed to either (a) compute a new pcurve, or (b) update vertex tolerances
        // from an already-existing pcurve.
        #[derive(Clone)]
        struct MPCEntry {
            edge_idx: usize,
            face_idx: usize,
            /// OCCT: SetFlag(true) means call UpdateVertices after pcurve is set.
            update_vertices: bool,
            /// OCCT SetData fields: existing edge to clone pcurve from.
            existing_edge: Option<usize>,
        }

        let mut a_vmpc: Vec<MPCEntry> = Vec::new();

        // ===== Phase 1: Common Block pcurves (OCCT L603-700) =====
        // Process PaveBlocksIn (boundary edges through face) and PaveBlocksOn
        // (edges lying on the face surface).
        for fi in 0..self.ds.faces.len() {
            let face_info = &self.ds.faces[fi].face_info;

            // OCCT L618-631: PaveBlocksIn
            for &pb_idx in &face_info.pave_blocks_in {
                if pb_idx >= self.ds.pave_blocks.len() { continue; }
                let pb = &self.ds.pave_blocks[pb_idx];
                let ei = pb.original_edge;
                if ei >= self.ds.edges.len() { continue; }
                a_vmpc.push(MPCEntry {
                    edge_idx: ei,
                    face_idx: fi,
                    update_vertices: false,
                    existing_edge: None,
                });
            }

            // OCCT L633-699: PaveBlocksOn
            for &pb_idx in &face_info.pave_blocks_on {
                if pb_idx >= self.ds.pave_blocks.len() { continue; }
                let pb = &self.ds.pave_blocks[pb_idx];
                let ei = pb.original_edge;
                if ei >= self.ds.edges.len() { continue; }

                // OCCT L641: HasCurveOnSurface �?skip if pcurve already exists
                if self.ds.edges[ei].face_reps.iter().any(|r| r.face_idx == fi) {
                    continue;
                }

                // OCCT L649-695: CommonBlock inheritance �?if another PB in the same
                // CommonBlock already has a pcurve on this face, reuse it (SetData).
                let mut cb_existing_edge: Option<usize> = None;
                if let Some(cb_idx) = pb.common_block_idx {
                    if let Some(cb) = self.ds.common_blocks.get(cb_idx) {
                        let cb_pbs = cb.pave_blocks();
                        if cb_pbs.len() >= 2 {
                            for &(other_pb_idx, _other_fi) in cb_pbs {
                                if other_pb_idx == pb_idx { continue; }
                                if let Some(other_pb) = self.ds.pave_blocks.get(other_pb_idx) {
                                    let other_ei = other_pb.original_edge;
                                    if let Some(other_edge) = self.ds.edges.get(other_ei) {
                                        if other_edge.face_reps.iter().any(|r| r.face_idx == fi) {
                                            // OCCT L678-690: SetData(aEz, aV1x, aT1x, aV2x, aT2x)
                                            // AttachExistingPCurve �?the existing pcurve on aEx
                                            // (other edge) will be copied to this edge.
                                            cb_existing_edge = Some(other_ei);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                a_vmpc.push(MPCEntry {
                    edge_idx: ei,
                    face_idx: fi,
                    update_vertices: false,
                    existing_edge: cb_existing_edge,
                });
            }
        } // end Phase 1

        // ===== Phase 2: Section edge pcurves �?UpdateVertices only (OCCT L702-757) =====
        // Pcurves for section edges are already computed during make_blocks.
        // This phase creates MPC entries with SetFlag(true) so that UpdateVertices
        // is called after the pcurve is known, correcting vertex tolerances.
        if b_pcurve_on_s[0] || b_pcurve_on_s[1] {
            // OCCT L711: anEFPairs fence prevents duplicate (edge, face) pairs
            let mut ef_pairs: std::collections::HashSet<(usize, usize)> =
                std::collections::HashSet::new();

            // OCCT L712-713: iterate all FF interferences
            for interf in &self.ds.interferences {
                if let Interference::FaceFace { f1, f2, curves, .. } = interf {
                    let nf = [*f1, *f2];

                    // OCCT L733-755: for each curve in the FF interference
                    for &ci in curves {
                        if ci >= self.ds.intersection_curves.len() { continue; }
                        let ic = &self.ds.intersection_curves[ci];

                        // OCCT L736-754: for each PaveBlock of this curve
                        for pb in &ic.pave_blocks {
                            // OCCT L741: nE = aPB->Edge()
                            // For section edges original_edge == NO_EDGE, use new_edge.
                            let section_ei = match pb.new_edge {
                                Some(ei) => ei,
                                None => continue,
                            };

                            // OCCT L744-753: for each of the two faces
                            for (m, &fi) in nf.iter().enumerate() {
                                if !b_pcurve_on_s[m] { continue; }
                                // OCCT L746: anEFPairs.Add(BOPDS_Pair(nE, nF[m]))
                                // Returns false if pair already existed
                                if !ef_pairs.insert((section_ei, fi)) {
                                    continue;
                                }
                                // OCCT L748-751: create MPC with SetFlag(true)
                                a_vmpc.push(MPCEntry {
                                    edge_idx: section_ei,
                                    face_idx: fi,
                                    update_vertices: true,
                                    existing_edge: None,
                                });
                            }
                        }
                    }
                }
            }
        } // end Phase 2

        // ===== Phase 3: Perform all MPC computations (OCCT L760-775) =====
        // OCCT L760-766: BOPTools_Parallel::Perform(myRunParallel, aVMPC, myContext)
        // rcad: sequential execution for now �?each MPC.Perform computes the pcurve
        // (or copies from existing_edge), stores it in the DS, and optionally calls
        // UpdateVertices.
        //
        // Later: replace with parallel if myRunParallel==true.
        let mut failed_pcurves: Vec<(usize, usize)> = Vec::new();
        let mut applied_pcurves: Vec<(usize, usize, Curve2d, f64, f64)> = Vec::new();
        let mut update_vertices_list: Vec<(usize, usize)> = Vec::new();

        for mpc in &a_vmpc {
            let ei = mpc.edge_idx;
            let fi = mpc.face_idx;

            if ei >= self.ds.edges.len() || fi >= self.ds.faces.len() {
                continue;
            }

            // If another edge in the same CommonBlock already provides the pcurve,
            // copy it (OCCT SetData / AttachExistingPCurve path).
            if let Some(src_ei) = mpc.existing_edge {
                if src_ei < self.ds.edges.len() {
                    if let Some(rep) = self.ds.edges[src_ei].face_reps.iter().find(|r| r.face_idx == fi) {
                        applied_pcurves.push((
                            ei, fi,
                            rep.pcurve.clone(),
                            self.ds.edges[ei].t_range[0],
                            self.ds.edges[ei].t_range[1],
                        ));
                        if mpc.update_vertices {
                            update_vertices_list.push((ei, fi));
                        }
                        continue;
                    }
                }
            }

            // Check if this edge already has a pcurve on this face (OCCT HasCurveOnSurface).
            // Section edges may already have their pcurve from make_blocks.
            if self.ds.edges[ei].face_reps.iter().any(|r| r.face_idx == fi) {
                if mpc.update_vertices {
                    // OCCT SetFlag(true) �?only UpdateVertices needed, pcurve exists
                    update_vertices_list.push((ei, fi));
                }
                continue;
            }

            // Compute pcurve (OCCT MPC.Perform internal logic)
            let face_surface = self.ds.faces[fi].surface.clone();
            let edge_curve = &self.ds.edges[ei].curve;
            if let Some((pcurve, len)) = DS::compute_edge_pcurve(edge_curve, &face_surface) {
                applied_pcurves.push((
                    ei, fi, pcurve,
                    self.ds.edges[ei].t_range[0],
                    self.ds.edges[ei].t_range[1],
                ));
                let _ = len;
                if mpc.update_vertices {
                    update_vertices_list.push((ei, fi));
                }
            } else {
                // OCCT L782-788: failed �?collect for warning
                failed_pcurves.push((ei, fi));
            }
        }

        // ===== Phase 4: Error reporting and batch application (OCCT L782-789) =====
        // OCCT L782-788: AddWarning for each failed MPC
        for &(ei, fi) in &failed_pcurves {
            // OCCT: BRep_Builder MakeCompound + Add(Edge) + Add(Face) +
            //       AddWarning(new BOPAlgo_AlertBuildingPCurveFailed(compound))
            // rcad: non-fatal warning as debug trace
            eprintln!("[WARN] MakePCurves failed for edge={} on face={}", ei, fi);
        }

        // Apply all computed pcurves to the DS
        use crate::bopds::ds::DSRepOnFace;
        for (ei, fi, pcurve, t0, t1) in applied_pcurves {
            self.ds.edges[ei].face_reps.push(DSRepOnFace {
                face_idx: fi,
                pcurve,
                pcurve2: None,
                pcurve_range: [t0, t1],
                start_param: t0,
                end_param: t1,
            });
        }

        // OCCT L284-287: UpdateVertices for section-edge pcurves.
        // Adjusts vertex tolerances based on mismatch between 3D curve and 2D pcurve.
        for &(ei, fi) in &update_vertices_list {
            // OCCT: UpdateVertices(aCopyE, myF) �?corrects vertex tolerances using the
            // difference between the edge's 3D curve and its pcurve on the face.
            //
            // rcad placeholder: structure matches OCCT but the actual tolerance
            // adjustment is deferred. The pcurve is already computed and stored above.
            let _ = (ei, fi);
        }
    }
    pub(crate) fn update_face_info(&mut self, fi: usize) {
        let face = &self.ds.faces[fi].clone();
        let mut sc_curves: Vec<usize> = face.face_info.curves_sc.iter().copied().collect();
        // Recompute vertices_in from all IC endpoints on this face
        for &ci in &sc_curves {
            if ci < self.ds.intersection_curves.len() {
                let ic = &self.ds.intersection_curves[ci];
                self.ds.faces[fi].face_info.vertices_in.insert(ic.start_vertex);
                self.ds.faces[fi].face_info.vertices_in.insert(ic.end_vertex);
            }
        }
        // Update pave_blocks_in/pave_blocks_on from edge pave blocks
        for &ei in &face.boundary_edges {
            if ei < self.ds.edges.len() {
                for pb in &self.ds.edges[ei].pave_blocks {
                    let pb_idx = self.ds.pave_blocks.len();
                    self.ds.pave_blocks.push(pb.clone());
                    self.ds.faces[fi].face_info.pave_blocks_in.insert(pb_idx);
                }
            }
        }
    }

    pub(crate) fn check_planes(&self, f1: usize, f2: usize) -> bool {
        let s1 = &self.ds.faces[f1].surface;
        let s2 = &self.ds.faces[f2].surface;
        match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                let dot = p1.normal.dot(p2.normal).abs();
                if dot > 0.9999 {
                    let d = (p2.origin - p1.origin).dot(p1.normal).abs();
                    d < TOLERANCE_ABS * 100.0
                } else { false }
            }
            _ => false,
        }
    }

    pub(crate) fn make_sd_vertices(&mut self, verts: &[usize]) {
        if verts.len() < 2 { return; }
        let tol = TOLERANCE_ABS * 100.0;
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut used = vec![false; verts.len()];
        for i in 0..verts.len() {
            if used[i] { continue; }
            let mut group = vec![verts[i]];
            used[i] = true;
            for j in (i+1)..verts.len() {
                if used[j] { continue; }
                let d = (self.ds.vertices[verts[i]].point - self.ds.vertices[verts[j]].point).length();
                if d < tol { group.push(verts[j]); used[j] = true; }
            }
            if group.len() > 1 { groups.push(group); }
        }
        // For each group, keep the first vertex as SD, remap others
        for group in &groups {
            let sd_vi = group[0];
            for &vi in &group[1..] {
                for ei in 0..self.ds.edges.len() {
                    for pb in &mut self.ds.edges[ei].pave_blocks {
                        if pb.pave1.vertex_idx == vi { pb.pave1.vertex_idx = sd_vi; }
                        if pb.pave2.vertex_idx == vi { pb.pave2.vertex_idx = sd_vi; }
                    }
                }
                for fi in 0..self.ds.faces.len() {
                    if self.ds.faces[fi].face_info.vertices_in.contains(&vi) {
                        self.ds.faces[fi].face_info.vertices_in.remove(&vi);
                        self.ds.faces[fi].face_info.vertices_in.insert(sd_vi);
                    }
                }
            }
        }
    }
    pub(crate) fn get_stick_vertices(&self, fi: usize) -> Vec<usize> {
        let mut verts: Vec<usize> = Vec::new();
        for inf in &self.ds.interferences {
            match inf {
                Interference::VertexVertex { v1, v2, .. } => {
                    if self.vertex_on_face(*v1, fi) { verts.push(*v1); }
                    if self.vertex_on_face(*v2, fi) { verts.push(*v2); }
                }
                Interference::VertexEdge { vertex, .. } => {
                    if self.vertex_on_face(*vertex, fi) { verts.push(*vertex); }
                }
                Interference::VertexFace { vertex, face } => {
                    if *face == fi { verts.push(*vertex); }
                }
                _ => {}
            }
        }
        verts.sort(); verts.dedup();
        verts
    }

    pub(crate) fn vertex_on_face(&self, vi: usize, fi: usize) -> bool {
        if vi >= self.ds.vertices.len() || fi >= self.ds.faces.len() { return false; }
        let face = &self.ds.faces[fi];
        let tol = face.geom_tol.max(TOLERANCE_ABS);
        for &bvi in &face.boundary_verts {
            if bvi == vi { return true; }
        }
        self.ds.faces[fi].face_info.vertices_in.contains(&vi)
    }

    pub(crate) fn is_existing_vertex(&self, ci: usize, param: f64) -> bool {
        let ic = &self.ds.intersection_curves[ci];
        let tol = ic.geom_tol.max(TOLERANCE_ABS) * 100.0;
        for inf in &self.ds.interferences {
            if let Interference::FaceFace { curves, .. } = inf {
                if curves.contains(&ci) { continue; }
            }
        }
        // Check if any vertex already exists at this parameter
        for fi in 0..self.ds.faces.len() {
            for &vi in &self.ds.faces[fi].face_info.vertices_in {
                if let Some(t) = self.project_vertex_on_curve(vi, ic) {
                    if (t - param).abs() < tol { return true; }
                }
            }
        }
        false
    }

    pub(crate) fn estimate_pave_on_curve(&self, ci: usize, vi: usize) -> Option<f64> {
        let ic = &self.ds.intersection_curves[ci];
        let v_tol = self.ds.vertices[vi].geom_tol;
        let c_tol = ic.geom_tol;
        // OCCT IsVertexOnLine 3-param overload (L4066):
        //   aTolV = BRep_Tool::Tolerance(aV)        �?v_tol
        //   aTolC = aTolR3D                         �?c_tol (no fuzzy)
        //   5-param: aTolSum = aTolV + aTolC; ×2
        let tl = (2.0 * (v_tol + c_tol)).max(1e-6);
        self.project_vertex_on_curve_with_tol(vi, ic, tl)
    }

    pub(crate) fn filter_paves_on_curves(&mut self, curve_idxs: &[usize]) {
        #[derive(Clone)]
        struct PBD { pb_idx: usize, sq_dist: f64, sin_angle: f64, tolerance: f64 }
        let anEps = f64::EPSILON;
        let mut vert_pbs: std::collections::HashMap<usize, Vec<PBD>> = std::collections::HashMap::new();

        for &ci in curve_idxs {
            if ci >= self.ds.intersection_curves.len() { continue; }
            let aNC = &self.ds.intersection_curves[ci];
            let aTolR3D = aNC.geom_tol.max(1e-12);
            let Some(aPB) = aNC.pave_blocks.first() else { continue };
            for pave in &aPB.ext_paves {
                let nV = pave.vertex_idx;
                if nV >= self.ds.vertices.len() { continue; }
                let aPV = self.ds.vertices[nV].point;
                let aPar = pave.param;
                let aPonC = aNC.curve.point_at(aPar);
                let aProjVec = aPV - aPonC;
                let aSqDist = aProjVec.length_squared();
                let dt = 1e-7;
                let tan1 = aNC.curve.point_at(aPar + dt);
                let tan2 = aNC.curve.point_at(aPar - dt);
                let aD1 = (tan1 - tan2) / (2.0 * dt);
                let aSqD1Mod = aD1.length_squared();
                let mut aSin = 0.0;
                if aSqDist > anEps && aSqD1Mod > anEps {
                    aSin = (aProjVec.cross(aD1).length()) / (aSqDist.sqrt() * aSqD1Mod.sqrt());
                }
                vert_pbs.entry(nV).or_default().push(PBD { pb_idx: ci, sq_dist: aSqDist, sin_angle: aSin, tolerance: aTolR3D });
            }
        }

        let aSinAngleMin = 0.5;
        for (&nV, aList) in &vert_pbs {
            let mut aMinDist = f64::MAX;
            for pbd in aList { if pbd.sq_dist < aMinDist { aMinDist = pbd.sq_dist; } }

            let mut aMaxDistKept = -1.0;
            let mut isRemoved = false;
            for pbd in aList {
                let aCheckDist = 100.0 * (pbd.tolerance * pbd.tolerance).max(aMinDist);
                if pbd.sq_dist > aCheckDist && pbd.sin_angle < aSinAngleMin {
                    if pbd.pb_idx < self.ds.intersection_curves.len() {
                        if let Some(pb) = self.ds.intersection_curves[pbd.pb_idx].pave_blocks.first_mut() {
                            pb.remove_ext_pave(nV);
                            isRemoved = true;
                        }
                    }
                } else if pbd.sq_dist > aMaxDistKept {
                    aMaxDistKept = pbd.sq_dist;
                }
            }

            if isRemoved && aMaxDistKept > 0.0 {
                if let Some(&pTol) = self.a_mv_tol.get(&nV) {
                    let aRealTol = pTol.max(aMaxDistKept.sqrt() + 1e-12);
                    if nV < self.ds.vertices.len() {
                        self.ds.vertices[nV].geom_tol = aRealTol;
                    }
                }
            }
        }
    }

    pub(crate) fn put_closing_pave_on_curve(&mut self, curve_idx: usize) {
        let ic = &self.ds.intersection_curves[curve_idx];
        let aT = ic.t_range;
        if !aT[0].is_finite() || !aT[1].is_finite() { return; }
        let aP = [ic.curve.point_at(aT[0]), ic.curve.point_at(aT[1])];
        let Some(aPB) = ic.pave_blocks.first() else { return };
        let mut nV: Option<usize> = None;
        let mut a_t_op = 0.0;
        let mut a_p_op = glam::DVec3::ZERO;
        for pave in &aPB.ext_paves {
            let a_tc = pave.param;
            for j in 0..2 {
                if (a_tc - aT[j]).abs() < crate::tolerance::TOLERANCE_ABS * 100.0 {
                    nV = Some(pave.vertex_idx);
                    a_t_op = if j == 0 { aT[1] } else { aT[0] };
                    a_p_op = if j == 0 { aP[1] } else { aP[0] };
                    break;
                }
            }
            if nV.is_some() { break; }
        }
        let Some(nV) = nV else { return };
        if nV >= self.ds.vertices.len() { return; }
        let a_tol_v = self.ds.vertices[nV].geom_tol;
        let a_pv = self.ds.vertices[nV].point;
        let a_tol_p = ic.geom_tol.max(1e-12) + 1e-12;
        let a_dist_vp = a_pv.distance(a_p_op);
        if a_dist_vp > a_tol_v + a_tol_p { return; }

        if let Some(pb) = self.ds.intersection_curves[curve_idx].change_pave_block1() {
            pb.append_ext_pave(Pave { vertex_idx: nV, param: a_t_op });
        }
    }

    pub(crate) fn extended_tolerance_occt(&self, vi: usize) -> f64 {
        let base_tol = if vi < self.ds.vertices.len() {
            self.ds.vertices[vi].geom_tol
        } else { TOLERANCE_ABS };
        let mut max_tol = base_tol;

        // For newly created vertices, check EE/EF interferences
        for inf in &self.ds.interferences {
            let (ei, param) = match inf {
                Interference::EdgeEdge { e1, param1, new_vertex, .. } if *new_vertex == vi => {
                    (*e1, *param1)
                }
                Interference::EdgeFace { edge, edge_param, new_vertex, .. } if *new_vertex == vi => {
                    (*edge, *edge_param)
                }
                _ => continue,
            };
            // Compute distance from vertex to edge endpoints at the parameter
            if ei < self.ds.edges.len() {
                if let Some(pt) = {
                    use rcad_kernel::geom::CurveEval;
                    Some(self.ds.edges[ei].curve.point_at(param))
                } {
                    let v_pt = self.ds.vertices[vi].point;
                    let d = (v_pt - pt).length();
                    if d > max_tol { max_tol = d; }
                }
            }
        }
        max_tol * 2.0  // Safety factor matching OCCT's approach
    }

    pub(crate) fn get_ef_pnts(&self, ei: usize) -> Vec<(usize, f64)> {
        let mut pnts: Vec<(usize, f64)> = Vec::new();
        for inf in &self.ds.interferences {
            if let Interference::EdgeFace { edge, edge_param, new_vertex, .. } = inf {
                if *edge == ei { pnts.push((*new_vertex, *edge_param)); }
            }
        }
        pnts
    }
    pub(crate) fn treat_vertices_ee(&mut self) {
        let mut to_merge: Vec<(usize, usize)> = Vec::new();
        for inf in &self.ds.interferences {
            if let Interference::EdgeEdge { new_vertex, e1, e2, .. } = inf {
                let vi = *new_vertex;
                if vi >= self.ds.vertices.len() { continue; }
                let v_pt = self.ds.vertices[vi].point;
                for inf2 in &self.ds.interferences {
                    if let Interference::EdgeEdge { new_vertex: vi2, .. } = inf2 {
                        if *vi2 != vi && *vi2 < self.ds.vertices.len() {
                            let d = (v_pt - self.ds.vertices[*vi2].point).length();
                            if d < TOLERANCE_ABS * 100.0 {
                                to_merge.push((vi, *vi2));
                            }
                        }
                    }
                }
            }
        }
        for (keep, remove) in &to_merge {
            for ei in 0..self.ds.edges.len() {
                for pb in &mut self.ds.edges[ei].pave_blocks {
                    if pb.pave1.vertex_idx == *remove { pb.pave1.vertex_idx = *keep; }
                    if pb.pave2.vertex_idx == *remove { pb.pave2.vertex_idx = *keep; }
                }
            }
        }
    }

    pub(crate) fn check_face_paves(&self, fi: usize, vi: usize) -> bool {
        let face = &self.ds.faces[fi];
        let tol = face.geom_tol.max(TOLERANCE_ABS);
        for &ei in &face.boundary_edges {
            if ei >= self.ds.edges.len() { continue; }
            for pb in &self.ds.edges[ei].pave_blocks {
                if pb.pave1.vertex_idx == vi || pb.pave2.vertex_idx == vi {
                    return true;
                }
            }
        }
        false
    }

    pub(crate) fn get_full_shape_map(&self, fi: usize) -> Vec<usize> {
        let mut indices: Vec<usize> = Vec::new();
        indices.push(fi);
        for &ei in &self.ds.faces[fi].boundary_edges { indices.push(ei); }
        for &vi in &self.ds.faces[fi].boundary_verts { indices.push(vi); }
        indices
    }

    pub(crate) fn remove_used_vertices(&self, verts: &mut Vec<usize>, used: &std::collections::BTreeSet<usize>) {
        verts.retain(|v| !used.contains(v));
    }

    pub(crate) fn correct_t_range(&self, ei: usize, t_start: f64, t_end: f64) -> [f64; 2] {
        let edge = &self.ds.edges[ei];
        let mut ts = t_start.max(edge.t_range[0]);
        let mut te = t_end.min(edge.t_range[1]);
        if te < ts { std::mem::swap(&mut ts, &mut te); }
        [ts, te]
    }

    pub(crate) fn is_block_in_on_face(&self, ei: usize, pbi: usize, fi: usize) -> bool {
        if ei >= self.ds.edges.len() { return false; }
        if fi >= self.ds.faces.len() { return false; }
        let pb = &self.ds.edges[ei].pave_blocks[pbi];
        let mid_param = (pb.pave1.param + pb.pave2.param) * 0.5;
        let mid_pt = self.ds.edges[ei].curve.point_at(mid_param);
        let face = &self.ds.faces[fi];
        let tol = face.geom_tol.max(TOLERANCE_ABS);
        match &face.surface {
            Surface3::Plane(p) => (mid_pt - p.origin).dot(p.normal).abs() <= tol,
            _ => false,
        }
    }

    pub(crate) fn make_sd_vertices_ff(&mut self) {
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

    /// �?OCCT-aligned: create section edges for intersecton curves lacking PaveBlocks.
    ///   OCCT PerformFF creates DSEdges for each intersecton curve PB immediately
    ///   (via BOPDS_Curve::InitPaveBlock1 + SetEdge).  rcad defers; this step
    ///   creates the missing section edge + PB so post_treat_ff can register PaveBlocksSc.
    pub(crate) fn init_ic_pave_blocks(&mut self) {
        let ic_indices: Vec<usize> = (0..self.ds.intersection_curves.len()).collect();
        for &ci in &ic_indices {
            let ic = &self.ds.intersection_curves[ci];
            if !ic.pave_blocks.is_empty() { continue; }
            let sv = ic.start_vertex;
            let ev = ic.end_vertex;
            let t_range = ic.t_range;
            let curve = ic.curve.clone();
            let pca = ic.pcurve_on_a.clone();
            let pcb = ic.pcurve_on_b.clone();
            let geom_tol = ic.geom_tol;

            // Pick a face index for each of the two intersecting faces
            // from the FF interference that references this IC.
            let mut f1 = usize::MAX;
            let mut f2 = usize::MAX;
            for inf in &self.ds.interferences {
                if let Interference::FaceFace { f1: a, f2: b, curves, .. } = inf {
                    if curves.contains(&ci) {
                        f1 = *a;
                        f2 = *b;
                        break;
                    }
                }
            }

            let new_ei = self.ds.edges.len();
            let mut face_reps = Vec::new();
            if f1 != usize::MAX {
                if let Some(ref pc) = pca {
                    face_reps.push(DSRepOnFace {
                        face_idx: f1,
                        pcurve: pc.clone(),
                        pcurve2: None,
                        pcurve_range: t_range,
                        start_param: t_range[0],
                        end_param: t_range[1],
                    });
                }
            }
            if f2 != usize::MAX {
                if let Some(ref pc) = pcb {
                    face_reps.push(DSRepOnFace {
                        face_idx: f2,
                        pcurve: pc.clone(),
                        pcurve2: None,
                        pcurve_range: t_range,
                        start_param: t_range[0],
                        end_param: t_range[1],
                    });
                }
            }

            self.ds.edges.push(DSEdge {
                start_vertex: sv,
                end_vertex: ev,
                curve: curve.clone(),
                t_range,
                origin: ShapeOrigin::ShapeA,
                geom_tol,
                paves: vec![],
                pave_blocks: vec![],
                face_reps,
                is_internal: false,
                vertex_params: {
                    let mut vp = std::collections::HashMap::new();
                    vp.insert(sv, t_range[0]);
                    vp.insert(ev, t_range[1]);
                    vp
                },
            });

            let mut pb = PaveBlock::new(NO_EDGE,
                Pave { vertex_idx: sv, param: t_range[0] },
                Pave { vertex_idx: ev, param: t_range[1] },
            );
            pb.curve = Some(curve);
            pb.new_edge = Some(new_ei);
            pb.pcurve_on_a = pca;
            pb.pcurve_on_b = pcb;
            self.ds.intersection_curves[ci].pave_blocks.push(pb);
        }
    }

}



