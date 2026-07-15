use super::*;

impl<'a> super::PaveFiller<'a> {
 /// OCCT BOPAlgo_PaveFiller_8.cxx L54-131: ProcessDE
 pub(crate) fn process_de(&mut self) {
 for an_edge_index in 0..self.ds.nb_source_shapes() {
  if an_edge_index >= self.ds.shape_info.len() { continue; }
  let an_edge_info = &self.ds.shape_info[an_edge_index];
  if an_edge_info.shape_type != rcad_kernel::topods::ShapeType::Edge { continue; }
  let flag = self.ds.edge_flag(an_edge_index);
  if flag <= 0 { continue; }
  // OCCT: HasFlag(nF) sets nF to the shape index directly.
  // rcad convention: the flag stores (face_idx + 1).
  let n_f = (flag - 1) as usize;
  // Map flat face index to shape_info index
  let n_f_si = if n_f < self.ds.face_shape_idx.len() { self.ds.face_shape_idx[n_f] } else { n_f };
  if n_f_si >= self.ds.shape_info.len() { continue; }
  let a_si_f = &self.ds.shape_info[n_f_si];
  // OCCT: nV = anEdgeInfo.SubShapes().First()
  let mut n_v = match an_edge_info.sub_shapes.first() {
   Some(&v) => v,
   None => continue,
  };
  // OCCT: if (HasShapeSD(nV, nVSD)) nV = nVSD
  if let Some(n_v_sd) = self.ds.has_shape_sd(n_v) {
   n_v = n_v_sd;
  }
  if a_si_f.shape_type == rcad_kernel::topods::ShapeType::Face {
   // OCCT: FindPaveBlocks(nV, nF, aLPBOut)
   let a_lpb_out = self.find_pave_blocks(n_v, n_f);
   if !a_lpb_out.is_empty() {
    // OCCT: aLPBD = myDS->ChangePaveBlocks(anEdgeIndex)
    // OCCT: aPBD = aLPBD.First()
    if an_edge_index < self.ds.edges.len() {
     let a_lpbd = &self.ds.edges[an_edge_index].pave_blocks;
     if !a_lpbd.is_empty() {
      let a_pbd = a_lpbd[0].clone();
      // OCCT: FillPaves(nV, anEdgeIndex, nF, aLPBOut, aPBD);
      self.fill_paves(n_v, an_edge_index, n_f, &a_lpb_out, &a_pbd);
      // OCCT: myDS->UpdatePaveBlock(aPBD);
      // rcad: PaveBlock update (shrunk range etc.) is applied during
      // SplitPaveBlocks, not here.
     }
    }
   }
   // OCCT: MakeSplitEdge(anEdgeIndex, nF);
   self.make_split_edge(an_edge_index, n_f);
  } else if a_si_f.shape_type == rcad_kernel::topods::ShapeType::Edge {
   // OCCT: EDGE branch — create a new degenerated edge for the vertex
   // Clone edge data to avoid borrow conflict with self.ds.push_edge
   let (a_curve, a_t_range, a_origin, a_geom_tol, a_is_geometric, a_location) = {
    if an_edge_index < self.ds.edges.len() {
     let e = &self.ds.edges[an_edge_index];
     (e.curve.clone(), e.t_range, e.origin, e.geom_tol, e.is_geometric, e.location)
    } else { continue }
   };
   if n_v >= self.ds.vertices.len() { continue; }
   // Create a new degenerated edge (start == end vertex)
   let de = DSEdge {
    start_vertex: n_v,
    end_vertex: n_v,
    curve: a_curve,
    t_range: a_t_range,
    origin: a_origin,
    geom_tol: a_geom_tol,
    paves: Vec::new(),
    pave_blocks: Vec::new(),
    face_reps: Vec::new(),
    is_internal: false,
    vertex_params: std::collections::HashMap::new(),
    face_tolerances: Vec::new(),
    is_geometric: a_is_geometric,
    location: a_location,
   };
   // OCCT: nEn = myDS->Append(aSI);
   let n_en = self.ds.push_edge(de, None);
   // OCCT: aPBD->SetEdge(nEn)
   if let Some(a_pbd) = self.ds.edges[an_edge_index].pave_blocks.first() {
    a_pbd.0.write().unwrap().new_edge = Some(n_en);
   }
  }
 }
 }

 pub(crate) fn fill_shrunk_data(&mut self) {
  // OCCT BOPAlgo_PaveFiller_9.cxx L65-138: FillShrunkData
  // Collect all edge PBs
  let mut all_pb: Vec<(usize, usize, usize, f64, f64)> = Vec::new();
  for ei in 0..self.ds.edges.len() {
   if self.ds.edge_has_flag(ei) || self.ds.is_edge_degenerated(ei) { continue; }
   let n_pbs = self.ds.edges[ei].pave_blocks.len();
   if n_pbs == 0 { continue; }
   for pi in 0..n_pbs {
    let pb = &self.ds.edges[ei].pave_blocks[pi];
    let r = pb.0.read().unwrap();
    if r.shrunk_range.is_some() { continue; }
    all_pb.push((ei, r.pave1.vertex_idx, r.pave2.vertex_idx, r.pave1.param, r.pave2.param));
   }
  }
  if all_pb.is_empty() { return; }
  let ec: Vec<Curve3> = self.ds.edges.iter().map(|e| e.curve.clone()).collect();
  let et: Vec<f64> = self.ds.edges.iter().map(|e| e.geom_tol).collect();
  let v_tols: Vec<f64> = (0..self.ds.vertex_origins.len()).map(|vi| self.ds.vertex_tolerance(vi)).collect();
  let results: Vec<(Option<[f64; 2]>, bool)> = all_pb.iter().map(|(ei, v1i, v2i, p1, p2)| {
   if self.ds.is_edge_degenerated(*ei) { return (None, false); }
   let mut sr = crate::inttools::shrunk_range::ShrunkRange::new();
   sr.set_data(*ei, [*p1, *p2], v_tols[*v1i], v_tols[*v2i], et[*ei]);
   sr.perform(&ec[*ei]);
   (sr.shrunk_range(), sr.is_splittable())
  }).collect();
  for ((ei, _v1i, _v2i, _p1, _p2), (range, splittable)) in all_pb.iter().zip(results.iter()) {
   if let Some(pb) = self.ds.edges.get(*ei).and_then(|e| e.pave_blocks.first()) {
    pb.0.write().unwrap().shrunk_range = *range;
    pb.0.write().unwrap().is_splittable = *splittable;
   }
  }
 }

 /// BOPAlgo_PaveFiller::AnalyzeShrunkData (PaveFiller_3.cxx L766-824).
 pub(crate) fn analyze_shrunk_data(&mut self, ei: usize, pb_idx: usize) {
     if pb_idx >= self.ds.pave_blocks.len() { return; }
     let (has_shrunk, is_splittable) = {
         let pb = &self.ds.pave_blocks[pb_idx];
         let r = pb.0.read().unwrap();
         (r.shrunk_range.is_some(), r.is_splittable)
     };
     if !has_shrunk || !is_splittable {
         let pb = &self.ds.pave_blocks[pb_idx];
         pb.0.write().unwrap().shrunk_range = Some([0.0, 0.0]);
         pb.0.write().unwrap().is_splittable = false;
     }
 }

 pub(crate) fn existing_pave_block(&self, ei: usize, vi: usize) -> bool {
 for pb in &self.ds.edges[ei].pave_blocks {
 if pb.0.read().unwrap().pave1.vertex_idx == vi || pb.0.read().unwrap().pave2.vertex_idx == vi { return true; }
 }
 false
 }

  /// SplitPaveBlocks (PaveFiller_2.cxx L419-626).
  /// Splits PaveBlocks of the given edges using extra paves.
  /// For each edge, reconstructs pave_blocks from the full pave set,
  /// computes shrunk data, and unifies vertices when no valid range.
    /// BOPAlgo_PaveFiller::SplitPaveBlocks (BOPAlgo_PaveFiller_2.cxx L419-626).
  pub(crate) fn split_pave_blocks(&mut self, edges: &std::collections::HashSet<usize>,
  add_interfs: bool) {
    // OCCT L422: fence map
    let mut a_mpairs: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    // OCCT L424-427: aMCBNewPB + aMVerticesToInitPB
    let mut a_mcb_new_pb: std::collections::HashMap<usize,
      Vec<(usize, PaveBlock)>> = std::collections::HashMap::new();
    let mut a_m_vertices_to_init_pb: Vec<usize> = Vec::new();

    // OCCT L432: for each edge in theMEdges
    for &n_e in edges {
      // OCCT L435-436: nE, aLPB = ChangePaveBlocks
      if n_e >= self.ds.edge_start_vertex.len() { continue; }
      if self.ds.edge_paves(n_e).len() < 2 { continue; }

      // OCCT L449: aCB = CommonBlock(aPB) — rcad: map old PB pairs to CB index
      let old_pair_to_cb: std::collections::HashMap<(usize, usize), usize> =
        self.ds.edge_pave_blocks(n_e).iter()
        .filter_map(|pb| {
          pb.0.read().unwrap().common_block_idx.map(|cb| {
            let v1 = pb.0.read().unwrap().pave1.vertex_idx;
            let v2 = pb.0.read().unwrap().pave2.vertex_idx;
            let key = if v1 <= v2 { (v1, v2) } else { (v2, v1) };
            (key, cb)
          })
        })
        .collect();

      // OCCT L452-453: aPB->Update(aLPBN) — split using ext_paves (rcad: edge.paves)
      let mut params: Vec<(usize, f64)> = self.ds.edge_paves(n_e).iter()
        .map(|p| (p.vertex_idx, p.param)).collect();
      params.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
      params.dedup_by(|a, b| (a.1 - b.1).abs() < TOLERANCE_ABS);

      // OCCT L457-524: for each new sub-PB in aLPBN
      let mut new_pbs: Vec<PaveBlock> = Vec::new();
      for w in params.windows(2) {
        let a_pb_n = PaveBlock::new(n_e,
          Pave { vertex_idx: w[0].0, param: w[0].1 },
          Pave { vertex_idx: w[1].0, param: w[1].1 });

        // OCCT L461: UpdatePaveBlockWithSDVertices (rcad: SD handled separately)

        // OCCT L462: FillShrunkData(aPBN)
        let v1_tol = self.ds.vertex_tolerance(w[0].0);
        let v2_tol = self.ds.vertex_tolerance(w[1].0);
        let mut a_sr = crate::inttools::shrunk_range::ShrunkRange::new();
        a_sr.set_data(n_e, [w[0].1, w[1].1], v1_tol, v2_tol,
                      self.ds.edge_tolerance(n_e));
        a_sr.perform(&self.ds.edge_curve(n_e).unwrap());

        // OCCT L464: bHasValidRange = HasShrunkData()
        let b_has_valid_range = a_sr.shrunk_range().is_some();
        // OCCT L468: bCheckDist = (bHasValidRange && !IsSplittable())
        let b_check_dist = b_has_valid_range && !a_sr.is_splittable();

        // OCCT L469: if (!bHasValidRange || bCheckDist)
        if !b_has_valid_range || b_check_dist {
          let (n_v1, n_v2) = (w[0].0, w[1].0);
          // OCCT L473-476: if (nV1 == nV2) continue
          if n_v1 == n_v2 { continue; }
          // OCCT L480-488: if (bCheckDist) { ComputeVV }
          let b_interfere = if b_check_dist {
            let tol = self.ds.vertex_tolerance(n_v1).max(self.ds.vertex_tolerance(n_v2))
                       .max(self.ds.fuzzy_tol);
            self.ds.vertex_point(n_v1).distance(self.ds.vertex_point(n_v2)) <= tol
          } else { false };
          // OCCT L491: if (!bHasValidRange) OR (interfering after ComputeVV)
          if !b_has_valid_range || b_interfere {
            let key = if n_v1 < n_v2 { (n_v1, n_v2) } else { (n_v2, n_v1) };
            // OCCT L495: if (aMPairs.Add(aPair))
            if a_mpairs.insert(key) {
              // OCCT L497-504: MakeSDVertices + aMVerticesToInitPB
              let n_v = n_v1.min(n_v2);
              self.ds.add_shape_sd(n_v1, n_v);
              self.ds.add_shape_sd(n_v2, n_v);
              if add_interfs {
                self.ds.interf_vv.push(InterferenceVV {
                  v1: n_v1, v2: n_v2, merged_vertex: n_v,
                });
              }
              a_m_vertices_to_init_pb.push(n_v1);
              a_m_vertices_to_init_pb.push(n_v2);
            }
            continue; // OCCT L506: skip this sub-PB (invalid range or merged)
          }
          // OCCT L509: falls through → valid PB, keep it (bCheckDist and NOT interfering)
        }
        // OCCT L511: aLPB.Append(aPBN)
        // OCCT L513-523: if CB, store to aMCBNewPB
        let (v1, v2) = (a_pb_n.pave1.vertex_idx, a_pb_n.pave2.vertex_idx);
        let new_key = if v1 <= v2 { (v1, v2) } else { (v2, v1) };
        if let Some(&cb_idx) = old_pair_to_cb.get(&new_key) {
          a_mcb_new_pb.entry(cb_idx).or_default().push((n_e, a_pb_n.clone()));
        }
        new_pbs.push(a_pb_n);
      }
      // OCCT L525-527: replace old PBs
      self.ds.edge_pave_blocks_mut(n_e).clone_from(
        &new_pbs.into_iter().map(crate::bopds::pave::SharedPB::new).collect::<Vec<_>>()
      );
    }

    // ===== OCCT L530-618: Make Common Blocks =====
    // (rcad: logic preserved from prior implementation)
    let mut a_m_inds: std::collections::HashMap<(usize, (usize, usize)), Vec<(usize, usize)>> =
      std::collections::HashMap::new();
    for (&cb_idx, pb_entries) in &a_mcb_new_pb {
      for &(ei, ref pb) in pb_entries {
        let (v1, v2) = (pb.pave1.vertex_idx, pb.pave2.vertex_idx);
        let key = if v1 <= v2 { (v1, v2) } else { (v2, v1) };
        if let Some(li) = self.ds.edges.get(ei)
          .and_then(|e| e.pave_blocks.iter().position(|p| {
            let (pv1, pv2) = (p.0.read().unwrap().pave1.vertex_idx,
                               p.0.read().unwrap().pave2.vertex_idx);
            (pv1 == v1 && pv2 == v2) || (pv1 == v2 && pv2 == v1)
          }))
        { a_m_inds.entry((cb_idx, key)).or_default().push((ei, li)); }
      }
    }
    // OCCT L559: for each group in aMInds
    for ((_, _vk), entries) in &a_m_inds {
      if entries.len() < 2 { continue; }
      let mut ue: Vec<usize> = entries.iter().map(|&(ei, _)| ei).collect();
      ue.sort(); ue.dedup();
      if ue.len() < 2 { continue; }
      let mut pfcb: Vec<(usize, usize)> = Vec::new();
      for &(ei, li) in entries {
        if let Some(pb) = self.ds.edges.get(ei).and_then(|e| e.pave_blocks.get(li)) {
          let gix = self.ds.pave_blocks.len();
          self.ds.pave_blocks.push(pb.clone());
          let oei = pb.0.read().unwrap().original_edge;
          for (fi, f) in self.ds.faces.iter().enumerate() {
            if f.boundary_edges.contains(&oei)
              || f.inner_boundary_edges.iter().any(|w| w.iter().any(|(e, _)| *e == oei))
            { pfcb.push((gix, fi)); }
          }
        }
      }
      pfcb.sort(); pfcb.dedup();
      if pfcb.len() < 2 { continue; }
      let mut cb = crate::bopds::common_block::CommonBlock::new();
      let mut uf: Vec<usize> = pfcb.iter().map(|&(_, fi)| fi).collect();
      uf.sort(); uf.dedup();
      for &(gix, fi) in &pfcb { cb.add_pave_block(gix, fi); }
      cb.set_faces(uf);
      let cb_idx = self.ds.common_blocks.len();
      self.ds.common_blocks.push(cb);
      for &(ei, li) in entries {
        if let Some(lp) = self.ds.edges.get_mut(ei).and_then(|e| e.pave_blocks.get_mut(li))
        { lp.0.write().unwrap().common_block_idx = Some(cb_idx); }
      }
    }

    // OCCT L620-625: InitPaveBlocksForVertex
    for &vi in &a_m_vertices_to_init_pb {
      self.ds.init_pave_blocks_for_vertex(vi);
    }

  }

 pub(crate) fn split_ics_at_periodic_boundary(&mut self) {
 let n_curves = self.ds.intersection_curves.len();
 let mut curve_faces = std::collections::HashMap::new();
 for ff in &self.ds.interf_ff {
 for &ci in &ff.curves {
 curve_faces.entry(ci).or_insert((ff.f1, ff.f2));
 }
 }
 for ci in 0..n_curves {
 let Some(&(fa, fb)) = curve_faces.get(&ci) else { continue };
 let t0 = self.ds.intersection_curves[ci].t_range[0];
 let t1 = self.ds.intersection_curves[ci].t_range[1];
 if (t1 - t0).abs() < TOLERANCE_ABS { continue; }
 // Clone the curve to avoid borrow conflict when push_vertex is needed later.
 let ic_curve = self.ds.intersection_curves[ci].curve.clone();
 for &(fi, use_a) in &[(fa, true), (fb, false)] {
 let (period, _, _) = match &self.ds.faces[fi].surface {
 Surface3::Sphere(_) => (std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
 Surface3::Cylinder(_) => (std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
 _ => continue,
 };
 let pcurve_owned = if use_a { self.ds.intersection_curves[ci].pcurve_on_a.clone() } else { self.ds.intersection_curves[ci].pcurve_on_b.clone() };
 let Some(pc) = pcurve_owned else { continue };
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
 let pt = ic_curve.point_at(t_cross);
 let v_new = match self.ds.find_vertex_near(pt, TOLERANCE_ABS * 1000.0) {
 Some(v) => v,
 None => {
 let vi = self.ds.vertices.len();
 self.ds.push_vertex(crate::bopds::ds::DSVertex {
 point: pt, geom_tol: TOLERANCE_ABS, origin: None,
 is_internal: false, location: 0,
 }, None);
 vi
 }
 };
 self.ds.face_info_mut(fa).vertices_in.insert(v_new);
 self.ds.face_info_mut(fb).vertices_in.insert(v_new);
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
 let aV_tol = self.ds.vertex_tolerance(nV);
 let mut aTolV = self.a_mv_tol.get(&nV).copied().unwrap_or(aV_tol);

 let mut aT: f64 = 0.0;
 let mut bIsVertexOnLine = self.is_vertex_on_line(nV, aTolV, curve_idx, aTolR3D + self.fuzzy_tolerance, &mut aT);
 if std::env::var("RCAD_DEBUG_MB").is_ok() {
  eprintln!("[PPC] nV={} aTolR3D={:.6e} tol_sum={:.6e} isOnLine={} aT={:.6e}",
   nV, aTolR3D, aTolR3D + self.fuzzy_tolerance, bIsVertexOnLine, aT);
 }

 if !bIsVertexOnLine && iCheckExtend != 0 && !self.verts_to_avoid_extension.contains(&nV) {
 let mut anExtraTol = aTolV;
 if self.extended_tolerance(nV, aMI, &mut anExtraTol, iCheckExtend) {
 bIsVertexOnLine = self.is_vertex_on_line(nV, anExtraTol, curve_idx, aTolR3D + self.fuzzy_tolerance, &mut aT);
 if bIsVertexOnLine {
 let aPOnC = ic_curve.point_at(aT);
 aTolV = aPOnC.distance(self.ds.vertex_point(nV));
 }
 }
 }

 if bIsVertexOnLine {
 if std::env::var("RCAD_DEBUG_MB").is_ok() {
  eprintln!("[PPC]  VERTEX_ON_LINE nV={} aT={:.6e}", nV, aT);
 }
 let aDTol = CONFUSION * 1e-5; // OCCT BOPTools_AlgoTools::DTolerance()
 let aPTol = Self::curve_parametric_tolerance(&ic_curve, aTolR3D.max(aTolV));

  let mut nVUsed = 0;
  // Clone the SharedPB outside the if-let to avoid keeping self.ds borrowed
  let shared_pb = self.ds.intersection_curves[curve_idx].pave_blocks.first().cloned();
  if let Some(spb) = shared_pb {
  let mut pb = spb.0.write().unwrap();
 let bExist = pb.contains_parameter(aT, aPTol, &mut nVUsed);
 if bExist {
 let pList = self.a_dmv_lv.entry(nVUsed).or_insert_with(|| {
 let mut list = Vec::new();
 list.push(nVUsed);
 if !self.a_mv_tol.contains_key(&nVUsed) {
 let aVUsed_tol = if nVUsed < self.ds.vertices.len() {
 self.ds.vertex_tolerance(nVUsed)
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
 let aP2 = self.ds.vertex_point(nV);
 let aDist = aP1.distance(aP2);
 if aTolV < aDist + aDTol {
 if !self.a_mv_tol.contains_key(&nV) {
 self.a_mv_tol.insert(nV, aTolV);
 }
 // Tolerance update (lock is dropped with spb scope)
 let new_tol = aDist + aDTol;
 self.ds.vertex_data_mut(nV).tolerance = new_tol;
 }
 }
 }
 return Some(aT);
 }
 if std::env::var("RCAD_DEBUG_MB").is_ok() {
  eprintln!("[PPC]  VERTEX_NOT_ON_LINE nV={}", nV);
 }
 None
 }

 /// BOPAlgo_PaveFiller::PutPavesOnCurve (L2404-2453).
 /// Takes pre-computed vertex sets instead of re-deriving from interferences.
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
 let c_box = crate::pave_filler::helpers::curve_bounding_box_simple(&ic.curve, 0.0);
 for &nV in theMVEF {
 self.put_pave_on_curve(nV, aTolR3D, curve_idx, aMI, 2);
 }
 for &nV in theMVOnIn {
 if theMVEF.contains(&nV) { continue; }
 if theMVCommon.contains(&nV) {
 self.put_pave_on_curve(nV, aTolR3D, curve_idx, aMI, 1);
 } else {
 if nV < self.ds.vertex_shape_idx.len() {
 let nv_si = self.ds.vertex_shape_idx[nV];
 if nv_si < self.ds.shape_info.len() {
 if let (Some(vmn), Some(vmx)) = (self.ds.shape_info[nv_si].box_min, self.ds.shape_info[nv_si].box_max) {
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
 }
 if nV < self.ds.vertex_shape_idx.len() && !self.ds.shape_info[self.ds.vertex_shape_idx[nV]].is_new {
 continue;
 }
 self.put_pave_on_curve(nV, aTolR3D, curve_idx, aMI, 1);
 }
 }
 }

 /// PutStickPavesOnCurve (BOPAlgo_PaveFiller_6.cxx L2780-2875).
 /// Processes stick vertices that are near the IC endpoints (rich criterion)
 /// and where both face normals are nearly opposite (crease criterion).
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
 if start_vertex < self.ds.vertices.len() && end_vertex < self.ds.vertices.len() {
 return;
 }
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
 let v_pt = self.ds.vertex_point(n_v);
 for m in 0..2 {
 if (m == 0 && start_vertex < self.ds.vertices.len()) ||
 (m == 1 && end_vertex < self.ds.vertices.len()) {
 continue;
 }
 let d2 = a_p[m].distance_squared(v_pt);
 if d2 > a_dt2 { continue; }
 let mut sc_pr = 1.0;
 if let (Some(s0), Some(s1)) = (&surf_both[0], &surf_both[1]) {
 let n0 = Self::estimate_surface_normal(s0, a_p[m]);
 let n1 = Self::estimate_surface_normal(s1, a_p[m]);
 sc_pr = n0.dot(n1);
 if sc_pr < 0.0 { sc_pr = -sc_pr; }
 sc_pr = 1.0 - sc_pr;
 }
 if sc_pr > a_d_sc_pr { continue; }
 let a_d = d2.sqrt();
 let a_tol_r3d_use = a_tol_r3d.max(self.ds.vertex_tolerance(n_v));
 self.put_pave_on_curve(n_v, a_d.min(a_tol_r3d_use), ci, aMI, 1);
 break;
 }
 }
 }

 /// PutEFPavesOnCurve (PaveFiller_6.cxx L2724-2778).
 /// For single-curve FF pairs, put EF-intersection vertices onto the curve.
 pub(crate) fn put_ef_paves_on_curve(
  &mut self,
  ci: usize,
  aMI: &std::collections::HashSet<usize>,
  aMVEF: &std::collections::HashSet<usize>,
 ) {
  let a_tol_r3d = {
   let ic = &self.ds.intersection_curves[ci];
   ic.geom_tol.max(ic.curve_extra.tangential_tol)
  };
  for &n_v in aMVEF {
   if n_v < self.ds.vertices.len() {
    self.put_pave_on_curve(n_v, a_tol_r3d, ci, aMI, 1);
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

 /// IntTools_Context::IsVertexOnLine (L786-992).
 /// Form-identical logic:
 ///   1. aTolSum = curve-type-adjusted (aTolV + aTolC)
 ///   2. First endpoint  ?local projection (Newton from aFirst),
 ///      then global fallback (closest_point_on_curve).
 ///   3. Last endpoint  ?same with bFirstValid shortcut.
 ///   4. Global projection  ?GeomAPI_ProjectPointOnCurve style
 ///      with NbPoints==0 fallback to bounded-curve endpoints.
 pub(crate) fn is_vertex_on_line(
 &self,
 nV: usize,
 aTolV: f64,
 curve_idx: usize,
 aTolC: f64,
 aT: &mut f64,
 ) -> bool {
 let ic = &self.ds.intersection_curves[curve_idx];
 let vp = self.ds.vertex_point(nV);
 // L800: aTolSum = aTolV + aTolC
 let mut aTolSum = aTolV + aTolC;

 // L802-819: curve-type-based tolerance adjustment
 let is_bspline_or_bezier = matches!(&ic.curve, Curve3::BSpline(_) | Curve3::Bezier(_));
 if is_bspline_or_bezier {
 aTolSum = 2.0 * aTolSum;
 if aTolSum < 1e-5 { aTolSum = 1e-5; }
 } else {
 aTolSum = 2.0 * aTolSum;
 if aTolSum < TOLERANCE_MESH_LEGACY { aTolSum = TOLERANCE_MESH_LEGACY; }
 }

 // L821-822: aFirst, aLast from curve parameter range
 let aFirst = ic.t_range[0];
 let aLast = ic.t_range[1];

 let mut bFirstValid = false;
 let mut aFirstDist = f64::MAX;

 // ── L824-883: First endpoint check ──────────────────────────────
 if !precision_is_infinite(aFirst) {
 let p_first = ic.curve.point_at(aFirst);
 aFirstDist = vp.distance(p_first);
 if aFirstDist < aTolSum {
 bFirstValid = true;
 *aT = aFirst;
 if aFirstDist > aTolV {
 // L840: local projection from aFirst (Newton refinement)
 let mut t_local = aFirst;
 let mut dist_local = aFirstDist;
 let clamp = |t: f64| t.clamp(aFirst, aLast);
 let dt = 1e-7;
 let max_iter = 20;
 let mut converged = false;
 for _ in 0..max_iter {
 let p = ic.curve.point_at(t_local);
 let diff = p - vp;
 let deriv = ic.curve.derivative_at(t_local);
 let deriv_sq = deriv.dot(deriv);
 if deriv_sq < 1e-20 { break; }
 let curv = (ic.curve.point_at(t_local + 2.0 * dt) - 2.0 * p
 + ic.curve.point_at(t_local - 2.0 * dt)) / (dt * dt);
 let denom = deriv_sq + diff.dot(curv);
 let delta = diff.dot(deriv) / if denom.abs() > 1e-20 { denom } else { deriv_sq };
 let new_t = clamp(t_local - delta);
 let new_dist = (ic.curve.point_at(new_t) - vp).length();
 if new_dist < dist_local {
 dist_local = new_dist;
 t_local = new_t;
 } else {
 break;
 }
 if delta.abs() < TOLERANCE_LINEAR_ULTRA_STRICT { converged = true; break; }
 }
 // L842-851: validate local projection result
 if converged || dist_local < aFirstDist {
 let mid = (aLast + aFirst) * 0.5;
 let p_proj = ic.curve.point_at(t_local);
 if (t_local > mid) || (dist_local > aTolSum)
 || (p_first.distance(p_proj) < TOLERANCE_ABS) {
 *aT = aFirst;
 } else {
 *aT = t_local;
 }
 } else {
 // L856-880: local projection failed  ?global fallback
 use rcad_kernel::projection::closest_point_on_curve;
 let proj = closest_point_on_curve(&ic.curve, vp, 64);
 let mid = (aLast + aFirst) * 0.5;
 if (proj.param > mid) || (proj.distance > aTolSum)
 || (p_first.distance(proj.point) < TOLERANCE_ABS) {
 *aT = aFirst;
 } else {
 *aT = proj.param;
 }
 }
 }
 }
 }

 // ── L886-951: Last endpoint check ───────────────────────────────
 if !precision_is_infinite(aLast) {
 let p_last = ic.curve.point_at(aLast);
 let d_last = vp.distance(p_last);
 // L890-893: bFirstValid shortcut
 if bFirstValid && aFirstDist < d_last {
 return true;
 }
 if d_last < aTolSum {
 *aT = aLast;
 if d_last > aTolV {
 // L901: local projection from aLast (Newton refinement)
 let mut t_local = aLast;
 let mut dist_local = d_last;
 let clamp = |t: f64| t.clamp(aFirst, aLast);
 let dt = 1e-7;
 let max_iter = 20;
 let mut converged = false;
 for _ in 0..max_iter {
 let p = ic.curve.point_at(t_local);
 let diff = p - vp;
 let deriv = ic.curve.derivative_at(t_local);
 let deriv_sq = deriv.dot(deriv);
 if deriv_sq < 1e-20 { break; }
 let curv = (ic.curve.point_at(t_local + 2.0 * dt) - 2.0 * p
 + ic.curve.point_at(t_local - 2.0 * dt)) / (dt * dt);
 let denom = deriv_sq + diff.dot(curv);
 let delta = diff.dot(deriv) / if denom.abs() > 1e-20 { denom } else { deriv_sq };
 let new_t = clamp(t_local - delta);
 let new_dist = (ic.curve.point_at(new_t) - vp).length();
 if new_dist < dist_local {
 dist_local = new_dist;
 t_local = new_t;
 } else {
 break;
 }
 if delta.abs() < TOLERANCE_LINEAR_ULTRA_STRICT { converged = true; break; }
 }
 // L905-911: validate local projection result
 if converged || dist_local < d_last {
 let mid = (aLast + aFirst) * 0.5;
 let p_proj = ic.curve.point_at(t_local);
 if (t_local < mid) || (dist_local > aTolSum)
 || (p_last.distance(p_proj) < TOLERANCE_ABS) {
 *aT = aLast;
 } else {
 *aT = t_local;
 }
 } else {
 // L916-940: local projection failed  ?global fallback
 use rcad_kernel::projection::closest_point_on_curve;
 let proj = closest_point_on_curve(&ic.curve, vp, 64);
 let mid = (aLast + aFirst) * 0.5;
 if (proj.param < mid) || (proj.distance > aTolSum)
 || (p_last.distance(proj.point) < TOLERANCE_ABS) {
 *aT = aLast;
 } else {
 *aT = proj.param;
 }
 }
 }
 return true;
 }
 } else if bFirstValid {
 return true;
 }

 // ── L953-992: General projection (GeomAPI_ProjectPointOnCurve) ──
 use rcad_kernel::projection::closest_point_on_curve;
 let proj = closest_point_on_curve(&ic.curve, vp, 64);
 let aNbProj = if proj.distance.is_finite() { 1 } else { 0 };
 if aNbProj == 0 {
 // L959-978: bounded curve fallback  ?check start/end points
 if matches!(&ic.curve, Curve3::BSpline(_) | Curve3::Bezier(_))
 || matches!(&ic.curve, Curve3::Line(_)) {
 let aPDist = vp.distance(ic.curve.point_at(aFirst));
 if aPDist < aTolSum {
 *aT = aFirst;
 return true;
 }
 let aPDist = vp.distance(ic.curve.point_at(aLast));
 if aPDist < aTolSum {
 *aT = aLast;
 return true;
 }
 }
 return false;
 }
 // L983-991: check distance vs tolerance
 let aDist = proj.distance;
 if aDist > aTolSum {
 return false;
 }
 *aT = proj.param;
 return true;
 }

 /// BOPAlgo_PaveFiller::ExtendedTolerance (PaveFiller_6.cxx L????).
 /// In-out aTolVExt: on input it is the current vertex tolerance (lower bound),
 /// on output it is set to the maximum extended distance found.
 pub(crate) fn extended_tolerance(
 &self,
 nV: usize,
 aMI: &std::collections::HashSet<usize>,
 aTolVExt: &mut f64,
 aType: i32,
 ) -> bool {
 if nV < self.ds.vertex_shape_idx.len() && !self.ds.shape_info[self.ds.vertex_shape_idx[nV]].is_new {
 return false;
 }
 let vp = self.ds.vertex_point(nV);
 let mut found = false;
 // Use input aTolVExt as initial threshold (OCCT in-out semantics)
 let mut max_ext = *aTolVExt;

 if aType == 0 || aType == 1 {
 for inf in &self.ds.interf_ee {
 if inf.new_vertex != nV { continue; }
 if !aMI.contains(&inf.e1) { continue; }
 if inf.e1 < self.ds.edges.len() {
 let p1 = crate::boptools::point_on_edge(&self.ds.edges[inf.e1], inf.param1);
 let p2 = crate::boptools::point_on_edge(&self.ds.edges[inf.e1], inf.param2);
 let d = vp.distance(p1).max(vp.distance(p2));
 if d > max_ext { max_ext = d; }
 found = true;
 }
 }
 }
 if aType == 0 || aType == 2 {
 for inf in &self.ds.interf_ef {
 if inf.new_vertex != nV { continue; }
 if !aMI.contains(&inf.edge) { continue; }
 if inf.edge < self.ds.edges.len() {
 let p1 = crate::boptools::point_on_edge(&self.ds.edges[inf.edge], inf.edge_param);
 let d = vp.distance(p1);
 if d > max_ext { max_ext = d; }
 found = true;
 }
 }
 }
 if found { *aTolVExt = max_ext; }
 found
 }

 pub(crate) fn project_vertex_on_curve(&self, vi: usize, ic: &IntersectionCurve) -> Option<f64> {
 use rcad_kernel::geom::CurveEval;
 let v_tol = self.ds.vertex_tolerance(vi);
 let c_tol = ic.geom_tol;
 let f_tol = self.ds.fuzzy_tol;
 // OCCT IsVertexOnLine L800: aTolSum = aTolV + aTolC
 // where aTolC = aTolR3D + myFuzzyValue (PutPaveOnCurve L2976)
 let raw_sum = v_tol + c_tol + f_tol;
 let tl = (2.0 * raw_sum).max(TOLERANCE_MESH_LEGACY);
 let result = self.project_vertex_on_curve_with_tol(vi, ic, tl);
 if result.is_some() { return result; }
 let ext_tol = self.extended_tolerance_occt(vi);
 if ext_tol > tl { self.project_vertex_on_curve_with_tol(vi, ic, ext_tol) } else { None }
 }

 pub(crate) fn project_vertex_on_curve_with_tol(&self, vi: usize, ic: &IntersectionCurve, tl: f64) -> Option<f64> {
 use rcad_kernel::geom::CurveEval;
 let vp = self.ds.vertex_point(vi);
 match &ic.curve {
 Curve3::Line(l) => { let v = vp - l.origin; let t = v.dot(l.direction);
 if (v - l.direction*t).length() <= tl { Some(t.clamp(ic.t_range[0], ic.t_range[1])) } else { None } }
 Curve3::Circle(c) => {
 let v = vp - c.center;
 let tl_cap = tl.min(c.radius.max(TOLERANCE_ABS) * 1e-4);
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
 // section face requires pcurves.
 if self.avoid_build_pcurve
 || (!self.section_attribute.pcurve_on_s1 && !self.section_attribute.pcurve_on_s2)
 {
 return;
 }
 let b_pcurve_on_s = [self.section_attribute.pcurve_on_s1, self.section_attribute.pcurve_on_s2];
 // rcad: the MPC (Make PCurve) concept is inlined  ?we store the data
 // needed to either (a) compute a new pcurve, or (b) update vertex tolerances
 // from an already-existing pcurve.
 #[derive(Clone)]
 struct MPCEntry {
 edge_idx: usize,
 face_idx: usize,
 /// OCCT: SetFlag(true) -> call UpdateVertices after pcurve is set.
 update_vertices: bool,
 /// OCCT SetData: existing edge whose pcurve on this face is reused.
 existing_edge: Option<usize>,
 }

 let mut a_vmpc: Vec<MPCEntry> = Vec::new();

 // ===== Phase 1: CommonBlock edges (OCCT L603-700) =====
 // Process PaveBlocksIn (boundary edges through face) and PaveBlocksOn
 // (edges lying on the face surface). If a CommonBlock peer has the pcurve,
 // SetData to copy it rather than recompute.
 for fi in 0..self.ds.faces.len() {
  let face_info = &self.ds.faces[fi].face_info;
  // OCCT L608-618: PaveBlocksIn
  for &pb_idx in &face_info.pave_blocks_in {
   if pb_idx >= self.ds.pave_blocks.len() { continue; }
   let pb = &self.ds.pave_blocks[pb_idx];
   let ei = pb.0.read().unwrap().original_edge;
   if ei >= self.ds.edges.len() { continue; }
   a_vmpc.push(MPCEntry {
    edge_idx: ei,
    face_idx: fi,
    update_vertices: false,
    existing_edge: None,
   });
  }
  // OCCT L620-700: PaveBlocksOn - check CB peers for existing pcurve
  for &pb_idx in &face_info.pave_blocks_on {
   if pb_idx >= self.ds.pave_blocks.len() { continue; }
   let pb = &self.ds.pave_blocks[pb_idx];
   let ei = pb.0.read().unwrap().original_edge;
   if ei >= self.ds.edges.len() { continue; }
   // OCCT L625: if pcurve already exists on this face -> skip
   if self.ds.edge_face_reps(ei).iter().any(|r| r.face_idx == fi) {
    continue;
   }
   // OCCT L630-676: check CommonBlock peer edges for existing pcurve
   let mut cb_existing_edge: Option<usize> = None;
   if let Some(cb_idx) = pb.0.read().unwrap().common_block_idx {
    if let Some(cb) = self.ds.common_blocks.get(cb_idx) {
     let cb_pbs = cb.pave_blocks();
     if cb_pbs.len() >= 2 {
      for &(other_pb_idx, _other_fi) in cb_pbs {
       if other_pb_idx == pb_idx { continue; }
       if let Some(other_pb) = self.ds.pave_blocks.get(other_pb_idx) {
        let other_ei = other_pb.0.read().unwrap().original_edge;
        if let Some(other_edge) = self.ds.edges.get(other_ei) {
         if other_edge.face_reps.iter().any(|r| r.face_idx == fi) {
          // OCCT L670: AttachExistingPCurve(aEx)
          cb_existing_edge = Some(other_ei);
          break;
         }
        }
       }
      }
     }
    }
   }
   // OCCT L678-692: push MPC entry (with or without existing_edge)
   a_vmpc.push(MPCEntry {
    edge_idx: ei,
    face_idx: fi,
    update_vertices: false,
    existing_edge: cb_existing_edge,
   });
  }
 }

 // ===== Phase 2: Section edges (OCCT L702-757) =====
 // Section edges are created from FF intersection curves. Their pcurves
 // are already stored in the edge's face_reps during make_blocks.
 // This phase creates MPC entries with update_vertices=true so that vertex
 // tolerances are corrected from the 2D<->3D deviation.
 // OCCT L706: iterate FF interferences
 if b_pcurve_on_s[0] || b_pcurve_on_s[1] {
  let mut ef_pairs: std::collections::HashSet<(usize, usize)> =
   std::collections::HashSet::new();
  for ff in &self.ds.interf_ff {
   let nf = [ff.f1, ff.f2];
   for &ci in &ff.curves {
    if ci >= self.ds.intersection_curves.len() { continue; }
    let ic = &self.ds.intersection_curves[ci];
    // OCCT L712: iterate PBs of the intersection curve
    for spb in &ic.pave_blocks {
     // OCCT L713: section edges have original_edge == NO_EDGE
     let section_ei = match spb.0.read().unwrap().new_edge {
      Some(ei) => ei,
      None => continue,
     };
     // OCCT L720-728: for each face of the FF pair
     for (m, &fi) in nf.iter().enumerate() {
      if !b_pcurve_on_s[m] { continue; }
      // OCCT L722: dedup pair (section_ei, fi)
      if !ef_pairs.insert((section_ei, fi)) {
       continue;
      }
      // OCCT L726: create MPC entry with SetFlag(true)
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

 // ===== Phase 3: Perform MPC computations (OCCT L760-775) =====
 // Each MPC entry either copies an existing pcurve (from CommonBlock peer)
 // or computes a new one via projection. The UpdateVertices list accumulates
 // edges whose vertex tolerances need correction from 2D<->3D deviation.
 let mut failed_pcurves: Vec<(usize, usize)> = Vec::new();
 let mut applied_pcurves: Vec<(usize, usize, Curve2d, f64, f64)> = Vec::new();
 let mut update_vertices_list: Vec<(usize, usize)> = Vec::new();

 for mpc in &a_vmpc {
 let ei = mpc.edge_idx;
 let fi = mpc.face_idx;
 if ei >= self.ds.edges.len() || fi >= self.ds.faces.len() {
 continue;
 }

 // OCCT L762: copy from existing edge in same CommonBlock
 if let Some(src_ei) = mpc.existing_edge {
 if src_ei < self.ds.edges.len() {
 if let Some(rep) = self.ds.edge_face_reps(src_ei).iter().find(|r| r.face_idx == fi) {
 applied_pcurves.push((
 ei, fi,
 rep.pcurve.clone(),
 self.ds.edge_range(ei)[0],
 self.ds.edge_range(ei)[1],
 ));
 if mpc.update_vertices {
 update_vertices_list.push((ei, fi));
 }
 continue;
 }
 }
 }

 // OCCT L768: HasCurveOnSurface check
 if self.ds.edge_face_reps(ei).iter().any(|r| r.face_idx == fi) {
 if mpc.update_vertices {
 update_vertices_list.push((ei, fi));
 }
 continue;
 }

 // OCCT L770: compute pcurve from edge curve + face surface
 let face_surface = self.ds.faces[fi].surface.clone();
 let edge_curve = &self.ds.edges[ei].curve;
 if let Some((pcurve, _len)) = DS::compute_edge_pcurve(edge_curve, &face_surface, None) {
 applied_pcurves.push((
 ei, fi, pcurve,
 self.ds.edge_range(ei)[0],
 self.ds.edge_range(ei)[1],
 ));
 if mpc.update_vertices {
 update_vertices_list.push((ei, fi));
 }
 } else {
 // OCCT L780: pcurve computation failed
 failed_pcurves.push((ei, fi));
 }
 }

 // ===== Phase 4: Apply results (OCCT L782-877) =====
 // OCCT L782-789: error reporting for failed pcurves
 for &(ei, fi) in &failed_pcurves {
 // OCCT: BRep_Builder MakeCompound + Add(Edge) + Add(Face) +
 //        AddWarning(new BOPAlgo_AlertBuildingPCurveFailed(compound))
 // rcad: non-fatal warning
 if std::env::var("RCAD_DEBUG_PF").is_ok() {
 eprintln!("[DBG] MakePCurves failed for edge={} on face={}", ei, fi);
 }
 }

 // OCCT L792-840: assign computed pcurves to edge face_reps
 use crate::bopds::ds::DSCurveRepOnFace;
 for (ei, fi, pcurve, t0, t1) in applied_pcurves {
 self.ds.edges[ei].face_reps.push(DSCurveRepOnFace {
 face_idx: fi,
 pcurve,
 pcurve2: None,
 pcurve_range: [t0, t1],
 start_param: t0,
 end_param: t1,
 });
 }

 // OCCT L842-877: UpdateVertices - correct vertex tolerances from
 // 2D pcurve deviation. For each vertex of the edge, compute
 // the 3D point from the pcurve and compare to the vertex position.
 for &(ei, fi) in &update_vertices_list {
  if ei >= self.ds.edges.len() || fi >= self.ds.faces.len() { continue; }
  let [a_t1, a_t2] = self.ds.edge_range(ei);
  let a_c3d = match self.ds.edge_curve(ei) {
   Some(c) => c.clone(),
   None => continue,
  };
  let face_reps: Vec<_> = self.ds.edge_face_reps(ei).iter().map(|r| r.clone()).collect();
  let face_surface = match self.ds.face_surface(fi) {
   Some(s) => s.clone(),
   None => continue,
  };
  let pc2d = face_reps.iter().find(|r| r.face_idx == fi).map(|r| &r.pcurve);
  let start_v = self.ds.edge_start_vertex_ds(ei);
  let end_v = self.ds.edge_end_vertex_ds(ei);
  if let Some(pc) = pc2d {
   for (j, a_t) in [a_t1, a_t2].iter().enumerate() {
    let a_p3d = a_c3d.point_at(*a_t);
    let vi = if j == 0 { start_v } else { end_v };
    if vi < self.ds.vertices.len() {
     let a_tol_v = self.ds.vertex_tolerance(vi);
     let uv = pc.point_at(*a_t);
     // OCCT: project 2D->3D using face surface
     let proj = match &face_surface {
      Surface3::Plane(p) => {
       p.origin + uv.x * p.u_dir + uv.y * p.v_dir
      }
      Surface3::Cylinder(c) => {
       let theta = uv.x / c.radius;
       let h = uv.y;
       let ref_dir = c.ref_dir.normalize_or_zero();
       let cross_dir = c.axis.cross(ref_dir).normalize_or_zero();
       c.origin + ref_dir * (c.radius * theta.cos())
           + cross_dir * (c.radius * theta.sin())
           + c.axis * h
      }
      _ => a_p3d,
     };
     let d2 = a_p3d.distance_squared(proj);
     if d2 > a_tol_v * a_tol_v {
      // OCCT L870: SetTolerance(aTolV + aCmp)
      let new_tol = d2.sqrt() + rcad_kernel::tolerance::COMPUTATIONAL;
      self.ds.vertex_data_mut(vi).tolerance = new_tol;
     }
    }
   }
  }
 }
 }
 pub(crate) fn update_face_info(&mut self, fi: usize) {
 let face = &self.ds.faces[fi].clone();
 let mut sc_curves: Vec<usize> = face.face_info.curves_sc.iter().copied().collect();
 // Recompute vertices_in from all IC endpoints on this face
 for &ci in &sc_curves {
 if ci < self.ds.intersection_curves.len() {
 let ic = &self.ds.intersection_curves[ci];
 let sv = ic.start_vertex;
 let ev = ic.end_vertex;
 self.ds.face_info_mut(fi).vertices_in.insert(sv);
 self.ds.face_info_mut(fi).vertices_in.insert(ev);
 }
 }
 // Update pave_blocks_in/pave_blocks_on from edge pave blocks
 for &ei in &face.boundary_edges {
 if ei < self.ds.edges.len() {
 let pbs = self.ds.edges[ei].pave_blocks.clone();
 for pb in &pbs {
 let pb_idx = self.ds.pave_blocks.len();
 self.ds.pave_blocks.push(pb.clone());
 self.ds.face_info_mut(fi).pave_blocks_in.insert(pb_idx);
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
  for spb in &mut self.ds.edges[ei].pave_blocks {
  let mut pb = spb.0.write().unwrap();
  if pb.pave1.vertex_idx == vi { pb.pave1.vertex_idx = sd_vi; }
  if pb.pave2.vertex_idx == vi { pb.pave2.vertex_idx = sd_vi; }
 }
 }
 for fi in 0..self.ds.faces.len() {
 if self.ds.face_info(fi).vertices_in.contains(&vi) {
 self.ds.face_info_mut(fi).vertices_in.remove(&vi);
 self.ds.face_info_mut(fi).vertices_in.insert(sd_vi);
 }
 }
 }
 }
 }
 pub(crate) fn get_stick_vertices(&self, fi: usize) -> Vec<usize> {
 let mut verts: Vec<usize> = Vec::new();
 for inf in &self.ds.interf_vv {
 if self.vertex_on_face(inf.v1, fi) { verts.push(inf.v1); }
 if self.vertex_on_face(inf.v2, fi) { verts.push(inf.v2); }
 }
 for inf in &self.ds.interf_ve {
 if self.vertex_on_face(inf.vertex, fi) { verts.push(inf.vertex); }
 }
 for inf in &self.ds.interf_vf {
 if inf.face == fi { verts.push(inf.vertex); }
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
 self.ds.face_info(fi).vertices_in.contains(&vi)
 }

 /// IsExistingVertex (PaveFiller_6.cxx L1200-1245).
 /// Checks if a 3D point `theP` coincides with any vertex in theMVOnIn
 /// within `theTol + vertex.geom_tol` distance.
 pub(crate) fn is_existing_vertex_at_point(
  &self,
  the_p: glam::DVec3,
  the_tol: f64,
  the_mv_on_in: &std::collections::HashSet<usize>,
 ) -> bool {
  for &n_v in the_mv_on_in {
   if n_v >= self.ds.vertices.len() { continue; }
   let a_pnt = self.ds.vertex_point(n_v);
   let a_tol_v = self.ds.vertex_tolerance(n_v);
   if a_pnt.distance(the_p) < the_tol + a_tol_v {
    return true;
   }
  }
  false
 }

 pub(crate) fn is_existing_vertex(&self, ci: usize, param: f64) -> bool {
 let ic = &self.ds.intersection_curves[ci];
 let tol = ic.geom_tol.max(TOLERANCE_ABS) * 100.0;
 // (FaceFace curve check removed  ?was a no-op loop over interferences)
 // Check if any vertex already exists at this parameter
 for fi in 0..self.ds.faces.len() {
 for &vi in &self.ds.face_info(fi).vertices_in {
 if let Some(t) = self.project_vertex_on_curve(vi, ic) {
 if (t - param).abs() < tol { return true; }
 }
 }
 }
 false
 }

 pub(crate) fn estimate_pave_on_curve(&self, ci: usize, vi: usize) -> Option<f64> {
 let ic = &self.ds.intersection_curves[ci];
 let v_tol = self.ds.vertex_tolerance(vi);
 let c_tol = ic.geom_tol;
 // OCCT IsVertexOnLine 3-param overload (L4066):
 // aTolV = BRep_Tool::Tolerance(aV) ?v_tol
 // aTolC = aTolR3D ?c_tol (no fuzzy)
 // 5-param: aTolSum = aTolV + aTolC;  2
 let tl = (2.0 * (v_tol + c_tol)).max(TOLERANCE_MESH_LEGACY);
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
 let aTolR3D = aNC.geom_tol.max(CONFUSION);
 let Some(aPB) = aNC.pave_blocks.first() else { continue };
 for pave in &aPB.0.read().unwrap().ext_paves {
 let nV = pave.vertex_idx;
 if nV >= self.ds.vertices.len() { continue; }
 let aPV = self.ds.vertex_point(nV);
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
  pb.0.write().unwrap().remove_ext_pave(nV);
 isRemoved = true;
 }
 }
 } else if pbd.sq_dist > aMaxDistKept {
 aMaxDistKept = pbd.sq_dist;
 }
 }

 if isRemoved && aMaxDistKept > 0.0 {
 if let Some(&pTol) = self.a_mv_tol.get(&nV) {
 let aRealTol = pTol.max(aMaxDistKept.sqrt() + CONFUSION);
 if nV < self.ds.vertices.len() {
 self.ds.vertex_data_mut(nV).tolerance = aRealTol;
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
 for pave in &aPB.0.read().unwrap().ext_paves {
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
 let a_tol_v = self.ds.vertex_tolerance(nV);
 let a_pv = self.ds.vertex_point(nV);
 let a_tol_p = ic.geom_tol.max(1e-12) + 1e-12;
 let a_dist_vp = a_pv.distance(a_p_op);
 if a_dist_vp > a_tol_v + a_tol_p { return; }

 if let Some(mut pb) = self.ds.intersection_curves[curve_idx].pave_blocks.first().map(|spb| spb.0.write().unwrap()) {
  pb.append_ext_pave(Pave { vertex_idx: nV, param: a_t_op });
 }
 }

 pub(crate) fn extended_tolerance_occt(&self, vi: usize) -> f64 {
 let base_tol = if vi < self.ds.vertices.len() {
 self.ds.vertex_tolerance(vi)
 } else { TOLERANCE_ABS };
 let mut max_tol = base_tol;

 // For newly created vertices, check EE/EF interferences
 for inf in &self.ds.interf_ee {
 if inf.new_vertex == vi {
 if inf.e1 < self.ds.edges.len() {
 if let Some(pt) = {
 use rcad_kernel::geom::CurveEval;
 Some(self.ds.edges[inf.e1].curve.point_at(inf.param1))
 } {
 let v_pt = self.ds.vertex_point(vi);
 let d = (v_pt - pt).length();
 if d > max_tol { max_tol = d; }
 }
 }
 }
 }
 for inf in &self.ds.interf_ef {
 if inf.new_vertex == vi {
 if inf.edge < self.ds.edges.len() {
 if let Some(pt) = {
 use rcad_kernel::geom::CurveEval;
 Some(self.ds.edges[inf.edge].curve.point_at(inf.edge_param))
 } {
 let v_pt = self.ds.vertex_point(vi);
 let d = (v_pt - pt).length();
 if d > max_tol { max_tol = d; }
 }
 }
 }
 }
 max_tol * 2.0  // Safety factor matching OCCT's approach
 }

 pub(crate) fn get_ef_pnts(&self, ei: usize) -> Vec<(usize, f64)> {
 let mut pnts: Vec<(usize, f64)> = Vec::new();
 for inf in &self.ds.interf_ef {
 if inf.edge == ei { pnts.push((inf.new_vertex, inf.edge_param)); }
 }
 pnts
 }
 // OCCT BOPAlgo_PaveFiller_4.cxx L305-390: TreatVerticesEE
 pub(crate) fn treat_vertices_ee(&mut self) {
  // OCCT L307-313: collect unique new vertex indices from EE interferences
  let mut a_liv: Vec<usize> = Vec::new();
  let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
  for aee in &self.ds.interf_ee {
   let n_v = aee.new_vertex;
   if a_mi.insert(n_v) {
    a_liv.push(n_v);
   }
  }
  if a_liv.is_empty() { return; }
  // OCCT L316-321: collect face indices
  let mut a_lif: Vec<usize> = Vec::new();
  for n_f in 0..self.ds.faces.len() {
   a_lif.push(n_f);
  }
  if a_lif.is_empty() { return; }
  // OCCT L324-387: cross-product iterate faces × EE vertices
  for &n_f in &a_lif {
   for &n_v in &a_liv {
    // OCCT L328: if (!aFI.VerticesOn().Contains(nV))
    if self.ds.face_info(n_f).vertices_on.contains(&n_v) { continue; }
    // OCCT L330-332: ComputeVF (aV, aF, aT1, aT2, dummy, myFuzzyValue)
    let tf = self.vf_tol(n_v, n_f);
    if let Ok(vf_result) = self.context.compute_vf(self.ds, n_v, n_f, tf) {
     // OCCT L334-335: aVF.SetIndices(nV, nF); aVF.SetUV(aT1, aT2);
     // OCCT L336: myDS->AddInterf(nV, nF);
     self.ds.try_add_interf(n_v, n_f);
     self.ds.interf_vf.push(InterferenceVF {
      vertex: n_v, face: n_f,
      u: vf_result.u, v: vf_result.v,
      index_new: None,
     });
     // OCCT L337: aFI.ChangeVerticesIn().Add(nV);
     self.ds.face_info_mut(n_f).vertices_in.insert(n_v);
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
 if pb.0.read().unwrap().pave1.vertex_idx == vi || pb.0.read().unwrap().pave2.vertex_idx == vi {
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
 let pb = &self.ds.edge_pave_blocks(ei)[pbi];
 let mid_param = (pb.0.read().unwrap().pave1.param + pb.0.read().unwrap().pave2.param) * 0.5;
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

 ///  create section edges for intersecton curves lacking PaveBlocks.
 /// OCCT PerformFF creates DSEdges for each intersecton curve PB immediately
 /// (via BOPDS_Curve::InitPaveBlock1 + SetEdge).  rcad defers; this step
 /// creates the missing section edge + PB so post_treat_ff can register PaveBlocksSc.
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
 if let Some(ff) = self.ds.interf_ff.iter().find(|ff| ff.curves.contains(&ci)) {
 f1 = ff.f1;
 f2 = ff.f2;
 }

 let new_ei = self.ds.edges.len();
 let mut face_reps = Vec::new();
 if f1 != usize::MAX {
 if let Some(ref pc) = pca {
 face_reps.push(DSCurveRepOnFace {
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
 face_reps.push(DSCurveRepOnFace {
 face_idx: f2,
 pcurve: pc.clone(),
 pcurve2: None,
 pcurve_range: t_range,
 start_param: t_range[0],
 end_param: t_range[1],
 });
 }
 }

 self.ds.push_edge(DSEdge {
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
  face_tolerances: Vec::new(),
  is_geometric: true,
  location: 0,
  }, None);

 let mut pb = PaveBlock::new(NO_EDGE,
 Pave { vertex_idx: sv, param: t_range[0] },
 Pave { vertex_idx: ev, param: t_range[1] },
 );
 pb.curve = Some(curve);
  pb.new_edge = Some(new_ei);
 pb.pcurve_on_a = pca;
 pb.pcurve_on_b = pcb;
 self.ds.intersection_curves[ci].pave_blocks.push(crate::bopds::pave::SharedPB::new(pb));
 }
 }

 /// BOPAlgo_PaveFiller::SplitEdge (PaveFiller_7.cxx).
 /// Splits an edge at parametric position, creating a new DS vertex and edge.
 pub(crate) fn split_edge(&mut self, n_e: usize, n_v1: usize, _a_t1: f64, n_v2: usize, _a_t2: f64) -> usize {
     if n_e >= self.ds.edges.len() { return n_e; }
     // Clone all needed data before mutating self.ds
     let (curve, t_range, origin, geom_tol, is_geometric, location) = {
         let e = &self.ds.edges[n_e];
         (e.curve.clone(), e.t_range, e.origin, e.geom_tol, e.is_geometric, e.location)
     };
     let de = DSEdge {
         start_vertex: n_v1,
         end_vertex: n_v2,
         curve,
         t_range,
         origin,
         geom_tol,
         paves: Vec::new(),
         pave_blocks: Vec::new(),
         face_reps: Vec::new(),
         is_internal: false,
         vertex_params: std::collections::HashMap::new(),
         face_tolerances: Vec::new(),
         is_geometric,
         location,
     };
     let new_ei = self.ds.push_edge(de, None);
     self.ds.edge_data_mut(new_ei).tolerance = geom_tol;
     new_ei
 }

 /// BOPAlgo_PaveFiller::FindPaveBlocks (PaveFiller_8.cxx).
 /// Find pave blocks that contain the given vertex on the given face.
 pub(crate) fn find_pave_blocks(&self, n_v: usize, n_f: usize) -> Vec<usize> {
     let mut result = Vec::new();
     if n_f >= self.ds.faces.len() { return result; }
     let fi = &self.ds.faces[n_f];
     for &pb_idx in fi.face_info.pave_blocks_on.iter()
         .chain(fi.face_info.pave_blocks_in.iter())
         .chain(fi.face_info.pave_blocks_sc.iter())
     {
         if pb_idx >= self.ds.pave_blocks.len() { continue; }
         let pb = &self.ds.pave_blocks[pb_idx];
         let (nv1, nv2) = pb.0.read().unwrap().indices();
         if nv1 == n_v || nv2 == n_v {
             result.push(pb_idx);
         }
     }
     result
 }

 /// BOPAlgo_PaveFiller::MakeSplitEdge (PaveFiller_8.cxx).
 /// Creates a new split edge for the vertex on face.
 pub(crate) fn make_split_edge(&mut self, n_v: usize, n_f: usize) {
     // In OCCT this creates a section edge for a vertex on a face.
     // For now, just ensure the vertex is registered in the face info.
     if n_v < self.ds.vertices.len() && n_f < self.ds.faces.len() {
         self.ds.face_info_mut(n_f).vertices_in.insert(n_v);
     }
 }

 /// OCCT BOPAlgo_PaveFiller_8.cxx L224-331: FillPaves.
 /// Adds intersection points between the degenerated edge's 2D curve
 /// and each passing edge's 2D curve as Extra paves.
 pub(crate) fn fill_paves(&mut self, _n_vd: usize, _n_ed: usize, _n_fd: usize,
     _a_lpb_out: &[usize], _a_pbd: &SharedPB) {
     // OCCT: Ensures the vertex is added as an Extra Pave on the PaveBlock
     // of the degenerated edge. Implementation uses Geom2dInt_GInter to
     // intersect the 2D curves and calls AddSplitPoint for valid results.
     // rcad: full implementation pending 2D curve intersection support.
 }

}



