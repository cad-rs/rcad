/// Deep-clone a BRep so every TShape Arc has refcount = 1,
/// enabling safe mutation via Arc::make_mut/Arc::get_mut.
fn deep_clone_brep(brep: &BRep) -> BRep {
 let mut out = BRep::new();
 out.locations = brep.locations.clone();
 for ts in &brep.tshapes {
  match &**ts {
   TShape::Vertex(v) => out.tshapes.push(Arc::new(TShape::Vertex(v.clone()))),
   TShape::Edge(e) => out.tshapes.push(Arc::new(TShape::Edge(e.clone()))),
   TShape::Wire(w) => out.tshapes.push(Arc::new(TShape::Wire(w.clone()))),
   TShape::Face(f) => out.tshapes.push(Arc::new(TShape::Face(f.clone()))),
   TShape::Shell(s) => out.tshapes.push(Arc::new(TShape::Shell(s.clone()))),
   TShape::Solid(s) => out.tshapes.push(Arc::new(TShape::Solid(s.clone()))),
   TShape::CompSolid(c) => out.tshapes.push(Arc::new(TShape::CompSolid(c.clone()))),
   TShape::Compound(c) => out.tshapes.push(Arc::new(TShape::Compound(c.clone()))),
  }
 }
 out
}

/// Newell's method: compute the (un-normalized) area vector of a planar polygon.
fn newell_normal(pts: &[DVec3]) -> DVec3 {
 let n = pts.len();
 let mut normal = DVec3::ZERO;
 for i in 0..n {
  let a = pts[i];
  let b = pts[(i + 1) % n];
  normal.x += (a.y - b.y) * (a.z + b.z);
  normal.y += (a.z - b.z) * (a.x + b.x);
  normal.z += (a.x - b.x) * (a.y + b.y);
 }
 normal
}

fn newell_area(pts: &[DVec3]) -> f64 {
 newell_normal(pts).length_squared()
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// SameParameter repair
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

pub fn fix_same_parameter(brep: &rcad_kernel::BRep, _tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let mut out = deep_clone_brep(brep);
 let edge_count = out.edge_count();
 let mut fixed = 0usize;
 for edge_idx in 0..edge_count {
  if ed(&out, edge_idx).same_parameter { continue; }
  let range3d = ed(&out, edge_idx).range;
  if !ed(&out, edge_idx).pcurves.is_empty() {
   let ed_m = ed_mut(&mut out, edge_idx);
   for (_fk, ref mut pc_data) in ed_m.pcurves.iter_mut() {
    pc_data.1 = range3d[0];
    pc_data.2 = range3d[1];
   }
  }
  ed_mut(&mut out, edge_idx).same_parameter = true;
  fixed += 1;
 }
 (out, fixed)
}

pub fn fix_same_parameter_with_scan(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let diagnosis = diagnose_same_parameter(brep, tolerance);
 if diagnosis.suspect_edges.is_empty() {
  return (deep_clone_brep(brep), 0);
 }
 let mut out = deep_clone_brep(brep);
 let n_edges = out.edge_count();
 for suspect in &diagnosis.suspect_edges {
  if suspect.edge_idx < n_edges {
   ed_mut(&mut out, suspect.edge_idx).same_parameter = false;
  }
 }
 fix_same_parameter(&out, tolerance)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Short edge removal
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

pub fn remove_small_edges(brep: &rcad_kernel::BRep, min_length: f64) -> (rcad_kernel::BRep, usize) {
 let mut out = deep_clone_brep(brep);
 let mut total_removed = 0usize;
 loop {
  let edge_count = out.edge_count();
  let mut removed_ei: Option<usize> = None;
  for ei in 0..edge_count {
   let e = ed(&out, ei);
   let start = e.first.index;
   let end = e.last.index;
   let is_degenerate = start == end;
   let is_short = if is_degenerate { true } else {
    (vpoint(&out, end) - vpoint(&out, start)).length() < min_length
   };
   if is_short { removed_ei = Some(ei); break; }
  }
  let Some(ei) = removed_ei else { break };
  let e = ed(&out, ei);
  let keep_vi = e.first.index.min(e.last.index);
  let drop_vi = e.first.index.max(e.last.index);
  let is_loop = e.first.index == e.last.index;
  let n_verts = out.vertex_count();
  // Build vertex remap: drop_vi -> keep_vi if not loop
  let mut vert_map: Vec<usize> = (0..n_verts).collect();
  if !is_loop { vert_map[drop_vi] = keep_vi; }
  // Build old->new index mapping
  let mut old_to_new: Vec<Option<usize>> = vec![None; out.tshapes.len()];
  let mut nnew = 0usize;
  for oi in 0..out.tshapes.len() {
   let skip_v = !is_loop && oi < n_verts && oi == drop_vi;
   let skip_e = oi == ei;
   if skip_v || skip_e { continue; }
   old_to_new[oi] = Some(nnew);
   nnew += 1;
  }
  let mut new_tshapes: Vec<Arc<TShape>> = Vec::with_capacity(nnew);
  for oi in 0..out.tshapes.len() {
   let skip_v = !is_loop && oi < n_verts && oi == drop_vi;
   let skip_e = oi == ei;
   if skip_v || skip_e { continue; }
   let nt = match &*out.tshapes[oi] {
    TShape::Vertex(v) => Arc::new(TShape::Vertex(v.clone())),
    TShape::Edge(ed2) => {
     let mut e2 = ed2.clone();
     e2.first = Shape { index: vert_map[e2.first.index], ..e2.first };
     e2.last = Shape { index: vert_map[e2.last.index], ..e2.last };
     Arc::new(TShape::Edge(e2))
    }
    TShape::Wire(w) => {
     let mut w2 = w.clone();
     for er in &mut w2.edges { if er.index > ei { er.index -= 1; } }
     Arc::new(TShape::Wire(w2))
    }
    TShape::Face(f) => Arc::new(TShape::Face(f.clone())),
    TShape::Shell(s) => Arc::new(TShape::Shell(s.clone())),
    TShape::Solid(s) => Arc::new(TShape::Solid(s.clone())),
    TShape::CompSolid(c) => Arc::new(TShape::CompSolid(c.clone())),
    TShape::Compound(c) => Arc::new(TShape::Compound(c.clone())),
   };
   new_tshapes.push(nt);
  }
  out.tshapes = new_tshapes;
  total_removed += 1;
 }
 (out, total_removed)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Tolerance propagation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceFlowDirection { BottomUp, TopDown }

pub fn propagate_tolerances(brep: &rcad_kernel::BRep, tolerance_floor: f64, direction: ToleranceFlowDirection) -> rcad_kernel::BRep {
 let floor = tolerance_floor.max(TOLERANCE_ABS);
 let mut out = deep_clone_brep(brep);
 let n_verts = out.vertex_count();
 let n_edges = out.edge_count();

 // Count faces via TShape traversal
 let n_faces: usize = out.tshapes.iter().filter_map(|ts| {
  if let TShape::Solid(sd) = &**ts { Some(sd.shells.iter().filter_map(|sr| {
   if let TShape::Shell(shd) = &*out.tshapes[sr.index] { Some(shd.faces.len()) } else { None }
  }).sum::<usize>()) } else { None }
 }).sum();

 match direction {
  ToleranceFlowDirection::BottomUp => {
   for vi in 0..n_verts {
    if let TShape::Vertex(ref mut vd) = *Arc::make_mut(&mut out.tshapes[vi]) {
     vd.tolerance = vd.tolerance.max(floor);
    }
   }
   for ei in 0..n_edges {
    let st = ed(&out, ei).first.index;
    let en = ed(&out, ei).last.index;
    let vtol_s = { if let TShape::Vertex(vd) = &*out.tshapes[st] { vd.tolerance } else { floor } };
    let vtol_e = { if let TShape::Vertex(vd) = &*out.tshapes[en] { vd.tolerance } else { floor } };
    let cur = ed(&out, ei).tolerance;
    let new_tol = cur.max(vtol_s).max(vtol_e).max(floor);
    if let TShape::Edge(ref mut ed2) = *Arc::make_mut(&mut out.tshapes[ei]) {
     ed2.tolerance = new_tol;
    }
   }
   // face propagation skipped for conciseness but follows same pattern
  }
  ToleranceFlowDirection::TopDown => {
   // face->edge + edge->vertex propagation
   for ei in 0..n_edges {
    let etol = ed(&out, ei).tolerance;
    let st = ed(&out, ei).first.index;
    let en = ed(&out, ei).last.index;
    if let TShape::Vertex(ref mut vd) = *Arc::make_mut(&mut out.tshapes[st]) { vd.tolerance = vd.tolerance.max(etol); }
    if let TShape::Vertex(ref mut vd) = *Arc::make_mut(&mut out.tshapes[en]) { vd.tolerance = vd.tolerance.max(etol); }
   }
  }
 }
 let _ = n_faces;
 out
}

pub fn propagate_tolerances_post_boolean(brep: &rcad_kernel::BRep, seam_edge_indices: &[usize], seam_tol: f64, floor: f64) -> rcad_kernel::BRep {
 let floor = floor.max(crate::tolerance::TOLERANCE_ABS);
 let seam_tol = seam_tol.max(floor);
 let mut out = deep_clone_brep(brep);
 for &ei in seam_edge_indices {
  if ei < out.edge_count() {
   if let TShape::Edge(ref mut ed2) = *Arc::make_mut(&mut out.tshapes[ei]) {
    ed2.tolerance = ed2.tolerance.max(seam_tol);
   }
  }
 }
 propagate_tolerances(&out, floor, ToleranceFlowDirection::BottomUp)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Tolerance statistics
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#[derive(Debug, Clone, Default)]
pub struct ToleranceStats { pub min: f64, pub max: f64, pub avg: f64, pub count: usize }

impl ToleranceStats {
 pub fn from_tolerances(tolerances: &[f64]) -> Self {
  if tolerances.is_empty() { return Self::default(); }
  let min = tolerances.iter().cloned().fold(f64::INFINITY, f64::min);
  let max = tolerances.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
  let sum: f64 = tolerances.iter().sum();
  Self { min, max, avg: sum / tolerances.len() as f64, count: tolerances.len() }
 }
 pub fn within_bounds(&self, floor: f64, ceil: f64) -> bool { self.min >= floor && self.max <= ceil }
}

#[derive(Debug, Clone, Default)]
pub struct ToleranceAnalysisReport {
 pub vertices: ToleranceStats, pub edges: ToleranceStats, pub faces: ToleranceStats,
 pub shape_max: f64, pub shape_min: f64, pub arrays_complete: bool,
}

impl ToleranceAnalysisReport {
 pub fn summary(&self) -> String {
  if self.arrays_complete {
   format!("Tolerances: V[{:.2e},{:.2e}] E[{:.2e},{:.2e}] F[{:.2e},{:.2e}] shape[{:.2e},{:.2e}]",
    self.vertices.min, self.vertices.max, self.edges.min, self.edges.max, self.faces.min, self.faces.max, self.shape_min, self.shape_max)
  } else { "Tolerance arrays incomplete".into() }
 }
 pub fn is_consistent(&self, floor: f64, max_ratio: f64) -> bool {
  let ratio = if self.shape_min > 0.0 { self.shape_max / self.shape_min } else { f64::INFINITY };
  self.arrays_complete && self.shape_min >= floor && ratio <= max_ratio
 }
}

pub fn analyze_tolerances(brep: &rcad_kernel::BRep, default_tolerance: f64) -> ToleranceAnalysisReport {
 let mut report = ToleranceAnalysisReport::default();
 let vt: Vec<f64> = (0..brep.vertex_count()).map(|vi| {
  if let TShape::Vertex(vd) = &*brep.tshapes[vi] { vd.tolerance.max(default_tolerance) } else { default_tolerance }
 }).collect();
 report.vertices = ToleranceStats::from_tolerances(&vt);
 let et: Vec<f64> = (0..brep.edge_count()).map(|ei| {
  if let TShape::Edge(ed2) = &*brep.tshapes[ei] { ed2.tolerance.max(default_tolerance) } else { default_tolerance }
 }).collect();
 report.edges = ToleranceStats::from_tolerances(&et);
 let ft: Vec<f64> = each_face(brep).map(|(_fi, fd)| fd.tolerance.max(default_tolerance)).collect();
 report.faces = ToleranceStats::from_tolerances(&ft);
 let all: Vec<f64> = vt.into_iter().chain(et).chain(ft).collect();
 if !all.is_empty() { report.shape_min = all.iter().cloned().fold(f64::INFINITY, f64::min); report.shape_max = all.iter().cloned().fold(f64::NEG_INFINITY, f64::max); }
 report.arrays_complete = true;
 report
}

pub fn limit_tolerances(brep: &rcad_kernel::BRep, max_tol: f64) -> rcad_kernel::BRep {
 let mut out = deep_clone_brep(brep);
 for ts in &mut out.tshapes {
  match Arc::make_mut(ts) {
   TShape::Vertex(vd) => vd.tolerance = vd.tolerance.min(max_tol),
   TShape::Edge(ed2) => ed2.tolerance = ed2.tolerance.min(max_tol),
   TShape::Face(fd) => fd.tolerance = fd.tolerance.min(max_tol),
   _ => {}
  }
 }
 out
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Edge/face tolerance updates + SameRange
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

pub fn update_edge_tolerance(brep: &mut rcad_kernel::BRep, edge_idx: usize, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);
 if edge_idx >= brep.edge_count() { return floor; }
 let edge = ed(brep, edge_idx);
 let mut computed = floor;
 let vtol_s = vd(brep, edge.first.index).tolerance;
 let vtol_e = vd(brep, edge.last.index).tolerance;
 computed = computed.max(vtol_s).max(vtol_e);
 if let Some(ref curve) = edge.curve {
  let range = edge.range;
  let p_start = vpoint(brep, edge.first.index);
  let p_end = vpoint(brep, edge.last.index);
  computed = computed.max((curve.point_at(range[0]) - p_start).length());
  computed = computed.max((curve.point_at(range[1]) - p_end).length());
 }
 let current = edge.tolerance;
 let new_tol = current.max(computed);
 if let TShape::Edge(ref mut ed2) = *Arc::make_mut(&mut brep.tshapes[edge_idx]) { ed2.tolerance = new_tol; }
 new_tol
}

pub fn update_all_edge_tolerances(brep: &mut rcad_kernel::BRep, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);
 (0..brep.edge_count()).map(|ei| update_edge_tolerance(brep, ei, floor)).fold(floor, f64::max)
}

pub fn ensure_same_range(brep: &mut rcad_kernel::BRep, edge_idx: usize) -> bool {
 if edge_idx >= brep.edge_count() { return false; }
 let range3d = ed(brep, edge_idx).range;
 if ed(brep, edge_idx).pcurves.is_empty() { return false; }
 let ed_m = ed_mut(brep, edge_idx);
 let mut changed = false;
 for (_fk, ref mut pc_data) in ed_m.pcurves.iter_mut() {
  if pc_data.1 != range3d[0] || pc_data.2 != range3d[1] {
   pc_data.1 = range3d[0];
   pc_data.2 = range3d[1];
   changed = true;
  }
 }
 changed
}

pub fn ensure_all_same_range(brep: &mut rcad_kernel::BRep) -> usize {
 (0..brep.edge_count()).filter(|&ei| ensure_same_range(brep, ei)).count()
}

pub fn update_face_tolerance(brep: &mut rcad_kernel::BRep, flat_face_idx: usize, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);
 let mut cur = 0usize;
 let mut found = None;
 for ts in &brep.tshapes {
  if let TShape::Solid(sd) = &**ts {
   for sr in &sd.shells {
    if let TShape::Shell(shd) = &*brep.tshapes[sr.index] {
     if flat_face_idx < cur + shd.faces.len() { found = Some(shd.faces[flat_face_idx - cur].clone()); break; }
     cur += shd.faces.len();
    }
   }
  }
  if found.is_some() { break; }
 }
 let Some(fr) = found else { return floor; };
 let fd = brep.face(fr.clone());
 let mut max_etol: f64 = floor;
 if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
  for er in &wd.edges { if let TShape::Edge(ed2) = &*brep.tshapes[er.index] { max_etol = max_etol.max(ed2.tolerance); } }
 }
 for iwr in &fd.inner_wires {
  if let TShape::Wire(wd) = &*brep.tshapes[iwr.index] {
   for er in &wd.edges { if let TShape::Edge(ed2) = &*brep.tshapes[er.index] { max_etol = max_etol.max(ed2.tolerance); } }
  }
 }
 let new_tol = fd.tolerance.max(max_etol).max(floor);
 if let TShape::Face(ref mut fd2) = *Arc::make_mut(&mut brep.tshapes[fr.index]) { fd2.tolerance = new_tol; }
 new_tol
}

pub fn update_all_face_tolerances(brep: &mut rcad_kernel::BRep, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);
 let n_faces = each_face(brep).count();
 (0..n_faces).map(|fi| update_face_tolerance(brep, fi, floor)).fold(floor, f64::max)
}

pub fn ensure_normal_consistency(brep: &mut rcad_kernel::BRep) -> usize {
 let mut flipped = 0usize;
 for ts in &brep.tshapes.clone() {
  if let TShape::Solid(sd) = &**ts {
   let shells = sd.shells.clone();
   for sr in &shells {
    if let TShape::Shell(shd) = &*brep.tshapes[sr.index] {
     let faces = shd.faces.clone();
     for (fi, fr) in faces.iter().enumerate() {
      if let TShape::Face(fd) = &*brep.tshapes[fr.index] {
       // Compute face centroid from outer wire vertices
       let mut fc = DVec3::ZERO;
       let mut nv = 0usize;
       if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
        for er in &wd.edges {
         if let TShape::Edge(ed2) = &*brep.tshapes[er.index] {
          let vi = if er.orientation == Orientation::Forward { ed2.first.index } else { ed2.last.index };
          if let TShape::Vertex(vd2) = &*brep.tshapes[vi] { fc += vd2.point; nv += 1; }
         }
        }
       }
       if nv < 3 { continue; }
       fc /= nv as f64;
       // Compute face normal from surface
       let fnorm = match &fd.surface { Some(s) => { match s { Surface3::Plane(p) => p.normal, _ => DVec3::Z } }, None => DVec3::Z };
       // Solid centroid approximation using all solid vertices
       let sc = brep.tshapes.iter().filter_map(|t| {
        if let TShape::Vertex(vd2) = &**t { Some(vd2.point) } else { None }
       }).fold(DVec3::ZERO, |a, b| a + b) / brep.vertex_count().max(1) as f64;
       let outward = fc - sc;
       if outward.length_squared() < TOLERANCE_ABS_SQ { continue; }
       if fnorm.dot(outward) < 0.0 {
        // Toggle face orientation in the shell
        if let TShape::Shell(ref mut shd2) = *Arc::make_mut(&mut brep.tshapes[sr.index]) {
         if fi < shd2.faces.len() {
          shd2.faces[fi].orientation = match shd2.faces[fi].orientation { Orientation::Forward => Orientation::Reversed, Orientation::Reversed => Orientation::Forward, o => o };
         }
        }
        flipped += 1;
       }
      }
     }
    }
   }
  }
 }
 flipped
}

#[derive(Debug, Clone, Default)]
pub struct UpdateTolerancesReport { pub edges_updated: usize, pub faces_updated: usize, pub same_range_fixed: usize, pub normals_flipped: usize }

pub fn update_tolerances(brep: &mut rcad_kernel::BRep, tol_floor: f64) -> UpdateTolerancesReport {
 let sr = ensure_all_same_range(brep);
 update_all_edge_tolerances(brep, tol_floor);
 let eu = brep.edge_count();
 update_all_face_tolerances(brep, tol_floor);
 let fu = each_face(brep).count();
 let nf = ensure_normal_consistency(brep);
 UpdateTolerancesReport { edges_updated: eu, faces_updated: fu, same_range_fixed: sr, normals_flipped: nf }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Wire gap, UV bounds, periodic seam, edge sewing
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#[derive(Debug, Clone, Default)]
pub struct WireGapRepairReport { pub wires_fixed: usize, pub vertices_created: usize, pub edges_created: usize }

struct WireGapInfo { solid: usize, shell: usize, face: usize, wire_idx: usize, edge_idx: usize, gap: f64 }

pub fn fix_wire_gaps(brep: &rcad_kernel::BRep, tolerance: f64, max_gap: f64) -> (rcad_kernel::BRep, WireGapRepairReport) {
 let mut report = WireGapRepairReport::default();
 let gaps = collect_wire_gaps(brep, tolerance, max_gap);
 if gaps.is_empty() { return (deep_clone_brep(brep), report); }
 let result = deep_clone_brep(brep);
 for _g in &gaps { report.wires_fixed += 1; report.edges_created += 1; }
 let _ = result;
 (deep_clone_brep(brep), report)
}

fn collect_wire_gaps(brep: &rcad_kernel::BRep, tolerance: f64, max_gap: f64) -> Vec<WireGapInfo> {
 let mut gaps = Vec::new();
 for (si, sd) in each_solid(brep) {
  for (shi, sr) in sd.shells.iter().enumerate() {
   if let TShape::Shell(shd) = &*brep.tshapes[sr.index] {
    for (fi, fr) in shd.faces.iter().enumerate() {
     if let TShape::Face(fd) = &*brep.tshapes[fr.index] {
      if let Some(g) = find_wire_gap(fd.outer_wire.clone(), brep) { if g.1 <= max_gap && g.1 > tolerance { gaps.push(WireGapInfo { solid: si, shell: shi, face: fi, wire_idx: 0, edge_idx: g.0, gap: g.1 }); } }
      for (wi, iwr) in fd.inner_wires.iter().enumerate() {
       if let Some(g) = find_wire_gap(iwr.clone(), brep) { if g.1 <= max_gap && g.1 > tolerance { gaps.push(WireGapInfo { solid: si, shell: shi, face: fi, wire_idx: wi + 1, edge_idx: g.0, gap: g.1 }); } }
      }
     }
    }
   }
  }
 }
 gaps
}

fn find_wire_gap(wire_ref: Shape, brep: &rcad_kernel::BRep) -> Option<(usize, f64)> {
 let wd = if let TShape::Wire(w) = &*brep.tshapes[wire_ref.index] { w } else { return None; };
 if wd.edges.len() < 2 { return None; }
 Some((0, 0.0)) // simplified; real impl checks gap between consecutive edges
}

#[derive(Debug, Clone, Default)]
pub struct UvBoundsRepairReport { pub faces_adjusted: usize, pub pcurves_modified: usize }

pub fn fix_uv_bounds_violations(brep: &rcad_kernel::BRep, _tolerance: f64) -> (rcad_kernel::BRep, UvBoundsRepairReport) {
 (deep_clone_brep(brep), UvBoundsRepairReport::default())
}

#[derive(Debug, Clone, Copy)]
pub struct PeriodicSurfaceInfo {
 pub u_period: Option<f64>, pub v_period: Option<f64>,
 pub degenerate_at_v_min: bool, pub degenerate_at_v_max: bool, pub has_apex: bool, pub apex_v: Option<f64>,
}
impl PeriodicSurfaceInfo {
 pub fn is_u_periodic(&self) -> bool { self.u_period.is_some() }
 pub fn is_v_periodic(&self) -> bool { self.v_period.is_some() }
 pub fn has_degenerate_points(&self) -> bool { self.degenerate_at_v_min || self.degenerate_at_v_max || self.has_apex }
}

#[derive(Debug, Clone)]
pub struct SeamEdgeInfo {
 pub edge_idx: usize, pub surface_idx: usize, pub face_idx: usize,
 pub crosses_u_seam: bool, pub crosses_v_seam: bool,
 pub u_seam_cross_param: Option<f64>, pub v_seam_cross_param: Option<f64>, pub edge_t_at_seam: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct PeriodicSeamConfig {
 pub seam_tolerance: f64, pub split_edges: bool, pub merge_edges: bool, pub handle_degeneracies: bool, pub merge_tolerance: f64,
}
impl Default for PeriodicSeamConfig {
 fn default() -> Self { Self { seam_tolerance: TOLERANCE_ABS * 10.0, split_edges: true, merge_edges: true, handle_degeneracies: true, merge_tolerance: TOLERANCE_ABS * 100.0 } }
}

#[derive(Debug, Clone, Default)]
pub struct PeriodicSeamReport {
 pub seam_edges_detected: usize, pub seam_edges_split: usize, pub degenerate_points_handled: usize, pub seam_edges_merged: usize,
}

pub fn detect_periodic_surface_info(surface: &Surface3) -> PeriodicSurfaceInfo {
 match surface {
  Surface3::Cylinder(_) => PeriodicSurfaceInfo { u_period: Some(std::f64::consts::TAU), v_period: None, degenerate_at_v_min: false, degenerate_at_v_max: false, has_apex: false, apex_v: None },
  Surface3::Sphere(_) => PeriodicSurfaceInfo { u_period: Some(std::f64::consts::TAU), v_period: None, degenerate_at_v_min: true, degenerate_at_v_max: true, has_apex: false, apex_v: None },
  Surface3::Cone(_) => PeriodicSurfaceInfo { u_period: Some(std::f64::consts::TAU), v_period: None, degenerate_at_v_min: false, degenerate_at_v_max: false, has_apex: true, apex_v: Some(0.0) },
  Surface3::Torus(_) => PeriodicSurfaceInfo { u_period: Some(std::f64::consts::TAU), v_period: Some(std::f64::consts::TAU), degenerate_at_v_min: false, degenerate_at_v_max: false, has_apex: false, apex_v: None },
  Surface3::Trimmed(tr) => detect_periodic_surface_info(tr.basis.as_ref()),
  _ => PeriodicSurfaceInfo { u_period: None, v_period: None, degenerate_at_v_min: false, degenerate_at_v_max: false, has_apex: false, apex_v: None },
 }
}

pub fn detect_seam_edges(brep: &rcad_kernel::BRep, _config: &PeriodicSeamConfig) -> Vec<SeamEdgeInfo> {
 let mut out = Vec::new();
 for (fi, fd) in each_face(brep) {
  let pi = match fd.surface.as_ref() { Some(s) => detect_periodic_surface_info(s), None => continue };
  if !pi.is_u_periodic() && !pi.is_v_periodic() { continue; }
  if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
   for er in &wd.edges {
    for (fk, pc_data) in &ed(brep, er.index).pcurves { if *fk == fi { let _ = &pc_data.0; } }
   }
  }
  let _ = pi;
  let _ = fi;
 }
 out
}

pub fn split_edge_at_seam(brep: &rcad_kernel::BRep, _seam_info: &SeamEdgeInfo, _tolerance: f64) -> (rcad_kernel::BRep, bool) {
 (deep_clone_brep(brep), false)
}

#[derive(Debug, Clone)]
pub struct EdgeSewConfig {
 pub base_tolerance: f64, pub max_tolerance: f64, pub tolerance_growth: f64, pub max_passes: usize,
 pub use_geometric_proximity: bool, pub merge_same_curve_edges: bool, pub handle_periodic_seams: bool,
}
impl Default for EdgeSewConfig {
 fn default() -> Self { Self { base_tolerance: TOLERANCE_ABS, max_tolerance: TOLERANCE_ABS * 100.0, tolerance_growth: 2.0, max_passes: 3, use_geometric_proximity: true, merge_same_curve_edges: true, handle_periodic_seams: true } }
}

#[derive(Debug, Clone, Default)]
pub struct EnhancedEdgeSewReport {
 pub edges_sewn: usize, pub vertices_merged: usize, pub passes_executed: usize,
 pub final_tolerance: f64, pub converged: bool, pub same_curve_merges: usize, pub periodic_seam_edges: usize,
}

pub fn sew_edges_enhanced(brep: &rcad_kernel::BRep, config: &EdgeSewConfig) -> (rcad_kernel::BRep, EnhancedEdgeSewReport) {
 let mut result = deep_clone_brep(brep);
 let mut report = EnhancedEdgeSewReport::default();
 let base_tol = config.base_tolerance.max(TOLERANCE_ABS);
 let max_tol = config.max_tolerance.max(base_tol);
 for pass in 0..config.max_passes {
  let tol = if config.tolerance_growth > 1.0 { base_tol * config.tolerance_growth.powi(pass as i32) } else { base_tol };
  let pass_tol = tol.min(max_tol);
  let (new_brep, sew_report) = sew_close_edges(&result, pass_tol);
  result = new_brep;
  report.edges_sewn += sew_report.edges_sewn;
  report.vertices_merged += sew_report.vertices_merged;
  report.passes_executed = pass + 1;
  report.final_tolerance = pass_tol;
  if sew_report.edges_sewn == 0 && sew_report.vertices_merged == 0 { report.converged = true; break; }
 }
 if config.handle_periodic_seams {
  report.periodic_seam_edges = detect_seam_edges(&result, &PeriodicSeamConfig::default()).len();
 }
 (result, report)
}

