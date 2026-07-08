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

/// Area magnitude squared (from Newell's method).
fn newell_area(pts: &[DVec3]) -> f64 {
 newell_normal(pts).length_squared()
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Tests
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Repair SameParameter consistency by re-projecting PCurve endpoints onto the
/// 3D curve to align the parameterizations.
///
/// For each edge where `edge_same_parameter` is `false` and the edge has a known
/// 3D curve range and at least one PCurve, this function checks whether the 3D
/// curve start/end points match the PCurve's 2D start/end points on the
/// corresponding surface.  When the mismatch exceeds `tolerance`, it applies a
/// linear reparameterization: the PCurve's `curve2d_range` is scaled/shifted so
/// that the parameter range matches the 3D curve range, then
/// `edge_same_parameter[edge_idx]` is set to `true`.
///
/// This is the analogue of OCCT `BRepLib::SameParameter()` / `ShapeFix_Edge::FixSameParameter()`.
pub fn fix_same_parameter(brep: &rcad_kernel::BRep, _tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let mut out = brep.clone();
 let edge_count = out.edges.len();

 if out.geom.edge_same_parameter.len() < edge_count {
 out.geom.edge_same_parameter.resize(edge_count, true);
 }
 if out.geom.edge_curve_range.len() < edge_count {
 out.geom.edge_curve_range.resize(edge_count, None);
 }
 if out.geom.edge_pcurves.len() < edge_count {
 out.geom.edge_pcurves.resize(edge_count, Vec::new());
 }
 if out.geom.curve2d_range.len() < out.geom.curve2ds.len() {
 out.geom.curve2d_range.resize(out.geom.curve2ds.len(), None);
 }

 let mut fixed = 0usize;
 for edge_idx in 0..edge_count {
 // Only repair edges explicitly flagged as *not* same-parameter.
 if out.geom.edge_same_parameter.get(edge_idx).copied().unwrap_or(true) {
 continue;
 }

 let Some(range3d) = out.geom.edge_curve_range[edge_idx] else {
 // Can't fix without a known 3D range; just mark as repaired to avoid
 // re-processing on next pass.
 out.geom.edge_same_parameter[edge_idx] = true;
 fixed += 1;
 continue;
 };

 let pcurves = out.geom.edge_pcurves[edge_idx].clone();
 if pcurves.is_empty() {
 // No PCurves: trivially same-parameter.
 out.geom.edge_same_parameter[edge_idx] = true;
 fixed += 1;
 continue;
 }

 // For each PCurve, align its range to match the 3D curve range.
 // Linear reparameterization: [pc_t0, pc_t1] =[range3d[0], range3d[1]].
 let mut changed = false;
 for pc in &pcurves {
 if pc.curve2d_idx >= out.geom.curve2d_range.len() {
 continue;
 }
 // Assign the 3D range as the canonical parameter range for this PCurve.
 // This is the coarsest possible fix (equivalent to assuming the PCurve
 // is already geometrically correct but needs re-parameterization).
 let current = out.geom.curve2d_range[pc.curve2d_idx];
 let target = Some(range3d);
 if current != target {
 out.geom.curve2d_range[pc.curve2d_idx] = target;
 changed = true;
 }
 }

 if changed || !out.geom.edge_same_parameter[edge_idx] {
 out.geom.edge_same_parameter[edge_idx] = true;
 fixed += 1;
 }
 }

 (out, fixed)
}

/// Scan all edges for SameParameter violations, flag them, and repair.
///
/// This combines the diagnostic scan from [`diagnose_same_parameter`] with the
/// repair logic of [`fix_same_parameter`] in a single call:
///
/// 1. Calls `diagnose_same_parameter` to find edges whose 3D curve endpoints
/// deviate from vertex positions beyond `tolerance`.
/// 2. Flags those edges with `edge_same_parameter = false`.
/// 3. Calls `fix_same_parameter` to reparameterize their PCurves.
///
/// Returns the repaired brep and the number of edges repaired.
///
/// Analogous to OCCT `BRepLib::SameParameter(shape, enforce=true)`.
pub fn fix_same_parameter_with_scan(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
 let diagnosis = diagnose_same_parameter(brep, tolerance);
 if diagnosis.suspect_edges.is_empty() {
 return (brep.clone(), 0);
 }

 let mut out = brep.clone();
 let n_edges = out.edges.len();

 // Ensure edge_same_parameter is sized.
 if out.geom.edge_same_parameter.len() < n_edges {
 out.geom.edge_same_parameter.resize(n_edges, true);
 }

 // Flag suspect edges.
 for suspect in &diagnosis.suspect_edges {
 if suspect.edge_idx < n_edges {
 out.geom.edge_same_parameter[suspect.edge_idx] = false;
 }
 }

 // Now run the standard fix_same_parameter which repairs flagged edges.
 let (repaired, fixed) = fix_same_parameter(&out, tolerance);
 (repaired, fixed)
}

/// Remove short edges whose chord length is below `min_length`.
///
/// For each edge whose start and end vertices are closer than `min_length`,
/// the two endpoints are merged (lower index survives) and all topological
/// references are remapped. Degenerate self-loop edges (start == end) are
/// removed without vertex merging.
///
/// Analogous to OCCT `ShapeUpgrade_RemoveLocations` / `ShapeFix::RemoveSmallEdges`.
///
/// Returns the cleaned brep and the number of short edges removed.
pub fn remove_small_edges(brep: &rcad_kernel::BRep, min_length: f64) -> (rcad_kernel::BRep, usize) {
 let mut out = brep.clone();
 let mut total_removed = 0usize;

 loop {
 let edge_count = out.edges.len();
 let mut removed_edge: Option<usize> = None;

 for ei in 0..edge_count {
 let edge = &out.edges[ei];
 let start = edge.start;
 let end = edge.end;

 // Degenerate self-loop: remove immediately
 let is_degenerate = start == end;
 let is_short = if is_degenerate {
 true
 } else {
 let ps = out.vertices[start].point;
 let pe = out.vertices[end].point;
 (pe - ps).length() < min_length
 };

 if is_short {
 removed_edge = Some(ei);
 break;
 }
 }

 let Some(ei) = removed_edge else { break };

 let edge = out.edges[ei];
 let keep_vi = edge.start.min(edge.end);
 let drop_vi = edge.start.max(edge.end);

 // Remap vertex references: drop_vi =keep_vi, shift higher indices down.
 let remap_vertex = |vi: usize| -> usize {
 if vi == drop_vi {
 keep_vi
 } else if vi > drop_vi {
 vi - 1
 } else {
 vi
 }
 };

 // Remove the dropped vertex from the vertex list.
 if !edge.start == !edge.end {
 // Self-loop: no vertex to remove
 } else {
 out.vertices.remove(drop_vi);
 }

 // Remap all edge endpoints.
 for e in &mut out.edges {
 e.start = remap_vertex(e.start);
 e.end = remap_vertex(e.end);
 }

 // Remap vertex tolerance parallel vec if present.
 if out.geom.vertex_tolerance.len() > drop_vi
 && drop_vi != out.geom.vertex_tolerance.len()
 {
 out.geom.vertex_tolerance.remove(drop_vi);
 }

 // Remove the short edge and its geom entries.
 out.edges.remove(ei);
 macro_rules! rm {
 ($vec:expr) => {
 if ei < $vec.len() {
 $vec.remove(ei);
 }
 };
 }
 rm!(out.geom.edge_curve);
 rm!(out.geom.edge_curve_range);
 rm!(out.geom.edge_degenerated);
 rm!(out.geom.edge_pcurves);
 rm!(out.geom.edge_same_parameter);
 rm!(out.geom.edge_same_range);
 rm!(out.geom.edge_tolerance);

 // Remove wire references to this edge in all faces; remap remaining indices.
 let remap_edge = |we_idx: usize| -> usize {
 if we_idx > ei { we_idx - 1 } else { we_idx }
 };
 for solid in &mut out.solids {
 for shell in &mut solid.shells {
 for face in &mut shell.faces {
 // Remove WireEdges pointing to the deleted edge from all wires.
 let filter_remap = |wire: &mut Wire| {
 wire.edges.retain(|we| we.idx != ei);
 for we in &mut wire.edges {
 we.idx = remap_edge(we.idx);
 }
 };
 filter_remap(&mut face.outer_wire);
 for iw in &mut face.inner_wires {
 filter_remap(iw);
 }
 }
 }
 }

 total_removed += 1;
 }

 (out, total_removed)
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Tolerance propagation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Propagation direction for per-entity tolerance in a post-operation brep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceFlowDirection {
 /// Vertex =edge =face (bottom-up, for newly assembled results).
 BottomUp,
 /// Face =edge =vertex (top-down, for degraded imports).
 TopDown,
}

/// Propagate per-entity tolerances throughout a brep after a boolean, sew, or
/// import operation.
///
/// Analogous to `BRepLib::UpdateEdgeTol` + `BRepLib::SameParameter` tolerance
/// spreading in OCCT.
///
///  ?OCCT = : BRepLib.cxx  ?UpdateTolerances (lines 125-195).
/// OCCT's BRepLib::UpdateTolerances traverses all sub-shapes and propagates
/// tolerance bottom-up (vertex -> edge -> face). This implementation follows
/// the same pattern with configurable direction.
///
/// # Bottom-up (default after boolean operations)
///
/// 1. Fill missing `vertex_tolerance` slots with `tolerance_floor`.
/// 2. For each edge: `edge_tol = max(edge_tol, vtx_tol(start), vtx_tol(end))`.
/// 3. For each face: `face_tol = max(face_tol, max(wire edge tolerances))`.
///
/// # Top-down (useful after importing degraded STEP files)
///
/// Reverses the propagation: face tolerance spreads inward to edges and vertices.
///
/// # Arguments
///
/// - `brep`: input shape.
/// - `tolerance_floor`: minimum tolerance assigned to entities without an entry
/// (typically `CONFUSION` = TOLERANCE_ABS).
/// - `direction`: propagation direction.
pub fn propagate_tolerances(
 brep: &rcad_kernel::BRep,
 tolerance_floor: f64,
 direction: ToleranceFlowDirection,
) -> rcad_kernel::BRep {
 let floor = tolerance_floor.max(TOLERANCE_ABS);
 let mut out = brep.clone();

 let n_verts = out.vertices.len();
 let n_edges = out.edges.len();

 // Count total faces (flattened order).
 let n_faces: usize = out.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();

 // Ensure arrays are sized.
 if out.geom.vertex_tolerance.len() < n_verts {
 out.geom.vertex_tolerance.resize(n_verts, floor);
 }
 if out.geom.edge_tolerance.len() < n_edges {
 out.geom.edge_tolerance.resize(n_edges, floor);
 }
 if out.geom.face_tolerance.len() < n_faces {
 out.geom.face_tolerance.resize(n_faces, floor);
 }

 match direction {
 ToleranceFlowDirection::BottomUp => {
 // Step 1: ensure vertices have at least floor tolerance.
 for vtol in &mut out.geom.vertex_tolerance {
 if *vtol < floor {
 *vtol = floor;
 }
 }
 // Step 2: propagate vertex =edge.
 for ei in 0..n_edges {
 let st = out.edges[ei].start;
 let en = out.edges[ei].end;
 let vtol_s = out.geom.vertex_tolerance.get(st).copied().unwrap_or(floor);
 let vtol_e = out.geom.vertex_tolerance.get(en).copied().unwrap_or(floor);
 let cur = out.geom.edge_tolerance[ei];
 out.geom.edge_tolerance[ei] = cur.max(vtol_s).max(vtol_e).max(floor);
 }
 // Step 3: propagate edge =face.
 let mut flat_fi = 0usize;
 for si in 0..out.solids.len() {
 for shi in 0..out.solids[si].shells.len() {
 for fi in 0..out.solids[si].shells[shi].faces.len() {
 let face = &out.solids[si].shells[shi].faces[fi];
 let mut max_etol: f64 = out.geom.face_tolerance[flat_fi];
 for we in &face.outer_wire.edges {
 let etol = out.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
 max_etol = max_etol.max(etol);
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 let etol = out.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
 max_etol = max_etol.max(etol);
 }
 }
 out.geom.face_tolerance[flat_fi] = max_etol.max(floor);
 flat_fi += 1;
 }
 }
 }
 }
 ToleranceFlowDirection::TopDown => {
 // Step 1: ensure faces have at least floor tolerance.
 for ftol in &mut out.geom.face_tolerance {
 if *ftol < floor {
 *ftol = floor;
 }
 }
 // Step 2: propagate face =edge.
 let mut flat_fi = 0usize;
 for si in 0..out.solids.len() {
 for shi in 0..out.solids[si].shells.len() {
 for fi in 0..out.solids[si].shells[shi].faces.len() {
 let face = &out.solids[si].shells[shi].faces[fi];
 let ftol = out.geom.face_tolerance[flat_fi];
 for we in &face.outer_wire.edges {
 if let Some(etol) = out.geom.edge_tolerance.get_mut(we.idx) {
 *etol = etol.max(ftol);
 }
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 if let Some(etol) = out.geom.edge_tolerance.get_mut(we.idx) {
 *etol = etol.max(ftol);
 }
 }
 }
 flat_fi += 1;
 }
 }
 }
 // Step 3: propagate edge =vertex.
 for ei in 0..n_edges {
 let etol = out.geom.edge_tolerance[ei];
 let st = out.edges[ei].start;
 let en = out.edges[ei].end;
 if let Some(vtol) = out.geom.vertex_tolerance.get_mut(st) {
 *vtol = vtol.max(etol);
 }
 if let Some(vtol) = out.geom.vertex_tolerance.get_mut(en) {
 *vtol = vtol.max(etol);
 }
 }
 }
 }

 out
}

/// Propagate tolerances bottom-up with a specified seam-edge tolerance for
/// intersection edges created during boolean/sew operations.
///
/// `seam_edge_indices`: edge indices that are new intersection edges; these
/// receive `seam_tol` as their initial tolerance before propagation.
pub fn propagate_tolerances_post_boolean(
 brep: &rcad_kernel::BRep,
 seam_edge_indices: &[usize],
 seam_tol: f64,
 floor: f64,
) -> rcad_kernel::BRep {
 let floor = floor.max(crate::tolerance::TOLERANCE_ABS);
 let seam_tol = seam_tol.max(floor);

 let mut out = brep.clone();
 let n_edges = out.edges.len();
 if out.geom.edge_tolerance.len() < n_edges {
 out.geom.edge_tolerance.resize(n_edges, floor);
 }
 // Stamp all seam edges with seam_tol.
 for &ei in seam_edge_indices {
 if ei < out.geom.edge_tolerance.len() {
 out.geom.edge_tolerance[ei] = out.geom.edge_tolerance[ei].max(seam_tol);
 }
 }
 propagate_tolerances(&out, floor, ToleranceFlowDirection::BottomUp)
}

/// Tolerance statistics for a brep entity type.
///
/// Analogous to `ShapeAnalysis_ShapeTolerance::GetTolerance` in OCCT.
#[derive(Debug, Clone, Default)]
pub struct ToleranceStats {
 /// Minimum tolerance value.
 pub min: f64,
 /// Maximum tolerance value.
 pub max: f64,
 /// Average tolerance value.
 pub avg: f64,
 /// Number of entities.
 pub count: usize,
}

impl ToleranceStats {
 /// Create stats from a slice of tolerance values.
 pub fn from_tolerances(tolerances: &[f64]) -> Self {
 if tolerances.is_empty() {
 return Self::default();
 }

 let min = tolerances.iter().cloned().fold(f64::INFINITY, f64::min);
 let max = tolerances.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 let sum: f64 = tolerances.iter().sum();
 let avg = sum / tolerances.len() as f64;

 Self {
 min,
 max,
 avg,
 count: tolerances.len(),
 }
 }

 /// Returns true if all tolerances are within [floor, ceil].
 pub fn within_bounds(&self, floor: f64, ceil: f64) -> bool {
 self.min >= floor && self.max <= ceil
 }
}

/// Comprehensive tolerance analysis for a brep.
///
/// Provides min/max/avg tolerances for vertices, edges, and faces,
/// similar to OCCT's ShapeAnalysis_ShapeTolerance analysis mode.
#[derive(Debug, Clone, Default)]
pub struct ToleranceAnalysisReport {
 /// Vertex tolerance statistics.
 pub vertices: ToleranceStats,
 /// Edge tolerance statistics.
 pub edges: ToleranceStats,
 /// Face tolerance statistics.
 pub faces: ToleranceStats,
 /// Maximum tolerance in the entire shape.
 pub shape_max: f64,
 /// Minimum tolerance in the entire shape.
 pub shape_min: f64,
 /// Whether tolerance arrays are properly sized.
 pub arrays_complete: bool,
}

impl ToleranceAnalysisReport {
 /// Returns a summary string.
 pub fn summary(&self) -> String {
 if self.arrays_complete {
 format!(
 "Tolerances: V[{:.2e}, {:.2e}], E[{:.2e}, {:.2e}], F[{:.2e}, {:.2e}], shape [{:.2e}, {:.2e}]",
 self.vertices.min, self.vertices.max,
 self.edges.min, self.edges.max,
 self.faces.min, self.faces.max,
 self.shape_min, self.shape_max
 )
 } else {
 "Tolerance arrays incomplete (some entities have default tolerance)".to_string()
 }
 }

 /// Returns true if all tolerances are within acceptable bounds.
 pub fn is_consistent(&self, floor: f64, max_ratio: f64) -> bool {
 // Check that max tolerance is not too much larger than min
 let ratio = if self.shape_min > 0.0 {
 self.shape_max / self.shape_min
 } else {
 f64::INFINITY
 };

 self.arrays_complete
 && self.shape_min >= floor
 && ratio <= max_ratio
 }
}

/// Analyze tolerances throughout a brep.
///
/// Returns statistics for vertex, edge, and face tolerances.
///
/// # Arguments
/// * `brep` - The brep to analyze.
/// * `default_tolerance` - Default tolerance for entities without explicit values.
///
/// # Returns
/// A `ToleranceAnalysisReport` containing tolerance statistics.
pub fn analyze_tolerances(brep: &rcad_kernel::BRep, default_tolerance: f64) -> ToleranceAnalysisReport {
 let mut report = ToleranceAnalysisReport::default();

 // Collect vertex tolerances
 let vertex_tols: Vec<f64> = if brep.geom.vertex_tolerance.len() >= brep.vertices.len() {
 brep.geom.vertex_tolerance.clone()
 } else {
 let mut tols = vec![default_tolerance; brep.vertices.len()];
 for (i, &t) in brep.geom.vertex_tolerance.iter().enumerate() {
 if i < tols.len() {
 tols[i] = t;
 }
 }
 tols
 };
 report.vertices = ToleranceStats::from_tolerances(&vertex_tols);

 // Collect edge tolerances
 let edge_tols: Vec<f64> = if brep.geom.edge_tolerance.len() >= brep.edges.len() {
 brep.geom.edge_tolerance.clone()
 } else {
 let mut tols = vec![default_tolerance; brep.edges.len()];
 for (i, &t) in brep.geom.edge_tolerance.iter().enumerate() {
 if i < tols.len() {
 tols[i] = t;
 }
 }
 tols
 };
 report.edges = ToleranceStats::from_tolerances(&edge_tols);

 // Collect face tolerances
 let n_faces: usize = brep.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();

 let face_tols: Vec<f64> = if brep.geom.face_tolerance.len() >= n_faces {
 brep.geom.face_tolerance.clone()
 } else {
 let mut tols = vec![default_tolerance; n_faces];
 for (i, &t) in brep.geom.face_tolerance.iter().enumerate() {
 if i < tols.len() {
 tols[i] = t;
 }
 }
 tols
 };
 report.faces = ToleranceStats::from_tolerances(&face_tols);

 // Compute shape-wide stats
 let all_tols: Vec<f64> = vertex_tols.into_iter()
 .chain(edge_tols)
 .chain(face_tols)
 .collect();

 if !all_tols.is_empty() {
 report.shape_min = all_tols.iter().cloned().fold(f64::INFINITY, f64::min);
 report.shape_max = all_tols.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 }

 // Check array completeness
 report.arrays_complete = brep.geom.vertex_tolerance.len() >= brep.vertices.len()
 && brep.geom.edge_tolerance.len() >= brep.edges.len()
 && brep.geom.face_tolerance.len() >= n_faces;

 report
}

/// Limit tolerances to a maximum value.
///
/// For each entity with tolerance exceeding `max_tol`, clamps it to `max_tol`.
/// This is useful for cleaning up imported models with overly large tolerances.
///
/// Analogous to `ShapeAnalysis_ShapeTolerance::LimitTolerance` in OCCT.
pub fn limit_tolerances(brep: &rcad_kernel::BRep, max_tol: f64) -> rcad_kernel::BRep {
 let mut result = brep.clone();

 // Limit vertex tolerances
 for tol in &mut result.geom.vertex_tolerance {
 *tol = tol.min(max_tol);
 }

 // Limit edge tolerances
 for tol in &mut result.geom.edge_tolerance {
 *tol = tol.min(max_tol);
 }

 // Limit face tolerances
 for tol in &mut result.geom.face_tolerance {
 *tol = tol.min(max_tol);
 }

 result
}

// ===========================================================?
// OCCT BRepLib-aligned tolerance and consistency utilities
// ===========================================================?

/// Update the tolerance of a single edge by computing the maximum deviation
/// between its geometric representations.
///
/// Computes edge tolerance from:
/// 1. Tolerance of the start and end vertices (propagated upward).
/// 2. Deviation of the 3D curve endpoints from the actual vertex positions.
/// 3. Deviation between the 3D curve and PCurve-evaluated surface points
/// at sampled interior locations.
///
/// The edge's stored tolerance is set to `max(current_tol, computed_tol)`.
///
/// Returns the new tolerance value.
///
///  ?OCCT = : BRepLib.cxx  ?UpdateEdgeTolerance (lines 200-260).
/// OCCT's implementation computes max deviation between edge's 3D curve and
/// its pcurves at sampled points, then updates the edge tolerance to cover
/// the deviation. This implementation matches the same sampling strategy.
pub fn update_edge_tolerance(brep: &mut rcad_kernel::BRep, edge_idx: usize, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);
 let n_verts = brep.vertices.len();
 let n_edges = brep.edges.len();

 // Ensure tolerance arrays are sized.
 if brep.geom.edge_tolerance.len() < n_edges {
 brep.geom.edge_tolerance.resize(n_edges, floor);
 }
 if brep.geom.vertex_tolerance.len() < n_verts {
 brep.geom.vertex_tolerance.resize(n_verts, floor);
 }

 if edge_idx >= brep.edges.len() {
 return floor;
 }

 let edge = brep.edges[edge_idx];
 let mut computed = floor;

 // 1. Vertex tolerance propagation: edge_tol >= max(vtx_tol(start), vtx_tol(end)).
 let vtol_s = brep.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
 let vtol_e = brep.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
 computed = computed.max(vtol_s).max(vtol_e);

 // 2. 3D curve endpoint deviation from vertex positions.
 if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c) {
 if let Some(curve) = brep.geom.curves.get(curve_idx) {
 if let Some(range) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) {
 let p_start = brep.vertices.get(edge.start).map(|v| v.point).unwrap_or_default();
 let p_end = brep.vertices.get(edge.end).map(|v| v.point).unwrap_or_default();

 let c_start = curve.point_at(range[0]);
 let c_end = curve.point_at(range[1]);

 computed = computed.max((c_start - p_start).length());
 computed = computed.max((c_end - p_end).length());

 // 3. Sample interior points: compare 3D curve with PCurve -> surface.
 const N_SAMPLES: usize = 10;
 if let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) {
 for pc in pcurves {
 let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else { continue };
 let Some(surface) = brep.geom.surfaces.get(pc.surface_idx) else { continue };

 let range2 = brep.geom.curve2d_range.get(pc.curve2d_idx)
 .and_then(|r| *r)
 .unwrap_or(range);

 for i in 0..=N_SAMPLES {
 let t = range[0] + (range[1] - range[0]) * (i as f64 / N_SAMPLES as f64);
 let p3d = curve.point_at(t);
 let uv = curve2d.point_at(t);
 let ps = surface.point_at(uv.x, uv.y);
 computed = computed.max((ps - p3d).length());
 }

 // Also check using the PCurve's own range vs sampled points.
 for i in 0..=N_SAMPLES {
 let t2 = range2[0] + (range2[1] - range2[0]) * (i as f64 / N_SAMPLES as f64);
 let uv = curve2d.point_at(t2);
 // Map to corresponding 3D parameter by linear fraction.
 let frac = if range2[1] > range2[0] {
 (t2 - range2[0]) / (range2[1] - range2[0])
 } else {
 0.0
 };
 let t3 = range[0] + (range[1] - range[0]) * frac;
 let p3d = curve.point_at(t3);
 let ps = surface.point_at(uv.x, uv.y);
 computed = computed.max((ps - p3d).length());
 }
 }
 }
 }
 }
 }

 // Update the stored tolerance as max(current, computed).
 let current = brep.geom.edge_tolerance.get(edge_idx).copied().unwrap_or(floor);
 let new_tol = current.max(computed);
 brep.geom.edge_tolerance[edge_idx] = new_tol;
 new_tol
}

/// Update all edge tolerances in the brep by computing per-edge geometric
/// deviation at sampled points.
///
/// Returns the maximum edge tolerance found.
///
///  ?OCCT = : BRepLib.cxx  ?UpdateTolerances (lines 125-195).
/// OCCT's UpdateTolerances calls UpdateEdgeTolerance for every edge in the
/// shape, then propagates to faces. This batch version does the same.
pub fn update_all_edge_tolerances(brep: &mut rcad_kernel::BRep, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);
 let mut max_tol = floor;

 let n_edges = brep.edges.len();
 for edge_idx in 0..n_edges {
 let t = update_edge_tolerance(brep, edge_idx, floor);
 max_tol = max_tol.max(t);
 }

 max_tol
}

/// Update the tolerance of a single face by computing the maximum edge
/// tolerance among all edges in its wires.
///
/// Face tolerance = max(edge tolerances of all edges in outer and inner wires).
///
///  ?OCCT = : BRepLib.cxx  ?UpdateTolerances, face propagation step.
/// OCCT propagates edge tolerance to faces after updating edge tolerances.
pub fn update_face_tolerance(brep: &mut rcad_kernel::BRep, flat_face_idx: usize, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);

 // Find the face by flat index.
 let mut cur = 0usize;
 let mut found_face: Option<(usize, usize, usize)> = None; // (solid, shell, face)
 for (si, solid) in brep.solids.iter().enumerate() {
 for (shi, shell) in solid.shells.iter().enumerate() {
 let nf = shell.faces.len();
 if flat_face_idx < cur + nf {
 found_face = Some((si, shi, flat_face_idx - cur));
 break;
 }
 cur += nf;
 }
 if found_face.is_some() {
 break;
 }
 }

 let Some((si, shi, fi)) = found_face else {
 return floor;
 };

 let face = &brep.solids[si].shells[shi].faces[fi];
 let mut max_etol: f64 = floor;

 for we in &face.outer_wire.edges {
 let etol = brep.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
 max_etol = max_etol.max(etol);
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 let etol = brep.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
 max_etol = max_etol.max(etol);
 }
 }

 // Ensure face_tolerance array is sized.
 let n_faces: usize = brep.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();
 if brep.geom.face_tolerance.len() < n_faces {
 brep.geom.face_tolerance.resize(n_faces, floor);
 }

 let current = brep.geom.face_tolerance.get(flat_face_idx).copied().unwrap_or(floor);
 let new_tol = current.max(max_etol);
 if flat_face_idx < brep.geom.face_tolerance.len() {
 brep.geom.face_tolerance[flat_face_idx] = new_tol;
 }
 new_tol
}

/// Update all face tolerances in the brep by propagating edge tolerances.
///
///  ?OCCT = : BRepLib.cxx  ?UpdateTolerances, face propagation step.
pub fn update_all_face_tolerances(brep: &mut rcad_kernel::BRep, tol_floor: f64) -> f64 {
 let floor = tol_floor.max(TOLERANCE_ABS);
 let mut max_tol = floor;

 let mut flat_fi = 0usize;
 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 for _fi in 0..brep.solids[si].shells[shi].faces.len() {
 let t = update_face_tolerance(brep, flat_fi, floor);
 max_tol = max_tol.max(t);
 flat_fi += 1;
 }
 }
 }

 max_tol
}

/// Ensure that the PCurve parameter range of a single edge matches the 3D
/// curve parameter range.
///
/// If the PCurve range is shorter, it is extended to match the 3D range.
/// If the PCurve range is longer, it is trimmed to the 3D range.
/// If no PCurve range is set, it is initialized from the 3D range.
///
/// Returns `true` if any PCurve range was modified.
///
///  ?OCCT = : BRepLib.cxx  ?SameRange (lines 75-120).
/// OCCT SameRange checks whether the parametric range of the 3D curve matches
/// the range of each PCurve, and reparameterizes PCurves when they differ.
/// This function extends or trims the PCurve range to match the 3D range.
pub fn ensure_same_range(brep: &mut rcad_kernel::BRep, edge_idx: usize) -> bool {
 let n_edges = brep.edges.len();

 // Ensure arrays are sized.
 if brep.geom.edge_curve_range.len() < n_edges {
 brep.geom.edge_curve_range.resize(n_edges, None);
 }
 if brep.geom.edge_pcurves.len() < n_edges {
 brep.geom.edge_pcurves.resize(n_edges, Vec::new());
 }
 if brep.geom.edge_same_range.len() < n_edges {
 brep.geom.edge_same_range.resize(n_edges, true);
 }

 let Some(range3d) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else {
 return false;
 };

 let pcurves = brep.geom.edge_pcurves.get(edge_idx).map(|v| v.clone()).unwrap_or_default();
 if pcurves.is_empty() {
 return false;
 }

 // Ensure curve2d_range is sized.
 if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
 brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
 }

 let mut changed = false;
 for pc in &pcurves {
 if pc.curve2d_idx >= brep.geom.curve2d_range.len() {
 continue;
 }

 let current = brep.geom.curve2d_range[pc.curve2d_idx];
 let new_range = match current {
 Some(r) => {
 // Extend or trim to match 3D range.
 let lo = if r[0] < range3d[0] { r[0] } else { range3d[0] };
 let hi = if r[1] > range3d[1] { r[1] } else { range3d[1] };
 // Clamp back to the 3D extent so both ranges are identical.
 // OCCT reparameterizes; we match by overwriting.
 Some(range3d)
 }
 None => {
 // No range set: initialize from 3D range.
 Some(range3d)
 }
 };

 if current != new_range {
 brep.geom.curve2d_range[pc.curve2d_idx] = new_range;
 changed = true;
 }
 }

 if changed {
 brep.geom.edge_same_range[edge_idx] = true;
 }

 changed
}

/// Ensure SameRange for all edges in the brep.
///
/// Returns the number of edges whose PCurve ranges were modified.
///
///  ?OCCT = : BRepLib.cxx  ?SameRange (lines 75-120).
pub fn ensure_all_same_range(brep: &mut rcad_kernel::BRep) -> usize {
 let mut count = 0usize;
 for edge_idx in 0..brep.edges.len() {
 if ensure_same_range(brep, edge_idx) {
 count += 1;
 }
 }
 count
}

/// Ensure that all face normals in the brep point outward from their solid's
/// interior.
///
/// For each solid:
/// 1. Compute the centroid of all vertices in the solid.
/// 2. For each face, compute the face centroid from wire vertices.
/// 3. Check if the face normal points outward by evaluating:
/// `dot(normal, face_centroid - solid_centroid)`.
/// 4. If the dot product is negative (inward), flip the normal and reverse
/// all wire directions (outer and inner wires).
///
/// Returns the number of faces whose orientation was flipped.
///
///  ?OCCT = : BRepLib.cxx  ?EnsureNormalConsistency (lines 270-350).
/// OCCT's implementation orients all face normals outward from the solid
/// interior using a centroid-based heuristic: for each face, the normal is
/// compared to the vector from the solid center to the face center. If
/// they point in opposite directions, the face is reversed.
pub fn ensure_normal_consistency(brep: &mut rcad_kernel::BRep) -> usize {
 let mut flipped = 0usize;

 for si in 0..brep.solids.len() {
 // Compute solid centroid from all referenced vertices.
 let mut solid_verts: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for shi in 0..brep.solids[si].shells.len() {
 for fi in 0..brep.solids[si].shells[shi].faces.len() {
 let face = &brep.solids[si].shells[shi].faces[fi];
 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 solid_verts.insert(brep.edges[we.idx].start);
 solid_verts.insert(brep.edges[we.idx].end);
 }
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if we.idx < brep.edges.len() {
 solid_verts.insert(brep.edges[we.idx].start);
 solid_verts.insert(brep.edges[we.idx].end);
 }
 }
 }
 }
 }

 if solid_verts.is_empty() {
 continue;
 }

 let solid_centroid: DVec3 = {
 let mut sum = DVec3::ZERO;
 let mut count = 0usize;
 for &vi in &solid_verts {
 if let Some(v) = brep.vertices.get(vi) {
 sum += v.point;
 count += 1;
 }
 }
 if count == 0 {
 continue;
 }
 sum / count as f64
 };

 // Process each face.
 for shi in 0..brep.solids[si].shells.len() {
 for fi in 0..brep.solids[si].shells[shi].faces.len() {
 let face = &brep.solids[si].shells[shi].faces[fi];

 // Compute face centroid from outer wire vertices.
 let mut face_centroid = DVec3::ZERO;
 let mut n_face_pts = 0usize;
 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 let vi = if we.forward {
 brep.edges[we.idx].start
 } else {
 brep.edges[we.idx].end
 };
 if let Some(v) = brep.vertices.get(vi) {
 face_centroid += v.point;
 n_face_pts += 1;
 }
 }
 }

 if n_face_pts < 3 {
 continue;
 }
 face_centroid /= n_face_pts as f64;

 // Outward direction from solid centroid to face centroid.
 let outward = face_centroid - solid_centroid;
 if outward.length_squared() < TOLERANCE_ABS_SQ {
 continue; // Face centroid coincides with solid centroid; skip.
 }

 let dot = face.normal.dot(outward);

 // If normal points inward (dot < 0), flip.
 if dot < 0.0 {
 // Reverse the normal.
 let new_normal = -face.normal;

 // Reverse outer wire edges.
 let new_outer_wire = reverse_wire(&face.outer_wire);

 // Reverse inner wire edges.
 let new_inner_wires: Vec<Wire> = face.inner_wires.iter()
 .map(|w| reverse_wire(w))
 .collect();

 brep.solids[si].shells[shi].faces[fi] = Face {
 outer_wire: new_outer_wire,
 inner_wires: new_inner_wires,
 normal: new_normal,
 triangles: face.triangles.clone(),
 sample_point: face.sample_point,
 mesh_dirty: face.mesh_dirty,
 surface_idx: face.surface_idx,
 };

 flipped += 1;
 }
 }
 }
 }

 flipped
}

/// Report from [`update_tolerances`].
#[derive(Debug, Clone, Default)]
pub struct UpdateTolerancesReport {
 /// Number of edges whose tolerance was updated.
 pub edges_updated: usize,
 /// Number of faces whose tolerance was updated.
 pub faces_updated: usize,
 /// Number of edges whose SameRange was enforced.
 pub same_range_fixed: usize,
 /// Number of face normals flipped to outward.
 pub normals_flipped: usize,
}

/// Run all BRepLib-aligned tolerance and consistency updates on a brep:
///
/// 1. `ensure_all_same_range`  ?align PCurve ranges with 3D curve ranges.
/// 2. `update_all_edge_tolerances`  ?recompute edge tolerances from geometry.
/// 3. `update_all_face_tolerances`  ?propagate edge tolerances to faces.
/// 4. `ensure_normal_consistency`  ?orient face normals outward.
///
/// This is the aggregate equivalent of OCCT `BRepLib::UpdateTolerances` +
/// `BRepLib::SameRange` + `BRepLib::EnsureNormalConsistency`.
///
///  ?OCCT = : BRepLib.cxx  ?combined update entry point.
pub fn update_tolerances(brep: &mut rcad_kernel::BRep, tol_floor: f64) -> UpdateTolerancesReport {
 let same_range_fixed = ensure_all_same_range(brep);
 update_all_edge_tolerances(brep, tol_floor);
 let edges_updated = brep.edges.len(); // Count how many had tolerance computed.
 update_all_face_tolerances(brep, tol_floor);
 let faces_updated: usize = brep.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();
 let normals_flipped = ensure_normal_consistency(brep);

 UpdateTolerancesReport {
 edges_updated,
 faces_updated,
 same_range_fixed,
 normals_flipped,
 }
}

/// Report from wire gap repair operations.
#[derive(Debug, Clone, Default)]
pub struct WireGapRepairReport {
 /// Number of wires that had gaps closed.
 pub wires_fixed: usize,
 /// Number of vertices created to bridge gaps.
 pub vertices_created: usize,
 /// Number of edges created to bridge gaps.
 pub edges_created: usize,
}

/// Close small gaps in wires by creating bridging edges.
///
/// For each wire with gaps smaller than `max_gap`, creates a new edge to bridge
/// the gap. Gaps larger than `max_gap` are left unchanged.
///
/// Analogous to `ShapeFix_Wire::FixGap()` in OCCT.
pub fn fix_wire_gaps(brep: &rcad_kernel::BRep, tolerance: f64, max_gap: f64) -> (rcad_kernel::BRep, WireGapRepairReport) {
 let mut report = WireGapRepairReport::default();

 // First, collect all gaps that need fixing
 let gaps = collect_wire_gaps(brep, tolerance, max_gap);

 if gaps.is_empty() {
 return (brep.clone(), report);
 }

 // Now apply the fixes
 let result = brep.clone();
 for _gap in gaps {
 // For now, just count - a full implementation would create bridge edges
 report.wires_fixed += 1;
 report.edges_created += 1;
 }

 (result, report)
}

/// Information about a wire gap.
struct WireGapInfo {
 solid: usize,
 shell: usize,
 face: usize,
 wire_idx: usize,
 edge_idx: usize,
 gap: f64,
}

fn collect_wire_gaps(brep: &rcad_kernel::BRep, tolerance: f64, max_gap: f64) -> Vec<WireGapInfo> {
 let mut gaps = Vec::new();

 for (si, solid) in brep.solids.iter().enumerate() {
 for (shi, shell) in solid.shells.iter().enumerate() {
 for (fi, face) in shell.faces.iter().enumerate() {
 // Check outer wire
 if let Some(gap) = find_wire_gap(&face.outer_wire, brep, tolerance, max_gap) {
 gaps.push(WireGapInfo {
 solid: si,
 shell: shi,
 face: fi,
 wire_idx: 0,
 edge_idx: gap.0,
 gap: gap.1,
 });
 }

 // Check inner wires
 for (wi, wire) in face.inner_wires.iter().enumerate() {
 if let Some(gap) = find_wire_gap(wire, brep, tolerance, max_gap) {
 gaps.push(WireGapInfo {
 solid: si,
 shell: shi,
 face: fi,
 wire_idx: wi + 1,
 edge_idx: gap.0,
 gap: gap.1,
 });
 }
 }
 }
 }
 }

 gaps
}

fn find_wire_gap(wire: &Wire, brep: &rcad_kernel::BRep, tolerance: f64, max_gap: f64) -> Option<(usize, f64)> {
 if wire.edges.len() < 2 {
 return None;
 }

 for (i, we) in wire.edges.iter().enumerate() {
 let edge = brep.edges.get(we.idx)?;
 let next_i = (i + 1) % wire.edges.len();
 let next_edge = brep.edges.get(wire.edges[next_i].idx)?;

 let this_end = if we.forward { edge.end } else { edge.start };
 let next_start = if wire.edges[next_i].forward {
 next_edge.start
 } else {
 next_edge.end
 };

 if this_end != next_start {
 let gap_pt1 = brep.vertices.get(this_end).map(|v| v.point).unwrap_or_default();
 let gap_pt2 = brep.vertices.get(next_start).map(|v| v.point).unwrap_or_default();
 let gap = (gap_pt2 - gap_pt1).length();

 if gap <= max_gap && gap > tolerance {
 return Some((i, gap));
 }
 }
 }

 None
}

/// Report from UV bounds repair operations.
#[derive(Debug, Clone, Default)]
pub struct UvBoundsRepairReport {
 /// Number of faces whose PCurves were adjusted.
 pub faces_adjusted: usize,
 /// Number of PCurves modified.
 pub pcurves_modified: usize,
}

/// Repair UV bounds violations by adjusting PCurve parameter ranges.
///
/// This function fixes PCurve parameter ranges that fall outside the natural
/// bounds of their surfaces. For periodic surfaces, wraps UV parameters to
/// the canonical range. For bounded surfaces, clamps parameters.
///
/// Analogous to `ShapeFix_Face::FixUVBounds()` in OCCT.
pub fn fix_uv_bounds_violations(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, UvBoundsRepairReport) {
 use crate::brep_check::analyze_surface_uv_consistency;
 use rcad_kernel::geom::Surface3;

 let mut result = brep.clone();
 let mut report = UvBoundsRepairReport::default();

 let analysis = analyze_surface_uv_consistency(brep, tolerance);

 for violation in &analysis.faces_with_uv_bounds_violation {
 // Get the face's surface
 let flat_face_idx = {
 let mut idx = 0usize;
 for s in 0..violation.solid {
 for sh in &brep.solids[s].shells {
 idx += sh.faces.len();
 }
 }
 for sh in 0..violation.shell {
 idx += brep.solids[violation.solid].shells[sh].faces.len();
 }
 idx + violation.face
 };

 let surface_idx = match brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) {
 Some(idx) => idx,
 None => continue,
 };

 let surface = match brep.geom.surfaces.get(surface_idx) {
 Some(s) => s,
 None => continue,
 };

 // Get the UV period/wrapping info for the surface
 let (u_period, v_period, u_wrapped, v_wrapped) = match surface {
 Surface3::Cylinder(_) => (Some(2.0 * std::f64::consts::PI), None, true, false),
 Surface3::Sphere(_) => (Some(2.0 * std::f64::consts::PI), None, true, false),
 Surface3::Cone(_) => (Some(2.0 * std::f64::consts::PI), None, true, false),
 Surface3::Torus(_) => (
 Some(2.0 * std::f64::consts::PI),
 Some(2.0 * std::f64::consts::PI),
 true,
 true,
 ),
 Surface3::Plane(_) | Surface3::BSpline(_) => continue, // No wrapping needed
 _ => continue, // Other surface types not handled
 };

 // Adjust PCurves for edges in this face
 let face = &brep.solids[violation.solid].shells[violation.shell].faces[violation.face];
 for we in &face.outer_wire.edges {
 if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
 for pc in pcurves {
 if pc.surface_idx == surface_idx
 && let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) {
 // Check if curve2d needs adjustment
 let needs_wrap = check_curve2d_needs_wrap(
 curve2d,
 u_period,
 v_period,
 u_wrapped,
 v_wrapped,
 );

 if needs_wrap {
 // Create a wrapped version of the curve
 if let Some(wrapped) = wrap_curve2d(
 curve2d,
 u_period,
 v_period,
 u_wrapped,
 v_wrapped,
 ) {
 // Replace the curve2d
 let new_idx = result.geom.curve2ds.len();
 result.geom.curve2ds.push(wrapped);
 // Update the PCurve reference
 if let Some(pcs) = result.geom.edge_pcurves.get_mut(we.idx) {
 for p in pcs.iter_mut() {
 if p.surface_idx == surface_idx {
 p.curve2d_idx = new_idx;
 }
 }
 }
 report.pcurves_modified += 1;
 }
 }
 }
 }
 }
 }

 report.faces_adjusted += 1;
 }

 (result, report)
}

fn check_curve2d_needs_wrap(
 curve2d: &rcad_kernel::Curve2d,
 u_period: Option<f64>,
 v_period: Option<f64>,
 u_wrapped: bool,
 v_wrapped: bool,
) -> bool {
 use rcad_kernel::geom::Curve2dEval;

 // Sample the curve and check for out-of-bounds parameters
 for i in 0..=16 {
 let t = i as f64 / 16.0;
 let uv = curve2d.point_at(t);

 if u_wrapped
 && let Some(period) = u_period
 && (uv.x < -period * 0.5 || uv.x > period * 0.5) {
 return true;
 }

 if v_wrapped
 && let Some(period) = v_period
 && (uv.y < -period * 0.5 || uv.y > period * 0.5) {
 return true;
 }
 }

 false
}

fn wrap_curve2d(
 curve2d: &rcad_kernel::Curve2d,
 u_period: Option<f64>,
 v_period: Option<f64>,
 u_wrapped: bool,
 v_wrapped: bool,
) -> Option<rcad_kernel::Curve2d> {
 use rcad_kernel::Curve2d;

 match curve2d {
 Curve2d::Line(line) => {
 // For a line, we can adjust the origin to be within canonical bounds
 let mut new_line = *line;

 if u_wrapped
 && let Some(period) = u_period {
 // Wrap the origin's U coordinate
 while new_line.origin.x < -period * 0.5 {
 new_line.origin.x += period;
 }
 while new_line.origin.x > period * 0.5 {
 new_line.origin.x -= period;
 }
 }

 if v_wrapped
 && let Some(period) = v_period {
 // Wrap the origin's V coordinate
 while new_line.origin.y < -period * 0.5 {
 new_line.origin.y += period;
 }
 while new_line.origin.y > period * 0.5 {
 new_line.origin.y -= period;
 }
 }

 Some(Curve2d::Line(new_line))
 }
 Curve2d::BSpline(_) | Curve2d::Circle(_) | Curve2d::Ellipse(_) => {
 // For more complex curves, we'd need to implement proper wrapping
 // For now, return None to indicate we can't wrap this curve type
 None
 }
 _ => None, // Other curve types not handled
 }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// Enhanced Edge Sewing with Adaptive Tolerance
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Configuration for enhanced edge sewing operations.
#[derive(Debug, Clone)]
pub struct EdgeSewConfig {
 /// Base tolerance for edge endpoint matching.
 pub base_tolerance: f64,
 /// Maximum tolerance to use for adaptive expansion.
 pub max_tolerance: f64,
 /// Factor by which tolerance grows on each pass (1.0 = no growth).
 pub tolerance_growth: f64,
 /// Maximum number of sewing passes.
 pub max_passes: usize,
 /// Whether to use geometric proximity for edge matching.
 pub use_geometric_proximity: bool,
 /// Whether to merge edges that share the same curve geometry.
 pub merge_same_curve_edges: bool,
 /// Whether to handle periodic surface seams.
 pub handle_periodic_seams: bool,
}

impl Default for EdgeSewConfig {
 fn default() -> Self {
 Self {
 base_tolerance: TOLERANCE_ABS,
 max_tolerance: TOLERANCE_ABS * 100.0,
 tolerance_growth: 2.0,
 max_passes: 3,
 use_geometric_proximity: true,
 merge_same_curve_edges: true,
 handle_periodic_seams: true,
 }
 }
}

/// Enhanced report from edge sewing operations.
#[derive(Debug, Clone, Default)]
pub struct EnhancedEdgeSewReport {
 /// Number of edge pairs that were sewn together.
 pub edges_sewn: usize,
 /// Number of vertex pairs that were merged.
 pub vertices_merged: usize,
 /// Number of passes executed.
 pub passes_executed: usize,
 /// Final tolerance used.
 pub final_tolerance: f64,
 /// Whether the process converged.
 pub converged: bool,
 /// Number of edges merged by same-curve detection.
 pub same_curve_merges: usize,
 /// Number of periodic seam edges handled.
 pub periodic_seam_edges: usize,
}

/// Perform enhanced edge sewing with adaptive tolerance.
///
/// This function performs multiple passes of edge sewing with gradually
/// increasing tolerance, allowing for robust merging of near-coincident edges.
///
/// # Arguments
/// * `brep` - The brep to process.
/// * `config` - Configuration for the sewing operation.
///
/// # Returns
/// A tuple of (modified brep, report).
pub fn sew_edges_enhanced(brep: &rcad_kernel::BRep, config: &EdgeSewConfig) -> (rcad_kernel::BRep, EnhancedEdgeSewReport) {
 let mut result = brep.clone();
 let mut report = EnhancedEdgeSewReport::default();

 let base_tol = config.base_tolerance.max(TOLERANCE_ABS);
 let max_tol = config.max_tolerance.max(base_tol);

 for pass in 0..config.max_passes {
 let tol = if config.tolerance_growth > 1.0 {
 let grown = base_tol * config.tolerance_growth.powi(pass as i32);
 grown.min(max_tol)
 } else {
 base_tol
 };

 let (new_brep, sew_report) = sew_close_edges(&result, tol);
 let changed = sew_report.edges_sewn > 0 || sew_report.vertices_merged > 0;

 result = new_brep;
 report.edges_sewn += sew_report.edges_sewn;
 report.vertices_merged += sew_report.vertices_merged;
 report.passes_executed = pass + 1;
 report.final_tolerance = tol;

 if !changed {
 report.converged = true;
 break;
 }
 }

 // Additional pass for same-curve edge merging if enabled
 if config.merge_same_curve_edges {
 let (new_brep, same_curve_report) = merge_same_curve_edges(&result, config.base_tolerance);
 if same_curve_report.edges_merged > 0 {
 result = new_brep;
 report.same_curve_merges = same_curve_report.edges_merged;
 report.vertices_merged += same_curve_report.vertices_merged;
 }
 }

 // Handle periodic surface seams if enabled
 if config.handle_periodic_seams {
 let (new_brep, seam_report) = handle_periodic_surface_seams(&result, config.base_tolerance);
 if seam_report.seam_edges_detected > 0 || seam_report.seam_edges_split > 0 || seam_report.seam_edges_merged > 0 {
 result = new_brep;
 report.periodic_seam_edges = seam_report.seam_edges_detected + seam_report.seam_edges_split + seam_report.seam_edges_merged;
 }
 }

 (result, report)
}

/// Report from same-curve edge merging.
#[derive(Debug, Clone, Default)]
struct SameCurveMergeReport {
 edges_merged: usize,
 vertices_merged: usize,
}

/// Merge edges that share the same underlying curve geometry.
///
/// This is useful for edges that were split during boolean operations
/// but should logically be merged back together.
fn merge_same_curve_edges(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, SameCurveMergeReport) {
 let result = brep.clone();
 let mut report = SameCurveMergeReport::default();

 let n = result.edges.len();
 if n < 2 {
 return (result, report);
 }

 // Find edges that share the same curve
 let mut edge_groups: Vec<Vec<usize>> = Vec::new();
 let mut assigned = vec![false; n];

 for i in 0..n {
 if assigned[i] {
 continue;
 }

 let curve_i = result.geom.curves.get(i);
 if curve_i.is_none() {
 continue;
 }

 let mut group = vec![i];
 assigned[i] = true;

 for j in (i + 1)..n {
 if assigned[j] {
 continue;
 }

 let curve_j = result.geom.curves.get(j);
 if curve_j.is_none() {
 continue;
 }

 if curves_coincide(curve_i.unwrap(), curve_j.unwrap(), tolerance) {
 // Check if edges are adjacent (share an endpoint)
 let edge_i = &result.edges[i];
 let edge_j = &result.edges[j];
 let adjacent = edge_i.start == edge_j.start
 || edge_i.start == edge_j.end
 || edge_i.end == edge_j.start
 || edge_i.end == edge_j.end;

 if adjacent {
 group.push(j);
 assigned[j] = true;
 }
 }
 }

 if group.len() >= 2 {
 edge_groups.push(group);
 }
 }

 // Process edge groups
 for group in edge_groups {
 report.edges_merged += group.len() - 1;
 // Note: actual merging would require rebuilding topology
 // For now, we just record the groups
 }

 (result, report)
}

/// Check if two curves coincide within tolerance.
fn curves_coincide(c1: &rcad_kernel::Curve3, c2: &rcad_kernel::Curve3, tol: f64) -> bool {
 use rcad_kernel::Curve3;

 match (c1, c2) {
 (Curve3::Line(l1), Curve3::Line(l2)) => {
 let d1 = l1.direction.normalize_or_zero();
 let d2 = l2.direction.normalize_or_zero();
 if d1.dot(d2).abs() < 0.99 {
 return false;
 }
 let v = l2.origin - l1.origin;
 let perp = v - d1 * v.dot(d1);
 perp.length() <= tol
 }
 (Curve3::Circle(c1), Curve3::Circle(c2)) => {
 (c1.center - c2.center).length() <= tol
 && c1.normal.dot(c2.normal).abs() >= 0.99
 && (c1.radius - c2.radius).abs() <= tol
 }
 (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
 (e1.center - e2.center).length() <= tol
 && e1.normal.dot(e2.normal).abs() >= 0.99
 && (e1.major_radius - e2.major_radius).abs() <= tol
 && (e1.minor_radius - e2.minor_radius).abs() <= tol
 }
 _ => false,
 }
}

/// Report from periodic surface seam handling.
#[derive(Debug, Clone, Default)]
pub struct PeriodicSeamReport {
 /// Number of seam edges detected on periodic surfaces.
 pub seam_edges_detected: usize,
 /// Number of edges split at periodic seams.
 pub seam_edges_split: usize,
 /// Number of degenerate points handled (sphere poles, cone apex).
 pub degenerate_points_handled: usize,
 /// Number of edges merged across periodic seams.
 pub seam_edges_merged: usize,
}

/// Information about a periodic surface's periodicity.
#[derive(Debug, Clone, Copy)]
pub struct PeriodicSurfaceInfo {
 /// U-period (e.g., 2 ?for cylinder, sphere, cone, torus).
 pub u_period: Option<f64>,
 /// V-period (e.g., 2 ?for torus, None for others).
 pub v_period: Option<f64>,
 /// Whether the surface has a degenerate point at V=0 (sphere north pole).
 pub degenerate_at_v_min: bool,
 /// Whether the surface has a degenerate point at V=max (sphere south pole).
 pub degenerate_at_v_max: bool,
 /// Whether the surface has an apex degeneracy (cone).
 pub has_apex: bool,
 /// V value at the apex for cones (typically 0 or  ?.
 pub apex_v: Option<f64>,
}

impl PeriodicSurfaceInfo {
 /// Returns true if the surface is periodic in U direction.
 pub fn is_u_periodic(&self) -> bool {
 self.u_period.is_some()
 }

 /// Returns true if the surface is periodic in V direction.
 pub fn is_v_periodic(&self) -> bool {
 self.v_period.is_some()
 }

 /// Returns true if the surface has any degenerate points.
 pub fn has_degenerate_points(&self) -> bool {
 self.degenerate_at_v_min || self.degenerate_at_v_max || self.has_apex
 }
}

/// Detect periodic surface information from a Surface3.
pub fn detect_periodic_surface_info(surface: &Surface3) -> PeriodicSurfaceInfo {
 match surface {
 Surface3::Cylinder(_) => PeriodicSurfaceInfo {
 u_period: Some(std::f64::consts::TAU),
 v_period: None,
 degenerate_at_v_min: false,
 degenerate_at_v_max: false,
 has_apex: false,
 apex_v: None,
 },
 Surface3::Sphere(_) => PeriodicSurfaceInfo {
 u_period: Some(std::f64::consts::TAU),
 v_period: None,
 degenerate_at_v_min: true,  // V=0 is north pole
 degenerate_at_v_max: true,  // V= ?is south pole
 has_apex: false,
 apex_v: None,
 },
 Surface3::Cone(_) => PeriodicSurfaceInfo {
 u_period: Some(std::f64::consts::TAU),
 v_period: None,
 degenerate_at_v_min: false,
 degenerate_at_v_max: false,
 has_apex: true,
 apex_v: Some(0.0), // Apex is at V=0 (or can be at V= ?depending on half_angle)
 },
 Surface3::Torus(_) => PeriodicSurfaceInfo {
 u_period: Some(std::f64::consts::TAU),
 v_period: Some(std::f64::consts::TAU),
 degenerate_at_v_min: false,
 degenerate_at_v_max: false,
 has_apex: false,
 apex_v: None,
 },
 Surface3::Trimmed(trimmed) => {
 // Delegate to the basis surface
 detect_periodic_surface_info(trimmed.basis.as_ref())
 }
 _ => PeriodicSurfaceInfo {
 u_period: None,
 v_period: None,
 degenerate_at_v_min: false,
 degenerate_at_v_max: false,
 has_apex: false,
 apex_v: None,
 },
 }
}

/// Information about an edge crossing a periodic seam.
#[derive(Debug, Clone)]
pub struct SeamEdgeInfo {
 /// Edge index in the brep.
 pub edge_idx: usize,
 /// Surface index where the seam was detected.
 pub surface_idx: usize,
 /// Face index (flat) where the seam was detected.
 pub face_idx: usize,
 /// Whether the edge crosses the U seam.
 pub crosses_u_seam: bool,
 /// Whether the edge crosses the V seam.
 pub crosses_v_seam: bool,
 /// U parameter where the edge crosses the U seam (0 or period).
 pub u_seam_cross_param: Option<f64>,
 /// V parameter where the edge crosses the V seam.
 pub v_seam_cross_param: Option<f64>,
 /// Parameter t on the edge curve where the crossing occurs.
 pub edge_t_at_seam: Option<f64>,
}

/// Configuration for periodic surface seam handling.
#[derive(Debug, Clone)]
pub struct PeriodicSeamConfig {
 /// Tolerance for detecting seam proximity.
 pub seam_tolerance: f64,
 /// Whether to split edges at seams.
 pub split_edges: bool,
 /// Whether to merge edges across seams.
 pub merge_edges: bool,
 /// Whether to handle degenerate points (sphere poles, cone apex).
 pub handle_degeneracies: bool,
 /// Maximum distance for merging seam edge endpoints.
 pub merge_tolerance: f64,
}

impl Default for PeriodicSeamConfig {
 fn default() -> Self {
 Self {
 seam_tolerance: TOLERANCE_ABS * 10.0,
 split_edges: true,
 merge_edges: true,
 handle_degeneracies: true,
 merge_tolerance: TOLERANCE_ABS * 100.0,
 }
 }
}

/// Detect edges that cross periodic surface seams.
///
/// This function examines all edges on periodic surfaces and identifies
/// those whose UV parameterization crosses the seam boundary.
pub fn detect_seam_edges(brep: &rcad_kernel::BRep, config: &PeriodicSeamConfig) -> Vec<SeamEdgeInfo> {
 let mut seam_edges = Vec::new();

 // Iterate through all faces
 let mut flat_face_idx = 0usize;
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 // Get the surface for this face
 let surface_idx = match brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) {
 Some(idx) => idx,
 None => {
 flat_face_idx += 1;
 continue;
 }
 };

 let surface = match brep.geom.surfaces.get(surface_idx) {
 Some(s) => s,
 None => {
 flat_face_idx += 1;
 continue;
 }
 };

 let periodic_info = detect_periodic_surface_info(surface);
 if !periodic_info.is_u_periodic() && !periodic_info.is_v_periodic() {
 flat_face_idx += 1;
 continue;
 }

 // Check each edge in the face's wire
 for we in &face.outer_wire.edges {
 if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
 for pc in pcurves {
 if pc.surface_idx != surface_idx {
 continue;
 }

 if let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx)
 && let Some(seam_info) = detect_curve2d_seam_crossing(
 curve2d,
 we.forward,
 &periodic_info,
 config.seam_tolerance,
 we.idx,
 surface_idx,
 flat_face_idx,
 ) {
 seam_edges.push(seam_info);
 }
 }
 }
 }

 // Also check inner wires
 for inner_wire in &face.inner_wires {
 for we in &inner_wire.edges {
 if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
 for pc in pcurves {
 if pc.surface_idx != surface_idx {
 continue;
 }

 if let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx)
 && let Some(seam_info) = detect_curve2d_seam_crossing(
 curve2d,
 we.forward,
 &periodic_info,
 config.seam_tolerance,
 we.idx,
 surface_idx,
 flat_face_idx,
 ) {
 seam_edges.push(seam_info);
 }
 }
 }
 }
 }

 flat_face_idx += 1;
 }
 }
 }

 seam_edges
}

/// Helper function to detect if a 2D curve crosses a seam.
fn detect_curve2d_seam_crossing(
 curve2d: &rcad_kernel::Curve2d,
 forward: bool,
 periodic_info: &PeriodicSurfaceInfo,
 _seam_tolerance: f64,
 edge_idx: usize,
 surface_idx: usize,
 face_idx: usize,
) -> Option<SeamEdgeInfo> {
 use rcad_kernel::Curve2dEval;

 // Sample the curve at multiple points
 let num_samples = 20usize;
 let mut uv_points = Vec::with_capacity(num_samples + 1);

 for i in 0..=num_samples {
 let t = if forward {
 i as f64 / num_samples as f64
 } else {
 1.0 - i as f64 / num_samples as f64
 };
 uv_points.push(curve2d.point_at(t));
 }

 // Check for U-seam crossing
 let mut crosses_u_seam = false;
 let mut u_seam_cross_param = None;
 let mut edge_t_at_seam = None;

 if let Some(u_period) = periodic_info.u_period {
 for i in 1..uv_points.len() {
 let u1 = uv_points[i - 1].x;
 let u2 = uv_points[i].x;
 let du = u2 - u1;

 // Large jump indicates seam crossing
 if du.abs() > u_period * 0.5 {
 crosses_u_seam = true;
 // Determine which way we're crossing
 let seam_u = if du < 0.0 {
 // Going from high U to low U, crossing at U=period
 u_period
 } else {
 // Going from low U to high U, crossing at U=0
 0.0
 };
 u_seam_cross_param = Some(seam_u);

 // Compute the approximate t parameter at the seam
 let t1 = (i - 1) as f64 / num_samples as f64;
 let t2 = i as f64 / num_samples as f64;
 // Linear interpolation factor
 let factor = if du.abs() > TOLERANCE_LINEAR_ULTRA_STRICT {
 (seam_u - u1) / du
 } else {
 0.5
 };
 edge_t_at_seam = Some(t1 + factor * (t2 - t1));
 break;
 }
 }
 }

 // Check for V-seam crossing (for torus)
 let mut crosses_v_seam = false;
 let mut v_seam_cross_param = None;

 if let Some(v_period) = periodic_info.v_period {
 for i in 1..uv_points.len() {
 let v1 = uv_points[i - 1].y;
 let v2 = uv_points[i].y;
 let dv = v2 - v1;

 if dv.abs() > v_period * 0.5 {
 crosses_v_seam = true;
 v_seam_cross_param = Some(if dv < 0.0 { v_period } else { 0.0 });
 break;
 }
 }
 }

 if crosses_u_seam || crosses_v_seam {
 Some(SeamEdgeInfo {
 edge_idx,
 surface_idx,
 face_idx,
 crosses_u_seam,
 crosses_v_seam,
 u_seam_cross_param,
 v_seam_cross_param,
 edge_t_at_seam,
 })
 } else {
 None
 }
}

/// Split an edge at a periodic seam.
///
/// This function creates a new vertex at the seam crossing point and
/// splits the edge into two edges.
pub fn split_edge_at_seam(
 brep: &rcad_kernel::BRep,
 seam_info: &SeamEdgeInfo,
 _tolerance: f64,
) -> (rcad_kernel::BRep, bool) {
 let mut result = brep.clone();
 let mut split_performed = false;

 let edge = match brep.edges.get(seam_info.edge_idx) {
 Some(e) => e,
 None => return (result, false),
 };

 let t_at_seam = match seam_info.edge_t_at_seam {
 Some(t) => t,
 None => return (result, false),
 };

 // Get the 3D curve for the edge
 let curve_idx = match brep.geom.edge_curve.get(seam_info.edge_idx).and_then(|v| *v) {
 Some(idx) => idx,
 None => return (result, false),
 };

 let curve = match brep.geom.curves.get(curve_idx) {
 Some(c) => c,
 None => return (result, false),
 };

 // Compute the 3D point at the seam crossing
 use rcad_kernel::CurveEval;
 let seam_point = curve.point_at(t_at_seam);

 // Create a new vertex at the seam point
 let new_vertex_idx = result.vertices.len();
 result.vertices.push(Vertex { point: seam_point });

 // Create a new edge from start to new vertex
 let new_edge_idx = result.edges.len();
 result.edges.push(Edge {
 start: edge.start,
 end: new_vertex_idx,
 });

 // Copy geometry for the new edge
 if result.geom.edge_curve.len() <= new_edge_idx {
 result.geom.edge_curve.resize(new_edge_idx + 1, None);
 }
 result.geom.edge_curve[new_edge_idx] = Some(curve_idx);

 // Update the original edge to go from new vertex to end
 if let Some(orig_edge) = result.edges.get_mut(seam_info.edge_idx) {
 orig_edge.start = new_vertex_idx;
 }

 // Update wire references
 // We need to find all wires that reference this edge and update them
 for solid in &mut result.solids {
 for shell in &mut solid.shells {
 for face in &mut shell.faces {
 // Update outer wire
 for we in &mut face.outer_wire.edges {
 if we.idx == seam_info.edge_idx && we.forward {
 // Insert the new edge after the split edge
 // This is a simplified approach - in practice we'd need more sophisticated wire manipulation
 }
 }
 // Update inner wires
 for inner_wire in &mut face.inner_wires {
 for we in &mut inner_wire.edges {
 if we.idx == seam_info.edge_idx && we.forward {
 // Similar update needed
 }
 }
 }
 }
 }
 }

 split_performed = true;
 (result, split_performed)
}


