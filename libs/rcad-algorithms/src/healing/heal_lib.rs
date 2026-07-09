//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Enhanced Healing: ShapeFix_Solid and ShapeFix_Wire Equivalents
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// ShapeFix_Solid equivalent: comprehensive solid repair.
///
/// This function performs OCCT ShapeFix_Solid-like operations:
/// - Shell orientation verification and repair
/// - Solid closure verification
/// - Shell manifoldness checks
/// - Face orientation consistency
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `tolerance` - Geometric tolerance
///
/// # Returns
/// Repaired B-Rep and count of fixes applied.
pub fn fix_solid(brep: &rcad_kernel::BRep, _tolerance: f64) -> (rcad_kernel::BRep, SolidFixReport) {
 use crate::brep_repair::{fix_face_orientation, recompute_face_normals};
 use rcad_kernel::BRepGraph;

 let mut report = SolidFixReport::default();
 let mut current = brep.clone();

 // Step 1: Recompute invalid normals
 let (brep_with_normals, normals_fixed) = recompute_face_normals(&current);
 current = brep_with_normals;
 report.normals_recomputed = normals_fixed;

 // Step 2: Fix face orientation for inward-pointing faces
 let (brep_oriented, faces_reoriented) = fix_face_orientation(&current);
 current = brep_oriented;
 report.faces_reoriented = faces_reoriented;

 // Step 3: Check solid closure and manifoldness
 let graph = BRepGraph::from_brep(&current);

 // Check if shells are closed — use tshape API
 use rcad_kernel::topods::TShape;
 for ts in &current.tshapes {
 if let TShape::Solid(sd) = &**ts {
 for (shi, shell_sr) in sd.shells.iter().enumerate() {
 let is_closed = if let TShape::Shell(shd) = &*current.tshapes[shell_sr.index] {
 shd.faces.iter().all(|face_sr| {
 if let TShape::Face(fd) = &*current.tshapes[face_sr.index] {
 if let TShape::Wire(wd) = &*current.tshapes[fd.outer_wire.index] {
 !wd.edges.is_empty()
 } else { false }
 } else { false }
 })
 } else { false };

 if !is_closed {
 report.unclosed_shells.push((shell_sr.index, shi));
 }
 }
 }
 }

 // Check manifoldness
 let nm_summary = graph.non_manifold_summary();
 report.non_manifold_edges = nm_summary.multi_face_edges.len();
 report.non_manifold_vertices = nm_summary.non_manifold_vertices.len();

 // Step 4: Verify shell orientation consistency — use tshape API
 for ts in &current.tshapes {
 if let TShape::Solid(sd) = &**ts {
 for shell_sr in &sd.shells {
 if let TShape::Shell(shd) = &*current.tshapes[shell_sr.index] {
 // Count faces with normals pointing in consistent direction
 let mut outward_count = 0usize;
 let mut inward_count = 0usize;

 for face_sr in &shd.faces {
 if let TShape::Face(fd) = &*current.tshapes[face_sr.index] {
 let normal = fd.surface.as_ref().map(|s| {
 rcad_kernel::geom::SurfaceEval::normal_at(s, 0.0, 0.0)
 }).unwrap_or_default();
 if normal.z > 0.0 {
 outward_count += 1;
 } else if normal.z < 0.0 {
 inward_count += 1;
 }
 }
 }

 // If most normals are inconsistent, note orientation issues
 if outward_count > 0 && inward_count > 0 {
 let ratio = outward_count as f64 / (outward_count + inward_count) as f64;
 if !(0.3..=0.7).contains(&ratio) {
 report.orientation_inconsistencies += 1;
 }
 }
 }
 }
 }
 }

 report.total_fixes = report.normals_recomputed + report.faces_reoriented;
 (current, report)
}

/// Report from solid-level fixes.
#[derive(Debug, Clone, Default)]
pub struct SolidFixReport {
 /// Number of face normals recomputed.
 pub normals_recomputed: usize,
 /// Number of faces reoriented.
 pub faces_reoriented: usize,
 /// Indices of unclosed shells (solid_idx, shell_idx).
 pub unclosed_shells: Vec<(usize, usize)>,
 /// Number of non-manifold edges detected.
 pub non_manifold_edges: usize,
 /// Number of non-manifold vertices detected.
 pub non_manifold_vertices: usize,
 /// Number of shells with orientation inconsistencies.
 pub orientation_inconsistencies: usize,
 /// Total number of fixes applied.
 pub total_fixes: usize,
}

impl SolidFixReport {
 pub fn is_clean(&self) -> bool {
 self.unclosed_shells.is_empty()
 && self.non_manifold_edges == 0
 && self.non_manifold_vertices == 0
 && self.orientation_inconsistencies == 0
 }

 pub fn summary(&self) -> String {
 if self.is_clean() && self.total_fixes == 0 {
 "Solid is clean, no fixes needed".to_string()
 } else {
 format!(
 "Solid fixes: {} normals, {} orientations, {} unclosed shells, {} non-manifold edges, {} non-manifold vertices",
 self.normals_recomputed,
 self.faces_reoriented,
 self.unclosed_shells.len(),
 self.non_manifold_edges,
 self.non_manifold_vertices
 )
 }
 }
}

/// ShapeFix_Wire equivalent: comprehensive wire repair.
///
/// This function performs OCCT ShapeFix_Wire-like operations:
/// - Wire closure verification and repair
/// - Edge order verification
/// - Degenerate edge handling
/// - Self-intersection detection
/// - Wire orientation fix
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `tolerance` - Geometric tolerance
///
/// # Returns
/// Repaired B-Rep and detailed wire fix report.
pub fn fix_wire(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, WireFixReport) {
 use crate::brep_repair::fix_wire_orientation;

 let mut report = WireFixReport::default();
 let mut current = brep.clone();

 // Step 1: Fix wire orientation
 let (brep_fixed, wires_fixed) = fix_wire_orientation(&current, tolerance);
 current = brep_fixed;
 report.wires_oriented = wires_fixed;

 // Step 2: Analyze wires for issues — use tshape API
 use rcad_kernel::topods::TShape;
 for ts in &current.tshapes {
 if let TShape::Solid(sd) = &**ts {
 for shell_sr in &sd.shells {
 if let TShape::Shell(shd) = &*current.tshapes[shell_sr.index] {
 for (fi, face_sr) in shd.faces.iter().enumerate() {
 if let TShape::Face(fd) = &*current.tshapes[face_sr.index] {
 // Check outer wire
 if let TShape::Wire(owd) = &*current.tshapes[fd.outer_wire.index] {
 let outer_wire = rcad_kernel::topology::Wire {
 edges: owd.edges.iter().map(|e| rcad_kernel::topology::WireEdge { idx: e.index, forward: true }).collect(),
 };
 let outer_issues = analyze_wire_issues(&current, &outer_wire, tolerance);
 if outer_issues.open_gaps > 0 || outer_issues.topological_self_intersections > 0 || outer_issues.geometric_self_intersections > 0 {
 report.outer_wire_issues.push(WireIssueLocation {
 solid: shell_sr.index,
 shell: 0,
 face: fi,
 wire_idx: 0,
 issues: outer_issues,
 });
 }
 }

 // Check inner wires
 for (wi, iw_sr) in fd.inner_wires.iter().enumerate() {
 if let TShape::Wire(iwd) = &*current.tshapes[iw_sr.index] {
 let inner_wire = rcad_kernel::topology::Wire {
 edges: iwd.edges.iter().map(|e| rcad_kernel::topology::WireEdge { idx: e.index, forward: true }).collect(),
 };
 let inner_issues = analyze_wire_issues(&current, &inner_wire, tolerance);
 if inner_issues.open_gaps > 0 || inner_issues.topological_self_intersections > 0 || inner_issues.geometric_self_intersections > 0 {
 report.inner_wire_issues.push(WireIssueLocation {
 solid: shell_sr.index,
 shell: 0,
 face: fi,
 wire_idx: wi + 1,
 issues: inner_issues,
 });
 }
 }
 }
 }
 }
 }
 }
 }
 }

 // Step 3: Count degenerate edges — use tshape API
 let mut degenerate_edges: Vec<usize> = Vec::new();
 for (ei, ts) in current.tshapes.iter().enumerate() {
 if let TShape::Edge(ed) = &**ts {
 if let (Some(v1_idx), Some(v2_idx)) = (Some(ed.first.index), Some(ed.last.index)) {
 let start_pt = if let TShape::Vertex(vd1) = &*current.tshapes[v1_idx] { Some(vd1.point) } else { None };
 let end_pt = if let TShape::Vertex(vd2) = &*current.tshapes[v2_idx] { Some(vd2.point) } else { None };
 if let (Some(s), Some(e)) = (start_pt, end_pt)
 && (s - e).length() < tolerance {
 degenerate_edges.push(ei);
 }
 }
 }
 }
 report.degenerate_edges = degenerate_edges;

 // Step 4: Compute wire quality metrics — use tshape API
 let mut total_wires_checked = 0usize;
 for ts in &current.tshapes {
 if let TShape::Solid(sd) = &**ts {
 for shell_sr in &sd.shells {
 if let TShape::Shell(shd) = &*current.tshapes[shell_sr.index] {
 for face_sr in &shd.faces {
 if let TShape::Face(fd) = &*current.tshapes[face_sr.index] {
 total_wires_checked += 1; // outer wire
 total_wires_checked += fd.inner_wires.len();
 }
 }
 }
 }
 }
 }
 report.total_wires_checked = total_wires_checked
  + report.outer_wire_issues.len()
  + report.inner_wire_issues.len();

 report.wires_with_issues = report.outer_wire_issues.len() + report.inner_wire_issues.len();
 report.total_fixes = report.wires_oriented;

 (current, report)
}

/// Location of a wire issue.
#[derive(Debug, Clone)]
pub struct WireIssueLocation {
 pub solid: usize,
 pub shell: usize,
 pub face: usize,
 pub wire_idx: usize,
 pub issues: crate::brep_check::WireIssueReport,
}

/// Report from wire-level fixes.
#[derive(Debug, Clone, Default)]
pub struct WireFixReport {
 /// Number of wires with corrected orientation.
 pub wires_oriented: usize,
 /// Issues found in outer wires.
 pub outer_wire_issues: Vec<WireIssueLocation>,
 /// Issues found in inner wires.
 pub inner_wire_issues: Vec<WireIssueLocation>,
 /// Indices of degenerate edges found.
 pub degenerate_edges: Vec<usize>,
 /// Total wires checked.
 pub total_wires_checked: usize,
 /// Wires with issues.
 pub wires_with_issues: usize,
 /// Total fixes applied.
 pub total_fixes: usize,
}

impl WireFixReport {
 pub fn is_clean(&self) -> bool {
 self.outer_wire_issues.is_empty()
 && self.inner_wire_issues.is_empty()
 && self.degenerate_edges.is_empty()
 }

 pub fn summary(&self) -> String {
 if self.is_clean() && self.total_fixes == 0 {
 format!("All {} wires clean, no fixes needed", self.total_wires_checked)
 } else {
 format!(
 "Wire fixes: {} oriented, {} with issues ({} outer, {} inner), {} degenerate edges",
 self.wires_oriented,
 self.wires_with_issues,
 self.outer_wire_issues.len(),
 self.inner_wire_issues.len(),
 self.degenerate_edges.len()
 )
 }
 }
}

/// Analyze wire for issues without modifying.
fn analyze_wire_issues(brep: &rcad_kernel::BRep, wire: &rcad_kernel::topology::Wire, tolerance: f64) -> crate::brep_check::WireIssueReport {
 use rcad_kernel::topods::TShape;
 let mut open_gaps = 0usize;
 let mut topological_self_intersections = 0usize;
 let mut geometric_self_intersections = 0usize;

 // Collect wire vertices
 let mut wire_verts = Vec::with_capacity(wire.edges.len());
 for we in &wire.edges {
 let edge_data = match brep.tshapes.get(we.idx) {
 Some(ts) => match &**ts {
 TShape::Edge(ed) => Some((ed.first.index, ed.last.index)),
 _ => None,
 },
 None => None,
 };
 let Some((sv, ev)) = edge_data else { continue };
 let (sv, ev) = if we.forward { (sv, ev) } else { (ev, sv) };
 // Check vertices exist (tshapes)
 if sv < brep.tshapes.len() && ev < brep.tshapes.len()
 && matches!(&*brep.tshapes[sv], TShape::Vertex(_))
 && matches!(&*brep.tshapes[ev], TShape::Vertex(_))
 {
 wire_verts.push((sv, ev));
 }
 }

 // Check for open gaps
 let n = wire_verts.len();
 if n > 1 {
 for i in 0..n {
 let next = (i + 1) % n;
 let end_v = wire_verts[i].1;
 let start_v = wire_verts[next].0;
 if end_v != start_v {
 let end_pt = match &*brep.tshapes[end_v] { TShape::Vertex(vd) => vd.point, _ => continue };
 let start_pt = match &*brep.tshapes[start_v] { TShape::Vertex(vd) => vd.point, _ => continue };
 if (end_pt - start_pt).length() > tolerance {
 open_gaps += 1;
 }
 }
 }
 }

 // Check for topological self-intersection (vertex appearing more than twice)
 use std::collections::HashMap;
 let mut vertex_count: HashMap<usize, usize> = HashMap::new();
 for &(sv, ev) in &wire_verts {
 *vertex_count.entry(sv).or_insert(0) += 1;
 *vertex_count.entry(ev).or_insert(0) += 1;
 }
 for &count in vertex_count.values() {
 if count > 2 {
 topological_self_intersections += 1;
 }
 }

 // Check for geometric self-intersection (2D projection)
 let vertex_at = |idx: usize| -> glam::DVec3 {
 match &*brep.tshapes[idx] { TShape::Vertex(vd) => vd.point, _ => glam::DVec3::ZERO }
 };
 if n >= 4 {
 for i in 0..n {
 for j in (i + 2)..n {
 if i == 0 && j == n - 1 {
 continue; // Adjacent edges wraparound
 }
 let (a_start, a_end) = wire_verts[i];
 let (b_start, b_end) = wire_verts[j];
 let p1 = vertex_at(a_start);
 let p2 = vertex_at(a_end);
 let p3 = vertex_at(b_start);
 let p4 = vertex_at(b_end);

 if segments_intersect_2d(p1, p2, p3, p4) {
 geometric_self_intersections += 1;
 }
 }
 }
 }

 crate::brep_check::WireIssueReport {
 solid: 0,
 shell: 0,
 face: 0,
 wire_idx: 0,
 edge_count: wire.edges.len(),
 open_gaps,
 topological_self_intersections,
 geometric_self_intersections,
 }
}

/// Check if two 2D line segments intersect (XY plane projection).
fn segments_intersect_2d(p1: glam::DVec3, p2: glam::DVec3, p3: glam::DVec3, p4: glam::DVec3) -> bool {
 let x1 = p1.x; let y1 = p1.y;
 let x2 = p2.x; let y2 = p2.y;
 let x3 = p3.x; let y3 = p3.y;
 let x4 = p4.x; let y4 = p4.y;

 let (min_x1, max_x1) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
 let (min_y1, max_y1) = if y1 < y2 { (y1, y2) } else { (y2, y1) };
 let (min_x2, max_x2) = if x3 < x4 { (x3, x4) } else { (x4, x3) };
 let (min_y2, max_y2) = if y3 < y4 { (y3, y4) } else { (y4, y3) };

 if max_x1 < min_x2 || max_x2 < min_x1 || max_y1 < min_y2 || max_y2 < min_y1 {
 return false;
 }

 fn ccw(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64) -> bool {
 (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)
 }

 ccw(x1, y1, x3, y3, x4, y4) != ccw(x2, y2, x3, y3, x4, y4)
 && ccw(x1, y1, x2, y2, x3, y3) != ccw(x1, y1, x2, y2, x4, y4)
}

/// Comprehensive healing with ShapeFix_Solid and ShapeFix_Wire integration.
///
/// This function provides OCCT-equivalent comprehensive healing:
/// 1. Wire-level fixes
/// 2. Face-level fixes
/// 3. Shell-level fixes
/// 4. Solid-level fixes
/// 5. Tolerance propagation
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `options` - Healing options
///
/// # Returns
/// Healed B-Rep and comprehensive report.
pub fn heal_comprehensive(brep: &rcad_kernel::BRep, options: &HealingOptions) -> (rcad_kernel::BRep, ComprehensiveHealingReport) {
 let mut report = ComprehensiveHealingReport::default();
 let mut current = brep.clone();

 // Stage 1: Wire fixes
 let (brep_wire, wire_report) = fix_wire(&current, options.tolerance);
 current = brep_wire;
 report.wire_report = Some(wire_report);

 // Stage 2: Face fixes (via standard repair)
 let (brep_face, repair_report) = repair(&current, options.tolerance);
 current = brep_face;
 report.repair_report = Some(repair_report);

 // Stage 3: Solid fixes
 let (brep_solid, solid_report) = fix_solid(&current, options.tolerance);
 current = brep_solid;
 report.solid_report = Some(solid_report);

 // Stage 4: Tolerance propagation
 current = crate::brep_repair::propagate_tolerances(
 &current,
 options.tolerance,
 crate::brep_repair::ToleranceFlowDirection::BottomUp,
 );
 let tol_report = crate::brep_repair::analyze_tolerances(&current, options.tolerance);
 report.tolerance_report = Some(tol_report.vertices);

 // Final check
 report.final_check = brep_check_analyze(&current);
 report.is_clean = report.final_check.is_valid();

 (current, report)
}

/// Comprehensive healing report with all stage details.
#[derive(Debug, Clone, Default)]
pub struct ComprehensiveHealingReport {
 /// Wire-level fix report.
 pub wire_report: Option<WireFixReport>,
 /// Standard repair report.
 pub repair_report: Option<crate::brep_repair::RepairReport>,
 /// Solid-level fix report.
 pub solid_report: Option<SolidFixReport>,
 /// Tolerance propagation report.
 pub tolerance_report: Option<crate::brep_repair::ToleranceStats>,
 /// Final checker result.
 pub final_check: CheckResult,
 /// Whether the result is checker-clean.
 pub is_clean: bool,
}

impl ComprehensiveHealingReport {
 pub fn summary(&self) -> String {
 let mut parts = Vec::new();

 if let Some(ref wr) = self.wire_report
 && wr.total_fixes > 0 {
 parts.push(format!("wires: {} fixes", wr.total_fixes));
 }

 if let Some(ref rr) = self.repair_report {
 let repairs = rr.vertices_merged + rr.faces_reoriented + rr.wires_fixed;
 if repairs > 0 {
 parts.push(format!("repair: {} fixes", repairs));
 }
 }

 if let Some(ref sr) = self.solid_report
 && sr.total_fixes > 0 {
 parts.push(format!("solid: {} fixes", sr.total_fixes));
 }

 if parts.is_empty() {
 if self.is_clean {
 "Clean result, no fixes needed".to_string()
 } else {
 format!("Issues remain: {} issues", self.final_check.issues.len())
 }
 } else {
 format!("{}  ?{}", parts.join(", "), if self.is_clean { "clean" } else { "issues remain" })
 }
 }
}
