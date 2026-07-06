//! BRep validity checker.
//!
//! Analogous to OCCT `BRepCheck_Analyzer`. Checks structural and geometric
//! consistency of a BRep without modifying it.
//!
//! # Checks performed
//!
//! - **C1 Wire closure**: every wire must form a closed chain  ?the end vertex of
//! each edge must equal the start vertex of the next edge.
//! - **C2 Face normal consistency**: each face's stored normal must not be a zero
//! vector.
//! - **C3 Degenerate face**: faces with fewer than 3 wire edges are degenerate.
//! - **C4 Edge index validity**: WireEdge indices must be within bounds of
//! `brep.edges`.
//! - **C5 Vertex index validity**: each edge's start/end indices must be within
//! bounds of `brep.vertices`.
//! - **C6 Manifold topology**: each edge must be shared by exactly 2 faces
//! (for closed manifold solids).
//! - **C7 Wire self-intersection**: a wire's edges must not share vertices
//! except at consecutive junctions (no figure-8 or self-touching wires).
//!
//! # Extended checks (OCCT BRepCheck_Analyzer equivalent)
//!
//! - **Surface continuity**: C0, C1, C2 continuity across adjacent faces
//! - **Curve-surface consistency**: 3D curve endpoints match surface evaluation
//! - **Edge-curve tolerance verification**: edge tolerance covers geometry deviation
//! - **Face-surface tolerance verification**: face tolerance covers surface deviation
//! - **Shell orientation consistency**: consistent normal orientation in shells
//! - **Solid closure verification**: all edges shared by exactly 2 faces
//! - **Wire orientation**: clockwise vs counter-clockwise validation
//! - **Nested wire validation**: inner loops properly contained within outer
//! - **Tolerance consistency**: adjacent faces have compatible tolerances
//! - **Vertex tolerance propagation**: vertices have appropriate tolerances
//! - **Aspect ratio checks**: face quality metrics
//! - **Degenerate geometry detection**: zero-length edges, collapsed faces
//! - **Sliver face detection**: very thin triangular faces
//! - **Small feature detection**: tiny faces, edges, vertices

use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::{topods, BRep};
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};

/// A single validity issue found during checking.
#[derive(Debug, Clone, PartialEq)]
pub enum CheckIssue {
 /// Wire is not closed: end vertex of edge `edge_idx` does not match start
 /// vertex of the next edge in the wire (solid `solid`, shell `shell`,
 /// face `face`, position `wire_pos`).
 OpenWire {
 solid: usize,
 shell: usize,
 face: usize,
 /// Index of the edge within the wire where the gap occurs.
 wire_pos: usize,
 },
 /// Face normal is a zero vector.
 ZeroNormal {
 solid: usize,
 shell: usize,
 face: usize,
 },
 /// Face outer wire has fewer than 3 edges.
 DegenerateFace {
 solid: usize,
 shell: usize,
 face: usize,
 },
 /// A WireEdge references an edge index that is out of bounds.
 InvalidEdgeIndex {
 solid: usize,
 shell: usize,
 face: usize,
 edge_idx: usize,
 },
 /// An edge references a vertex index that is out of bounds.
 InvalidVertexIndex { edge: usize, vertex_idx: usize },
 /// An edge is shared by more or fewer than 2 faces (non-manifold).
 NonManifoldEdge { edge_idx: usize, face_count: usize },
 /// A wire has self-intersecting topology: a vertex appears more than
 /// twice (once as start, once as end) in the same wire.
 SelfIntersectingWire {
 solid: usize,
 shell: usize,
 face: usize,
 wire_idx: usize,
 vertex: usize,
 },
 /// A wire's outer boundary edges intersect each other geometrically
 /// (non-adjacent edges in the wire cross in 3D space).
 ///
 /// This catches cases where a face wire forms a figure-eight or butterfly
 /// polygon rather than a simple closed loop.
 GeometricSelfIntersection {
 solid: usize,
 shell: usize,
 face: usize,
 /// Index of one of the crossing edges within the outer wire.
 edge_a: usize,
 /// Index of the other crossing edge within the outer wire.
 edge_b: usize,
 },
 //  € € Geometry validation issues  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
 /// Surface continuity violation between adjacent faces.
 SurfaceContinuityViolation {
 solid: usize,
 face_a: usize,
 face_b: usize,
 shared_edge: usize,
 /// Expected continuity (0=C0, 1=C1, 2=C2)
 expected: u8,
 /// Actual continuity achieved
 actual: u8,
 /// Gap or angle deviation at the junction
 deviation: f64,
 },
 /// Curve-surface consistency violation: 3D curve doesn't match surface evaluation.
 CurveSurfaceMismatch {
 edge: usize,
 surface: usize,
 /// Maximum deviation between 3D curve and surface curve
 max_deviation: f64,
 },
 /// Edge tolerance insufficient to cover geometry deviation.
 EdgeToleranceViolation {
 edge: usize,
 stored_tolerance: f64,
 required_tolerance: f64,
 },
 /// Face tolerance insufficient to cover surface deviation.
 FaceToleranceViolation {
 solid: usize,
 shell: usize,
 face: usize,
 stored_tolerance: f64,
 required_tolerance: f64,
 },
 //  € € Topology validation issues  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
 /// Shell has inconsistent orientation (mixed inward/outward normals).
 ShellOrientationInconsistent {
 solid: usize,
 shell: usize,
 faces_with_inverted_normals: usize,
 },
 /// Solid is not closed (has boundary edges).
 SolidNotClosed {
 solid: usize,
 boundary_edge_count: usize,
 },
 /// Wire orientation is incorrect for its role (outer vs inner).
 WireOrientationIncorrect {
 solid: usize,
 shell: usize,
 face: usize,
 wire_idx: usize,
 /// true = should be CCW (outer), false = should be CW (inner)
 expected_ccw: bool,
 actual_ccw: bool,
 },
 /// Inner wire is not properly contained within outer wire.
 NestedWireViolation {
 solid: usize,
 shell: usize,
 face: usize,
 inner_wire_idx: usize,
 /// Number of inner wire vertices outside outer wire boundary
 vertices_outside: usize,
 },
 //  € € Tolerance issues  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
 /// Adjacent faces have inconsistent tolerances.
 ToleranceInconsistency {
 edge: usize,
 face_a: usize,
 face_b: usize,
 tolerance_a: f64,
 tolerance_b: f64,
 ratio: f64,
 },
 /// Vertex tolerance doesn't cover incident edge endpoints.
 VertexToleranceViolation {
 vertex: usize,
 stored_tolerance: f64,
 required_tolerance: f64,
 },
 //  € € Quality metric issues  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
 /// Face has poor aspect ratio.
 PoorAspectRatio {
 solid: usize,
 shell: usize,
 face: usize,
 aspect_ratio: f64,
 },
 /// Edge has near-zero length.
 DegenerateEdge {
 edge: usize,
 length: f64,
 },
 /// Face is a sliver (very thin).
 SliverFace {
 solid: usize,
 shell: usize,
 face: usize,
 area: f64,
 min_dimension: f64,
 },
 /// Small feature detected (tiny face or edge).
 SmallFeature {
 solid: usize,
 shell: usize,
 face: usize,
 feature_type: SmallFeatureType,
 size: f64,
 },
}

/// Type of small feature detected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SmallFeatureType {
 TinyFace,
 TinyEdge,
 TinyVertexGap,
}

impl std::fmt::Display for CheckIssue {
 fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
 match self {
 CheckIssue::OpenWire {
 solid,
 shell,
 face,
 wire_pos,
 } => write!(
 f,
 "OpenWire: solid={solid} shell={shell} face={face} at wire pos {wire_pos}"
 ),
 CheckIssue::ZeroNormal { solid, shell, face } => {
 write!(f, "ZeroNormal: solid={solid} shell={shell} face={face}")
 }
 CheckIssue::DegenerateFace { solid, shell, face } => {
 write!(f, "DegenerateFace: solid={solid} shell={shell} face={face}")
 }
 CheckIssue::InvalidEdgeIndex {
 solid,
 shell,
 face,
 edge_idx,
 } => write!(
 f,
 "InvalidEdgeIndex: solid={solid} shell={shell} face={face} edge={edge_idx}"
 ),
 CheckIssue::InvalidVertexIndex { edge, vertex_idx } => {
 write!(f, "InvalidVertexIndex: edge={edge} vertex={vertex_idx}")
 }
 CheckIssue::NonManifoldEdge { edge_idx, face_count } => {
 write!(f, "NonManifoldEdge: edge={edge_idx} shared by {face_count} faces (expected 2)")
 }
 CheckIssue::SelfIntersectingWire {
 solid,
 shell,
 face,
 wire_idx,
 vertex,
 } => write!(
 f,
 "SelfIntersectingWire: solid={solid} shell={shell} face={face} wire={wire_idx} vertex={vertex}"
 ),
 CheckIssue::GeometricSelfIntersection {
 solid,
 shell,
 face,
 edge_a,
 edge_b,
 } => write!(
 f,
 "GeometricSelfIntersection: solid={solid} shell={shell} face={face} edges {edge_a} and {edge_b} cross"
 ),
 // Geometry validation
 CheckIssue::SurfaceContinuityViolation {
 solid,
 face_a,
 face_b,
 shared_edge,
 expected,
 actual,
 deviation,
 } => write!(
 f,
 "SurfaceContinuityViolation: solid={solid} faces {face_a}/{face_b} edge={shared_edge} expected C{expected} got C{actual} deviation={deviation:.6e}"
 ),
 CheckIssue::CurveSurfaceMismatch { edge, surface, max_deviation } => {
 write!(f, "CurveSurfaceMismatch: edge={edge} surface={surface} deviation={max_deviation:.6e}")
 }
 CheckIssue::EdgeToleranceViolation { edge, stored_tolerance, required_tolerance } => {
 write!(f, "EdgeToleranceViolation: edge={edge} stored={stored_tolerance:.6e} required={required_tolerance:.6e}")
 }
 CheckIssue::FaceToleranceViolation { solid, shell, face, stored_tolerance, required_tolerance } => {
 write!(f, "FaceToleranceViolation: solid={solid} shell={shell} face={face} stored={stored_tolerance:.6e} required={required_tolerance:.6e}")
 }
 // Topology validation
 CheckIssue::ShellOrientationInconsistent { solid, shell, faces_with_inverted_normals } => {
 write!(f, "ShellOrientationInconsistent: solid={solid} shell={shell} {faces_with_inverted_normals} inverted faces")
 }
 CheckIssue::SolidNotClosed { solid, boundary_edge_count } => {
 write!(f, "SolidNotClosed: solid={solid} {boundary_edge_count} boundary edges")
 }
 CheckIssue::WireOrientationIncorrect { solid, shell, face, wire_idx, expected_ccw, actual_ccw } => {
 let expected = if *expected_ccw { "CCW" } else { "CW" };
 let actual = if *actual_ccw { "CCW" } else { "CW" };
 write!(f, "WireOrientationIncorrect: solid={solid} shell={shell} face={face} wire={wire_idx} expected={expected} got={actual}")
 }
 CheckIssue::NestedWireViolation { solid, shell, face, inner_wire_idx, vertices_outside } => {
 write!(f, "NestedWireViolation: solid={solid} shell={shell} face={face} inner_wire={inner_wire_idx} {vertices_outside} vertices outside")
 }
 // Tolerance issues
 CheckIssue::ToleranceInconsistency { edge, face_a, face_b, tolerance_a, tolerance_b, ratio } => {
 write!(f, "ToleranceInconsistency: edge={edge} faces {face_a}/{face_b} tol_a={tolerance_a:.6e} tol_b={tolerance_b:.6e} ratio={ratio:.2}")
 }
 CheckIssue::VertexToleranceViolation { vertex, stored_tolerance, required_tolerance } => {
 write!(f, "VertexToleranceViolation: vertex={vertex} stored={stored_tolerance:.6e} required={required_tolerance:.6e}")
 }
 // Quality metrics
 CheckIssue::PoorAspectRatio { solid, shell, face, aspect_ratio } => {
 write!(f, "PoorAspectRatio: solid={solid} shell={shell} face={face} ratio={aspect_ratio:.2}")
 }
 CheckIssue::DegenerateEdge { edge, length } => {
 write!(f, "DegenerateEdge: edge={edge} length={length:.6e}")
 }
 CheckIssue::SliverFace { solid, shell, face, area, min_dimension } => {
 write!(f, "SliverFace: solid={solid} shell={shell} face={face} area={area:.6e} min_dim={min_dimension:.6e}")
 }
 CheckIssue::SmallFeature { solid, shell, face, feature_type, size } => {
 let type_str = match feature_type {
 SmallFeatureType::TinyFace => "TinyFace",
 SmallFeatureType::TinyEdge => "TinyEdge",
 SmallFeatureType::TinyVertexGap => "TinyVertexGap",
 };
 write!(f, "SmallFeature: solid={solid} shell={shell} face={face} type={type_str} size={size:.6e}")
 }
 }
 }
}

/// Result of a BRep validity check.
#[derive(Debug, Clone)]
#[derive(Default)]
pub struct CheckResult {
 pub issues: Vec<CheckIssue>,
}

impl CheckResult {
 /// Returns `true` if no issues were found.
 pub fn is_valid(&self) -> bool {
 self.issues.is_empty()
 }
}

/// Check the validity of a BRep and return a `CheckResult` with any issues found.
///
/// Analogous to OCCT `BRepCheck_Analyzer::Perform()`.
pub fn brep_check_analyze(brep: &BRep) -> CheckResult {
 let mut issues = Vec::new();
 let n_edges = brep.edges.len();
 let n_verts = brep.vertices.len();

 // C5: edge vertex bounds
 for (eidx, edge) in brep.edges.iter().enumerate() {
 if edge.start >= n_verts {
 issues.push(CheckIssue::InvalidVertexIndex {
 edge: eidx,
 vertex_idx: edge.start,
 });
 }
 if edge.end >= n_verts {
 issues.push(CheckIssue::InvalidVertexIndex {
 edge: eidx,
 vertex_idx: edge.end,
 });
 }
 }

 // C6: manifold check  ?each edge must be shared by exactly 2 faces.
 // Count how many faces reference each edge across all solids/shells/faces.
 let mut edge_face_count: Vec<usize> = vec![0; n_edges];
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 // Count edges in outer wire
 for we in &face.outer_wire.edges {
 if we.idx < n_edges {
 edge_face_count[we.idx] += 1;
 }
 }
 // Count edges in inner wires
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if we.idx < n_edges {
 edge_face_count[we.idx] += 1;
 }
 }
 }
 }
 }
 }
 for (eidx, &count) in edge_face_count.iter().enumerate() {
 if count != 2 {
 issues.push(CheckIssue::NonManifoldEdge {
 edge_idx: eidx,
 face_count: count,
 });
 }
 }

 for (si, solid) in brep.solids.iter().enumerate() {
 for (shi, shell) in solid.shells.iter().enumerate() {
 for (fi, face) in shell.faces.iter().enumerate() {
 let wire = &face.outer_wire;

 // C2: zero normal
 if face.normal == DVec3::ZERO {
 issues.push(CheckIssue::ZeroNormal {
 solid: si,
 shell: shi,
 face: fi,
 });
 }

 // C3: degenerate face
 if wire.edges.len() < 3 {
 issues.push(CheckIssue::DegenerateFace {
 solid: si,
 shell: shi,
 face: fi,
 });
 // Can't check wire closure for degenerate face
 continue;
 }

 // C4: edge index bounds + collect start/end vertices for wire closure check
 let mut valid = true;
 let mut wire_verts: Vec<(usize, usize)> = Vec::new(); // (start_vidx, end_vidx)
 for we in &wire.edges {
 if we.idx >= n_edges {
 issues.push(CheckIssue::InvalidEdgeIndex {
 solid: si,
 shell: shi,
 face: fi,
 edge_idx: we.idx,
 });
 valid = false;
 } else {
 let edge = &brep.edges[we.idx];
 let (sv, ev) = if we.forward {
 (edge.start, edge.end)
 } else {
 (edge.end, edge.start)
 };
 wire_verts.push((sv, ev));
 }
 }

 if !valid {
 continue;
 }

 // C1: wire closure  ?end of edge[i] must match start of edge[i+1]
 // OCCT BRepCheck_Wire::Closed equivalent (BRepCheck_Wire.cxx lines ~60-95)
 let n = wire_verts.len();
 for i in 0..n {
 let next = (i + 1) % n;
 let end_v = wire_verts[i].1;
 let start_v = wire_verts[next].0;
 if end_v != start_v {
 // Tolerance check: allow same position even if different vertex objects
 let end_pt = brep.vertices[end_v].point;
 let start_pt = brep.vertices[start_v].point;
 if (end_pt - start_pt).length() > TOLERANCE_MESH_LEGACY {
 issues.push(CheckIssue::OpenWire {
 solid: si,
 shell: shi,
 face: fi,
 wire_pos: i,
 });
 }
 }
 }

 // C7: wire self-intersection  ?each vertex should appear at most
 // OCCT BRepCheck_Wire::SelfIntersection equivalent (BRepCheck_Wire.cxx lines ~100-145)
 // twice in the wire (once as start of an edge, once as end of another).
 check_wire_self_intersection(
 &wire_verts,
 &brep.vertices,
 si, shi, fi, 0, // outer wire index = 0
 &mut issues,
 );

 // C8: geometric self-intersection  ?check if non-adjacent edges of
 // the outer wire cross each other in 3D space (projects to 2D via
 // the face plane for planar faces).
 check_geometric_self_intersection(
 &wire_verts,
 &brep.vertices,
 si, shi, fi,
 &mut issues,
 );

 // Check inner wires too
 for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
 if inner_wire.edges.len() < 2 {
 continue; // too few edges to self-intersect
 }
 let mut inner_verts: Vec<(usize, usize)> = Vec::new();
 let mut inner_valid = true;
 for we in &inner_wire.edges {
 if we.idx >= n_edges {
 issues.push(CheckIssue::InvalidEdgeIndex {
 solid: si,
 shell: shi,
 face: fi,
 edge_idx: we.idx,
 });
 inner_valid = false;
 } else {
 let edge = &brep.edges[we.idx];
 let (sv, ev) = if we.forward {
 (edge.start, edge.end)
 } else {
 (edge.end, edge.start)
 };
 inner_verts.push((sv, ev));
 }
 }
 if !inner_valid {
 continue;
 }

 // Inner wire closure check
 let n_inner = inner_verts.len();
 for i in 0..n_inner {
 let next = (i + 1) % n_inner;
 let end_v = inner_verts[i].1;
 let start_v = inner_verts[next].0;
 if end_v != start_v {
 let end_pt = brep.vertices[end_v].point;
 let start_pt = brep.vertices[start_v].point;
 if (end_pt - start_pt).length() > TOLERANCE_MESH_LEGACY {
 issues.push(CheckIssue::OpenWire {
 solid: si,
 shell: shi,
 face: fi,
 wire_pos: i,
 });
 }
 }
 }

 // Inner wire self-intersection
 check_wire_self_intersection(
 &inner_verts,
 &brep.vertices,
 si, shi, fi,
 wi + 1, // inner wire indices start after outer
 &mut issues,
 );
 }
 }
 }
 }

 CheckResult { issues }
}

/// Convenience short alias for [`brep_check_analyze`].
pub fn check_brep(brep: &BRep) -> CheckResult { brep_check_analyze(brep) }

/// Check a single wire for self-intersecting topology.
///
/// A valid wire wire should have each vertex appear at most twice across
/// all edge endpoints: once as the start of some edge and once as the end
/// of another edge. If a vertex appears 3+ times, the wire self-intersects.
///
/// Aligned with OCCT BRepCheck_Wire::SelfIntersection concept
/// (BRepCheck_Wire.cxx lines ~100-145).
fn check_wire_self_intersection(
 wire_verts: &[(usize, usize)],
 vertices: &[rcad_kernel::topology::Vertex],
 solid: usize,
 shell: usize,
 face: usize,
 wire_idx: usize,
 issues: &mut Vec<CheckIssue>,
) {
 use std::collections::HashMap;
 let mut vertex_count: HashMap<usize, usize> = HashMap::new();
 for &(sv, ev) in wire_verts {
 *vertex_count.entry(sv).or_insert(0) += 1;
 *vertex_count.entry(ev).or_insert(0) += 1;
 }
 // In a closed wire, each vertex should appear exactly twice (once as start,
 // once as end). Allow tolerance for vertices at the same position.
 for (&vidx, &count) in &vertex_count {
 if count > 2 {
 // Check if it's actually a geometric self-intersection (different
 // positions) or just the same vertex referenced multiple times
 // (which could be valid for some topologies).
 if vidx < vertices.len() {
 issues.push(CheckIssue::SelfIntersectingWire {
 solid,
 shell,
 face,
 wire_idx,
 vertex: vidx,
 });
 }
 }
 }
}

/// Check whether non-adjacent edges of a wire intersect geometrically.
///
/// Projects the wire edge endpoints onto the face's 2D plane (using any two
/// non-collinear edges to form a local basis) and runs 2D segment intersection
/// tests on all non-adjacent edge pairs. Adjacent edges share an endpoint and
/// therefore trivially "intersect" at that endpoint  ?they are excluded.
///
/// This check only tests the start/end vertices; curved edges (circles, BSplines)
/// are approximated by their chord.
fn check_geometric_self_intersection(
 wire_verts: &[(usize, usize)],
 vertices: &[rcad_kernel::topology::Vertex],
 solid: usize,
 shell: usize,
 face: usize,
 issues: &mut Vec<CheckIssue>,
) {
 let n = wire_verts.len();
 if n < 4 {
 return; // Need at least 4 edges for non-adjacent crossing to exist
 }

 // Build the 3D segment list: each edge is (p0, p1).
 let segs: Vec<(DVec3, DVec3)> = wire_verts
 .iter()
 .map(|&(sv, ev)| {
 let p0 = vertices.get(sv).map(|v| v.point).unwrap_or(DVec3::ZERO);
 let p1 = vertices.get(ev).map(|v| v.point).unwrap_or(DVec3::ZERO);
 (p0, p1)
 })
 .collect();

 // Project to 2D using a local basis from the first non-degenerate edge.
 let (origin, axis_u, axis_v) = {
 let mut found = None;
 for i in 0..n {
 let d = segs[i].1 - segs[i].0;
 if d.length() > TOLERANCE_LEN_MIN {
 let u = d.normalize();
 // Pick an arbitrary perpendicular
 let tmp = if u.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let v = u.cross(tmp).normalize();
 found = Some((segs[i].0, u, v));
 break;
 }
 }
 match found {
 Some(b) => b,
 None => return, // All edges are degenerate; skip
 }
 };

 let project = |p: DVec3| -> [f64; 2] {
 let d = p - origin;
 [d.dot(axis_u), d.dot(axis_v)]
 };

 let seg2d: Vec<([f64; 2], [f64; 2])> =
 segs.iter().map(|&(p0, p1)| (project(p0), project(p1))).collect();

 // Check all non-adjacent pairs: i and j are non-adjacent when |i - j| > 1
 // (mod n, so also check the wrap-around adjacency).
 for i in 0..n {
 for j in (i + 2)..n {
 // Skip adjacent pair at the wrap-around (edge n-1 and edge 0 share a vertex)
 if i == 0 && j == n - 1 {
 continue;
 }
 if segments_2d_properly_intersect(seg2d[i].0, seg2d[i].1, seg2d[j].0, seg2d[j].1) {
 issues.push(CheckIssue::GeometricSelfIntersection {
 solid,
 shell,
 face,
 edge_a: i,
 edge_b: j,
 });
 return; // Report only the first crossing per face
 }
 }
 }
}

/// Returns `true` if the open segment p1 2 properly intersects segment p3 4.
/// Returns `false` if they only share an endpoint (T-intersection) or don't cross.
fn segments_2d_properly_intersect(
 p1: [f64; 2],
 p2: [f64; 2],
 p3: [f64; 2],
 p4: [f64; 2],
) -> bool {
 let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
 let d2 = [p4[0] - p3[0], p4[1] - p3[1]];

 let cross = d1[0] * d2[1] - d1[1] * d2[0];

 if cross.abs() < TOLERANCE_LEN_MIN {
 return false; // Parallel or collinear
 }

 let dx = p3[0] - p1[0];
 let dy = p3[1] - p1[1];
 let t = (dx * d2[1] - dy * d2[0]) / cross;
 let s = (dx * d1[1] - dy * d1[0]) / cross;

 // Proper interior intersection: t and s must be strictly in (0, 1)
 let eps = TOLERANCE_COORD_SUB;
 t > eps && t < 1.0 - eps && s > eps && s < 1.0 - eps
}

fn count_geometric_self_intersections(
 wire_verts: &[(usize, usize)],
 v_points: &[DVec3],
) -> usize {
 let n = wire_verts.len();
 if n < 4 {
 return 0;
 }

 let segs: Vec<(DVec3, DVec3)> = wire_verts
 .iter()
 .map(|&(sv, ev)| {
 let p0 = v_points.get(sv).copied().unwrap_or(DVec3::ZERO);
 let p1 = v_points.get(ev).copied().unwrap_or(DVec3::ZERO);
 (p0, p1)
 })
 .collect();

 let (origin, axis_u, axis_v) = {
 let mut found = None;
 for i in 0..n {
 let d = segs[i].1 - segs[i].0;
 if d.length() > TOLERANCE_LEN_MIN {
 let u = d.normalize();
 let tmp = if u.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let v = u.cross(tmp).normalize();
 found = Some((segs[i].0, u, v));
 break;
 }
 }
 match found {
 Some(b) => b,
 None => return 0,
 }
 };

 let project = |p: DVec3| -> [f64; 2] {
 let d = p - origin;
 [d.dot(axis_u), d.dot(axis_v)]
 };

 let seg2d: Vec<([f64; 2], [f64; 2])> =
 segs.iter().map(|&(p0, p1)| (project(p0), project(p1))).collect();

 let mut count = 0usize;
 for i in 0..n {
 for j in (i + 2)..n {
 if i == 0 && j == n - 1 {
 continue;
 }
 if segments_2d_properly_intersect(seg2d[i].0, seg2d[i].1, seg2d[j].0, seg2d[j].1) {
 count += 1;
 }
 }
 }
 count
}
// ----- OCCT BRepCheck alignment: Shell/Wire/Face validation -----
//
// Functions in this section are aligned with OCCT's BRepCheck classes:
// BRepCheck_Shell.cxx  (Shell closure - each edge shared by exactly 2 faces)
// BRepCheck_Wire.cxx (Wire closure + self-intersection)
// BRepCheck_Face.cxx (Wire-on-surface check)
//
// OCCT source: $OCCT_SRC/src/BRepCheck/
// -----------------------------------------------------------------

/// Find a face by its flat (global) index across all solids/shells.
///
/// Returns `(solid_idx, shell_idx, face_idx, &Face)` or `None` if the index
/// is out of range.
fn find_face_by_flat_idx<'a>(
 brep: &'a BRep,
 flat_idx: usize,
) -> Option<(usize, usize, usize, &'a rcad_kernel::topology::Face)> {
 let mut idx = 0usize;
 for (si, solid) in brep.solids.iter().enumerate() {
 for (shi, shell) in solid.shells.iter().enumerate() {
 for (fi, face) in shell.faces.iter().enumerate() {
 if idx == flat_idx {
 return Some((si, shi, fi, face));
 }
 idx += 1;
 }
 }
 }
 None
}

/// BRepCheck_Wire::Closed equivalent.
///
/// Checks that every wire belonging to the face at `face_idx` (flat index)
/// forms a closed loop: the end vertex of each edge matches the start vertex
/// of the next edge, and the last edge wraps back to the first.
///
/// Aligned with OCCT BRepCheck_Wire::Closed (BRepCheck_Wire.cxx lines ~60-95)
pub fn check_wire_closed(brep: &BRep, face_idx: usize) -> bool {
 let (_, _, _, face) = match find_face_by_flat_idx(brep, face_idx) {
 Some(f) => f,
 None => return false,
 };
 // Check outer wire
 if !check_single_wire_closed(brep, &face.outer_wire) {
 return false;
 }
 // Check all inner wires (holes)
 for wire in &face.inner_wires {
 if !check_single_wire_closed(brep, wire) {
 return false;
 }
 }
 true
}

/// Internal helper: checks closure of a single wire.
///
/// For each consecutive edge pair (including wrap-around), verifies that the
/// end vertex of edge[i] matches the start vertex of edge[i+1]. Falls back to
/// vertex position tolerance when vertex indices differ.
fn check_single_wire_closed(brep: &BRep, wire: &rcad_kernel::topology::Wire) -> bool {
 let n = wire.edges.len();
 if n == 0 {
 return true;
 }
 if n == 1 {
 let we = &wire.edges[0];
 if let Some(edge) = brep.edges.get(we.idx) {
 let (sv, ev) = if we.forward { (edge.start, edge.end) } else { (edge.end, edge.start) };
 if sv == ev { return true; }
 let s_pt = brep.vertices.get(sv).map(|v| v.point).unwrap_or_default();
 let e_pt = brep.vertices.get(ev).map(|v| v.point).unwrap_or_default();
 return (s_pt - e_pt).length() <= TOLERANCE_MESH_LEGACY;
 }
 return false;
 }
 for i in 0..n {
 let next = (i + 1) % n;
 let we_cur = &wire.edges[i];
 let we_next = &wire.edges[next];
 let edge_cur = match brep.edges.get(we_cur.idx) { Some(e) => e, None => return false };
 let edge_next = match brep.edges.get(we_next.idx) { Some(e) => e, None => return false };
 let (_, ev) = if we_cur.forward { (edge_cur.start, edge_cur.end) } else { (edge_cur.end, edge_cur.start) };
 let (sv, _) = if we_next.forward { (edge_next.start, edge_next.end) } else { (edge_next.end, edge_next.start) };
 if ev != sv {
 let end_pt = brep.vertices.get(ev).map(|v| v.point).unwrap_or_default();
 let start_pt = brep.vertices.get(sv).map(|v| v.point).unwrap_or_default();
 if (end_pt - start_pt).length() > TOLERANCE_MESH_LEGACY {
 return false;
 }
 }
 }
 true
}

/// BRepCheck_Wire::SelfIntersection equivalent.
///
/// Checks if any wire of the face at `face_idx` has edges that intersect
/// each other by sharing vertices at non-consecutive positions (i.e. a vertex
/// appears in more than two edge endpoint positions in the same wire).
///
/// Returns a list of `(edge_idx_in_wire_a, edge_idx_in_wire_b)` pairs for
/// each topological self-intersection found in any wire of the face.
///
/// Aligned with OCCT BRepCheck_Wire::SelfIntersection (BRepCheck_Wire.cxx lines ~100-145)
pub fn check_wire_self_intersection_pairs(
 brep: &BRep,
 face_idx: usize,
) -> Vec<(usize, usize)> {
 let (_, _, _, face) = match find_face_by_flat_idx(brep, face_idx) {
 Some(f) => f,
 None => return Vec::new(),
 };
 let mut result = Vec::new();
 result.extend(check_single_wire_self_intersection_pairs(brep, &face.outer_wire));
 for wire in &face.inner_wires {
 result.extend(check_single_wire_self_intersection_pairs(brep, wire));
 }
 result
}

/// Internal helper: finds self-intersecting edge pairs in a single wire.
fn check_single_wire_self_intersection_pairs(
 brep: &BRep,
 wire: &rcad_kernel::topology::Wire,
) -> Vec<(usize, usize)> {
 use std::collections::HashMap;
 let n = wire.edges.len();
 if n < 4 { return Vec::new(); }
 let mut vertex_occurrences: HashMap<usize, Vec<(usize, bool)>> = HashMap::new();
 for (i, we) in wire.edges.iter().enumerate() {
 if let Some(edge) = brep.edges.get(we.idx) {
 let (sv, ev) = if we.forward { (edge.start, edge.end) } else { (edge.end, edge.start) };
 vertex_occurrences.entry(sv).or_default().push((i, true));
 vertex_occurrences.entry(ev).or_default().push((i, false));
 }
 }
 let mut pairs = Vec::new();
 for (&_vidx, occurrences) in &vertex_occurrences {
 if occurrences.len() <= 2 { continue; }
 let edge_positions: Vec<usize> = occurrences.iter().map(|(pos, _)| *pos).collect();
 for a in 0..edge_positions.len() {
 for b in (a + 1)..edge_positions.len() {
 let ea = edge_positions[a];
 let eb = edge_positions[b];
 let diff = if ea > eb { ea - eb } else { eb - ea };
 let is_adjacent = diff == 1 || (ea == 0 && eb == n - 1) || (eb == 0 && ea == n - 1);
 if is_adjacent { continue; }
 let pair = if ea < eb { (ea, eb) } else { (eb, ea) };
 if !pairs.contains(&pair) { pairs.push(pair); }
 }
 }
 }
 pairs
}

/// BRepCheck_Face::Intersection equivalent (wire-on-surface check).
///
/// Checks that every edge in the face's wires lies on the face surface
/// within the given tolerance. For each edge with a 3D curve, samples 3
/// points along the curve (beginning, middle, end), projects them through
/// the face's PCurve + surface, and verifies that the deviation from the
/// 3D curve is within tolerance.
///
/// Returns `true` if all edge curves stay on the surface within tolerance.
///
/// Aligned with OCCT BRepCheck_Face::Intersection (BRepCheck_Face.cxx lines ~70-130)
pub fn check_face_wire_on_surface(
 brep: &BRep,
 face_idx: usize,
 tolerance: f64,
) -> bool {
 let (_, _, _, face) = match find_face_by_flat_idx(brep, face_idx) {
 Some(f) => f,
 None => return false,
 };
 let surface_idx = match brep.geom.face_surface.get(face_idx).and_then(|v| *v) {
 Some(idx) => idx,
 None => return true,
 };
 let surface = match brep.geom.surfaces.get(surface_idx) {
 Some(s) => s,
 None => return true,
 };
 for we in &face.outer_wire.edges {
 if !check_edge_on_surface(brep, we.idx, surface_idx, surface, tolerance) {
 return false;
 }
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if !check_edge_on_surface(brep, we.idx, surface_idx, surface, tolerance) {
 return false;
 }
 }
 }
 true
}

/// Check that a single edge's 3D curve lies on the given surface within tolerance.
fn check_edge_on_surface(
 brep: &BRep,
 edge_idx: usize,
 surface_idx: usize,
 surface: &rcad_kernel::geom::Surface3,
 tolerance: f64,
) -> bool {
 let curve_idx = match brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) {
 Some(idx) => idx, None => return true,
 };
 let curve = match brep.geom.curves.get(curve_idx) {
 Some(c) => c, None => return true,
 };
 let range = match brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) {
 Some(r) => r, None => return true,
 };
 let pcurves = match brep.geom.edge_pcurves.get(edge_idx) {
 Some(pc) => pc, None => return true,
 };
 let pc = match pcurves.iter().find(|pc| pc.surface_idx == surface_idx) {
 Some(pc) => pc, None => return true,
 };
 let curve2d = match brep.geom.curve2ds.get(pc.curve2d_idx) {
 Some(c) => c, None => return true,
 };
 let range2d = brep.geom.curve2d_range.get(pc.curve2d_idx).and_then(|r| *r).unwrap_or(range);
 let sample_ts = [0.0, 0.5, 1.0];
 for &t_frac in &sample_ts {
 let t3 = range[0] + t_frac * (range[1] - range[0]);
 let p3d = curve.point_at(t3);
 let t2 = range2d[0] + t_frac * (range2d[1] - range2d[0]);
 let uv = curve2d.point_at(t2);
 let p_surf = surface.point_at(uv.x, uv.y);
 if (p3d - p_surf).length() > tolerance {
 return false;
 }
 }
 true
}


//  € € SameParameter diagnosis  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// A single edge whose 3D curve endpoints deviate from the vertex positions.
///
/// Analogous to the diagnostic output of `BRepCheck_Edge::SameParameter()` in OCCT.
#[derive(Debug, Clone)]
pub struct SuspectEdge {
 /// Index of the edge in `BRep.edges`.
 pub edge_idx: usize,
 /// Distance from `curve.point_at(t1)` to the start vertex position.
 pub start_gap: f64,
 /// Distance from `curve.point_at(t2)` to the end vertex position.
 pub end_gap: f64,
}

/// Report from [`diagnose_same_parameter`].
#[derive(Debug, Clone, Default)]
pub struct SameParameterDiagnosis {
 /// Edges whose curve endpoints deviate from vertex positions beyond `tolerance`.
 pub suspect_edges: Vec<SuspectEdge>,
}

impl SameParameterDiagnosis {
 /// Returns `true` if no SameParameter violations were found.
 pub fn is_clean(&self) -> bool {
 self.suspect_edges.is_empty()
 }
}

/// A single edge whose PCurve ranges deviate from the 3D edge range.
///
/// Analogous to the diagnostic output of `BRepCheck_Edge::SameRange()` in OCCT.
#[derive(Debug, Clone)]
pub struct SuspectSameRangeEdge {
 /// Index of the edge in `BRep.edges`.
 pub edge_idx: usize,
 /// Number of attached PCurves that do not match the 3D edge range.
 pub mismatched_pcurves: usize,
 /// Maximum endpoint mismatch magnitude among attached PCurves.
 pub max_delta: f64,
}

/// Report from [`diagnose_same_range`].
#[derive(Debug, Clone, Default)]
pub struct SameRangeDiagnosis {
 /// Edges whose PCurve ranges deviate from their 3D edge range beyond tolerance.
 pub suspect_edges: Vec<SuspectSameRangeEdge>,
}

impl SameRangeDiagnosis {
 /// Returns `true` if no SameRange violations were found.
 pub fn is_clean(&self) -> bool {
 self.suspect_edges.is_empty()
 }
}

/// A single edge-PCurve pair whose UV-evaluated surface endpoints do not match
/// the edge's 3D curve endpoints.
#[derive(Debug, Clone)]
pub struct SuspectFaceSurfaceEdge {
 pub edge_idx: usize,
 pub pcurve_pos: usize,
 pub surface_idx: usize,
 pub start_gap: f64,
 pub end_gap: f64,
 pub max_gap: f64,
}

/// Report from [`diagnose_face_surface_consistency`].
#[derive(Debug, Clone, Default)]
pub struct FaceSurfaceConsistencyDiagnosis {
 pub suspect_edges: Vec<SuspectFaceSurfaceEdge>,
}

impl FaceSurfaceConsistencyDiagnosis {
 pub fn is_clean(&self) -> bool {
 self.suspect_edges.is_empty()
 }
}

/// Diagnose face-on-surface consistency via PCurves.
///
/// For each edge with a 3D curve range and attached PCurves, evaluates:
/// - 3D endpoints `C3(t1), C3(t2)`
/// - UV endpoints `C2(t1), C2(t2)`
/// - Surface points `S(u1,v1), S(u2,v2)`
///
/// If the surface points deviate from the 3D edge endpoints beyond `tolerance`,
/// the edge-PCurve pair is reported as inconsistent.
pub fn diagnose_face_surface_consistency(
 brep: &BRep,
 tolerance: f64,
) -> FaceSurfaceConsistencyDiagnosis {
 use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};

 let mut suspect_edges = Vec::new();

 for edge_idx in 0..brep.edges.len() {
 let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) else {
 continue;
 };
 let Some(curve3) = brep.geom.curves.get(curve_idx) else {
 continue;
 };
 let Some(range3) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else {
 continue;
 };
 let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
 continue;
 };

 for (pcurve_pos, pc) in pcurves.iter().enumerate() {
 let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else {
 continue;
 };
 let Some(surface) = brep.geom.surfaces.get(pc.surface_idx) else {
 continue;
 };

 let range2 = brep
 .geom
 .curve2d_range
 .get(pc.curve2d_idx)
 .and_then(|r| *r)
 .unwrap_or(range3);

 let p3_start = curve3.point_at(range3[0]);
 let p3_end = curve3.point_at(range3[1]);

 let uv_start = curve2d.point_at(range2[0]);
 let uv_end = curve2d.point_at(range2[1]);

 let ps_start = surface.point_at(uv_start.x, uv_start.y);
 let ps_end = surface.point_at(uv_end.x, uv_end.y);

 let start_gap = (ps_start - p3_start).length();
 let end_gap = (ps_end - p3_end).length();
 let max_gap = start_gap.max(end_gap);

 if max_gap > tolerance {
 suspect_edges.push(SuspectFaceSurfaceEdge {
 edge_idx,
 pcurve_pos,
 surface_idx: pc.surface_idx,
 start_gap,
 end_gap,
 max_gap,
 });
 }
 }
 }

 FaceSurfaceConsistencyDiagnosis { suspect_edges }
}

/// Per-wire diagnostics for gap and self-intersection analysis.
#[derive(Debug, Clone, Default)]
pub struct WireIssueReport {
 pub solid: usize,
 pub shell: usize,
 pub face: usize,
 /// 0 = outer wire, 1..N = inner wire index + 1.
 pub wire_idx: usize,
 pub edge_count: usize,
 pub open_gaps: usize,
 pub topological_self_intersections: usize,
 pub geometric_self_intersections: usize,
}

/// Aggregated report from [`analyze_wire_issues`].
#[derive(Debug, Clone, Default)]
pub struct WireAnalysisReport {
 pub wires: Vec<WireIssueReport>,
 pub total_open_gaps: usize,
 pub total_topological_self_intersections: usize,
 pub total_geometric_self_intersections: usize,
}

impl WireAnalysisReport {
 /// Returns true when no gap or self-intersection issue was found.
 pub fn is_clean(&self) -> bool {
 self.total_open_gaps == 0
 && self.total_topological_self_intersections == 0
 && self.total_geometric_self_intersections == 0
 }
}

/// Analyze all face wires for gap and self-intersection issues.
///
/// This is a structured counterpart to checker issues C1/C7/C8 and is useful
/// for import diagnostics and healing reports.
/// Analyze wire issues in a BRep by examining topology directly.
/// Works on topods::BRep — no old BRep bridge needed.
pub fn analyze_wire_issues(brep: &topods::BRep, tolerance: f64) -> WireAnalysisReport {
 let tshapes = &brep.tshapes;

 let mut report = WireAnalysisReport::default();

 // Build vertex point lookup: tshape index → point
 let v_points: Vec<DVec3> = tshapes.iter().filter_map(|ts| {
 if let topods::TShape::Vertex(vd) = &**ts {
 Some(vd.point)
 } else {
 None
 }
 }).collect();

 let mut si = 0usize;
 for ts in tshapes {
 let topods::TShape::Solid(sd) = &**ts else { continue };

 let mut shi = 0usize;
 for shell_sr in &sd.shells {
 let topods::TShape::Shell(shd) = &*tshapes[shell_sr.index] else { continue };

 let mut fi = 0usize;
 for face_sr in &shd.faces {
 let topods::TShape::Face(fd) = &*tshapes[face_sr.index] else { continue };

 // Collect wires: [outer, inner1, inner2, ...]
 let mut all_wires: Vec<(usize, &topods::ShapeRef)> = Vec::new();
 all_wires.push((0, &fd.outer_wire));
 for (wi, inner) in fd.inner_wires.iter().enumerate() {
 all_wires.push((wi + 1, inner));
 }

 for (wire_idx, wire_sr) in all_wires {
 let topods::TShape::Wire(wd) = &*tshapes[wire_sr.index] else { continue };

 let mut wire_verts = Vec::with_capacity(wd.edges.len());
 let mut valid = true;

 for we in &wd.edges {
 if we.index >= tshapes.len() {
 valid = false;
 break;
 }
 let topods::TShape::Edge(ed) = &*tshapes[we.index] else { valid = false; break };
 let (sv, ev) = if we.orientation.is_forward() {
 (ed.first.index, ed.last.index)
 } else {
 (ed.last.index, ed.first.index)
 };
 if sv >= v_points.len() || ev >= v_points.len() {
 valid = false;
 break;
 }
 wire_verts.push((sv, ev));
 }

 if !valid {
 continue;
 }

 let mut open_gaps = 0usize;
 let n = wire_verts.len();
 if n > 1 {
 for i in 0..n {
 let next = (i + 1) % n;
 let end_v = wire_verts[i].1;
 let start_v = wire_verts[next].0;
 if end_v != start_v {
 let end_pt = v_points[end_v];
 let start_pt = v_points[start_v];
 if (end_pt - start_pt).length() > tolerance {
 open_gaps += 1;
 }
 }
 }
 }

 let mut vertex_count: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 for &(sv, ev) in &wire_verts {
 *vertex_count.entry(sv).or_insert(0) += 1;
 *vertex_count.entry(ev).or_insert(0) += 1;
 }
 let topological_self_intersections =
 vertex_count.values().filter(|&&c| c > 2).count();
 let geometric_self_intersections =
 count_geometric_self_intersections(&wire_verts, &v_points);

 if open_gaps > 0
 || topological_self_intersections > 0
 || geometric_self_intersections > 0
 {
 report.wires.push(WireIssueReport {
 solid: si,
 shell: shi,
 face: fi,
 wire_idx,
 edge_count: wire_verts.len(),
 open_gaps,
 topological_self_intersections,
 geometric_self_intersections,
 });
 }

 report.total_open_gaps += open_gaps;
 report.total_topological_self_intersections += topological_self_intersections;
 report.total_geometric_self_intersections += geometric_self_intersections;
 }
 fi += 1;
 }
 shi += 1;
 }
 si += 1;
 }

 report
}

/// Legacy internal helper: old flat BRep for tests that build flat structures.
#[cfg(test)]
pub(crate) fn analyze_wire_issues_flat(brep: &BRep, tolerance: f64) -> WireAnalysisReport {
 // Reuse the topods implementation by converting
 analyze_wire_issues(&brep.to_topods(), tolerance)
}

/// Scan all edges in `brep` for SameParameter violations.
///
/// For each edge that has a known 3D curve and curve range `[t1, t2]`, evaluates
/// the curve at `t1` and `t2` and compares the resulting points to the
/// corresponding vertex positions.  Edges whose endpoint gap exceeds `tolerance`
/// are reported as suspects and their `edge_same_parameter` flag is updated to
/// `false` so that a subsequent [`fix_same_parameter`] pass can repair them.
///
/// Analogous to `BRepLib::SameParameter()` diagnosis step in OCCT.
///
/// Returns a non-mutating diagnosis; call `fix_same_parameter_with_scan` to
/// also repair the flagged edges.
pub fn diagnose_same_parameter(brep: &BRep, tolerance: f64) -> SameParameterDiagnosis {
 use rcad_kernel::geom::CurveEval;

 let mut suspects = Vec::new();
 let n_edges = brep.edges.len();
 let n_verts = brep.vertices.len();

 for edge_idx in 0..n_edges {
 let curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c);
 let Some(ci) = curve_idx else { continue };
 let Some(curve) = brep.geom.curves.get(ci) else { continue };
 let Some(range) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else { continue };

 let edge = &brep.edges[edge_idx];
 if edge.start >= n_verts || edge.end >= n_verts {
 continue;
 }

 let p_start = brep.vertices[edge.start].point;
 let p_end = brep.vertices[edge.end].point;

 let eval_start = curve.point_at(range[0]);
 let eval_end = curve.point_at(range[1]);

 let start_gap = (eval_start - p_start).length();
 let end_gap = (eval_end - p_end).length();

 if start_gap > tolerance || end_gap > tolerance {
 suspects.push(SuspectEdge {
 edge_idx,
 start_gap,
 end_gap,
 });
 }
 }

 SameParameterDiagnosis { suspect_edges: suspects }
}

/// Scan all edges in `brep` for SameRange violations.
///
/// For each edge with known 3D range and attached PCurves, checks whether each
/// referenced `curve2d_range` matches the edge's `[t1, t2]` within `tolerance`.
/// Missing 2D ranges are also treated as violations.
pub fn diagnose_same_range(brep: &BRep, tolerance: f64) -> SameRangeDiagnosis {
 let mut suspects = Vec::new();
 let n_edges = brep.edges.len();

 for edge_idx in 0..n_edges {
 let Some(range3d) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else {
 continue;
 };
 let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
 continue;
 };
 if pcurves.is_empty() {
 continue;
 }

 let mut mismatched_pcurves = 0usize;
 let mut max_delta = 0.0f64;

 for pc in pcurves {
 let Some(range2d) = brep.geom.curve2d_range.get(pc.curve2d_idx).and_then(|r| *r) else {
 mismatched_pcurves += 1;
 max_delta = max_delta.max(tolerance);
 continue;
 };

 let d0 = (range2d[0] - range3d[0]).abs();
 let d1 = (range2d[1] - range3d[1]).abs();
 let d = d0.max(d1);
 if d > tolerance {
 mismatched_pcurves += 1;
 max_delta = max_delta.max(d);
 }
 }

 let flagged_false = brep
 .geom
 .edge_same_range
 .get(edge_idx)
 .is_some_and(|v| !*v);

 if mismatched_pcurves > 0 || flagged_false {
 suspects.push(SuspectSameRangeEdge {
 edge_idx,
 mismatched_pcurves,
 max_delta,
 });
 }
 }

 SameRangeDiagnosis { suspect_edges: suspects }
}

//  € € Shell topology analysis  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Topology analysis report for a BRep's shell structure.
///
/// Analogous to `ShapeAnalysis_Shell` in OCCT's TKShHealing.
#[derive(Debug, Clone, Default)]
pub struct ShellTopologyReport {
 /// `true` if every edge is referenced by exactly 2 faces (no free edges).
 pub is_closed: bool,
 /// `true` if no edge is referenced by more than 2 faces.
 pub is_manifold: bool,
 /// Number of edges referenced by exactly 1 face (free / open edges).
 pub open_edge_count: usize,
 /// Number of edges referenced by 3 or more faces (non-manifold edges).
 pub non_manifold_edge_count: usize,
 /// Number of vertices not referenced by any edge in the BRep.
 pub isolated_vertex_count: usize,
 /// Total edge count in the BRep.
 pub total_edges: usize,
 /// Total face count across all solids/shells.
 pub total_faces: usize,
}

/// Analyze the shell topology of `brep`.
///
/// Counts free edges, non-manifold edges, and isolated vertices to determine
/// whether the BRep represents a closed manifold solid.
///
/// Analogous to `ShapeAnalysis_Shell::LoadUnorientedEdges()` + checks in OCCT.
pub fn analyze_shell_topology(brep: &BRep) -> ShellTopologyReport {
 let total_edges = brep.edges.len();

 // Count how many faces each edge is referenced by.
 let mut edge_face_count: Vec<usize> = vec![0; total_edges];
 let mut total_faces = 0usize;

 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 total_faces += 1;
 for we in &face.outer_wire.edges {
 if we.idx < total_edges {
 edge_face_count[we.idx] += 1;
 }
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 if we.idx < total_edges {
 edge_face_count[we.idx] += 1;
 }
 }
 }
 }
 }
 }

 let open_edge_count = edge_face_count.iter().filter(|&&c| c == 1).count();
 let non_manifold_edge_count = edge_face_count.iter().filter(|&&c| c > 2).count();

 // Count isolated vertices: those not referenced by any edge as start or end.
 let mut vertex_used = vec![false; brep.vertices.len()];
 for edge in &brep.edges {
 if edge.start < vertex_used.len() {
 vertex_used[edge.start] = true;
 }
 if edge.end < vertex_used.len() {
 vertex_used[edge.end] = true;
 }
 }
 let isolated_vertex_count = vertex_used.iter().filter(|&&used| !used).count();

 ShellTopologyReport {
 is_closed: open_edge_count == 0,
 is_manifold: non_manifold_edge_count == 0,
 open_edge_count,
 non_manifold_edge_count,
 isolated_vertex_count,
 total_edges,
 total_faces,
 }
}

//  € € Euler characteristic analysis  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Euler characteristic and topological genus for a single solid.
///
/// For a closed orientable 2-manifold of genus *g*: = V  ?E + F = 2  ?2g.
///
/// | Shape | | genus |
/// | sphere  | 2  | 0 |
/// | torus | 0  | 1 |
/// | 2-torus | -2 | 2 |
#[derive(Debug, Clone)]
pub struct EulerAnalysis {
 pub solid_idx: usize,
 /// Unique vertices referenced by this solid's face boundaries.
 pub vertices: usize,
 /// Unique edges referenced by this solid's face boundaries.
 pub edges: usize,
 /// Total faces across all shells of this solid.
 pub faces: usize,
 /// Euler characteristic: V  ?E + F.
 pub euler_number: i64,
 /// `true` if every edge of this solid is referenced by exactly 2 faces
 /// (no free boundary edges  ?closed shell).
 pub is_closed: bool,
 /// Topological genus, computed as `(2  ?euler_number) / 2`.
 /// `None` when `!is_closed` or when the result is not a non-negative integer.
 pub genus: Option<i64>,
}

/// Compute per-solid Euler analysis for every solid in `brep`.
///
/// Analogous to the topological analysis portion of `BRepCheck_Analyzer`.
pub fn euler_analysis(brep: &BRep) -> Vec<EulerAnalysis> {
 let mut results = Vec::with_capacity(brep.solids.len());

 for (si, solid) in brep.solids.iter().enumerate() {
 let mut unique_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut face_count = 0usize;

 for shell in &solid.shells {
 for face in &shell.faces {
 face_count += 1;
 for we in &face.outer_wire.edges {
 unique_edges.insert(we.idx);
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 unique_edges.insert(we.idx);
 }
 }
 }
 }

 // Collect unique vertex indices from the unique edges.
 let mut unique_verts: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for &ei in &unique_edges {
 if let Some(e) = brep.edges.get(ei) {
 unique_verts.insert(e.start);
 unique_verts.insert(e.end);
 }
 }

 let v = unique_verts.len();
 let e = unique_edges.len();
 let f = face_count;
 let euler_number = v as i64 - e as i64 + f as i64;

 // Determine if the solid is closed (every edge shared by exactly 2 faces).
 let mut edge_face_count: std::collections::HashMap<usize, usize> =
 std::collections::HashMap::new();
 for shell in &solid.shells {
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 *edge_face_count.entry(we.idx).or_insert(0) += 1;
 }
 for wire in &face.inner_wires {
 for we in &wire.edges {
 *edge_face_count.entry(we.idx).or_insert(0) += 1;
 }
 }
 }
 }
 let is_closed = edge_face_count.values().all(|&c| c == 2);

 // Genus only makes sense on a closed manifold.
 let genus = if is_closed {
 let g = (2 - euler_number) / 2;
 // Valid genus: (2  ? ) must be even and non-negative.
 if (2 - euler_number) % 2 == 0 && g >= 0 {
 Some(g)
 } else {
 None
 }
 } else {
 None
 };

 results.push(EulerAnalysis {
 solid_idx: si,
 vertices: v,
 edges: e,
 faces: f,
 euler_number,
 is_closed,
 genus,
 });
 }

 results
}

//  € € Orientation consistency analysis  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// A face whose stored normal appears to point inward rather than outward.
///
/// Determined geometrically: the dot product of `face.normal` with the vector
/// from the solid's vertex centroid to the face's centroid is negative.
#[derive(Debug, Clone)]
pub struct OrientationIssue {
 /// Solid index.
 pub solid_idx: usize,
 /// Flat face index (counting across all solids/shells in traversal order).
 pub face_idx: usize,
 /// Dot product of `face.normal` with the outward radial direction.
 /// Negative values indicate an inward-pointing normal.
 pub dot_product: f64,
}

/// Report from [`check_orientation_consistency`].
#[derive(Debug, Clone, Default)]
pub struct OrientationReport {
 /// `true` iff all face normals appear to point outward from the solid interior.
 pub is_consistent: bool,
 /// Faces whose stored normals appear to point inward.
 pub issues: Vec<OrientationIssue>,
 /// Number of faces with outward-pointing normals.
 pub consistent_face_count: usize,
 /// Number of faces with inward-pointing normals.
 pub inconsistent_face_count: usize,
}

/// Check that every face's stored normal points **outward** from the solid interior.
///
/// Uses a geometric heuristic: the solid's interior centroid is approximated as
/// the average of all its vertices.  For each face the face centroid is computed
/// from the outer-wire corner vertices.  A face normal is considered outward when
/// `face.normal (face_centroid  ?solid_centroid) > 0`.
///
/// This correctly handles all primitives (box, sphere, cylinder, cone, torus)
/// whose normals are computed analytically during construction.  For shapes
/// produced by Boolean operations or user-constructed BReps, the heuristic may
/// give false positives on highly non-convex solids  ?treat the report as an
/// advisory rather than a hard constraint.
///
/// Analogous to `BRepCheck_Shell::Orientation()` in OCCT.
pub fn check_orientation_consistency(brep: &BRep) -> OrientationReport {
 use glam::DVec3;

 let mut issues = Vec::new();
 let mut consistent_face_count = 0usize;
 let mut inconsistent_face_count = 0usize;
 let mut flat_face_idx = 0usize;

 for (si, solid) in brep.solids.iter().enumerate() {
 // Compute the solid's interior centroid from all vertices that appear
 // in ANY of this solid's face wires.
 let mut solid_verts = std::collections::HashSet::new();
 for shell in &solid.shells {
 for face in &shell.faces {
 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 solid_verts.insert(brep.edges[we.idx].start);
 solid_verts.insert(brep.edges[we.idx].end);
 }
 }
 }
 }
 if solid_verts.is_empty() {
 // No geometry  ?skip.
 for shell in &solid.shells {
 flat_face_idx += shell.faces.len();
 }
 continue;
 }
 let solid_centroid: DVec3 = {
 let sum: DVec3 = solid_verts
 .iter()
 .filter(|&&vi| vi < brep.vertices.len())
 .map(|&vi| brep.vertices[vi].point)
 .sum();
 sum / solid_verts.len() as f64
 };

 for shell in &solid.shells {
 for face in &shell.faces {
 // Face centroid: average of first vertices of each outer-wire edge.
 let mut face_centroid = DVec3::ZERO;
 let mut n = 0usize;
 for we in &face.outer_wire.edges {
 if we.idx < brep.edges.len() {
 let vi = if we.forward {
 brep.edges[we.idx].start
 } else {
 brep.edges[we.idx].end
 };
 if vi < brep.vertices.len() {
 face_centroid += brep.vertices[vi].point;
 n += 1;
 }
 }
 }
 if n == 0 {
 flat_face_idx += 1;
 continue;
 }
 face_centroid /= n as f64;

 let outward = face_centroid - solid_centroid;
 let dot = face.normal.dot(outward);
 if dot >= 0.0 {
 consistent_face_count += 1;
 } else {
 inconsistent_face_count += 1;
 issues.push(OrientationIssue {
 solid_idx: si,
 face_idx: flat_face_idx,
 dot_product: dot,
 });
 }
 flat_face_idx += 1;
 }
 }
 }

 OrientationReport {
 is_consistent: issues.is_empty(),
 issues,
 consistent_face_count,
 inconsistent_face_count,
 }
}

//  € € Comprehensive richer validity analysis  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Aggregated validity report combining all available checks.
///
/// This is the RCAD equivalent of OCCT's `BRepCheck_Analyzer` + `ShapeAnalysis`
/// combined output, giving a single entry-point for full BRep validation.
#[derive(Debug, Clone)]
pub struct RicherValidityReport {
 /// Basic structural check result.
 pub check_result: CheckResult,
 /// Shell topology (free edges, non-manifold edges, isolation).
 pub shell_topology: ShellTopologyReport,
 /// Per-solid Euler characteristic and genus.
 pub euler: Vec<EulerAnalysis>,
 /// Face wire orientation consistency across shared edges.
 pub orientation: OrientationReport,
 /// `true` iff BRep passes all structural checks and orientation is consistent.
 pub is_fully_valid: bool,
}

impl RicherValidityReport {
 /// A short human-readable summary of this report.
 pub fn summary(&self) -> String {
 let issues = self.check_result.issues.len();
 let euler_issues = self.euler.iter().filter(|e| e.genus.is_none_or(|g| g < 0)).count();
 let orient_issues = self.orientation.inconsistent_face_count;
 if self.is_fully_valid {
 format!(
 "valid: {} solids, closed={}, manifold={}",
 self.euler.len(),
 self.shell_topology.is_closed,
 self.shell_topology.is_manifold,
 )
 } else {
 format!(
 "INVALID: {} structural issue(s), {} orientation inconsistency/ies, {} genus anomaly/ies",
 issues, orient_issues, euler_issues,
 )
 }
 }
}

/// Run all available validity checks on `brep` and return a consolidated report.
///
/// This is the preferred entry point for comprehensive BRep validation:
/// it combines the basic `brep_check_analyze()`, shell topology, Euler analysis, and
/// orientation consistency into a single call.
pub fn richer_validity_analysis(brep: &BRep) -> RicherValidityReport {
 let check_result = brep_check_analyze(brep);
 let shell_topology = analyze_shell_topology(brep);
 let euler = euler_analysis(brep);
 let orientation = check_orientation_consistency(brep);

 let is_fully_valid = check_result.is_valid() && orientation.is_consistent;

 RicherValidityReport {
 check_result,
 shell_topology,
 euler,
 orientation,
 is_fully_valid,
 }
}

//  € € Surface UV Analysis (ShapeAnalysis_Surface equivalent)  € € € € € € € € € € € € € € € € € € € € € € €

/// Report from surface UV domain analysis.
///
/// Analogous to OCCT's `ShapeAnalysis_Surface` which checks UV bounds,
/// periodicity, and parameter space consistency.
#[derive(Debug, Clone, Default)]
pub struct SurfaceAnalysisReport {
 /// Number of faces analyzed.
 pub faces_analyzed: usize,
 /// Faces whose PCurve parameter ranges violate expected surface bounds.
 pub faces_with_uv_bounds_violation: Vec<UvBoundsViolation>,
 /// Total number of issues detected.
 pub total_issues: usize,
}

impl SurfaceAnalysisReport {
 pub fn is_clean(&self) -> bool {
 self.total_issues == 0
 }

 pub fn summary(&self) -> String {
 if self.is_clean() {
 format!("{} faces analyzed, no UV issues", self.faces_analyzed)
 } else {
 format!(
 "{} faces analyzed, {} UV bounds violations",
 self.faces_analyzed,
 self.faces_with_uv_bounds_violation.len()
 )
 }
 }
}

/// UV bounds violation for a single face.
#[derive(Debug, Clone)]
pub struct UvBoundsViolation {
 pub solid: usize,
 pub shell: usize,
 pub face: usize,
 /// Surface type (Plane, Cylinder, etc.)
 pub surface_type: String,
 /// Expected UV bounds [u_min, u_max, v_min, v_max] for the surface.
 pub expected_bounds: [f64; 4],
 /// Actual UV bounds of the face's PCurves.
 pub actual_bounds: [f64; 4],
 /// Maximum violation distance.
 pub violation: f64,
}
include!("e1.rs");
include!("tests_inc.rs");
