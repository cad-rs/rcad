//! Surface tessellation (adaptive chord error in UV / world space).
//!
//! **Phase C note:** several internal checks use [`TOLERANCE_MESH_LEGACY`] and related constants as
//! **dimensionless UV / angle slacks** (parameter-domain noise, multi-turn heuristics), not as a
//! substitute for world-space BRep pairing. For mesh esh segment chaining tied to operand
//! topology, use [`crate::tolerance::tessellation_merge_linear_from_brep`] /
//! [`crate::tolerance::tessellation_merge_linear_from_two_breps`] and
//! [`crate::section::intersect_triangle_soups_adaptive`].
//!

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::topods;
use rcad_kernel::topology;
use rcad_kernel::geom::{any_perpendicular, Curve3, CurveEval, Surface3, SurfaceEval};
use std::collections::HashMap;

use crate::offset::project_point_to_surface_uv;

/// High-quality tessellation of a parametric surface.
#[derive(Debug, Clone)]
pub struct SurfaceMesh {
 /// Tessellated node positions (world coordinates).
 pub nodes: Vec<DVec3>,
 /// Triangle corner indices (three node indices per triangle).
 pub triangles: Vec<[usize; 3]>,
 /// Per-node shading normals.
 pub normals: Vec<DVec3>,
 /// When `true` the mesh data is out of date with respect to the source
 /// geometry and must be recomputed before use.
 ///
 /// `triangulate_surface` always returns a clean mesh (`dirty = false`).
 /// Callers that cache a `SurfaceMesh` should call [`SurfaceMesh::invalidate`]
 /// whenever the source geometry changes.
 pub dirty: bool,
}

impl SurfaceMesh {
 /// Mark this mesh as stale.  The next render or query should recompute it.
 pub fn invalidate(&mut self) {
 self.dirty = true;
 }

 /// Returns `true` if the mesh data is up-to-date with the source geometry.
 pub fn is_clean(&self) -> bool {
 !self.dirty
 }
}

/// Parameters controlling surface tessellation.
#[derive(Debug, Clone)]
pub struct TessellationParams {
 // --- Basic controls ---
 /// Maximum chordal deviation (midpoint of a triangle edge to the true surface).
 /// Smaller values yield finer meshes; typical range `0.001`= 0.1`.
 pub chord_tolerance: f64,
 /// Maximum angle error (radians) between adjacent triangle normals before splitting.
 pub angle_tolerance: f64,
 /// Minimum UV step; stops runaway refinement.
 pub min_step: f64,
 /// Maximum UV step.
 pub max_step: f64,

 // --- Size limits ---
 /// Minimum triangle size in world units; smaller triangles are not split further (`0.0` = off).
 pub min_triangle_size: f64,
 /// Maximum triangle size in world units; larger patches are split (`f64::MAX` = off).
 pub max_triangle_size: f64,

 // --- Quality ---
 /// Enable adaptive refinement (curvature-driven density). Default `true`.
 pub adaptive_refinement: bool,
 /// Prefer finer mesh in high-curvature regions. Default `true`.
 pub curvature_sensitive: bool,
 /// Soft quality cap on triangle aspect ratio (diagnostics). Default `20.0`.
 pub max_aspect_ratio: f64,

 // --- Boundary / seams ---
 /// Preserve boundary nodes when welding samples. Default `true`.
 pub boundary_preservation: bool,
 /// Preserve seam edges (special handling when welding). Default `true`.
 pub seam_preservation: bool,

 // --- Performance ---
 /// Maximum recursion depth (safety cap). Default `8`.
 pub max_depth: usize,
 /// If `true`, callers may tessellate faces in parallel (not used everywhere). Default `false`.
 pub parallel: bool,
}

impl Default for TessellationParams {
 fn default() -> Self {
 Self {
 chord_tolerance: 0.01,
 angle_tolerance: 0.1,  // ~5.7 degrees
 min_step: TOLERANCE_RETRY_LADDER_COARSE,
 max_step: 0.5,
 min_triangle_size: 0.0,
 max_triangle_size: f64::MAX,
 adaptive_refinement: true,
 curvature_sensitive: true,
 max_aspect_ratio: 20.0,
 boundary_preservation: true,
 seam_preservation: true,
 max_depth: 8,
 parallel: false,
 }
 }
}

impl TessellationParams {
 /// Fast preview preset (interactive viewing, favors speed).
 pub fn preview() -> Self {
 Self {
 chord_tolerance: 0.1,
 angle_tolerance: 0.3,  // ~17 degrees
 min_step: TOLERANCE_ADAPTIVE_MAX,
 max_step: 1.0,
 min_triangle_size: 0.01,
 max_triangle_size: f64::MAX,
 adaptive_refinement: false,
 curvature_sensitive: false,
 max_aspect_ratio: 30.0,
 boundary_preservation: true,
 seam_preservation: false,
 max_depth: 4,
 parallel: true,
 }
 }

 /// Balanced default preset for general use.
 pub fn standard() -> Self {
 Self {
 chord_tolerance: 0.01,
 angle_tolerance: 0.1,  // ~5.7 degrees
 min_step: TOLERANCE_RETRY_LADDER_COARSE,
 max_step: 0.5,
 min_triangle_size: 0.0,
 max_triangle_size: f64::MAX,
 adaptive_refinement: true,
 curvature_sensitive: true,
 max_aspect_ratio: 20.0,
 boundary_preservation: true,
 seam_preservation: true,
 max_depth: 8,
 parallel: false,
 }
 }

 /// High-quality preset for rendering and visualization.
 pub fn high_quality() -> Self {
 Self {
 chord_tolerance: 0.001,
 angle_tolerance: 0.05,  // ~2.9 degrees
 min_step: TOLERANCE_RETRY_LADDER_MID,
 max_step: 0.2,
 min_triangle_size: 0.0,
 max_triangle_size: f64::MAX,
 adaptive_refinement: true,
 curvature_sensitive: true,
 max_aspect_ratio: 10.0,
 boundary_preservation: true,
 seam_preservation: true,
 max_depth: 12,
 parallel: false,
 }
 }

 /// Export-oriented preset (STL/OBJ, reasonable file size).
 pub fn export() -> Self {
 Self {
 chord_tolerance: 0.005,
 angle_tolerance: 0.08,  // ~4.6 degrees
 min_step: TOLERANCE_RETRY_LADDER_COARSE,
 max_step: 0.3,
 min_triangle_size: 0.0,
 max_triangle_size: f64::MAX,
 adaptive_refinement: true,
 curvature_sensitive: true,
 max_aspect_ratio: 15.0,
 boundary_preservation: true,
 seam_preservation: true,
 max_depth: 10,
 parallel: false,
 }
 }

 /// Analysis-oriented preset (FEA/CFD-style mesh quality).
 pub fn analysis() -> Self {
 Self {
 chord_tolerance: 0.0005,
 angle_tolerance: 0.03,  // ~1.7 degrees
 min_step: TOLERANCE_MESH_LEGACY,
 max_step: 0.1,
 min_triangle_size: 0.0,
 max_triangle_size: f64::MAX,
 adaptive_refinement: true,
 curvature_sensitive: true,
 max_aspect_ratio: 5.0,
 boundary_preservation: true,
 seam_preservation: true,
 max_depth: 15,
 parallel: false,
 }
 }

 /// Scale tolerances from a nominal triangle count heuristic; returns a new `Self`.
 pub fn with_target_triangle_count(&self, target_count: usize) -> Self {
 let factor = (target_count as f64 / 1000.0).powf(1.0 / 3.0).max(0.1).min(10.0);
 Self {
 chord_tolerance: self.chord_tolerance * factor,
 angle_tolerance: self.angle_tolerance * factor,
 min_step: self.min_step,
 max_step: self.max_step / factor,
 min_triangle_size: self.min_triangle_size,
 max_triangle_size: self.max_triangle_size,
 adaptive_refinement: self.adaptive_refinement,
 curvature_sensitive: self.curvature_sensitive,
 max_aspect_ratio: self.max_aspect_ratio,
 boundary_preservation: self.boundary_preservation,
 seam_preservation: self.seam_preservation,
 max_depth: self.max_depth,
 parallel: self.parallel,
 }
 }
}

/// Adaptive chord-error tessellation of a parametric surface.
///
/// Algorithm (UV domain):
/// 1. Start from a uniform quad grid over `[u_min, u_max]  ?[v_min, v_max]`.
/// 2. For each quad, measure chord error (linear patch vs true surface).
/// 3. Recursively split quads that exceed `params.chord_tolerance` (and related checks).
/// 4. Emit two triangles per leaf quad.
///
/// # Arguments
/// - `surface`: surface to tessellate
/// - `u_range`: `[u_min, u_max]`
/// - `v_range`: `[v_min, v_max]`
/// - `params`: tessellation controls
pub fn triangulate_surface(
 surface: &Surface3,
 u_range: [f64; 2],
 v_range: [f64; 2],
 params: &TessellationParams,
) -> SurfaceMesh {
 let mut nodes: Vec<DVec3> = Vec::new();
 let mut normals: Vec<DVec3> = Vec::new();
 let mut triangles: Vec<[usize; 3]> = Vec::new();

 // Initial UV grid resolution (at least 2 ? quads)
 let initial_steps = 4usize;
 let [u0, u1] = u_range;
 let [v0, v1] = v_range;
 let du = (u1 - u0) / initial_steps as f64;
 let dv = (v1 - v0) / initial_steps as f64;

 // Adaptive refinement for each seed quad
 for i in 0..initial_steps {
 for j in 0..initial_steps {
 let ua = u0 + i as f64 * du;
 let ub = ua + du;
 let va = v0 + j as f64 * dv;
 let vb = va + dv;

 subdivide_quad(
 surface,
 [ua, ub],
 [va, vb],
 params,
 0,
 &mut nodes,
 &mut normals,
 &mut triangles,
 );
 }
 }

 weld_surface_mesh_nodes(SurfaceMesh {
 nodes,
 triangles,
 normals,
 dirty: false,
 })
}

fn weld_surface_mesh_nodes(mesh: SurfaceMesh) -> SurfaceMesh {
 const WELD_TOLERANCE: f64 = TOLERANCE_COORD_SUB;

 let mut remap = vec![0usize; mesh.nodes.len()];
 let mut welded_nodes: Vec<DVec3> = Vec::new();
 let mut welded_normals: Vec<DVec3> = Vec::new();
 let mut normal_counts = Vec::new();
 let mut buckets: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
 let scale = 1.0 / WELD_TOLERANCE;

 for (index, point) in mesh.nodes.iter().enumerate() {
 let key = [
 (point.x * scale).round() as i64,
 (point.y * scale).round() as i64,
 (point.z * scale).round() as i64,
 ];

 let mut matched = None;
 if let Some(candidates) = buckets.get(&key) {
 for &candidate in candidates {
 if (welded_nodes[candidate] - *point).length_squared() <= WELD_TOLERANCE * WELD_TOLERANCE {
 matched = Some(candidate);
 break;
 }
 }
 }

 let target = if let Some(existing) = matched {
 existing
 } else {
 let new_index = welded_nodes.len();
 welded_nodes.push(*point);
 welded_normals.push(DVec3::ZERO);
 normal_counts.push(0usize);
 buckets.entry(key).or_default().push(new_index);
 new_index
 };

 remap[index] = target;
 if let Some(normal) = mesh.normals.get(index) {
 welded_normals[target] += *normal;
 normal_counts[target] += 1;
 }
 }

 let welded_triangles: Vec<[usize; 3]> = mesh
 .triangles
 .iter()
 .filter_map(|&[a, b, c]| {
 let ra = remap[a];
 let rb = remap[b];
 let rc = remap[c];
 if ra == rb || rb == rc || rc == ra {
 None
 } else {
 Some([ra, rb, rc])
 }
 })
 .collect();

 let welded_normals: Vec<DVec3> = welded_normals
 .into_iter()
 .zip(normal_counts)
 .map(|(normal, count)| {
 if count == 0 {
 DVec3::ZERO
 } else {
 normal.normalize_or_zero()
 }
 })
 .collect();

 SurfaceMesh {
 nodes: welded_nodes,
 triangles: welded_triangles,
 normals: welded_normals,
 dirty: mesh.dirty,
 }
}

/// Recursively refine one UV-space quad.
fn subdivide_quad(
 surface: &Surface3,
 u_range: [f64; 2],
 v_range: [f64; 2],
 params: &TessellationParams,
 depth: usize,
 nodes: &mut Vec<DVec3>,
 normals: &mut Vec<DVec3>,
 triangles: &mut Vec<[usize; 3]>,
) {
 let [u0, u1] = u_range;
 let [v0, v1] = v_range;

 // Corner evaluations
 let p00 = surface.point_at(u0, v0);
 let p10 = surface.point_at(u1, v0);
 let p01 = surface.point_at(u0, v1);
 let p11 = surface.point_at(u1, v1);

 let um = (u0 + u1) * 0.5;
 let vm = (v0 + v1) * 0.5;

 // Decide whether to subdivide further
 let should_subdivide = if depth < params.max_depth {
 let step_u = u1 - u0;
 let step_v = v1 - v0;

 // Stop if UV steps are already at the minimum
 if step_u < params.min_step * 2.0 && step_v < params.min_step * 2.0 {
 false
 } else {
 // Chord error: triangle centroids vs surface
 let chord_exceeded = check_chord_tolerance(surface, p00, p10, p11, p01, um, vm, params.chord_tolerance);

 // Normal variation across the quad
 let angle_exceeded = depth < params.max_depth / 2 && check_angle_tolerance(surface, u0, u1, v0, v1, params.angle_tolerance);

 // World-space size cap
 let size_exceeded = if params.max_triangle_size < f64::MAX {
 let diag = (p11 - p00).length();
 diag > params.max_triangle_size
 } else {
 false
 };

 // Combine criteria when adaptive refinement is on
 if params.adaptive_refinement {
 chord_exceeded || angle_exceeded || size_exceeded
 } else {
 size_exceeded
 }
 }
 } else {
 false
 };

 if should_subdivide {
 // Split into four child quads
 subdivide_quad(surface, [u0, um], [v0, vm], params, depth + 1, nodes, normals, triangles);
 subdivide_quad(surface, [um, u1], [v0, vm], params, depth + 1, nodes, normals, triangles);
 subdivide_quad(surface, [u0, um], [vm, v1], params, depth + 1, nodes, normals, triangles);
 subdivide_quad(surface, [um, u1], [vm, v1], params, depth + 1, nodes, normals, triangles);
 } else {
 // Emit two triangles
 let n = nodes.len();

 // Corner normals
 let n00 = surface.normal_at(u0, v0);
 let n10 = surface.normal_at(u1, v0);
 let n01 = surface.normal_at(u0, v1);
 let n11 = surface.normal_at(u1, v1);

 // Drop degenerate corner data
 let valid = [p00, p10, p01, p11].iter().all(|p| p.is_finite());
 if !valid {
 return;
 }

 nodes.extend_from_slice(&[p00, p10, p11, p01]);
 normals.extend_from_slice(&[n00, n10, n11, n01]);

 // Pick the shorter diagonal for better triangle shape
 let d0 = (p11 - p00).length_squared();
 let d1 = (p10 - p01).length_squared();
 if d0 <= d1 {
 // Diagonal p00= 11
 triangles.push([n, n + 1, n + 2]);
 triangles.push([n, n + 2, n + 3]);
 } else {
 // Diagonal p10= 01
 triangles.push([n, n + 1, n + 3]);
 triangles.push([n + 1, n + 2, n + 3]);
 }
 }
}

/// Whether chord error exceeds tolerance (triangle centroids vs surface).
fn check_chord_tolerance(
 surface: &Surface3,
 p00: DVec3, p10: DVec3, p11: DVec3, p01: DVec3,
 um: f64, vm: f64,
 tolerance: f64,
) -> bool {
 // Centroid of triangle (p00, p10, p11)
 let c1 = (p00 + p10 + p11) / 3.0;
 // Centroid of triangle (p00, p11, p01)
 let c2 = (p00 + p11 + p01) / 3.0;

 // True surface point at UV center
 let surf_mid = surface.point_at(um, vm);

 // Bilinear patch midpoint vs surface
 let interp_mid = (p00 + p10 + p11 + p01) / 4.0;
 let chord_err = (surf_mid - interp_mid).length();

 // Also sample chord error at triangle centroids
 let t1_u = (c1 - p00).length() / (p11 - p00).length().max(TOLERANCE_LINEAR_ULTRA_STRICT);
 let _ = t1_u; // reserved for finer UV mapping at centroids
 let chord1 = (surface.point_at(um, vm) - c1).length();
 let chord2 = (surface.point_at(um, vm) - c2).length();

 chord_err > tolerance || chord1 > tolerance || chord2 > tolerance
}

/// Whether normal variation across the quad exceeds `tolerance` (radians).
fn check_angle_tolerance(
 surface: &Surface3,
 u0: f64, u1: f64, v0: f64, v1: f64,
 tolerance: f64,
) -> bool {
 let n00 = surface.normal_at(u0, v0);
 let n11 = surface.normal_at(u1, v1);
 let n10 = surface.normal_at(u1, v0);
 let n01 = surface.normal_at(u0, v1);

 // Angle between normals at adjacent corners
 for (a, b) in [(n00, n10), (n00, n01), (n11, n10), (n11, n01)] {
 let la = a.length();
 let lb = b.length();
 if la < 0.5 || lb < 0.5 {
 continue;
 }
 let cos_a = (a.dot(b) / (la * lb)).clamp(-1.0, 1.0);
 let angle = cos_a.acos();
 if angle > tolerance {
 return true;
 }
 }
 false
}

/// Ear-clipping triangulation for a simple polygon in 3D.
/// Projects to 2D using the given normal, then runs ear-clipping.
pub fn triangulate_polygon(nodes: &[DVec3], normal: DVec3) -> Vec<[usize; 3]> {
 let n = nodes.len();
 if n < 3 {
 return vec![];
 }
 if n == 3 {
 return vec![[0, 1, 2]];
 }
 if n == 4 {
 return vec![[0, 1, 2], [0, 2, 3]];
 }

 let (u_axis, v_axis) = local_basis(normal);
 let pts_2d: Vec<[f64; 2]> = nodes
 .iter()
 .map(|p| [p.dot(u_axis), p.dot(v_axis)])
 .collect();

 ear_clip(&pts_2d)
}

/// Triangulate an outer ring with optional hole rings.
///
/// Returns local indices into `outer + holes` point order.
pub(crate) fn triangulate_polygon_with_holes(
 outer: &[DVec3],
 holes: &[Vec<DVec3>],
 normal: DVec3,
) -> Vec<[usize; 3]> {
 if outer.len() < 3 {
 return Vec::new();
 }
 let (u_axis, v_axis) = local_basis(normal);
 let mut flat: Vec<f64> = Vec::new();
 let mut hole_starts: Vec<usize> = Vec::new();
 let mut vertex_count = 0usize;

 for p in outer {
 flat.push(p.dot(u_axis));
 flat.push(p.dot(v_axis));
 vertex_count += 1;
 }
 for hole in holes {
 if hole.len() < 3 {
 continue;
 }
 hole_starts.push(vertex_count);
 for p in hole {
 flat.push(p.dot(u_axis));
 flat.push(p.dot(v_axis));
 vertex_count += 1;
 }
 }

 let coords: Vec<[f64; 2]> = flat
 .chunks_exact(2)
 .map(|c| [c[0], c[1]])
 .collect();
 let mut indices: Vec<usize> = Vec::new();
 {
 let mut ear = earcut::Earcut::new();
 ear.earcut(coords, &hole_starts, &mut indices);
 }
 let mut tris = Vec::new();
 for tri in indices.chunks_exact(3) {
 tris.push([tri[0], tri[1], tri[2]]);
 }
 tris
}

fn estimate_polygon_normal(points: &[DVec3]) -> Option<DVec3> {
 if points.len() < 3 {
 return None;
 }
 // Newell-style robust polygon normal estimate from boundary points.
 let mut n = DVec3::ZERO;
 for i in 0..points.len() {
 let a = points[i];
 let b = points[(i + 1) % points.len()];
 n.x += (a.y - b.y) * (a.z + b.z);
 n.y += (a.z - b.z) * (a.x + b.x);
 n.z += (a.x - b.x) * (a.y + b.y);
 }
 let len2 = n.length_squared();
 if len2 <= TOLERANCE_METRIC_SQ_NEAR_ZERO {
 None
 } else {
 Some(n / len2.sqrt())
 }
}


fn local_basis(normal: DVec3) -> (DVec3, DVec3) {
 let ref_dir = if normal.x.abs() < 0.9 {
 DVec3::X
 } else {
 DVec3::Y
 };
 let u = normal.cross(ref_dir).normalize();
 let v = normal.cross(u).normalize();
 (u, v)
}

/// Build an ordered 3D polygon from a wire.
///
/// Curved edges are sampled using their analytic 3D curve + edge range.
/// Straight or missing-geometry edges contribute only their end vertex.
fn sample_wire_polygon_points(brep: &rcad_kernel::BRep, wire: &rcad_kernel::topology::Wire) -> Vec<DVec3> {
 use std::collections::HashSet;
 let mut pts: Vec<DVec3> = Vec::new();
 let two_pi = 2.0 * std::f64::consts::PI;
 let mut seen_edge_indices: HashSet<usize> = HashSet::new();

 for we in &wire.edges {
 // Some imported wires repeat seam/periodic edges to force loop closure.
 // Keep one geometric contribution per topological edge to avoid
 // self-overlapping polygon chains that destabilize ear clipping.
 if !seen_edge_indices.insert(we.idx) {
 continue;
 }
 let Some(topods::TShape::Edge(ed)) = brep.tshapes.get(we.idx).map(|ts| ts.as_ref()) else {
 continue;
 };
 let topo_closed = ed.first.index == ed.last.index;

 let start_idx = if we.forward { ed.first.index } else { ed.last.index };
 let end_idx = if we.forward { ed.last.index } else { ed.first.index };

 let p_start = match brep.vertex_point(start_idx) {
 Some(pt) => pt,
 None => continue,
 };
 let p_end = match brep.vertex_point(end_idx) {
 Some(pt) => pt,
 None => continue,
 };
 let edge_tol = ed.tolerance.max(0.0);

 let mut sampled = false;
 if let Some(curve) = &ed.curve
 && !matches!(curve, Curve3::Line(_)) {
 let [r0, r1] = ed.range;

 let mut t0 = r0;
 let mut t1 = r1;
 if matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_)) {
 // Robust unit disambiguation by endpoint fit against topological
 // edge vertices; avoids misclassifying valid unwrapped radians.
 let range_rad = [t0, t1];
 let range_deg = [t0.to_radians(), t1.to_radians()];
 let err_for = |r: [f64; 2]| -> f64 {
 let a = curve.point_at(r[0]);
 let b = curve.point_at(r[1]);
 let e_direct = (a - p_start).length() + (b - p_end).length();
 let e_swapped = (a - p_end).length() + (b - p_start).length();
 e_direct.min(e_swapped)
 };
 let err_rad = err_for(range_rad);
 let err_deg = err_for(range_deg);
 let max_abs = t0.abs().max(t1.abs());
 let span_abs = (t1 - t0).abs();
 let seam_like = (p_start - p_end).length() <= TOLERANCE_ABS;
 let degree_likely_by_magnitude =
 max_abs > two_pi + TOLERANCE_COORD_SUB && max_abs <= 360.0 + TOLERANCE_MESH_LEGACY;
 let tie =
 (err_deg - err_rad).abs() <= TOLERANCE_LINEAR_RELAX_8 * (1.0 + err_rad.abs().max(err_deg.abs()));
 if err_deg + TOLERANCE_COORD_SUB < err_rad
 || (tie && degree_likely_by_magnitude)
 || (seam_like && degree_likely_by_magnitude && span_abs > two_pi * 1.5)
 {
 t0 = range_deg[0];
 t1 = range_deg[1];
 }
 }
 if !we.forward {
 std::mem::swap(&mut t0, &mut t1);
 }
 let near_full_turn = ((t1 - t0).abs() - two_pi).abs() <= TOLERANCE_ADAPTIVE_MAX;
 let choose_arc_delta = |a0: f64, a1: f64, span_hint: f64| -> f64 {
 let mut minor = a1 - a0;
 if minor > std::f64::consts::PI {
 minor -= two_pi;
 } else if minor < -std::f64::consts::PI {
 minor += two_pi;
 }
 let major = if minor >= 0.0 {
 minor - two_pi
 } else {
 minor + two_pi
 };
 if !span_hint.is_finite() || span_hint.abs() <= TOLERANCE_LEN_MIN {
 return minor;
 }
 let score = |cand: f64| -> f64 {
 let mag = (cand.abs() - span_hint.abs()).abs();
 let sign_penalty = if cand.signum() == span_hint.signum() {
 0.0
 } else {
 two_pi
 };
 mag + sign_penalty
 };
 if score(major) < score(minor) {
 major
 } else {
 minor
 }
 };

 // Repair clearly wrong full-period range on circular/elliptic edges.
 match curve {
 Curve3::Circle(c) => {
 let wrap_2pi = |t: f64| -> f64 {
 let mut out = t % two_pi;
 if out < 0.0 {
 out += two_pi;
 }
 out
 };
 let span_hint = t1 - t0;
 let multi_turn = span_hint.abs() > two_pi + TOLERANCE_MESH_LEGACY;
 // Seam edge: same geometric vertex at start/end.
 let seam_tol = (TOLERANCE_ABS * c.radius.max(1.0)).max(edge_tol * 5.0).max(TOLERANCE_LINEAR_RELAX_8);
 let seam = topo_closed || (p_start - p_end).length() <= seam_tol;
 if seam {
 // Only force full turn when source trim already indicates
 // a near-full period.
 if topo_closed {
 let sign = if span_hint < 0.0 { -1.0 } else { 1.0 };
 t0 = 0.0;
 t1 = sign * two_pi;
 } else if near_full_turn {
 let sign = if span_hint < 0.0 { -1.0 } else { 1.0 };
 t0 = 0.0;
 t1 = sign * two_pi;
 }
 } else {
 if near_full_turn {
 let x_ax = rcad_kernel::geom::any_perpendicular(c.normal);
 let y_ax = c.normal.cross(x_ax);
 let v0 = p_start - c.center;
 let v1 = p_end - c.center;
 let a0 = wrap_2pi(v0.dot(y_ax).atan2(v0.dot(x_ax)));
 let a1 = wrap_2pi(v1.dot(y_ax).atan2(v1.dot(x_ax)));
 let reliable_hint = span_hint.signum() * TOLERANCE_ADAPTIVE_MAX;
 let mut dt = choose_arc_delta(a0, a1, reliable_hint);
 if dt.abs() < TOLERANCE_MESH_LEGACY && span_hint.is_finite() && span_hint.abs() > TOLERANCE_ADAPTIVE_MAX {
 let sign = if span_hint < 0.0 { -1.0 } else { 1.0 };
 dt = sign * span_hint.abs().clamp(TOLERANCE_ADAPTIVE_MAX, two_pi - TOLERANCE_MESH_LEGACY);
 }
 let chord = (p_end - p_start).length();
 let chord_tol = TOLERANCE_LINEAR_RELAX_8 * c.radius.max(1.0);
 if chord > chord_tol && dt.abs() < TOLERANCE_ADAPTIVE_MAX {
 let ratio = (chord / (2.0 * c.radius.max(TOLERANCE_LEN_MIN))).clamp(0.0, 1.0);
 let minor = (2.0 * ratio.asin()).clamp(TOLERANCE_MESH_LEGACY, std::f64::consts::PI);
 let sign = if span_hint < 0.0 { -1.0 } else { 1.0 };
 dt = sign * minor;
 }
 t0 = a0;
 t1 = a0 + dt;
 } else if multi_turn {
 // Non-closed edges with >2  span are usually unit/trim artifacts.
 // Rebuild a single-turn arc from endpoints to avoid rosette sampling.
 let x_ax = rcad_kernel::geom::any_perpendicular(c.normal);
 let y_ax = c.normal.cross(x_ax);
 let v0 = p_start - c.center;
 let v1 = p_end - c.center;
 let a0 = wrap_2pi(v0.dot(y_ax).atan2(v0.dot(x_ax)));
 let a1 = wrap_2pi(v1.dot(y_ax).atan2(v1.dot(x_ax)));
 let reliable_hint = span_hint.signum() * two_pi;
 let dt = choose_arc_delta(a0, a1, reliable_hint);
 t0 = a0;
 t1 = a0 + dt;
 }
 }
 }
 Curve3::Ellipse(e) => {
 let wrap_2pi = |t: f64| -> f64 {
 let mut out = t % two_pi;
 if out < 0.0 {
 out += two_pi;
 }
 out
 };
 let span_hint = t1 - t0;
 let multi_turn = span_hint.abs() > two_pi + TOLERANCE_MESH_LEGACY;
 let seam_tol = (TOLERANCE_ABS * e.major_radius.max(e.minor_radius).max(1.0))
 .max(edge_tol * 5.0)
 .max(TOLERANCE_LINEAR_RELAX_8);
 let seam = topo_closed || (p_start - p_end).length() <= seam_tol;
 if seam {
 // Only force full turn when source trim already indicates
 // a near-full period.
 if topo_closed {
 let sign = if span_hint < 0.0 { -1.0 } else { 1.0 };
 t0 = 0.0;
 t1 = sign * two_pi;
 } else if near_full_turn {
 let sign = if span_hint < 0.0 { -1.0 } else { 1.0 };
 t0 = 0.0;
 t1 = sign * two_pi;
 }
 } else {
 if near_full_turn {
 let x_ax = e.major_dir.normalize();
 let y_ax = e.normal.cross(x_ax).normalize();
 let v0 = p_start - e.center;
 let v1 = p_end - e.center;
 let a0 = wrap_2pi((v0.dot(y_ax) / e.minor_radius).atan2(v0.dot(x_ax) / e.major_radius));
 let a1 = wrap_2pi((v1.dot(y_ax) / e.minor_radius).atan2(v1.dot(x_ax) / e.major_radius));
 let reliable_hint = span_hint.signum() * TOLERANCE_ADAPTIVE_MAX;
 let mut dt = choose_arc_delta(a0, a1, reliable_hint);
 if dt.abs() < TOLERANCE_MESH_LEGACY && span_hint.is_finite() && span_hint.abs() > TOLERANCE_ADAPTIVE_MAX {
 let sign = if span_hint < 0.0 { -1.0 } else { 1.0 };
 dt = sign * span_hint.abs().clamp(TOLERANCE_ADAPTIVE_MAX, two_pi - TOLERANCE_MESH_LEGACY);
 }
 t0 = a0;
 t1 = a0 + dt;
 } else if multi_turn {
 // Non-closed edges with >2  span are usually unit/trim artifacts.
 // Rebuild a single-turn arc from endpoints to avoid rosette sampling.
 let x_ax = e.major_dir.normalize();
 let y_ax = e.normal.cross(x_ax).normalize();
 let v0 = p_start - e.center;
 let v1 = p_end - e.center;
 let a0 = wrap_2pi((v0.dot(y_ax) / e.minor_radius).atan2(v0.dot(x_ax) / e.major_radius));
 let a1 = wrap_2pi((v1.dot(y_ax) / e.minor_radius).atan2(v1.dot(x_ax) / e.major_radius));
 let reliable_hint = span_hint.signum() * two_pi;
 let dt = choose_arc_delta(a0, a1, reliable_hint);
 t0 = a0;
 t1 = a0 + dt;
 }
 }
 }
 _ => {}
 }

 let span = (t1 - t0).abs();
 if span > TOLERANCE_LEN_MIN {
 let n_segs = match curve {
 Curve3::Circle(_) => {
 let segs = (span / (2.0 * std::f64::consts::PI) * 64.0).ceil() as usize;
 segs.clamp(4, 64)
 }
 Curve3::Ellipse(_) => 24,
 _ => 16,
 };
 if pts.is_empty() {
 pts.push(p_start);
 }
 for i in 1..=n_segs {
 let t = t0 + (t1 - t0) * (i as f64 / n_segs as f64);
 pts.push(curve.point_at(t));
 }
 // Keep sampled chain anchored to topological edge endpoints.
 if let Some(last) = pts.last_mut() {
 *last = p_end;
 }
 sampled = true;
 }
 }

 if !sampled {
 if pts.is_empty() {
 pts.push(p_start);
 }
 pts.push(p_end);
 }
 }

 // Drop duplicated closing point if present.
 if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length() < TOLERANCE_COORD_SUB {
 pts.pop();
 }

 pts
}

fn ear_clip(pts: &[[f64; 2]]) -> Vec<[usize; 3]> {
 let n = pts.len();
 let mut indices: Vec<usize> = (0..n).collect();
 let mut triangles = Vec::new();

 // Ensure CCW winding
 let area = signed_area_2d(pts, &indices);
 if area < 0.0 {
 indices.reverse();
 }

 let mut remaining = indices;
 while remaining.len() > 3 {
 let len = remaining.len();
 let mut ear_found = false;

 for i in 0..len {
 let prev = if i == 0 { len - 1 } else { i - 1 };
 let next = if i == len - 1 { 0 } else { i + 1 };

 let a = remaining[prev];
 let b = remaining[i];
 let c = remaining[next];

 // Check convexity (left turn)
 if cross_2d(pts[a], pts[b], pts[c]) <= 0.0 {
 continue;
 }

 // Check no other mesh node inside this triangle
 let mut contains_other = false;
 for j in 0..len {
 if j == prev || j == i || j == next {
 continue;
 }
 if point_in_triangle_2d(pts[remaining[j]], pts[a], pts[b], pts[c]) {
 contains_other = true;
 break;
 }
 }

 if !contains_other {
 triangles.push([a, b, c]);
 remaining.remove(i);
 ear_found = true;
 break;
 }
 }

 if !ear_found {
 // Degenerate polygon =emit remaining as fan
 for i in 1..remaining.len() - 1 {
 triangles.push([remaining[0], remaining[i], remaining[i + 1]]);
 }
 break;
 }
 }

 if remaining.len() == 3 {
 triangles.push([remaining[0], remaining[1], remaining[2]]);
 }

 triangles
}

fn signed_area_2d(pts: &[[f64; 2]], indices: &[usize]) -> f64 {
 let n = indices.len();
 let mut area = 0.0;
 for i in 0..n {
 let j = (i + 1) % n;
 let a = pts[indices[i]];
 let b = pts[indices[j]];
 area += a[0] * b[1] - b[0] * a[1];
 }
 area * 0.5
}

fn cross_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
 (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn point_in_triangle_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
 let d1 = cross_2d(a, b, p);
 let d2 = cross_2d(b, c, p);
 let d3 = cross_2d(c, a, p);

 let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
 let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

 !(has_neg && has_pos)
}


/// Tessellate all faces of a BRep in-place, writing triangle indices into
/// `face.triangles`.
///
/// Analogous to OCCT `BRepMesh_IncrementalMesh`.
///
/// For each face:
/// - If the face has an associated `Surface3` in `brep.geom`, the surface is
/// sampled adaptively using `triangulate_surface` with the given `params`.
/// The resulting world-space vertices are appended to `brep.vertices` and
/// the triangle indices are stored in `face.triangles`.
/// - Faces without a surface entry fall back to fan-triangulation of the
/// outer wire vertices (same as the existing rendering path).
///
/// Faces whose [`Face::mesh_dirty`] flag is `false` (clean) are **skipped**
/// unless their `triangles` is empty =allowing incremental updates when only
/// part of the model changes.  To force a full retessellation call
/// [`BRep::invalidate_mesh`] first.
///
/// After tessellating a face its `mesh_dirty` flag is set to `false`.
pub fn mesh_brep(brep: &mut rcad_kernel::BRep, params: &TessellationParams) {
 let mut face_flat_idx = 0usize;

 // Collect all Face TShape indices in flat order.
 let face_tsi: Vec<usize> = brep.tshapes.iter().enumerate()
 .filter_map(|(i, ts)| if matches!(ts.as_ref(), topods::TShape::Face(_)) { Some(i) } else { None })
 .collect();

 // Collect vertex points to add (to avoid borrow conflicts with immutable tshapes iteration).
 let mut new_verts: Vec<DVec3> = Vec::new();

 {
 for ts in &brep.tshapes {
 let topods::TShape::Solid(sd) = ts.as_ref() else { continue };
 for sr in &sd.shells {
 let Some(shd) = (if sr.index < brep.tshapes.len() {
 match &*brep.tshapes[sr.index] { topods::TShape::Shell(sh) => Some(sh), _ => None }
 } else { None }) else { continue };
 for _fr in &shd.faces {
 let tsi = face_tsi.get(face_flat_idx).copied().unwrap_or(usize::MAX);
 let Some(face_data) = (if tsi != usize::MAX {
 match &*brep.tshapes[tsi] { topods::TShape::Face(fd) => Some(fd), _ => None }
 } else { None }) else {
 face_flat_idx += 1;
 continue;
 };

 // Resolve surface and UV domain from TFaceData.
 let surf_and_domain: Option<(Surface3, [f64; 4])> = face_data.surface.clone().map(|surf| {
 let domain = face_data.uv_domain.unwrap_or_else(|| surf.default_domain());
 (surf, domain)
 });

 let mut filled = false;

 if let Some((surf, domain)) = surf_and_domain {
 let [u0, u1, v0, v1] = domain;

 // Clamp infinite domains using vertex projections.
 let (u0, u1, v0, v1) = clamp_domain_to_vertices(
 brep, face_flat_idx, &surf, u0, u1, v0, v1,
 );

 let is_plane = matches!(surf, Surface3::Plane(_));
 let uv_domain_min = uv_polyline_trim_closed_len_sq(u1 - u0, v1 - v0).sqrt();
 let plane_analytic_ok = !is_plane
 || (u1 - u0).abs() * (v1 - v0).abs() >= uv_domain_min * uv_domain_min;

 let domain_ok = (u1 - u0).abs() >= uv_domain_min
 && (v1 - v0).abs() >= uv_domain_min;
 let has_inner_wires = !face_data.inner_wires.is_empty();

 if domain_ok && plane_analytic_ok && !has_inner_wires && !is_plane {
 let mesh = triangulate_surface(&surf, [u0, u1], [v0, v1], params);
 if !mesh.triangles.is_empty() {
 new_verts.extend_from_slice(&mesh.nodes);
 filled = true;
 }
 }
 }

 if !filled {
 // Wire-based triangulation fallback.
 let outer_we: Vec<topology::WireEdge> = (if face_data.outer_wire.index < brep.tshapes.len() {
 match &*brep.tshapes[face_data.outer_wire.index] {
 topods::TShape::Wire(wd) => wd.edges.iter().map(|e| topology::WireEdge::new(e.index, e.orientation.is_forward())).collect(),
 _ => Vec::new(),
 }
 } else { Vec::new() });
 let outer_wire = topology::Wire { edges: outer_we };

 let inner_wires: Vec<topology::Wire> = face_data.inner_wires.iter().map(|iw_sr| {
 let edges = (if iw_sr.index < brep.tshapes.len() {
 match &*brep.tshapes[iw_sr.index] {
 topods::TShape::Wire(wd) => wd.edges.iter().map(|e| topology::WireEdge::new(e.index, e.orientation.is_forward())).collect(),
 _ => Vec::new(),
 }
 } else { Vec::new() });
 topology::Wire { edges }
 }).collect();

 let outer_pts = sample_wire_polygon_points(brep, &outer_wire);
 let hole_pts: Vec<Vec<DVec3>> = inner_wires.iter()
 .map(|wire| sample_wire_polygon_points(brep, wire))
 .filter(|pts| pts.len() >= 3)
 .collect();
 if outer_pts.len() >= 3 {
 let mut poly_pts = outer_pts.clone();
 for hole in &hole_pts {
 poly_pts.extend_from_slice(hole);
 }
 let local_tris = if hole_pts.is_empty() {
 triangulate_polygon(&outer_pts, DVec3::Z)
 } else {
 triangulate_polygon_with_holes(&outer_pts, &hole_pts, DVec3::Z)
 };
 if std::env::var("RCAD_DEBUG_FACE_TRI").is_ok() && !hole_pts.is_empty() {
 eprintln!(
 "[rcad-tri][debug] face_flat={} outer_pts={} holes={} hole_pts_total={} tris={}",
 face_flat_idx, outer_pts.len(), hole_pts.len(),
 hole_pts.iter().map(|h| h.len()).sum::<usize>(), local_tris.len()
 );
 }
 if !local_tris.is_empty() {
 new_verts.extend_from_slice(&poly_pts);
 }
 }
 }

 face_flat_idx += 1;
 }
 }
 }
 }

 // Add collected vertices to the BRep.
 for pt in &new_verts {
 brep.add_tvertex(*pt);
 }
}

/// Absolute span (radians on periodic axes) below which a finite stored trim is treated
/// as degenerate (STEP often ships near-zero boxes =one-strip tessellation).
const DEGENERATE_TRIM_ABS_MIN: f64 = 0.08;
/// Minimum span as a fraction of a full period / natural range before we prefer a hull from wire vertices.
const DEGENERATE_TRIM_REL: f64 = 1.0 / 64.0;

/// Canonical periodic / bounded metadata for analytic surfaces used in [`clamp_domain_to_vertices`].
#[derive(Clone, Copy)]
struct CanonicalUvAxes {
 /// Period in `u` when `u` is periodic (e.g. azimuth / longitude), else `None`.
 u_period: Option<f64>,
 /// Period in `v` when `v` is periodic (torus minor angle), else `None`.
 v_period: Option<f64>,
 /// Natural finite bounds for a non-periodic `v` (e.g. sphere colatitude in `[0,  `).
 v_natural: Option<(f64, f64)>,
}

fn canonical_uv_axes(surf: &Surface3) -> CanonicalUvAxes {
 use std::f64::consts::{PI, TAU};
 match surf {
 Surface3::Cylinder(_) => CanonicalUvAxes {
 u_period: Some(TAU),
 v_period: None,
 v_natural: None,
 },
 Surface3::Cone(_) => CanonicalUvAxes {
 u_period: Some(TAU),
 v_period: None,
 v_natural: None,
 },
 Surface3::Sphere(_) => CanonicalUvAxes {
 u_period: Some(TAU),
 v_period: None,
 v_natural: Some((0.0, PI)),
 },
 Surface3::Torus(_) => CanonicalUvAxes {
 u_period: Some(TAU),
 v_period: Some(TAU),
 v_natural: None,
 },
 Surface3::Trimmed(t) => canonical_uv_axes(t.basis.as_ref()),
 Surface3::Offset(o) => canonical_uv_axes(o.basis.as_ref()),
 _ => CanonicalUvAxes {
 u_period: None,
 v_period: None,
 v_natural: None,
 },
 }
}

fn span_too_small_for_axis(span: f64, period: Option<f64>, natural: Option<(f64, f64)>) -> bool {
 if span < TOLERANCE_LEN_MIN {
 return true;
 }
 if let Some(p) = period {
 let thr = DEGENERATE_TRIM_ABS_MIN.max(p * DEGENERATE_TRIM_REL);
 return span < thr;
 }
 if let Some((lo, hi)) = natural
 && lo.is_finite() && hi.is_finite() {
 let range = (hi - lo).abs();
 if range > TOLERANCE_LEN_MIN {
 let thr = DEGENERATE_TRIM_ABS_MIN.max(range * DEGENERATE_TRIM_REL);
 return span < thr;
 }
 }
 false
}

/// Unwrap a sequence of angles in wire order so min/max reflect a contiguous patch (same idea as STEP seam handling).
fn unwrap_1d_periodic_chain(vals: &mut [f64], period: f64) {
 if vals.len() < 2 {
 return;
 }
 let mut offset = 0.0;
 let mut previous = vals[0];
 for v in vals.iter_mut().skip(1) {
 let raw = *v + offset;
 let delta = raw - previous;
 if delta > period * 0.5 {
 offset -= period;
 } else if delta < -period * 0.5 {
 offset += period;
 }
 *v += offset;
 previous = *v;
 }
}

/// UV bounding box from boundary sample nodes (with margins), using the same projection as offset/PCurve code.
fn hull_uv_box_from_wire(surf: &Surface3, pts: &[DVec3]) -> Option<(f64, f64, f64, f64)> {
 if pts.is_empty() {
 return None;
 }
 let ax = canonical_uv_axes(surf);
 let mut uv_pairs: Vec<[f64; 2]> = Vec::with_capacity(pts.len());
 for &p in pts {
 uv_pairs.push(project_point_to_surface_uv(p, surf, None)?);
 }
 let mut us: Vec<f64> = uv_pairs.iter().map(|a| a[0]).collect();
 let mut vs: Vec<f64> = uv_pairs.iter().map(|a| a[1]).collect();
 if let Some(per) = ax.u_period {
 unwrap_1d_periodic_chain(&mut us, per);
 }
 if let Some(per) = ax.v_period {
 unwrap_1d_periodic_chain(&mut vs, per);
 }
 let u_min = us.iter().copied().fold(f64::INFINITY, f64::min);
 let u_max = us.iter().copied().fold(f64::NEG_INFINITY, f64::max);
 let v_min = vs.iter().copied().fold(f64::INFINITY, f64::min);
 let v_max = vs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
 if !u_min.is_finite() || !v_min.is_finite() {
 return None;
 }
 let uv_floor = uv_polyline_trim_closed_len_sq(u_max - u_min, v_max - v_min).sqrt();
 let mu = (u_max - u_min).abs() * 0.05 + uv_floor;
 let mv = (v_max - v_min).abs() * 0.05 + uv_floor;
 let nu0 = u_min - mu;
 let nu1 = u_max + mu;
 let mut nv0 = v_min - mv;
 let mut nv1 = v_max + mv;
 if let Some((lo, hi)) = ax.v_natural {
 nv0 = nv0.max(lo);
 nv1 = nv1.min(hi);
 }
 Some((nu0, nu1, nv0, nv1))
}

fn fallback_infinite_domain(u0: f64, u1: f64, v0: f64, v1: f64) -> (f64, f64, f64, f64) {
 let eff_u0 = if u0.is_finite() { u0 } else { -10.0 };
 let eff_u1 = if u1.is_finite() { u1 } else { 10.0 };
 let eff_v0 = if v0.is_finite() { v0 } else { -10.0 };
 let eff_v1 = if v1.is_finite() { v1 } else { 10.0 };
 (eff_u0, eff_u1, eff_v0, eff_v1)
}

/// Clamp a potentially infinite UV domain to the range spanned by the face's
/// wire vertices projected onto the surface parameters.
fn clamp_domain_to_vertices(
 brep: &rcad_kernel::BRep,
 face_flat_idx: usize,
 surf: &Surface3,
 u0: f64, u1: f64, v0: f64, v1: f64,
) -> (f64, f64, f64, f64) {
 let ax = canonical_uv_axes(surf);
 let need_u = !u0.is_finite() || !u1.is_finite();
 let need_v = !v0.is_finite() || !v1.is_finite();
 let du = if u0.is_finite() && u1.is_finite() {
 (u1 - u0).abs()
 } else {
 f64::NAN
 };
 let dv = if v0.is_finite() && v1.is_finite() {
 (v1 - v0).abs()
 } else {
 f64::NAN
 };

 let u_want_hull =
 need_u || (du.is_finite() && span_too_small_for_axis(du, ax.u_period, None));
 let v_want_hull = need_v
 || (dv.is_finite() && span_too_small_for_axis(dv, ax.v_period, ax.v_natural));

 if !u_want_hull && !v_want_hull {
 return (u0, u1, v0, v1);
 }

 // Collect sampled wire points for this face so curved boundaries
 // (e.g. cylinder/cone circles) contribute to UV hull estimation.
 // Find Face at face_flat_idx in TShape order.
 let face_tsi: Vec<usize> = brep.tshapes.iter().enumerate()
 .filter_map(|(i, ts)| if matches!(ts.as_ref(), topods::TShape::Face(_)) { Some(i) } else { None })
 .collect();
 let Some(&face_tsi) = face_tsi.get(face_flat_idx) else {
 return (u0, u1, v0, v1);
 };
 let topods::TShape::Face(fd) = &*brep.tshapes[face_tsi] else {
 return (u0, u1, v0, v1);
 };

 // Build a dangling topology::Wire from the outer_wire ShapeRef for sampling.
 let outer_we: Vec<topology::WireEdge> = (if fd.outer_wire.index < brep.tshapes.len() {
 match &*brep.tshapes[fd.outer_wire.index] {
 topods::TShape::Wire(wd) => wd.edges.iter().map(|e| topology::WireEdge::new(e.index, e.orientation.is_forward())).collect(),
 _ => Vec::new(),
 }
 } else { Vec::new() });
 let outer_wire_w = topology::Wire { edges: outer_we };

 let mut pts = sample_wire_polygon_points(brep, &outer_wire_w);
 if pts.is_empty() {
 // Fallback to topological endpoints when wire sampling is unavailable.
 pts = outer_wire_w
 .edges
 .iter()
 .filter_map(|we| {
 brep.tshapes.get(we.idx).and_then(|ts| {
 if let topods::TShape::Edge(ed) = ts.as_ref() {
 let vi = if we.forward { ed.first.index } else { ed.last.index };
 brep.vertex_point(vi)
 } else { None }
 })
 })
 .collect();
 }

 if pts.is_empty() {
 return (u0, u1, v0, v1);
 }

 match surf {
 Surface3::Plane(plane) => {
 // Project vertices onto the plane's local UV frame.
 let u_ax = any_perpendicular(plane.normal);
 let v_ax = plane.normal.cross(u_ax).normalize_or_zero();
 let us: Vec<f64> = pts.iter().map(|&p| (p - plane.origin).dot(u_ax)).collect();
 let vs: Vec<f64> = pts.iter().map(|&p| (p - plane.origin).dot(v_ax)).collect();
 let pu0 = us.iter().cloned().fold(f64::INFINITY, f64::min);
 let pu1 = us.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 let pv0 = vs.iter().cloned().fold(f64::INFINITY, f64::min);
 let pv1 = vs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
 let uv_floor = uv_polyline_trim_closed_len_sq(pu1 - pu0, pv1 - pv0).sqrt();
 let mu = (pu1 - pu0).abs() * 0.05 + uv_floor;
 let mv = (pv1 - pv0).abs() * 0.05 + uv_floor;
 (pu0 - mu, pu1 + mu, pv0 - mv, pv1 + mv)
 }
 Surface3::Cylinder(_)
 | Surface3::Cone(_)
 | Surface3::Sphere(_)
 | Surface3::Torus(_)
 | Surface3::Trimmed(_)
 | Surface3::Offset(_) => {
 if let Some((nu0, nu1, nv0, nv1)) = hull_uv_box_from_wire(surf, &pts) {
 let out_u0 = if u_want_hull { nu0 } else { u0 };
 let out_u1 = if u_want_hull { nu1 } else { u1 };
 let out_v0 = if v_want_hull { nv0 } else { v0 };
 let out_v1 = if v_want_hull { nv1 } else { v1 };
 (out_u0, out_u1, out_v0, out_v1)
 } else {
 fallback_infinite_domain(u0, u1, v0, v1)
 }
 }
 _ => fallback_infinite_domain(u0, u1, v0, v1),
 }
}

// ============================================================================
// Mesh quality metrics
// ============================================================================

/// Scalar summaries for a triangle mesh (aspect ratio, edge lengths, areas).
#[derive(Debug, Clone, Default)]
pub struct MeshQualityMetrics {
 /// Number of triangles.
 pub triangle_count: usize,
 /// Number of mesh nodes referenced by the mesh.
 pub node_count: usize,
 /// Mean edge aspect ratio (max/min edge per triangle).
 pub average_aspect_ratio: f64,
 /// Worst-case aspect ratio among triangles.
 pub max_aspect_ratio: f64,
 /// Count of triangles with aspect ratio above the diagnostic threshold.
 pub poor_aspect_ratio_count: usize,
 /// Mean edge length (all edges of all triangles).
 pub average_edge_length: f64,
 /// Standard deviation of edge lengths.
 pub edge_length_stddev: f64,
 /// Shortest edge in the mesh.
 pub min_edge_length: f64,
 /// Longest edge in the mesh.
 pub max_edge_length: f64,
 /// Mean triangle area.
 pub average_area: f64,
 /// Standard deviation of triangle areas.
 pub area_stddev: f64,
 /// Smallest triangle area.
 pub min_area: f64,
 /// Largest triangle area.
 pub max_area: f64,
 /// Triangles with near-zero area.
 pub degenerate_count: usize,
}

impl MeshQualityMetrics {
 /// Heuristic = ood mesh=check against a maximum allowed aspect ratio.
 pub fn is_good(&self, max_aspect_ratio: f64) -> bool {
 self.degenerate_count == 0
 && self.max_aspect_ratio <= max_aspect_ratio
 && (self.triangle_count <= 10 || self.poor_aspect_ratio_count < self.triangle_count / 10)
 }

 /// Composite score in `[0, 1]` from aspect ratio, degenerates, and edge-length CV.
 pub fn quality_score(&self) -> f64 {
 if self.triangle_count == 0 {
 return 0.0;
 }

 let aspect_score = if self.max_aspect_ratio > 0.0 {
 (10.0 / self.max_aspect_ratio).min(1.0)
 } else {
 0.0
 };

 let degenerate_ratio = self.degenerate_count as f64 / self.triangle_count as f64;
 let degenerate_score = 1.0 - degenerate_ratio;

 let uniformity_score = if self.average_edge_length > 0.0 {
 let cv = self.edge_length_stddev / self.average_edge_length;
 (1.0 - cv).max(0.0)
 } else {
 0.0
 };

 (aspect_score * 0.4 + degenerate_score * 0.4 + uniformity_score * 0.2).clamp(0.0, 1.0)
 }
}

/// Compute [`MeshQualityMetrics`] for raw node/triangle buffers.
pub fn compute_mesh_quality(nodes: &[DVec3], triangles: &[[usize; 3]]) -> MeshQualityMetrics {
 if nodes.is_empty() || triangles.is_empty() {
 return MeshQualityMetrics::default();
 }

 let mut metrics = MeshQualityMetrics {
 triangle_count: triangles.len(),
 node_count: nodes.len(),
 ..Default::default()
 };

 let mut aspect_ratios: Vec<f64> = Vec::with_capacity(triangles.len());
 let mut areas: Vec<f64> = Vec::with_capacity(triangles.len());
 let mut edge_lengths: Vec<f64> = Vec::new();

 for &tri in triangles {
 let [i0, i1, i2] = tri;
 if i0 >= nodes.len() || i1 >= nodes.len() || i2 >= nodes.len() {
 continue;
 }

 let p0 = nodes[i0];
 let p1 = nodes[i1];
 let p2 = nodes[i2];

 // Edge lengths
 let e0 = (p1 - p0).length();
 let e1 = (p2 - p1).length();
 let e2 = (p0 - p2).length();

 edge_lengths.push(e0);
 edge_lengths.push(e1);
 edge_lengths.push(e2);

 // Area via Heron's formula
 let s = (e0 + e1 + e2) * 0.5;
 let area = if s > 0.0 {
 (s * (s - e0) * (s - e1) * (s - e2)).sqrt()
 } else {
 0.0
 };

 areas.push(area);

 // Degenerate (near-collinear) triangles
 if area < TOLERANCE_LEN_MIN {
 metrics.degenerate_count += 1;
 }

 // Aspect ratio = longest / shortest edge
 let max_edge = e0.max(e1).max(e2);
 let min_edge = e0.min(e1).min(e2);
 let aspect_ratio = if min_edge > TOLERANCE_LEN_MIN {
 max_edge / min_edge
 } else {
 f64::INFINITY
 };
 aspect_ratios.push(aspect_ratio);
 }

 // Aggregate statistics
 if !aspect_ratios.is_empty() {
 metrics.max_aspect_ratio = aspect_ratios.iter().cloned().fold(0.0, f64::max);
 metrics.average_aspect_ratio = aspect_ratios.iter().sum::<f64>() / aspect_ratios.len() as f64;
 metrics.poor_aspect_ratio_count = aspect_ratios.iter().filter(|&&ar| ar > 20.0).count();
 }

 if !edge_lengths.is_empty() {
 metrics.min_edge_length = edge_lengths.iter().cloned().fold(f64::INFINITY, f64::min);
 metrics.max_edge_length = edge_lengths.iter().cloned().fold(0.0, f64::max);
 metrics.average_edge_length = edge_lengths.iter().sum::<f64>() / edge_lengths.len() as f64;

 let variance = edge_lengths.iter()
 .map(|&l| (l - metrics.average_edge_length).powi(2))
 .sum::<f64>() / edge_lengths.len() as f64;
 metrics.edge_length_stddev = variance.sqrt();
 }

 if !areas.is_empty() {
 metrics.min_area = areas.iter().cloned().fold(f64::INFINITY, f64::min);
 metrics.max_area = areas.iter().cloned().fold(0.0, f64::max);
 metrics.average_area = areas.iter().sum::<f64>() / areas.len() as f64;

 let variance = areas.iter()
 .map(|&a| (a - metrics.average_area).powi(2))
 .sum::<f64>() / areas.len() as f64;
 metrics.area_stddev = variance.sqrt();
 }

 metrics
}

impl SurfaceMesh {
 /// Convenience wrapper around [`compute_mesh_quality`].
 pub fn compute_quality(&self) -> MeshQualityMetrics {
 compute_mesh_quality(&self.nodes, &self.triangles)
 }
}

// ============================================================================
// Adaptive mesh subdivision
// ============================================================================

/// Midpoint / 4-to-1 split rules driven by normal variation or edge length.
#[derive(Debug, Clone)]
pub struct AdaptiveSubdivider {
 /// Split an edge when the angle between endpoint normals exceeds this (radians).
 pub curvature_threshold: f64,
 /// Split an edge when its length exceeds this (world units).
 pub distance_threshold: f64,
 /// Maximum recursion depth for uniform splits (reserved for future use).
 pub max_subdivision_levels: usize,
 /// Reserved flag for boundary-aware splitting.
 pub preserve_boundary: bool,
}

impl Default for AdaptiveSubdivider {
 fn default() -> Self {
 Self {
 curvature_threshold: 0.1,  // ~5.7 degrees
 distance_threshold: 0.1,
 max_subdivision_levels: 3,
 preserve_boundary: true,
 }
 }
}

impl AdaptiveSubdivider {
 /// Default-configured subdivider.
 pub fn new() -> Self {
 Self::default()
 }

 /// Builder: set [`Self::curvature_threshold`].
 pub fn with_curvature_threshold(mut self, threshold: f64) -> Self {
 self.curvature_threshold = threshold;
 self
 }

 /// Builder: set [`Self::distance_threshold`].
 pub fn with_distance_threshold(mut self, threshold: f64) -> Self {
 self.distance_threshold = threshold;
 self
 }

 /// Builder: set [`Self::max_subdivision_levels`].
 pub fn with_max_levels(mut self, levels: usize) -> Self {
 self.max_subdivision_levels = levels;
 self
 }

 /// Split triangles whose edges exceed the normal-difference threshold.
 pub fn subdivide_by_curvature(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
 if mesh.triangles.is_empty() || mesh.normals.is_empty() {
 return mesh.clone();
 }

 let mut nodes = mesh.nodes.clone();
 let mut normals = mesh.normals.clone();
 let mut triangles = Vec::new();

 // Canonical edge -> new midpoint vertex index
 let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();

 for &tri in &mesh.triangles {
 let [i0, i1, i2] = tri;

 // Per-edge normal change
 let n0 = normals.get(i0).copied().unwrap_or(DVec3::ZERO);
 let n1 = normals.get(i1).copied().unwrap_or(DVec3::ZERO);
 let n2 = normals.get(i2).copied().unwrap_or(DVec3::ZERO);

 let split_01 = self.should_split_by_curvature(n0, n1);
 let split_12 = self.should_split_by_curvature(n1, n2);
 let split_20 = self.should_split_by_curvature(n2, n0);

 if split_01 || split_12 || split_20 {
 self.subdivide_triangle(
 tri,
 &mut nodes,
 &mut normals,
 &mut triangles,
 &mut edge_midpoints,
 );
 } else {
 triangles.push(tri);
 }
 }

 SurfaceMesh {
 nodes,
 triangles,
 normals,
 dirty: false,
 }
 }

 /// Split triangles whose edges exceed [`Self::distance_threshold`].
 pub fn subdivide_by_distance(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
 if mesh.triangles.is_empty() {
 return mesh.clone();
 }

 let mut nodes = mesh.nodes.clone();
 let mut normals = mesh.normals.clone();
 let mut triangles = Vec::new();

 let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();

 for &tri in &mesh.triangles {
 let [i0, i1, i2] = tri;

 let p0 = nodes[i0];
 let p1 = nodes[i1];
 let p2 = nodes[i2];

 let split_01 = self.should_split_by_distance(p0, p1);
 let split_12 = self.should_split_by_distance(p1, p2);
 let split_20 = self.should_split_by_distance(p2, p0);

 if split_01 || split_12 || split_20 {
 self.subdivide_triangle(
 tri,
 &mut nodes,
 &mut normals,
 &mut triangles,
 &mut edge_midpoints,
 );
 } else {
 triangles.push(tri);
 }
 }

 SurfaceMesh {
 nodes,
 triangles,
 normals,
 dirty: false,
 }
 }

 fn should_split_by_curvature(&self, n0: DVec3, n1: DVec3) -> bool {
 let len0 = n0.length();
 let len1 = n1.length();
 if len0 < 0.5 || len1 < 0.5 {
 return false;
 }
 let cos_angle = (n0.dot(n1) / (len0 * len1)).clamp(-1.0, 1.0);
 let angle = cos_angle.acos();
 angle > self.curvature_threshold
 }

 fn should_split_by_distance(&self, p0: DVec3, p1: DVec3) -> bool {
 (p1 - p0).length() > self.distance_threshold
 }

 fn subdivide_triangle(
 &self,
 tri: [usize; 3],
 nodes: &mut Vec<DVec3>,
 normals: &mut Vec<DVec3>,
 triangles: &mut Vec<[usize; 3]>,
 edge_midpoints: &mut HashMap<(usize, usize), usize>,
 ) {
 let [i0, i1, i2] = tri;

 let p0 = nodes[i0];
 let p1 = nodes[i1];
 let p2 = nodes[i2];

 let n0 = normals.get(i0).copied().unwrap_or(DVec3::ZERO);
 let n1 = normals.get(i1).copied().unwrap_or(DVec3::ZERO);
 let n2 = normals.get(i2).copied().unwrap_or(DVec3::ZERO);

 // Reuse midpoints on shared edges
 let m01 = self.get_or_create_midpoint(i0, i1, p0, p1, n0, n1, nodes, normals, edge_midpoints);
 let m12 = self.get_or_create_midpoint(i1, i2, p1, p2, n1, n2, nodes, normals, edge_midpoints);
 let m20 = self.get_or_create_midpoint(i2, i0, p2, p0, n2, n0, nodes, normals, edge_midpoints);

 // Four-way split =four triangles
 triangles.push([i0, m01, m20]);
 triangles.push([m01, i1, m12]);
 triangles.push([m20, m12, i2]);
 triangles.push([m01, m12, m20]);
 }

 fn get_or_create_midpoint(
 &self,
 i0: usize,
 i1: usize,
 p0: DVec3,
 p1: DVec3,
 n0: DVec3,
 n1: DVec3,
 nodes: &mut Vec<DVec3>,
 normals: &mut Vec<DVec3>,
 edge_midpoints: &mut HashMap<(usize, usize), usize>,
 ) -> usize {
 let key = if i0 < i1 { (i0, i1) } else { (i1, i0) };

 if let Some(&idx) = edge_midpoints.get(&key) {
 return idx;
 }

 let mid_point = (p0 + p1) * 0.5;
 let mid_normal = (n0 + n1).normalize_or_zero();

 let idx = nodes.len();
 nodes.push(mid_point);
 normals.push(mid_normal);
 edge_midpoints.insert(key, idx);
 idx
 }
}

// ============================================================================
// Boundary-aware tessellation helpers
// ============================================================================

/// Sharp crease between two triangles, expressed as an undirected edge.
#[derive(Debug, Clone)]
pub struct FeatureEdge {
 /// Edge start node index.
 pub start_node: usize,
 /// Edge end node index.
 pub end_node: usize,
 /// Dihedral angle between adjacent faces (radians).
 pub feature_angle: f64,
}

/// Detects crease edges and optionally protects them during welding.
#[derive(Debug, Clone)]
pub struct BoundarySensitiveTessellator {
 /// Dihedral angle above which an edge is treated as a feature.
 pub feature_angle_threshold: f64,
 /// Manually supplied or auto-detected crease edges.
 pub feature_edges: Vec<FeatureEdge>,
 /// When `true`, [`Self::detect_feature_edges`] fills `feature_edges`.
 pub auto_detect_features: bool,
}

impl Default for BoundarySensitiveTessellator {
 fn default() -> Self {
 Self {
 feature_angle_threshold: 0.52,  // ~30 degrees
 feature_edges: Vec::new(),
 auto_detect_features: true,
 }
 }
}

impl BoundarySensitiveTessellator {
 /// Default tessellator with ~30 ?crease threshold.
 pub fn new() -> Self {
 Self::default()
 }

 /// Builder: override [`Self::feature_angle_threshold`].
 pub fn with_feature_angle(mut self, angle: f64) -> Self {
 self.feature_angle_threshold = angle;
 self
 }

 /// Append a user-defined crease edge.
 pub fn add_feature_edge(mut self, start: usize, end: usize, angle: f64) -> Self {
 self.feature_edges.push(FeatureEdge {
 start_node: start,
 end_node: end,
 feature_angle: angle,
 });
 self
 }

 /// Populate `feature_edges` from mesh dihedral angles (internal edges only).
 pub fn detect_feature_edges(&mut self, nodes: &[DVec3], triangles: &[[usize; 3]], _normals: &[DVec3]) {
 if !self.auto_detect_features {
 return;
 }

 self.feature_edges.clear();

 // Edge -> incident triangle indices
 let mut edge_to_tris: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
 for (tri_idx, &tri) in triangles.iter().enumerate() {
 let edges = [
 (tri[0].min(tri[1]), tri[0].max(tri[1])),
 (tri[1].min(tri[2]), tri[1].max(tri[2])),
 (tri[2].min(tri[0]), tri[2].max(tri[0])),
 ];
 for edge in edges {
 edge_to_tris.entry(edge).or_default().push(tri_idx);
 }
 }

 // Two-manifold interior edges only
 for (edge, tri_indices) in &edge_to_tris {
 if tri_indices.len() == 2 {
 let tri0 = &triangles[tri_indices[0]];
 let tri1 = &triangles[tri_indices[1]];

 // Face normals
 let n0 = compute_triangle_normal(nodes, tri0);
 let n1 = compute_triangle_normal(nodes, tri1);

 // Dihedral angle
 let cos_angle = n0.dot(n1).clamp(-1.0, 1.0);
 let angle = cos_angle.acos();

 if angle > self.feature_angle_threshold {
 self.feature_edges.push(FeatureEdge {
 start_node: edge.0,
 end_node: edge.1,
 feature_angle: angle,
 });
 }
 }
 }
 }

 /// Re-weld nodes while pinning crease endpoints from `feature_edges`.
 pub fn preserve_feature_edges(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
 if self.feature_edges.is_empty() {
 return mesh.clone();
 }

 // Nodes incident on a crease must not snap to neighbors
 let feature_nodes: std::collections::HashSet<usize> = self.feature_edges.iter()
 .flat_map(|e| [e.start_node, e.end_node])
 .collect();

 // Custom weld pass with exclusions
 let mut result = mesh.clone();
 result = weld_surface_mesh_nodes_with_exclusion(&result, &feature_nodes);
 result
 }
}

fn compute_triangle_normal(nodes: &[DVec3], tri: &[usize; 3]) -> DVec3 {
 if tri[0] >= nodes.len() || tri[1] >= nodes.len() || tri[2] >= nodes.len() {
 return DVec3::Z;
 }
 let p0 = nodes[tri[0]];
 let p1 = nodes[tri[1]];
 let p2 = nodes[tri[2]];
 (p1 - p0).cross(p2 - p0).normalize_or_zero()
}

include!("extra.rs");

