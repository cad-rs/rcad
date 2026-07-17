//! Bounding Volume Hierarchy (BVH) for spatial acceleration.
//!
//! Built with the SAH (Surface Area Heuristic) to speed up:
//! - Ray picking
//! - `min_distance` between shapes
//! - Gap / overlap detection (`detect_gaps_overlaps`)
//!
//! Analogous to OCCT `BVH_Tree` / `BVH_Builder`.


use glam::DVec3;
use rcad_kernel::topods::{self, TShape};
use rcad_kernel::geom::{Surface3, SurfaceEval};

use crate::tolerance::*;

/// Axis-aligned bounding box (AABB).
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
 pub min: DVec3,
 pub max: DVec3,
 pub gap: f64,
}

impl Aabb {
 /// Empty AABB (`min > max`, contains no points).
 pub fn empty() -> Self {
 Self {
 min: DVec3::splat(f64::INFINITY),
 max: DVec3::splat(f64::NEG_INFINITY),
 gap: 0.0,
 }
 }

 /// Build an AABB that encloses all given points.
 pub fn from_points(pts: &[DVec3]) -> Self {
 let mut aabb = Self::empty();
 for &p in pts {
 aabb.expand_point(p);
 }
 aabb
 }

 /// Expand to include one point.
 pub fn expand_point(&mut self, p: DVec3) {
 self.min = self.min.min(p);
 self.max = self.max.max(p);
 }

 /// Expand to include another AABB.
 pub fn expand_aabb(&mut self, other: &Aabb) {
 self.min = self.min.min(other.min);
 self.max = self.max.max(other.max);
 }

 /// AABB center.
 pub fn center(&self) -> DVec3 {
 (self.min + self.max) * 0.5
 }

 /// Surface area (for SAH cost).
 pub fn surface_area(&self) -> f64 {
 let d = self.max - self.min;
 if d.x < 0.0 || d.y < 0.0 || d.z < 0.0 {
 return 0.0; // empty AABB
 }
 2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
 }

 /// Whether this AABB intersects another.
 pub fn intersects(&self, other: &Aabb) -> bool {
 // OCCT Bnd_Box::IsOut
 self.min.x - self.gap <= other.max.x + other.gap
  && self.max.x + self.gap >= other.min.x - other.gap
  && self.min.y - self.gap <= other.max.y + other.gap
  && self.max.y + self.gap >= other.min.y - other.gap
  && self.min.z - self.gap <= other.max.z + other.gap
  && self.max.z + self.gap >= other.min.z - other.gap
 }

 /// Ray–AABB intersection; returns entry parameter `t` along the ray (forward hits only).
 pub fn ray_intersect(&self, origin: DVec3, inv_dir: DVec3) -> Option<f64> {
 let t1 = (self.min - origin) * inv_dir;
 let t2 = (self.max - origin) * inv_dir;

 let t_min = t1.min(t2);
 let t_max = t1.max(t2);

 let t_enter = t_min.x.max(t_min.y).max(t_min.z);
 let t_exit = t_max.x.min(t_max.y).min(t_max.z);

 if t_exit >= t_enter.max(0.0) {
 Some(t_enter.max(0.0))
 } else {
 None
 }
 }

 /// Squared minimum distance from a point to this AABB.
 pub fn point_dist_sq(&self, p: DVec3) -> f64 {
 let clamped = p.clamp(self.min, self.max);
 (p - clamped).length_squared()
 }

 /// Whether a point lies inside the AABB (inclusive of the boundary).
 pub fn contains_point(&self, p: DVec3) -> bool {
 p.cmpge(self.min).all() && p.cmple(self.max).all()
 }
}

/// BVH node (internal or leaf).
#[derive(Debug, Clone)]
enum BvhNode {
 /// Leaf: holds a range of face indices.
 Leaf {
 aabb: Aabb,
 /// Face index range `[start, end)` into `Bvh.face_indices`.
 start: usize,
 end: usize,
 },
 /// Internal node: indices of left/right children in `Bvh.nodes`.
 Internal {
 aabb: Aabb,
 left: usize,
 right: usize,
 },
}

impl BvhNode {
 fn aabb(&self) -> &Aabb {
 match self {
 BvhNode::Leaf { aabb, .. } => aabb,
 BvhNode::Internal { aabb, .. } => aabb,
 }
 }
}

/// BVH over the faces of one `rcad_kernel::BRep`.
///
/// After construction, supports accelerated ray casts, nearest-face queries, etc.
pub struct Bvh {
 /// Node array (root at index `0`).
 nodes: Vec<BvhNode>,
 /// Face index array (leaves reference `[start, end)` ranges).
 face_indices: Vec<usize>,
 /// Per-face AABBs (indexed by original face index).
 face_aabbs: Vec<Aabb>,
 /// Per-face centers (for SAH splits).
 face_centers: Vec<DVec3>,
}

/// Maximum number of faces per leaf.
const MAX_LEAF_SIZE: usize = 4;

/// Number of SAH bucket boundaries sampled along each axis.
const SAH_BUCKETS: usize = 8;

impl Bvh {
 /// Face count (same order as [`Bvh::build`] source iteration).
 pub fn face_count(&self) -> usize {
 self.face_aabbs.len()
 }

 /// Conservative bounds for `face_index` as used in [`Bvh::candidate_pairs`].
 pub fn face_aabb(&self, face_index: usize) -> Option<&Aabb> {
 self.face_aabbs.get(face_index)
 }

 /// Build a BVH over all faces of `brep`.
 ///
 /// Sampling: each face uses boundary vertices plus a small grid of interior samples
 /// on curved patches so the AABB conservatively covers the face.
 pub fn build(brep: &rcad_kernel::BRep) -> Self {
 // Collect face tshape indices.
 let face_ts_indices: Vec<usize> = brep.tshapes.iter().enumerate()
 .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_)))
 .map(|(fi, _)| fi)
 .collect();
 let n_faces = face_ts_indices.len();

 let mut face_aabbs = Vec::with_capacity(n_faces);
 let mut face_centers = Vec::with_capacity(n_faces);

 for &fi in &face_ts_indices {
 let ts = &brep.tshapes[fi];
 let TShape::Face(fd) = ts.as_ref() else { continue };
 let mut aabb = Aabb::empty();

 // Seed AABB from boundary vertices via outer wire edges.
 if let Some(wts) = brep.tshapes.get(fd.outer_wire.index) {
 if let TShape::Wire(wd) = wts.as_ref() {
 for er in &wd.edges {
 if let Some(ets) = brep.tshapes.get(er.index) {
 if let TShape::Edge(ed) = ets.as_ref() {
 if let Some(p0) = brep.vertex_point(ed.first.index) {
 aabb.expand_point(p0);
 }
 if let Some(p1) = brep.vertex_point(ed.last.index) {
 aabb.expand_point(p1);
 }
 }
 }
 }
 }
 }

 // OCCT : rcad_kernel::BRepBndLib::Add (rcad_kernel::BRepBndLib.cxx L83+) — AABB。
 // OCCT BndLib_AddSurface (BndLib_AddSurface.cxx L275-306):
 // Sphere: ± (L299-306); Cylinder/Cone/Torus: BndLib (L287-297)
 // Plane: 4 (L278-286); BSpline/Bezier:  /  (L315-316)
 if let Some(surf) = &fd.surface {
 match surf {
 Surface3::Sphere(s) => {
 let r = s.radius.abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
 aabb.expand_point(s.center - DVec3::splat(r));
 aabb.expand_point(s.center + DVec3::splat(r));
 }
 Surface3::Cylinder(c) => {
 let domain = surf.default_domain();
 let [_, _, v0, v1] = domain;
 let ax = c.axis.normalize_or_zero();
 let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let u_dir = ax.cross(perp).normalize_or_zero();
 let v_dir = ax.cross(u_dir).normalize_or_zero();
 let r = c.radius.abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
 for &vh in &[v0, v1] {
 for k in 0..8 {
 let a = std::f64::consts::TAU * k as f64 / 8.0;
 let p = c.origin + ax * vh + u_dir * r * a.cos() + v_dir * r * a.sin();
 aabb.expand_point(p);
 }
 }
 }
 Surface3::Cone(c) => {
 let domain = surf.default_domain();
 let [_, _, v0, v1] = domain;
 let ax = c.axis.normalize_or_zero();
 let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let u_dir = ax.cross(perp).normalize_or_zero();
 let v_dir = ax.cross(u_dir).normalize_or_zero();
 for &vh in &[v0, v1] {
 let r_at = (c.radius + vh * c.half_angle_rad.tan()).abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
 let center = c.apex + ax * vh;
 for k in 0..8 {
 let a = std::f64::consts::TAU * k as f64 / 8.0;
 let p = center + u_dir * r_at * a.cos() + v_dir * r_at * a.sin();
 aabb.expand_point(p);
 }
 }
 }
 Surface3::Torus(t) => {
 let r_out = t.major_radius.abs() + t.minor_radius.abs() + TOLERANCE_LINEAR_ULTRA_STRICT;
 let ax = t.axis.normalize_or_zero();
 let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
 let u_dir = ax.cross(perp).normalize_or_zero();
 let v_dir = ax.cross(u_dir).normalize_or_zero();
 for k in 0..8 {
 let a = std::f64::consts::TAU * k as f64 / 8.0;
 let c = t.center + u_dir * t.major_radius * a.cos() + v_dir * t.major_radius * a.sin();
 aabb.expand_point(c + ax * t.minor_radius);
 aabb.expand_point(c - ax * t.minor_radius);
 }
 }
 // Plane  : + UV  
 _ => {
 let domain = surf.default_domain();
 let [u0, u1, v0, v1] = domain;
 for i in 0..=2 {
 for j in 0..=2 {
 let u = u0 + (u1 - u0) * i as f64 / 2.0;
 let v = v0 + (v1 - v0) * j as f64 / 2.0;
 let p = surf.point_at(u, v);
 if p.is_finite() { aabb.expand_point(p); }
 }
 }
 }
 }
 }

 // Degenerate faces: nudge AABB to non-zero extent
 let size = aabb.max - aabb.min;
 if size.x < TOLERANCE_LINEAR_ULTRA_STRICT {
 aabb.min.x -= TOLERANCE_LINEAR_ULTRA_STRICT;
 aabb.max.x += TOLERANCE_LINEAR_ULTRA_STRICT;
 }
 if size.y < TOLERANCE_LINEAR_ULTRA_STRICT {
 aabb.min.y -= TOLERANCE_LINEAR_ULTRA_STRICT;
 aabb.max.y += TOLERANCE_LINEAR_ULTRA_STRICT;
 }
 if size.z < TOLERANCE_LINEAR_ULTRA_STRICT {
 aabb.min.z -= TOLERANCE_LINEAR_ULTRA_STRICT;
 aabb.max.z += TOLERANCE_LINEAR_ULTRA_STRICT;
 }

 let center = aabb.center();
 face_aabbs.push(aabb);
 face_centers.push(center);
 }

 // Rebuild ordered face index array (sequential, matching the collected face order).
 let ordered_indices: Vec<usize> = (0..n_faces).collect();
 let mut bvh = Bvh {
 nodes: Vec::new(),
 face_indices: ordered_indices,
 face_aabbs,
 face_centers,
 };

 if n_faces > 0 {
 bvh.build_recursive(0, n_faces);
 }

 bvh
 }

 /// Recursively build nodes; returns the new node index in `nodes`.
 fn build_recursive(&mut self, start: usize, end: usize) -> usize {
 let count = end - start;

 // Union AABB for the current face range
 let mut aabb = Aabb::empty();
 for i in start..end {
 aabb.expand_aabb(&self.face_aabbs[self.face_indices[i]]);
 }

 // Leaf when few enough faces
 if count <= MAX_LEAF_SIZE {
 let node_idx = self.nodes.len();
 self.nodes.push(BvhNode::Leaf { aabb, start, end });
 return node_idx;
 }

 // SAH: pick axis and split plane
 let (split_axis, split_pos) = self.sah_split(start, end, &aabb);

 // Partition `face_indices` in place
 let mid = self.partition(start, end, split_axis, split_pos);

 // Avoid degenerate splits (all faces on one side)
 let mid = if mid == start || mid == end {
 (start + end) / 2
 } else {
 mid
 };

 // Placeholder internal node before recursing into children
 let node_idx = self.nodes.len();
 self.nodes.push(BvhNode::Internal { aabb: Aabb::empty(), left: 0, right: 0 });

 let left = self.build_recursive(start, mid);
 let right = self.build_recursive(mid, end);

 // Fill in internal node AABB and child links
 self.nodes[node_idx] = BvhNode::Internal { aabb, left, right };

 node_idx
 }

 /// SAH split: returns `(axis 0/1/2, split coordinate)`.
 fn sah_split(&self, start: usize, end: usize, parent_aabb: &Aabb) -> (usize, f64) {
 let parent_sa = parent_aabb.surface_area().max(TOLERANCE_LEN_SQ_DIV_SAFE);
 let mut best_cost = f64::INFINITY;
 let mut best_axis = 0usize;
 let mut best_pos = 0.0f64;

 for axis in 0..3usize {
 let axis_min = match axis {
 0 => parent_aabb.min.x,
 1 => parent_aabb.min.y,
 _ => parent_aabb.min.z,
 };
 let axis_max = match axis {
 0 => parent_aabb.max.x,
 1 => parent_aabb.max.y,
 _ => parent_aabb.max.z,
 };
 let span = axis_max - axis_min;
 if span < TOLERANCE_FLOAT_LOOSE {
 continue;
 }

 for b in 1..SAH_BUCKETS {
 let split = axis_min + span * b as f64 / SAH_BUCKETS as f64;

 let mut left_aabb = Aabb::empty();
 let mut right_aabb = Aabb::empty();
 let mut left_count = 0usize;
 let mut right_count = 0usize;

 for i in start..end {
 let fi = self.face_indices[i];
 let center_val = match axis {
 0 => self.face_centers[fi].x,
 1 => self.face_centers[fi].y,
 _ => self.face_centers[fi].z,
 };
 if center_val < split {
 left_aabb.expand_aabb(&self.face_aabbs[fi]);
 left_count += 1;
 } else {
 right_aabb.expand_aabb(&self.face_aabbs[fi]);
 right_count += 1;
 }
 }

 if left_count == 0 || right_count == 0 {
 continue;
 }

 let cost = (left_count as f64 * left_aabb.surface_area()
 + right_count as f64 * right_aabb.surface_area())
 / parent_sa;

 if cost < best_cost {
 best_cost = cost;
 best_axis = axis;
 best_pos = split;
 }
 }
 }

 // Fallback: split at midpoint along the longest AABB axis
 if best_cost.is_infinite() {
 let d = parent_aabb.max - parent_aabb.min;
 best_axis = if d.x >= d.y && d.x >= d.z { 0 } else if d.y >= d.z { 1 } else { 2 };
 best_pos = parent_aabb.center()[best_axis];
 }

 (best_axis, best_pos)
 }

 /// In-place partition of `face_indices[start..end]`; returns the split index.
 fn partition(&mut self, start: usize, end: usize, axis: usize, split_pos: f64) -> usize {
 let mut mid = start;
 for i in start..end {
 let fi = self.face_indices[i];
 let center_val = match axis {
 0 => self.face_centers[fi].x,
 1 => self.face_centers[fi].y,
 _ => self.face_centers[fi].z,
 };
 if center_val < split_pos {
 self.face_indices.swap(i, mid);
 mid += 1;
 }
 }
 mid
 }

 // ──────────────────────────────────────────────────────────────────────────
 // Query API
 // ──────────────────────────────────────────────────────────────────────────

 /// Ray cast: first face hit and ray parameter `t`.
 ///
 /// `origin` is the ray origin; `dir` is the direction (need not be unit).
 pub fn ray_cast(&self, origin: DVec3, dir: DVec3) -> Option<(usize, f64)> {
 if self.nodes.is_empty() {
 return None;
 }
 let inv_dir = DVec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
 let mut best: Option<(usize, f64)> = None;
 self.ray_cast_node(0, origin, inv_dir, &mut best);
 best
 }

 fn ray_cast_node(
 &self,
 node_idx: usize,
 origin: DVec3,
 inv_dir: DVec3,
 best: &mut Option<(usize, f64)>,
 ) {
 let node = &self.nodes[node_idx];
 let t_aabb = node.aabb().ray_intersect(origin, inv_dir);
 let t_hit = match t_aabb {
 None => return,
 Some(t) => t,
 };

 // Prune if AABB entry is farther than the best hit so far
 if let Some((_, best_t)) = best
 && t_hit > *best_t {
 return;
 }

 match node {
 BvhNode::Leaf { start, end, .. } => {
 for i in *start..*end {
 let fi = self.face_indices[i];
 // Coarse test: face AABB only (exact ray–face test is caller’s job)
 if let Some(t) = self.face_aabbs[fi].ray_intersect(origin, inv_dir) {
 let update = best.is_none_or(|(_, bt)| t < bt);
 if update {
 *best = Some((fi, t));
 }
 }
 }
 }
 BvhNode::Internal { left, right, .. } => {
 self.ray_cast_node(*left, origin, inv_dir, best);
 self.ray_cast_node(*right, origin, inv_dir, best);
 }
 }
 }

 /// All face indices whose AABB intersects `query`.
 ///
 /// Used to cull face pairs for gap/overlap detection.
 pub fn query_aabb(&self, query: &Aabb) -> Vec<usize> {
 let mut result = Vec::new();
 if !self.nodes.is_empty() {
 self.query_aabb_node(0, query, &mut result);
 }
 result
 }

 fn query_aabb_node(&self, node_idx: usize, query: &Aabb, result: &mut Vec<usize>) {
 let node = &self.nodes[node_idx];
 if !node.aabb().intersects(query) {
 return;
 }
 match node {
 BvhNode::Leaf { start, end, .. } => {
 for i in *start..*end {
 let fi = self.face_indices[i];
 if self.face_aabbs[fi].intersects(query) {
 result.push(fi);
 }
 }
 }
 BvhNode::Internal { left, right, .. } => {
 self.query_aabb_node(*left, query, result);
 self.query_aabb_node(*right, query, result);
 }
 }
 }

 /// Up to `max_k` nearest faces to `point` (approximate, sorted by AABB distance).
 ///
 /// `max_dist` is the search radius (faces farther away are omitted).
 pub fn nearest_faces(&self, point: DVec3, max_dist: f64, max_k: usize) -> Vec<(usize, f64)> {
 let mut candidates: Vec<(usize, f64)> = Vec::new();
 if self.nodes.is_empty() {
 return candidates;
 }
 let max_dist_sq = max_dist * max_dist;
 self.nearest_faces_node(0, point, max_dist_sq, &mut candidates);
 candidates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
 candidates.truncate(max_k);
 candidates
 }

 fn nearest_faces_node(
 &self,
 node_idx: usize,
 point: DVec3,
 max_dist_sq: f64,
 result: &mut Vec<(usize, f64)>,
 ) {
 let node = &self.nodes[node_idx];
 let d_sq = node.aabb().point_dist_sq(point);
 if d_sq > max_dist_sq {
 return;
 }
 match node {
 BvhNode::Leaf { start, end, .. } => {
 for i in *start..*end {
 let fi = self.face_indices[i];
 let face_d_sq = self.face_aabbs[fi].point_dist_sq(point);
 if face_d_sq <= max_dist_sq {
 result.push((fi, face_d_sq.sqrt()));
 }
 }
 }
 BvhNode::Internal { left, right, .. } => {
 self.nearest_faces_node(*left, point, max_dist_sq, result);
 self.nearest_faces_node(*right, point, max_dist_sq, result);
 }
 }
 }

 /// Candidate face pairs between two BVHs whose AABBs may intersect.
 ///
 /// Culls pairs before Boolean operations instead of an O(n²) brute force.
 pub fn candidate_pairs(bvh_a: &Bvh, bvh_b: &Bvh) -> Vec<(usize, usize)> {
 let mut pairs = Vec::new();
 if bvh_a.nodes.is_empty() || bvh_b.nodes.is_empty() {
 return pairs;
 }
 Self::candidate_pairs_node(bvh_a, 0, bvh_b, 0, &mut pairs);
 pairs
 }

 fn candidate_pairs_node(
 bvh_a: &Bvh,
 node_a: usize,
 bvh_b: &Bvh,
 node_b: usize,
 pairs: &mut Vec<(usize, usize)>,
 ) {
 let na = &bvh_a.nodes[node_a];
 let nb = &bvh_b.nodes[node_b];

 if !na.aabb().intersects(nb.aabb()) {
 return;
 }

 match (na, nb) {
 (BvhNode::Leaf { start: sa, end: ea, .. }, BvhNode::Leaf { start: sb, end: eb, .. }) => {
 for ia in *sa..*ea {
 for ib in *sb..*eb {
 let fa = bvh_a.face_indices[ia];
 let fb = bvh_b.face_indices[ib];
 if bvh_a.face_aabbs[fa].intersects(&bvh_b.face_aabbs[fb]) {
 pairs.push((fa, fb));
 }
 }
 }
 }
 (BvhNode::Internal { left: la, right: ra, .. }, _) => {
 Self::candidate_pairs_node(bvh_a, *la, bvh_b, node_b, pairs);
 Self::candidate_pairs_node(bvh_a, *ra, bvh_b, node_b, pairs);
 }
 (_, BvhNode::Internal { left: lb, right: rb, .. }) => {
 Self::candidate_pairs_node(bvh_a, node_a, bvh_b, *lb, pairs);
 Self::candidate_pairs_node(bvh_a, node_a, bvh_b, *rb, pairs);
 }
 }
 }

 /// Debug / profiling statistics for this BVH.
 pub fn stats(&self) -> BvhStats {
 let mut stats = BvhStats::default();
 if !self.nodes.is_empty() {
 self.stats_node(0, 0, &mut stats);
 }
 stats
 }

 fn stats_node(&self, node_idx: usize, depth: usize, stats: &mut BvhStats) {
 stats.node_count += 1;
 stats.max_depth = stats.max_depth.max(depth);
 match &self.nodes[node_idx] {
 BvhNode::Leaf { start, end, .. } => {
 stats.leaf_count += 1;
 stats.total_leaf_faces += end - start;
 stats.max_leaf_faces = stats.max_leaf_faces.max(end - start);
 }
 BvhNode::Internal { left, right, .. } => {
 self.stats_node(*left, depth + 1, stats);
 self.stats_node(*right, depth + 1, stats);
 }
 }
 }
}

/// Aggregated BVH statistics.
#[derive(Debug, Default)]
pub struct BvhStats {
 pub node_count: usize,
 pub leaf_count: usize,
 pub max_depth: usize,
 pub total_leaf_faces: usize,
 pub max_leaf_faces: usize,
}

impl BvhStats {
 pub fn avg_leaf_faces(&self) -> f64 {
 if self.leaf_count == 0 {
 0.0
 } else {
 self.total_leaf_faces as f64 / self.leaf_count as f64
 }
 }
}

// ============================================================================
// DS-level BVH for spatial acceleration in PaveFiller (VE, EE, VF, EF)
// ============================================================================

/// A simplified BVH over DS entities (vertices or edges) for pair culling.
///
/// OCCT `BOPDS_Iterator` builds BVH trees over all DS sub-shapes to
/// avoid O(n²) pair enumeration.  `DsBvh` provides the same culling
/// BVH for DS entity pair filtering.
/// Uses median-split builder (same algorithm as OCCT BVH_SpatialMedianBuilder).
pub struct DsBvh {
 nodes: Vec<BvhNode>,
 indices: Vec<usize>,
 aabbs: Vec<Aabb>,
}

impl DsBvh {
 /// Build a BVH from a list of entity AABBs.
 /// `indices` and `aabbs` must have the same length.
 pub fn build(indices: Vec<usize>, aabbs: Vec<Aabb>) -> Self {
 let n = indices.len();
 let mut order: Vec<usize> = (0..n).collect();
 let mut nodes = Vec::new();
 let mut sorted_indices = Vec::with_capacity(n);
 let mut sorted_aabbs = Vec::with_capacity(n);

 if n > 0 {
 Self::build_rec(&mut order, &aabbs, 0, n, &mut nodes,
 &mut sorted_indices, &mut sorted_aabbs, 0);
 }

 Self { nodes, indices: sorted_indices, aabbs: sorted_aabbs }
 }

 // BVH_SpatialMedianBuilder — median-split construction.
 fn build_rec(
 order: &mut [usize],
 aabbs: &[Aabb],
 start: usize,
 end: usize,
 nodes: &mut Vec<BvhNode>,
 out_indices: &mut Vec<usize>,
 out_aabbs: &mut Vec<Aabb>,
 depth: usize,
 ) -> usize {
 let mut node_aabb = Aabb::empty();
 for &oi in &order[start..end] {
 node_aabb.expand_aabb(&aabbs[oi]);
 }

 let count = end - start;
 if count <= 4 {
 let idx = nodes.len();
 let leaf_start = out_indices.len();
 for &oi in &order[start..end] {
 out_indices.push(oi);
 out_aabbs.push(aabbs[oi]);
 }
 nodes.push(BvhNode::Leaf {
 aabb: node_aabb,
 start: leaf_start,
 end: out_indices.len(),
 });
 return idx;
 }

 // Median split along the largest axis
 let axis = {
 let size = node_aabb.max - node_aabb.min;
 if size.x >= size.y && size.x >= size.z { 0 }
 else if size.y >= size.z { 1 }
 else { 2 }
 };

 order[start..end].sort_by(|&a, &b| {
 let ca = &aabbs[a];
 let cb = &aabbs[b];
 let va = [ca.min.x + ca.max.x, ca.min.y + ca.max.y, ca.min.z + ca.max.z][axis];
 let vb = [cb.min.x + cb.max.x, cb.min.y + cb.max.y, cb.min.z + cb.max.z][axis];
 va.partial_cmp(&vb).unwrap_or(std::cmp::Ordering::Equal)
 });

 let mid = start + count / 2;
 let left = Self::build_rec(order, aabbs, start, mid, nodes, out_indices, out_aabbs, depth + 1);
 let right = Self::build_rec(order, aabbs, mid, end, nodes, out_indices, out_aabbs, depth + 1);

 let idx = nodes.len();
 nodes.push(BvhNode::Internal { aabb: node_aabb, left, right });
 idx
 }

 /// Candidate entity pairs between two DS BVHs whose AABBs overlap.
 /// dual-tree BVH traversal (matching IntPatch_BVHTraversal).
 pub fn candidate_pairs(bvh_a: &DsBvh, bvh_b: &DsBvh) -> Vec<(usize, usize)> {
 let mut pairs = Vec::new();
 if bvh_a.nodes.is_empty() || bvh_b.nodes.is_empty() {
 return pairs;
 }
 Self::candidate_pairs_node(bvh_a, 0, bvh_b, 0, &mut pairs);
 pairs
 }

 fn candidate_pairs_node(
 bvh_a: &DsBvh, node_a: usize,
 bvh_b: &DsBvh, node_b: usize,
 pairs: &mut Vec<(usize, usize)>,
 ) {
 let na = &bvh_a.nodes[node_a];
 let nb = &bvh_b.nodes[node_b];

 if !na.aabb().intersects(nb.aabb()) {
 return;
 }

 match (na, nb) {
 (BvhNode::Leaf { start: sa, end: ea, .. }, BvhNode::Leaf { start: sb, end: eb, .. }) => {
 for ia in *sa..*ea {
 for ib in *sb..*eb {
 let ia_idx = bvh_a.indices[ia];
 let ib_idx = bvh_b.indices[ib];
 if bvh_a.aabbs[ia].intersects(&bvh_b.aabbs[ib]) {
 pairs.push((ia_idx, ib_idx));
 }
 }
 }
 }
 (BvhNode::Internal { left: la, right: ra, .. }, _) => {
 Self::candidate_pairs_node(bvh_a, *la, bvh_b, node_b, pairs);
 Self::candidate_pairs_node(bvh_a, *ra, bvh_b, node_b, pairs);
 }
 (_, BvhNode::Internal { left: lb, right: rb, .. }) => {
 Self::candidate_pairs_node(bvh_a, node_a, bvh_b, *lb, pairs);
 Self::candidate_pairs_node(bvh_a, node_a, bvh_b, *rb, pairs);
 }
 }
 }

 /// Query all items whose AABB overlaps the query AABB.
 pub fn query_aabb(&self, query: &Aabb) -> Vec<usize> {
 if self.nodes.is_empty() { return vec![]; }
 let mut results = Vec::new();
 self.query_aabb_node(0, query, &mut results);
 results
 }

 fn query_aabb_node(&self, node_idx: usize, query: &Aabb, results: &mut Vec<usize>) {
 let node = &self.nodes[node_idx];
 if !node.aabb().intersects(query) { return; }
 match node {
 BvhNode::Leaf { start, end, .. } => {
 for i in *start..*end {
 if self.aabbs[i].intersects(query) {
 results.push(self.indices[i]);
 }
 }
 }
 BvhNode::Internal { left, right, .. } => {
 self.query_aabb_node(*left, query, results);
 self.query_aabb_node(*right, query, results);
 }
 }
 }
}


