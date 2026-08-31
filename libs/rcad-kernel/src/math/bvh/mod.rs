//! OCCT BVH: Bounding Volume Hierarchy for spatial acceleration.
//!
//! Analogous to OCCT `BVH_Tree` / `BVH_Builder`.
//!
//! Built with the SAH (Surface Area Heuristic) to speed up:
//! - Face pair culling for Boolean operations
//! - Ray picking
//! - Nearest-face queries

pub mod bounding_box;

pub use bounding_box::Aabb;

use glam::DVec3;

use crate::geom::{Surface3, SurfaceEval};
use crate::math::bnd::BndBox;
use crate::topods::{self, TShape};

const TOL_LINEAR_ULTRA_STRICT: f64 = 1e-10;
const TOL_FLOAT_LOOSE: f64 = 1e-14;
const TOL_LEN_SQ_DIV_SAFE: f64 = 1e-30;

/// BVH node (internal or leaf).
#[derive(Debug, Clone)]
pub(crate) enum BvhNode {
    Leaf {
        aabb: BndBox,
        start: usize,
        end: usize,
    },
    Internal {
        aabb: BndBox,
        left: usize,
        right: usize,
    },
}

impl BvhNode {
    pub(crate) fn aabb(&self) -> &BndBox {
        match self {
            BvhNode::Leaf { aabb, .. } => aabb,
            BvhNode::Internal { aabb, .. } => aabb,
        }
    }
}

/// BVH over the faces of one BRep.
pub struct Bvh {
    nodes: Vec<BvhNode>,
    face_indices: Vec<usize>,
    face_aabbs: Vec<BndBox>,
    face_centers: Vec<DVec3>,
}

const MAX_LEAF_SIZE: usize = 4;
const SAH_BUCKETS: usize = 8;

impl Bvh {
    pub fn face_count(&self) -> usize { self.face_aabbs.len() }

    pub fn face_aabb(&self, face_index: usize) -> Option<&BndBox> {
        self.face_aabbs.get(face_index)
    }

    /// Build a BVH over all faces of `brep`.
    pub fn build(brep: &crate::BRep) -> Self {
        let face_ts_indices: Vec<usize> = brep.tshapes.iter().enumerate()
            .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_)))
            .map(|(fi, _)| fi).collect();
        let n_faces = face_ts_indices.len();

        let mut face_aabbs = Vec::with_capacity(n_faces);
        let mut face_centers = Vec::with_capacity(n_faces);

        for &fi in &face_ts_indices {
            let ts = &brep.tshapes[fi];
            let TShape::Face(fd) = ts.as_ref() else { continue; };
            let mut aabb = BndBox::new();

            if let Some(wts) = brep.tshapes.get(fd.outer_wire.index) {
                if let TShape::Wire(wd) = wts.as_ref() {
                    for er in &wd.edges {
                        if let Some(ets) = brep.tshapes.get(er.index) {
                            if let TShape::Edge(ed) = ets.as_ref() {
                                if let Some(p0) = brep.vertex_point(ed.first.index) {
                                    aabb.add_point(p0);
                                }
                                if let Some(p1) = brep.vertex_point(ed.last.index) {
                                    aabb.add_point(p1);
                                }
                            }
                        }
                    }
                }
            }

            if let Some(surf) = &fd.surface {
                match surf {
                    Surface3::Sphere(s) => {
                        let r = s.radius.abs() + TOL_LINEAR_ULTRA_STRICT;
                        aabb.add_point(s.center - DVec3::splat(r));
                        aabb.add_point(s.center + DVec3::splat(r));
                    }
                    Surface3::Cylinder(c) => {
                        let [_, _, v0, v1] = surf.default_domain();
                        let ax = c.axis.normalize_or_zero();
                        let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                        let u_dir = ax.cross(perp).normalize_or_zero();
                        let v_dir = ax.cross(u_dir).normalize_or_zero();
                        let r = c.radius.abs() + TOL_LINEAR_ULTRA_STRICT;
                        for &vh in &[v0, v1] {
                            for k in 0..8 {
                                let a = std::f64::consts::TAU * k as f64 / 8.0;
                                let p = c.origin + ax * vh + u_dir * r * a.cos() + v_dir * r * a.sin();
                                aabb.add_point(p);
                            }
                        }
                    }
                    Surface3::Cone(c) => {
                        let [_, _, v0, v1] = surf.default_domain();
                        let ax = c.axis.normalize_or_zero();
                        let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                        let u_dir = ax.cross(perp).normalize_or_zero();
                        let v_dir = ax.cross(u_dir).normalize_or_zero();
                        for &vh in &[v0, v1] {
                            let r_at = (c.radius + vh * c.half_angle_rad.tan()).abs() + TOL_LINEAR_ULTRA_STRICT;
                            let center = c.apex + ax * vh;
                            for k in 0..8 {
                                let a = std::f64::consts::TAU * k as f64 / 8.0;
                                let p = center + u_dir * r_at * a.cos() + v_dir * r_at * a.sin();
                                aabb.add_point(p);
                            }
                        }
                    }
                    Surface3::Torus(t) => {
                        let r_out = t.major_radius.abs() + t.minor_radius.abs() + TOL_LINEAR_ULTRA_STRICT;
                        let ax = t.axis.normalize_or_zero();
                        let perp = if ax.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                        let u_dir = ax.cross(perp).normalize_or_zero();
                        let v_dir = ax.cross(u_dir).normalize_or_zero();
                        for k in 0..8 {
                            let a = std::f64::consts::TAU * k as f64 / 8.0;
                            let c = t.center + u_dir * t.major_radius * a.cos() + v_dir * t.major_radius * a.sin();
                            aabb.add_point(c + ax * t.minor_radius);
                            aabb.add_point(c - ax * t.minor_radius);
                        }
                    }
                    _ => {
                        let [u0, u1, v0, v1] = surf.default_domain();
                        for i in 0..=2 {
                            for j in 0..=2 {
                                let u = u0 + (u1 - u0) * i as f64 / 2.0;
                                let v = v0 + (v1 - v0) * j as f64 / 2.0;
                                let p = surf.point_at(u, v);
                                if p.is_finite() { aabb.add_point(p); }
                            }
                        }
                    }
                }
            }

            // Degenerate faces: nudge AABB to non-zero extent
            if aabb.dx() < TOL_LINEAR_ULTRA_STRICT {
                aabb.add_point(aabb.raw_min() - DVec3::splat(TOL_LINEAR_ULTRA_STRICT));
                aabb.add_point(aabb.raw_max() + DVec3::splat(TOL_LINEAR_ULTRA_STRICT));
            }

            let center = aabb.center();
            face_aabbs.push(aabb);
            face_centers.push(center);
        }

        let ordered_indices: Vec<usize> = (0..n_faces).collect();
        let mut bvh = Bvh { nodes: Vec::new(), face_indices: ordered_indices, face_aabbs, face_centers };
        if n_faces > 0 { bvh.build_recursive(0, n_faces); }
        bvh
    }

    fn build_recursive(&mut self, start: usize, end: usize) -> usize {
        let count = end - start;
        let mut aabb = BndBox::new();
        for i in start..end { aabb.add_box(&self.face_aabbs[self.face_indices[i]]); }

        if count <= MAX_LEAF_SIZE {
            let node_idx = self.nodes.len();
            self.nodes.push(BvhNode::Leaf { aabb, start, end });
            return node_idx;
        }

        let (split_axis, split_pos) = self.sah_split(start, end, &aabb);
        let mid = self.partition(start, end, split_axis, split_pos);
        let mid = if mid == start || mid == end { (start + end) / 2 } else { mid };

        let node_idx = self.nodes.len();
        self.nodes.push(BvhNode::Internal { aabb: BndBox::new(), left: 0, right: 0 });

        let left = self.build_recursive(start, mid);
        let right = self.build_recursive(mid, end);
        self.nodes[node_idx] = BvhNode::Internal { aabb, left, right };

        node_idx
    }

    fn sah_split(&self, start: usize, end: usize, parent_aabb: &BndBox) -> (usize, f64) {
        let parent_sa = parent_aabb.surface_area().max(TOL_LEN_SQ_DIV_SAFE);
        let mut best_cost = f64::INFINITY;
        let mut best_axis = 0usize;
        let mut best_pos = 0.0;

        let pmin = parent_aabb.raw_min();
        let pmax = parent_aabb.raw_max();

        for axis in 0..3usize {
            let axis_min = pmin[axis];
            let axis_max = pmax[axis];
            let span = axis_max - axis_min;
            if span < TOL_FLOAT_LOOSE { continue; }

            for b in 1..SAH_BUCKETS {
                let split = axis_min + span * b as f64 / SAH_BUCKETS as f64;
                let mut left_aabb = BndBox::new();
                let mut right_aabb = BndBox::new();
                let mut left_count = 0usize;
                let mut right_count = 0usize;

                for i in start..end {
                    let fi = self.face_indices[i];
                    let cv = self.face_centers[fi][axis];
                    if cv < split { left_aabb.add_box(&self.face_aabbs[fi]); left_count += 1; }
                    else { right_aabb.add_box(&self.face_aabbs[fi]); right_count += 1; }
                }
                if left_count == 0 || right_count == 0 { continue; }

                let cost = (left_count as f64 * left_aabb.surface_area()
                    + right_count as f64 * right_aabb.surface_area()) / parent_sa;
                if cost < best_cost { best_cost = cost; best_axis = axis; best_pos = split; }
            }
        }

        if best_cost.is_infinite() {
            let d = pmax - pmin;
            best_axis = if d.x >= d.y && d.x >= d.z { 0 } else if d.y >= d.z { 1 } else { 2 };
            best_pos = parent_aabb.center()[best_axis];
        }
        (best_axis, best_pos)
    }

    fn partition(&mut self, start: usize, end: usize, axis: usize, split_pos: f64) -> usize {
        let mut mid = start;
        for i in start..end {
            let fi = self.face_indices[i];
            if self.face_centers[fi][axis] < split_pos { self.face_indices.swap(i, mid); mid += 1; }
        }
        mid
    }

    // ── Query API ─────────────────────────────────────────────────────────

    /// Ray cast: first face hit and ray parameter `t`.
    pub fn ray_cast(&self, origin: DVec3, dir: DVec3) -> Option<(usize, f64)> {
        if self.nodes.is_empty() { return None; }
        let inv_dir = DVec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
        let mut best: Option<(usize, f64)> = None;
        self.ray_cast_node(0, origin, inv_dir, &mut best);
        best
    }

    fn ray_cast_node(&self, node_idx: usize, origin: DVec3, inv_dir: DVec3, best: &mut Option<(usize, f64)>) {
        let node = &self.nodes[node_idx];
        let t_hit = match node.aabb().ray_intersect(origin, inv_dir) { None => return, Some(t) => t };
        if let Some((_, bt)) = best && t_hit > *bt { return; }
        match node {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    let fi = self.face_indices[i];
                    if let Some(t) = self.face_aabbs[fi].ray_intersect(origin, inv_dir) {
                        if best.is_none_or(|(_, bt)| t < bt) { *best = Some((fi, t)); }
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
    pub fn query_aabb(&self, query: &BndBox) -> Vec<usize> {
        let mut result = Vec::new();
        if !self.nodes.is_empty() { self.query_aabb_node(0, query, &mut result); }
        result
    }

    fn query_aabb_node(&self, node_idx: usize, query: &BndBox, result: &mut Vec<usize>) {
        let node = &self.nodes[node_idx];
        if !node.aabb().intersects(query) { return; }
        match node {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    let fi = self.face_indices[i];
                    if self.face_aabbs[fi].intersects(query) { result.push(fi); }
                }
            }
            BvhNode::Internal { left, right, .. } => {
                self.query_aabb_node(*left, query, result);
                self.query_aabb_node(*right, query, result);
            }
        }
    }

    /// Up to `max_k` nearest faces to `point`.
    pub fn nearest_faces(&self, point: DVec3, max_dist: f64, max_k: usize) -> Vec<(usize, f64)> {
        let mut candidates = Vec::new();
        if self.nodes.is_empty() { return candidates; }
        let max_dist_sq = max_dist * max_dist;
        self.nearest_faces_node(0, point, max_dist_sq, &mut candidates);
        candidates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        candidates.truncate(max_k);
        candidates
    }

    fn nearest_faces_node(&self, node_idx: usize, point: DVec3, max_dist_sq: f64, result: &mut Vec<(usize, f64)>) {
        let node = &self.nodes[node_idx];
        if node.aabb().point_dist_sq(point) > max_dist_sq { return; }
        match node {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    let fi = self.face_indices[i];
                    let d_sq = self.face_aabbs[fi].point_dist_sq(point);
                    if d_sq <= max_dist_sq { result.push((fi, d_sq.sqrt())); }
                }
            }
            BvhNode::Internal { left, right, .. } => {
                self.nearest_faces_node(*left, point, max_dist_sq, result);
                self.nearest_faces_node(*right, point, max_dist_sq, result);
            }
        }
    }

    /// Candidate face pairs between two BVHs.
    pub fn candidate_pairs(bvh_a: &Bvh, bvh_b: &Bvh) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if bvh_a.nodes.is_empty() || bvh_b.nodes.is_empty() { return pairs; }
        Self::candidate_pairs_node(bvh_a, 0, bvh_b, 0, &mut pairs);
        pairs
    }

    fn candidate_pairs_node(bvh_a: &Bvh, node_a: usize, bvh_b: &Bvh, node_b: usize, pairs: &mut Vec<(usize, usize)>) {
        let na = &bvh_a.nodes[node_a];
        let nb = &bvh_b.nodes[node_b];
        if !na.aabb().intersects(nb.aabb()) { return; }
        match (na, nb) {
            (BvhNode::Leaf { start: sa, end: ea, .. }, BvhNode::Leaf { start: sb, end: eb, .. }) => {
                for ia in *sa..*ea {
                    for ib in *sb..*eb {
                        let fa = bvh_a.face_indices[ia];
                        let fb = bvh_b.face_indices[ib];
                        if bvh_a.face_aabbs[fa].intersects(&bvh_b.face_aabbs[fb]) { pairs.push((fa, fb)); }
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

    pub fn stats(&self) -> BvhStats {
        let mut stats = BvhStats::default();
        if !self.nodes.is_empty() { self.stats_node(0, 0, &mut stats); }
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
        if self.leaf_count == 0 { 0.0 } else { self.total_leaf_faces as f64 / self.leaf_count as f64 }
    }
}
