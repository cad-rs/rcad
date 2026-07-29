// OCCT IntTools_EdgeEdge — edge-edge intersection.
// OCCT IntTools_EdgeEdge.cxx / .lxx
use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};

#[derive(Debug, Clone)]
pub struct CommonPrt {
    pub is_edge: bool,
    pub range1: [f64; 2],
    pub ranges2: Vec<[f64; 2]>,
    pub vertex_param1: f64,
    pub vertex_param2: f64,
    pub bounding_point1: DVec3,
    pub bounding_point2: DVec3,
}

pub struct EdgeEdgeIntersector {
    curve1: Curve3, curve2: Curve3,
    range1: [f64; 2], range2: [f64; 2],
    common_parts: Vec<CommonPrt>, done: bool,
    fuzzy_value: f64,
    quick_coincidence_check: bool,
}
impl EdgeEdgeIntersector {
    pub fn new() -> Self {
        EdgeEdgeIntersector {
            curve1: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            curve2: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            range1: [0.0, 1.0], range2: [0.0, 1.0],
            common_parts: Vec::new(), done: false,
            fuzzy_value: 1e-7,
            quick_coincidence_check: false,
        }
    }
    pub fn use_quick_coincidence_check(&mut self, b: bool) { self.quick_coincidence_check = b; }
    pub fn set_edges(&mut self, ei1: usize, r1: [f64; 2], ei2: usize, r2: [f64; 2], ds: &crate::bop::ds::DS) -> &mut Self {
        if let Some(c) = ds.edge_curve(ei1) { self.curve1 = c.clone(); }
        if let Some(c) = ds.edge_curve(ei2) { self.curve2 = c.clone(); }
        self.range1 = r1; self.range2 = r2; self
    }
    pub fn set_fuzzy_value(&mut self, f: f64) { self.fuzzy_value = f.max(1e-7); }

    /// OCCT IntTools_EdgeEdge::Perform (IntTools_EdgeEdge.cxx L185-243).
    pub fn perform(&mut self) {
        self.common_parts.clear();
        // OCCT L198-201: Line-Line analytic
        if matches!(&self.curve1, Curve3::Line(_)) && matches!(&self.curve2, Curve3::Line(_)) {
            self.compute_line_line(); self.done = true; return;
        }
        // OCCT L204-215: Quick coincident check
        if self.quick_coincidence_check && self.is_coincident() {
            self.common_parts.push(CommonPrt {
                is_edge: true,
                range1: self.range1, ranges2: vec![self.range2],
                vertex_param1: (self.range1[0] + self.range1[1]) * 0.5,
                vertex_param2: (self.range2[0] + self.range2[1]) * 0.5,
                bounding_point1: self.curve1.point_at((self.range1[0] + self.range1[1]) * 0.5),
                bounding_point2: self.curve2.point_at((self.range2[0] + self.range2[1]) * 0.5),
            });
            self.done = true; return;
        }
        // OCCT L237-242: FindSolutions + MergeSolutions
        self.find_and_merge_solutions();
        self.done = true;
    }

    /// OCCT IntTools_EdgeEdge::IsCoincident (L247-285): 24pt sampling.
    fn is_coincident(&self) -> bool {
        let a_nb_seg = 23usize;
        let t11 = self.range1[0]; let t12 = self.range1[1];
        let t21 = self.range2[0]; let t22 = self.range2[1];
        let dt = (t12 - t11) / a_nb_seg as f64;
        let mut i_cnt = 0;
        for i in 0..=a_nb_seg {
            let a_t1 = t11 + i as f64 * dt;
            let p1 = self.curve1.point_at(a_t1);
            let mut best_d = f64::MAX;
            for j in 0..=a_nb_seg {
                let a_t2 = t21 + (t22 - t21) * j as f64 / a_nb_seg as f64;
                let d = (p1 - self.curve2.point_at(a_t2)).length();
                if d < best_d { best_d = d; }
            }
            if best_d < self.fuzzy_value { i_cnt += 1; }
        }
        i_cnt as f64 / (a_nb_seg + 1) as f64 > 0.5
    }

    /// OCCT FindSolutions + MergeSolutions: adaptive sampling with local refinement.
    fn find_and_merge_solutions(&mut self) {
        let t1 = self.range1[0]; let t2 = self.range1[1];
        let t3 = self.range2[0]; let t4 = self.range2[1];

        // Coarse pass: find candidate regions
        let coarse = 23usize;
        let dt1 = (t2 - t1) / coarse as f64;
        let dt2 = (t4 - t3) / coarse as f64;
        let mut candidates: Vec<(f64, f64, f64)> = Vec::new(); // (t1, t2, distance)
        for i in 0..=coarse {
            let a_t1 = t1 + i as f64 * dt1;
            let p1 = self.curve1.point_at(a_t1);
            let mut best_d = f64::MAX; let mut best_t2 = t3;
            for j in 0..=coarse {
                let a_t2 = t3 + j as f64 * dt2;
                let d = (p1 - self.curve2.point_at(a_t2)).length();
                if d < best_d { best_d = d; best_t2 = a_t2; }
            }
            if best_d < self.fuzzy_value * 10.0 {
                candidates.push((a_t1, best_t2, best_d));
            }
        }
        if candidates.is_empty() { return; }

        // Local refinement: sub-sample around each candidate
        let fine = 10usize;
        for (base_t1, base_t2, _base_d) in &candidates {
            let mut best_d = f64::MAX; let mut best_t1 = *base_t1; let mut best_t2 = *base_t2;
            let refine_range = dt1.max(dt2) * 0.5;
            for i in 0..=fine {
                let a_t1 = (*base_t1 - refine_range) + 2.0 * refine_range * i as f64 / fine as f64;
                if a_t1 < t1 || a_t1 > t2 { continue; }
                let p1 = self.curve1.point_at(a_t1);
                for j in 0..=fine {
                    let a_t2 = (*base_t2 - refine_range) + 2.0 * refine_range * j as f64 / fine as f64;
                    if a_t2 < t3 || a_t2 > t4 { continue; }
                    let d = (p1 - self.curve2.point_at(a_t2)).length();
                    if d < best_d { best_d = d; best_t1 = a_t1; best_t2 = a_t2; }
                }
            }
            if best_d <= self.fuzzy_value {
                self.common_parts.push(CommonPrt {
                    is_edge: false,
                    range1: [best_t1, best_t1], ranges2: vec![[best_t2, best_t2]],
                    vertex_param1: best_t1, vertex_param2: best_t2,
                    bounding_point1: self.curve1.point_at(best_t1),
                    bounding_point2: self.curve2.point_at(best_t2),
                });
            }
        }
    }

    /// OCCT IntTools_EdgeEdge::ComputeLineLine (L300+).
    fn compute_line_line(&mut self) {
        let (l1, l2) = match (&self.curve1, &self.curve2) {
            (Curve3::Line(l1), Curve3::Line(l2)) => (l1, l2),
            _ => return,
        };
        let o1 = l1.origin; let d1 = l1.direction;
        let o2 = l2.origin; let d2 = l2.direction;
        let r = o2 - o1;
        let d1xd2 = d1.cross(d2);
        let denom = d1xd2.length_squared();
        if denom < 1e-30 {
            let dist = r.cross(d1).length() / d1.length().max(1e-30);
            if dist <= self.fuzzy_value {
                let s1 = self.range2[0]; let s2 = self.range2[1];
                let d1_sq = d1.length_squared().max(1e-30);
                let proj_s1 = (r.dot(d1) + s1 * d2.dot(d1)) / d1_sq;
                let proj_s2 = (r.dot(d1) + s2 * d2.dot(d1)) / d1_sq;
                let ov_start = self.range1[0].max(proj_s1.min(proj_s2));
                let ov_end = self.range1[1].min(proj_s1.max(proj_s2));
                if ov_end > ov_start + 1e-12 {
                    self.common_parts.push(CommonPrt {
                        is_edge: true, range1: [ov_start, ov_end],
                        ranges2: vec![[s1, s2]],
                        vertex_param1: ov_start, vertex_param2: s1,
                        bounding_point1: self.curve1.point_at(ov_start),
                        bounding_point2: self.curve2.point_at(s1),
                    });
                }
            }
            return;
        }
        let t = r.cross(d2).dot(d1xd2) / denom;
        let s = r.cross(d1).dot(d1xd2) / denom;
        if t >= self.range1[0] - 1e-12 && t <= self.range1[1] + 1e-12
            && s >= self.range2[0] - 1e-12 && s <= self.range2[1] + 1e-12 {
            let p = o1 + d1 * t;
            self.common_parts.push(CommonPrt {
                is_edge: false, range1: [t, t], ranges2: vec![[s, s]],
                vertex_param1: t, vertex_param2: s,
                bounding_point1: p, bounding_point2: p,
            });
        }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn common_parts(&self) -> &[CommonPrt] { &self.common_parts }
}
