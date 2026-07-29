// OCCT IntTools_EdgeEdge ?edge-edge intersection.
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
}
impl EdgeEdgeIntersector {
    pub fn new() -> Self {
        EdgeEdgeIntersector {
            curve1: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            curve2: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            range1: [0.0, 1.0], range2: [0.0, 1.0],
            common_parts: Vec::new(), done: false,
            fuzzy_value: 1e-7,
        }
    }
    pub fn set_edges(&mut self, ei1: usize, r1: [f64; 2], ei2: usize, r2: [f64; 2], ds: &crate::bop::ds::DS) -> &mut Self {
        if let Some(c) = ds.edge_curve(ei1) { self.curve1 = c.clone(); }
        if let Some(c) = ds.edge_curve(ei2) { self.curve2 = c.clone(); }
        self.range1 = r1; self.range2 = r2; self
    }
    pub fn set_fuzzy_value(&mut self, f: f64) {
        self.fuzzy_value = f;
    }
    pub fn perform(&mut self) {
        self.common_parts.clear();
        let n_samples = 64usize;
        let pts1: Vec<(f64, DVec3)> = (0..=n_samples).map(|i| {
            let t = self.range1[0] + (self.range1[1] - self.range1[0]) * i as f64 / n_samples as f64;
            (t, self.curve1.point_at(t))
        }).collect();
        let pts2: Vec<(f64, DVec3)> = (0..=n_samples).map(|i| {
            let t = self.range2[0] + (self.range2[1] - self.range2[0]) * i as f64 / n_samples as f64;
            (t, self.curve2.point_at(t))
        }).collect();
        let mut best_dist = f64::MAX;
        let mut best = (0usize, 0usize);
        for (i, (_, p1)) in pts1.iter().enumerate() {
            for (j, (_, p2)) in pts2.iter().enumerate() {
                let d = (p1 - p2).length();
                if d < best_dist { best_dist = d; best = (i, j); }
            }
        }
        if best_dist <= 1e-7 {
            let (i, j) = best;
            let (t1, p1) = pts1[i];
            let (t2, p2) = pts2[j];
            self.common_parts.push(CommonPrt {
                is_edge: false, range1: [t1, t1], ranges2: vec![[t2, t2]],
                vertex_param1: t1, vertex_param2: t2,
                bounding_point1: p1, bounding_point2: p2,
            });
        }
        self.done = true;
    }
    pub fn is_done(&self) -> bool { self.done }
    pub fn common_parts(&self) -> &[CommonPrt] { &self.common_parts }
}
