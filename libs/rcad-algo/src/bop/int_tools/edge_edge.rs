// OCCT IntTools_EdgeEdge — edge-edge intersection.
//
// Finds intersection points or overlapping segments between two 3D curves.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};

/// Common part between two edges (vertex touch or edge overlap).
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

/// IntTools_EdgeEdge — edge-edge intersection engine.
pub struct EdgeEdge {
    curve1: Curve3,
    curve2: Curve3,
    range1: [f64; 2],
    range2: [f64; 2],
    tol1: f64,
    tol2: f64,
    fuzzy: f64,
    common_parts: Vec<CommonPrt>,
    done: bool,
}

impl EdgeEdge {
    pub fn new() -> Self {
        EdgeEdge {
            curve1: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            curve2: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            range1: [0.0, 1.0],
            range2: [0.0, 1.0],
            tol1: 1e-7, tol2: 1e-7,
            fuzzy: 0.0,
            common_parts: Vec::new(),
            done: false,
        }
    }

    pub fn set_curves(&mut self, c1: Curve3, c2: Curve3) { self.curve1 = c1; self.curve2 = c2; }
    pub fn set_ranges(&mut self, r1: [f64; 2], r2: [f64; 2]) { self.range1 = r1; self.range2 = r2; }
    pub fn set_tolerances(&mut self, t1: f64, t2: f64) { self.tol1 = t1; self.tol2 = t2; }
    pub fn set_fuzzy(&mut self, f: f64) { self.fuzzy = f; }
    pub fn is_done(&self) -> bool { self.done }
    pub fn common_parts(&self) -> &[CommonPrt] { &self.common_parts }

    /// Perform intersection.
    pub fn perform(&mut self) {
        self.common_parts.clear();
        let res1 = self.curve_resolution(&self.curve1, self.range1);
        let res2 = self.curve_resolution(&self.curve2, self.range2);
        let tol = self.tol1.max(self.tol2) + self.fuzzy;
        let n_samples = 64usize;

        // Sample both curves and find close point pairs
        let pts1: Vec<(f64, DVec3)> = (0..=n_samples).map(|i| {
            let t = self.range1[0] + (self.range1[1] - self.range1[0]) * i as f64 / n_samples as f64;
            (t, self.curve1.point_at(t))
        }).collect();
        let pts2: Vec<(f64, DVec3)> = (0..=n_samples).map(|i| {
            let t = self.range2[0] + (self.range2[1] - self.range2[0]) * i as f64 / n_samples as f64;
            (t, self.curve2.point_at(t))
        }).collect();

        // Find closest pairs
        let mut best_dist = f64::MAX;
        let mut best_pair = (0usize, 0usize);
        for (i, (t1, p1)) in pts1.iter().enumerate() {
            for (j, (t2, p2)) in pts2.iter().enumerate() {
                let d = (p1 - p2).length();
                if d < best_dist {
                    best_dist = d;
                    best_pair = (i, j);
                }
            }
        }

        // If close enough, record as vertex-type common part
        if best_dist <= tol {
            let (i, j) = best_pair;
            let (t1, p1) = pts1[i];
            let (t2, p2) = pts2[j];
            let res = res1.max(res2).max(tol);
            let t1_start = (t1 - res).max(self.range1[0]);
            let t1_end = (t1 + res).min(self.range1[1]);
            let t2_start = (t2 - res).max(self.range2[0]);
            let t2_end = (t2 + res).min(self.range2[1]);
            // Check if overlapping (edge-type common part)
            let overlap_len1 = t1_end - t1_start;
            let overlap_len2 = t2_end - t2_start;
            if overlap_len1 > res * 2.0 && overlap_len2 > res * 2.0 {
                self.common_parts.push(CommonPrt {
                    is_edge: true,
                    range1: [t1_start, t1_end],
                    ranges2: vec![[t2_start, t2_end]],
                    vertex_param1: t1, vertex_param2: t2,
                    bounding_point1: p1, bounding_point2: p2,
                });
            } else {
                self.common_parts.push(CommonPrt {
                    is_edge: false,
                    range1: [t1, t1],
                    ranges2: vec![[t2, t2]],
                    vertex_param1: t1, vertex_param2: t2,
                    bounding_point1: p1, bounding_point2: p2,
                });
            }
        }
        self.done = true;
    }

    fn curve_resolution(&self, curve: &Curve3, range: [f64; 2]) -> f64 {
        let n = 10usize;
        let dt = (range[1] - range[0]) / n as f64;
        let mut max_step = 0.0;
        let mut t = range[0];
        let mut prev = curve.point_at(t);
        for _ in 0..n {
            t += dt;
            let p = curve.point_at(t);
            let d = (p - prev).length();
            if d > max_step { max_step = d; }
            prev = p;
        }
        max_step
    }
}

impl Default for EdgeEdge { fn default() -> Self { Self::new() } }

/// Check if a curve type is Line (for fast path).
pub fn is_line_curve(curve: &Curve3) -> bool { matches!(curve, Curve3::Line(_)) }
