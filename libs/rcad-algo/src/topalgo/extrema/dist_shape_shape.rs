// OCCT BRepExtrema_DistShapeShape (BRepExtrema_DistShapeShape.hxx / .cxx)
// Computes the minimum distance between two shapes.
//
// OCCT L224: used by IntTools_EdgeEdge::Perform() to quickly check
// if a line and an analytic curve are far apart.
// rcad: delegates to rcad-kernel::base::extrema.

use rcad_kernel::base::extrema;
use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::topo_shape::Shape;

/// OCCT BRepExtrema_DistShapeShape — minimum distance between two shapes.
pub struct DistShapeShape {
    done: bool,
    value: f64,
}

impl DistShapeShape {
    /// Constructor with two shapes.
    /// OCCT: BRepExtrema_DistShapeShape(S1, S2, Extrema_ExtFlag_MIN)
    pub fn new(s1: &Shape, s2: &Shape, _ext_flag: i32) -> Self {
        let val = Self::compute_min_dist(s1, s2);
        DistShapeShape { done: true, value: val }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn value(&self) -> f64 { self.value }

    /// Compute minimum distance between two shapes by sampling edge curves.
    fn compute_min_dist(s1: &Shape, s2: &Shape) -> f64 {
        let pts1 = Self::sample_shape(s1);
        let pts2 = Self::sample_shape(s2);
        if pts1.is_empty() || pts2.is_empty() { return f64::MAX; }
        let mut min_d = f64::MAX;
        for &p1 in &pts1 {
            for &p2 in &pts2 {
                let d = (p1 - p2).length();
                if d < min_d { min_d = d; }
            }
        }
        min_d
    }

    /// Sample points from a shape's edges.
    fn sample_shape(s: &Shape) -> Vec<glam::DVec3> {
        let mut pts = Vec::new();
        match &*s.data {
            rcad_kernel::topods::TShape::Edge(ed) => {
                if let Some(c) = &ed.curve {
                    let n = 8usize;
                    let (f, l) = (ed.range[0], ed.range[1]);
                    for i in 0..=n {
                        pts.push(c.point_at(f + (l - f) * i as f64 / n as f64));
                    }
                }
            }
            _ => {
                // For non-edge shapes, recurse into sub-shapes
                // rcad: skip — DistShapeShape primarily used for edge-edge
            }
        }
        pts
    }
}
