// OCCT IntTools_CommonPrt
use glam::DVec3;
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
impl CommonPrt {
    pub fn new_vertex(t1: f64, t2: f64, p: DVec3, r1: [f64; 2], r2: Vec<[f64; 2]>) -> Self {
        CommonPrt { is_edge: false, range1: r1, ranges2: r2, vertex_param1: t1, vertex_param2: t2, bounding_point1: p, bounding_point2: p }
    }
    pub fn new_edge(r1: [f64; 2], r2: Vec<[f64; 2]>, p1: DVec3, p2: DVec3) -> Self {
        let vp2 = r2.first().map_or(0.0, |r| r[0]);
        CommonPrt { is_edge: true, range1: r1, ranges2: r2, vertex_param1: r1[0], vertex_param2: vp2, bounding_point1: p1, bounding_point2: p2 }
    }
}
