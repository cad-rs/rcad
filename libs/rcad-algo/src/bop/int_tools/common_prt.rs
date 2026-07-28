// OCCT IntTools_CommonPrt — common part between two edges.
use glam::DVec3;

/// IntTools_CommonPrt: describes a common part (vertex or edge) between two edges.
#[derive(Debug, Clone)]
pub struct CommonPrt {
    /// true = EDGE type (overlap), false = VERTEX type (point touch)
    pub is_edge: bool,
    /// Parameter range on the first edge
    pub range1: [f64; 2],
    /// Parameter ranges on the second edge (1-to-N mapping)
    pub ranges2: Vec<[f64; 2]>,
    /// Parameter on edge1 for vertex-type (OCCT: VertexParameter1)
    pub vertex_param1: f64,
    /// Parameter on edge2 for vertex-type
    pub vertex_param2: f64,
    /// 3D bounding point (OCCT: BoundingPoint1/2)
    pub bounding_point1: DVec3,
    pub bounding_point2: DVec3,
    /// All-null flag — true when both ranges are degenerate
    pub all_null: bool,
}

impl CommonPrt {
    /// OCCT: default constructor
    pub fn new() -> Self {
        CommonPrt {
            is_edge: false,
            range1: [0.0, 0.0],
            ranges2: Vec::new(),
            vertex_param1: 0.0,
            vertex_param2: 0.0,
            bounding_point1: DVec3::ZERO,
            bounding_point2: DVec3::ZERO,
            all_null: false,
        }
    }

    pub fn set_type(&mut self, is_edge: bool) { self.is_edge = is_edge; }
    pub fn is_edge_type(&self) -> bool { self.is_edge }
    pub fn set_ranges1(&mut self, r: [f64; 2]) { self.range1 = r; }
    pub fn set_ranges2(&mut self, r: Vec<[f64; 2]>) { self.ranges2 = r; }
    pub fn set_vertex_parameter1(&mut self, p: f64) { self.vertex_param1 = p; }
    pub fn vertex_parameter1(&self) -> f64 { self.vertex_param1 }
    pub fn set_vertex_parameter2(&mut self, p: f64) { self.vertex_param2 = p; }
    pub fn vertex_parameter2(&self) -> f64 { self.vertex_param2 }
    pub fn set_bounding_points(&mut self, p1: DVec3, p2: DVec3) { self.bounding_point1 = p1; self.bounding_point2 = p2; }
    pub fn bounding_point1(&self) -> DVec3 { self.bounding_point1 }
    pub fn bounding_point2(&self) -> DVec3 { self.bounding_point2 }
    pub fn set_all_null(&mut self, v: bool) { self.all_null = v; }
    pub fn is_all_null(&self) -> bool { self.all_null }

    pub fn new_vertex(t1: f64, t2: f64, p: DVec3, range1: [f64; 2], ranges2: Vec<[f64; 2]>) -> Self {
        CommonPrt {
            is_edge: false,
            range1,
            ranges2,
            vertex_param1: t1,
            vertex_param2: t2,
            bounding_point1: p,
            bounding_point2: p,
            all_null: false,
        }
    }

    pub fn new_edge(range1: [f64; 2], ranges2: Vec<[f64; 2]>, p1: DVec3, p2: DVec3) -> Self {
        CommonPrt {
            is_edge: true,
            range1,
            ranges2,
            vertex_param1: range1[0],
            vertex_param2: ranges2.first().map_or(0.0, |r| r[0]),
            bounding_point1: p1,
            bounding_point2: p2,
            all_null: false,
        }
    }
}

impl Default for CommonPrt { fn default() -> Self { Self::new() } }
