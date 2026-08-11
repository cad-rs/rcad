// OCCT IntTools_CommonPrt
use glam::DVec3;

/// OCCT TopAbs_ShapeEnum — the type of the common part (myType).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CommonPrtType {
    Vertex,
    Edge,
}

#[derive(Debug, Clone)]
pub struct CommonPrt {
    /// OCCT myType — TopAbs_EDGE or TopAbs_VERTEX.
    pub type_: CommonPrtType,
    /// OCCT myRange1 (edge 1 parameter range).
    pub range1: [f64; 2],
    /// OCCT myRanges2 (edge 2 parameter ranges).
    pub ranges2: Vec<[f64; 2]>,
    /// OCCT myVertPar1 / myVertPar2.
    pub vertex_param1: f64,
    pub vertex_param2: f64,
    /// OCCT myPnt1 / myPnt2 (BoundingPoints).
    pub bounding_point1: DVec3,
    pub bounding_point2: DVec3,
    /// OCCT myAllNullFlag.
    pub all_null_flag: bool,
}

impl CommonPrt {
    pub fn new() -> Self {
        CommonPrt {
            type_: CommonPrtType::Vertex,
            range1: [0.0, 0.0],
            ranges2: Vec::new(),
            vertex_param1: 0.0,
            vertex_param2: 0.0,
            bounding_point1: DVec3::ZERO,
            bounding_point2: DVec3::ZERO,
            all_null_flag: false,
        }
    }

    /// OCCT SetType(TopAbs_ShapeEnum).
    pub fn set_type(&mut self, t: CommonPrtType) {
        self.type_ = t;
    }
    /// OCCT Type().
    pub fn get_type(&self) -> CommonPrtType {
        self.type_
    }
    /// OCCT SetRange1(f, l).
    pub fn set_range1(&mut self, f: f64, l: f64) {
        self.range1 = [f, l];
    }
    /// OCCT SetVertexParameter1(t).
    pub fn set_vertex_parameter1(&mut self, t: f64) {
        self.vertex_param1 = t;
    }
    /// OCCT SetBoundingPoints(p1, p2).
    pub fn set_bounding_points(&mut self, p1: DVec3, p2: DVec3) {
        self.bounding_point1 = p1;
        self.bounding_point2 = p2;
    }
    /// OCCT SetAllNullFlag(f).
    pub fn set_all_null_flag(&mut self, f: bool) {
        self.all_null_flag = f;
    }
    /// OCCT AllNullFlag().
    pub fn all_null_flag(&self) -> bool {
        self.all_null_flag
    }

    /// rcad convenience constructor — a VERTEX common part.
    pub fn new_vertex(t1: f64, t2: f64, p: DVec3, r1: [f64; 2], r2: Vec<[f64; 2]>) -> Self {
        CommonPrt {
            type_: CommonPrtType::Vertex,
            range1: r1,
            ranges2: r2,
            vertex_param1: t1,
            vertex_param2: t2,
            bounding_point1: p,
            bounding_point2: p,
            all_null_flag: false,
        }
    }

    /// rcad convenience constructor — an EDGE common part.
    pub fn new_edge(r1: [f64; 2], r2: Vec<[f64; 2]>, p1: DVec3, p2: DVec3) -> Self {
        let vp2 = r2.first().map_or(0.0, |r| r[0]);
        CommonPrt {
            type_: CommonPrtType::Edge,
            range1: r1,
            ranges2: r2,
            vertex_param1: r1[0],
            vertex_param2: vp2,
            bounding_point1: p1,
            bounding_point2: p2,
            all_null_flag: false,
        }
    }
}
