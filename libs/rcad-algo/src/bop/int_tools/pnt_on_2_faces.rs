// OCCT IntTools_PntOn2Faces
use glam::DVec3;
#[derive(Debug, Clone)]
pub struct PntOn2Faces { pub pnt: DVec3, pub uv1: (f64, f64), pub uv2: (f64, f64), pub index: i32 }
impl PntOn2Faces {
    pub fn new() -> Self { PntOn2Faces { pnt: DVec3::ZERO, uv1: (0.0, 0.0), uv2: (0.0, 0.0), index: -1 } }
    pub fn set_pnt(&mut self, p: DVec3) { self.pnt = p; }
    pub fn pnt(&self) -> DVec3 { self.pnt }
    pub fn set_uv1(&mut self, u: f64, v: f64) { self.uv1 = (u, v); }
    pub fn uv1(&self) -> (f64, f64) { self.uv1 }
    pub fn set_uv2(&mut self, u: f64, v: f64) { self.uv2 = (u, v); }
    pub fn uv2(&self) -> (f64, f64) { self.uv2 }
    pub fn set_index(&mut self, i: i32) { self.index = i; }
    pub fn index(&self) -> i32 { self.index }
}
impl Default for PntOn2Faces { fn default() -> Self { Self::new() } }

/// OCCT IntSurf_PntOn2S — a point given by its UV parameters on two surfaces.
/// Used by BOPAlgo_PaveFiller::GetEFPnts to pass Edge-Face intersection points
/// into IntTools_FaceFace (SetList) as starting points for intersection curves.
#[derive(Debug, Clone)]
pub struct PntOn2S {
    pub uv1: (f64, f64),
    pub uv2: (f64, f64),
}

impl PntOn2S {
    pub fn new() -> Self {
        PntOn2S { uv1: (0.0, 0.0), uv2: (0.0, 0.0) }
    }
    pub fn set_value(&mut self, u1: f64, v1: f64, u2: f64, v2: f64) {
        self.uv1 = (u1, v1);
        self.uv2 = (u2, v2);
    }
    pub fn parameters(&self) -> (f64, f64, f64, f64) {
        (self.uv1.0, self.uv1.1, self.uv2.0, self.uv2.1)
    }
}
impl Default for PntOn2S { fn default() -> Self { Self::new() } }
