// OCCT IntTools_PntOnFace — a point on a face with UV coordinates.
use glam::{DVec2, DVec3};
#[derive(Debug, Clone)]
pub struct PntOnFace {
    pub uv: DVec2,
    pub pnt: DVec3,
    pub index: i32,
}
impl PntOnFace {
    pub fn new() -> Self { PntOnFace { uv: DVec2::ZERO, pnt: DVec3::ZERO, index: -1 } }
    pub fn set_uv(&mut self, u: f64, v: f64) { self.uv = DVec2::new(u, v); }
    pub fn uv(&self) -> DVec2 { self.uv }
    pub fn set_pnt(&mut self, p: DVec3) { self.pnt = p; }
    pub fn pnt(&self) -> DVec3 { self.pnt }
    pub fn set_index(&mut self, i: i32) { self.index = i; }
    pub fn index(&self) -> i32 { self.index }
}
impl Default for PntOnFace { fn default() -> Self { Self::new() } }
