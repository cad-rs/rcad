// OCCT IntTools_Curve
use rcad_kernel::geom::Curve3;
#[derive(Debug, Clone)]
pub struct Curve {
    pub curve: Curve3,
    pub tolerance: f64,
    pub has_forward: bool,
    pub has_reversed: bool,
}
impl Curve {
    pub fn new() -> Self { Curve { curve: Curve3::Line(rcad_kernel::geom::Line3 { origin: glam::DVec3::ZERO, direction: glam::DVec3::X }), tolerance: 0.0, has_forward: false, has_reversed: false } }
    pub fn set_curve(&mut self, c: Curve3) { self.curve = c; }
    pub fn curve(&self) -> &Curve3 { &self.curve }
    pub fn set_tolerance(&mut self, t: f64) { self.tolerance = t; }
    pub fn tolerance(&self) -> f64 { self.tolerance }
    pub fn set_has_forward(&mut self, v: bool) { self.has_forward = v; }
    pub fn has_forward(&self) -> bool { self.has_forward }
    pub fn set_has_reversed(&mut self, v: bool) { self.has_reversed = v; }
    pub fn has_reversed(&self) -> bool { self.has_reversed }
}
impl Default for Curve { fn default() -> Self { Self::new() } }
