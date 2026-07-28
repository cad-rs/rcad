// OCCT IntTools_Range — parameter range.
#[derive(Debug, Clone, Copy)]
pub struct Range {
    pub first: f64,
    pub last: f64,
}
impl Range {
    pub fn new() -> Self { Range { first: 0.0, last: 0.0 } }
    pub fn new_range(f: f64, l: f64) -> Self { Range { first: f, last: l } }
    pub fn set_first(&mut self, f: f64) { self.first = f; }
    pub fn first(&self) -> f64 { self.first }
    pub fn set_last(&mut self, l: f64) { self.last = l; }
    pub fn last(&self) -> f64 { self.last }
    pub fn range(&self) -> (f64, f64) { (self.first, self.last) }
}
impl Default for Range { fn default() -> Self { Self::new() } }
