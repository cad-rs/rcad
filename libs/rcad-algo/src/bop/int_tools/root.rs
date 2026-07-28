// OCCT IntTools_Root
#[derive(Debug, Clone)]
pub struct Root {
    pub root: f64, pub first: f64, pub last: f64,
    pub is_infinite: bool, pub is_valid: bool,
}
impl Root {
    pub fn new() -> Self { Root { root: 0.0, first: 0.0, last: 0.0, is_infinite: false, is_valid: false } }
    pub fn root(&self) -> f64 { self.root }
    pub fn set_root(&mut self, r: f64) { self.root = r; }
    pub fn interval(&self) -> (f64, f64) { (self.first, self.last) }
    pub fn set_interval(&mut self, f: f64, l: f64) { self.first = f; self.last = l; }
    pub fn is_infinite(&self) -> bool { self.is_infinite }
    pub fn set_infinite(&mut self, v: bool) { self.is_infinite = v; }
    pub fn is_valid(&self) -> bool { self.is_valid }
    pub fn set_valid(&mut self, v: bool) { self.is_valid = v; }
}
impl Default for Root { fn default() -> Self { Self::new() } }
