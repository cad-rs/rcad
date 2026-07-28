// OCCT BOPDS_Iterator 鈥?BVH-based pair enumeration.
use crate::bop::ds::DS;
use crate::bop::tools::bvh::Aabb;
/// BOPDS_Iterator 鈥?iterates pairs of interfering shapes.
pub struct BOPDS_Iterator<'a> {
    ds: &'a DS,
    fuzzy_tol: f64,
    my_lists: Vec<Vec<(usize, usize)>>,
    current_list: Vec<(usize, usize)>,
    current_pos: usize,
}
impl<'a> BOPDS_Iterator<'a> {
    pub fn new(ds: &'a DS, fuzzy_tol: f64) -> Self {
        let n = 10; // NbInterfTypes
        let mut my_lists = Vec::with_capacity(n);
        for _ in 0..n { my_lists.push(Vec::new()); }
        BOPDS_Iterator {
            ds, fuzzy_tol,
            my_lists,
            current_list: Vec::new(),
            current_pos: 0,
        }
    }
    pub fn prepare(&mut self) {
        // TODO: populate pair lists based on DS shape AABBs
    }
    pub fn initialize(&mut self, t1: u32, t2: u32) {
        let idx = Self::type_to_integer(t1, t2);
        self.current_list = if idx < self.my_lists.len() {
            self.my_lists[idx].clone()
        } else {
            Vec::new()
        };
        self.current_pos = 0;
    }
    pub fn more(&self) -> bool { self.current_pos < self.current_list.len() }
    pub fn next(&mut self) { self.current_pos += 1; }
    pub fn value(&self) -> (usize, usize) { self.current_list[self.current_pos] }
    pub fn pairs(&self, t1: u32, t2: u32) -> &[(usize, usize)] {
        let idx = Self::type_to_integer(t1, t2);
        if idx < self.my_lists.len() { &self.my_lists[idx] } else { &[] }
    }
    fn type_to_integer(t1: u32, t2: u32) -> usize {
        match (t1, t2) {
            (0, 0) => 0, // VV
            (0, 1) => 1, // VE
            (1, 1) => 2, // EE
            (0, 2) => 3, // VF
            (1, 2) => 4, // EF
            (2, 2) => 5, // FF
            (0, 3) => 6, // VZ
            (1, 3) => 7, // EZ
            (2, 3) => 8, // FZ
            (3, 3) => 9, // ZZ
            _ => 0,
        }
    }
}