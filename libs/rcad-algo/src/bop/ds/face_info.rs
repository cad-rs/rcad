// OCCT BOPDS_FaceInfo 1:1 translation.
use indexmap::IndexSet;

/// BOPDS_FaceInfo — per-face intersection state.
#[derive(Debug, Clone, Default)]
pub struct FaceInfo {
    pub face_index: usize,
    pub pave_blocks_in: IndexSet<usize>,
    pub pave_blocks_on: IndexSet<usize>,
    pub pave_blocks_sc: IndexSet<usize>,
    pub curves_sc: IndexSet<usize>,
    pub vertices_in: IndexSet<usize>,
    pub vertices_on: IndexSet<usize>,
    pub vertices_sc: IndexSet<usize>,
}

impl FaceInfo {
    pub fn clear(&mut self) {
        self.face_index = 0;
        self.pave_blocks_in.clear();
        self.pave_blocks_on.clear();
        self.pave_blocks_sc.clear();
        self.curves_sc.clear();
        self.vertices_in.clear();
        self.vertices_on.clear();
        self.vertices_sc.clear();
    }
    pub fn set_index(&mut self, i: usize) { self.face_index = i; }
    pub fn index(&self) -> usize { self.face_index }
    pub fn pave_blocks_in(&self) -> &IndexSet<usize> { &self.pave_blocks_in }
    pub fn change_pave_blocks_in(&mut self) -> &mut IndexSet<usize> { &mut self.pave_blocks_in }
    pub fn vertices_in(&self) -> &IndexSet<usize> { &self.vertices_in }
    pub fn change_vertices_in(&mut self) -> &mut IndexSet<usize> { &mut self.vertices_in }
    pub fn pave_blocks_on(&self) -> &IndexSet<usize> { &self.pave_blocks_on }
    pub fn change_pave_blocks_on(&mut self) -> &mut IndexSet<usize> { &mut self.pave_blocks_on }
    pub fn vertices_on(&self) -> &IndexSet<usize> { &self.vertices_on }
    pub fn change_vertices_on(&mut self) -> &mut IndexSet<usize> { &mut self.vertices_on }
    pub fn pave_blocks_sc(&self) -> &IndexSet<usize> { &self.pave_blocks_sc }
    pub fn change_pave_blocks_sc(&mut self) -> &mut IndexSet<usize> { &mut self.pave_blocks_sc }
    pub fn vertices_sc(&self) -> &IndexSet<usize> { &self.vertices_sc }
    pub fn change_vertices_sc(&mut self) -> &mut IndexSet<usize> { &mut self.vertices_sc }
    pub fn curves_sc_only(&self) -> Vec<usize> { self.curves_sc.iter().copied().collect() }
    pub fn has_any_interference(&self) -> bool {
        !self.pave_blocks_in.is_empty() || !self.pave_blocks_on.is_empty()
            || !self.pave_blocks_sc.is_empty() || !self.curves_sc.is_empty()
            || !self.vertices_in.is_empty() || !self.vertices_sc.is_empty()
            || !self.vertices_on.is_empty()
    }
}
