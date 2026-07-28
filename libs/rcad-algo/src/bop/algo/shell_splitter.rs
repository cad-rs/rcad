// OCCT BOPAlgo_ShellSplitter ?shell partitioning.
use crate::bop::ds::DS;

pub struct ShellSplitter {
    ds: DS,
}
impl ShellSplitter {
    pub fn new(ds: DS) -> Self { ShellSplitter { ds } }
    pub fn perform(&mut self) {}
}

/// Build connexity blocks from connected faces.
pub fn make_connexity_blocks(_faces: &[usize], _ds: &DS, _out: &mut Vec<Vec<usize>>) {}

