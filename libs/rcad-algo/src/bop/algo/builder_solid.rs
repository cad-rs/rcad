// OCCT BOPAlgo_BuilderSolid — solid construction from face set.
use crate::bop::ds::DS;

pub struct BuilderSolid {
    ds: DS,
}
impl BuilderSolid {
    pub fn new(ds: DS) -> Self { BuilderSolid { ds } }
    pub fn perform(&mut self) {}
}
