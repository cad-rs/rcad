// OCCT BOPAlgo_BuilderFace — face splitting.
use crate::bop::ds::DS;

pub struct BuilderFace {
    ds: DS,
}
impl BuilderFace {
    pub fn new(ds: DS) -> Self { BuilderFace { ds } }
    pub fn perform(&mut self) {}
}
