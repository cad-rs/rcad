//! OCCT BRepAlgoAPI_BuilderOperation 閳?boolean operation wrapper.
use crate::bop::algo::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::bop::ds::DS;
pub struct BuilderOperation {
    pub shape_a: rcad_kernel::BRep,
    pub shape_b: rcad_kernel::BRep,
    pub op_type: BooleanOpType,
}
impl BuilderOperation {
    pub fn new(a: rcad_kernel::BRep, b: rcad_kernel::BRep, op: BooleanOpType) -> Self {
        BuilderOperation { shape_a: a, shape_b: b, op_type: op }
    }
    pub fn perform(&mut self) -> Result<rcad_kernel::BRep, BooleanError> {
        let mut ds = DS::new();
        // TODO: populate DS from BReps
        let mut builder = BooleanBuilder::new(&ds, self.op_type);
        builder.build().map_err(|_| BooleanError::InvalidOperation)
    }
}