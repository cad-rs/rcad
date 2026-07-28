// OCCT BOPAlgo_Builder — shape construction from DS.

pub use crate::bop::algo::BooleanOpType;
use crate::bop::algo::{GlueEnum, Report};
use crate::bop::ds::DS;

/// BOPAlgo_Builder — result builder for boolean operations.
pub struct BooleanBuilder<'a> {
    ds: &'a DS,
    my_report: Report,
    my_operation: BooleanOpType,
    my_glue: GlueEnum,
}

impl<'a> BooleanBuilder<'a> {
    pub fn new(ds: &'a DS, op: BooleanOpType) -> Self {
        BooleanBuilder {
            ds, my_report: Report::new(),
            my_operation: op, my_glue: GlueEnum::GlueOff,
        }
    }
    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }
    pub fn report(&self) -> &Report { &self.my_report }
    pub fn build(&mut self) -> Result<rcad_kernel::BRep, ()> {
        Ok(rcad_kernel::BRep::new())
    }
}

/// Project a 3D point to UV on a surface (placeholder).
pub fn world_to_uv(_p: glam::DVec3, _surf: &rcad_kernel::geom::Surface3) -> glam::DVec2 { glam::DVec2::ZERO }

/// Boolean operation error type.
#[derive(Debug, Clone)]
pub enum BooleanError {
    InvalidOperation,
    TooFewArguments,
    NoFiller,
    BOPNotAllowed,
    BOPNotSet,
    EmptyShape,
    EmptyInput,
    DegenerateResult,
    NumericalFailure(&'static str),
    InvalidResult(&'static str),
}
