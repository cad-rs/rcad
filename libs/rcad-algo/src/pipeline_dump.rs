//! Pipeline dump — minimal stub for boolean pipeline debugging.
pub struct DumpCtx;
impl DumpCtx {
    pub fn new(_grid: &str, _case: &str) -> Self { Self }
    pub fn new_with_module(_grid: &str, _case: &str, _module: &str) -> Self { Self }
    pub fn snapshot(&self, _stage: &str, _ds: &crate::bop::ds::DS, _brep: Option<&rcad_kernel::topods::BRep>) {}
}
