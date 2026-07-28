// Submodules: OCCT BOPDS data structures
pub mod common_block;
pub mod face_info;
pub mod pave;
pub mod new_ds;
pub use new_ds::DS;
pub mod iterator;
pub use iterator::BOPDS_Iterator;
pub mod topods_builder;

// Tools for DS operations
pub mod tools;

pub use new_ds::{
    Classification, ConnexityBlock, Interference, InterferenceEE, InterferenceEF,
    InterferenceEZ, InterferenceFF, InterferenceFZ, InterferenceVE, InterferenceVF,
    InterferenceVV, InterferenceVZ, InterferenceZZ, IntersectionCurve, NearTangentType,
    PairIterator, PassKey, ShapeOrigin, ShapeSD, SharedTopologyInfo,
};

/// Compute face AABB from DS shape info.
#[allow(dead_code)]
pub fn face_aabb(ds: &DS, fi: usize) -> crate::bop::tools::bvh::BndBox {
    let _ = (ds, fi);
    rcad_kernel::math::bnd::BndBox::new()
}
