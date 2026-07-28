//! rcad-algo: OCCT TKBO — boolean operation algorithms.
//!
//! Depends only on rcad-kernel (TKMath + TKGeomBase) and rcad-brep (TKBRep).

pub mod bop;

// Re-export boolean operation API at top level
pub use crate::bop::brep_algo_api::{common, cut, fuse};
pub use crate::bop::algo::BooleanOpType;

// Re-export from rcad-kernel
pub use rcad_kernel::base::bnd_lib;
pub use rcad_kernel::base::extrema;
pub use rcad_kernel::math::math_poly::solve_quartic;
pub use rcad_kernel::math::lin::inverse_3x3;
pub use rcad_kernel::math::opt::golden_section_max;
