//! rcad-algo: OCCT TKBO — boolean operation algorithms (new DS-based).
//!
//! Depends only on rcad-kernel (TKMath + TKGeomBase) and rcad-brep (TKBRep).

pub mod bop;
pub mod tolerance;
pub mod bnd_box;
pub mod classify;
pub mod history;
pub mod triangulate;
pub mod pipeline_dump;
pub mod brep_check;
pub mod geom2d_api;
pub mod geom_populate;
pub mod math_utils;
pub use rcad_kernel::math::math_poly::solve_quartic;
pub use rcad_kernel::math::lin::inverse_3x3;
pub use rcad_kernel::math::opt::golden_section_max;

// Re-export from rcad-kernel (TKMath/TKGeomBase equivalents)
pub use rcad_kernel::base::bnd_lib;
pub use rcad_kernel::base::extrema;
