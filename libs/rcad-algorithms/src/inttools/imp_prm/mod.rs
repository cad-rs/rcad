//! ✅ Placeholder: IntPatch_ImpPrmIntersection — analytic-parametric surface intersection.
//!
//! OCCT IntPatch_ImpPrmIntersection.hxx / .cxx (3892 lines).
//! See intersection.rs for the main Perform method.
//!
//! WIP: supporting types SOnBounds, SearchInside, IWalking are placeholders.
//! Full 1:1 translation is in progress.

pub mod arc_function;
pub mod surf_function;
pub mod s_on_bounds;
pub mod search_inside;
pub mod i_walking;
pub mod decompose;
pub mod intersection;

pub use intersection::ImpPrmIntersection;
