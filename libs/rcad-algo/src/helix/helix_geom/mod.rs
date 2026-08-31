//! OCCT HelixGeom package (TKHelix) — helix curve construction.
//!
//! Ported classes:
//! - [`helix_curve::HelixCurve`]      (HelixGeom_HelixCurve)
//! - [`tools`]                        (HelixGeom_Tools)
//! - [`builder_approx_curve`]         (HelixGeom_BuilderApproxCurve)
//! - [`builder_helix_gen::BuilderHelixGen`]     (HelixGeom_BuilderHelixGen)
//! - [`builder_helix_coil::BuilderHelixCoil`]   (HelixGeom_BuilderHelixCoil)
//! - [`builder_helix::BuilderHelix`]  (HelixGeom_BuilderHelix)

pub mod builder_approx_curve;
pub mod builder_helix;
pub mod builder_helix_coil;
pub mod builder_helix_gen;
pub mod helix_curve;
pub mod tools;

pub use builder_helix::BuilderHelix;
pub use builder_helix_coil::BuilderHelixCoil;
pub use builder_helix_gen::BuilderHelixGen;
pub use helix_curve::HelixCurve;
