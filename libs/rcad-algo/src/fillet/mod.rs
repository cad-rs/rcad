//! OCCT TKFillet — modules named 1:1 after the OCCT TKFillet packages:
//! `chfi_ds` (ChFiDS), `chfi3d` (ChFi3d), `brep_fillet_api` (BRepFilletAPI).
//!
//! `fillet` holds the legacy `make_fillet_edge` compatibility helper
//! (blend on a single edge) pending absorption into the aligned pipeline.

pub mod brep_fillet_api;
pub mod chfi3d;
pub mod chfi3d_builder_0;
pub mod chfi_ds;
pub mod chfi_kpart;
pub mod fillet;
pub mod topopebrepds;

pub use fillet::make_fillet_edge;
