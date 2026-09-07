//! OCCT TKFillet — modules named 1:1 after the OCCT TKFillet packages:
//! `chfi2d*` (ChFi2d), `chfi_ds` (ChFiDS), `chfi3d*` (ChFi3d),
//! `brep_fillet_api` (BRepFilletAPI).
//!
//! `fillet` holds the legacy `make_fillet_edge` compatibility helper
//! (blend on a single edge) pending absorption into the aligned pipeline.

pub mod brep_blend;
pub mod brep_blend_curv_point_rad_inv;
pub mod brep_blend_func;
pub mod brep_blend_func_inv;
pub mod brep_blend_function;
pub mod brep_blend_point;
pub mod brep_blend_ruled;
pub mod brep_fillet_api;
pub mod chfi2d;
pub mod chfi2d_ana_fillet_algo;
pub mod chfi2d_builder;
pub mod chfi2d_builder_0;
pub mod chfi2d_chamfer_api;
pub mod chfi2d_fillet_algo;
pub mod chfi2d_fillet_api;
pub mod chfi3d;
pub mod chfi3d_builder_0;
pub mod chfi3d_builder_0_filds;
pub mod chfi3d_builder_c1;
pub mod chfi3d_builder_spkp;
pub mod chfi3d_ds;
pub mod chfi_ds;
pub mod chfi_ds_spine;
pub mod chfi_ds_stripe;
pub mod chfi_ds_surfdata;
pub mod chfi_kpart;
pub mod fillet;
pub mod hbuilder;

pub use fillet::make_fillet_edge;
