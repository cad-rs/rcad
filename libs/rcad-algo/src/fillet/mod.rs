//! OCCT TKFillet — ChFi3d / BRepFilletAPI fillet algorithms.
//!
//! `occt` holds the OCCT-aligned 1:1 translation of TKFillet (ChFiDS data
//! structures, ChFi3d_Builder/FilBuilder/ChBuilder, BRepFilletAPI_MakeFillet
//! / MakeChamfer); `fillet` holds the legacy `make_fillet_edge` compat
//! helper pending absorption into the aligned pipeline.

pub mod fillet;
pub mod occt;

pub use fillet::make_fillet_edge;
