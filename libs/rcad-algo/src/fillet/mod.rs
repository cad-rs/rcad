//! OCCT TKFillet — ChFi3d / BRepFilletAPI fillet algorithms.
//!
//! Current content migrated from `algo_ext`: the legacy
//! `make_fillet_edge` compatibility helper (blend on a single edge).

pub mod fillet;

pub use fillet::make_fillet_edge;
