//! OCCT ShapeExtend package (TKShHealing) — extended data structures for
//! shape analysis and fixing.

pub mod status;
pub mod wire_data;

#[cfg(test)]
mod wire_data_tests;

pub use status::{ShapeExtendStatus, decode_status, encode_status};
pub use wire_data::WireData;
