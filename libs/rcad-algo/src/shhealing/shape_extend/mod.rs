//! OCCT ShapeExtend package (TKShHealing) — extended data structures for
//! shape analysis and fixing.

pub mod basic_msg_registrator;
pub mod complex_curve;
pub mod composite_surface;
pub mod explorer;
pub mod msg;
pub mod msg_registrator;
pub mod status;
pub mod wire_data;

#[cfg(test)]
mod wire_data_tests;

pub use basic_msg_registrator::BasicMsgRegistrator;
pub use complex_curve::ShapeExtendComplexCurve;
pub use composite_surface::{ShapeExtendCompositeSurface, SurfaceResD1, SurfaceResD2, SurfaceResD3};
pub use explorer::ShapeExtendExplorer;
pub use msg::MessageMsg;
pub use msg_registrator::MsgRegistrator;
pub use status::{
    ShapeExtendParametrisation, ShapeExtendStatus, decode_status, encode_status, init,
};
pub use wire_data::WireData;
