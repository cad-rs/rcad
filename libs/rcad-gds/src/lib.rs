//! GDSII (GDS) file format support for RCAD.
//!
//! Provides import and export capabilities for GDSII layout files,
//! with 2D-to-3D extrusion support.

pub mod convert;
pub mod error;
pub mod layer_config;
pub mod reader;
pub mod types;
pub mod writer;

pub use error::GdsError;
pub use layer_config::{LayerConfig, LayerSettings};
pub use reader::GdsReader;
pub use types::{
    GdsBoundary, GdsLibrary, GdsPath, GdsReference, GdsStructure, GdsText, GdsUnits, Transform2D,
};
pub use writer::GdsWriter;
