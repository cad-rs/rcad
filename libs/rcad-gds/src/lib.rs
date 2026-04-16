//! GDSII (GDS) file format support for RCAD.
//!
//! Provides import and export capabilities for GDSII layout files,
//! with 2D-to-3D extrusion support.

pub mod error;
pub mod types;
pub mod layer_config;
pub mod reader;
pub mod writer;
pub mod convert;

pub use error::GdsError;
pub use types::{GdsLibrary, GdsStructure, GdsBoundary, GdsPath, GdsText, GdsReference, GdsUnits};
pub use layer_config::{LayerConfig, LayerSettings};
pub use reader::GdsReader;
pub use writer::GdsWriter;
