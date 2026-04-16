//! Unified I/O for layout formats (GDS, OASIS).

pub mod error;
pub mod format;
pub mod traits;
pub mod detection;
pub mod layer_config;

pub use format::LayoutFormat;
pub use detection::detect_format;
pub use error::{IoError, Result};

// Re-export from sub-crates (GDS)
pub use rcad_gds::{GdsLibrary, GdsStructure, GdsBoundary, GdsPath, GdsText, GdsReference, GdsUnits, Transform2D};
pub use rcad_gds::{GdsReader, GdsWriter};
pub use rcad_gds::GdsError;

// Re-export from sub-crates (OASIS)
pub use rcad_oas::{OasLibrary, OasCell, OasPolygon, OasPath, OasText, OasPlacement};
pub use rcad_oas::{OasReader, OasWriter};
pub use rcad_oas::OasError;

// Re-export unified layer config
pub use layer_config::UnifiedLayerConfig;
