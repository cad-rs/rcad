use glam::DVec2;
use serde::{Deserialize, Serialize};

/// GDS units information.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct GdsUnits {
    /// User units (e.g., 1e-6 for microns)
    pub user_unit: f64,
    /// Meters per database unit
    pub meter_unit: f64,
}

impl Default for GdsUnits {
    fn default() -> Self {
        Self {
            user_unit: 1e-6,  // microns
            meter_unit: 1e-9, // nanometers
        }
    }
}

/// GDS library (top-level container).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GdsLibrary {
    pub name: String,
    pub units: GdsUnits,
    pub structures: std::collections::HashMap<String, GdsStructure>,
}

/// GDS structure (cell).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GdsStructure {
    pub name: String,
    pub boundaries: Vec<GdsBoundary>,
    pub paths: Vec<GdsPath>,
    pub texts: Vec<GdsText>,
    pub references: Vec<GdsReference>,
}

/// GDS boundary (closed polygon).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GdsBoundary {
    pub layer: i16,
    pub datatype: i16,
    pub points: Vec<DVec2>,
}

/// GDS path (wire with width).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EndCapType {
    Flush,
    Round,
    Square,
}

impl Default for EndCapType {
    fn default() -> Self {
        Self::Flush
    }
}

/// GDS path (wire).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GdsPath {
    pub layer: i16,
    pub datatype: i16,
    pub width: f64,
    pub points: Vec<DVec2>,
    #[serde(default)]
    pub end_cap: EndCapType,
}

/// GDS text annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GdsText {
    pub layer: i16,
    pub text_type: i16,
    pub position: DVec2,
    pub content: String,
}

/// 2D transform for cell references.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Transform2D {
    pub translation: DVec2,
    pub rotation: f64,  // radians
    pub reflection: bool,
    pub magnification: f64,
}

impl Default for Transform2D {
    fn default() -> Self {
        Self {
            translation: DVec2::ZERO,
            rotation: 0.0,
            reflection: false,
            magnification: 1.0,
        }
    }
}

impl Transform2D {
    pub fn identity() -> Self {
        Self::default()
    }

    pub fn from_translation(x: f64, y: f64) -> Self {
        Self {
            translation: DVec2::new(x, y),
            ..Default::default()
        }
    }

    pub fn transform_point(&self, p: DVec2) -> DVec2 {
        let mut result = p * self.magnification;
        if self.reflection {
            result.x = -result.x;
        }
        let (sin, cos) = self.rotation.sin_cos();
        let rotated = DVec2::new(
            result.x * cos - result.y * sin,
            result.x * sin + result.y * cos,
        );
        rotated + self.translation
    }
}

/// Array parameters for AREF.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ArrayParams {
    pub columns: u16,
    pub rows: u16,
    pub column_offset: DVec2,
    pub row_offset: DVec2,
}

/// GDS cell reference (AREF or SREF).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GdsReference {
    pub cell_name: String,
    pub transform: Transform2D,
    pub array: Option<ArrayParams>,
}
