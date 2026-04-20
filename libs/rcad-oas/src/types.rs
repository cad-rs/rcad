use glam::DVec2;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// OASIS library (top-level container).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct OasLibrary {
    pub name: Option<String>,
    pub cells: HashMap<String, OasCell>,
}

/// OASIS cell.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct OasCell {
    pub name: String,
    pub polygons: Vec<OasPolygon>,
    pub paths: Vec<OasPath>,
    pub texts: Vec<OasText>,
    pub placements: Vec<OasPlacement>,
}

/// OASIS polygon.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OasPolygon {
    pub layer: i32,
    pub datatype: i32,
    pub points: Vec<DVec2>,
}

/// OASIS path.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OasPath {
    pub layer: i32,
    pub datatype: i32,
    pub width: f64,
    pub points: Vec<DVec2>,
    #[serde(default)]
    pub start_extension: f64,
    #[serde(default)]
    pub end_extension: f64,
}

/// OASIS text.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OasText {
    pub layer: i32,
    pub text_type: i32,
    pub position: DVec2,
    pub content: String,
}

/// OASIS cell placement.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OasPlacement {
    pub cell_name: String,
    pub x: f64,
    pub y: f64,
    pub rotation: f64,
    pub reflection: bool,
    pub magnification: f64,
    pub array: Option<(u32, u32, f64, f64, f64, f64)>,
}
