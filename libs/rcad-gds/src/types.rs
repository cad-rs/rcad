use glam::DVec2;
use serde::{Deserialize, Serialize};

/// GDS units information.
///
/// In GDSII, units are specified as:
/// - user_unit: size of one user unit in meters (e.g., 1e-6 for microns)
/// - meter_unit: size of one database unit in meters (e.g., 1e-9 for nanometers)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct GdsUnits {
    /// Size of one user unit in meters (e.g., 1e-6 for microns)
    pub user_unit: f64,
    /// Size of one database unit in meters (e.g., 1e-9 for nanometers)
    pub meter_unit: f64,
}

impl Default for GdsUnits {
    fn default() -> Self {
        Self {
            user_unit: 1e-6,  // 1 micron
            meter_unit: 1e-9, // 1 nanometer
        }
    }
}

/// GDS library (top-level container).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum EndCapType {
    #[default]
    Flush,
    Round,
    Square,
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
    pub rotation: f64, // radians
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
            result.y = -result.y; // Reflect across X-axis (flip Y)
        }
        let (sin, cos) = self.rotation.sin_cos();
        let rotated = DVec2::new(
            result.x * cos - result.y * sin,
            result.x * sin + result.y * cos,
        );
        rotated + self.translation
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: DVec2, b: DVec2) {
        let tolerance = 1e-10;
        assert!((a.x - b.x).abs() < tolerance, "x: {} != {}", a.x, b.x);
        assert!((a.y - b.y).abs() < tolerance, "y: {} != {}", a.y, b.y);
    }

    #[test]
    fn test_identity_transform() {
        let transform = Transform2D::identity();
        let point = DVec2::new(5.0, 10.0);
        approx_eq(transform.transform_point(point), point);
    }

    #[test]
    fn test_translation_only() {
        let transform = Transform2D::from_translation(100.0, 200.0);
        let point = DVec2::new(5.0, 10.0);
        approx_eq(transform.transform_point(point), DVec2::new(105.0, 210.0));
    }

    #[test]
    fn test_rotation_90_degrees() {
        let transform = Transform2D {
            translation: DVec2::ZERO,
            rotation: std::f64::consts::FRAC_PI_2, // 90 degrees
            reflection: false,
            magnification: 1.0,
        };
        let point = DVec2::new(1.0, 0.0);
        let result = transform.transform_point(point);
        // 90 deg rotation: (1, 0) -> (0, 1)
        approx_eq(result, DVec2::new(0.0, 1.0));
    }

    #[test]
    fn test_reflection_across_x_axis() {
        let transform = Transform2D {
            translation: DVec2::ZERO,
            rotation: 0.0,
            reflection: true,
            magnification: 1.0,
        };
        // Point (2, 3) reflected across X-axis should be (2, -3)
        let point = DVec2::new(2.0, 3.0);
        let result = transform.transform_point(point);
        approx_eq(result, DVec2::new(2.0, -3.0));

        // Point with negative Y
        let point2 = DVec2::new(5.0, -4.0);
        let result2 = transform.transform_point(point2);
        approx_eq(result2, DVec2::new(5.0, 4.0));
    }

    #[test]
    fn test_combined_reflection_and_translation() {
        let transform = Transform2D {
            translation: DVec2::new(10.0, 20.0),
            rotation: 0.0,
            reflection: true,
            magnification: 1.0,
        };
        // Point (2, 3) reflected -> (2, -3), then translated -> (12, 17)
        let point = DVec2::new(2.0, 3.0);
        let result = transform.transform_point(point);
        approx_eq(result, DVec2::new(12.0, 17.0));
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
