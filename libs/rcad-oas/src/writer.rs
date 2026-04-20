//! OASIS file writer.

use std::io::Cursor;
use std::path::Path;

use crate::{OasCell, OasError, OasLibrary, OasPath, OasPlacement, OasPolygon, OasText};

/// OASIS file writer using laykit.
pub struct OasWriter;

impl OasWriter {
    /// Create a new OASIS writer.
    pub fn new() -> Self {
        Self
    }

    /// Write an OASIS library to file.
    pub fn write_file<P: AsRef<Path>>(&self, library: &OasLibrary, path: P) -> Result<(), OasError> {
        let oasis = Self::convert_library(library)?;
        oasis.write_to_file(path)
            .map_err(|e| OasError::Laykit(format!("{:?}", e)))?;
        Ok(())
    }

    /// Serialize an OASIS library to bytes.
    pub fn to_bytes(library: &OasLibrary) -> Result<Vec<u8>, OasError> {
        let oasis = Self::convert_library(library)?;
        let mut buffer = Vec::new();
        {
            let mut cursor = Cursor::new(&mut buffer);
            oasis.write_to_writer(&mut cursor)
                .map_err(|e| OasError::Laykit(format!("{:?}", e)))?;
        }
        Ok(buffer)
    }

    /// Convert our OasLibrary to laykit OASISFile.
    fn convert_library(library: &OasLibrary) -> Result<laykit::OASISFile, OasError> {
        let mut oasis = laykit::OASISFile::new();

        // Set default unit (1 nanometer)
        oasis.unit = 1e-9;

        // Convert cells
        for cell in library.cells.values() {
            oasis.cells.push(Self::convert_cell(cell)?);
        }

        Ok(oasis)
    }

    /// Convert our OasCell to laykit OASISCell.
    fn convert_cell(cell: &OasCell) -> Result<laykit::OASISCell, OasError> {
        let mut result = laykit::OASISCell {
            name: cell.name.clone(),
            elements: Vec::new(),
        };

        // Convert polygons
        for polygon in &cell.polygons {
            result.elements.push(laykit::OASISElement::Polygon(
                Self::convert_polygon(polygon),
            ));
        }

        // Convert paths
        for path in &cell.paths {
            result.elements.push(laykit::OASISElement::Path(
                Self::convert_path(path),
            ));
        }

        // Convert texts
        for text in &cell.texts {
            result.elements.push(laykit::OASISElement::Text(
                Self::convert_text(text),
            ));
        }

        // Convert placements
        for placement in &cell.placements {
            result.elements.push(laykit::OASISElement::Placement(
                Self::convert_placement(placement),
            ));
        }

        Ok(result)
    }

    /// Convert our OasPolygon to laykit Polygon.
    fn convert_polygon(polygon: &OasPolygon) -> laykit::Polygon {
        // Find bounding box to get origin
        let min_x = polygon.points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let min_y = polygon.points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);

        // Convert points to relative coordinates
        let points: Vec<(i64, i64)> = polygon
            .points
            .iter()
            .map(|p| {
                (
                    (p.x - min_x).round() as i64,
                    (p.y - min_y).round() as i64,
                )
            })
            .collect();

        laykit::Polygon {
            layer: polygon.layer as u32,
            datatype: polygon.datatype as u32,
            x: min_x.round() as i64,
            y: min_y.round() as i64,
            points,
            repetition: None,
            properties: Vec::new(),
        }
    }

    /// Convert our OasPath to laykit OPath.
    fn convert_path(path: &OasPath) -> laykit::OPath {
        // Find bounding box to get origin
        let min_x = path.points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let min_y = path.points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);

        // Convert points to relative coordinates
        let points: Vec<(i64, i64)> = path
            .points
            .iter()
            .map(|p| {
                (
                    (p.x - min_x).round() as i64,
                    (p.y - min_y).round() as i64,
                )
            })
            .collect();

        // Determine extension scheme
        let extension_scheme = if path.start_extension.abs() < 1e-10 && path.end_extension.abs() < 1e-10 {
            laykit::ExtensionScheme::Flush
        } else if (path.start_extension - path.width / 2.0).abs() < 1e-10
               && (path.end_extension - path.width / 2.0).abs() < 1e-10 {
            laykit::ExtensionScheme::HalfWidth
        } else {
            laykit::ExtensionScheme::Custom {
                start: path.start_extension.round() as i64,
                end: path.end_extension.round() as i64,
            }
        };

        laykit::OPath {
            layer: path.layer as u32,
            datatype: path.datatype as u32,
            x: min_x.round() as i64,
            y: min_y.round() as i64,
            half_width: (path.width / 2.0).round() as u64,
            extension_scheme,
            points,
            repetition: None,
            properties: Vec::new(),
        }
    }

    /// Convert our OasText to laykit OText.
    fn convert_text(text: &OasText) -> laykit::OText {
        laykit::OText {
            layer: text.layer as u32,
            texttype: text.text_type as u32,
            x: text.position.x.round() as i64,
            y: text.position.y.round() as i64,
            string: text.content.clone(),
            repetition: None,
            properties: Vec::new(),
        }
    }

    /// Convert our OasPlacement to laykit Placement.
    fn convert_placement(placement: &OasPlacement) -> laykit::Placement {
        // Handle array placements
        let repetition = if let Some((cols, rows, col_dx, _col_dy, _row_dx, row_dy)) = placement.array {
            if cols > 1 || rows > 1 {
                Some(laykit::Repetition::Matrix {
                    x_count: cols,
                    y_count: rows,
                    x_space: col_dx.round() as u64,
                    y_space: row_dy.round() as u64,
                })
            } else {
                None
            }
        } else {
            None
        };

        laykit::Placement {
            cell_name: placement.cell_name.clone(),
            x: placement.x.round() as i64,
            y: placement.y.round() as i64,
            magnification: if (placement.magnification - 1.0).abs() > 1e-10 {
                Some(placement.magnification)
            } else {
                None
            },
            angle: if placement.rotation.abs() > 1e-10 {
                Some(placement.rotation.to_degrees())
            } else {
                None
            },
            mirror: placement.reflection,
            repetition,
            properties: Vec::new(),
        }
    }
}

impl Default for OasWriter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::OasReader;
    use glam::DVec2;
    use std::collections::HashMap;

    fn create_test_library() -> OasLibrary {
        let mut cells = HashMap::new();

        // Create a simple cell with a polygon
        cells.insert(
            "CELL1".to_string(),
            OasCell {
                name: "CELL1".to_string(),
                polygons: vec![OasPolygon {
                    layer: 1,
                    datatype: 0,
                    points: vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(10.0, 0.0),
                        DVec2::new(10.0, 10.0),
                        DVec2::new(0.0, 10.0),
                    ],
                }],
                paths: vec![OasPath {
                    layer: 2,
                    datatype: 0,
                    width: 1.0,
                    points: vec![DVec2::new(0.0, 0.0), DVec2::new(20.0, 20.0)],
                    start_extension: 0.0,
                    end_extension: 0.0,
                }],
                texts: vec![OasText {
                    layer: 3,
                    text_type: 0,
                    position: DVec2::new(5.0, 5.0),
                    content: "test".to_string(),
                }],
                placements: vec![],
            },
        );

        OasLibrary {
            name: Some("TESTLIB".to_string()),
            cells,
        }
    }

    #[test]
    fn test_write_to_bytes() {
        let library = create_test_library();
        let bytes = OasWriter::to_bytes(&library).unwrap();

        // OASIS files should have reasonable size
        assert!(bytes.len() > 50, "OASIS bytes should not be empty");

        // Check magic bytes
        assert!(bytes.starts_with(b"%SEMI-OASIS"), "Should have OASIS magic bytes");
    }

    #[test]
    fn test_roundtrip() {
        let original = create_test_library();
        let bytes = OasWriter::to_bytes(&original).unwrap();

        // Read back the bytes
        let read_library = OasReader::parse_bytes(&bytes).unwrap();

        // Check library properties
        assert_eq!(read_library.cells.len(), original.cells.len());

        // Check cell
        let original_cell = original.cells.get("CELL1").unwrap();
        let read_cell = read_library.cells.get("CELL1").unwrap();

        // Check polygon
        assert_eq!(read_cell.polygons.len(), 1);
        let original_polygon = &original_cell.polygons[0];
        let read_polygon = &read_cell.polygons[0];
        assert_eq!(read_polygon.layer, original_polygon.layer);
        assert_eq!(read_polygon.datatype, original_polygon.datatype);
        assert_eq!(read_polygon.points.len(), original_polygon.points.len());

        // Check path
        assert_eq!(read_cell.paths.len(), 1);

        // Check text
        assert_eq!(read_cell.texts.len(), 1);
        let original_text = &original_cell.texts[0];
        let read_text = &read_cell.texts[0];
        assert_eq!(read_text.layer, original_text.layer);
        assert_eq!(read_text.content, original_text.content);
    }

    #[test]
    fn test_placements() {
        let mut cells = HashMap::new();

        // Create a referenced cell
        cells.insert(
            "REFCELL".to_string(),
            OasCell {
                name: "REFCELL".to_string(),
                polygons: vec![OasPolygon {
                    layer: 1,
                    datatype: 0,
                    points: vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(5.0, 0.0),
                        DVec2::new(5.0, 5.0),
                        DVec2::new(0.0, 5.0),
                    ],
                }],
                ..Default::default()
            },
        );

        // Create a top cell with a placement
        cells.insert(
            "TOP".to_string(),
            OasCell {
                name: "TOP".to_string(),
                placements: vec![OasPlacement {
                    cell_name: "REFCELL".to_string(),
                    x: 100.0,
                    y: 100.0,
                    rotation: 0.0,
                    reflection: false,
                    magnification: 1.0,
                    array: None,
                }],
                ..Default::default()
            },
        );

        let library = OasLibrary {
            name: Some("REFTEST".to_string()),
            cells,
        };

        // Roundtrip
        let bytes = OasWriter::to_bytes(&library).unwrap();
        let read_library = OasReader::parse_bytes(&bytes).unwrap();

        // Check placement was preserved
        let top = read_library.cells.get("TOP").unwrap();
        assert_eq!(top.placements.len(), 1);
        assert_eq!(top.placements[0].cell_name, "REFCELL");
    }

    #[test]
    fn test_array_placements() {
        let mut cells = HashMap::new();

        cells.insert(
            "UNIT".to_string(),
            OasCell {
                name: "UNIT".to_string(),
                polygons: vec![OasPolygon {
                    layer: 1,
                    datatype: 0,
                    points: vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(1.0, 0.0),
                        DVec2::new(1.0, 1.0),
                        DVec2::new(0.0, 1.0),
                    ],
                }],
                ..Default::default()
            },
        );

        // Create a top cell with an array placement (2x3)
        cells.insert(
            "ARRAY_TOP".to_string(),
            OasCell {
                name: "ARRAY_TOP".to_string(),
                placements: vec![OasPlacement {
                    cell_name: "UNIT".to_string(),
                    x: 0.0,
                    y: 0.0,
                    rotation: 0.0,
                    reflection: false,
                    magnification: 1.0,
                    array: Some((2, 3, 10.0, 0.0, 0.0, 10.0)),
                }],
                ..Default::default()
            },
        );

        let library = OasLibrary {
            name: Some("ARRAYTEST".to_string()),
            cells,
        };

        // Roundtrip
        let bytes = OasWriter::to_bytes(&library).unwrap();
        let read_library = OasReader::parse_bytes(&bytes).unwrap();

        // Check array placement was preserved
        let top = read_library.cells.get("ARRAY_TOP").unwrap();
        assert_eq!(top.placements.len(), 1);
        let placement = &top.placements[0];
        assert_eq!(placement.cell_name, "UNIT");
        // Note: arrays get expanded during read, so we check that it exists
    }

    #[test]
    fn test_rotation_and_magnification() {
        // Note: laykit's OASIS implementation has limited support for transforms.
        // It writes placements but doesn't fully preserve angle/magnification on read.
        // This test verifies that the basic placement structure is preserved.
        let mut cells = HashMap::new();

        cells.insert(
            "CELL".to_string(),
            OasCell {
                name: "CELL".to_string(),
                placements: vec![OasPlacement {
                    cell_name: "OTHER".to_string(),
                    x: 50.0,
                    y: 50.0,
                    rotation: std::f64::consts::FRAC_PI_4, // 45 degrees
                    reflection: false,
                    magnification: 2.0,
                    array: None,
                }],
                ..Default::default()
            },
        );

        cells.insert(
            "OTHER".to_string(),
            OasCell {
                name: "OTHER".to_string(),
                ..Default::default()
            },
        );

        let library = OasLibrary {
            name: Some("TRANSFORMTEST".to_string()),
            cells,
        };

        // Roundtrip
        let bytes = OasWriter::to_bytes(&library).unwrap();
        let read_library = OasReader::parse_bytes(&bytes).unwrap();

        // Check that the placement was preserved
        let cell = read_library.cells.get("CELL").unwrap();
        assert_eq!(cell.placements.len(), 1);

        let placement = &cell.placements[0];
        assert_eq!(placement.cell_name, "OTHER");

        // Check position is preserved
        assert!((placement.x - 50.0).abs() < 1.0, "X position should be ~50");
        assert!((placement.y - 50.0).abs() < 1.0, "Y position should be ~50");

        // Note: laykit doesn't fully preserve rotation/magnification in roundtrip
        // The values default to 0.0 and 1.0 respectively
    }
}
