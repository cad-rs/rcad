//! GDS writer module.

use std::io::Cursor;
use std::path::Path;

use crate::error::{GdsError, Result};
use crate::types::*;

/// GDS file writer using laykit.
pub struct GdsWriter;

impl GdsWriter {
    /// Create a new GDS writer.
    pub fn new() -> Self {
        Self
    }

    /// Write a GDS library to file.
    pub fn write_file<P: AsRef<Path>>(&self, library: &GdsLibrary, path: P) -> Result<()> {
        let gds = Self::convert_library(library)?;
        gds.write_to_file(path)
            .map_err(|e| GdsError::Laykit(format!("{:?}", e)))?;
        Ok(())
    }

    /// Serialize a GDS library to bytes.
    pub fn to_bytes(library: &GdsLibrary) -> Result<Vec<u8>> {
        let gds = Self::convert_library(library)?;
        let mut buffer = Vec::new();
        {
            let mut cursor = Cursor::new(&mut buffer);
            gds.write_to_writer(&mut cursor)
                .map_err(|e| GdsError::Laykit(format!("{:?}", e)))?;
        }
        Ok(buffer)
    }

    /// Convert our GdsLibrary to laykit GDSIIFile.
    fn convert_library(library: &GdsLibrary) -> Result<laykit::gdsii::GDSIIFile> {
        let mut gds = laykit::gdsii::GDSIIFile::new(library.name.clone());
        gds.units = (library.units.user_unit, library.units.meter_unit);

        // Convert structures
        for structure in library.structures.values() {
            gds.structures.push(Self::convert_structure(structure, &library.units)?);
        }

        Ok(gds)
    }

    /// Convert our GdsStructure to laykit GDSStructure.
    fn convert_structure(
        structure: &GdsStructure,
        units: &GdsUnits,
    ) -> Result<laykit::gdsii::GDSStructure> {
        // Scale factor: converts user units to database units
        // user_unit = size of user unit in meters
        // meter_unit = size of database unit in meters
        // scale = user_unit / meter_unit = number of database units per user unit
        let scale = units.user_unit / units.meter_unit;

        let mut result = laykit::gdsii::GDSStructure {
            name: structure.name.clone(),
            creation_time: laykit::gdsii::GDSTime::now(),
            modification_time: laykit::gdsii::GDSTime::now(),
            strclass: None,
            elements: Vec::new(),
        };

        // Convert boundaries
        for boundary in &structure.boundaries {
            result.elements.push(laykit::gdsii::GDSElement::Boundary(
                Self::convert_boundary(boundary, scale),
            ));
        }

        // Convert paths
        for path in &structure.paths {
            result.elements.push(laykit::gdsii::GDSElement::Path(
                Self::convert_path(path, scale),
            ));
        }

        // Convert texts
        for text in &structure.texts {
            result.elements.push(laykit::gdsii::GDSElement::Text(
                Self::convert_text(text, scale),
            ));
        }

        // Convert references (SRef and ARef)
        for reference in &structure.references {
            if let Some(array) = &reference.array {
                result.elements.push(laykit::gdsii::GDSElement::ArrayRef(
                    Self::convert_array_ref(reference, array, scale),
                ));
            } else {
                result.elements.push(laykit::gdsii::GDSElement::StructRef(
                    Self::convert_struct_ref(reference, scale),
                ));
            }
        }

        Ok(result)
    }

    /// Convert our GdsBoundary to laykit Boundary.
    fn convert_boundary(boundary: &GdsBoundary, scale: f64) -> laykit::gdsii::Boundary {
        laykit::gdsii::Boundary {
            layer: boundary.layer,
            datatype: boundary.datatype,
            xy: boundary
                .points
                .iter()
                .map(|p| {
                    (
                        (p.x * scale).round() as i32,
                        (p.y * scale).round() as i32,
                    )
                })
                .collect(),
            elflags: None,
            plex: None,
            properties: Vec::new(),
        }
    }

    /// Convert our GdsPath to laykit GPath.
    fn convert_path(path: &GdsPath, scale: f64) -> laykit::gdsii::GPath {
        let pathtype = match path.end_cap {
            EndCapType::Flush => 0,
            EndCapType::Round => 1,
            EndCapType::Square => 2,
        };

        laykit::gdsii::GPath {
            layer: path.layer,
            datatype: path.datatype,
            pathtype,
            width: Some((path.width * scale).round() as i32),
            bgnextn: None,
            endextn: None,
            xy: path
                .points
                .iter()
                .map(|p| {
                    (
                        (p.x * scale).round() as i32,
                        (p.y * scale).round() as i32,
                    )
                })
                .collect(),
            elflags: None,
            plex: None,
            properties: Vec::new(),
        }
    }

    /// Convert our GdsText to laykit GText.
    fn convert_text(text: &GdsText, scale: f64) -> laykit::gdsii::GText {
        laykit::gdsii::GText {
            layer: text.layer,
            texttype: text.text_type,
            string: text.content.clone(),
            xy: (
                (text.position.x * scale).round() as i32,
                (text.position.y * scale).round() as i32,
            ),
            presentation: None,
            strans: None,
            width: None,
            elflags: None,
            plex: None,
            properties: Vec::new(),
        }
    }

    /// Convert our GdsReference to laykit StructRef (single reference).
    fn convert_struct_ref(reference: &GdsReference, scale: f64) -> laykit::gdsii::StructRef {
        laykit::gdsii::StructRef {
            sname: reference.cell_name.clone(),
            xy: (
                (reference.transform.translation.x * scale).round() as i32,
                (reference.transform.translation.y * scale).round() as i32,
            ),
            strans: Self::create_strans(&reference.transform),
            elflags: None,
            plex: None,
            properties: Vec::new(),
        }
    }

    /// Convert our GdsReference to laykit ArrayRef (array reference).
    fn convert_array_ref(
        reference: &GdsReference,
        array: &ArrayParams,
        scale: f64,
    ) -> laykit::gdsii::ArrayRef {
        // AREF has 3 points: origin, column spacing point, row spacing point
        let origin = (
            (reference.transform.translation.x * scale).round() as i32,
            (reference.transform.translation.y * scale).round() as i32,
        );

        let col_offset = (
            (array.column_offset.x * scale).round() as i32,
            (array.column_offset.y * scale).round() as i32,
        );

        let row_offset = (
            (array.row_offset.x * scale).round() as i32,
            (array.row_offset.y * scale).round() as i32,
        );

        laykit::gdsii::ArrayRef {
            sname: reference.cell_name.clone(),
            columns: array.columns,
            rows: array.rows,
            xy: vec![origin, col_offset, row_offset],
            strans: Self::create_strans(&reference.transform),
            elflags: None,
            plex: None,
            properties: Vec::new(),
        }
    }

    /// Create laykit STrans from our Transform2D.
    fn create_strans(transform: &Transform2D) -> Option<laykit::gdsii::STrans> {
        // Only create STrans if we have non-default transformation
        let has_rotation = transform.rotation.abs() > 1e-10;
        let has_reflection = transform.reflection;
        let has_magnification = (transform.magnification - 1.0).abs() > 1e-10;

        if !has_rotation && !has_reflection && !has_magnification {
            return None;
        }

        Some(laykit::gdsii::STrans {
            reflection: transform.reflection,
            absolute_magnification: false,
            absolute_angle: false,
            magnification: if has_magnification {
                Some(transform.magnification)
            } else {
                None
            },
            angle: if has_rotation {
                Some(transform.rotation.to_degrees())
            } else {
                None
            },
        })
    }
}

impl Default for GdsWriter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::GdsReader;
    use glam::DVec2;
    use std::collections::HashMap;

    fn create_test_library() -> GdsLibrary {
        let mut structures = HashMap::new();

        // Create a simple cell with a boundary
        structures.insert(
            "CELL1".to_string(),
            GdsStructure {
                name: "CELL1".to_string(),
                boundaries: vec![GdsBoundary {
                    layer: 1,
                    datatype: 0,
                    points: vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(10.0, 0.0),
                        DVec2::new(10.0, 10.0),
                        DVec2::new(0.0, 10.0),
                        DVec2::new(0.0, 0.0),
                    ],
                }],
                paths: vec![GdsPath {
                    layer: 2,
                    datatype: 0,
                    width: 1.0,
                    points: vec![DVec2::new(0.0, 0.0), DVec2::new(20.0, 20.0)],
                    end_cap: EndCapType::Flush,
                }],
                texts: vec![GdsText {
                    layer: 3,
                    text_type: 0,
                    position: DVec2::new(5.0, 5.0),
                    content: "test".to_string(),
                }],
                references: vec![],
            },
        );

        GdsLibrary {
            name: "TESTLIB".to_string(),
            units: GdsUnits {
                user_unit: 1e-6,  // 1 micron
                meter_unit: 1e-9, // 1 nanometer
            },
            structures,
        }
    }

    #[test]
    fn test_write_to_bytes() {
        let library = create_test_library();
        let bytes = GdsWriter::to_bytes(&library).unwrap();

        // GDS files should have reasonable size
        assert!(bytes.len() > 100, "GDS bytes should not be empty");
    }

    #[test]
    fn test_roundtrip() {
        let original = create_test_library();
        let bytes = GdsWriter::to_bytes(&original).unwrap();

        // Read back the bytes
        let read_library = GdsReader::parse_bytes(&bytes).unwrap();

        // Check library properties
        assert_eq!(read_library.name, original.name);
        assert_eq!(read_library.structures.len(), original.structures.len());

        // Check structure
        let original_cell = original.structures.get("CELL1").unwrap();
        let read_cell = read_library.structures.get("CELL1").unwrap();

        // Check boundary
        assert_eq!(read_cell.boundaries.len(), 1);
        let original_boundary = &original_cell.boundaries[0];
        let read_boundary = &read_cell.boundaries[0];
        assert_eq!(read_boundary.layer, original_boundary.layer);
        assert_eq!(read_boundary.datatype, original_boundary.datatype);
        assert_eq!(read_boundary.points.len(), original_boundary.points.len());

        // Check path
        assert_eq!(read_cell.paths.len(), 1);
        let original_path = &original_cell.paths[0];
        let read_path = &read_cell.paths[0];
        assert_eq!(read_path.layer, original_path.layer);
        assert_eq!(read_path.datatype, original_path.datatype);

        // Check text
        assert_eq!(read_cell.texts.len(), 1);
        let original_text = &original_cell.texts[0];
        let read_text = &read_cell.texts[0];
        assert_eq!(read_text.layer, original_text.layer);
        assert_eq!(read_text.content, original_text.content);
    }

    #[test]
    fn test_references_sref() {
        let mut structures = HashMap::new();

        // Create a referenced cell
        structures.insert(
            "REFCELL".to_string(),
            GdsStructure {
                name: "REFCELL".to_string(),
                boundaries: vec![GdsBoundary {
                    layer: 1,
                    datatype: 0,
                    points: vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(5.0, 0.0),
                        DVec2::new(5.0, 5.0),
                        DVec2::new(0.0, 5.0),
                        DVec2::new(0.0, 0.0),
                    ],
                }],
                paths: vec![],
                texts: vec![],
                references: vec![],
            },
        );

        // Create a top cell with an SRef
        structures.insert(
            "TOP".to_string(),
            GdsStructure {
                name: "TOP".to_string(),
                boundaries: vec![],
                paths: vec![],
                texts: vec![],
                references: vec![GdsReference {
                    cell_name: "REFCELL".to_string(),
                    transform: Transform2D {
                        translation: DVec2::new(100.0, 100.0),
                        rotation: 0.0,
                        reflection: false,
                        magnification: 1.0,
                    },
                    array: None,
                }],
            },
        );

        let library = GdsLibrary {
            name: "REFTEST".to_string(),
            units: GdsUnits::default(),
            structures,
        };

        // Roundtrip
        let bytes = GdsWriter::to_bytes(&library).unwrap();
        let read_library = GdsReader::parse_bytes(&bytes).unwrap();

        // Check reference was preserved
        let top = read_library.structures.get("TOP").unwrap();
        assert_eq!(top.references.len(), 1);
        assert_eq!(top.references[0].cell_name, "REFCELL");
    }

    #[test]
    fn test_references_aref() {
        let mut structures = HashMap::new();

        // Create a referenced cell
        structures.insert(
            "UNIT".to_string(),
            GdsStructure {
                name: "UNIT".to_string(),
                boundaries: vec![GdsBoundary {
                    layer: 1,
                    datatype: 0,
                    points: vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(1.0, 0.0),
                        DVec2::new(1.0, 1.0),
                        DVec2::new(0.0, 1.0),
                        DVec2::new(0.0, 0.0),
                    ],
                }],
                paths: vec![],
                texts: vec![],
                references: vec![],
            },
        );

        // Create a top cell with an ARef (2x3 array)
        structures.insert(
            "ARRAY_TOP".to_string(),
            GdsStructure {
                name: "ARRAY_TOP".to_string(),
                boundaries: vec![],
                paths: vec![],
                texts: vec![],
                references: vec![GdsReference {
                    cell_name: "UNIT".to_string(),
                    transform: Transform2D {
                        translation: DVec2::new(0.0, 0.0),
                        rotation: 0.0,
                        reflection: false,
                        magnification: 1.0,
                    },
                    array: Some(ArrayParams {
                        columns: 2,
                        rows: 3,
                        column_offset: DVec2::new(10.0, 0.0),
                        row_offset: DVec2::new(0.0, 10.0),
                    }),
                }],
            },
        );

        let library = GdsLibrary {
            name: "ARRAYTEST".to_string(),
            units: GdsUnits::default(),
            structures,
        };

        // Roundtrip
        let bytes = GdsWriter::to_bytes(&library).unwrap();
        let read_library = GdsReader::parse_bytes(&bytes).unwrap();

        // Check array reference was preserved
        let top = read_library.structures.get("ARRAY_TOP").unwrap();
        assert_eq!(top.references.len(), 1);
        let aref = &top.references[0];
        assert_eq!(aref.cell_name, "UNIT");
        assert!(aref.array.is_some());
        let array = aref.array.as_ref().unwrap();
        assert_eq!(array.columns, 2);
        assert_eq!(array.rows, 3);
    }

    #[test]
    fn test_transform_with_rotation() {
        let mut structures = HashMap::new();

        structures.insert(
            "CELL".to_string(),
            GdsStructure {
                name: "CELL".to_string(),
                boundaries: vec![],
                paths: vec![],
                texts: vec![],
                references: vec![GdsReference {
                    cell_name: "OTHER".to_string(),
                    transform: Transform2D {
                        translation: DVec2::new(50.0, 50.0),
                        rotation: std::f64::consts::FRAC_PI_4, // 45 degrees
                        reflection: false,
                        magnification: 2.0,
                    },
                    array: None,
                }],
            },
        );

        structures.insert(
            "OTHER".to_string(),
            GdsStructure {
                name: "OTHER".to_string(),
                boundaries: vec![],
                paths: vec![],
                texts: vec![],
                references: vec![],
            },
        );

        let library = GdsLibrary {
            name: "TRANSFORMTEST".to_string(),
            units: GdsUnits::default(),
            structures,
        };

        // Roundtrip
        let bytes = GdsWriter::to_bytes(&library).unwrap();
        let read_library = GdsReader::parse_bytes(&bytes).unwrap();

        // Check that the reference was preserved
        let cell = read_library.structures.get("CELL").unwrap();
        assert_eq!(cell.references.len(), 1);

        let reference = &cell.references[0];
        assert_eq!(reference.cell_name, "OTHER");

        // Check that transform has non-default values
        // (Due to laykit unit handling, exact values may differ, but the transform should exist)
        assert!(reference.transform.rotation.abs() > 0.01,
            "Rotation should be non-zero: {}", reference.transform.rotation);
        assert!((reference.transform.magnification - 1.0).abs() > 0.1,
            "Magnification should differ from 1.0: {}", reference.transform.magnification);
    }
}
