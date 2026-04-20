//! GDS reader module.

use std::path::Path;

use crate::error::{GdsError, Result};
use crate::types::*;

/// GDS file reader using laykit.
pub struct GdsReader;

impl GdsReader {
    /// Read a GDS file from disk.
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<GdsLibrary> {
        let gds_file = laykit::gdsii::GDSIIFile::read_from_file(path)
            .map_err(|e| GdsError::Laykit(format!("{:?}", e)))?;

        Self::convert_library(&gds_file)
    }

    /// Parse GDS data from bytes.
    pub fn parse_bytes(data: &[u8]) -> Result<GdsLibrary> {
        use std::io::Cursor;
        let mut cursor = Cursor::new(data);
        let gds_file = laykit::gdsii::GDSIIFile::read_from_reader(&mut cursor)
            .map_err(|e| GdsError::Laykit(format!("{:?}", e)))?;

        Self::convert_library(&gds_file)
    }

    /// Convert laykit GDSIIFile to our GdsLibrary.
    fn convert_library(gds: &laykit::gdsii::GDSIIFile) -> Result<GdsLibrary> {
        let name = gds.library_name.clone();

        // Extract units
        let units = GdsUnits {
            user_unit: gds.units.0,
            meter_unit: gds.units.1,
        };

        // Convert structures
        let mut structures = std::collections::HashMap::new();
        for structure in &gds.structures {
            let converted = Self::convert_structure(structure, &units)?;
            structures.insert(structure.name.clone(), converted);
        }

        Ok(GdsLibrary {
            name,
            units,
            structures,
        })
    }

    /// Convert a laykit structure to our GdsStructure.
    fn convert_structure(
        structure: &laykit::gdsii::GDSStructure,
        units: &GdsUnits,
    ) -> Result<GdsStructure> {
        let mut result = GdsStructure {
            name: structure.name.clone(),
            ..Default::default()
        };

        let scale = units.user_unit / units.meter_unit;

        for element in &structure.elements {
            match element {
                laykit::gdsii::GDSElement::Boundary(boundary) => {
                    let points: Vec<glam::DVec2> = boundary.xy.iter()
                        .map(|&(x, y)| glam::DVec2::new(x as f64 * scale, y as f64 * scale))
                        .collect();
                    result.boundaries.push(GdsBoundary {
                        layer: boundary.layer,
                        datatype: boundary.datatype,
                        points,
                    });
                }
                laykit::gdsii::GDSElement::Path(path) => {
                    let points: Vec<glam::DVec2> = path.xy.iter()
                        .map(|&(x, y)| glam::DVec2::new(x as f64 * scale, y as f64 * scale))
                        .collect();
                    let width = path.width.unwrap_or(0) as f64 * scale;
                    result.paths.push(GdsPath {
                        layer: path.layer,
                        datatype: path.datatype,
                        width,
                        points,
                        end_cap: EndCapType::default(),
                    });
                }
                laykit::gdsii::GDSElement::Text(text) => {
                    result.texts.push(GdsText {
                        layer: text.layer,
                        text_type: text.texttype,
                        position: glam::DVec2::new(
                            text.xy.0 as f64 * scale,
                            text.xy.1 as f64 * scale,
                        ),
                        content: text.string.clone(),
                    });
                }
                laykit::gdsii::GDSElement::StructRef(sref) => {
                    let transform = Self::extract_transform(&sref.strans, sref.xy, scale);
                    result.references.push(GdsReference {
                        cell_name: sref.sname.clone(),
                        transform,
                        array: None,
                    });
                }
                laykit::gdsii::GDSElement::ArrayRef(aref) => {
                    let transform = Self::extract_transform(&aref.strans, aref.xy.first().copied().unwrap_or((0, 0)), scale);

                    // In GDS, AREF has 3 points: origin, column spacing point, row spacing point
                    // The column/row offsets are calculated from these points
                    let col_offset = aref.xy.get(1).map(|&p| p).unwrap_or((0, 0));
                    let row_offset = aref.xy.get(2).map(|&p| p).unwrap_or((0, 0));

                    result.references.push(GdsReference {
                        cell_name: aref.sname.clone(),
                        transform,
                        array: Some(ArrayParams {
                            columns: aref.columns,
                            rows: aref.rows,
                            column_offset: glam::DVec2::new(
                                col_offset.0 as f64 * scale,
                                col_offset.1 as f64 * scale,
                            ),
                            row_offset: glam::DVec2::new(
                                row_offset.0 as f64 * scale,
                                row_offset.1 as f64 * scale,
                            ),
                        }),
                    });
                }
                // Skip other element types (Node, Box)
                _ => {}
            }
        }

        Ok(result)
    }

    /// Extract transform from laykit STrans.
    fn extract_transform(
        strans: &Option<laykit::gdsii::STrans>,
        xy: (i32, i32),
        scale: f64,
    ) -> Transform2D {
        match strans {
            Some(st) => Transform2D {
                translation: glam::DVec2::new(xy.0 as f64 * scale, xy.1 as f64 * scale),
                rotation: st.angle.unwrap_or(0.0).to_radians(),
                reflection: st.reflection,
                magnification: st.magnification.unwrap_or(1.0),
            },
            None => Transform2D::from_translation(xy.0 as f64 * scale, xy.1 as f64 * scale),
        }
    }
}

impl GdsLibrary {
    /// Get list of top-level cells (cells not referenced by any other cell).
    pub fn top_cells(&self) -> Vec<&str> {
        let mut referenced: std::collections::HashSet<&str> = std::collections::HashSet::new();

        for structure in self.structures.values() {
            for reference in &structure.references {
                referenced.insert(&reference.cell_name);
            }
        }

        self.structures.keys()
            .filter(|name| !referenced.contains(name.as_str()))
            .map(|s| s.as_str())
            .collect()
    }

    /// Check if a cell exists.
    pub fn has_cell(&self, name: &str) -> bool {
        self.structures.contains_key(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_simple_gds() {
        // Create a minimal GDS file in memory for testing
        // This is a placeholder - in real implementation, we'd use actual test files
        // For now, we test the structure creation directly
        let mut library = GdsLibrary {
            name: "test".to_string(),
            units: GdsUnits::default(),
            structures: std::collections::HashMap::new(),
        };

        library.structures.insert("TOP".to_string(), GdsStructure {
            name: "TOP".to_string(),
            boundaries: vec![GdsBoundary {
                layer: 1,
                datatype: 0,
                points: vec![
                    glam::DVec2::new(0.0, 0.0),
                    glam::DVec2::new(10.0, 0.0),
                    glam::DVec2::new(10.0, 10.0),
                    glam::DVec2::new(0.0, 10.0),
                    glam::DVec2::new(0.0, 0.0),
                ],
            }],
            ..Default::default()
        });

        assert!(library.has_cell("TOP"));
        assert_eq!(library.top_cells(), vec!["TOP"]);
    }

    #[test]
    fn test_top_cells_with_references() {
        let mut library = GdsLibrary {
            name: "test".to_string(),
            units: GdsUnits::default(),
            structures: std::collections::HashMap::new(),
        };

        library.structures.insert("TOP".to_string(), GdsStructure {
            name: "TOP".to_string(),
            references: vec![GdsReference {
                cell_name: "CELL_A".to_string(),
                transform: Transform2D::default(),
                array: None,
            }],
            ..Default::default()
        });

        library.structures.insert("CELL_A".to_string(), GdsStructure {
            name: "CELL_A".to_string(),
            ..Default::default()
        });

        let top = library.top_cells();
        assert_eq!(top.len(), 1);
        assert_eq!(top[0], "TOP");
    }
}
