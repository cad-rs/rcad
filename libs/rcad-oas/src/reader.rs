//! OASIS file reader.

use std::collections::HashMap;
use std::io::Cursor;
use std::path::Path;

use crate::{OasCell, OasError, OasLibrary, OasPath, OasPlacement, OasPolygon, OasText};

/// OASIS file reader using laykit.
pub struct OasReader;

impl OasReader {
    /// Read an OASIS file from disk.
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<OasLibrary, OasError> {
        let oasis = laykit::OASISFile::read_from_file(path)
            .map_err(|e| OasError::Laykit(format!("{:?}", e)))?;

        Self::convert_library(&oasis)
    }

    /// Parse OASIS data from bytes.
    pub fn parse_bytes(data: &[u8]) -> Result<OasLibrary, OasError> {
        let mut cursor = Cursor::new(data);
        let oasis = laykit::OASISFile::read_from_reader(&mut cursor)
            .map_err(|e| OasError::Laykit(format!("{:?}", e)))?;

        Self::convert_library(&oasis)
    }

    /// Convert laykit OASISFile to our OasLibrary.
    fn convert_library(oasis: &laykit::OASISFile) -> Result<OasLibrary, OasError> {
        let mut cells = HashMap::new();

        for cell in &oasis.cells {
            let converted = Self::convert_cell(cell)?;
            cells.insert(converted.name.clone(), converted);
        }

        Ok(OasLibrary {
            name: None, // OASIS doesn't have a library name like GDSII
            cells,
        })
    }

    /// Convert a laykit OASISCell to our OasCell.
    fn convert_cell(cell: &laykit::OASISCell) -> Result<OasCell, OasError> {
        let mut result = OasCell {
            name: cell.name.clone(),
            ..Default::default()
        };

        for element in &cell.elements {
            match element {
                laykit::OASISElement::Polygon(polygon) => {
                    // Expand repetitions
                    for expanded in Self::expand_polygon(polygon) {
                        result.polygons.push(expanded);
                    }
                }
                laykit::OASISElement::Path(path) => {
                    for expanded in Self::expand_path(path) {
                        result.paths.push(expanded);
                    }
                }
                laykit::OASISElement::Text(text) => {
                    for expanded in Self::expand_text(text) {
                        result.texts.push(expanded);
                    }
                }
                laykit::OASISElement::Placement(placement) => {
                    for expanded in Self::expand_placement(placement) {
                        result.placements.push(expanded);
                    }
                }
                laykit::OASISElement::Rectangle(rect) => {
                    // Convert rectangle to polygon
                    for expanded in Self::expand_rectangle(rect) {
                        result.polygons.push(expanded);
                    }
                }
                laykit::OASISElement::Trapezoid(trap) => {
                    // Convert trapezoid to polygon
                    for expanded in Self::expand_trapezoid(trap) {
                        result.polygons.push(expanded);
                    }
                }
                laykit::OASISElement::CTrapezoid(ctrap) => {
                    // Convert ctrapezoid to polygon
                    for expanded in Self::expand_ctrapezoid(ctrap) {
                        result.polygons.push(expanded);
                    }
                }
                laykit::OASISElement::Circle(circle) => {
                    // Convert circle to polygon (approximation)
                    for expanded in Self::expand_circle(circle) {
                        result.polygons.push(expanded);
                    }
                }
            }
        }

        Ok(result)
    }

    /// Expand polygon with repetition.
    fn expand_polygon(polygon: &laykit::Polygon) -> Vec<OasPolygon> {
        let mut result = Vec::new();

        let base_points: Vec<glam::DVec2> = polygon
            .points
            .iter()
            .map(|&(x, y)| {
                glam::DVec2::new(polygon.x as f64 + x as f64, polygon.y as f64 + y as f64)
            })
            .collect();

        let offsets = Self::get_repetition_offsets(&polygon.repetition);

        for (dx, dy) in offsets {
            let points: Vec<glam::DVec2> = base_points
                .iter()
                .map(|p| glam::DVec2::new(p.x + dx, p.y + dy))
                .collect();

            result.push(OasPolygon {
                layer: polygon.layer as i32,
                datatype: polygon.datatype as i32,
                points,
            });
        }

        result
    }

    /// Expand path with repetition.
    fn expand_path(path: &laykit::OPath) -> Vec<OasPath> {
        let mut result = Vec::new();

        let base_points: Vec<glam::DVec2> = path
            .points
            .iter()
            .map(|&(x, y)| glam::DVec2::new(path.x as f64 + x as f64, path.y as f64 + y as f64))
            .collect();

        let width = (path.half_width * 2) as f64;
        let (start_ext, end_ext) = match &path.extension_scheme {
            laykit::ExtensionScheme::Flush => (0.0, 0.0),
            laykit::ExtensionScheme::HalfWidth => (width / 2.0, width / 2.0),
            laykit::ExtensionScheme::Custom { start, end } => (*start as f64, *end as f64),
        };

        let offsets = Self::get_repetition_offsets(&path.repetition);

        for (dx, dy) in offsets {
            let points: Vec<glam::DVec2> = base_points
                .iter()
                .map(|p| glam::DVec2::new(p.x + dx, p.y + dy))
                .collect();

            result.push(OasPath {
                layer: path.layer as i32,
                datatype: path.datatype as i32,
                width,
                points,
                start_extension: start_ext,
                end_extension: end_ext,
            });
        }

        result
    }

    /// Expand text with repetition.
    fn expand_text(text: &laykit::OText) -> Vec<OasText> {
        let mut result = Vec::new();

        let base_pos = glam::DVec2::new(text.x as f64, text.y as f64);
        let offsets = Self::get_repetition_offsets(&text.repetition);

        for (dx, dy) in offsets {
            result.push(OasText {
                layer: text.layer as i32,
                text_type: text.texttype as i32,
                position: glam::DVec2::new(base_pos.x + dx, base_pos.y + dy),
                content: text.string.clone(),
            });
        }

        result
    }

    /// Expand placement with repetition.
    fn expand_placement(placement: &laykit::Placement) -> Vec<OasPlacement> {
        let mut result = Vec::new();

        let offsets = Self::get_repetition_offsets(&placement.repetition);

        for (dx, dy) in offsets {
            let x = placement.x as f64 + dx;
            let y = placement.y as f64 + dy;

            // Convert repetition to array format if it's a matrix
            let array = if let Some(laykit::Repetition::Matrix {
                x_count,
                y_count,
                x_space,
                y_space,
            }) = &placement.repetition
            {
                if *x_count > 1 || *y_count > 1 {
                    Some((
                        *x_count,
                        *y_count,
                        *x_space as f64,
                        0.0,
                        0.0,
                        *y_space as f64,
                    ))
                } else {
                    None
                }
            } else {
                None
            };

            result.push(OasPlacement {
                cell_name: placement.cell_name.clone(),
                x,
                y,
                rotation: placement.angle.unwrap_or(0.0).to_radians(),
                reflection: placement.mirror,
                magnification: placement.magnification.unwrap_or(1.0),
                array,
            });
        }

        result
    }

    /// Expand rectangle to polygon(s).
    fn expand_rectangle(rect: &laykit::Rectangle) -> Vec<OasPolygon> {
        let mut result = Vec::new();

        let offsets = Self::get_repetition_offsets(&rect.repetition);

        for (dx, dy) in offsets {
            let x = rect.x as f64 + dx;
            let y = rect.y as f64 + dy;
            let w = rect.width as f64;
            let h = rect.height as f64;

            let points = vec![
                glam::DVec2::new(x, y),
                glam::DVec2::new(x + w, y),
                glam::DVec2::new(x + w, y + h),
                glam::DVec2::new(x, y + h),
            ];

            result.push(OasPolygon {
                layer: rect.layer as i32,
                datatype: rect.datatype as i32,
                points,
            });
        }

        result
    }

    /// Expand trapezoid to polygon(s).
    fn expand_trapezoid(trap: &laykit::Trapezoid) -> Vec<OasPolygon> {
        let mut result = Vec::new();

        let offsets = Self::get_repetition_offsets(&trap.repetition);

        for (dx, dy) in offsets {
            let x = trap.x as f64 + dx;
            let y = trap.y as f64 + dy;
            let w = trap.width as f64;
            let h = trap.height as f64;
            let da = trap.delta_a as f64;
            let db = trap.delta_b as f64;

            let points = if trap.orientation {
                // Horizontal trapezoid
                vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y),
                    glam::DVec2::new(x + w + db, y + h),
                    glam::DVec2::new(x + da, y + h),
                ]
            } else {
                // Vertical trapezoid
                vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y + da),
                    glam::DVec2::new(x + w, y + h + da),
                    glam::DVec2::new(x, y + h),
                ]
            };

            result.push(OasPolygon {
                layer: trap.layer as i32,
                datatype: trap.datatype as i32,
                points,
            });
        }

        result
    }

    /// Expand ctrapezoid to polygon(s).
    fn expand_ctrapezoid(ctrap: &laykit::CTrapezoid) -> Vec<OasPolygon> {
        let mut result = Vec::new();

        let offsets = Self::get_repetition_offsets(&ctrap.repetition);

        for (dx, dy) in offsets {
            let x = ctrap.x as f64 + dx;
            let y = ctrap.y as f64 + dy;
            let w = ctrap.width as f64;
            let h = ctrap.height as f64;

            // CTRAPEZOID types 0-7 are half-width triangles, 8-23 are various trapezoids
            let points = match ctrap.trap_type {
                0 => vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y),
                    glam::DVec2::new(x + w, y + h),
                ],
                1 => vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y),
                    glam::DVec2::new(x, y + h),
                ],
                2 => vec![
                    glam::DVec2::new(x, y + h),
                    glam::DVec2::new(x + w, y),
                    glam::DVec2::new(x + w, y + h),
                ],
                3 => vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x, y + h),
                    glam::DVec2::new(x + w, y + h),
                ],
                4 => vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y + h),
                    glam::DVec2::new(x, y + h),
                ],
                5 => vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y),
                    glam::DVec2::new(x, y + h),
                ],
                6 => vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y + h),
                    glam::DVec2::new(x + w, y),
                ],
                7 => vec![
                    glam::DVec2::new(x, y + h),
                    glam::DVec2::new(x + w, y),
                    glam::DVec2::new(x, y),
                ],
                // For types 8+, we default to a simple rectangle approximation
                _ => vec![
                    glam::DVec2::new(x, y),
                    glam::DVec2::new(x + w, y),
                    glam::DVec2::new(x + w, y + h),
                    glam::DVec2::new(x, y + h),
                ],
            };

            result.push(OasPolygon {
                layer: ctrap.layer as i32,
                datatype: ctrap.datatype as i32,
                points,
            });
        }

        result
    }

    /// Expand circle to polygon approximation.
    fn expand_circle(circle: &laykit::Circle) -> Vec<OasPolygon> {
        let mut result = Vec::new();

        let offsets = Self::get_repetition_offsets(&circle.repetition);

        // Approximate circle with 32 segments
        let num_segments = 32;
        let radius = circle.radius as f64;

        for (dx, dy) in offsets {
            let cx = circle.x as f64 + dx;
            let cy = circle.y as f64 + dy;

            let points: Vec<glam::DVec2> = (0..num_segments)
                .map(|i| {
                    let angle = 2.0 * std::f64::consts::PI * i as f64 / num_segments as f64;
                    glam::DVec2::new(cx + radius * angle.cos(), cy + radius * angle.sin())
                })
                .collect();

            result.push(OasPolygon {
                layer: circle.layer as i32,
                datatype: circle.datatype as i32,
                points,
            });
        }

        result
    }

    /// Get offsets from repetition pattern.
    fn get_repetition_offsets(repetition: &Option<laykit::Repetition>) -> Vec<(f64, f64)> {
        match repetition {
            None => vec![(0.0, 0.0)],
            Some(laykit::Repetition::ReusePrevious) => vec![(0.0, 0.0)], // Would need to track previous
            Some(laykit::Repetition::Matrix {
                x_count,
                y_count,
                x_space,
                y_space,
            }) => {
                let mut offsets = Vec::new();
                for iy in 0..*y_count {
                    for ix in 0..*x_count {
                        offsets.push((ix as f64 * *x_space as f64, iy as f64 * *y_space as f64));
                    }
                }
                offsets
            }
            Some(laykit::Repetition::Arbitrary {
                x_displacements,
                y_displacements,
            }) => {
                let mut offsets = vec![(0.0, 0.0)];
                let mut x_acc = 0i64;
                let mut y_acc = 0i64;

                for &dx in x_displacements {
                    x_acc += dx;
                    offsets.push((x_acc as f64, y_acc as f64));
                }

                for &dy in y_displacements {
                    y_acc += dy;
                    offsets.push((x_acc as f64, y_acc as f64));
                }

                offsets
            }
            Some(laykit::Repetition::Grid { count, grid_space }) => {
                let mut offsets = Vec::new();
                for i in 0..*count {
                    offsets.push((i as f64 * *grid_space as f64, 0.0));
                }
                offsets
            }
        }
    }
}

impl OasLibrary {
    /// Get list of top-level cells (cells not referenced by any other cell).
    pub fn top_cells(&self) -> Vec<&str> {
        let mut referenced: std::collections::HashSet<&str> = std::collections::HashSet::new();

        for cell in self.cells.values() {
            for placement in &cell.placements {
                referenced.insert(&placement.cell_name);
            }
        }

        self.cells
            .keys()
            .filter(|name| !referenced.contains(name.as_str()))
            .map(|s| s.as_str())
            .collect()
    }

    /// Check if a cell exists.
    pub fn has_cell(&self, name: &str) -> bool {
        self.cells.contains_key(name)
    }
}

impl Default for OasReader {
    fn default() -> Self {
        Self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_library_structure() {
        let mut library = OasLibrary {
            name: Some("test".to_string()),
            cells: HashMap::new(),
        };

        library.cells.insert(
            "TOP".to_string(),
            OasCell {
                name: "TOP".to_string(),
                polygons: vec![OasPolygon {
                    layer: 1,
                    datatype: 0,
                    points: vec![
                        glam::DVec2::new(0.0, 0.0),
                        glam::DVec2::new(10.0, 0.0),
                        glam::DVec2::new(10.0, 10.0),
                        glam::DVec2::new(0.0, 10.0),
                    ],
                }],
                ..Default::default()
            },
        );

        assert!(library.has_cell("TOP"));
        assert_eq!(library.top_cells(), vec!["TOP"]);
    }

    #[test]
    fn test_top_cells_with_references() {
        let mut library = OasLibrary {
            name: Some("test".to_string()),
            cells: HashMap::new(),
        };

        library.cells.insert(
            "TOP".to_string(),
            OasCell {
                name: "TOP".to_string(),
                placements: vec![OasPlacement {
                    cell_name: "CELL_A".to_string(),
                    x: 0.0,
                    y: 0.0,
                    rotation: 0.0,
                    reflection: false,
                    magnification: 1.0,
                    array: None,
                }],
                ..Default::default()
            },
        );

        library.cells.insert(
            "CELL_A".to_string(),
            OasCell {
                name: "CELL_A".to_string(),
                ..Default::default()
            },
        );

        let top = library.top_cells();
        assert_eq!(top.len(), 1);
        assert_eq!(top[0], "TOP");
    }
}
