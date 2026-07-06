//! Conversion between OASIS geometry and RCAD kernel types.

use rcad_kernel::{BRep, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

use crate::{OasCell, OasError, OasLibrary, OasPath, OasPolygon};
use crate::layer_config::LayerConfig;

/// Build a flat-faced BRep from a polygon at a given Z offset.
fn polygon_to_brep_face(
    points: &[glam::DVec2],
    z_offset: f64,
) -> Result<BRep, OasError> {
    if points.len() < 3 {
        return Err(OasError::GeometryError(
            "Polygon has fewer than 3 points".to_string(),
        ));
    }

    let mut brep = BRep::default();
    let n_edges = points.len();

    // Vertices at z_offset
    let vi_start = brep.vertices.len();
    for pt in points.iter().take(n_edges) {
        brep.vertices.push(Vertex {
            point: glam::DVec3::new(pt.x, pt.y, z_offset),
        });
    }

    // Edges
    let ei_start = brep.edges.len();
    for i in 0..n_edges {
        brep.edges.push(Edge {
            start: vi_start + i,
            end: vi_start + (i + 1) % n_edges,
        });
    }

    // Wire
    let wire = Wire {
        edges: (0..n_edges)
            .map(|i| WireEdge::fwd(ei_start + i))
            .collect(),
    };

    let face = Face {
        outer_wire: wire,
        inner_wires: Vec::new(),
        normal: glam::DVec3::Z,
        triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };

        let shell = Shell { faces: vec![face] };
        brep.solids.push(Solid { shells: vec![shell] });

    Ok(brep)
}

/// Merge all vertices, edges, and solids from `src` into `dst`,
/// offsetting vertex and edge indices.
fn merge_into(dst: &mut BRep, src: &BRep) {
    let _v_off = dst.vertices.len();
    let e_off = dst.edges.len();

    dst.vertices.extend_from_slice(&src.vertices);
    dst.edges.extend_from_slice(&src.edges);

    for solid in &src.solids {
        let mut new_solid = Solid { shells: Vec::new() };
        for shell in &solid.shells {
            let mut new_shell = Shell { faces: Vec::new() };
            for face in &shell.faces {
                let mut new_face = face.clone();
                let outer = Wire {
                    edges: new_face.outer_wire.edges.iter()
                        .map(|we| WireEdge {
                            idx: we.idx + e_off,
                            forward: we.forward,
                        })
                        .collect(),
                };
                let inner: Vec<Wire> = new_face.inner_wires.iter().map(|w| Wire {
                    edges: w.edges.iter()
                        .map(|we| WireEdge {
                            idx: we.idx + e_off,
                            forward: we.forward,
                        })
                        .collect(),
                }).collect();
                new_face.outer_wire = outer;
                new_face.inner_wires = inner;
                new_shell.faces.push(new_face);
            }
            new_solid.shells.push(new_shell);
        }
        dst.solids.push(new_solid);
    }
}

/// Convert OASIS geometry to RCAD kernel shapes.
pub struct OasConverter {
    config: LayerConfig,
}

impl OasConverter {
    /// Create a new converter with layer configuration.
    pub fn new(config: LayerConfig) -> Self {
        Self { config }
    }

    /// Convert an OASIS polygon to a BRep face.
    pub fn polygon_to_face(&self, polygon: &OasPolygon) -> Result<Face, OasError> {
        if polygon.points.len() < 3 {
            return Err(OasError::GeometryError("Polygon has fewer than 3 points".to_string()));
        }

        let layer_settings = self.config.get(polygon.layer);

        // Create vertices for this polygon (used for reference)
        let _face_vertices: Vec<Vertex> = polygon
            .points
            .iter()
            .map(|p| Vertex {
                point: glam::DVec3::new(p.x, p.y, layer_settings.z_offset),
            })
            .collect();

        // Create edges
        let mut face_edges: Vec<Edge> = Vec::new();
        let mut wire_edges: Vec<WireEdge> = Vec::new();

        for i in 0..polygon.points.len() {
            let next_idx = (i + 1) % polygon.points.len();
            face_edges.push(Edge {
                start: i,
                end: next_idx,
            });
            wire_edges.push(WireEdge::fwd(face_edges.len() - 1));
        }

        let wire = Wire {
            edges: wire_edges,
        };

        let face = Face {
            outer_wire: wire,
            inner_wires: Vec::new(),
            normal: glam::DVec3::Z,
            triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };

        Ok(face)
    }

    /// Convert an OASIS path to a BRep face (as polygon).
    pub fn path_to_face(&self, path: &OasPath) -> Result<Face, OasError> {
        if path.points.len() < 2 {
            return Err(OasError::GeometryError("Path has fewer than 2 points".to_string()));
        }

        // Convert path to polygon by offsetting the centerline
        let polygon = Self::path_to_polygon(path);

        self.polygon_to_face(&polygon)
    }

    /// Convert a path to a polygon by offsetting the centerline.
    fn path_to_polygon(path: &OasPath) -> OasPolygon {
        let half_width = path.width / 2.0;

        if path.points.len() == 2 {
            // Simple case: straight line
            let p0 = path.points[0];
            let p1 = path.points[1];

            let dx = p1.x - p0.x;
            let dy = p1.y - p0.y;
            let len = (dx * dx + dy * dy).sqrt();

            if len < 1e-10 {
                return OasPolygon {
                    layer: path.layer,
                    datatype: path.datatype,
                    points: vec![p0],
                };
            }

            // Perpendicular unit vector
            let nx = -dy / len;
            let ny = dx / len;

            // Calculate corners with extensions
            let ext_start = if path.start_extension > 0.0 {
                path.start_extension
            } else {
                0.0
            };
            let ext_end = if path.end_extension > 0.0 {
                path.end_extension
            } else {
                0.0
            };

            // Direction unit vector
            let dir_x = dx / len;
            let dir_y = dy / len;

            let points = vec![
                glam::DVec2::new(
                    p0.x - nx * half_width - dir_x * ext_start,
                    p0.y - ny * half_width - dir_y * ext_start,
                ),
                glam::DVec2::new(
                    p1.x - nx * half_width + dir_x * ext_end,
                    p1.y - ny * half_width + dir_y * ext_end,
                ),
                glam::DVec2::new(
                    p1.x + nx * half_width + dir_x * ext_end,
                    p1.y + ny * half_width + dir_y * ext_end,
                ),
                glam::DVec2::new(
                    p0.x + nx * half_width - dir_x * ext_start,
                    p0.y + ny * half_width - dir_y * ext_start,
                ),
            ];

            OasPolygon {
                layer: path.layer,
                datatype: path.datatype,
                points,
            }
        } else {
            // Multi-segment path: compute mitered corners
            let mut left_points = Vec::new();
            let mut right_points = Vec::new();

            for i in 0..path.points.len() - 1 {
                let p0 = path.points[i];
                let p1 = path.points[i + 1];

                let dx = p1.x - p0.x;
                let dy = p1.y - p0.y;
                let len = (dx * dx + dy * dy).sqrt();

                if len < 1e-10 {
                    continue;
                }

                let nx = -dy / len;
                let ny = dx / len;

                if i == 0 {
                    left_points.push(glam::DVec2::new(p0.x + nx * half_width, p0.y + ny * half_width));
                    right_points.push(glam::DVec2::new(p0.x - nx * half_width, p0.y - ny * half_width));
                }

                left_points.push(glam::DVec2::new(p1.x + nx * half_width, p1.y + ny * half_width));
                right_points.push(glam::DVec2::new(p1.x - nx * half_width, p1.y - ny * half_width));
            }

            // Combine left and reversed right points
            let mut points = left_points;
            right_points.reverse();
            points.extend(right_points);

            OasPolygon {
                layer: path.layer,
                datatype: path.datatype,
                points,
            }
        }
    }

    /// Convert an OASIS cell to a BRep compound with true 3D extrusion.
    ///
    /// Each polygon/path is first built as a flat face at `z_offset`,
    /// then extruded along +Z by the layer thickness. If thickness 鈮?0
    /// the face is kept as a 2D sheet.
    pub fn cell_to_brep(&self, cell: &OasCell) -> Result<BRep, OasError> {
        let mut result = BRep::default();

        // Process polygons
        for polygon in &cell.polygons {
            let ls = self.config.get(polygon.layer);
            let flat = polygon_to_brep_face(&polygon.points, ls.z_offset)?;

            if ls.thickness > 0.0 {
                match rcad_modeling::builder::ops::extrude(
                    &flat, 0, glam::DVec3::Z, ls.thickness,
                ) {
                    Ok(extruded) => merge_into(&mut result, &extruded),
                    Err(_) => merge_into(&mut result, &flat),
                }
            } else {
                merge_into(&mut result, &flat);
            }
        }

        // Process paths (convert to polygons first)
        for path in &cell.paths {
            let polygon = Self::path_to_polygon(path);
            let ls = self.config.get(path.layer);
            let flat = polygon_to_brep_face(&polygon.points, ls.z_offset)?;

            if ls.thickness > 0.0 {
                match rcad_modeling::builder::ops::extrude(
                    &flat, 0, glam::DVec3::Z, ls.thickness,
                ) {
                    Ok(extruded) => merge_into(&mut result, &extruded),
                    Err(_) => merge_into(&mut result, &flat),
                }
            } else {
                merge_into(&mut result, &flat);
            }
        }

        Ok(result)
    }
}

impl Default for OasConverter {
    fn default() -> Self {
        Self::new(LayerConfig::default())
    }
}

impl OasLibrary {
    /// Convert to BRep with layer configuration.
    pub fn to_brep(&self, cell_name: &str, config: &LayerConfig) -> Result<BRep, OasError> {
        let cell = self
            .cells
            .get(cell_name)
            .ok_or_else(|| OasError::CellNotFound(cell_name.to_string()))?;

        let converter = OasConverter::new(config.clone());
        converter.cell_to_brep(cell)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LayerSettings;
    use glam::DVec2;

    fn create_test_polygon() -> OasPolygon {
        OasPolygon {
            layer: 1,
            datatype: 0,
            points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(10.0, 0.0),
                DVec2::new(10.0, 10.0),
                DVec2::new(0.0, 10.0),
            ],
        }
    }

    fn create_test_path() -> OasPath {
        OasPath {
            layer: 2,
            datatype: 0,
            width: 2.0,
            points: vec![DVec2::new(0.0, 0.0), DVec2::new(20.0, 0.0)],
            start_extension: 0.0,
            end_extension: 0.0,
        }
    }

    #[test]
    fn test_polygon_to_face() {
        let polygon = create_test_polygon();
        let converter = OasConverter::default();
        let face = converter.polygon_to_face(&polygon).unwrap();

        // Face should have 4 edges (for a quad)
        assert_eq!(face.outer_wire.edges.len(), 4);
    }

    #[test]
    fn test_polygon_to_face_too_few_points() {
        let polygon = OasPolygon {
            layer: 1,
            datatype: 0,
            points: vec![DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0)],
        };
        let converter = OasConverter::default();
        let result = converter.polygon_to_face(&polygon);

        assert!(result.is_err());
    }

    #[test]
    fn test_path_to_face() {
        let path = create_test_path();
        let converter = OasConverter::default();
        let face = converter.path_to_face(&path).unwrap();

        // Path should be converted to a polygon with 4 edges
        assert_eq!(face.outer_wire.edges.len(), 4);
    }

    #[test]
    fn test_cell_to_brep() {
        let cell = OasCell {
            name: "TEST".to_string(),
            polygons: vec![create_test_polygon()],
            paths: vec![create_test_path()],
            texts: vec![],
            placements: vec![],
        };

        let converter = OasConverter::default();
        let brep = converter.cell_to_brep(&cell).unwrap();

        // Should have 2 solids (one per polygon/path)
        assert_eq!(brep.solids.len(), 2);
        assert!(!brep.vertices.is_empty());
    }

    #[test]
    fn test_library_to_brep() {
        let mut library = OasLibrary {
            name: Some("test".to_string()),
            cells: std::collections::HashMap::new(),
        };

        library.cells.insert("TOP".to_string(), OasCell {
            name: "TOP".to_string(),
            polygons: vec![create_test_polygon()],
            ..Default::default()
        });

        let config = LayerConfig::new().with_layer(1, LayerSettings::new(5.0));

        let brep = library.to_brep("TOP", &config).unwrap();
        assert!(!brep.solids.is_empty());
        assert!(!brep.vertices.is_empty());
    }

    #[test]
    fn test_library_to_brep_cell_not_found() {
        let library = OasLibrary {
            name: Some("test".to_string()),
            cells: std::collections::HashMap::new(),
        };
        let config = LayerConfig::default();

        let result = library.to_brep("NONEXISTENT", &config);
        assert!(matches!(result, Err(OasError::CellNotFound(_))));
    }

    #[test]
    fn test_layer_config_thickness() {
        let polygon = create_test_polygon();
        let config = LayerConfig::new().with_layer(1, LayerSettings {
            thickness: 5.0,
            z_offset: 10.0,
            color: None,
            name: None,
        });

        let converter = OasConverter::new(config);
        let face = converter.polygon_to_face(&polygon).unwrap();

        // Verify the face was created (z_offset is used in vertex creation)
        assert_eq!(face.outer_wire.edges.len(), 4);
    }
}
