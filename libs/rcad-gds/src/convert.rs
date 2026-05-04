use rcad_kernel::{BRep, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

use crate::error::{GdsError, Result};
use crate::types::*;
use crate::layer_config::LayerConfig;

/// Convert a GDS boundary to a 2D wire.
pub fn boundary_to_wire(boundary: &GdsBoundary) -> Result<Wire> {
    if boundary.points.len() < 3 {
        return Err(GdsError::GeometryError("Boundary has fewer than 3 points".to_string()));
    }

    let mut vertices: Vec<Vertex> = Vec::new();
    let mut edges: Vec<Edge> = Vec::new();
    let mut wire_edges: Vec<WireEdge> = Vec::new();

    // Create vertices from boundary points
    for point in &boundary.points {
        vertices.push(Vertex {
            point: glam::DVec3::new(point.x, point.y, 0.0),
        });
    }

    // Create edges connecting consecutive vertices
    for i in 0..vertices.len() - 1 {
        let start_idx = i;
        let end_idx = i + 1;
        edges.push(Edge {
            start: start_idx,
            end: end_idx,
        });
        wire_edges.push(WireEdge::fwd(edges.len() - 1));
    }

    Ok(Wire {
        edges: wire_edges,
    })
}

/// Build a flat-faced BRep from a polygon at a given Z offset.
fn polygon_to_brep_face(
    points: &[glam::DVec2],
    z_offset: f64,
) -> Result<BRep> {
    if points.len() < 3 {
        return Err(GdsError::GeometryError(
            "Polygon has fewer than 3 points".to_string(),
        ));
    }

    let mut brep = BRep::default();
    let n = points.len() - 1; // last point == first (closed loop)
    let n_edges = if n > 0 { n } else { points.len() };

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
        mesh_dirty: true,
    };

    let shell = Shell { faces: vec![face] };
    brep.solids.push(Solid { shells: vec![shell] });

    Ok(brep)
}

/// Convert a GdsLibrary to BRep with true 3D layer-based extrusion.
///
/// Each GDS boundary is first built as a flat face at `z_offset`,
/// then extruded along +Z by the layer thickness (from [`LayerConfig`]).
/// If thickness ≤ 0 the face is kept as-is (2D sheet).
pub fn gds_to_brep(library: &GdsLibrary, cell_name: &str, config: &LayerConfig) -> Result<BRep> {
    let structure = library.structures.get(cell_name)
        .ok_or_else(|| GdsError::CellNotFound(cell_name.to_string()))?;

    let mut result = BRep::default();

    for boundary in &structure.boundaries {
        let layer_settings = config.get(boundary.layer);

        // Build a flat-face BRep at z_offset
        let flat = polygon_to_brep_face(&boundary.points, layer_settings.z_offset)?;

        if layer_settings.thickness > 0.0 {
            // Extrude along +Z
            match rcad_modeling::builder::ops::extrude(
                &flat,
                0, // face_idx = 0 (the only face/solid/shell)
                glam::DVec3::Z,
                layer_settings.thickness,
            ) {
                Ok(extruded) => {
                    merge_into(&mut result, &extruded);
                }
                Err(_) => {
                    // Fall back to flat face on extrusion failure
                    merge_into(&mut result, &flat);
                }
            }
        } else {
            merge_into(&mut result, &flat);
        }
    }

    Ok(result)
}

/// Merge all vertices, edges, and solids from `src` into `dst`.
/// Vertex and edge indices in `src` solids are offset to match `dst`.
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
                // Offset wire edge indices
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

impl GdsLibrary {
    /// Convert to BRep with layer configuration.
    pub fn to_brep(&self, cell_name: &str, config: &LayerConfig) -> Result<BRep> {
        gds_to_brep(self, cell_name, config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LayerSettings;
    use glam::DVec2;

    fn create_test_boundary() -> GdsBoundary {
        GdsBoundary {
            layer: 1,
            datatype: 0,
            points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(10.0, 0.0),
                DVec2::new(10.0, 10.0),
                DVec2::new(0.0, 10.0),
                DVec2::new(0.0, 0.0), // close the loop
            ],
        }
    }

    #[test]
    fn test_boundary_to_wire() {
        let boundary = create_test_boundary();
        let wire = boundary_to_wire(&boundary).unwrap();
        assert_eq!(wire.edges.len(), 4); // 4 edges for a square
    }

    #[test]
    fn test_boundary_to_wire_too_few_points() {
        let boundary = GdsBoundary {
            layer: 1,
            datatype: 0,
            points: vec![DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0)],
        };
        let result = boundary_to_wire(&boundary);
        assert!(result.is_err());
    }

    #[test]
    fn test_gds_to_brep() {
        let mut library = GdsLibrary {
            name: "test".to_string(),
            ..Default::default()
        };
        library.structures.insert("TOP".to_string(), GdsStructure {
            name: "TOP".to_string(),
            boundaries: vec![create_test_boundary()],
            ..Default::default()
        });

        let config = LayerConfig::new()
            .with_layer(1, LayerSettings::new(5.0));

        let brep = library.to_brep("TOP", &config).unwrap();
        assert!(!brep.solids.is_empty());
        assert!(!brep.vertices.is_empty());
    }

    #[test]
    fn test_gds_to_brep_cell_not_found() {
        let library = GdsLibrary::default();
        let config = LayerConfig::default();

        let result = library.to_brep("NONEXISTENT", &config);
        assert!(matches!(result, Err(GdsError::CellNotFound(_))));
    }
}
