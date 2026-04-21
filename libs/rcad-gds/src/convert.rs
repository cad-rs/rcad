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

/// Convert a GdsLibrary to BRep with layer-based extrusion.
pub fn gds_to_brep(library: &GdsLibrary, cell_name: &str, config: &LayerConfig) -> Result<BRep> {
    let structure = library.structures.get(cell_name)
        .ok_or_else(|| GdsError::CellNotFound(cell_name.to_string()))?;

    let mut brep = BRep::default();

    // Process boundaries
    for boundary in &structure.boundaries {
        let layer_settings = config.get(boundary.layer);

        // Create a face from the boundary
        let _wire = boundary_to_wire(boundary)?;

        // Create vertices for this boundary
        let mut face_vertices: Vec<Vertex> = boundary.points.iter()
            .map(|p| Vertex {
                point: glam::DVec3::new(p.x, p.y, layer_settings.z_offset),
            })
            .collect();

        // Create edges
        let mut face_edges: Vec<Edge> = Vec::new();
        let mut wire_edges: Vec<WireEdge> = Vec::new();

        for i in 0..boundary.points.len() - 1 {
            face_edges.push(Edge {
                start: i,
                end: i + 1,
            });
            wire_edges.push(WireEdge::fwd(face_edges.len() - 1));
        }

        let _vertex_start = brep.vertices.len();
        let edge_start = brep.edges.len();

        brep.vertices.append(&mut face_vertices);
        brep.edges.append(&mut face_edges);

        // Adjust wire edge indices
        let wire = Wire {
            edges: wire_edges.iter()
                .map(|we| WireEdge::fwd(we.idx + edge_start))
                .collect(),
        };

        let face = Face {
            outer_wire: wire,
            inner_wires: Vec::new(),
            normal: glam::DVec3::Z,
            triangles: Vec::new(),
            mesh_dirty: true,
        };

        // Create a solid from this face (simplified)
        let shell = Shell {
            faces: vec![face],
        };
        let solid = Solid {
            shells: vec![shell],
        };
        brep.solids.push(solid);
    }

    Ok(brep)
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
