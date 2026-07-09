use rcad_kernel::{topods, BRep, Edge, Vertex, Wire, WireEdge};

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
 let mut vrefs = Vec::with_capacity(n_edges);
 for pt in points.iter().take(n_edges) {
  let sr = brep.add_tvertex(glam::DVec3::new(pt.x, pt.y, z_offset));
  vrefs.push(sr);
 }

 // Edges (line segments)
 let mut erefs = Vec::with_capacity(n_edges);
 for i in 0..n_edges {
  let start = vrefs[i];
  let end = vrefs[(i + 1) % n_edges];
  let sr = brep.add_tedge(None, start, end, [0.0, 1.0]);
  erefs.push(sr);
 }

 // Wire
 let wire_sr = brep.add_twire(erefs);

 // Face (planar, no surface geometry attached)
 let face_sr = brep.add_tface(None, wire_sr, vec![], None, None, vec![], false);

 // Shell
 let shell_sr = brep.add_tshell(vec![face_sr]);

 // Solid
 brep.add_tsolid(vec![shell_sr]);

 Ok(brep)
}

/// Convert a GdsLibrary to BRep with true 3D layer-based extrusion.
///
/// Each GDS boundary is first built as a flat face at `z_offset`,
/// then extruded along +Z by the layer thickness (from [`LayerConfig`]).
/// If thickness is 0 the face is kept as-is (2D sheet).
pub fn gds_to_brep(library: &GdsLibrary, cell_name: &str, config: &LayerConfig) -> Result<BRep> {
 let structure = library.structures.get(cell_name)
 .ok_or_else(|| GdsError::CellNotFound(cell_name.to_string()))?;

 let mut result = BRep::default();

 for boundary in &structure.boundaries {
  let layer_settings = config.get(boundary.layer);

  // Build a flat-face BRep at z_offset
  let flat = polygon_to_brep_face(&boundary.points, layer_settings.z_offset)?;

  if layer_settings.thickness > 0.0 {
  // Find the face index for extrusion
  let face_idx = flat.tshapes.iter().position(|ts| {
   matches!(&**ts, topods::TShape::Face(_))
  }).ok_or_else(|| {
   GdsError::GeometryError("No face in flat polygon".to_string())
  })?;

  // Extrude along +Z
  match rcad_modeling::builder::ops::extrude(
   &flat,
   face_idx,
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

/// Merge all TShapes from `src` into `dst`.
/// ShapeRef ptr_ids stay valid because Arc references are shared.
fn merge_into(dst: &mut BRep, src: &BRep) {
 for ts in &src.tshapes {
  dst.tshapes.push(std::sync::Arc::clone(ts));
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
  assert!(brep.has_solids());
  assert!(brep.vertex_count() > 0);
 }

 #[test]
 fn test_gds_to_brep_cell_not_found() {
  let library = GdsLibrary::default();
  let config = LayerConfig::default();

  let result = library.to_brep("NONEXISTENT", &config);
  assert!(matches!(result, Err(GdsError::CellNotFound(_))));
 }
}
