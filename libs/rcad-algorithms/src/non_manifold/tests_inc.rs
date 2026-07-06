#[cfg(test)]
mod tests {
 use crate::non_manifold::*;
 use rcad_kernel::{PrimitiveSolid, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
 use glam::DVec3;

 fn unit_box() -> BRep {
 BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 })
 }

 /// Build a minimal non-manifold BRep where edge 0 is shared by 3 faces.
 fn non_manifold_tripod() -> BRep {
 let vertices = vec![
 Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0
 Vertex { point: DVec3::new(1.0, 0.0, 0.0) }, // 1
 Vertex { point: DVec3::new(0.0, 1.0, 0.0) }, // 2
 Vertex { point: DVec3::new(0.0, 0.0, 1.0) }, // 3
 Vertex { point: DVec3::new(0.0, -1.0, 0.0) }, // 4
 ];

 let edges = vec![
 Edge { start: 0, end: 1 }, // shared by 3 faces
 Edge { start: 1, end: 2 },
 Edge { start: 2, end: 0 },
 Edge { start: 1, end: 3 },
 Edge { start: 3, end: 0 },
 Edge { start: 1, end: 4 },
 Edge { start: 4, end: 0 },
 ];

 let f0 = Face {
 outer_wire: Wire {
 edges: vec![
 WireEdge::new(0, true),
 WireEdge::new(1, true),
 WireEdge::new(2, true),
 ],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 let f1 = Face {
 outer_wire: Wire {
 edges: vec![
 WireEdge::new(0, true),
 WireEdge::new(3, true),
 WireEdge::new(4, true),
 ],
 },
 inner_wires: vec![],
 normal: DVec3::Y,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 let f2 = Face {
 outer_wire: Wire {
 edges: vec![
 WireEdge::new(0, true),
 WireEdge::new(5, true),
 WireEdge::new(6, true),
 ],
 },
 inner_wires: vec![],
 normal: -DVec3::Y,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 BRep {
 vertices,
 edges,
 solids: vec![Solid {
 shells: vec![Shell {
 faces: vec![f0, f1, f2],
 }],
 }],
 geom: Default::default(),
 compound: None,
 compsolid: None,
 }
 }

 /// Build a bow-tie vertex configuration (two edge fans meeting at a vertex).
 fn bow_tie_vertex() -> BRep {
 // Two separate triangles sharing only vertex 0
 let vertices = vec![
 Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0 - bow-tie vertex
 Vertex { point: DVec3::new(1.0, 0.0, 0.0) }, // 1
 Vertex { point: DVec3::new(0.0, 1.0, 0.0) }, // 2
 Vertex { point: DVec3::new(-1.0, 0.0, 0.0) }, // 3
 Vertex { point: DVec3::new(0.0, -1.0, 0.0) }, // 4
 ];

 let edges = vec![
 Edge { start: 0, end: 1 }, // triangle 1
 Edge { start: 1, end: 2 },
 Edge { start: 2, end: 0 },
 Edge { start: 0, end: 3 }, // triangle 2
 Edge { start: 3, end: 4 },
 Edge { start: 4, end: 0 },
 ];

 let f0 = Face {
 outer_wire: Wire {
 edges: vec![
 WireEdge::new(0, true),
 WireEdge::new(1, true),
 WireEdge::new(2, true),
 ],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 let f1 = Face {
 outer_wire: Wire {
 edges: vec![
 WireEdge::new(3, true),
 WireEdge::new(4, true),
 WireEdge::new(5, true),
 ],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 BRep {
 vertices,
 edges,
 solids: vec![Solid {
 shells: vec![Shell {
 faces: vec![f0, f1],
 }],
 }],
 geom: Default::default(),
 compound: None,
 compsolid: None,
 }
 }

 #[test]
 fn test_is_manifold_for_box() {
 let brep = unit_box();
 assert!(is_manifold(&brep));
 }

 #[test]
 fn test_is_manifold_for_tripod() {
 let brep = non_manifold_tripod();
 assert!(!is_manifold(&brep));
 }

 #[test]
 fn test_non_manifold_edges_for_box() {
 let brep = unit_box();
 let nm_edges = non_manifold_edges(&brep);
 assert!(nm_edges.is_empty());
 }

 #[test]
 fn test_non_manifold_edges_for_tripod() {
 let brep = non_manifold_tripod();
 let nm_edges = non_manifold_edges(&brep);
 // Edge 0 is multi-face, edges 1-6 are boundary
 assert_eq!(nm_edges.len(), 7); // 1 multi-face + 6 boundary
 assert!(nm_edges.contains(&0));
 }

 #[test]
 fn test_multi_face_edges_for_tripod() {
 let brep = non_manifold_tripod();
 let multi = multi_face_edges(&brep);
 assert_eq!(multi, vec![0]);
 }

 #[test]
 fn test_non_manifold_vertices_for_tripod() {
 let brep = non_manifold_tripod();
 let verts = non_manifold_vertices(&brep);
 assert_eq!(verts, vec![0, 1]); // endpoints of edge 0
 }

 #[test]
 fn test_analyze_non_manifold_for_box() {
 let brep = unit_box();
 let report = analyze_non_manifold(&brep);
 assert!(report.is_manifold);
 assert!(report.is_closed);
 assert_eq!(report.multi_face_edge_count, 0);
 assert_eq!(report.boundary_edge_count, 0);
 assert!(report.is_clean());
 }

 #[test]
 fn test_analyze_non_manifold_for_tripod() {
 let brep = non_manifold_tripod();
 let report = analyze_non_manifold(&brep);
 assert!(!report.is_manifold);
 assert!(!report.is_closed);
 assert_eq!(report.multi_face_edge_count, 1);
 assert_eq!(report.boundary_edge_count, 6);
 assert_eq!(report.non_manifold_vertex_count, 2);
 assert!(!report.is_clean());
 }

 #[test]
 fn test_split_non_manifold_edges_for_box() {
 let brep = unit_box();
 let (result, report) = split_non_manifold_edges(&brep);
 assert!(is_manifold(&result));
 assert_eq!(report.edges_split, 0);
 }

 #[test]
 fn test_split_non_manifold_edges_for_tripod() {
 let brep = non_manifold_tripod();
 let (result, report) = split_non_manifold_edges(&brep);

 // After splitting, the multi-face edge should be resolved
 assert!(report.edges_split > 0);
 assert!(report.new_edges_created > 0);

 // Verify the mapping
 assert!(report.edge_mapping.contains_key(&0));
 }

 #[test]
 fn test_make_manifold_for_box() {
 let brep = unit_box();
 let (result, report) = make_manifold(&brep).expect("should succeed");
 assert!(report.was_already_manifold);
 assert!(report.is_manifold);
 }

 #[test]
 fn test_make_manifold_for_tripod() {
 let brep = non_manifold_tripod();
 let (result, report) = make_manifold(&brep).expect("should succeed");
 assert!(!report.was_already_manifold);
 // After splitting, boundary edges remain, so not fully manifold in the closed sense
 // but the multi-face edge should be resolved
 }

 #[test]
 fn test_non_manifold_traversal() {
 let brep = non_manifold_tripod();
 let graph = BRepGraph::from_brep(&brep);

 // Test non_manifold_adjacent_faces
 let adj = graph.non_manifold_adjacent_faces(0);
 // Face 0 shares edge 0 with faces 1 and 2
 assert!(adj.contains(&1));
 assert!(adj.contains(&2));

 // Test manifold_regions
 let regions = graph.manifold_regions();
 // With a multi-face edge, faces should still be connected via that edge
 // but our manifold_regions skips non-manifold edges
 assert!(!regions.is_empty());

 // Test non_manifold_edge_info
 let info = graph.non_manifold_edge_info();
 assert_eq!(info.len(), 1);
 assert_eq!(info[0].0, 0); // edge 0
 assert_eq!(info[0].1.len(), 3); // 3 adjacent faces
 }

 #[test]
 fn test_boundary_edges_for_tripod() {
 let brep = non_manifold_tripod();
 let bounds = boundary_edges(&brep);
 // Edges 1-6 are boundary edges (1 face each)
 assert_eq!(bounds.len(), 6);
 }

 #[test]
 fn test_orphan_edges() {
 let brep = unit_box();
 let orphans = orphan_edges(&brep);
 assert!(orphans.is_empty());
 }

 //  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
 // Tests for new functionality
 //  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

 #[test]
 fn test_detect_non_manifold_topology_box() {
 let brep = unit_box();
 let report = detect_non_manifold_topology(&brep);

 assert!(report.is_manifold);
 assert!(report.is_closed);
 assert!(report.edge_details.is_empty());
 assert!(report.vertex_details.is_empty());
 assert_eq!(report.counts.manifold_edges, brep.edges.len());
 }

 #[test]
 fn test_detect_non_manifold_topology_tripod() {
 let brep = non_manifold_tripod();
 let report = detect_non_manifold_topology(&brep);

 assert!(!report.is_manifold);
 assert!(!report.is_closed);
 assert_eq!(report.counts.multi_face_edges, 1);
 assert_eq!(report.counts.boundary_edges, 6);

 // Check edge details
 assert_eq!(report.edge_details.len(), 7);
 let multi_face_edge = report.edge_details.iter()
 .find(|e| e.edge_type == NonManifoldEdgeType::MultiFace);
 assert!(multi_face_edge.is_some());
 let edge = multi_face_edge.unwrap();
 assert_eq!(edge.adjacent_face_count, 3);

 // Check vertex details
 assert!(report.vertex_details.iter().any(|v| v.vertex_type == NonManifoldVertexType::MultiFaceJunction));
 }

 #[test]
 fn test_detect_non_manifold_topology_bow_tie() {
 let brep = bow_tie_vertex();
 let report = detect_non_manifold_topology(&brep);

 // Bow-tie has a vertex where two edge fans meet
 assert!(report.vertex_details.iter().any(|v| {
 matches!(v.vertex_type, NonManifoldVertexType::BowTie) || v.fan_count > 1
 }));
 }

 #[test]
 fn test_is_non_manifold() {
 let box_brep = unit_box();
 let tripod = non_manifold_tripod();

 assert!(!is_non_manifold(&box_brep));
 assert!(is_non_manifold(&tripod));
 }

 #[test]
 fn test_count_non_manifold_entities() {
 let brep = non_manifold_tripod();
 let counts = count_non_manifold_entities(&brep);

 assert_eq!(counts.multi_face_edges, 1);
 assert_eq!(counts.boundary_edges, 6);
 assert!(!counts.is_manifold());
 }

 #[test]
 fn test_convert_to_manifold_box() {
 let brep = unit_box();
 let (result, report) = convert_to_manifold(&brep);

 assert!(report.was_already_manifold);
 assert!(report.is_manifold);
 assert_eq!(report.edges_split, 0);
 assert_eq!(report.vertices_duplicated, 0);
 }

 #[test]
 fn test_convert_to_manifold_tripod() {
 let brep = non_manifold_tripod();
 let (result, report) = convert_to_manifold(&brep);

 assert!(!report.was_already_manifold);
 assert!(report.edges_split > 0);
 assert!(report.new_edges_created > 0);

 // Check that we have split details
 assert!(!report.edge_split_details.is_empty());
 }

 #[test]
 fn test_convert_to_manifold_with_options_conservative() {
 let brep = non_manifold_tripod();
 let options = ManifoldConversionOptions::conservative();
 let (result, report) = convert_to_manifold_with_options(&brep, options);

 // Conservative options should still split edges
 assert!(report.edges_split > 0);
 // But not duplicate vertices
 assert_eq!(report.vertices_duplicated, 0);
 }

 #[test]
 fn test_convert_to_manifold_with_options_aggressive() {
 let brep = non_manifold_tripod();
 let options = ManifoldConversionOptions::aggressive();
 let (result, report) = convert_to_manifold_with_options(&brep, options);

 assert!(report.edges_split > 0);
 }

 #[test]
 fn test_manifold_conversion_report_summary() {
 let brep = unit_box();
 let (_, report) = convert_to_manifold(&brep);
 assert_eq!(report.summary(), "Already manifold");

 let tripod = non_manifold_tripod();
 let (_, report) = convert_to_manifold(&tripod);
 assert!(report.summary().contains("edges split"));
 }

 #[test]
 fn test_non_manifold_sewing_options() {
 let strict = NonManifoldSewingOptions::strict(TOLERANCE_MESH_LEGACY);
 assert_eq!(strict.non_manifold_mode, NonManifoldSewingMode::StrictManifold);

 let allow = NonManifoldSewingOptions::allow_non_manifold(TOLERANCE_MESH_LEGACY);
 assert_eq!(allow.non_manifold_mode, NonManifoldSewingMode::AllowNonManifold);

 let create = NonManifoldSewingOptions::create_non_manifold(TOLERANCE_MESH_LEGACY);
 assert_eq!(create.non_manifold_mode, NonManifoldSewingMode::CreateNonManifold);
 }

 #[test]
 fn test_sew_non_manifold_aware_strict() {
 let brep = unit_box();
 let options = NonManifoldSewingOptions::strict(TOLERANCE_MESH_LEGACY);
 let (result, report) = sew_non_manifold_aware(&brep, &options);

 // A closed box has no free edges to sew
 assert!(report.is_successful());
 assert!(report.is_manifold);
 }

 #[test]
 fn test_sew_non_manifold_aware_allow_non_manifold() {
 let brep = unit_box();
 let options = NonManifoldSewingOptions::allow_non_manifold(TOLERANCE_MESH_LEGACY);
 let (result, report) = sew_non_manifold_aware(&brep, &options);

 assert!(report.is_successful());
 assert!(report.is_manifold);
 }

 #[test]
 fn test_make_connected_non_manifold_aware() {
 let brep = unit_box();
 let options = NonManifoldMakeConnectedOptions::default();
 let (result, report) = make_connected_non_manifold_aware(&brep, &options);

 assert!(report.is_manifold);
 assert_eq!(report.vertices_merged, 0); // No duplicates to merge
 }

 #[test]
 fn test_non_manifold_edge_type_classification() {
 let brep = non_manifold_tripod();
 let report = detect_non_manifold_topology(&brep);

 // Check that edges are classified correctly
 let multi_face = report.edge_details.iter()
 .filter(|e| e.edge_type == NonManifoldEdgeType::MultiFace)
 .count();
 assert_eq!(multi_face, 1);

 let boundary = report.edge_details.iter()
 .filter(|e| e.edge_type == NonManifoldEdgeType::Boundary)
 .count();
 assert_eq!(boundary, 6);
 }

 #[test]
 fn test_non_manifold_vertex_type_classification() {
 let brep = non_manifold_tripod();
 let report = detect_non_manifold_topology(&brep);

 // Vertices 0 and 1 are on the multi-face edge
 let multi_face_junctions = report.vertex_details.iter()
 .filter(|v| v.vertex_type == NonManifoldVertexType::MultiFaceJunction)
 .count();
 assert!(multi_face_junctions >= 2);
 }

 #[test]
 fn test_non_manifold_counts_methods() {
 let mut counts = NonManifoldCounts::default();
 counts.multi_face_edges = 2;
 counts.boundary_edges = 3;
 counts.orphan_edges = 1;

 assert_eq!(counts.non_manifold_edge_count(), 6);
 assert!(!counts.is_manifold());
 }

 #[test]
 fn test_detailed_non_manifold_report_clean() {
 let brep = unit_box();
 let report = detect_non_manifold_topology(&brep);

 assert!(report.is_clean());
 // Verify basic report structure
 assert!(report.counts.total_edges > 0);
 }

 #[test]
 fn test_detailed_non_manifold_report_non_clean() {
 let brep = non_manifold_tripod();
 let report = detect_non_manifold_topology(&brep);

 assert!(!report.is_clean());
 assert!(report.counts.multi_face_edges > 0);
 }

 #[test]
 fn test_manifold_region_count() {
 let brep = unit_box();
 let report = detect_non_manifold_topology(&brep);

 // A single box should have 1 manifold region
 assert_eq!(report.manifold_region_count, 1);
 }

 #[test]
 fn test_edge_fan_computation() {
 let brep = bow_tie_vertex();
 let graph = BRepGraph::from_brep(&brep);

 // Vertex 0 is the bow-tie center
 let incident_edges: Vec<usize> = graph.vertex_adjacent_edges(0).to_vec();
 let fans = compute_edge_fans(0, &incident_edges, &graph, &brep);

 // Should have 2 separate fans (one for each triangle)
 assert_eq!(fans.len(), 2);
 }
}
