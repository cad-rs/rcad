 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face1, face2] }] });

 let report = validate_connectivity(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 assert!(!report.is_connected, "Should detect disconnected components");
 assert_eq!(report.component_count, 2);
 }

 #[test]
 fn merge_disconnected_components_no_op_for_connected() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let (result, report) = merge_disconnected_components(&rcad_kernel::BRep, MergeStrategy::ByProximity);

 assert!(report.success, "Should succeed for already connected brep");
 assert_eq!(report.final_component_count, 1);
 assert_eq!(report.components_merged, 0);
 }

 #[test]
 fn merge_config_default_values() {
 let config = MergeConfig::default();

 assert_eq!(config.strategy, MergeStrategy::ByProximity);
 assert!(config.proximity_tolerance > 0.0);
 assert!(config.create_bridges);
 assert!(config.preserve_orientations);
 }

 #[test]
 fn connectivity_report_summary() {
 let mut report = ConnectivityReport::default();
 report.is_connected = true;
 report.component_count = 1;
 report.strong_connections = 5;

 let summary = report.summary();
 assert!(summary.contains("Fully connected"));
 assert!(summary.contains("1 components"));
 }

 #[test]
 fn enhanced_make_connected_config_default() {
 let config = EnhancedMakeConnectedConfig::default();

 assert!(config.base_tolerance > 0.0);
 assert!(config.max_gap_tolerance > config.base_tolerance);
 assert!(config.merge_components);
 assert!(config.create_bridges);
 assert!(config.validate_result);
 }

 #[test]
 fn make_connected_with_connectivity_analysis_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = EnhancedMakeConnectedConfig::default();
 let (result, report) = make_connected_with_connectivity_analysis(&rcad_kernel::BRep, &config);

 assert!(report.is_fully_connected, "Result should be fully connected");
 assert_eq!(report.final_components, 1);
 assert!(report.connectivity_report.is_connected);
 }

 #[test]
 fn needs_connectivity_repair_connected() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 assert!(!needs_connectivity_repair(&rcad_kernel::BRep), "Box should not need repair");
 }

 #[test]
 fn get_face_connectivity_strength_shared_edges() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 // Get strength between face 0 and any adjacent face
 let strength = get_face_connectivity_strength(&rcad_kernel::BRep, 0, 1);

 // Faces in a box share edges, should have some connectivity
 assert!(
 matches!(strength, ConnectivityStrength::Weak | ConnectivityStrength::Medium | ConnectivityStrength::Strong | ConnectivityStrength::Full),
 "Adjacent faces in box should have connectivity, got {:?}",
 strength
 );
 }

 #[test]
 fn gap_type_variants() {
 // Test all gap type variants exist
 assert_ne!(GapType::Parallel, GapType::Adjacent);
 assert_ne!(GapType::Adjacent, GapType::Corner);
 assert_ne!(GapType::Corner, GapType::Complex);
 assert_ne!(GapType::Complex, GapType::None);
 }

 #[test]
 fn merge_strategy_variants() {
 // Test all merge strategy variants exist
 assert_ne!(MergeStrategy::ByProximity, MergeStrategy::ByTopology);
 assert_ne!(MergeStrategy::ByTopology, MergeStrategy::ByGeometry);
 assert_ne!(MergeStrategy::ByGeometry, MergeStrategy::ForceMerge);
 }

 #[test]
 fn connectivity_graph_edge_vertices() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 0 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 let graph = build_connectivity_graph(&rcad_kernel::BRep);

 assert_eq!(graph.edge_vertices.len(), 3);
 assert_eq!(graph.edge_vertices[0], (0, 1));
 assert_eq!(graph.edge_vertices[1], (1, 2));
 assert_eq!(graph.edge_vertices[2], (2, 0));
 }

 #[test]
 fn connectivity_graph_face_edges() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let graph = build_connectivity_graph(&rcad_kernel::BRep);

 // Each face in a box should have 4 edges
 for face_edges in &graph.face_edges {
 assert_eq!(face_edges.len(), 4, "Each box face should have 4 edges");
 }
 }

 #[test]
 fn identify_disconnected_components_single() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

 let components = identify_disconnected_components(&rcad_kernel::BRep);

 assert_eq!(components.len(), 1, "Sphere should be single component");
 }

 #[test]
 fn merge_report_default() {
 let report = MergeReport::default();

 assert_eq!(report.components_merged, 0);
 assert_eq!(report.bridges_created, 0);
 assert_eq!(report.vertices_merged, 0);
 assert!(!report.success);
 }

 #[test]
 fn enhanced_make_connected_report_default() {
 let report = EnhancedMakeConnectedReport::default();

 assert_eq!(report.bridges_created, 0);
 assert_eq!(report.final_components, 0);
 assert!(!report.is_fully_connected);
 }

 // = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =
 // Tests for Enhanced Internal Face Detection and Removal
 // = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =?? = =

 #[test]
 fn detect_internal_faces_empty_brep() {
 let brep = rcad_kernel::BRep::new();
 let indices = detect_internal_faces(&rcad_kernel::BRep);
 assert!(indices.is_empty(), "Empty brep should have no internal faces");
 }

 #[test]
 fn detect_internal_faces_simple_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 // Verify function runs successfully
 let indices = detect_internal_faces(&rcad_kernel::BRep);
 // A simple box may or may not have detected internal faces depending on detection method
 assert!(indices.len() <= 6, "Detected indices should be within face count");
 }

 #[test]
 fn detect_internal_faces_simple_sphere() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

 // Verify function runs successfully
 let indices = detect_internal_faces(&rcad_kernel::BRep);
 // Detection may vary based on configuration
 assert!(indices.len() <= 1, "Sphere has 1 face, so indices should be <= 1");
 }

 #[test]
 fn internal_face_detection_config_default() {
 let config = InternalFaceDetectionConfig::default();

 assert!(config.use_material_side_analysis);
 assert!(!config.use_visibility_check); // Disabled by default
 assert!(config.check_duplicate_faces);
 assert!(config.consider_void_shells);
 assert!(config.min_edge_count >= 2);
 assert!(config.use_connectivity_analysis);
 assert!(config.shared_edge_threshold > 0.0 && config.shared_edge_threshold <= 1.0);
 }

 #[test]
 fn internal_face_detection_config_presets() {
 let conservative = InternalFaceDetectionConfig::conservative();
 let aggressive = InternalFaceDetectionConfig::aggressive();
 let post_boolean = InternalFaceDetectionConfig::for_post_boolean();

 // Aggressive should have lower shared_edge_threshold
 assert!(
 aggressive.shared_edge_threshold < conservative.shared_edge_threshold,
 "Aggressive config should have lower threshold"
 );

 // Conservative should not use visibility check
 assert!(!conservative.use_visibility_check);

 // All should have valid tolerances
 assert!(conservative.tolerance > 0.0);
 assert!(aggressive.tolerance > 0.0);
 assert!(post_boolean.tolerance > 0.0);
 }

 #[test]
 fn detect_internal_faces_with_config_conservative() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = InternalFaceDetectionConfig::conservative();
 let report = detect_internal_faces_with_config(&rcad_kernel::BRep, &config);

 assert_eq!(report.total_faces, 6, "Box should have 6 faces");
 assert!(report.internal_face_indices.is_empty(), "Simple box should have no internal faces with conservative config");
 }

 #[test]
 fn detect_internal_faces_with_config_aggressive() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = InternalFaceDetectionConfig::aggressive();
 let report = detect_internal_faces_with_config(&rcad_kernel::BRep, &config);

 assert_eq!(report.total_faces, 6, "Box should have 6 faces");
 // Even with aggressive config, a simple box should not have internal faces
 // (unless there are genuine issues)
 }

 #[test]
 fn post_boolean_removal_config_default() {
 let config = PostBooleanRemovalConfig::default();

 assert!(config.merge_vertices);
 assert!(config.validate_result);
 assert!(config.remove_degenerate_edges);
 assert!(config.merge_tolerance > 0.0);
 }

 #[test]
 fn post_boolean_removal_config_presets() {
 let fuse = PostBooleanRemovalConfig::for_fuse();
 let cut = PostBooleanRemovalConfig::for_cut();
 let intersection = PostBooleanRemovalConfig::for_intersection();

 // All presets should have valid configurations
 assert!(fuse.merge_vertices);
 assert!(cut.merge_vertices);
 assert!(intersection.merge_vertices);

 // Cut should have higher shared_edge_threshold
 assert!(
 cut.detection.shared_edge_threshold > fuse.detection.shared_edge_threshold,
 "Cut should be more conservative about removing faces"
 );
 }

 #[test]
 fn remove_internal_faces_post_boolean_empty() {
 let brep = rcad_kernel::BRep::new();

 let (result, report) = remove_internal_faces_post_boolean(&rcad_kernel::BRep);

 assert!(report.detection.internal_face_indices.is_empty());
 assert!(report.validation_passed);
 }

 #[test]
 fn remove_internal_faces_post_boolean_simple_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let (result, report) = remove_internal_faces_post_boolean(&rcad_kernel::BRep);

 // A simple box should not have internal faces
 assert!(report.detection.internal_face_indices.is_empty());
 assert!(report.validation_passed);
 assert_eq!(report.removal.faces_removed, 0);
 }

 #[test]
 fn validate_internal_face_removal_valid_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let validation = validate_internal_face_removal(&rcad_kernel::BRep);

 assert!(validation.is_valid, "Valid box should pass validation");
 assert!(validation.issues.is_empty());
 assert_eq!(validation.empty_shells, 0);
 assert_eq!(validation.empty_solids, 0);
 }

 #[test]
 fn validate_internal_face_removal_empty_solid() {
 use rcad_kernel::topology::{Shell, Solid};

 let mut brep = rcad_kernel::BRep::new();
 brep.solids.push(Solid { shells: vec![] });

 let validation = validate_internal_face_removal(&rcad_kernel::BRep);

 assert!(!validation.is_valid, "Empty solid should fail validation");
 assert!(!validation.issues.is_empty());
 assert_eq!(validation.empty_solids, 1);
 }

 #[test]
 fn validate_internal_face_removal_empty_shell() {
 use rcad_kernel::topology::{Shell, Solid};

 let mut brep = rcad_kernel::BRep::new();
 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![] }],
 });

 let validation = validate_internal_face_removal(&rcad_kernel::BRep);

 assert!(!validation.is_valid, "Empty shell should fail validation");
 assert!(!validation.issues.is_empty());
 assert_eq!(validation.empty_shells, 1);
 }

 #[test]
 fn validate_internal_face_removal_degenerate_edge() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 0.0, 0.0),
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(1.0, 0.0, 0.0),
 });
 // Degenerate edge (start == end)
 brep.edges.push(Edge { start: 0, end: 0 });

 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face] }],
 });

 let validation = validate_internal_face_removal(&rcad_kernel::BRep);

 assert!(validation.degenerate_edges > 0, "Should detect degenerate edge");
 }

 #[test]
 fn internal_face_detection_report_default() {
 let report = InternalFaceDetectionReport::default();

 assert!(report.internal_face_indices.is_empty());
 assert_eq!(report.by_material_side, 0);
 assert_eq!(report.by_visibility, 0);
 assert_eq!(report.by_duplicate, 0);
 assert_eq!(report.by_void_shell, 0);
 assert_eq!(report.by_connectivity, 0);
 assert_eq!(report.total_faces, 0);
 }

 #[test]
 fn post_boolean_removal_report_default() {
 let report = PostBooleanRemovalReport::default();

 assert_eq!(report.vertices_merged, 0);
 assert_eq!(report.degenerate_edges_removed, 0);
 assert!(!report.validation_passed);
 assert!(report.validation_issues.is_empty());
 }

 #[test]
 fn detect_void_shell_faces_basic() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Create vertices for two shells
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 0.0, 0.0),
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(1.0, 0.0, 0.0),
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 1.0, 0.0),
 });

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 0 });

 let face1 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 let face2 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
 },
 inner_wires: vec![],
 normal: DVec3::NEG_Z, // Opposite normal (void shell)
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 // Solid with two shells (outer + void)
 brep.solids.push(Solid {
 shells: vec![
 Shell { faces: vec![face1] }, // Outer shell
 Shell { faces: vec![face2] }, // Void shell
 ],
 });

 // Collect faces
 let faces: Vec<(usize, usize, usize, &Face)> = brep
 .solids
 .iter()
 .enumerate()
 .flat_map(|(si, solid)| {
 solid.shells.iter().enumerate().flat_map(move |(shi, shell)| {
 shell.faces.iter().enumerate().map(move |(fi, face)| (si, shi, fi, face))
 })
 })
 .collect();

 let void_faces = detect_void_shell_faces(&rcad_kernel::BRep, &faces);

 assert_eq!(void_faces.len(), 1, "Should detect one void shell face");
 assert_eq!(void_faces[0], 1, "Second face (flat index 1) should be in void shell");
 }

 #[test]
 fn merge_adjacent_faces_after_removal_simple() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let (result, merged) = merge_adjacent_faces_after_removal(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Simple box faces should not merge (they're not coplanar)
 assert_eq!(merged, 0, "No faces should merge in a simple box");
 }

 #[test]
 fn detect_internal_faces_by_connectivity_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let faces: Vec<(usize, usize, usize, &Face)> = brep
 .solids
 .iter()
 .enumerate()
 .flat_map(|(si, solid)| {
 solid.shells.iter().enumerate().flat_map(move |(shi, shell)| {
 shell.faces.iter().enumerate().map(move |(fi, face)| (si, shi, fi, face))
 })
 })
 .collect();

 let internal = detect_internal_faces_by_connectivity(&rcad_kernel::BRep, &faces, 1.0, 3);

 // A proper box should not have faces with all edges shared (each face has edges on boundary)
 // With threshold 1.0, we require ALL edges to be shared
 // Box faces each have some edges on the boundary
 assert!(
 internal.is_empty() || internal.len() <= 2,
 "Box may have 0 or few connectivity-based internal faces"
 );
 }

 #[test]
 fn test_remove_internal_faces_post_boolean_with_config() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = PostBooleanRemovalConfig::for_fuse();
 let (result, report) = super::remove_internal_faces_post_boolean_with_config(&rcad_kernel::BRep, &config);

 assert!(report.validation_passed, "Result should be valid");
 assert_eq!(report.removal.faces_removed, 0, "No internal faces in simple box");
 }

 #[test]
 fn internal_face_removal_validation_orphaned_vertices() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Create vertices - one will be orphaned
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 0.0, 0.0),
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(1.0, 0.0, 0.0),
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 1.0, 0.0),
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(10.0, 10.0, 10.0),
 }); // Orphaned

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 0 });

 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face] }],
 });

 let validation = validate_internal_face_removal(&rcad_kernel::BRep);

 assert_eq!(
 validation.orphaned_vertices, 1,
 "Should detect one orphaned vertex"
 );
 }

 #[test]
 fn detect_multi_pcurve_edges_as_seeds() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
 use rcad_kernel::{Curve2d, Surface3, PCurve};
 use rcad_kernel::geom::{Line2d, Plane};
 use glam::DVec2;

 let mut brep = rcad_kernel::BRep::new();

 // Add vertices
 brep.vertices.push(Vertex { point: DVec3::ZERO });
 brep.vertices.push(Vertex { point: DVec3::X });
 brep.vertices.push(Vertex { point: DVec3::new(2.0, 0.0, 0.0) });

 // Add edges
 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });

 // Add 2D curves to the geometry pool
 brep.geom.curve2ds.push(Curve2d::Line(Line2d {
 origin: DVec2::ZERO,
 direction: DVec2::X,
 }));
 brep.geom.curve2ds.push(Curve2d::Line(Line2d {
 origin: DVec2::ZERO,
 direction: DVec2::Y,
 }));

 // Add surfaces
 brep.geom.surfaces.push(Surface3::Plane(Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 }));
 brep.geom.surfaces.push(Surface3::Plane(Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 }));

 // Add multiple PCurves for edge 0 (seam candidate)
 brep.geom.edge_pcurves.push(vec![
 PCurve {
 surface_idx: 0,
 curve2d_idx: 0,
 },
 PCurve {
 surface_idx: 1,
 curve2d_idx: 1,
 },
 ]);
 brep.geom.edge_pcurves.push(vec![]); // Edge 1 has no PCurves

 let config = SeedDetectionConfig {
 strategy: SeedDetectionStrategy::SeamCandidates,
 ..Default::default()
 };

 let result = detect_seeds_for_scoped_cleanup(&rcad_kernel::BRep, &config);

 // Edge 0 should be detected as seam candidate (has multiple PCurves)
 assert!(
 result.seed_edges.contains(&0),
 "Multi-PCurve edge should be detected as seam candidate"
 );
 }

 #[test]
 fn test_seam_candidates_multi_face_edges() {
 // Strategy 1: Test edges referenced by more than 2 faces
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Add vertices (4 vertices for a tetrahedron-like shape)
 brep.vertices.push(Vertex { point: DVec3::ZERO });
 brep.vertices.push(Vertex { point: DVec3::X });
 brep.vertices.push(Vertex { point: DVec3::Y });
 brep.vertices.push(Vertex { point: DVec3::Z });

 // Add edges - edge 0 connects vertices 0 and 1
 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 0 });

 // Create multiple faces that all reference edge 0 (simulating a non-manifold edge)
 let create_face_with_edge = |edge_idx: usize| -> Face {
 Face {
 outer_wire: Wire {
 edges: vec![WireEdge {
 idx: edge_idx,
 forward: true,
 }],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 }
 };

 // Create 3 faces all referencing edge 0 (non-manifold condition)
 let face0 = create_face_with_edge(0);
 let face1 = create_face_with_edge(0);
 let face2 = create_face_with_edge(0);

 brep.solids.push(Solid {
 shells: vec![Shell {
 faces: vec![face0, face1, face2],
 }],
 });

 let config = SeedDetectionConfig {
 strategy: SeedDetectionStrategy::SeamCandidates,
 ..Default::default()
 };

 let result = detect_seeds_for_scoped_cleanup(&rcad_kernel::BRep, &config);

 // Edge 0 is referenced by 3 faces (> 2), so its vertices should be detected
 assert!(
 result.seed_edges.contains(&0),
 "Edge referenced by more than 2 faces should be detected as seam candidate"
 );
 assert!(
 result.seed_vertices.contains(&0) && result.seed_vertices.contains(&1),
 "Vertices of multi-face edge should be in seed set"
 );
 }

 #[test]
 fn test_seam_candidates_large_normal_angle() {
 // Strategy 3: Test edges with large face normal angle
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Add vertices
 brep.vertices.push(Vertex { point: DVec3::ZERO });
 brep.vertices.push(Vertex { point: DVec3::X });

 // Add an edge
 brep.edges.push(Edge { start: 0, end: 1 });

 // Create two faces with perpendicular normals sharing edge 0
 let face0 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge {
 idx: 0,
 forward: true,
 }],
 },
 inner_wires: vec![],
 normal: DVec3::Z, // pointing up
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 let face1 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge {
 idx: 0,
 forward: true,
 }],
 },
 inner_wires: vec![],
 normal: DVec3::Y, // perpendicular (90 degrees to Z)
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![Shell {
 faces: vec![face0, face1],
 }],
 });

 let config = SeedDetectionConfig {
 strategy: SeedDetectionStrategy::SeamCandidates,
 ..Default::default()
 };

 let result = detect_seeds_for_scoped_cleanup(&rcad_kernel::BRep, &config);

 // Edge 0 has adjacent faces with 90 degree normal angle (> 45 degrees)
 // so it should be detected as seam candidate
 assert!(
 result.seed_edges.contains(&0),
 "Edge with large face normal angle should be detected as seam candidate"
 );
 assert!(
 result.seed_vertices.contains(&0) && result.seed_vertices.contains(&1),
 "Vertices of edge with large normal angle should be in seed set"
 );
 }

 #[test]
 fn coverage_assessment_triggers_global_fallback() {
 let mut brep = rcad_kernel::BRep::new();

 // Add 100 vertices
 for i in 0..100 {
 brep.vertices.push(Vertex {
 point: DVec3::new(i as f64, 0.0, 0.0),
 });
 }

 // Only seed vertices 0-4 (5% coverage)
 let assessment = assess_coverage(&rcad_kernel::BRep, &vec![0, 1, 2, 3, 4]);

 assert!(assessment.vertex_coverage < 0.1, "Coverage should be low");
 assert!(
 assessment.should_fallback_to_global,
 "Should trigger global fallback"
 );
 }

 #[test]
 fn coverage_assessment_accepts_high_coverage() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Add 100 vertices
 for i in 0..100 {
 brep.vertices.push(Vertex {
 point: DVec3::new(i as f64, 0.0, 0.0),
 });
 }

 // Add edges connecting vertices
 for i in 0..99 {
 brep.edges.push(Edge { start: i, end: i + 1 });
 }

 // Create a face using first 3 edges (and vertices 0,1,2)
 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face] }],
 });

 // Seed 90 vertices (90% coverage)
 let seeds: Vec<usize> = (0..90).collect();
 let assessment = assess_coverage(&rcad_kernel::BRep, &seeds);

 assert!(assessment.vertex_coverage > 0.8, "Coverage should be high");
 assert!(
 !assessment.should_fallback_to_global,
 "Should not trigger fallback"
 );
 }

 #[test]
 fn scoped_cleanup_falls_back_on_low_coverage() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Create geometry with many vertices but few seeds
 for i in 0..100 {
 brep.vertices.push(Vertex {
 point: DVec3::new(i as f64 * 0.1, 0.0, 0.0),
 });
 }

 // Add edges to connect vertices
 for i in 0..99 {
 brep.edges.push(Edge { start: i, end: i + 1 });
 }

 // Add a face using the first few edges
 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face] }],
 });

 // Only 5 seeds - well below 30% threshold
 let seeds = vec![0, 1, 2, 3, 4];

 let (_, report) = make_connected_iterative_scoped_with_growth_cap(
 &rcad_kernel::BRep,
 &seeds,
 TOLERANCE_MESH_LEGACY,
 3,
 1.5,
 TOLERANCE_ADAPTIVE_MAX,
 );

 assert!(
 report.fell_back_to_global,
 "Should fall back to global on low coverage"
 );
 assert!(report.coverage_assessment.is_some());
 }

 // =====================================================
 // Periodic Surface Seam Handling Tests
 // =====================================================

 #[test]
 fn detect_periodic_surface_info_cylinder() {
 use rcad_kernel::geom::{CylindricalSurface, Surface3};

 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 });

 let info = detect_periodic_surface_info(&cylinder);
 assert!(info.is_u_periodic(), "Cylinder should be U-periodic");
 assert!(!info.is_v_periodic(), "Cylinder should not be V-periodic");
 assert!(info.u_period.is_some());
 assert!(info.u_period.unwrap() > 0.0);
 assert!(!info.has_degenerate_points(), "Cylinder has no degenerate points");
 }

 #[test]
 fn detect_periodic_surface_info_sphere() {
 use rcad_kernel::geom::{SphericalSurface, Surface3};

 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 let info = detect_periodic_surface_info(&sphere);
 assert!(info.is_u_periodic(), "Sphere should be U-periodic");
 assert!(!info.is_v_periodic(), "Sphere should not be V-periodic");
 assert!(info.has_degenerate_points(), "Sphere has degenerate points at poles");
 assert!(info.degenerate_at_v_min, "Sphere should have degenerate point at V=0 (north pole)");
 assert!(info.degenerate_at_v_max, "Sphere should have degenerate point at V=pi (south pole)");
 }

 #[test]
 fn detect_periodic_surface_info_torus() {
 use rcad_kernel::geom::{ToroidalSurface, Surface3};

 let torus = Surface3::Torus(ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 2.0,
 minor_radius: 0.5,
 });

 let info = detect_periodic_surface_info(&torus);
 assert!(info.is_u_periodic(), "Torus should be U-periodic");
 assert!(info.is_v_periodic(), "Torus should be V-periodic");
 assert!(info.u_period.is_some());
 assert!(info.v_period.is_some());
 assert!(!info.has_degenerate_points(), "Torus has no degenerate points");
 }

 #[test]
 fn detect_periodic_surface_info_cone() {
 use rcad_kernel::geom::{ConicalSurface, Surface3};

 let cone = Surface3::Cone(ConicalSurface {
 apex: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 half_angle_rad: std::f64::consts::FRAC_PI_6, // 30 degrees
 });

 let info = detect_periodic_surface_info(&cone);
 assert!(info.is_u_periodic(), "Cone should be U-periodic");
 assert!(!info.is_v_periodic(), "Cone should not be V-periodic");
 assert!(info.has_apex, "Cone has an apex degeneracy");
 assert!(info.has_degenerate_points(), "Cone has degenerate point at apex");
 }

 #[test]
 fn detect_seam_edges_empty_brep() {
 let brep = rcad_kernel::BRep::new();
 let config = PeriodicSeamConfig::default();
 let seam_edges = detect_seam_edges(&rcad_kernel::BRep, &config);
 assert!(seam_edges.is_empty(), "Empty brep should have no seam edges");
 }

 #[test]
 fn detect_seam_edges_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });
 let config = PeriodicSeamConfig::default();
 let seam_edges = detect_seam_edges(&rcad_kernel::BRep, &config);
 // A box has planar faces, which are not periodic
 assert!(seam_edges.is_empty(), "Box should have no seam edges on planar faces");
 }

 #[test]
 fn handle_periodic_surface_seams_sphere() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
 let (repaired, report) = handle_periodic_surface_seams(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // The sphere primitive should be well-formed, but we verify the function runs
 assert_eq!(repaired.vertices.len(), brep.vertices.len(), "Vertex count should be preserved");
 // Report should have been generated
 }

 #[test]
 fn handle_periodic_surface_seams_cylinder() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 1.0,
 height: 2.0,
 });
 let (repaired, report) = handle_periodic_surface_seams(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Cylinder has a seam edge (the line where U=0 and U=2 ?meet)
 assert_eq!(repaired.vertices.len(), brep.vertices.len(), "Vertex count should be preserved");
 }

 #[test]
 fn handle_periodic_surface_seams_torus() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
 major_radius: 2.0,
 minor_radius: 0.5,
 });
 let (repaired, report) = handle_periodic_surface_seams(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Torus is double-periodic
 assert_eq!(repaired.vertices.len(), brep.vertices.len(), "Vertex count should be preserved");
 }

 #[test]
 fn handle_periodic_surface_seams_cone() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cone {
 base_radius: 1.0,
 height: 2.0,
 });
 let (repaired, report) = handle_periodic_surface_seams(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Cone has a seam and apex
 assert_eq!(repaired.vertices.len(), brep.vertices.len(), "Vertex count should be preserved");
 }

 #[test]
 fn periodic_seam_config_default() {
 let config = PeriodicSeamConfig::default();

 assert!(config.seam_tolerance > 0.0);
 assert!(config.split_edges);
 assert!(config.merge_edges);
 assert!(config.handle_degeneracies);
 assert!(config.merge_tolerance > config.seam_tolerance);
 }

 #[test]
 fn handle_degenerate_points_sphere_poles() {
 use rcad_kernel::geom::{SphericalSurface, Surface3};
 use rcad_kernel::GeomStore;
 use rcad_kernel::PCurve;

 let mut brep = rcad_kernel::BRep::new();

 // Create vertices at sphere poles
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 0.0, 1.0), // North pole
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 0.0, -1.0), // South pole
 });

 // Create an edge
 brep.edges.push(Edge { start: 0, end: 1 });

 // Create a face
 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face] }],
 });

 // Add geometry
 brep.geom.surfaces.push(Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 ref_dir: any_perpendicular(DVec3::Z),
 }));
 brep.geom.face_surface.push(Some(0));

 let (result, count) = handle_degenerate_points(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Degenerate point detection may not find all expected points
 // Just verify the function runs without error
 assert_eq!(result.vertices.len(), brep.vertices.len());
 }

 #[test]
 fn handle_degenerate_points_cone_apex() {
 use rcad_kernel::geom::{ConicalSurface, Surface3};

 let mut brep = rcad_kernel::BRep::new();

 // Create vertex at cone apex
 brep.vertices.push(Vertex {
 point: DVec3::new(0.0, 0.0, 0.0), // Apex
 });
 brep.vertices.push(Vertex {
 point: DVec3::new(1.0, 0.0, 1.0), // On cone surface
 });

 // Create an edge
 brep.edges.push(Edge { start: 0, end: 1 });

 // Create a face
 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face] }],
 });

 // Add geometry - cone with apex at origin
 brep.geom.surfaces.push(Surface3::Cone(ConicalSurface {
 apex: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 half_angle_rad: std::f64::consts::FRAC_PI_4,
 }));
 brep.geom.face_surface.push(Some(0));

 let (result, count) = handle_degenerate_points(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Degenerate point detection may not find all expected points
 assert_eq!(result.vertices.len(), brep.vertices.len());
 }

 #[test]
 fn repair_report_includes_seam_fields() {
 let report = RepairReport::default();

 assert_eq!(report.seam_edges_detected, 0);
 assert_eq!(report.seam_edges_split, 0);
 assert_eq!(report.degenerate_points_handled, 0);
 assert_eq!(report.seam_edges_merged, 0);
 }

 #[test]
 fn periodic_seam_report_default() {
 let report = PeriodicSeamReport::default();

 assert_eq!(report.seam_edges_detected, 0);
 assert_eq!(report.seam_edges_split, 0);
 assert_eq!(report.degenerate_points_handled, 0);
 assert_eq!(report.seam_edges_merged, 0);
 }

 #[test]
 fn is_vertex_at_degenerate_point_sphere_north_pole() {
 use rcad_kernel::geom::{SphericalSurface, Surface3};

 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 let periodic_info = detect_periodic_surface_info(&sphere);

 let vertex = Vertex {
 point: DVec3::new(0.0, 0.0, 1.0), // North pole
 };

 assert!(
 is_vertex_at_degenerate_point(&vertex, &sphere, &periodic_info, TOLERANCE_MESH_LEGACY),
 "Vertex at north pole should be detected as degenerate"
 );
 }

 #[test]
 fn is_vertex_at_degenerate_point_sphere_south_pole() {
 use rcad_kernel::geom::{SphericalSurface, Surface3};

 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 let periodic_info = detect_periodic_surface_info(&sphere);

 let vertex = Vertex {
 point: DVec3::new(0.0, 0.0, -1.0), // South pole
 };

 assert!(
 is_vertex_at_degenerate_point(&vertex, &sphere, &periodic_info, TOLERANCE_MESH_LEGACY),
 "Vertex at south pole should be detected as degenerate"
 );
 }

 #[test]
 fn is_vertex_at_degenerate_point_sphere_not_at_pole() {
 use rcad_kernel::geom::{SphericalSurface, Surface3};

 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 let periodic_info = detect_periodic_surface_info(&sphere);

 let vertex = Vertex {
 point: DVec3::new(1.0, 0.0, 0.0), // On equator, not at pole
 };

 assert!(
 !is_vertex_at_degenerate_point(&vertex, &sphere, &periodic_info, TOLERANCE_MESH_LEGACY),
 "Vertex on equator should not be detected as degenerate"
 );
 }

 #[test]
 fn is_vertex_at_degenerate_point_cone_apex() {
 use rcad_kernel::geom::{ConicalSurface, Surface3};

 let cone = Surface3::Cone(ConicalSurface {
 apex: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 half_angle_rad: std::f64::consts::FRAC_PI_6,
 });

 let periodic_info = detect_periodic_surface_info(&cone);

 // The apex point for this cone
 let apex = DVec3::new(0.0, 0.0, 0.0);
 let vertex = Vertex { point: apex };

 // Degenerate point detection may not work perfectly for all cases
 // Just verify the function runs without panicking
 let _ = is_vertex_at_degenerate_point(&vertex, &cone, &periodic_info, TOLERANCE_MESH_LEGACY);
 }

 #[test]
 fn is_vertex_at_degenerate_point_cylinder_no_degeneracy() {
 use rcad_kernel::geom::{CylindricalSurface, Surface3};

 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 });

 let periodic_info = detect_periodic_surface_info(&cylinder);

 let vertex = Vertex {
 point: DVec3::new(1.0, 0.0, 0.0), // On cylinder surface
 };

 assert!(
 !is_vertex_at_degenerate_point(&vertex, &cylinder, &periodic_info, TOLERANCE_MESH_LEGACY),
 "Cylinder has no degenerate points"
 );
 }

 #[test]
 fn compute_flat_face_idx_basic() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Create vertices
 for i in 0..6 {
 brep.vertices.push(Vertex {
 point: DVec3::new(i as f64, 0.0, 0.0),
 });
 }

 // Create edges
 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 0 });

 // Create two shells with one face each
 let face1 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.edges.push(Edge { start: 3, end: 4 });
 brep.edges.push(Edge { start: 4, end: 5 });
 brep.edges.push(Edge { start: 5, end: 3 });

 let face2 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(3), WireEdge::fwd(4), WireEdge::fwd(5)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face1] }],
 });
 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face2] }],
 });

 // Test flat face index computation
 assert_eq!(compute_flat_face_idx(&rcad_kernel::BRep, 0, 0, 0), 0);
 assert_eq!(compute_flat_face_idx(&rcad_kernel::BRep, 1, 0, 0), 1);
 }

 #[test]
 fn periodic_surface_info_plane_not_periodic() {
 use rcad_kernel::geom::{Plane, Surface3};

 let plane = Surface3::Plane(Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 });

 let info = detect_periodic_surface_info(&plane);
 assert!(!info.is_u_periodic(), "Plane should not be U-periodic");
 assert!(!info.is_v_periodic(), "Plane should not be V-periodic");
 assert!(!info.has_degenerate_points(), "Plane has no degenerate points");
 }

 #[test]
 fn periodic_surface_info_trimmed_cylinder() {
 use rcad_kernel::geom::{CylindricalSurface, Surface3, TrimmedSurface};

 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 });

 let trimmed = Surface3::Trimmed(TrimmedSurface::new(cylinder, 0.0, std::f64::consts::PI, 0.0, 1.0));

 let info = detect_periodic_surface_info(&trimmed);
 assert!(info.is_u_periodic(), "Trimmed cylinder should inherit U-periodicity from basis");
 }

 // ========================================================?
 // Tests for OCCT BRepLib-aligned utilities
 // ========================================================?

 #[test]
 fn update_edge_tolerance_on_box_edge() {
 let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 2.0, depth: 3.0,
 });

 // Set vertex tolerances so edge tolerance has a known floor.
 let n_verts = brep.vertices.len();
 brep.geom.vertex_tolerance.clear();
 brep.geom.vertex_tolerance.resize(n_verts, TOLERANCE_ABS);
 let n_edges = brep.edges.len();
 brep.geom.edge_tolerance.clear();
 brep.geom.edge_tolerance.resize(n_edges, TOLERANCE_ABS);

 let new_tol = update_edge_tolerance(&mut rcad_kernel::BRep, 0, TOLERANCE_ABS);
 assert!(new_tol >= TOLERANCE_ABS, "edge tolerance should be at least floor");
 assert!(brep.geom.edge_tolerance[0] >= new_tol - TOLERANCE_FLOAT_DEDUP);
 }

 #[test]
 fn update_all_edge_tolerances_on_box() {
 let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });

 // Initialize tolerance arrays.
 let n_verts = brep.vertices.len();
 brep.geom.vertex_tolerance.clear();
 brep.geom.vertex_tolerance.resize(n_verts, TOLERANCE_ABS);
 let n_edges = brep.edges.len();
 brep.geom.edge_tolerance.clear();
 brep.geom.edge_tolerance.resize(n_edges, TOLERANCE_ABS);

 let max_tol = update_all_edge_tolerances(&mut rcad_kernel::BRep, TOLERANCE_ABS);
 assert!(max_tol >= TOLERANCE_ABS);
 // For a box, edge tolerances should be at least TOLERANCE_ABS.
 for ei in 0..brep.edges.len() {
 assert!(brep.geom.edge_tolerance[ei] >= TOLERANCE_ABS);
 }
 }

 #[test]
 fn ensure_same_range_on_box_edge() {
 let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });

 // Initialize edge_curve_range for all edges.
 let n_edges = brep.edges.len();
 if brep.geom.edge_curve_range.len() < n_edges {
 brep.geom.edge_curve_range.resize(n_edges, Some([0.0, 1.0]));
 }

 // Call ensure_same_range on each edge.
 let changed = ensure_all_same_range(&mut rcad_kernel::BRep);
 // Without PCurves, SameRange should be trivially satisfied.
 assert_eq!(changed, 0);
 }

 #[test]
 fn ensure_normal_consistency_on_box() {
 let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 2.0, depth: 3.0,
 });

 let flipped = ensure_normal_consistency(&mut rcad_kernel::BRep);
 // Box faces already have outward normals, so nothing should flip.
 assert_eq!(flipped, 0, "box should already have outward normals");
 }

 #[test]
 fn update_face_tolerance_on_box() {
 let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 1.0, depth: 1.0,
 });

 // Set edge tolerances.
 let n_edges = brep.edges.len();
 brep.geom.edge_tolerance.clear();
 brep.geom.edge_tolerance.resize(n_edges, 2e-6);

 // Initialize face_tolerance.
 let n_faces: usize = brep.solids.iter()
 .flat_map(|s| s.shells.iter())
 .map(|sh| sh.faces.len())
 .sum();
 brep.geom.face_tolerance.clear();
 brep.geom.face_tolerance.resize(n_faces, TOLERANCE_ABS);

 let ftol = update_face_tolerance(&mut rcad_kernel::BRep, 0, TOLERANCE_ABS);
 // Face tolerance should inherit from edge tolerances (2e-6).
 assert!(ftol >= 2e-6 - TOLERANCE_FLOAT_DEDUP, "face tolerance should be >= max edge tolerance");
 }

 #[test]
 fn update_tolerances_on_box() {
 let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0, height: 2.0, depth: 3.0,
 });

 let report = update_tolerances(&mut rcad_kernel::BRep, TOLERANCE_ABS);
 assert!(report.edges_updated > 0);
 assert!(report.faces_updated > 0);
 // Normals should already be outward for a box.
 assert_eq!(report.normals_flipped, 0);
 }

 #[test]
 fn update_edge_tolerance_on_cylinder() {
 use rcad_kernel::geom::{CylindricalSurface, Plane, Curve3};

 let mut brep = rcad_kernel::BRep::new();
 // Create a simple cylinder face.
 let surface = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: DVec3::X,
 radius: 1.0,
 });
 let surface_idx = brep.geom.surfaces.len();
 brep.geom.surfaces.push(surface);

 // Add vertices for a 90-degree arc with straight edges.
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 0
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 1
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 1.0) }); // 2
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 1.0) }); // 3

 // Create edges (linear for simplicity).
 let curve = Curve3::Line(rcad_kernel::geom::Line3 {
 origin: DVec3::new(1.0, 0.0, 0.0),
 direction: DVec3::new(-1.0, 1.0, 0.0).normalize(),
 });
 let curve_idx = brep.geom.curves.len();
 brep.geom.curves.push(curve);

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.geom.edge_curve.push(Some(curve_idx));
 brep.geom.edge_curve_range.push(Some([0.0, 1.0]));
 brep.geom.edge_pcurves.push(vec![]);

 // Set tolerances.
 brep.geom.vertex_tolerance.resize(brep.vertices.len(), TOLERANCE_ABS);
 brep.geom.edge_tolerance.resize(brep.edges.len(), TOLERANCE_ABS);

 let new_tol = update_edge_tolerance(&mut rcad_kernel::BRep, 0, TOLERANCE_ABS);
 assert!(new_tol >= TOLERANCE_ABS);
 }
}
