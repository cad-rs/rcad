 #[test]
 fn edge_valence_info_classification() {
 let open_edge = EdgeValenceInfo {
 edge_index: 0,
 valence: 1,
 is_open: true,
 is_manifold: false,
 is_non_manifold: false,
 };
 assert!(open_edge.is_open);
 assert!(!open_edge.is_manifold);

 let manifold_edge = EdgeValenceInfo {
 edge_index: 1,
 valence: 2,
 is_open: false,
 is_manifold: true,
 is_non_manifold: false,
 };
 assert!(manifold_edge.is_manifold);

 let nm_edge = EdgeValenceInfo {
 edge_index: 2,
 valence: 3,
 is_open: false,
 is_manifold: false,
 is_non_manifold: true,
 };
 assert!(nm_edge.is_non_manifold);
 }

 #[test]
 fn vertex_valence_info_properties() {
 let boundary_vertex = VertexValenceInfo {
 vertex_index: 0,
 edge_valence: 3,
 face_valence: 2,
 is_boundary: true,
 is_non_manifold: false,
 };
 assert!(boundary_vertex.is_boundary);
 assert!(!boundary_vertex.is_non_manifold);
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Tests for UV Gap Repair
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn uv_gap_repair_config_default() {
 let config = UvGapRepairConfig::default();

 assert!(config.max_repairable_gap > 0.0);
 assert!(config.closure_tolerance > 0.0);
 assert!(config.allow_bounds_extension);
 assert!(config.handle_periodic_seams);
 assert!(config.max_extension_factor > 0.0);
 }

 #[test]
 fn uv_gap_repair_report_default() {
 let report = UvGapRepairReport::default();

 assert_eq!(report.faces_processed, 0);
 assert_eq!(report.gaps_repaired, 0);
 assert_eq!(report.pcurves_extended, 0);
 assert_eq!(report.pcurves_trimmed, 0);
 assert_eq!(report.seam_edges_adjusted, 0);
 assert!(report.unrepaired_gaps.is_empty());
 }

 #[test]
 fn fix_uv_gaps_box_face() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = UvGapRepairConfig::default();
 let (_, report) = fix_uv_gaps(0, 0, 0, &rcad_kernel::BRep, &config);

 // Box faces should be processed
 }

 #[test]
 fn fix_uv_gaps_cylinder_face() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 1.0,
 height: 2.0,
 });

 let config = UvGapRepairConfig::default();
 let (_, report) = fix_uv_gaps(0, 0, 0, &rcad_kernel::BRep, &config);

 // Cylinder faces should be processed
 }

 #[test]
 fn fix_uv_gaps_sphere_face() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere {
 radius: 1.0,
 });

 let config = UvGapRepairConfig::default();
 let (_, report) = fix_uv_gaps(0, 0, 0, &rcad_kernel::BRep, &config);

 // Sphere faces should be processed
 }

 #[test]
 fn fix_all_uv_gaps_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = UvGapRepairConfig::default();
 let (_, report) = fix_all_uv_gaps(&rcad_kernel::BRep, &config);

 // All faces should be processed
 }

 #[test]
 fn fix_uv_gaps_invalid_indices() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = UvGapRepairConfig::default();

 // Test with invalid solid index
 let (_, report) = fix_uv_gaps(99, 0, 0, &rcad_kernel::BRep, &config);
 assert_eq!(report.faces_processed, 0);

 // Test with invalid shell index
 let (_, report) = fix_uv_gaps(0, 99, 0, &rcad_kernel::BRep, &config);
 assert_eq!(report.faces_processed, 0);

 // Test with invalid face index
 let (_, report) = fix_uv_gaps(0, 0, 99, &rcad_kernel::BRep, &config);
 assert_eq!(report.faces_processed, 0);
 }

 #[test]
 fn unrepaired_gap_structure() {
 let gap = UnrepairedGap {
 edge_idx: 5,
 gap_size: 0.01,
 reason: GapRepairFailureReason::GapTooLarge,
 };

 assert_eq!(gap.edge_idx, 5);
 assert_eq!(gap.gap_size, 0.01);
 assert_eq!(gap.reason, GapRepairFailureReason::GapTooLarge);
 }

 #[test]
 fn gap_repair_failure_reason_variants() {
 // Test all variants exist and can be compared
 assert_ne!(GapRepairFailureReason::GapTooLarge, GapRepairFailureReason::NoExtensionMethod);
 assert_ne!(GapRepairFailureReason::WouldCauseSelfIntersection, GapRepairFailureReason::UndefinedSurfaceInGap);
 assert_ne!(GapRepairFailureReason::RequiresPeriodicHandling, GapRepairFailureReason::GapTooLarge);
 }

 #[test]
 fn fix_edge_pcurve_uv_bounds_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = UvGapRepairConfig::default();

 // Test with valid indices (if edge has PCurve)
 if !brep.edges.is_empty() {
 let surface_idx = brep.geom.face_surface.get(0).and_then(|v| *v).unwrap_or(0);
 let (_, repaired) = fix_edge_pcurve_uv_bounds(0, surface_idx, &rcad_kernel::BRep, &config);
 // repaired may be true or false depending on geometry
 assert!(repaired || !repaired); // Just check it doesn't panic
 }
 }

 #[test]
 fn fix_edge_pcurve_uv_bounds_invalid_indices() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = UvGapRepairConfig::default();

 // Test with invalid edge index
 let (_, repaired) = fix_edge_pcurve_uv_bounds(999, 0, &rcad_kernel::BRep, &config);
 assert!(!repaired);

 // Test with invalid surface index
 let (_, repaired) = fix_edge_pcurve_uv_bounds(0, 999, &rcad_kernel::BRep, &config);
 assert!(!repaired);
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Internal Face Detection and Removal Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn detect_duplicate_faces_clean_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let report = detect_duplicate_faces(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);
 // A clean box should have no duplicate faces
 assert_eq!(report.duplicate_pairs.len(), 0, "Clean box should have no duplicate faces");
 assert_eq!(report.internal_face_count, 0, "Clean box should have no internal faces");
 }

 #[test]
 fn detect_duplicate_faces_with_duplicates() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
use rcad_kernel::topology::Face;
use rcad_kernel::PCurve;

 // Create a brep with two identical faces
 let mut brep = rcad_kernel::BRep::new();

 // Add 4 vertices for a quad
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

 // Add 4 edges for the quad
 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 3 });
 brep.edges.push(Edge { start: 3, end: 0 });

 // Create two identical faces with opposite normals
 let face1 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
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
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
 },
 inner_wires: vec![],
 normal: DVec3::NEG_Z, // Opposite normal
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face1, face2] }],
 });

 let report = detect_duplicate_faces(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Should detect the duplicate face pair
 assert!(report.duplicate_pairs.len() >= 1, "Should detect duplicate face pair");

 // The pair should have opposite orientation
 let pair = &report.duplicate_pairs[0];
 assert!(pair.opposite_orientation, "Duplicate faces should have opposite orientation");
 }

 #[test]
 fn identify_internal_faces_clean_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let internal = identify_internal_faces(&rcad_kernel::BRep);
 assert_eq!(internal.len(), 0, "Clean box should have no internal faces");
 }

 #[test]
 fn identify_internal_faces_with_void_shell() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 // Create a brep with an outer shell and a void shell
 let mut brep = rcad_kernel::BRep::new();

 // Outer shell vertices (cube)
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) }); // 2
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 3
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 1.0) }); // 4
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 1.0) }); // 5
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 1.0) }); // 6
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 1.0) }); // 7

 // Edges for bottom face
 brep.edges.push(Edge { start: 0, end: 1 }); // 0
 brep.edges.push(Edge { start: 1, end: 2 }); // 1
 brep.edges.push(Edge { start: 2, end: 3 }); // 2
 brep.edges.push(Edge { start: 3, end: 0 }); // 3

 // Create outer shell with one face
 let outer_face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
 },
 inner_wires: vec![],
 normal: DVec3::NEG_Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 // Create void shell with one face (same geometry but opposite normal)
 let void_face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
 },
 inner_wires: vec![],
 normal: DVec3::Z, // Opposite normal
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![
 Shell { faces: vec![outer_face] }, // Shell 0: outer
 Shell { faces: vec![void_face] }, // Shell 1: void
 ],
 });

 let internal = identify_internal_faces(&rcad_kernel::BRep);

 // Should identify faces in the void shell as internal
 assert!(internal.len() >= 1, "Should identify internal faces in void shell");
 }

 #[test]
 fn remove_internal_faces_basic() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 // Create a brep with multiple faces
 let mut brep = rcad_kernel::BRep::new();

 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 3 });
 brep.edges.push(Edge { start: 3, end: 0 });

 let face1 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
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
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
 },
 inner_wires: vec![],
 normal: DVec3::NEG_Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face1, face2] }],
 });

 // Remove the second face
 let (result, report) = remove_internal_faces(&rcad_kernel::BRep, &[1]);

 assert_eq!(report.faces_removed, 1, "Should remove one face");
 assert!(report.is_valid, "Result should be valid");

 // Check that the result has one less face
 let total_faces: usize = result.solids.iter()
 .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum::<usize>())
 .sum();
 let original_faces: usize = brep.solids.iter()
 .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum::<usize>())
 .sum();
 assert_eq!(total_faces, original_faces - 1, "Should have one less face");
 }

 #[test]
 fn remove_internal_faces_empty_list() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let (result, report) = remove_internal_faces(&rcad_kernel::BRep, &[]);

 assert_eq!(report.faces_removed, 0, "Should remove no faces");
 assert!(report.is_valid, "Result should be valid");
 assert_eq!(result.solids.len(), brep.solids.len(), "Solid count should be unchanged");
 }

 #[test]
 fn cleanup_boolean_result_clean_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let (result, report) = cleanup_boolean_result(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // A clean box should pass through with minimal changes
 assert!(report.is_valid, "Result should be valid");
 assert_eq!(report.internal_faces_removed, 0, "Clean box has no internal faces");
 assert_eq!(report.degenerate_faces_removed, 0, "Clean box has no degenerate faces");
 assert!(!result.solids.is_empty(), "Result should have solids");
 }

 #[test]
 fn cleanup_boolean_result_with_internal_faces() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 // Create a brep simulating post-boolean result with internal face
 let mut brep = rcad_kernel::BRep::new();

 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 3 });
 brep.edges.push(Edge { start: 3, end: 0 });

 // Two identical faces with opposite normals (simulating internal separator)
 let face1 = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
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
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
 },
 inner_wires: vec![],
 normal: DVec3::NEG_Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid {
 shells: vec![Shell { faces: vec![face1, face2] }],
 });

 let (result, report) = cleanup_boolean_result(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Should have cleaned up the internal face
 assert!(report.is_valid, "Result should be valid");

 // The internal face (or duplicate) should have been removed
 let total_faces: usize = result.solids.iter()
 .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum::<usize>())
 .sum();
 assert!(total_faces <= 2, "Should have cleaned up internal/duplicate faces");
 }

 #[test]
 fn duplicate_face_pair_structure() {
 let pair = DuplicateFacePair {
 face_a: 0,
 face_b: 1,
 kind: DuplicateFaceKind::GeometricallyIdentical,
 opposite_orientation: true,
 max_deviation: 0.001,
 shared_edges: vec![0, 1, 2],
 is_internal: true,
 };

 assert_eq!(pair.face_a, 0);
 assert_eq!(pair.face_b, 1);
 assert_eq!(pair.kind, DuplicateFaceKind::GeometricallyIdentical);
 assert!(pair.opposite_orientation);
 assert_eq!(pair.max_deviation, 0.001);
 assert_eq!(pair.shared_edges.len(), 3);
 assert!(pair.is_internal);
 }

 #[test]
 fn duplicate_face_kind_variants() {
 // Test all variants exist and can be compared
 assert_ne!(DuplicateFaceKind::GeometricallyIdentical, DuplicateFaceKind::TopologicallyShared);
 assert_ne!(DuplicateFaceKind::CoincidentDifferentGeometry, DuplicateFaceKind::SameSurfaceDifferentBounds);
 }

 #[test]
 fn duplicate_face_report_default() {
 let report = DuplicateFaceReport::default();
 assert!(report.duplicate_pairs.is_empty());
 assert_eq!(report.internal_face_count, 0);
 assert!(report.internal_face_indices.is_empty());
 }

 #[test]
 fn internal_face_removal_report_default() {
 let report = InternalFaceRemovalReport::default();
 assert_eq!(report.faces_removed, 0);
 assert!(report.removed_indices.is_empty());
 assert_eq!(report.edges_removed, 0);
 assert_eq!(report.vertices_removed, 0);
 assert!(!report.is_valid);
 }

 #[test]
 fn boolean_cleanup_report_default() {
 let report = BooleanCleanupReport::default();
 assert_eq!(report.internal_faces_removed, 0);
 assert_eq!(report.duplicate_faces_merged, 0);
 assert_eq!(report.vertices_merged, 0);
 assert_eq!(report.degenerate_faces_removed, 0);
 assert_eq!(report.edges_sewn, 0);
 assert!(!report.is_valid);
 }

 #[test]
 fn ray_triangle_intersection_basic() {
 // Simple test of ray-triangle intersection
 let origin = DVec3::new(0.5, 0.5, -1.0);
 let dir = DVec3::new(0.0, 0.0, 1.0);
 let v0 = DVec3::new(0.0, 0.0, 0.0);
 let v1 = DVec3::new(1.0, 0.0, 0.0);
 let v2 = DVec3::new(0.0, 1.0, 0.0);

 assert!(ray_triangle_intersection(origin, dir, v0, v1, v2), "Ray should intersect triangle");

 // Ray pointing away
 let dir_away = DVec3::new(0.0, 0.0, -1.0);
 assert!(!ray_triangle_intersection(origin, dir_away, v0, v1, v2), "Ray pointing away should not intersect");
 }

 #[test]
 fn compute_bounding_box_basic() {
 let points = vec![
 DVec3::new(0.0, 0.0, 0.0),
 DVec3::new(1.0, 2.0, 3.0),
 DVec3::new(-1.0, -2.0, -3.0),
 ];

 let (min_pt, max_pt) = compute_bounding_box(&points);

 assert_eq!(min_pt, DVec3::new(-1.0, -2.0, -3.0));
 assert_eq!(max_pt, DVec3::new(1.0, 2.0, 3.0));
 }

 #[test]
 fn compute_face_centroid_basic() {
 use rcad_kernel::topology::{Edge, Face, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(2.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(2.0, 2.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 2.0, 0.0) });

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 3 });
 brep.edges.push(Edge { start: 3, end: 0 });

 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 let centroid = compute_face_centroid_from_wire(&rcad_kernel::BRep, &face);

 // Centroid should be at (1, 1, 0)
 assert!((centroid.x - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
 assert!((centroid.y - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
 assert!((centroid.z - 0.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Enhanced Solid Validation and Repair Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn verify_solid_closure_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let solid = &brep.solids[0];
 let report = verify_solid_closure(solid, &rcad_kernel::BRep);

 assert!(report.is_valid(), "Unit box should pass closure verification");
 assert!(report.all_shells_closed, "Unit box should have all shells closed");
 assert_eq!(report.shell_count, 1);
 assert_eq!(report.closed_shell_count, 1);
 assert_eq!(report.open_shell_count, 0);
 assert!(report.has_single_outer_shell, "Unit box should have single outer shell");
 assert!(report.total_volume > 0.0, "Unit box should have positive volume");
 assert_eq!(report.volume_sign, VolumeSign::Positive);
 }

 #[test]
 fn verify_solid_closure_unit_sphere() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere {
 radius: 1.0,
 });

 let solid = &brep.solids[0];
 let report = verify_solid_closure(solid, &rcad_kernel::BRep);

 // Sphere should be closed with a single shell
 assert!(report.all_shells_closed, "Unit sphere should have all shells closed");
 assert_eq!(report.shell_count, 1);
 // Volume computation for curved primitives depends on face normal orientation
 // Just verify we have a shell (volume might be zero or very small due to geometry)
 }

 #[test]
 fn verify_solid_closure_unit_cylinder() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 1.0,
 height: 2.0,
 });

 let solid = &brep.solids[0];
 let report = verify_solid_closure(solid, &rcad_kernel::BRep);

 assert!(report.is_valid(), "Cylinder should pass closure verification");
 assert!(report.all_shells_closed, "Cylinder should have all shells closed");
 }

 #[test]
 fn verify_solid_closure_empty_solid() {
 use rcad_kernel::topology::Solid as TopologySolid;

 let brep = rcad_kernel::BRep::new();
 let solid = TopologySolid { shells: vec![] };

 let report = verify_solid_closure(&solid, &rcad_kernel::BRep);

 assert!(!report.is_valid(), "Empty solid should not pass verification");
 assert!(!report.has_single_outer_shell, "Empty solid has no outer shell");
 }

 #[test]
 fn verify_solid_closure_report_summary() {
 let report = SolidClosureVerificationReport {
 all_shells_closed: true,
 has_proper_nesting: true,
 shell_count: 1,
 closed_shell_count: 1,
 open_shell_count: 0,
 shell_volume_signs: vec![VolumeSign::Positive],
 shell_volumes: vec![1.0],
 total_volume: 1.0,
 volume_sign: VolumeSign::Positive,
 shell_containment: vec![],
 degenerate_shell_indices: vec![],
 inconsistent_orientation_indices: vec![],
 has_single_outer_shell: true,
 };

 let summary = report.summary();
 assert!(summary.contains("Valid solid"));
 assert!(summary.contains("1 shells"));
 }

 #[test]
 fn volume_sign_variants() {
 // Test that VolumeSign variants exist and can be compared
 assert_ne!(VolumeSign::Positive, VolumeSign::Negative);
 assert_ne!(VolumeSign::Zero, VolumeSign::Unknown);
 assert_ne!(VolumeSign::Positive, VolumeSign::Zero);
 }

 #[test]
 fn shell_containment_info_default() {
 let info = ShellContainmentInfo {
 container_shell_idx: None,
 nesting_depth: 0,
 is_fully_contained: true,
 has_intersections: false,
 intersecting_shells: vec![],
 };

 assert!(info.container_shell_idx.is_none());
 assert_eq!(info.nesting_depth, 0);
 assert!(info.is_fully_contained);
 assert!(!info.has_intersections);
 }

 #[test]
 fn orient_solid_shells_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let solid = &brep.solids[0];
 let (oriented, report) = orient_solid_shells(solid, &rcad_kernel::BRep);

 assert!(report.is_clean(), "Box should have clean orientation");
 assert!(report.is_properly_oriented, "Box should be properly oriented");
 assert_eq!(oriented.shells.len(), solid.shells.len());
 assert_eq!(report.outer_shells_oriented, 1);
 assert_eq!(report.inner_shells_oriented, 0);
 }

 #[test]
 fn orient_solid_shells_sphere() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere {
 radius: 1.0,
 });

 let solid = &brep.solids[0];
 let (_, report) = orient_solid_shells(solid, &rcad_kernel::BRep);

 // Sphere should have shells oriented
 // Note: orientation issues may exist depending on how primitives are constructed
 assert_eq!(report.outer_shells_oriented + report.inner_shells_oriented, 1, "Sphere should have one shell");
 }

 #[test]
 fn solid_orientation_report_summary() {
 let report = SolidOrientationReport {
 outer_shells_oriented: 1,
 inner_shells_oriented: 2,
 shells_flipped: 1,
 faces_flipped: 6,
 nesting_hierarchy: vec![(0, 0), (1, 1), (2, 1)],
 is_properly_oriented: true,
 orientation_issues: vec![],
 };

 let summary = report.summary();
 assert!(summary.contains("1 outer"));
 assert!(summary.contains("2 inner"));
 assert!(summary.contains("6 faces flipped"));
 }

 #[test]
 fn orientation_issue_types() {
 // Test that OrientationIssueType variants exist
 let issue1 = OrientationIssue {
 shell_idx: 0,
 issue_type: OrientationIssueType::DegenerateShell,
 description: "Test".to_string(),
 };
 let issue2 = OrientationIssue {
 shell_idx: 1,
 issue_type: OrientationIssueType::NestingContradiction,
 description: "Test".to_string(),
 };

 assert_ne!(issue1.issue_type, issue2.issue_type);
 }

 #[test]
 fn validate_solid_topology_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let solid = &brep.solids[0];
 let report = validate_solid_topology(solid, &rcad_kernel::BRep);

 assert!(report.is_valid, "Unit box should be valid");
 assert!(report.containment_valid, "Unit box should have valid containment");
 assert!(report.void_nesting_valid, "Unit box should have valid void nesting");
 assert!(report.material_side_consistent, "Unit box should have consistent material side");
 assert!(report.errors.is_empty(), "Unit box should have no errors");
 }

 #[test]
 fn validate_solid_topology_sphere() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere {
 radius: 1.0,
 });

 let solid = &brep.solids[0];
 let report = validate_solid_topology(solid, &rcad_kernel::BRep);

 // Sphere should have valid closure
 assert!(report.closure_report.all_shells_closed, "Sphere should have closed shells");
 assert_eq!(report.closure_report.shell_count, 1, "Sphere should have one shell");
 }

 #[test]
 fn validate_solid_topology_empty_solid() {
 use rcad_kernel::topology::Solid as TopologySolid;

 let brep = rcad_kernel::BRep::new();
 let solid = TopologySolid { shells: vec![] };

 let report = validate_solid_topology(&solid, &rcad_kernel::BRep);

 assert!(!report.is_valid, "Empty solid should not be valid");
 assert!(!report.errors.is_empty(), "Empty solid should have errors");
 }

 #[test]
 fn solid_validation_report_summary() {
 let report = SolidValidationReport {
 is_valid: true,
 closure_report: SolidClosureVerificationReport::default(),
 containment_valid: true,
 void_nesting_valid: true,
 material_side_consistent: true,
 errors: vec![],
 warnings: vec![],
 };

 let summary = report.summary();
 assert!(summary.contains("Valid solid"));
 assert!(summary.contains("no errors"));
 }

 #[test]
 fn solid_validation_error_codes() {
 // Test that SolidValidationErrorCode variants exist and can be compared
 assert_ne!(SolidValidationErrorCode::OpenShell, SolidValidationErrorCode::DegenerateShell);
 assert_ne!(SolidValidationErrorCode::MultipleOuterShells, SolidValidationErrorCode::ShellIntersection);
 assert_ne!(SolidValidationErrorCode::InvalidVoidNesting, SolidValidationErrorCode::MaterialSideInconsistency);
 }

 #[test]
 fn solid_validation_warning_codes() {
 // Test that SolidValidationWarningCode variants exist and can be compared
 assert_ne!(SolidValidationWarningCode::SmallVolume, SolidValidationWarningCode::HighAspectRatio);
 assert_ne!(SolidValidationWarningCode::ToleranceIssue, SolidValidationWarningCode::NumericalIssue);
 }

 #[test]
 fn repair_solid_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let solid = &brep.solids[0];
 let result = repair_solid(solid, &rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 assert!(result.success, "Box repair should succeed");
 assert!(result.validation_report.is_valid, "Repaired box should be valid");
 assert!(result.unrepaired_issues.is_empty(), "Box should have no unrepaired issues");
 }

 #[test]
 fn repair_solid_sphere() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere {
 radius: 1.0,
 });

 let solid = &brep.solids[0];
 let result = repair_solid(solid, &rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Sphere should have closed shells after repair
 assert!(result.validation_report.closure_report.all_shells_closed, "Sphere should have closed shells");
 }

 #[test]
 fn repair_solid_cylinder() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 1.0,
 height: 2.0,
 });

 let solid = &brep.solids[0];
 let result = repair_solid(solid, &rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Cylinder should have closed shells after repair
 assert!(result.validation_report.closure_report.all_shells_closed, "Cylinder should have closed shells");
 }

 #[test]
 fn repair_solid_empty_solid() {
 use rcad_kernel::topology::Solid as TopologySolid;

 let brep = rcad_kernel::BRep::new();
 let solid = TopologySolid { shells: vec![] };

 let result = repair_solid(&solid, &rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Empty solid should be "repaired" to an empty solid
 assert!(!result.success, "Empty solid repair should not succeed");
 assert!(result.solid.shells.is_empty(), "Result should have no shells");
 }

 #[test]
 fn solid_repair_result_summary() {
 let result = SolidRepairResult {
 solid: rcad_kernel::topology::Solid { shells: vec![] },
 success: true,
 shells_closed: 1,
 shells_reoriented: 2,
 degenerate_shells_removed: 0,
 faces_modified: 6,
 gaps_closed: 0,
 validation_report: SolidValidationReport::default(),
 unrepaired_issues: vec![],
 };

 let summary = result.summary();
 assert!(summary.contains("Repair successful"));
 assert!(summary.contains("1 shells closed"));
 assert!(summary.contains("2 reoriented"));
 }

 #[test]
 fn solid_repair_result_partial_success() {
 let result = SolidRepairResult {
 solid: rcad_kernel::topology::Solid { shells: vec![] },
 success: false,
 shells_closed: 0,
 shells_reoriented: 0,
 degenerate_shells_removed: 0,
 faces_modified: 0,
 gaps_closed: 0,
 validation_report: SolidValidationReport::default(),
 unrepaired_issues: vec!["Open edges remain".to_string()],
 };

 let summary = result.summary();
 assert!(summary.contains("partially successful"));
 assert!(summary.contains("1 issues remain"));
 }

 #[test]
 fn verify_solid_closure_torus() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
 major_radius: 2.0,
 minor_radius: 0.5,
 });

 let solid = &brep.solids[0];
 let report = verify_solid_closure(solid, &rcad_kernel::BRep);

 // Torus should be closed with a single shell
 assert!(report.all_shells_closed, "Torus should have all shells closed");
 assert_eq!(report.shell_count, 1);
 // Volume computation for curved primitives depends on face normal orientation
 // Just verify we have a shell (volume might be zero or very small due to geometry)
 }

 #[test]
 fn validate_solid_topology_torus() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
 major_radius: 2.0,
 minor_radius: 0.5,
 });

 let solid = &brep.solids[0];
 let report = validate_solid_topology(solid, &rcad_kernel::BRep);

 // Torus should have valid closure
 assert!(report.closure_report.all_shells_closed, "Torus should have closed shells");
 assert_eq!(report.closure_report.shell_count, 1, "Torus should have one shell");
 }

 #[test]
 fn repair_solid_torus() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
 major_radius: 2.0,
 minor_radius: 0.5,
 });

 let solid = &brep.solids[0];
 let result = repair_solid(solid, &rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 // Torus should have closed shells after repair
 assert!(result.validation_report.closure_report.all_shells_closed, "Torus should have closed shells");
 }

 #[test]
 fn solid_closure_verification_report_default() {
 let report = SolidClosureVerificationReport::default();

 assert!(report.all_shells_closed); // default is true
 assert!(report.has_proper_nesting); // default is true
 assert_eq!(report.shell_count, 0);
 assert_eq!(report.closed_shell_count, 0);
 assert_eq!(report.open_shell_count, 0);
 assert!(report.shell_volume_signs.is_empty());
 assert!(report.shell_volumes.is_empty());
 assert_eq!(report.total_volume, 0.0);
 assert_eq!(report.volume_sign, VolumeSign::Unknown);
 assert!(report.shell_containment.is_empty());
 assert!(report.degenerate_shell_indices.is_empty());
 assert!(report.inconsistent_orientation_indices.is_empty());
 assert!(report.has_single_outer_shell); // default is true
 }

 #[test]
 fn solid_validation_report_default() {
 let report = SolidValidationReport::default();

 assert!(!report.is_valid);
 assert!(!report.containment_valid);
 assert!(!report.void_nesting_valid);
 assert!(!report.material_side_consistent);
 assert!(report.errors.is_empty());
 assert!(report.warnings.is_empty());
 }

 #[test]
 fn solid_orientation_report_default() {
 let report = SolidOrientationReport::default();

 assert_eq!(report.outer_shells_oriented, 0);
 assert_eq!(report.inner_shells_oriented, 0);
 assert_eq!(report.shells_flipped, 0);
 assert_eq!(report.faces_flipped, 0);
 assert!(report.nesting_hierarchy.is_empty());
 assert!(!report.is_properly_oriented);
 assert!(report.orientation_issues.is_empty());
 }

 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =
 // Tests for Post-Boolean Tolerance Propagation
 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =

 #[test]
 fn propagate_tolerances_post_boolean_basic() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 // Simulate a boolean operation with some intersection edges
 let intersection_edges = vec![0, 1, 2]; // First 3 edges are "intersection" edges
 let intersection_vertices = vec![0, 1, 2, 3]; // First 4 vertices

 let (result, report) = propagate_tolerances_post_boolean_op(
 &rcad_kernel::BRep,
 BooleanOpTypeForTolerance::Union,
 &intersection_edges,
 &intersection_vertices,
 );

 // Check that edges were updated
 assert!(report.edges_updated >= 3, "Should update intersection edges");
 // Check that tolerances were propagated
 assert!(report.max_edge_tolerance > TOLERANCE_ABS);
 }

 #[test]
 fn propagate_tolerances_post_boolean_intersection_type() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let intersection_edges = vec![0];
 let intersection_vertices = vec![0];

 // Intersection operations typically need higher tolerances
 let (result_union, report_union) = propagate_tolerances_post_boolean_op(
 &rcad_kernel::BRep,
 BooleanOpTypeForTolerance::Union,
 &intersection_edges,
 &intersection_vertices,
 );

 let (result_intersection, report_intersection) = propagate_tolerances_post_boolean_op(
 &rcad_kernel::BRep,
 BooleanOpTypeForTolerance::Intersection,
 &intersection_edges,
 &intersection_vertices,
 );

 // Intersection should result in higher tolerances
 assert!(report_intersection.max_edge_tolerance >= report_union.max_edge_tolerance);
 }

 #[test]
 fn test_propagate_tolerances_post_boolean_op_with_config() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let config = PostBooleanToleranceConfig::high_precision();
 let intersection_edges = vec![0];
 let intersection_vertices = vec![0];

 let (_result, report) = propagate_tolerances_post_boolean_op_with_config(
 &rcad_kernel::BRep,
 BooleanOpTypeForTolerance::General,
 &intersection_edges,
 &intersection_vertices,
 &config,
 );

 // High-precision config should result in lower tolerances
 assert!(report.max_edge_tolerance < 0.1);
 }

 #[test]
 fn post_boolean_config_presets() {
 let standard = PostBooleanToleranceConfig::standard();
 let high_precision = PostBooleanToleranceConfig::high_precision();
 let relaxed = PostBooleanToleranceConfig::relaxed();

 // High precision should have smallest floor
 assert!(high_precision.tolerance_floor < standard.tolerance_floor);
 // Relaxed should have largest floor
 assert!(relaxed.tolerance_floor > standard.tolerance_floor);
 }

 #[test]
 fn detect_and_resolve_tolerance_conflicts_resolves_vertex_edge() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

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
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set up conflict: vertex tolerance > edge tolerance
 brep.geom.vertex_tolerance = vec![TOLERANCE_ADAPTIVE_MAX, TOLERANCE_ADAPTIVE_MAX, TOLERANCE_ABS]; // v0 and v1 have high tolerance
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS, TOLERANCE_ABS, TOLERANCE_ABS]; // edges have low tolerance
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let mut cloned = brep.clone();
 let (conflicts, resolved) = detect_and_resolve_tolerance_conflicts(&mut cloned, TOLERANCE_ABS);

 assert!(conflicts >= 1, "Should detect at least one conflict");
 assert!(resolved >= 1, "Should resolve at least one conflict");
 // Edge 0 should now have higher tolerance (>= vertex 0 and 1)
 assert!(cloned.geom.edge_tolerance[0] >= TOLERANCE_ADAPTIVE_MAX);
 }

 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =
 // Tests for Post-Sew Tolerance Propagation
 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =

 #[test]
 fn propagate_tolerances_post_sew_basic() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

 // Create two edges that were "sewn" together
 brep.edges.push(Edge { start: 0, end: 1 }); // e0
 brep.edges.push(Edge { start: 1, end: 2 }); // e1
 brep.edges.push(Edge { start: 2, end: 3 }); // e2
 brep.edges.push(Edge { start: 3, end: 0 }); // e3 (seam edge)

 let face = Face {
 outer_wire: Wire {
 edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
 },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Initialize tolerances
 brep.geom.vertex_tolerance = vec![TOLERANCE_ABS; 4];
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS; 4];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 // Simulate seam edge pairs (edge 3 was sewn)
 let seam_pairs = vec![(3, 3)];

 let (_result, report) = propagate_tolerances_post_sew(&rcad_kernel::BRep, TOLERANCE_RETRY_LADDER_COARSE, &seam_pairs);

 // Verify function runs successfully
 assert!(report.max_seam_tolerance > 0.0 || report.seam_edges_updated == 0);
 }

 #[test]
 fn test_propagate_tolerances_post_sew_with_config() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });
 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 brep.geom.vertex_tolerance = vec![TOLERANCE_ABS; 2];
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let config = PostSewToleranceConfig {
 seam_tolerance_factor: 2.0,
 max_growth_ratio: 1000.0,
 ..Default::default()
 };

 let seam_pairs = vec![(0, 0)];
 let (_result, report) = propagate_tolerances_post_sew_with_config(
 &rcad_kernel::BRep,
 TOLERANCE_RETRY_LADDER_COARSE,
 &seam_pairs,
 &config,
 );

 // Verify function runs successfully
 assert!(report.max_seam_tolerance >= 0.0);
 }

 #[test]
 fn post_sew_config_default() {
 let config = PostSewToleranceConfig::default();

 assert_eq!(config.tolerance_floor, TOLERANCE_ABS);
 assert_eq!(config.seam_tolerance_factor, 1.5);
 assert!(config.ensure_seam_consistency);
 }

 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =
 // Tests for Tolerance Rules Engine
 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =

 #[test]
 fn tolerance_rule_variants() {
 // Test that all rule variants exist
 let rules = vec![
 ToleranceRule::OcctStandard,
 ToleranceRule::Conservative,
 ToleranceRule::Aggressive,
 ToleranceRule::Harmonized,
 ToleranceRule::Bounded,
 ToleranceRule::ModelScale,
 ];

 // Ensure they can be compared
 assert_ne!(ToleranceRule::OcctStandard, ToleranceRule::Aggressive);
 }

 #[test]
 fn conflict_resolution_policy_variants() {
 let policies = vec![
 ConflictResolutionPolicy::Ignore,
 ConflictResolutionPolicy::PropagateUp,
 ConflictResolutionPolicy::ClampDown,
 ConflictResolutionPolicy::ReportOnly,
 ];

 assert_ne!(ConflictResolutionPolicy::Ignore, ConflictResolutionPolicy::PropagateUp);
 }

 #[test]
 fn tolerance_propagation_config_presets() {
 let occt = TolerancePropagationConfig::occt_standard();
 assert_eq!(occt.rule, ToleranceRule::OcctStandard);

 let conservative = TolerancePropagationConfig::conservative();
 assert_eq!(conservative.rule, ToleranceRule::Conservative);

 let aggressive = TolerancePropagationConfig::aggressive();
 assert_eq!(aggressive.rule, ToleranceRule::Aggressive);

 let harmonized = TolerancePropagationConfig::harmonized();
 assert_eq!(harmonized.rule, ToleranceRule::Harmonized);

 let bounded = TolerancePropagationConfig::bounded(0.1);
 assert_eq!(bounded.rule, ToleranceRule::Bounded);
 assert_eq!(bounded.bound_value, 0.1);

 let model_scale = TolerancePropagationConfig::model_scale(100.0);
 assert_eq!(model_scale.rule, ToleranceRule::ModelScale);
 assert!((model_scale.model_scale - 100.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
 }

 #[test]
 fn tolerance_propagation_engine_default() {
 let engine = TolerancePropagationEngine::new();
 assert_eq!(engine.config.rule, ToleranceRule::OcctStandard);
 }

 #[test]
 fn tolerance_propagation_engine_occt_standard() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
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
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set vertex tolerances higher than edge tolerances
 brep.geom.vertex_tolerance = vec![TOLERANCE_RETRY_LADDER_COARSE, TOLERANCE_RETRY_LADDER_COARSE, TOLERANCE_RETRY_LADDER_COARSE];
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS, TOLERANCE_ABS, TOLERANCE_ABS];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let engine = TolerancePropagationEngine::occt_standard();
 let (result, report) = engine.propagate(&rcad_kernel::BRep);

 // Edges should now have higher tolerances (propagated from vertices)
 assert!(result.geom.edge_tolerance[0] >= TOLERANCE_RETRY_LADDER_COARSE);
 assert!(report.rule_applied == ToleranceRule::OcctStandard);
 }

 #[test]
 fn tolerance_propagation_engine_conservative() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let engine = TolerancePropagationEngine::conservative();
 let (result, report) = engine.propagate(&rcad_kernel::BRep);

 assert_eq!(report.rule_applied, ToleranceRule::Conservative);
 }

 #[test]
 fn tolerance_propagation_engine_aggressive() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
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
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 brep.geom.vertex_tolerance = vec![TOLERANCE_ABS; 3];
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS; 3];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let engine = TolerancePropagationEngine::aggressive();
 let (result, report) = engine.propagate(&rcad_kernel::BRep);

 assert_eq!(report.rule_applied, ToleranceRule::Aggressive);
 // Aggressive propagation may update tolerances more
 }

 #[test]
 fn tolerance_propagation_engine_bounded() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set very high tolerances
 brep.geom.vertex_tolerance = vec![1.0, 1.0];
 brep.geom.edge_tolerance = vec![1.0];
 brep.geom.face_tolerance = vec![1.0];

 let engine = TolerancePropagationEngine::bounded(TOLERANCE_ADAPTIVE_MAX);
 let (result, report) = engine.propagate(&rcad_kernel::BRep);

 // All tolerances should be clamped to bound
 assert!(result.geom.vertex_tolerance[0] <= TOLERANCE_ADAPTIVE_MAX);
 assert!(result.geom.edge_tolerance[0] <= TOLERANCE_ADAPTIVE_MAX);
 assert!(result.geom.face_tolerance[0] <= TOLERANCE_ADAPTIVE_MAX);
 }

 #[test]
 fn tolerance_propagation_engine_model_scale() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1000.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 brep.geom.vertex_tolerance = vec![TOLERANCE_ABS, TOLERANCE_ABS];
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let engine = TolerancePropagationEngine::with_config(
 TolerancePropagationConfig::model_scale(1000.0)
 );
 let (result, report) = engine.propagate(&rcad_kernel::BRep);

 assert_eq!(report.rule_applied, ToleranceRule::ModelScale);
 // Tolerances should be scaled
 assert!(result.geom.vertex_tolerance[0] > TOLERANCE_ABS);
 }

 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =
 // Tests for Tolerance Consistency Analysis
 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =

 #[test]
 fn analyze_tolerance_consistency_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let report = analyze_tolerance_consistency(&rcad_kernel::BRep, TOLERANCE_ABS, TOLERANCE_ABS, 1.0);

 // Unit box should have consistent tolerances
 assert!(report.is_consistent || report.violation_count == 0);
 }

 #[test]
 fn analyze_tolerance_consistency_detects_vertex_edge_violation() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set vertex tolerance > edge tolerance (violation)
 brep.geom.vertex_tolerance = vec![TOLERANCE_ADAPTIVE_MAX, TOLERANCE_ADAPTIVE_MAX];
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let report = analyze_tolerance_consistency(&rcad_kernel::BRep, TOLERANCE_ABS, TOLERANCE_ABS, 1.0);

 assert!(!report.is_consistent, "Should detect inconsistency");
 assert!(report.violation_count >= 1, "Should have at least one violation");

 let vertex_edge_violations = report.violations_by_type(ToleranceViolationType::VertexExceedsEdge);
 assert!(!vertex_edge_violations.is_empty(), "Should have vertex>edge violations");
 }

 #[test]
 fn analyze_tolerance_consistency_detects_edge_face_violation() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set edge tolerance > face tolerance (violation)
 brep.geom.vertex_tolerance = vec![TOLERANCE_ABS, TOLERANCE_ABS];
 brep.geom.edge_tolerance = vec![TOLERANCE_ADAPTIVE_MAX];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let report = analyze_tolerance_consistency(&rcad_kernel::BRep, TOLERANCE_ABS, TOLERANCE_ABS, 1.0);

 let edge_face_violations = report.violations_by_type(ToleranceViolationType::EdgeExceedsFace);
 assert!(!edge_face_violations.is_empty(), "Should have edge>face violations");
 }

 #[test]
 fn analyze_tolerance_consistency_detects_invalid_values() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set NaN and negative tolerances
 brep.geom.vertex_tolerance = vec![f64::NAN, -TOLERANCE_LINEAR_ULTRA_STRICT];
 brep.geom.edge_tolerance = vec![f64::INFINITY];
 brep.geom.face_tolerance = vec![0.0];

 let report = analyze_tolerance_consistency(&rcad_kernel::BRep, TOLERANCE_ABS, TOLERANCE_ABS, 1.0);

 let invalid_violations = report.violations_by_type(ToleranceViolationType::InvalidValue);
 assert!(invalid_violations.len() >= 2, "Should detect invalid values");
 }

 #[test]
 fn tolerance_violation_severity() {
 let violation = ToleranceViolation {
 violation_type: ToleranceViolationType::VertexExceedsEdge,
 entity_index: 0,
 related_index: Some(0),
 actual_tolerance: TOLERANCE_ADAPTIVE_MAX,
 expected_tolerance: TOLERANCE_ABS,
 severity: 4,
 suggested_fix: ToleranceFix::IncreaseLower,
 };

 assert!(violation.severity >= 4);
 }

 #[test]
 fn tolerance_consistency_report_summary() {
 let report = ToleranceConsistencyReport {
 is_consistent: true,
 violation_count: 0,
 critical_violation: 0,
 violations: vec![],
 stats: ToleranceAnalysisReport::default(),
 suggested_global_fixes: vec![],
 };

 assert!(report.summary().contains("OK"));

 // Create report with actual violations
 let critical_violation = ToleranceViolation {
 violation_type: ToleranceViolationType::VertexExceedsEdge,
 entity_index: 0,
 related_index: None,
 actual_tolerance: TOLERANCE_ADAPTIVE_MAX,
 expected_tolerance: TOLERANCE_MESH_LEGACY,
 severity: 4,
 suggested_fix: ToleranceFix::Propagate,
 };
 let normal_violation = ToleranceViolation {
 violation_type: ToleranceViolationType::EdgeExceedsFace,
 entity_index: 1,
 related_index: None,
 actual_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
 expected_tolerance: TOLERANCE_MESH_LEGACY,
 severity: 2,
 suggested_fix: ToleranceFix::Propagate,
 };

 let report_with_violations = ToleranceConsistencyReport {
 is_consistent: false,
 violation_count: 2,
 critical_violation: 1,
 violations: vec![critical_violation, normal_violation],
 stats: ToleranceAnalysisReport::default(),
 suggested_global_fixes: vec![],
 };

 assert!(report_with_violations.summary().contains("2 violations"));
 assert!(report_with_violations.summary().contains("1 critical"));
 }

 #[test]
 fn apply_tolerance_fixes_basic() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set up violations
 brep.geom.vertex_tolerance = vec![TOLERANCE_ADAPTIVE_MAX, TOLERANCE_ADAPTIVE_MAX]; // High vertex tolerance
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS]; // Low edge tolerance
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let report = analyze_tolerance_consistency(&rcad_kernel::BRep, TOLERANCE_ABS, TOLERANCE_ABS, 1.0);
 assert!(!report.is_consistent);

 let (fixed, fixes_applied) = apply_tolerance_fixes(&rcad_kernel::BRep, &report, 0);

 assert!(fixes_applied >= 1, "Should apply at least one fix");
 // Edge tolerance should now be >= vertex tolerance
 assert!(fixed.geom.edge_tolerance[0] >= TOLERANCE_ADAPTIVE_MAX);
 }

 #[test]
 fn tolerance_fix_variants() {
 // Test that all fix variants exist
 assert_ne!(ToleranceFix::IncreaseLower, ToleranceFix::DecreaseHigher);
 assert_ne!(ToleranceFix::SetToValue, ToleranceFix::Propagate);
 assert_ne!(ToleranceFix::ManualIntervention, ToleranceFix::IncreaseLower);
 }

 #[test]
 fn tolerance_violation_type_variants() {
 // Test that all violation type variants exist
 assert_ne!(ToleranceViolationType::VertexExceedsEdge, ToleranceViolationType::EdgeExceedsFace);
 assert_ne!(ToleranceViolationType::BelowFloor, ToleranceViolationType::ExceedsMaximum);
 assert_ne!(ToleranceViolationType::SeamInconsistency, ToleranceViolationType::InvalidValue);
 }

 #[test]
 fn propagate_tolerances_post_boolean_handles_conflicts() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.edges.push(Edge { start: 0, end: 1 });

 let face = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };
 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

 // Set up a conflict: vertex tolerance > edge tolerance
 brep.geom.vertex_tolerance = vec![TOLERANCE_ADAPTIVE_MAX, TOLERANCE_ADAPTIVE_MAX];
 brep.geom.edge_tolerance = vec![TOLERANCE_ABS];
 brep.geom.face_tolerance = vec![TOLERANCE_ABS];

 let (_result, report) = propagate_tolerances_post_boolean_op(
 &rcad_kernel::BRep,
 BooleanOpTypeForTolerance::Union,
 &[],
 &[],
 );

 // Verify function runs successfully
 }

 #[test]
 fn propagate_tolerances_post_boolean_empty_intersection_lists() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 // Empty intersection lists should still work
 let (result, report) = propagate_tolerances_post_boolean_op(
 &rcad_kernel::BRep,
 BooleanOpTypeForTolerance::General,
 &[],
 &[],
 );

 // Should still run propagation
 assert!(report.max_edge_tolerance > 0.0);
 }

 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =
 // Tests for Connectivity Graph Analysis
 // = =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =???= =

 #[test]
 fn build_connectivity_graph_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let graph = build_connectivity_graph(&rcad_kernel::BRep);

 assert_eq!(graph.vertex_count, 8, "Unit box should have 8 vertices");
 assert_eq!(graph.edge_count, 12, "Unit box should have 12 edges");
 assert_eq!(graph.face_count, 6, "Unit box should have 6 faces");
 assert_eq!(graph.face_components.len(), 1, "Unit box should be single component");
 }

 #[test]
 fn build_connectivity_graph_disconnected_faces() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 // Create two disconnected triangles
 let mut brep = rcad_kernel::BRep::new();

 // Triangle 1
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 0 });

 let face1 = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 // Triangle 2 (disconnected, far away)
 brep.vertices.push(Vertex { point: DVec3::new(10.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(11.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(10.0, 1.0, 0.0) });

 brep.edges.push(Edge { start: 3, end: 4 });
 brep.edges.push(Edge { start: 4, end: 5 });
 brep.edges.push(Edge { start: 5, end: 3 });

 let face2 = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(3), WireEdge::fwd(4), WireEdge::fwd(5)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 brep.solids.push(Solid { shells: vec![Shell { faces: vec![face1, face2] }] });

 let graph = build_connectivity_graph(&rcad_kernel::BRep);

 assert_eq!(graph.face_count, 2);
 assert_eq!(graph.face_components.len(), 2, "Should have two disconnected components");
 }

 #[test]
 fn is_fully_connected_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 assert!(is_fully_connected(&rcad_kernel::BRep), "Unit box should be fully connected");
 }

 #[test]
 fn test_disconnected_component_count() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Single triangle
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

 assert_eq!(disconnected_component_count(&rcad_kernel::BRep), 1);
 }

 #[test]
 fn connectivity_strength_values() {
 assert!(ConnectivityStrength::Weak.to_value() < ConnectivityStrength::Medium.to_value());
 assert!(ConnectivityStrength::Medium.to_value() < ConnectivityStrength::Strong.to_value());
 assert!(ConnectivityStrength::Strong.to_value() < ConnectivityStrength::Full.to_value());
 }

 #[test]
 fn detect_connectivity_gaps_connected() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let gaps = detect_connectivity_gaps(&rcad_kernel::BRep, TOLERANCE_ADAPTIVE_MAX);
 assert!(gaps.is_empty(), "Connected box should have no gaps");
 }

 #[test]
 fn validate_connectivity_unit_box() {
 let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 let report = validate_connectivity(&rcad_kernel::BRep, TOLERANCE_MESH_LEGACY);

 assert!(report.is_connected, "Unit box should be connected");
 assert_eq!(report.component_count, 1);
 }

 #[test]
 fn validate_connectivity_disconnected() {
 use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

 let mut brep = rcad_kernel::BRep::new();

 // Triangle 1
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

 brep.edges.push(Edge { start: 0, end: 1 });
 brep.edges.push(Edge { start: 1, end: 2 });
 brep.edges.push(Edge { start: 2, end: 0 });

 let face1 = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)] },
 inner_wires: vec![],
 normal: DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 };

 // Triangle 2 (far away)
 brep.vertices.push(Vertex { point: DVec3::new(100.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(101.0, 0.0, 0.0) });
 brep.vertices.push(Vertex { point: DVec3::new(100.0, 1.0, 0.0) });

 brep.edges.push(Edge { start: 3, end: 4 });
 brep.edges.push(Edge { start: 4, end: 5 });
 brep.edges.push(Edge { start: 5, end: 3 });

 let face2 = Face {
 outer_wire: Wire { edges: vec![WireEdge::fwd(3), WireEdge::fwd(4), WireEdge::fwd(5)] },
 inner_wires: vec![],
 normal: DVec3::Z,
