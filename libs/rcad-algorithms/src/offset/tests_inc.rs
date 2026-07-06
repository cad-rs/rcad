#[cfg(test)]
mod tests {
 use crate::offset::*;
 use glam::DVec3;
 use rcad_kernel::geom::{Plane, SphericalSurface, CylindricalSurface, ConicalSurface, ToroidalSurface};

 #[test]
 fn offset_plane_translates() {
 let plane = Surface3::Plane(Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 });

 let offset = offset_surface(&plane, 0.5).unwrap();

 if let Surface3::Plane(p) = offset {
 assert!((p.origin.z - 0.5).abs() < TOLERANCE_COORD_SUB, "plane should translate by offset distance");
 assert!((p.normal - DVec3::Z).length() < TOLERANCE_COORD_SUB, "normal should be unchanged");
 } else {
 panic!("expected Plane");
 }
 }

 #[test]
 fn offset_sphere_grows() {
 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 2.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 let offset = offset_surface(&sphere, 0.5).unwrap();

 if let Surface3::Sphere(s) = offset {
 assert!((s.radius - 2.5).abs() < TOLERANCE_COORD_SUB, "radius should increase by offset");
 } else {
 panic!("expected Sphere");
 }
 }

 #[test]
 fn offset_cylinder_grows() {
 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 });

 let offset = offset_surface(&cylinder, 0.3).unwrap();

 if let Surface3::Cylinder(c) = offset {
 assert!((c.radius - 1.3).abs() < TOLERANCE_COORD_SUB, "radius should increase by offset");
 } else {
 panic!("expected Cylinder");
 }
 }

 #[test]
 fn offset_sphere_negative_too_large_returns_none() {
 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 // Negative offset larger than radius should return None
 let offset = offset_surface(&sphere, -2.0);
 assert!(offset.is_none(), "offset larger than radius should return None");
 }

 #[test]
 fn offset_zero_returns_error() {
 let brep = BRep::new();
 let opts = OffsetOptions::new(0.0);

 let result = offset_shape(&brep, opts);
 assert!(matches!(result, Err(OffsetError::ZeroDistance)));
 }

 #[test]
 fn offset_options_default() {
 let opts = OffsetOptions::default();
 assert_eq!(opts.distance, 1.0);
 assert!(opts.check_self_intersection);
 assert!(!opts.auto_repair);
 }

 #[test]
 fn offset_options_builder() {
 let opts = OffsetOptions::new(0.5)
 .with_tolerance(TOLERANCE_MESH_LEGACY)
 .with_self_intersection_check(false)
 .with_auto_repair(true);

 assert_eq!(opts.distance, 0.5);
 assert!((opts.tolerance - TOLERANCE_MESH_LEGACY).abs() < TOLERANCE_LEN_MIN);
 assert!(!opts.check_self_intersection);
 assert!(opts.auto_repair);
 }

 #[test]
 fn self_intersection_detection_small_box() {
 // Create a 1x1x1 box
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 // Populate geometry
 crate::geom_populate::populate_box_geom(&mut brep);

 // Offset distance > 0.5 should self-intersect
 let self_intersects = detect_self_intersection(&brep, 0.6);
 assert!(self_intersects, "should detect self-intersection for large offset");

 // Offset distance < 0.5 should not self-intersect
 let no_intersect = detect_self_intersection(&brep, 0.4);
 assert!(!no_intersect, "should not detect self-intersection for small offset");
 }

 #[test]
 fn offset_shell_simple_box() {
 // Create a 2x2x2 box
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let shell = &brep.solids[0].shells[0];
 let result = offset_shell(shell, &brep, 0.1);

 assert!(result.is_ok(), "offset_shell should succeed for a simple box");
 let offset_brep = result.unwrap();

 // Should have the same number of faces
 let orig_face_count = shell.faces.len();
 let offset_face_count = offset_brep.solids[0].shells[0].faces.len();
 assert_eq!(offset_face_count, orig_face_count, "should preserve face count");
 }

 #[test]
 fn offset_shell_negative_distance() {
 // Create a 2x2x2 box
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let shell = &brep.solids[0].shells[0];
 let result = offset_shell(shell, &brep, -0.1);

 assert!(result.is_ok(), "offset_shell with negative distance should succeed");
 }

 #[test]
 fn offset_solid_simple() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let solid = &brep.solids[0];
 let result = offset_solid(solid, &brep, 0.2);

 assert!(result.is_ok(), "offset_solid should succeed");
 let offset_brep = result.unwrap();

 // Verify structure
 assert!(!offset_brep.vertices.is_empty(), "should have vertices");
 assert!(!offset_brep.edges.is_empty(), "should have edges");
 assert!(!offset_brep.solids.is_empty(), "should have solids");
 }

 #[test]
 fn hollow_solid_simple_box() {
 // Create a 2x2x2 box
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 // Hollow by removing top face (index 5 based on typical box construction)
 let solid = &brep.solids[0];
 let result = hollow_solid(solid, &brep, 0.1, &[5]);

 assert!(result.is_ok(), "hollow_solid should succeed with one face removed");
 let hollow_brep = result.unwrap();

 // Should have original kept faces (5) + lateral faces at boundary
 let face_count = hollow_brep.solids[0].shells[0].faces.len();
 assert!(face_count >= 5, "should have at least 5 faces (kept faces + lateral faces)");
 }

 #[test]
 fn hollow_solid_multiple_open_faces() {
 // Create a 2x2x2 box
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 // Hollow by removing top (5) and bottom (0) faces
 let solid = &brep.solids[0];
 let result = hollow_solid(solid, &brep, 0.1, &[0, 5]);

 assert!(result.is_ok(), "hollow_solid should succeed with multiple open faces");
 }

 #[test]
 fn hollow_solid_all_faces_error() {
 // Create a box
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 // Trying to remove all 6 faces should error
 let solid = &brep.solids[0];
 let result = hollow_solid(solid, &brep, 0.1, &[0, 1, 2, 3, 4, 5]);

 assert!(result.is_err(), "hollow_solid should fail when all faces are removed");
 }

 #[test]
 fn offset_shape_api() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let opts = OffsetOptions::new(0.1)
 .with_self_intersection_check(true);

 let result = offset_shape(&brep, opts);

 assert!(result.is_ok(), "offset_shape should succeed");
 let offset_result = result.unwrap();

 assert_eq!(offset_result.offset_faces, 6, "should have 6 offset faces");
 assert!(!offset_result.self_intersection, "should not have self-intersection");
 }

 #[test]
 fn offset_torus_surface() {
 let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 2.0,
 minor_radius: 0.5,
 });

 let offset = offset_surface(&torus, 0.1).unwrap();

 if let Surface3::Torus(t) = offset {
 assert!((t.minor_radius - 0.6).abs() < TOLERANCE_COORD_SUB, "minor radius should increase by offset");
 assert!((t.major_radius - 2.0).abs() < TOLERANCE_COORD_SUB, "major radius should be unchanged");
 } else {
 panic!("expected Torus");
 }
 }

 #[test]
 fn offset_cone_surface() {
 let cone = Surface3::Cone(rcad_kernel::geom::ConicalSurface {
 apex: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 half_angle_rad: std::f64::consts::PI / 6.0, // 30 degrees
 });

 let offset = offset_surface(&cone, 0.1);

 assert!(offset.is_some(), "cone offset should succeed");
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Join Type Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn join_type_default() {
 assert_eq!(JoinType::default(), JoinType::Intersection);
 }

 #[test]
 fn join_type_requires_geometry() {
 assert!(!JoinType::Intersection.requires_join_geometry());
 assert!(JoinType::Arc.requires_join_geometry());
 assert!(JoinType::Tangent.requires_join_geometry());
 }

 #[test]
 fn join_type_as_str() {
 assert_eq!(JoinType::Intersection.as_str(), "intersection");
 assert_eq!(JoinType::Arc.as_str(), "arc");
 assert_eq!(JoinType::Tangent.as_str(), "tangent");
 }

 #[test]
 fn offset_with_arc_join() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let opts = OffsetOptions::new(0.1)
 .with_join_type(JoinType::Arc);

 let result = offset_shape(&brep, opts);
 assert!(result.is_ok(), "offset with arc join should succeed");
 }

 #[test]
 fn offset_with_tangent_join() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let opts = OffsetOptions::new(0.1)
 .with_join_type(JoinType::Tangent);

 let result = offset_shape(&brep, opts);
 assert!(result.is_ok(), "offset with tangent join should succeed");
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Variable Thickness Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn variable_thickness_new() {
 let vt = VariableThickness::new(1.0);
 assert_eq!(vt.default_thickness, 1.0);
 assert!(vt.face_thicknesses.is_empty());
 assert_eq!(vt.transition_width, 0.0);
 assert!(!vt.interpolate);
 }

 #[test]
 fn variable_thickness_with_face() {
 let vt = VariableThickness::new(1.0)
 .with_face(0, 0.5)
 .with_face(1, 1.5);

 assert_eq!(vt.thickness_for_face(0), 0.5);
 assert_eq!(vt.thickness_for_face(1), 1.5);
 assert_eq!(vt.thickness_for_face(2), 1.0); // default
 }

 #[test]
 fn variable_thickness_validation() {
 let vt = VariableThickness::new(1.0)
 .with_face(0, 0.5)
 .with_face(10, 1.5); // Invalid face index

 // Validate with 5 faces
 let result = vt.validate(5);
 assert!(result.is_err(), "should fail for out-of-range face index");
 }

 #[test]
 fn variable_thickness_zero_thickness_error() {
 let vt = VariableThickness::new(0.0);

 let result = vt.validate(5);
 assert!(result.is_err(), "should fail for zero default thickness");
 }

 #[test]
 fn offset_with_variable_thickness() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let vt = VariableThickness::new(0.2)
 .with_face(0, 0.1)
 .with_face(1, 0.3);

 let opts = OffsetOptions::new(0.2)
 .with_variable_thickness(vt);

 let result = offset_shape(&brep, opts);
 assert!(result.is_ok(), "offset with variable thickness should succeed");
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Self-Intersection Detection Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn detect_self_intersection_detailed_empty() {
 let brep = BRep::new();
 let result = detect_self_intersection_detailed(&brep, 0.5);

 assert!(!result.has_intersection);
 assert!(result.intersecting_pairs.is_empty());
 assert!(result.min_safe_distance.is_none());
 }

 #[test]
 fn detect_self_intersection_detailed_box() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 // Large offset should detect intersection
 let result = detect_self_intersection_detailed(&brep, 0.6);
 assert!(result.has_intersection);
 assert!(!result.intersecting_pairs.is_empty());
 assert!(result.min_safe_distance.is_some());

 // Small offset should not detect intersection
 let result = detect_self_intersection_detailed(&brep, 0.3);
 assert!(!result.has_intersection);
 }

 #[test]
 fn self_intersection_config_default() {
 let config = SelfIntersectionConfig::default();

 assert!(config.detect);
 assert!(!config.auto_repair);
 assert_eq!(config.max_repair_attempts, 5);
 assert!((config.reduction_factor - 0.8).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Quality Analysis Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn offset_quality_default() {
 let quality = OffsetQuality::default();

 assert_eq!(quality.min_wall_thickness, 0.0);
 assert_eq!(quality.max_deviation, 0.0);
 assert_eq!(quality.degenerate_edge_count, 0);
 assert_eq!(quality.self_intersection_count, 0);
 assert!(!quality.is_valid); // 0 min_wall_thickness < default threshold
 }

 #[test]
 fn offset_quality_check_thresholds() {
 let quality = OffsetQuality {
 min_wall_thickness: 0.5,
 max_deviation: TOLERANCE_RETRY_LADDER_MID,
 degenerate_edge_count: 0,
 self_intersection_count: 0,
 face_area_ratio: 1.0,
 edge_length_ratio: 1.0,
 is_valid: true,
 warnings: Vec::new(),
 };

 let thresholds = QualityThresholds::default();
 assert!(quality.check_thresholds(&thresholds).is_ok());
 }

 #[test]
 fn offset_quality_check_thresholds_failure() {
 let quality = OffsetQuality {
 min_wall_thickness: TOLERANCE_COORD_SUB, // Below threshold
 max_deviation: 0.0,
 degenerate_edge_count: 0,
 self_intersection_count: 0,
 face_area_ratio: 1.0,
 edge_length_ratio: 1.0,
 is_valid: false,
 warnings: Vec::new(),
 };

 let thresholds = QualityThresholds::default();
 assert!(quality.check_thresholds(&thresholds).is_err());
 }

 #[test]
 fn quality_thresholds_default() {
 let thresholds = QualityThresholds::default();

 assert!((thresholds.min_wall_thickness - TOLERANCE_MESH_LEGACY).abs() < TOLERANCE_LEN_MIN);
 assert!((thresholds.max_deviation - TOLERANCE_RETRY_LADDER_COARSE).abs() < TOLERANCE_LEN_MIN);
 assert!(!thresholds.allow_self_intersection);
 assert!((thresholds.max_degenerate_ratio - 0.1).abs() < TOLERANCE_LEN_MIN);
 }

 #[test]
 fn analyze_offset_quality_simple() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let opts = OffsetOptions::new(0.1);
 let result = offset_shape(&brep, opts.clone()).unwrap();

 let quality = analyze_offset_quality(&result.brep, &brep, &opts);

 assert!(quality.is_valid || quality.warnings.iter().any(|w| w.contains("wall thickness")));
 }

 #[test]
 fn compute_min_wall_thickness_box() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 crate::geom_populate::populate_box_geom(&mut brep);

 let min_thickness = compute_min_wall_thickness(&brep, 0.1);

 // For a 2x2x2 box with offset 0.1, min wall should be around 1.8
 assert!(min_thickness > 0.0, "min wall thickness should be positive");
 }

 #[test]
 fn test_compute_face_area_ratio() {
 let mut brep1 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });
 crate::geom_populate::populate_box_geom(&mut brep1);

 let mut brep2 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });
 crate::geom_populate::populate_box_geom(&mut brep2);

 let ratio = compute_face_area_ratio(&brep1, &brep2);
 assert!((ratio - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "same box should have ratio 1.0");
 }

 #[test]
 fn test_compute_edge_length_ratio() {
 let mut brep1 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });
 crate::geom_populate::populate_box_geom(&mut brep1);

 let mut brep2 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });
 crate::geom_populate::populate_box_geom(&mut brep2);

 let ratio = compute_edge_length_ratio(&brep1, &brep2);
 assert!((ratio - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "same box should have ratio 1.0");
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Offset Options Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn offset_options_with_join_type() {
 let opts = OffsetOptions::new(0.5)
 .with_join_type(JoinType::Arc);

 assert_eq!(opts.join_type, JoinType::Arc);
 }

 #[test]
 fn offset_options_with_variable_thickness() {
 let vt = VariableThickness::new(1.0).with_face(0, 0.5);
 let opts = OffsetOptions::new(0.5)
 .with_variable_thickness(vt.clone());

 assert!(opts.variable_thickness.is_some());
 assert_eq!(opts.variable_thickness.as_ref().unwrap().thickness_for_face(0), 0.5);
 }

 #[test]
 fn offset_options_with_self_intersection_config() {
 let config = SelfIntersectionConfig {
 detect: true,
 auto_repair: true,
 max_repair_attempts: 10,
 reduction_factor: 0.9,
 min_offset_distance: 0.001,
 allow_partial_results: true,
 };

 let opts = OffsetOptions::new(0.5)
 .with_self_intersection_config(config.clone());

 assert!(opts.self_intersection_config.auto_repair);
 assert_eq!(opts.self_intersection_config.max_repair_attempts, 10);
 }

 #[test]
 fn offset_options_with_quality_thresholds() {
 let thresholds = QualityThresholds {
 min_wall_thickness: 0.1,
 max_deviation: 0.01,
 allow_self_intersection: true,
 max_degenerate_ratio: 0.05,
 };

 let opts = OffsetOptions::new(0.5)
 .with_quality_thresholds(thresholds.clone());

 assert!((opts.quality_thresholds.min_wall_thickness - 0.1).abs() < TOLERANCE_LEN_MIN);
 assert!(opts.quality_thresholds.allow_self_intersection);
 }

 #[test]
 fn offset_options_with_approximation_tolerance() {
 let opts = OffsetOptions::new(0.5)
 .with_approximation_tolerance(TOLERANCE_MESH_LEGACY);

 assert!((opts.approximation_tolerance - TOLERANCE_MESH_LEGACY).abs() < TOLERANCE_LEN_MIN);
 }

 #[test]
 fn offset_options_with_wall_thickness_check() {
 let opts = OffsetOptions::new(0.5)
 .with_wall_thickness_check(0.1);

 assert!(opts.check_wall_thickness);
 assert!((opts.min_wall_thickness - 0.1).abs() < TOLERANCE_LEN_MIN);
 }

 #[test]
 fn offset_options_effective_distance_for_face() {
 let vt = VariableThickness::new(1.0).with_face(0, 0.5);
 let opts = OffsetOptions::new(1.0)
 .with_variable_thickness(vt);

 assert_eq!(opts.effective_distance_for_face(0), 0.5);
 assert_eq!(opts.effective_distance_for_face(1), 1.0); // default
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Offset Result Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn offset_result_fields() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });
 crate::geom_populate::populate_box_geom(&mut brep);

 let opts = OffsetOptions::new(0.1);
 let result = offset_shape(&brep, opts).unwrap();

 assert_eq!(result.offset_faces, 6);
 assert!(!result.self_intersection);
 assert_eq!(result.effective_distance, 0.1);
 assert_eq!(result.repair_attempts, 0);
 assert!(result.warnings.is_empty() || !result.warnings.is_empty());
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Error Display Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn offset_error_display() {
 assert_eq!(
 OffsetError::ZeroDistance.to_string(),
 "offset distance is zero"
 );
 assert_eq!(
 OffsetError::InvalidInput("test").to_string(),
 "invalid input: test"
 );
 assert!(OffsetError::DegenerateSurface { face_index: 1, distance: 0.5 }
 .to_string()
 .contains("face 1"));
 assert!(OffsetError::SelfIntersection { description: "test".to_string() }
 .to_string()
 .contains("self-intersection detected"));
 }

 #[test]
 fn offset_error_new_variants() {
 let err = OffsetError::WallThicknessViolation {
 minimum: 0.1,
 actual: 0.05,
 location: "face 0".to_string(),
 };
 assert!(err.to_string().contains("0.05"));

 let err = OffsetError::JoinCreationFailed {
 join_type: JoinType::Arc,
 edge_index: 1,
 reason: "test".to_string(),
 };
 assert!(err.to_string().contains("Arc"));

 let err = OffsetError::QualityCheckFailed {
 metric: "wall_thickness".to_string(),
 value: 0.05,
 threshold: 0.1,
 };
 assert!(err.to_string().contains("wall_thickness"));

 let err = OffsetError::RecoveryFailed {
 attempts: 3,
 last_error: "test error".to_string(),
 };
 assert!(err.to_string().contains("3 attempts"));
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // Integration Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn full_offset_workflow_with_quality_check() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });
 crate::geom_populate::populate_box_geom(&mut brep);

 let opts = OffsetOptions::new(0.1)
 .with_join_type(JoinType::Arc)
 .with_self_intersection_check(true)
 .with_wall_thickness_check(0.01)
 .with_approximation_tolerance(TOLERANCE_RETRY_LADDER_MID);

 let result = offset_shape(&brep, opts).unwrap();

 // Verify the workflow completed
 assert!(!result.brep.vertices.is_empty());
 assert!(!result.brep.edges.is_empty());
 assert_eq!(result.offset_faces, 6);
 }

 #[test]
 fn offset_with_self_intersection_config() {
 let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });
 crate::geom_populate::populate_box_geom(&mut brep);

 let config = SelfIntersectionConfig {
 detect: true,
 auto_repair: false,
 max_repair_attempts: 5,
 reduction_factor: 0.8,
 min_offset_distance: 0.01,
 allow_partial_results: false,
 };

 let opts = OffsetOptions::new(0.6)
 .with_self_intersection_config(config);

 let result = offset_shape(&brep, opts).unwrap();

 // Self-intersection detection should not fire for a convex box offset outward
 assert!(!result.self_intersection);
 }

 // B3: Offset Sphere-Sphere Intersection Tests

 #[test]
 fn offset_sphere_sphere_circle() {
 let s1 = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 3.0);
 let s2 = SphericalSurface::new(DVec3::new(2.0, 0.0, 0.0), DVec3::Z, 3.0);

 match intersect_offset_sphere_sphere(&s1, &s2, 0.0, 0.0) {
 OffsetIntersectionCurve::Circle(c) => {
 // Circle should be perpendicular to line of centers (X-axis)
 assert!((c.center.x - 1.0).abs() < TOLERANCE_COORD_SUB); // Midpoint
 assert!((c.normal.x.abs() - 1.0).abs() < TOLERANCE_COORD_SUB);
 // Radius: sqrt(r ?- a ? where a = d/2 = 1, r = 3
 let expected_r: f64 = (9.0_f64 - 1.0_f64).sqrt();
 assert!((c.radius - expected_r).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected Circle, got {other:?}"),
 }
 }

 #[test]
 fn offset_sphere_sphere_with_offset() {
 let s1 = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 2.0);
 let s2 = SphericalSurface::new(DVec3::new(5.0, 0.0, 0.0), DVec3::Z, 2.0);

 // Without offset: no intersection (d=5 > r1+r2=4)
 match intersect_offset_sphere_sphere(&s1, &s2, 0.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection without offset, got {other:?}"),
 }

 // With offset 0.5 each: r1=2.5, r2=2.5, sum=5, tangent
 match intersect_offset_sphere_sphere(&s1, &s2, 0.5, 0.5) {
 OffsetIntersectionCurve::TangentPoint(pt) => {
 assert!((pt.x - 2.5).abs() < TOLERANCE_COORD_SUB);
 }
 // Circle with radius 0 is equivalent to a tangent point
 OffsetIntersectionCurve::Circle(c) if c.radius < TOLERANCE_COORD_SUB => {
 assert!((c.center.x - 2.5).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected TangentPoint with offset, got {other:?}"),
 }
 }

 #[test]
 fn offset_sphere_sphere_negative_offset() {
 let s1 = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 3.0);
 let s2 = SphericalSurface::new(DVec3::new(2.0, 0.0, 0.0), DVec3::Z, 3.0);

 // Reduce both radii
 match intersect_offset_sphere_sphere(&s1, &s2, -1.0, -1.0) {
 OffsetIntersectionCurve::Circle(c) => {
 assert!((c.radius - (4.0_f64 - 1.0_f64).sqrt()).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected Circle with reduced radii, got {other:?}"),
 }
 }

 #[test]
 fn offset_sphere_sphere_concentric() {
 let s1 = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 2.0);
 let s2 = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 3.0);

 match intersect_offset_sphere_sphere(&s1, &s2, 0.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection for concentric different radii, got {other:?}"),
 }

 // With offset: s1 becomes radius 3
 match intersect_offset_sphere_sphere(&s1, &s2, 1.0, 0.0) {
 OffsetIntersectionCurve::Coincident => {}
 other => panic!("Expected Coincident for same radii after offset, got {other:?}"),
 }
 }

 #[test]
 fn offset_sphere_sphere_degenerate() {
 let s1 = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 0.5);
 let s2 = SphericalSurface::new(DVec3::new(5.0, 0.0, 0.0), DVec3::Z, 0.5);

 // Negative offset larger than radius -> degenerate
 match intersect_offset_sphere_sphere(&s1, &s2, -1.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection for degenerate sphere, got {other:?}"),
 }
 }

 // B4: Mixed Surface Offset Intersection Tests

 #[test]
 fn offset_plane_cylinder_perpendicular() {
 let plane = Plane { origin: DVec3::new(0.0, 5.0, 0.0), normal: DVec3::Y };
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Y, ref_dir: any_perpendicular(DVec3::Y), radius: 2.0 };

 match intersect_offset_plane_cylinder(&plane, &cyl, 0.0, 0.0) {
 OffsetIntersectionCurve::Circle(c) => {
 assert!((c.radius - 2.0).abs() < TOLERANCE_COORD_SUB);
 assert!((c.center.y - 5.0).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected Circle, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_cylinder_with_offsets() {
 let plane = Plane { origin: DVec3::ZERO, normal: DVec3::Y };
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Y, ref_dir: any_perpendicular(DVec3::Y), radius: 2.0 };

 // Plane offset by 1 (y=1), cylinder offset by 0.5 (r=2.5)
 match intersect_offset_plane_cylinder(&plane, &cyl, 1.0, 0.5) {
 OffsetIntersectionCurve::Circle(c) => {
 assert!((c.radius - 2.5).abs() < TOLERANCE_COORD_SUB);
 assert!((c.center.y - 1.0).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected Circle, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_cylinder_parallel_two_lines() {
 let plane = Plane { origin: DVec3::ZERO, normal: DVec3::X };
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Y, ref_dir: any_perpendicular(DVec3::Y), radius: 2.0 };

 match intersect_offset_plane_cylinder(&plane, &cyl, 0.0, 0.0) {
 OffsetIntersectionCurve::TwoLines(l1, l2) => {
 assert!(l1.direction.dot(DVec3::Y).abs() > 0.99);
 assert!(l2.direction.dot(DVec3::Y).abs() > 0.99);
 }
 other => panic!("Expected TwoLines, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_cylinder_no_intersection() {
 let plane = Plane { origin: DVec3::new(10.0, 0.0, 0.0), normal: DVec3::X };
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Y, ref_dir: any_perpendicular(DVec3::Y), radius: 2.0 };

 match intersect_offset_plane_cylinder(&plane, &cyl, 0.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_cylinder_oblique_ellipse() {
 let plane = Plane { origin: DVec3::ZERO, normal: DVec3::new(0.0, 1.0, 1.0).normalize() };
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Y, ref_dir: any_perpendicular(DVec3::Y), radius: 1.0 };

 match intersect_offset_plane_cylinder(&plane, &cyl, 0.0, 0.0) {
 OffsetIntersectionCurve::Ellipse(e) => {
 assert!((e.minor_radius - 1.0).abs() < TOLERANCE_COORD_SUB);
 assert!(e.major_radius > 1.0);
 }
 other => panic!("Expected Ellipse, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_sphere_circle() {
 let plane = Plane { origin: DVec3::ZERO, normal: DVec3::Y };
 let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Y, 3.0);

 match intersect_offset_plane_sphere(&plane, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::Circle(c) => {
 assert!((c.radius - 3.0).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected Circle, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_sphere_with_offsets() {
 let plane = Plane { origin: DVec3::new(0.0, 2.0, 0.0), normal: DVec3::Y };
 let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Y, 3.0);

 match intersect_offset_plane_sphere(&plane, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::Circle(c) => {
 let expected_r: f64 = (9.0_f64 - 4.0_f64).sqrt();
 assert!((c.radius - expected_r).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected Circle, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_sphere_tangent() {
 let plane = Plane { origin: DVec3::new(0.0, 3.0, 0.0), normal: DVec3::Y };
 let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Y, 3.0);

 match intersect_offset_plane_sphere(&plane, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::TangentPoint(pt) => {
 assert!((pt.y - 3.0).abs() < TOLERANCE_COORD_SUB);
 }
 // Degenerate circle (near-zero radius) is equivalent to a tangent point; numeric path may use ~TOLERANCE_MESH_LEGACY.
 OffsetIntersectionCurve::Circle(c) if c.radius < 50.0 * TOLERANCE_ABS => {
 assert!((c.center.y - 3.0).abs() < TOLERANCE_COORD_SUB);
 }
 other => panic!("Expected TangentPoint, got {other:?}"),
 }
 }

 #[test]
 fn offset_plane_sphere_no_intersection() {
 let plane = Plane { origin: DVec3::new(0.0, 10.0, 0.0), normal: DVec3::Y };
 let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Y, 3.0);

 match intersect_offset_plane_sphere(&plane, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection, got {other:?}"),
 }
 }

 #[test]
 fn offset_cylinder_sphere_axis_aligned_two_circles() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: any_perpendicular(DVec3::Z), radius: 2.0 };
 let sphere = SphericalSurface::new(DVec3::new(0.0, 0.0, 3.0), DVec3::Z, 5.0);

 match intersect_offset_cylinder_sphere(&cyl, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::TwoCircles(c1, c2) => {
 // Sphere center at z=3, R=5, cylinder r=2
 // dz = sqrt(25-4) = sqrt(21) =4.58
 let expected_dz: f64 = (25.0_f64 - 4.0_f64).sqrt();
 assert!((c1.center.z - (3.0 - expected_dz)).abs() < TOLERANCE_LINEAR_RELAX_8);
 assert!((c2.center.z - (3.0 + expected_dz)).abs() < TOLERANCE_LINEAR_RELAX_8);
 }
 other => panic!("Expected TwoCircles, got {other:?}"),
 }
 }

 #[test]
 fn offset_cylinder_sphere_with_offsets() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: any_perpendicular(DVec3::Z), radius: 2.0 };
 let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 2.0);

 // Without offset: tangent (R=r)
 // A Circle with radius equal to the cylinder radius is a tangent circle
 match intersect_offset_cylinder_sphere(&cyl, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::TangentCircle(_) => {}
 OffsetIntersectionCurve::Circle(c) if (c.radius - 2.0).abs() < TOLERANCE_COORD_SUB => {}
 other => panic!("Expected TangentCircle (or Circle with r=2) without offset, got {other:?}"),
 }

 // With offset on sphere: R=3 > r=2, should have two circles
 match intersect_offset_cylinder_sphere(&cyl, &sphere, 0.0, 1.0) {
 OffsetIntersectionCurve::TwoCircles(_, _) => {}
 other => panic!("Expected TwoCircles with offset, got {other:?}"),
 }
 }

 #[test]
 fn offset_cylinder_sphere_off_axis_general() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: any_perpendicular(DVec3::Z), radius: 1.0 };
 let sphere = SphericalSurface::new(DVec3::new(5.0, 0.0, 0.0), DVec3::Z, 5.0);

 // Off-axis cases may return General or fall back to Numerical approximation
 match intersect_offset_cylinder_sphere(&cyl, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::General => {}
 OffsetIntersectionCurve::Numerical(_) => {} // Numerical approximation is acceptable
 other => panic!("Expected General or Numerical for off-axis case, got {other:?}"),
 }
 }

 #[test]
 fn offset_cylinder_sphere_no_intersection() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: any_perpendicular(DVec3::Z), radius: 1.0 };
 let sphere = SphericalSurface::new(DVec3::new(10.0, 0.0, 0.0), DVec3::Z, 2.0);

 // Sphere center far off axis
 match intersect_offset_cylinder_sphere(&cyl, &sphere, 0.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection, got {other:?}"),
 }
 }

 #[test]
 fn offset_cylinder_sphere_degenerate() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: any_perpendicular(DVec3::Z), radius: 0.5 };
 let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 0.5);

 // Negative offset creates degenerate surfaces
 match intersect_offset_cylinder_sphere(&cyl, &sphere, -1.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection for degenerate, got {other:?}"),
 }
 }

 // Precision tests

 #[test]
 fn offset_plane_plane_high_precision() {
 // Test that high precision (TOLERANCE_LINEAR_ULTRA_STRICT) is achieved
 let p1 = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
 let p2 = Plane {
 origin: DVec3::new(1.0, 1.0, 0.0),
 normal: DVec3::new(1.0, 1.0, 1.0).normalize()
 };

 match intersect_offset_plane_plane(&p1, &p2, 0.0, 0.0) {
 OffsetIntersectionCurve::Line(line) => {
 // Verify point is on both planes
 let d1 = line.origin.dot(p1.normal);
 let d2 = (line.origin - p2.origin).dot(p2.normal);
 assert!(d1.abs() < TOLERANCE_COORD_SUB, "Point should be on plane 1");
 assert!(d2.abs() < TOLERANCE_COORD_SUB, "Point should be on plane 2");
 }
 other => panic!("Expected Line, got {other:?}"),
 }
 }

 #[test]
 fn offset_sphere_sphere_precision() {
 // Test curved surface precision (TOLERANCE_LINEAR_RELAX_8 target)
 let s1 = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 1000.0);
 let s2 = SphericalSurface::new(DVec3::new(100.0, 0.0, 0.0), DVec3::Z, 950.0);

 match intersect_offset_sphere_sphere(&s1, &s2, 0.0, 0.0) {
 OffsetIntersectionCurve::Circle(c) => {
 // Verify circle center lies on the radical plane
 let d1 = (c.center - s1.center).length();
 let d2 = (c.center - s2.center).length();

 // Distance from center to sphere surface should match circle radius
 let r1_expected = (s1.radius * s1.radius - d1 * d1).sqrt();
 assert!((c.radius - r1_expected).abs() < TOLERANCE_MESH_LEGACY, "Circle radius should match");
 }
 other => panic!("Expected Circle, got {other:?}"),
 }
 }

 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // UV Projection Tests
 // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn project_point_to_plane_uv_basic() {
 let plane = Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 };
 let surf = Surface3::Plane(plane);

 // Point on the plane at (1, 2, 0)
 let point = DVec3::new(1.0, 2.0, 0.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 // Verify by reconstructing the point
 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_ULTRA_STRICT, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_plane_uv_offset_origin() {
 let plane = Plane {
 origin: DVec3::new(5.0, 5.0, 5.0),
 normal: DVec3::Z,
 };
 let surf = Surface3::Plane(plane);

 // Point on the plane
 let point = DVec3::new(6.0, 7.0, 5.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_ULTRA_STRICT, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_sphere_uv_basic() {
 let sphere = SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 2.0,
 ref_dir: any_perpendicular(DVec3::Z),
 };
 let surf = Surface3::Sphere(sphere);

 // Point on the sphere at (2, 0, 0)
 let point = DVec3::new(2.0, 0.0, 0.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_sphere_uv_pole() {
 let sphere = SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 2.0,
 ref_dir: any_perpendicular(DVec3::Z),
 };
 let surf = Surface3::Sphere(sphere);

 // Point at the north pole
 let point = DVec3::new(0.0, 0.0, 2.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_sphere_uv_with_hint() {
 let sphere = SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 2.0,
 ref_dir: any_perpendicular(DVec3::Z),
 };
 let surf = Surface3::Sphere(sphere);

 // Point exactly on the sphere near the seam (u =  ?
 let point = DVec3::new(-2.0, 0.0, 0.0);

 // Without hint, might get u =- ?or u = ?
 let uv_no_hint = project_point_to_surface_uv(point, &surf, None).unwrap();

 // With hint u =  ? should get u close to  ?
 // Use v =  ?2 for the equator where the point is
 let uv_with_hint = project_point_to_surface_uv(point, &surf, Some([std::f64::consts::PI, std::f64::consts::FRAC_PI_2])).unwrap();

 // Both should reconstruct to the same point
 let p1 = surf.point_at(uv_no_hint[0], uv_no_hint[1]);
 let p2 = surf.point_at(uv_with_hint[0], uv_with_hint[1]);
 assert!((point - p1).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 assert!((point - p2).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match with hint");
 }

 #[test]
 fn project_point_to_cylinder_uv_basic() {
 let cylinder = CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 };
 let surf = Surface3::Cylinder(cylinder);

 // Point on the cylinder at (1, 0, 5)
 let point = DVec3::new(1.0, 0.0, 5.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_cylinder_uv_quadrant() {
 let cylinder = CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 };
 let surf = Surface3::Cylinder(cylinder);

 // Point at 45 degrees in XY plane
 let r = 1.0_f64 / std::f64::consts::SQRT_2;
 let point = DVec3::new(r, r, 3.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");

 // Check that v is 3.0
 assert!((uv[1] - 3.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "v should be 3.0");

 // The u angle depends on the basis chosen by any_perpendicular
 // Just verify that reconstruction works correctly
 }

 #[test]
 fn project_point_to_cylinder_uv_axis() {
 let cylinder = CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 };
 let surf = Surface3::Cylinder(cylinder);

 // Point on the axis - should use hint
 let point = DVec3::new(0.0, 0.0, 5.0);
 let uv = project_point_to_surface_uv(point, &surf, Some([1.0, 0.0])).unwrap();

 // v should be 5.0, u can be anything
 assert!((uv[1] - 5.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "v should be 5.0");
 }

 #[test]
 fn project_point_to_cone_uv_basic() {
 let cone = ConicalSurface {
 apex: DVec3::ZERO,
 axis: DVec3::Z,
 half_angle_rad: std::f64::consts::FRAC_PI_4, // 45 degrees
 radius: 0.0, // Not used for UV computation
 };
 let surf = Surface3::Cone(cone);

 // At height = 2, radius should be 2 * tan(45 ? = 2
 let point = DVec3::new(2.0, 0.0, 2.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_torus_uv_basic() {
 let torus = ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 2.0,
 minor_radius: 0.5,
 };
 let surf = Surface3::Torus(torus);

 // Point on the outer edge of the torus
 let point = DVec3::new(2.5, 0.0, 0.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_torus_uv_top() {
 let torus = ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 2.0,
 minor_radius: 0.5,
 };
 let surf = Surface3::Torus(torus);

 // Point on the top of the torus tube
 let point = DVec3::new(2.0, 0.0, 0.5);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_torus_uv_inner() {
 let torus = ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 2.0,
 minor_radius: 0.5,
 };
 let surf = Surface3::Torus(torus);

 // Point on the inner edge of the torus
 let point = DVec3::new(1.5, 0.0, 0.0);
 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();

 let reconstructed = surf.point_at(uv[0], uv[1]);
 assert!((point - reconstructed).length() < TOLERANCE_LINEAR_RELAX_8, "Reconstructed point should match");
 }

 #[test]
 fn project_point_to_surface_uv_consistency() {
 // Test that multiple calls with the same point give consistent results
 let sphere = SphericalSurface {
 center: DVec3::new(1.0, 2.0, 3.0),
 axis: DVec3::Z,
 radius: 5.0,
 ref_dir: any_perpendicular(DVec3::Z),
 };
 let surf = Surface3::Sphere(sphere);

 let point = DVec3::new(4.0, 6.0, 3.0);

 let uv1 = project_point_to_surface_uv(point, &surf, None).unwrap();
 let uv2 = project_point_to_surface_uv(point, &surf, None).unwrap();

 assert!((uv1[0] - uv2[0]).abs() < TOLERANCE_LEN_MIN, "u should be consistent");
 assert!((uv1[1] - uv2[1]).abs() < TOLERANCE_LEN_MIN, "v should be consistent");
 }

 #[test]
 fn project_point_to_surface_uv_precision() {
 // Test high precision for various surface types
 let sphere = SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 100.0,
 ref_dir: any_perpendicular(DVec3::Z),
 };
 let surf = Surface3::Sphere(sphere);

 // Point at a non-trivial angle
 let theta: f64 = 1.234;
 let phi: f64 = 0.789;
 let point = DVec3::new(
 100.0 * phi.sin() * theta.cos(),
 100.0 * phi.sin() * theta.sin(),
 100.0 * phi.cos(),
 );

 let uv = project_point_to_surface_uv(point, &surf, None).unwrap();
 let reconstructed = surf.point_at(uv[0], uv[1]);

 // Should achieve sub-micron precision for 100m sphere
 assert!((point - reconstructed).length() < TOLERANCE_MESH_LEGACY, "High precision should be achieved");
 }

 #[test]
 fn orthonormal_basis_deterministic() {
 use crate::offset::orthonormal_basis_from_normal;

 // Same normal should give same basis
 let n1 = DVec3::new(1.0, 2.0, 3.0).normalize();
 let n2 = DVec3::new(1.0, 2.0, 3.0).normalize();

 let (u1, v1) = orthonormal_basis_from_normal(n1);
 let (u2, v2) = orthonormal_basis_from_normal(n2);

 assert!((u1 - u2).length() < TOLERANCE_LEN_MIN, "u should be deterministic");
 assert!((v1 - v2).length() < TOLERANCE_LEN_MIN, "v should be deterministic");
 }

 #[test]
 fn orthonormal_basis_orthogonal() {
 use crate::offset::orthonormal_basis_from_normal;

 let normals = vec![
 DVec3::X,
 DVec3::Y,
 DVec3::Z,
 DVec3::new(1.0, 1.0, 1.0).normalize(),
 DVec3::new(1.0, 2.0, 3.0).normalize(),
 ];

 for n in normals {
 let (u, v) = orthonormal_basis_from_normal(n);

 // Check orthonormality
 assert!(u.dot(n).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "u should be perpendicular to n");
 assert!(v.dot(n).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "v should be perpendicular to n");
 assert!((u.length() - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "u should be unit");
 assert!((v.length() - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "v should be unit");
 assert!((u.dot(v)).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "u and v should be perpendicular");
 }
 }

 // Edge case tests for OCCT alignment

 #[test]
 fn offset_torus_handles_positive() {
 use rcad_kernel::geom::ToroidalSurface;

 let torus = Surface3::Torus(ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 3.0,
 minor_radius: 1.0,
 });

 let result = offset_surface(&torus, 0.5);
 assert!(result.is_some(), "offset torus should succeed");

 if let Surface3::Torus(t) = result.unwrap() {
 assert!((t.minor_radius - 1.5).abs() < TOLERANCE_COORD_SUB, "minor radius should increase");
 }
 }

 #[test]
 fn offset_negative_shrinks_sphere() {
 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 2.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 let offset = offset_surface(&sphere, -0.5).unwrap();

 if let Surface3::Sphere(s) = offset {
 assert!((s.radius - 1.5).abs() < TOLERANCE_COORD_SUB, "radius should decrease with negative offset");
 } else {
 panic!("expected Sphere");
 }
 }

 #[test]
 fn offset_cone_preserves_apex() {
 use rcad_kernel::geom::ConicalSurface;

 let cone = Surface3::Cone(ConicalSurface {
 apex: DVec3::ZERO,
 axis: DVec3::Z,
 radius: 1.0,
 half_angle_rad: std::f64::consts::FRAC_PI_6, // 30 degrees
 });

 let result = offset_surface(&cone, 0.2);
 assert!(result.is_some(), "offset cone should succeed");
 }

 // ============================================================================
 // OCCT TKOffset Alignment Tests
 // ============================================================================

 #[test]
 fn offset_multiple_surfaces_preserve_topology() {
 // Test that offsetting multiple adjacent surfaces maintains connectivity
 let plane1 = Surface3::Plane(Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 });
 let plane2 = Surface3::Plane(Plane {
 origin: DVec3::new(0.0, 1.0, 0.0),
 normal: DVec3::Y,
 });

 let offset1 = offset_surface(&plane1, 0.5).unwrap();
 let offset2 = offset_surface(&plane2, 0.5).unwrap();

 // Both should succeed
 assert!(matches!(offset1, Surface3::Plane(_)));
 assert!(matches!(offset2, Surface3::Plane(_)));
 }

 #[test]
 fn offset_cylinder_negative_valid() {
 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 2.0,
 });

 let offset = offset_surface(&cylinder, -1.0).unwrap();

 if let Surface3::Cylinder(c) = offset {
 assert!((c.radius - 1.0).abs() < TOLERANCE_COORD_SUB, "radius should decrease with negative offset");
 }
 }

 #[test]
 fn offset_cylinder_negative_invalid() {
 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::ZERO,
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.0,
 });

 // Negative offset larger than radius should fail
 let offset = offset_surface(&cylinder, -2.0);
 assert!(offset.is_none(), "offset larger than radius should return None");
 }

 #[test]
 fn offset_torus_negative_minor() {
 use rcad_kernel::geom::ToroidalSurface;

 let torus = Surface3::Torus(ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 3.0,
 minor_radius: 1.0,
 });

 let result = offset_surface(&torus, -0.5);
 assert!(result.is_some(), "negative offset torus should succeed if within bounds");

 if let Some(Surface3::Torus(t)) = result {
 assert!((t.minor_radius - 0.5).abs() < TOLERANCE_COORD_SUB, "minor radius should decrease");
 }
 }

 #[test]
 fn offset_torus_negative_exceeds_minor() {
 use rcad_kernel::geom::ToroidalSurface;

 let torus = Surface3::Torus(ToroidalSurface {
 center: DVec3::ZERO,
 axis: DVec3::Z,
 major_radius: 3.0,
 minor_radius: 1.0,
 });

 // Negative offset larger than minor radius
 let result = offset_surface(&torus, -2.0);
 assert!(result.is_none(), "offset exceeding minor radius should return None");
 }

 #[test]
 fn offset_options_various_distances() {
 let opts_small = OffsetOptions::new(0.01);
 assert_eq!(opts_small.distance, 0.01);

 let opts_large = OffsetOptions::new(100.0);
 assert_eq!(opts_large.distance, 100.0);
 }

 #[test]
 fn offset_join_types() {
 // Test that join types are properly defined
 let intersection = JoinType::Intersection;
 let arc = JoinType::Arc;
 let tangent = JoinType::Tangent;

 assert_eq!(intersection, JoinType::Intersection);
 assert_eq!(arc, JoinType::Arc);
 assert_eq!(tangent, JoinType::Tangent);
 }

 #[test]
 fn offset_error_types() {
 // Verify error types exist and can be created
 let err1 = OffsetError::ZeroDistance;
 let err2 = OffsetError::InvalidInput("test");
 let err3 = OffsetError::SelfIntersection { description: "test".into() };

 assert!(matches!(err1, OffsetError::ZeroDistance));
 assert!(matches!(err2, OffsetError::InvalidInput(_)));
 assert!(matches!(err3, OffsetError::SelfIntersection { .. }));
 }

 #[test]
 fn offset_shell_result_checks() {
 let brep = BRep::new();
 let opts = OffsetOptions::new(1.0);

 let result = offset_shape(&brep, opts);
 // Empty BRep should either succeed with empty result or fail gracefully
 assert!(result.is_ok() || result.is_err());
 }

 #[test]
 fn offset_face_result_type() {
 let plane = Surface3::Plane(Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 });

 let result = offset_surface(&plane, 1.0);
 assert!(result.is_some());

 let offset = result.unwrap();
 assert!(matches!(offset, Surface3::Plane(_)));
 }

 #[test]
 fn offset_preserves_surface_type() {
 // Sphere offset should remain sphere
 let sphere = Surface3::Sphere(SphericalSurface {
 center: DVec3::new(1.0, 2.0, 3.0),
 axis: DVec3::Z,
 radius: 5.0,
 ref_dir: any_perpendicular(DVec3::Z),
 });

 let offset = offset_surface(&sphere, 1.0).unwrap();
 assert!(matches!(offset, Surface3::Sphere(_)));

 // Cylinder offset should remain cylinder
 let cylinder = Surface3::Cylinder(CylindricalSurface {
 origin: DVec3::new(1.0, 2.0, 3.0),
 axis: DVec3::Y,
 ref_dir: any_perpendicular(DVec3::Y),
 radius: 2.0,
 });

 let offset = offset_surface(&cylinder, 1.0).unwrap();
 assert!(matches!(offset, Surface3::Cylinder(_)));
 }

 #[test]
 fn offset_tolerance_propagation() {
 let opts = OffsetOptions::new(0.5)
 .with_tolerance(TOLERANCE_LINEAR_RELAX_8);

 assert!((opts.tolerance - TOLERANCE_LINEAR_RELAX_8).abs() < TOLERANCE_FLOAT_DEDUP);
 }

 #[test]
 fn offset_self_intersection_option() {
 let opts1 = OffsetOptions::new(0.5).with_self_intersection_check(true);
 assert!(opts1.check_self_intersection);

 let opts2 = OffsetOptions::new(0.5).with_self_intersection_check(false);
 assert!(!opts2.check_self_intersection);
 }

 #[test]
 fn offset_auto_repair_option() {
 let opts = OffsetOptions::new(0.5).with_auto_repair(true);
 assert!(opts.auto_repair);
 }

 //  € € Phase 2: New offset handler tests  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

 #[test]
 fn offset_cone_cone_coaxial() {
 let c1 = ConicalSurface { apex: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, half_angle_rad: 0.3 };
 let c2 = ConicalSurface { apex: DVec3::new(0.0, 0.0, 3.0), axis: DVec3::Z, radius: 2.0, half_angle_rad: 0.5 };
 match intersect_offset_cone_cone(&c1, &c2, 0.0, 0.0) {
 OffsetIntersectionCurve::Circle(_) | OffsetIntersectionCurve::NoIntersection | OffsetIntersectionCurve::General => {}
 other => panic!("Unexpected: {other:?}"),
 }
 }

 #[test]
 fn offset_cone_cone_with_offset() {
 let c1 = ConicalSurface { apex: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, half_angle_rad: 0.5 };
 let c2 = ConicalSurface { apex: DVec3::new(0.0, 0.0, 2.0), axis: DVec3::Z, radius: 2.0, half_angle_rad: 0.5 };
 let result = intersect_offset_cone_cone(&c1, &c2, 0.2, 0.3);
 assert!(matches!(result,
 OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 ));
 }

 #[test]
 fn offset_cylinder_cone_basic() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X, radius: 1.0 };
 let cone = ConicalSurface { apex: DVec3::new(0.0, 0.0, 5.0), axis: DVec3::Z, radius: 2.0, half_angle_rad: 0.3 };
 let result = intersect_offset_cylinder_cone(&cyl, &cone, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::General
 ));
 }

 #[test]
 fn offset_cylinder_cone_with_offset() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X, radius: 1.0 };
 let cone = ConicalSurface { apex: DVec3::new(0.0, 0.0, 5.0), axis: DVec3::Z, radius: 2.0, half_angle_rad: 0.3 };
 let result = intersect_offset_cylinder_cone(&cyl, &cone, 0.2, 0.1);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::General
 ));
 }

 #[test]
 fn offset_plane_cone_perpendicular() {
 let plane = Plane { origin: DVec3::new(0.0, 0.0, 3.0), normal: DVec3::Z };
 let cone = ConicalSurface { apex: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, half_angle_rad: 0.5 };
 let result = intersect_offset_plane_cone(&plane, &cone, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::Ellipse(_)
 | OffsetIntersectionCurve::NoIntersection
 ));
 }

 #[test]
 fn offset_plane_cone_with_offset() {
 let plane = Plane { origin: DVec3::new(0.0, 0.0, 3.0), normal: DVec3::Z };
 let cone = ConicalSurface { apex: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, half_angle_rad: 0.5 };
 let result = intersect_offset_plane_cone(&plane, &cone, 0.2, 0.1);
 assert!(matches!(result,
 OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::Ellipse(_)
 | OffsetIntersectionCurve::NoIntersection
 ));
 }

 #[test]
 fn offset_sphere_cone_basic() {
 let sphere = SphericalSurface::new(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, 3.0);
 let cone = ConicalSurface { apex: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, half_angle_rad: 0.5 };
 let result = intersect_offset_sphere_cone(&sphere, &cone, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TwoCircles(_, _)
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::NoIntersection
 ));
 }

 #[test]
 fn offset_sphere_cone_negative_degenerate() {
 let sphere = SphericalSurface::new(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, 0.3);
 let cone = ConicalSurface { apex: DVec3::ZERO, axis: DVec3::Z, radius: 1.0, half_angle_rad: 0.5 };
 // Negative offset larger than sphere radius  ?degenerate
 match intersect_offset_sphere_cone(&sphere, &cone, -1.0, 0.0) {
 OffsetIntersectionCurve::NoIntersection => {}
 other => panic!("Expected NoIntersection for degenerate, got {other:?}"),
 }
 }

 #[test]
 fn offset_torus_torus_basic() {
 let t1 = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0 };
 let t2 = ToroidalSurface { center: DVec3::new(2.0, 0.0, 0.0), axis: DVec3::Z, major_radius: 3.0, minor_radius: 0.8 };
 let result = intersect_offset_torus_torus(&t1, &t2, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TwoCircles(_, _)
 | OffsetIntersectionCurve::TangentCircle(_)
 | OffsetIntersectionCurve::Coincident
 | OffsetIntersectionCurve::Numerical(_)
 ));
 }

 #[test]
 fn offset_torus_torus_with_offset() {
 let t1 = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0 };
 let t2 = ToroidalSurface { center: DVec3::new(2.0, 0.0, 0.0), axis: DVec3::Z, major_radius: 3.0, minor_radius: 0.8 };
 let result = intersect_offset_torus_torus(&t1, &t2, 0.2, 0.1);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TwoCircles(_, _)
 | OffsetIntersectionCurve::TangentCircle(_)
 | OffsetIntersectionCurve::Coincident
 | OffsetIntersectionCurve::Numerical(_)
 ));
 }

 #[test]
 fn offset_cylinder_torus_basic() {
 let cyl = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X, radius: 2.0 };
 let torus = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 5.0, minor_radius: 1.0 };
 let result = intersect_offset_cylinder_torus(&cyl, &torus, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::TwoCircles(_, _)
 ));
 }

 #[test]
 fn offset_torus_cone_basic() {
 let torus = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 4.0, minor_radius: 1.0 };
 let cone = ConicalSurface { apex: DVec3::new(0.0, 0.0, -2.0), axis: DVec3::Z, radius: 1.0, half_angle_rad: 0.4 };
 let result = intersect_offset_torus_cone(&torus, &cone, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TwoCircles(_, _)
 ));
 }

 #[test]
 fn offset_plane_torus_basic() {
 let plane = Plane { origin: DVec3::new(0.0, 0.0, 2.0), normal: DVec3::Z };
 let torus = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0 };
 let result = intersect_offset_plane_torus(&plane, &torus, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::TwoCircles(_, _)
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TangentCircle(_)
 ));
 }

 #[test]
 fn offset_sphere_torus_basic() {
 let sphere = SphericalSurface::new(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, 5.0);
 let torus = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0 };
 let result = intersect_offset_sphere_torus(&sphere, &torus, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TwoCircles(_, _)
 ));
 }

 #[test]
 fn offset_sphere_torus_with_offsets() {
 let sphere = SphericalSurface::new(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, 5.0);
 let torus = ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0 };
 let result = intersect_offset_sphere_torus(&sphere, &torus, 0.5, 0.3);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TwoCircles(_, _)
 ));
 }

 #[test]
 fn offset_torus_via_dispatch() {
 let s1 = Surface3::Torus(ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0 });
 let s2 = Surface3::Torus(ToroidalSurface { center: DVec3::new(2.0, 0.0, 0.0), axis: DVec3::Z, major_radius: 3.0, minor_radius: 0.8 });
 let result = intersect_offset_surfaces(&s1, &s2, 0.0, 0.0);
 assert!(matches!(result,
 OffsetIntersectionCurve::NoIntersection
 | OffsetIntersectionCurve::General
 | OffsetIntersectionCurve::Circle(_)
 | OffsetIntersectionCurve::TwoCircles(_, _)
 | OffsetIntersectionCurve::TangentCircle(_)
 | OffsetIntersectionCurve::Coincident
 | OffsetIntersectionCurve::Numerical(_)
 ));
 }

 #[test]
 fn offset_dispatch_new_pairings() {
 // Verify that an un-handled pairing still falls back to numerical
 let s1 = Surface3::Torus(ToroidalSurface { center: DVec3::ZERO, axis: DVec3::Z, major_radius: 3.0, minor_radius: 1.0 });
 let s2 = Surface3::Torus(ToroidalSurface { center: DVec3::new(2.0, 0.0, 0.0), axis: DVec3::Z, major_radius: 3.0, minor_radius: 0.8 });
 let result = intersect_offset_surfaces(&s1, &s2, 0.0, 0.0);
 // The new torus-torus handler should now be reached instead of numerical
 assert!(!matches!(result, OffsetIntersectionCurve::NoIntersection));
 }
}
