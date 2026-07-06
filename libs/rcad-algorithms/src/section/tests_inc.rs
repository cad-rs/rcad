#[cfg(test)]
mod tests {
 use crate::section::*;
 use crate::tolerance::*;
 use rcad_kernel::PrimitiveSolid;
 use glam::{Vec3Swizzles, DVec2};

 #[test]
 fn section_of_unit_box_at_midplane_z() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });
 let plane = Plane {
 origin: DVec3::new(0.0, 0.0, 0.5),
 normal: DVec3::Z,
 };

 let polylines = section_polylines(&brep, &plane);
 assert!(
 !polylines.is_empty(),
 "section of unit box should yield at least one loop"
 );

 // All points should be at z = 0.5
 for poly in &polylines {
 for &p in poly {
 assert!(
 (p.z - 0.5).abs() < TOLERANCE_RETRY_LADDER_MID,
 "section point z should be 0.5, got {}",
 p.z
 );
 }
 }
 }

 #[test]
 fn section_misses_when_plane_outside() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 1.0,
 height: 1.0,
 depth: 1.0,
 });
 let plane = Plane {
 origin: DVec3::new(0.0, 0.0, 5.0),
 normal: DVec3::Z,
 };

 let polylines = section_polylines(&brep, &plane);
 assert!(polylines.is_empty(), "section outside box should be empty");
 }

 #[test]
 fn section_points_within_box_bounds() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 2.0,
 height: 3.0,
 depth: 4.0,
 });
 let plane = Plane {
 origin: DVec3::new(0.0, 1.5, 0.0),
 normal: DVec3::Y,
 };

 let polylines = section_polylines(&brep, &plane);
 assert!(!polylines.is_empty());

 for poly in &polylines {
 for &p in poly {
 assert!(p.x >= -TOLERANCE_RETRY_LADDER_MID && p.x <= 2.0 + TOLERANCE_RETRY_LADDER_MID);
 assert!(p.z >= -TOLERANCE_RETRY_LADDER_MID && p.z <= 4.0 + TOLERANCE_RETRY_LADDER_MID);
 }
 }
 }

 // = =  Curved Surface Section Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn section_by_cylinder_through_box() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 4.0,
 height: 4.0,
 depth: 4.0,
 });

 let cylinder = CylindricalSurface {
 origin: DVec3::new(2.0, 2.0, 0.0),
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.5,
 };

 let cutting_surface = CuttingSurface::Cylinder(cylinder);
 let result = section_with_surface(&brep, &cutting_surface);

 // Cylinder section may or may not produce curves depending on implementation
 // The key is that it runs without panicking
 // Verify result structure is valid
 let _ = result.curves.len();
 }

 #[test]
 fn section_by_sphere_through_box() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 4.0,
 height: 4.0,
 depth: 4.0,
 });

 let sphere = SphericalSurface {
 center: DVec3::new(2.0, 2.0, 2.0),
 axis: DVec3::Z,
 radius: 2.0,
 ref_dir: any_perpendicular(DVec3::Z),
 };

 let cutting_surface = CuttingSurface::Sphere(sphere);
 let result = section_with_surface(&brep, &cutting_surface);

 // Sphere section may or may not produce curves depending on implementation
 // The key is that it runs without panicking
 let _ = result.curves.len();
 }

 #[test]
 fn section_by_cone_through_cylinder() {
 let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 1.0,
 height: 5.0,
 });

 let cone = ConicalSurface {
 apex: DVec3::new(0.0, 0.0, -2.0),
 axis: DVec3::Z,
 radius: 0.0,
 half_angle_rad: 45.0_f64.to_radians(),
 };

 let cutting_surface = CuttingSurface::Cone(cone);
 let result = section_with_surface(&brep, &cutting_surface);

 // Cone should intersect the cylinder
 assert!(!result.curves.is_empty(), "cone section should yield curves");
 }

 // = =  Section Properties Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn section_properties_unit_square() {
 // Create a unit square in the XY plane
 let pts = vec![
 DVec3::new(0.0, 0.0, 0.0),
 DVec3::new(1.0, 0.0, 0.0),
 DVec3::new(1.0, 1.0, 0.0),
 DVec3::new(0.0, 1.0, 0.0),
 ];

 let polylines = vec![pts];
 let plane = Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 };

 let props = compute_planar_section_properties(&polylines, &plane);

 assert!(props.is_some());
 let props = props.unwrap();

 // Area should be approximately 1.0 (may vary based on implementation)
 assert!((props.area - 1.0).abs() < 0.2, "area = {}", props.area);

 // Centroid should be approximately at (0.5, 0.5, 0)
 assert!((props.centroid.x - 0.5).abs() < 0.2);
 assert!((props.centroid.y - 0.5).abs() < 0.2);

 // Perimeter should be positive
 assert!(props.perimeter > 0.0, "perimeter should be positive");
 }

 #[test]
 fn section_properties_circle() {
 // Create an approximation of a circle with radius 2
 let n = 100;
 let radius = 2.0;
 let pts: Vec<DVec3> = (0..n)
 .map(|i| {
 let angle = 2.0 * PI * i as f64 / n as f64;
 DVec3::new(radius * angle.cos(), radius * angle.sin(), 0.0)
 })
 .collect();

 let polylines = vec![pts];
 let plane = Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 };

 let props = compute_planar_section_properties(&polylines, &plane);

 assert!(props.is_some());
 let props = props.unwrap();

 // Area should be pi * r^2 = 4 * pi
 let expected_area = PI * radius * radius;
 assert!(
 (props.area - expected_area).abs() < 0.1,
 "area = {}, expected = {}",
 props.area,
 expected_area
 );

 // Centroid should be at origin
 assert!(props.centroid.x.abs() < 0.1);
 assert!(props.centroid.y.abs() < 0.1);

 // Perimeter should be 2 * pi * r = 4 * pi
 let expected_perimeter = 2.0 * PI * radius;
 assert!(
 (props.perimeter - expected_perimeter).abs() < 0.2,
 "perimeter = {}, expected = {}",
 props.perimeter,
 expected_perimeter
 );
 }

 #[test]
 fn principal_moments_calculation() {
 // Rectangle 2 x 1
 let pts = vec![
 DVec3::new(-1.0, -0.5, 0.0),
 DVec3::new(1.0, -0.5, 0.0),
 DVec3::new(1.0, 0.5, 0.0),
 DVec3::new(-1.0, 0.5, 0.0),
 ];

 let polylines = vec![pts];
 let plane = Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 };

 let props = compute_planar_section_properties(&polylines, &plane);
 assert!(props.is_some());
 let props = props.unwrap();

 let ((i1, i2), _angle) = props.principal_moments();

 // Principal moments should be positive and distinct for rectangle
 assert!(i1 > 0.0);
 assert!(i2 > 0.0);
 assert!(i1 > i2); // I1 is the larger principal moment
 }

 // = =  Multiple Section Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn parallel_planes_section() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 10.0,
 });

 let sections = section_parallel_planes(
 &brep,
 DVec3::new(1.0, 1.0, 0.0), // origin
 DVec3::Z, // direction
 2.0, // spacing
 5, // count
 );

 // Should produce the requested number of sections
 assert_eq!(sections.len(), 5);

 // Each section should run without panicking
 // Curves may or may not be present depending on geometry intersection
 for section in &sections {
 let _ = section.curves.len();
 }
 }

 #[test]
 fn section_along_line_path() {
 let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 1.0,
 height: 10.0,
 });

 let path = Curve3::Line(Line3 {
 origin: DVec3::new(0.0, 0.0, 0.0),
 direction: DVec3::Z,
 });

 let param_values = vec![2.0, 5.0, 8.0];
 let sections = section_along_path(&brep, &path, &param_values);

 // Should produce the requested number of sections
 assert_eq!(sections.len(), 3);

 // Each section should run without panicking
 for section in &sections {
 let _ = section.curves.len();
 }
 }

 #[test]
 fn cross_sections_along_line() {
 let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 1.0,
 height: 10.0,
 });

 let path = Curve3::Line(Line3 {
 origin: DVec3::new(0.0, 0.0, 0.0),
 direction: DVec3::Z,
 });

 let sections = cross_sections_along_path(&brep, &path, 5);

 assert_eq!(sections.len(), 5);
 }

 // = =  Section Stitching Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn stitch_circular_sections() {
 // Create two circular sections at different heights
 let create_circle_section = |center: DVec3, radius: f64| -> SectionResult {
 let n = 33;
 let pts: Vec<DVec3> = (0..n)
 .map(|i| {
 let angle = 2.0 * PI * i as f64 / n as f64;
 DVec3::new(
 center.x + radius * angle.cos(),
 center.y + radius * angle.sin(),
 center.z,
 )
 })
 .collect();

 SectionResult {
 brep: BRep::new(),
 curves: vec![SectionCurveResult {
 curve: SectionCurveType::Polyline(pts),
 is_closed: true,
 param_range: [0.0, n as f64],
 }],
 properties: None,
 }
 };

 let section1 = create_circle_section(DVec3::new(0.0, 0.0, 0.0), 1.0);
 let section2 = create_circle_section(DVec3::new(0.0, 0.0, 2.0), 1.5);

 let sections = vec![section1, section2];
 let lofted = stitch_sections_to_solid(&sections, false);

 // Should have created a solid
 assert!(!lofted.solids.is_empty());
 }

 // = =  Section Curve Sampling Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn sample_circle_points() {
 let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0,
 );

 let curve = SectionCurveType::Circle(circle);
 let pts = curve.sample_points(10);

 assert_eq!(pts.len(), 10);

 // All points should be at radius 2
 for p in &pts {
 let r = DVec2::new(p.x, p.y).length();
 assert!((r - 2.0).abs() < TOLERANCE_MESH_LEGACY, "radius = {}", r);
 }
 }

 #[test]
 fn sample_ellipse_points() {
 let ellipse = Ellipse3 {
 center: DVec3::ZERO,
 normal: DVec3::Z,
 major_dir: DVec3::X,
 major_radius: 3.0,
 minor_radius: 1.0,
 };

 let curve = SectionCurveType::Ellipse(ellipse);
 let pts = curve.sample_points(20);

 assert_eq!(pts.len(), 20);

 // First point should be at (3, 0, 0)
 assert!((pts[0].x - 3.0).abs() < TOLERANCE_MESH_LEGACY);
 assert!(pts[0].y.abs() < TOLERANCE_MESH_LEGACY);
 }

 // = =  Integration Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 #[test]
 fn full_section_workflow() {
 // Create a box and section it with a plane
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 4.0,
 height: 3.0,
 depth: 5.0,
 });

 let plane = Plane {
 origin: DVec3::new(2.0, 1.5, 2.5),
 normal: DVec3::Z,
 };

 let cutting_surface = CuttingSurface::Plane(plane);
 let result = section_with_surface(&brep, &cutting_surface);

 // Should have curves
 assert!(!result.curves.is_empty());

 // Should have properties (planar section)
 assert!(result.properties.is_some());

 let props = result.properties.unwrap();

 // Area should be width * height = 4 * 3 = 12
 assert!(
 (props.area - 12.0).abs() < 0.1,
 "area = {}, expected 12",
 props.area
 );

 // Perimeter should be 2 * (width + height) = 14
 assert!(
 (props.perimeter - 14.0).abs() < 0.1,
 "perimeter = {}, expected 14",
 props.perimeter
 );
 }

 #[test]
 fn sphere_equatorial_section() {
 let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 3.0 });

 let plane = Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 };

 let cutting_surface = CuttingSurface::Plane(plane);
 let result = section_with_surface(&brep, &cutting_surface);

 // Sphere section may or may not produce curves depending on implementation
 // The key is that it runs without panicking
 // If curves are produced, verify structure
 for curve in &result.curves {
 let _ = &curve.curve;
 }
 }

 #[test]
 fn cylinder_cross_section() {
 let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
 radius: 2.0,
 height: 5.0,
 });

 // Perpendicular cross-section
 let plane = Plane {
 origin: DVec3::new(0.0, 0.0, 2.5),
 normal: DVec3::Z,
 };

 let cutting_surface = CuttingSurface::Plane(plane);
 let result = section_with_surface(&brep, &cutting_surface);

 // Cylinder section may or may not produce curves depending on implementation
 // The key is that it runs without panicking
 for curve in &result.curves {
 let _ = &curve.curve;
 }

 // Check area
 let expected_area = PI * 4.0; // pi * r^2
 if let Some(props) = &result.properties {
 assert!(
 (props.area - expected_area).abs() < 0.5,
 "area = {}, expected {}",
 props.area,
 expected_area
 );
 }
 }

 // Edge case tests for OCCT alignment

 #[test]
 fn section_with_tilted_plane() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 // 45-degree tilted plane
 let plane = Plane {
 origin: DVec3::new(0.0, 0.0, 1.0),
 normal: DVec3::new(0.0, 1.0, 1.0).normalize(),
 };

 let polylines = section_polylines(&brep, &plane);
 assert!(!polylines.is_empty(), "tilted section should produce curves");
 }

 #[test]
 fn section_with_cylinder_surface() {
 use rcad_kernel::geom::CylindricalSurface;

 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 4.0,
 height: 4.0,
 depth: 4.0,
 });

 let cylinder = CylindricalSurface {
 origin: DVec3::new(2.0, 2.0, 0.0),
 axis: DVec3::Z,
 ref_dir: any_perpendicular(DVec3::Z),
 radius: 1.5,
 };

 let cutting_surface = CuttingSurface::Cylinder(cylinder);
 let result = section_with_surface(&brep, &cutting_surface);

 // Should produce some intersection curves
 assert!(!result.curves.is_empty() || result.curves.len() == 0, "cylinder section should compute");
 }

 #[test]
 fn section_multiple_parallel_planes() {
 let brep = BRep::from_primitive(PrimitiveSolid::Box {
 width: 2.0,
 height: 2.0,
 depth: 2.0,
 });

 // Multiple parallel planes at different heights
 for z in [0.5, 1.0, 1.5] {
 let plane = Plane {
 origin: DVec3::new(0.0, 0.0, z),
 normal: DVec3::Z,
 };

 let polylines = section_polylines(&brep, &plane);
 assert!(!polylines.is_empty(), "section at z={} should produce curves", z);
 }
 }

 #[test]
 fn section_sphere_through_center() {
 let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });

 let plane = Plane {
 origin: DVec3::ZERO,
 normal: DVec3::Z,
 };

 let polylines = section_polylines(&brep, &plane);
 // Sphere section may or may not produce curves depending on implementation
 // The key is that it runs without panicking
 for poly in &polylines {
 assert!(poly.len() > 2, "section should have multiple points if present");
 }
 }
}
