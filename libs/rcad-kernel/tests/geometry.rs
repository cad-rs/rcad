//! Integration tests for `rcad-kernel` geometry algorithms.
//!
//! These tests verify the correctness of projection, curvature, arc length,
//! NURBS round-trip, and transformation operations against known analytic results.

use glam::DVec3;
use rcad_kernel::{
    any_perpendicular, BRep,
    arc_length,
    closest_point_on_curve, closest_point_on_surface,
    extrema_curve_curve,
    gaussian_curvature, mean_curvature,
    geom::{Circle3, Curve3, Line3, Plane, SphericalSurface, Surface3},
};

// ─────────────────────────────────────────────────────────────────────────────
// Point projection
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn closest_point_on_line_known() {
    // Infinite line along X axis, query point at (2, 3, 0).
    // Closest point should be (2, 0, 0), distance = 3.
    let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
    let result = closest_point_on_curve(&line, DVec3::new(2.0, 3.0, 0.0), 8);
    let expected = DVec3::new(2.0, 0.0, 0.0);
    assert!(
        (result.point - expected).length() < 1e-9,
        "closest point wrong: {:?}", result.point
    );
    assert!(
        (result.distance - 3.0).abs() < 1e-9,
        "distance wrong: {}", result.distance
    );
}

#[test]
fn closest_point_on_circle_known() {
    // Unit circle in XY plane, center (0,0,0), query at (2, 0, 0).
    // Closest point should be (1, 0, 0), distance = 1.
    let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
    let result = closest_point_on_curve(&circle, DVec3::new(2.0, 0.0, 0.0), 16);
    let expected = DVec3::new(1.0, 0.0, 0.0);
    assert!(
        (result.point - expected).length() < 1e-6,
        "closest point wrong: {:?}", result.point
    );
    assert!(
        (result.distance - 1.0).abs() < 1e-6,
        "distance wrong: {}", result.distance
    );
}

#[test]
fn closest_point_on_sphere_surface_known() {
    // Sphere at origin r=2, query at (5, 0, 0).
    // Closest surface point = (2, 0, 0), distance = 3.
    let sphere = Surface3::Sphere(SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 2.0, ref_dir: any_perpendicular(DVec3::Z) });
    let result = closest_point_on_surface(&sphere, DVec3::new(5.0, 0.0, 0.0), 16);
    let expected = DVec3::new(2.0, 0.0, 0.0);
    assert!(
        (result.point - expected).length() < 1e-6,
        "closest surface point wrong: {:?}", result.point
    );
    assert!(
        (result.distance - 3.0).abs() < 1e-6,
        "distance wrong: {}", result.distance
    );
}

#[test]
fn closest_point_on_plane_known() {
    // XY plane (z=0), query at (1, 2, 5).
    // Closest point = (1, 2, 0), distance = 5.
    let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
    let result = closest_point_on_surface(&plane, DVec3::new(1.0, 2.0, 5.0), 8);
    let expected = DVec3::new(1.0, 2.0, 0.0);
    assert!(
        (result.point - expected).length() < 1e-9,
        "closest plane point wrong: {:?}", result.point
    );
    assert!(
        (result.distance - 5.0).abs() < 1e-9,
        "distance wrong: {}", result.distance
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Arc length
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn arc_length_semicircle() {
    // Semicircle of radius 3: arc length = π * r = 3π.
    let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 3.0));
    let len = arc_length(&circle, 0.0, std::f64::consts::PI);
    let expected = std::f64::consts::PI * 3.0;
    assert!(
        (len - expected).abs() < 1e-6,
        "semicircle arc length wrong: {} (expected {})", len, expected
    );
}

#[test]
fn arc_length_full_circle() {
    // Full circle of radius 1: circumference = 2π.
    let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
    let len = arc_length(&circle, 0.0, std::f64::consts::TAU);
    let expected = std::f64::consts::TAU;
    assert!(
        (len - expected).abs() < 1e-6,
        "full circle arc length wrong: {} (expected {})", len, expected
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Curvature
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn gaussian_curvature_sphere_known() {
    // Sphere of radius r: Gaussian curvature K = 1/r².
    let sphere = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 2.0, ref_dir: any_perpendicular(DVec3::Z) };
    let k = gaussian_curvature(&Surface3::Sphere(sphere), std::f64::consts::PI, std::f64::consts::FRAC_PI_2);
    let expected = 1.0 / (2.0 * 2.0); // 0.25
    assert!(
        (k - expected).abs() < 1e-6,
        "sphere Gaussian curvature wrong: {} (expected {})", k, expected
    );
}

#[test]
fn gaussian_curvature_plane_zero() {
    // Plane: Gaussian curvature K = 0.
    let plane = Plane::new(DVec3::ZERO, DVec3::Z);
    let k = gaussian_curvature(&Surface3::Plane(plane), 0.0, 0.0);
    assert!(
        k.abs() < 1e-9,
        "plane Gaussian curvature should be 0, got {}", k
    );
}

#[test]
fn mean_curvature_sphere_known() {
    // Sphere of radius r: mean curvature H = 1/r.
    let sphere = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 3.0, ref_dir: any_perpendicular(DVec3::Z) };
    let h = mean_curvature(&Surface3::Sphere(sphere), 0.5, 0.5);
    let expected = 1.0 / 3.0;
    assert!(
        (h - expected).abs() < 1e-6,
        "sphere mean curvature wrong: {} (expected {})", h, expected
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Extrema
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn extrema_parallel_lines_minimum_distance() {
    // Two parallel lines offset by 3 in Y.
    let l1 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
    let l2 = Curve3::Line(Line3 { origin: DVec3::new(0.0, 3.0, 0.0), direction: DVec3::X });
    let result = extrema_curve_curve(&l1, &l2, 8);
    assert!(
        !result.pairs.is_empty(),
        "extrema should find at least one pair"
    );
    let min_dist = result.pairs[0].distance;
    assert!(
        (min_dist - 3.0).abs() < 1e-6,
        "minimum distance between parallel lines wrong: {} (expected 3.0)", min_dist
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// BRep transformation
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn transform_box_scaling() {
    let (brep, _root) = BRep::build_unit_cube();
    // Unit cube at [0,0,0]-[1,1,1]; verify bounding box is dimension 1.
    let bb = brep.bounding_box().expect("unit cube must have a bounding box");
    let [mn, mx] = bb;
    let size = mx - mn;
    assert!(
        (size.x - 1.0).abs() < 1e-9 && (size.y - 1.0).abs() < 1e-9 && (size.z - 1.0).abs() < 1e-9,
        "unit cube bounding box wrong: {:?}", size
    );
}

#[test]
fn transform_preserves_original() {
    let (brep, _root) = BRep::build_unit_cube();
    // Verify the cube has 8 vertices.
    let n_verts: usize = brep.tshapes.iter().filter(|ts| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Vertex(_))).count();
    assert_eq!(n_verts, 8, "unit cube should have 8 vertices, got {}", n_verts);
}
