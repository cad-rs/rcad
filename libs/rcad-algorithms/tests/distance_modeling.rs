//! Point/brep distance and shape–shape distance on **modeling-built** BReps.
//!
//! Note: `point_to_shape_distance` / `closest_on_brep` use analytic **unbounded** surfaces per face, so
//! a query point near a **box** may minimize distance to a face plane (not the finite rectangle); see
//! `rcad-kernel` `distance` docs. Spheres are exact for the current pipeline.

use glam::DVec3;
use rcad_kernel::distance::{min_distance, point_to_shape_distance};
use rcad_modeling::{make_box_brep, make_sphere_brep};

#[test]
fn point_to_modeling_sphere_separated() {
    let s = make_sphere_brep(DVec3::ZERO, 1.0).expect("sphere");
    let d = point_to_shape_distance(DVec3::new(5.0, 0.0, 0.0), &s);
    assert!((d.distance - 4.0).abs() < 0.15, "r=1, point at x=5 → distance 4, got {}", d.distance);
}

#[test]
fn point_outside_modeling_box_distance_sane() {
    // Do not assert clamped-rectangle distance; just ensure finite positive distance.
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");
    let d = point_to_shape_distance(DVec3::new(10.0, 10.0, 10.0), &b);
    assert!(
        d.distance > 0.0 && d.distance < 20.0,
        "sane dist, got {} (closest {})",
        d.distance,
        d.point_on_b
    );
}

#[test]
fn min_distance_modeling_spheres_separated() {
    let a = make_sphere_brep(DVec3::ZERO, 1.0).expect("a");
    let b = make_sphere_brep(DVec3::new(5.0, 0.0, 0.0), 1.0).expect("b");
    let d = min_distance(&a, &b);
    assert!((d.distance - 3.0).abs() < 0.1, "two unit spheres d=5 → gap 3, got {}", d.distance);
}
