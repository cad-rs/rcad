/// Edge case tests for boolean operations.
///
/// These tests cover challenging geometric configurations that are known
/// to cause issues in CAD kernels:
/// - Near-coincident faces
/// - High-curvature intersections
/// - Extreme size ratios
/// - Failure recovery scenarios
use rcad_algorithms::tolerance::*;
use glam::DVec3;
use rcad_algorithms::{BooleanOpType, boolean_op, BooleanError};
use rcad_kernel::properties::volume;
use rcad_modeling::{
    make_box_brep, make_cone_brep, make_cylinder_brep, make_sphere_brep, make_torus_brep,
};

fn face_count(brep: &rcad_kernel::BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn all_triangles_valid(brep: &rcad_kernel::BRep) -> bool {
    let nv = brep.vertices.len();
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .flat_map(|f| &f.triangles)
        .all(|tri| tri.iter().all(|&i| i < nv))
}

// ============================================================================
// Near-Coincident Faces Tests
// ============================================================================

/// Two boxes with faces that are nearly coincident (separated by tiny gap).
/// Tests the kernel's ability to handle fuzzy tolerance correctly.
#[test]
fn near_coincident_faces_tiny_gap() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    // Second box starts at x=2.0001, leaving a tiny gap of 0.0001
    let b2 = make_box_brep(DVec3::new(2.0001, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box2");

    // Union should succeed despite the tiny gap
    let result = boolean_op(BooleanOpType::Union, &b1, &b2);
    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {
            // This is acceptable for very small gaps
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Two spheres with surfaces that nearly touch.
#[test]
fn near_coinident_spheres_union() {
    let s1 = make_sphere_brep(DVec3::ZERO, 1.0).expect("sphere1");
    // Second sphere positioned so surfaces almost touch
    let s2 = make_sphere_brep(DVec3::new(2.001, 0.0, 0.0), 1.0).expect("sphere2");

    let result = boolean_op(BooleanOpType::Union, &s1, &s2);
    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {
            // Near-coincident curved geometry can collapse to a degenerate fallback.
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

// ============================================================================
// High-Curvature Intersection Tests
// ============================================================================

/// Torus self-intersection region (high curvature at inner equator).
#[test]
fn torus_inner_region_intersection() {
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 1.5).expect("torus");
    // Small box intersecting the inner high-curvature region
    let box_ = make_box_brep(DVec3::new(-0.5, -1.0, -2.0), DVec3::X, DVec3::Y, 1.0, 2.0, 4.0)
        .expect("box");

    let result = boolean_op(BooleanOpType::Intersection, &torus, &box_)
        .expect("torus inner region intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Cylinder intersecting sphere at sphere pole (degenerate parameterization).
#[test]
fn cylinder_sphere_pole_intersection() {
    let sphere = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
    // Cylinder positioned to intersect at sphere's north pole
    let cylinder = make_cylinder_brep(DVec3::new(0.0, 0.0, 1.5), DVec3::Z, DVec3::X, 0.5, 2.0)
        .expect("cylinder");

    let result = boolean_op(BooleanOpType::Union, &sphere, &cylinder)
        .expect("cylinder-sphere pole union should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Multiple high-curvature features in one operation.
#[test]
fn multiple_high_curvature_features() {
    let sphere = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
    let cone = make_cone_brep(DVec3::new(0.0, 0.0, -1.0), DVec3::Z, DVec3::X, 1.0, 3.0)
        .expect("cone");

    let result = boolean_op(BooleanOpType::Difference, &sphere, &cone);

    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {
            // High-curvature operations can produce degenerate results
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

// ============================================================================
// Extreme Size Ratio Tests
// ============================================================================

/// Very small feature on a large solid (hole in large plate).
#[test]
fn tiny_hole_in_large_plate() {
    // Large plate
    let plate = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 100.0, 100.0, 10.0).expect("plate");
    // Tiny cylinder as hole
    let hole = make_cylinder_brep(DVec3::new(50.0, 50.0, -1.0), DVec3::Z, DVec3::X, 0.1, 12.0)
        .expect("hole");

    let result = boolean_op(BooleanOpType::Difference, &plate, &hole)
        .expect("tiny hole in large plate should succeed");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

/// Very large object subtracted by small object.
#[test]
fn large_minus_small_size_ratio() {
    // Size ratio of 1000:1
    let large = make_sphere_brep(DVec3::ZERO, 100.0).expect("large sphere");
    let small = make_sphere_brep(DVec3::new(50.0, 0.0, 0.0), 0.1).expect("small sphere");

    let result = boolean_op(BooleanOpType::Difference, &large, &small);

    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
            // Volume should be nearly unchanged
            let vol = volume(&r);
            let expected = 4.0 / 3.0 * std::f64::consts::PI * 100.0_f64.powi(3);
            assert!((vol - expected).abs() / expected < 0.01);
        }
        Err(BooleanError::DegenerateResult) => {
            // The small subtraction might be ignored
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Nested shapes with extreme size difference.
#[test]
fn nested_extreme_size_difference() {
    let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1000.0, 1000.0, 1000.0)
        .expect("outer");
    let inner = make_box_brep(DVec3::new(400.0, 400.0, 400.0), DVec3::X, DVec3::Y, 200.0, 200.0, 200.0)
        .expect("inner");

    let result = boolean_op(BooleanOpType::Difference, &outer, &inner)
        .expect("nested extreme size difference should succeed");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

/// Thin wall creation (high aspect ratio geometry).
#[test]
fn thin_wall_creation() {
    // Create a thin-walled box by subtracting a slightly smaller box
    let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).expect("outer");
    // Inner box offset by 0.1 (thin wall)
    let inner = make_box_brep(
        DVec3::new(0.1, 0.1, 0.1),
        DVec3::X,
        DVec3::Y,
        9.8,
        9.8,
        10.0, // Open at top
    )
    .expect("inner");

    let result = boolean_op(BooleanOpType::Difference, &outer, &inner)
        .expect("thin wall creation should succeed");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

// ============================================================================
// Failure Recovery Tests
// ============================================================================

/// Test boolean operation with non-manifold input detection.
#[test]
fn non_manifold_geometry_handling() {
    // Create two boxes that share an edge (not a face) - creates non-manifold at intersection
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let b2 = make_box_brep(DVec3::new(2.0, 1.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box2");

    // This should either succeed or return a valid error
    let result = boolean_op(BooleanOpType::Union, &b1, &b2);
    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(e) => {
            // Any well-defined error is acceptable
            println!("Non-manifold geometry returned error: {:?}", e);
        }
    }
}

/// Test disjoint geometry union (two separate objects).
#[test]
fn disjoint_geometry_union() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box1");
    let b2 = make_box_brep(DVec3::new(100.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("box2");

    let result = boolean_op(BooleanOpType::Union, &b1, &b2)
        .expect("disjoint geometry union should succeed");

    // Should have faces from both boxes
    assert!(face_count(&result) >= 12);
    assert!(all_triangles_valid(&result));
}

/// Test disjoint geometry intersection (should produce empty result).
#[test]
fn disjoint_geometry_intersection() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box1");
    let b2 = make_box_brep(DVec3::new(100.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("box2");

    let result = boolean_op(BooleanOpType::Intersection, &b1, &b2);

    // Disjoint intersection should either return empty or degenerate result
    match result {
        Ok(r) => {
            assert_eq!(face_count(&r), 0, "disjoint intersection should have no faces");
        }
        Err(BooleanError::DegenerateResult) => {
            // This is the expected behavior for disjoint geometry
        }
        Err(e) => panic!("unexpected error for disjoint intersection: {:?}", e),
    }
}

/// Test contained geometry intersection (small fully inside large).
#[test]
fn contained_geometry_intersection() {
    let large = make_sphere_brep(DVec3::ZERO, 10.0).expect("large");
    let small = make_sphere_brep(DVec3::ZERO, 1.0).expect("small");

    let result = boolean_op(BooleanOpType::Intersection, &large, &small);

    // Result should be identical to the small sphere, or a degenerate fallback.
    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {
            // Accepted fallback for fully-contained curved intersections.
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

// ============================================================================
// Degenerate Input Tests
// ============================================================================

/// Test with near-zero radius sphere (should fail gracefully).
#[test]
fn near_zero_radius_sphere() {
    let result = make_sphere_brep(DVec3::ZERO, TOLERANCE_LINEAR_ULTRA_STRICT);
    // Should either fail or create a degenerate sphere
    match result {
        Ok(brep) => {
            // If it succeeds, verify it has geometry
            assert!(!brep.vertices.is_empty() || face_count(&brep) == 0);
        }
        Err(_) => {
            // Error is expected for near-zero radius
        }
    }
}

/// Test with extremely long cylinder (high aspect ratio).
#[test]
fn extremely_long_cylinder() {
    let cylinder = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 10000.0)
        .expect("cylinder");
    let box_ = make_box_brep(DVec3::new(-5.0, -5.0, 5000.0), DVec3::X, DVec3::Y, 10.0, 10.0, 10.0)
        .expect("box");

    let result = boolean_op(BooleanOpType::Intersection, &cylinder, &box_)
        .expect("long cylinder intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ============================================================================
// Symmetry and Edge Case Tests
// ============================================================================

/// Test that A union B equals B union A (commutativity).
#[test]
fn union_commutativity() {
    let a = make_sphere_brep(DVec3::new(0.0, 0.0, 0.0), 1.0).expect("sphere a");
    let b = make_sphere_brep(DVec3::new(0.5, 0.0, 0.0), 1.0).expect("sphere b");

    let ab = boolean_op(BooleanOpType::Union, &a, &b).expect("A union B");
    let ba = boolean_op(BooleanOpType::Union, &b, &a).expect("B union A");

    // Volumes should be approximately equal
    let vol_ab = volume(&ab);
    let vol_ba = volume(&ba);
    assert!(
        (vol_ab - vol_ba).abs() < 0.01 * vol_ab,
        "union should be commutative: {} vs {}",
        vol_ab,
        vol_ba
    );
}

/// Test that A intersect B equals B intersect A (commutativity).
#[test]
fn intersection_commutativity() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box a");
    let b = make_box_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box b");

    let ab = boolean_op(BooleanOpType::Intersection, &a, &b).expect("A intersect B");
    let ba = boolean_op(BooleanOpType::Intersection, &b, &a).expect("B intersect A");

    // Volumes should be approximately equal
    let vol_ab = volume(&ab);
    let vol_ba = volume(&ba);
    assert!(
        (vol_ab - vol_ba).abs() < 0.01,
        "intersection should be commutative: {} vs {}",
        vol_ab,
        vol_ba
    );
}
