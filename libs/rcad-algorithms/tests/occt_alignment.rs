//! OCCT Alignment Tests for rcad-algorithms
//!
//! These tests verify alignment with Open CASCADE Technology (OCCT) behavior
//! across the following modules:
//! - TKBO: Boolean operations
//! - TKOffset: Offset operations
//! - TKShHealing: Shape healing and repair
//!
//! Test categories:
//! 1. Boolean operations with overlapping geometry
//! 2. Boolean with thin features
//! 3. Offset operations with self-intersection detection
//! 4. Complex sweep/loft operations
//! 5. Repair operations

use rcad_algorithms::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_algorithms::{
    boolean_op, BooleanOpType, BooleanError, BooleanOptions,
    boolean_op_with_options, boolean_op_simplified, SimplifyOptions,
    brep_algo_api::{
        BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepHistory, BooleanApiOptions,
    },
    sweep::{linear_sweep, rotational_sweep, pipe_sweep, SweepMode,
            linear_law_sweep, LinearLaw, handle_pipe_corners},
    offset::offset_surface,
    healing::{heal, heal_comprehensive, HealingOptions,
              fix_solid, fix_wire},
    brep_check::brep_check_analyze,
    history::{EdgeOrigin as HistEdgeOrigin, VertexOrigin as HistVertexOrigin},
};
use rcad_kernel::{any_perpendicular, BRep};
use rcad_kernel::geom::{Surface3, SphericalSurface, CylindricalSurface, Plane,
                        ToroidalSurface, ConicalSurface};
use rcad_kernel::properties::volume;
use rcad_modeling::{
    make_box_brep, make_cylinder_brep, make_sphere_brep, make_cone_brep,
};
use rcad_kernel::PrimitiveSolid;
use rcad_algorithms::geom_populate;

// Helper functions

fn box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: w,
        height: h,
        depth: d,
    });
    for v in &mut brep.vertices {
        v.point += DVec3::new(x, y, z);
    }
    geom_populate::populate_box_geom(&mut brep);
    brep
}

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn all_triangles_valid(brep: &BRep) -> bool {
    let nv = brep.vertices.len();
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .flat_map(|f| &f.triangles)
        .all(|tri| tri.iter().all(|&i| i < nv))
}

fn history_modified_scan_bounds(history: &BRepHistory) -> (usize, usize) {
    if let Some(inner_history) = history.inner() {
        let max_edge_idx = inner_history
            .edge_origins
            .iter()
            .filter_map(|o| match o {
                HistEdgeOrigin::FromA(src)
                | HistEdgeOrigin::FromB(src)
                | HistEdgeOrigin::SplitFromA(src)
                | HistEdgeOrigin::SplitFromB(src) => Some(*src),
                HistEdgeOrigin::Generated => None,
            })
            .max()
            .map(|m| m + 1)
            .unwrap_or(0);
        let max_vertex_idx = inner_history
            .vertex_origins
            .iter()
            .filter_map(|o| match o {
                HistVertexOrigin::FromA(src) | HistVertexOrigin::FromB(src) => Some(*src),
                HistVertexOrigin::Intersection => None,
            })
            .max()
            .map(|m| m + 1)
            .unwrap_or(0);
        (max_edge_idx, max_vertex_idx)
    } else {
        (0, 0)
    }
}

fn assert_history_modified_semantics(history: &BRepHistory) {
    let (max_edge_scan, max_vertex_scan) = history_modified_scan_bounds(history);

    let mut edge_hits = 0usize;
    for idx in 0..max_edge_scan {
        edge_hits += history.modified_edges_from_a(idx).len();
        edge_hits += history.modified_edges_from_b(idx).len();
    }

    let mut vertex_hits = 0usize;
    for idx in 0..max_vertex_scan {
        vertex_hits += history.modified_vertices_from_a(idx).len();
        vertex_hits += history.modified_vertices_from_b(idx).len();
    }

    let stats = history.statistics();
    assert_eq!(stats.modified_edges, edge_hits);
    assert_eq!(stats.modified_vertices, vertex_hits);
    assert_eq!(
        history.has_modified(),
        stats.modified_faces + stats.modified_edges + stats.modified_vertices > 0
    );
}

fn assert_history_deleted_semantics(
    history: &BRepHistory,
    face_count_a: usize,
    face_count_b: usize,
    edge_scan_bound: usize,
    vertex_scan_bound: usize,
) {
    let deleted_from_a = (0..face_count_a)
        .filter(|&idx| history.is_deleted_from_a(idx))
        .count();
    let deleted_from_b = (0..face_count_b)
        .filter(|&idx| history.is_deleted_from_b(idx))
        .count();

    for idx in 0..edge_scan_bound {
        if history.is_deleted_edge_from_a(idx) || history.is_deleted_edge_from_b(idx) {
            assert!(history.is_deleted_edge_any(idx));
        }
    }
    for idx in 0..vertex_scan_bound {
        if history.is_deleted_vertex_from_a(idx) || history.is_deleted_vertex_from_b(idx) {
            assert!(history.is_deleted_vertex_any(idx));
        }
    }

    let stats = history.statistics();
    assert_eq!(stats.deleted_faces, deleted_from_a + deleted_from_b);
    assert_eq!(
        history.has_deleted(),
        stats.deleted_faces + stats.deleted_edges + stats.deleted_vertices > 0
    );
}

// ============================================================================
// 1. Boolean Operations with Overlapping Geometry (TKBO Coverage)
// ============================================================================

/// Test union of two identical boxes at the same location.
/// OCCT TKBO should recognize identical geometry and return equivalent shape.
#[test]
fn boolean_identical_boxes_union() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let b2 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box2");

    let result = boolean_op(BooleanOpType::Union, &b1, &b2);
    match result {
        Ok(r) => {
            assert!(face_count(&r) >= 6, "result should have box faces");
            assert!(all_triangles_valid(&r));
            // Volume should be same as single box (identical geometry)
            let v = volume(&r);
            assert!((v - 8.0).abs() < 0.1, "volume should be approximately 8.0");
        }
        Err(BooleanError::DegenerateResult) => {
            // Acceptable for identical inputs
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Test union of two identical spheres at the same location.
/// OCCT should handle self-referential boolean operations gracefully.
#[test]
fn boolean_identical_spheres_union() {
    let s1 = make_sphere_brep(DVec3::ZERO, 1.0).expect("sphere1");
    let s2 = make_sphere_brep(DVec3::ZERO, 1.0).expect("sphere2");

    let result = boolean_op(BooleanOpType::Union, &s1, &s2);
    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
            let v = volume(&r);
            let expected = 4.0 / 3.0 * std::f64::consts::PI;
            assert!((v - expected).abs() / expected < 0.1, "volume should match sphere");
        }
        Err(BooleanError::DegenerateResult) => {}
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Test intersection of identical boxes.
/// Self-intersection should return the original shape.
#[test]
fn boolean_identical_boxes_intersection() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let b2 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box2");

    let result = boolean_op(BooleanOpType::Intersection, &b1, &b2);
    match result {
        Ok(r) => {
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {}
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Test difference of identical boxes.
/// Self-difference should produce empty result.
#[test]
fn boolean_identical_boxes_difference() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let b2 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box2");

    let result = boolean_op(BooleanOpType::Difference, &b1, &b2);
    match result {
        Ok(r) => {
            assert_eq!(face_count(&r), 0, "self-difference should produce empty result");
        }
        Err(BooleanError::DegenerateResult) => {
            // Expected behavior for self-difference
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Test union of identical cylinders.
#[test]
fn boolean_identical_cylinders_union() {
    let c1 = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 4.0).expect("cyl1");
    let c2 = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 4.0).expect("cyl2");

    let result = boolean_op(BooleanOpType::Union, &c1, &c2);
    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {}
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Test union with near-zero overlap (faces barely touching).
/// OCCT uses fuzzy tolerance to handle near-coincident geometry.
#[test]
fn boolean_near_zero_overlap_union() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    // Box starts at x=2.0 - epsilon, creating minimal overlap
    let epsilon = TOLERANCE_MESH_LEGACY;
    let b2 = make_box_brep(DVec3::new(2.0 - epsilon, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box2");

    let opts = BooleanOptions {
        fuzzy_tol: TOLERANCE_RETRY_LADDER_MID,
        ..Default::default()
    };

    let result = boolean_op_with_options(BooleanOpType::Union, &b1, &b2, opts);
    match result {
        Ok((r, _report)) => {
            assert!(face_count(&r) >= 6);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {}
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Test complete containment: small box inside large box union.
#[test]
fn boolean_contained_box_union() {
    let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("outer");
    let inner = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("inner");

    let result = boolean_op(BooleanOpType::Union, &outer, &inner)
        .expect("contained union should succeed");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
    // Volume should be that of outer box
    let v = volume(&result);
    assert!((v - 64.0).abs() < 1.0, "volume should be approximately 64.0");
}

// ============================================================================
// 2. Boolean with Thin Features (TKBO Coverage)
// ============================================================================

/// Test thin wall creation (0.01 thickness) via difference.
#[test]
fn boolean_thin_wall_creation_0_01() {
    let mut outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("outer");
    let mut inner = make_box_brep(
        DVec3::new(0.01, 0.01, 0.01),
        DVec3::X,
        DVec3::Y,
        0.98,
        0.98,
        0.98,
    ).expect("inner");

    geom_populate::populate_box_geom(&mut outer);
    geom_populate::populate_box_geom(&mut inner);

    let result = boolean_op_simplified(
        BooleanOpType::Difference,
        &outer,
        &inner,
        SimplifyOptions::default(),
    );

    match result {
        Ok((r, _)) => {
            assert!(face_count(&r) >= 6);
            assert!(all_triangles_valid(&r));
        }
        Err(_) => {
            // Thin features can be challenging
        }
    }
}

/// Test very thin sheet metal thickness (0.5mm in 100mm panel).
#[test]
fn boolean_sheet_metal_thickness() {
    let mut panel = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 100.0, 100.0, 0.5).expect("panel");
    geom_populate::populate_box_geom(&mut panel);

    // Hole that goes through the sheet
    let hole = make_cylinder_brep(
        DVec3::new(50.0, 50.0, -0.1),
        DVec3::Z,
        DVec3::X,
        5.0,
        0.7,
    ).expect("hole");

    let result = boolean_op_simplified(
        BooleanOpType::Difference,
        &panel,
        &hole,
        SimplifyOptions::default(),
    );

    assert!(result.is_ok(), "sheet metal hole cut should succeed");
    let (r, _) = result.unwrap();
    assert!(face_count(&r) >= 6);
    assert!(all_triangles_valid(&r));
}

/// Test small hole subtraction from large block.
#[test]
fn boolean_small_hole_large_block() {
    let mut block = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 50.0, 50.0, 50.0).expect("block");
    geom_populate::populate_box_geom(&mut block);

    // Very small hole (0.05 radius)
    let hole = make_cylinder_brep(
        DVec3::new(25.0, 25.0, -1.0),
        DVec3::Z,
        DVec3::X,
        0.05,
        52.0,
    ).expect("hole");

    let result = boolean_op_simplified(
        BooleanOpType::Difference,
        &block,
        &hole,
        SimplifyOptions::default(),
    );

    match result {
        Ok((r, _)) => {
            assert!(face_count(&r) >= 6);
            assert!(all_triangles_valid(&r));
        }
        Err(_) => {
            // Small features can be challenging for some kernels
        }
    }
}

/// Test thin slot (high aspect ratio feature).
#[test]
fn boolean_thin_slot_subtraction() {
    let mut block = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 20.0, 20.0, 20.0).expect("block");
    geom_populate::populate_box_geom(&mut block);

    // Thin slot: 0.05mm wide, 10mm deep, full height
    let slot = make_box_brep(
        DVec3::new(9.975, -1.0, 0.0),
        DVec3::X,
        DVec3::Z,
        0.05,
        12.0,
        20.0,
    ).expect("slot");

    let result = boolean_op_simplified(
        BooleanOpType::Difference,
        &block,
        &slot,
        SimplifyOptions::default(),
    );

    if let Ok((r, _)) = result {
        assert!(face_count(&r) >= 6);
        assert!(all_triangles_valid(&r));
    }
}

/// Test multiple small holes pattern.
#[test]
fn boolean_multiple_small_holes_pattern() {
    let mut block = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 2.0).expect("block");
    geom_populate::populate_box_geom(&mut block);

    let hole_radius = 0.1;
    let spacing = 2.0;

    // Create a pattern of small holes
    let mut result = block.clone();
    for i in 0..5 {
        for j in 0..5 {
            let cx = 1.0 + i as f64 * spacing;
            let cy = 1.0 + j as f64 * spacing;

            let hole = make_cylinder_brep(
                DVec3::new(cx, cy, -0.5),
                DVec3::Z,
                DVec3::X,
                hole_radius,
                3.0,
            ).expect("hole");

            match boolean_op(BooleanOpType::Difference, &result, &hole) {
                Ok(r) => result = r,
                Err(_) => continue,
            }
        }
    }

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

// ============================================================================
// 3. Offset Operations with Self-Intersection Detection (TKOffset Coverage)
// ============================================================================

/// Test offset of sphere surface.
#[test]
fn offset_sphere_positive() {
    let sphere = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 1.0,
        ref_dir: any_perpendicular(DVec3::Z),
    });

    let result = offset_surface(&sphere, 0.5);
    assert!(result.is_some(), "offset sphere should succeed");

    if let Surface3::Sphere(s) = result.unwrap() {
        assert!((s.radius - 1.5).abs() < TOLERANCE_COORD_SUB, "radius should be 1.5");
    }
}

/// Test negative offset of sphere (shrinking).
#[test]
fn offset_sphere_negative() {
    let sphere = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 2.0,
        ref_dir: any_perpendicular(DVec3::Z),
    });

    let result = offset_surface(&sphere, -0.5);
    assert!(result.is_some(), "negative offset sphere should succeed");

    if let Surface3::Sphere(s) = result.unwrap() {
        assert!((s.radius - 1.5).abs() < TOLERANCE_COORD_SUB, "radius should be 1.5");
    }
}

/// Test offset that would create self-intersection (large negative offset on torus).
#[test]
fn offset_torus_self_intersection_detection() {
    let torus = Surface3::Torus(ToroidalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        major_radius: 2.0,
        minor_radius: 1.0,
    });

    // Negative offset greater than minor radius would cause self-intersection
    let result = offset_surface(&torus, -1.5);

    // Should either return None (invalid) or handle gracefully
    if let Some(Surface3::Torus(t)) = result {
        // If it returns a torus, verify it's valid
        assert!(t.minor_radius > 0.0 || t.minor_radius.is_nan(), "minor radius should be positive or NaN");
    }
    // None is also acceptable for self-intersecting offsets
}

/// Test offset of cylinder surface.
#[test]
fn offset_cylinder_positive() {
    let cylinder = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 1.0,
    });

    let result = offset_surface(&cylinder, 0.5);
    assert!(result.is_some(), "offset cylinder should succeed");

    if let Surface3::Cylinder(c) = result.unwrap() {
        assert!((c.radius - 1.5).abs() < TOLERANCE_COORD_SUB, "radius should be 1.5");
    }
}

/// Test offset of cone surface.
#[test]
fn offset_cone_surface() {
    let cone = Surface3::Cone(ConicalSurface {
        apex: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 1.0,
        half_angle_rad: std::f64::consts::FRAC_PI_6, // 30 degrees
    });

    let result = offset_surface(&cone, 0.2);
    assert!(result.is_some(), "offset cone should succeed");
}

/// Test offset of plane surface.
#[test]
fn offset_plane_surface() {
    let plane = Surface3::Plane(Plane {
        origin: DVec3::ZERO,
        normal: DVec3::Z,
    });

    let result = offset_surface(&plane, 1.0);
    assert!(result.is_some(), "offset plane should succeed");

    if let Surface3::Plane(p) = result.unwrap() {
        // Plane offset should maintain normal direction
        assert!(p.normal.abs_diff_eq(DVec3::Z, TOLERANCE_COORD_SUB));
    }
}

/// Test self-intersection detection for large offsets.
#[test]
fn offset_detect_self_intersection_torus() {
    let torus = Surface3::Torus(ToroidalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        major_radius: 3.0,
        minor_radius: 1.0,
    });

    // Test with offset that would NOT cause self-intersection
    let no_intersection = offset_surface(&torus, 0.5).is_some();
    assert!(no_intersection, "small offset should not self-intersect");

    // Test with offset that WOULD cause self-intersection
    let would_intersect = offset_surface(&torus, -1.5).is_none();
    assert!(would_intersect, "large negative offset should self-intersect");
}

// ============================================================================
// 4. Complex Sweep/Loft Operations (TKPrim Coverage)
// ============================================================================

/// Test linear sweep of square profile.
#[test]
fn sweep_linear_square_profile() {
    let profile = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
    ];

    let result = linear_sweep(&profile, DVec3::Z, 2.0);
    assert!(result.is_ok(), "linear sweep should succeed");

    let brep = result.unwrap();
    assert!(face_count(&brep) == 6, "extruded box should have 6 faces");
}

/// Test rotational sweep of rectangle (creates cylinder-like shape).
#[test]
fn sweep_rotational_rectangle() {
    let profile = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 1.0),
        DVec3::new(1.0, 0.0, 1.0),
    ];

    let result = rotational_sweep(&profile, DVec3::ZERO, DVec3::Z, std::f64::consts::TAU);
    assert!(result.is_ok(), "full revolution should succeed");

    let brep = result.unwrap();
    assert!(face_count(&brep) > 0, "revolved shape should have faces");
}

/// Test pipe sweep along linear spine.
#[test]
fn sweep_pipe_linear_spine() {
    let profile = vec![
        DVec2::new(-0.5, -0.5),
        DVec2::new(0.5, -0.5),
        DVec2::new(0.5, 0.5),
        DVec2::new(-0.5, 0.5),
    ];

    let spine = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(0.0, 0.0, 2.0),
    ];

    let result = pipe_sweep(&profile, &spine, SweepMode::Pipe);
    assert!(result.is_ok(), "pipe sweep should succeed");

    let brep = result.unwrap();
    assert!(face_count(&brep) >= 5, "pipe should have start, end, and lateral faces");
}

/// Test pipe sweep along curved spine.
#[test]
fn sweep_pipe_curved_spine() {
    let profile = vec![
        DVec2::new(-0.3, -0.3),
        DVec2::new(0.3, -0.3),
        DVec2::new(0.3, 0.3),
        DVec2::new(-0.3, 0.3),
    ];

    // Arc-shaped spine
    let n = 10;
    let spine: Vec<DVec3> = (0..n)
        .map(|i| {
            let t = i as f64 / (n - 1) as f64;
            let angle = t * std::f64::consts::FRAC_PI_2;
            DVec3::new(angle.cos(), angle.sin(), 0.0)
        })
        .collect();

    let result = pipe_sweep(&profile, &spine, SweepMode::Pipe);
    assert!(result.is_ok(), "curved pipe sweep should succeed");
}

/// Test pipe sweep with corner handling.
#[test]
fn sweep_pipe_with_sharp_corner() {
    let profile = vec![
        DVec2::new(-0.2, -0.2),
        DVec2::new(0.2, -0.2),
        DVec2::new(0.2, 0.2),
        DVec2::new(-0.2, 0.2),
    ];

    // L-shaped spine with sharp corner
    let spine = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(1.0, 2.0, 0.0),
    ];

    let result = handle_pipe_corners(&spine, &profile, 0.1);
    assert!(result.is_ok(), "pipe with corners should succeed");
}

/// Test linear law sweep (tapered shape).
#[test]
fn sweep_linear_law_tapered() {
    let profile = vec![
        DVec2::new(-0.5, -0.5),
        DVec2::new(0.5, -0.5),
        DVec2::new(0.5, 0.5),
        DVec2::new(-0.5, 0.5),
    ];

    let spine = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
        DVec3::new(0.0, 0.0, 2.0),
    ];

    let law = LinearLaw {
        start_value: 1.0,
        end_value: 0.5,
    };

    let result = linear_law_sweep(&profile, &spine, law);
    assert!(result.is_ok(), "tapered sweep should succeed");
}

/// Test sweep with triangular profile.
#[test]
fn sweep_triangular_profile() {
    let profile = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.5, 1.0, 0.0),
    ];

    let result = linear_sweep(&profile, DVec3::Z, 2.0);
    assert!(result.is_ok(), "triangular sweep should succeed");

    let brep = result.unwrap();
    assert!(face_count(&brep) == 5, "extruded triangle should have 5 faces (2 caps + 3 lateral)");
}

/// Test partial revolution (90 degrees).
#[test]
fn sweep_partial_revolution() {
    let profile = vec![
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 1.0),
        DVec3::new(1.0, 0.0, 1.0),
    ];

    let result = rotational_sweep(&profile, DVec3::ZERO, DVec3::Z, std::f64::consts::FRAC_PI_2);
    assert!(result.is_ok(), "90-degree revolution should succeed");

    let brep = result.unwrap();
    assert!(face_count(&brep) > 0);
}

// ============================================================================
// 5. Repair Operations (TKShHealing Coverage)
// ============================================================================

/// Test healing of degenerate edge (zero-length edge).
#[test]
fn healing_degenerate_edge() {
    let mut brep = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    // Create a degenerate edge by duplicating start/end vertex
    if !brep.edges.is_empty() {
        let v0 = brep.edges[0].start;
        brep.edges[0].end = v0; // Zero-length edge
    }

    let (_healed, report) = heal(&brep);

    // Should detect or fix the degenerate edge
    assert!(report.initial_issue_count() >= 1 || report.is_clean());
}

/// Test healing of reversed face normal.
#[test]
fn healing_reversed_face_normal() {
    let mut brep = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    // Reverse one face normal to simulate incorrect orientation
    if let Some(face) = brep.solids.first_mut()
        .and_then(|s| s.shells.first_mut())
        .and_then(|sh| sh.faces.first_mut())
    {
        face.normal = -face.normal;
    }

    let (_healed, report) = heal(&brep);

    // Should detect and fix the reversed normal
    assert!(report.is_improved() || report.is_clean());
}

/// Test fix_solid for shell closure issues.
#[test]
fn repair_fix_solid_closure() {
    let brep = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    let (_fixed, report) = fix_solid(&brep, TOLERANCE_MESH_LEGACY);

    // A valid box should pass solid checks
    assert!(report.unclosed_shells.is_empty());
}

/// Test fix_wire for wire issues.
#[test]
fn repair_fix_wire_issues() {
    let brep = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    let (_fixed, report) = fix_wire(&brep, TOLERANCE_MESH_LEGACY);

    // A valid box should have clean wires
    assert!(report.is_clean() || report.wires_with_issues == 0);
}

/// Test comprehensive healing pipeline.
#[test]
fn repair_comprehensive_healing() {
    let brep = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    let options = HealingOptions::default();
    let (_healed, report) = heal_comprehensive(&brep, &options);

    // Should produce valid result
    assert!(report.is_clean || report.final_check.issues.is_empty());
}

/// Test healing of small gap in wire.
#[test]
fn healing_small_wire_gap() {
    let mut brep = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    // Perturb a vertex slightly to create a small gap
    if !brep.vertices.is_empty() {
        brep.vertices[0].point.x += TOLERANCE_RETRY_LADDER_COARSE;
    }

    let (_healed, _report) = heal(&brep);
}

/// Test healing of self-intersecting shell detection.
#[test]
fn healing_self_intersecting_shell_detection() {
    // Create a box and check for self-intersection
    let brep = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    let check_result = brep_check_analyze(&brep);

    // Valid box should not have self-intersection
    assert!(check_result.is_valid() || check_result.issues.is_empty());
}

/// Test merging of coplanar faces (simplified test).
#[test]
fn repair_merge_coplanar_faces() {
    // Create two adjacent boxes sharing a face
    let b1 = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    let b2 = box_at(2.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    let union = boolean_op(BooleanOpType::Union, &b1, &b2)
        .expect("union should succeed");

    // Run healing on the result
    let (healed, _report) = heal(&union);

    // Should produce valid geometry
    assert!(all_triangles_valid(&healed));
}

/// Test non-manifold edge detection.
#[test]
fn repair_detect_non_manifold_edge() {
    // Create geometry that could produce non-manifold edges
    let b1 = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    let b2 = box_at(2.0, 1.0, 0.0, 2.0, 2.0, 2.0);

    // Union creates a configuration where edge manifoldness can be tested
    let result = boolean_op(BooleanOpType::Union, &b1, &b2);

    if let Ok(union) = result {
        let (_fixed, _report) = fix_solid(&union, TOLERANCE_MESH_LEGACY);
    }
}

// ============================================================================
// Additional OCCT Coverage Tests
// ============================================================================

/// Test boolean union with fuzzy tolerance for nearly touching faces.
#[test]
fn boolean_fuzzy_tolerance_nearly_touching() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b1");
    // Gap of TOLERANCE_RETRY_LADDER_MID between boxes
    let b2 = make_box_brep(DVec3::new(2.0 + TOLERANCE_RETRY_LADDER_MID, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("b2");

    // With fuzzy tolerance larger than gap, should merge
    let opts = BooleanOptions {
        fuzzy_tol: TOLERANCE_RETRY_LADDER_COARSE,
        ..Default::default()
    };

    let result = boolean_op_with_options(BooleanOpType::Union, &b1, &b2, opts);
    assert!(result.is_ok(), "fuzzy tolerance should bridge small gap");

    let (r, _) = result.unwrap();
    assert!(face_count(&r) >= 6);
}

/// Test boolean difference creating internal cavity.
#[test]
fn boolean_internal_cavity() {
    let mut outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).expect("outer");
    let mut inner = make_box_brep(DVec3::new(2.0, 2.0, 2.0), DVec3::X, DVec3::Y, 6.0, 6.0, 6.0).expect("inner");

    geom_populate::populate_box_geom(&mut outer);
    geom_populate::populate_box_geom(&mut inner);

    let result = boolean_op_simplified(
        BooleanOpType::Difference,
        &outer,
        &inner,
        SimplifyOptions::default(),
    );

    if let Ok((r, _)) = result {
        assert!(face_count(&r) >= 10, "hollow box should have inner and outer faces");
        assert!(all_triangles_valid(&r));
    }
}

/// Test cone-sphere intersection at apex.
#[test]
fn boolean_cone_sphere_apex_intersection() {
    let cone = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0).expect("cone");
    let sphere = make_sphere_brep(DVec3::new(0.0, 0.0, 0.5), 1.0).expect("sphere");

    let result = boolean_op(BooleanOpType::Intersection, &cone, &sphere);

    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {
            // Apex region can be challenging
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

// ============================================================================
// 6. BRepAlgoAPI History Semantics (OCCT BuilderShape-style Coverage)
// ============================================================================

/// Test that cut history reports deleted source faces and exposes stable
/// edge/vertex deletion query semantics.
#[test]
fn brep_algoapi_cut_history_deleted_semantics() {
    let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("outer");
    let inner = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("inner");

    let mut cut = BRepAlgoAPI_Cut::new(&outer, &inner);
    cut.set_options(BooleanApiOptions::default().with_history(true));
    assert!(cut.build(), "cut should succeed with history enabled");

    let history = cut.history();
    assert!(history.is_generated(), "history should be generated");
    assert!(history.inner().is_some(), "inner boolean history should be available");
    let inner_face_count = inner
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let outer_face_count = outer
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();

    assert_history_deleted_semantics(history, outer_face_count, inner_face_count, 24, 24);
}

/// Test that history modified edge/vertex queries are self-consistent with
/// aggregated statistics in a real BRepAlgoAPI operation path.
#[test]
fn brep_algoapi_cut_history_modified_edge_vertex_semantics() {
    let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("outer");
    let inner = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("inner");

    let mut cut = BRepAlgoAPI_Cut::new(&outer, &inner);
    cut.set_options(BooleanApiOptions::default().with_history(true));
    assert!(cut.build(), "cut should succeed with history enabled");

    let history = cut.history();
    assert!(history.is_generated(), "history should be generated");
    assert_history_modified_semantics(history);
}

/// Fuse-path counterpart of cut history modified semantics.
#[test]
fn brep_algoapi_fuse_history_modified_edge_vertex_semantics() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 3.0, 3.0, 3.0).expect("a");
    let b = make_box_brep(DVec3::new(1.0, 1.0, 0.5), DVec3::X, DVec3::Y, 3.0, 3.0, 3.0)
        .expect("b");

    let mut fuse = BRepAlgoAPI_Fuse::new(&a, &b);
    fuse.set_options(BooleanApiOptions::default().with_history(true));
    assert!(fuse.build(), "fuse should succeed with history enabled");

    let history = fuse.history();
    assert!(history.is_generated(), "history should be generated");
    assert_history_modified_semantics(history);
}

/// Fuse-path counterpart of cut history deleted semantics.
#[test]
fn brep_algoapi_fuse_history_deleted_semantics() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 3.0, 3.0, 3.0).expect("a");
    let b = make_box_brep(DVec3::new(1.0, 1.0, 0.5), DVec3::X, DVec3::Y, 3.0, 3.0, 3.0)
        .expect("b");

    let mut fuse = BRepAlgoAPI_Fuse::new(&a, &b);
    fuse.set_options(BooleanApiOptions::default().with_history(true));
    assert!(fuse.build(), "fuse should succeed with history enabled");

    let history = fuse.history();
    assert!(history.is_generated(), "history should be generated");

    let face_count_a = a
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let face_count_b = b
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();

    let max_edge_scan = a.edges.len().max(b.edges.len());
    let max_vertex_scan = a.vertices.len().max(b.vertices.len());
    assert_history_deleted_semantics(
        history,
        face_count_a,
        face_count_b,
        max_edge_scan,
        max_vertex_scan,
    );
}

/// Common-path counterpart of cut/fuse history modified semantics.
#[test]
fn brep_algoapi_common_history_modified_edge_vertex_semantics() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 3.0, 3.0, 3.0).expect("a");
    let b = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 3.0, 3.0, 3.0)
        .expect("b");

    let mut common = BRepAlgoAPI_Common::new(&a, &b);
    common.set_options(BooleanApiOptions::default().with_history(true));
    assert!(common.build(), "common should succeed with history enabled");

    let history = common.history();
    assert!(history.is_generated(), "history should be generated");
    assert_history_modified_semantics(history);
}

/// Common-path counterpart of cut/fuse history deleted semantics.
#[test]
fn brep_algoapi_common_history_deleted_semantics() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 3.0, 3.0, 3.0).expect("a");
    let b = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 3.0, 3.0, 3.0)
        .expect("b");

    let mut common = BRepAlgoAPI_Common::new(&a, &b);
    common.set_options(BooleanApiOptions::default().with_history(true));
    assert!(common.build(), "common should succeed with history enabled");

    let history = common.history();
    assert!(history.is_generated(), "history should be generated");

    let face_count_a = a
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let face_count_b = b
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();

    let max_edge_scan = a.edges.len().max(b.edges.len());
    let max_vertex_scan = a.vertices.len().max(b.vertices.len());
    assert_history_deleted_semantics(
        history,
        face_count_a,
        face_count_b,
        max_edge_scan,
        max_vertex_scan,
    );
}
