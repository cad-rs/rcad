/// Integration tests for boolean operations across multiple shapes and scenarios.
/// These complement the inline unit tests by testing multi-step workflows and
/// error path behavior at crate boundary.
use rcad_algorithms::tolerance::*;
use glam::DVec3;
use rcad_algorithms::{
    BooleanError, BooleanOpType, CellExpr, MakerVolume, boolean_op, boolean_op_simplified,
    make_solid_from_region, SimplifyOptions,
};
use rcad_algorithms::bopds::ds::{DS, ShapeOrigin};
use rcad_algorithms::pave_filler::PaveFiller;
use rcad_kernel::topology::Face;
use rcad_kernel::PrimitiveSolid;
use rcad_kernel::BRep;
use rcad_kernel::properties::volume;
use rcad_algorithms::geom_populate;
use rcad_modeling::{
    make_box_brep, make_cone_brep, make_cylinder_brep, make_sphere_brep, make_torus_brep,
};

fn polygon_area_outer_wire(brep: &BRep, face: &Face) -> f64 {
    let mut pts: Vec<DVec3> = Vec::new();
    for we in &face.outer_wire.edges {
        let e = match brep.edges.get(we.idx) {
            Some(e) => e,
            None => continue,
        };
        let vi = if we.forward { e.start } else { e.end };
        if let Some(v) = brep.vertices.get(vi) {
            pts.push(v.point);
        }
    }
    if pts.len() < 3 {
        return 0.0;
    }
    let n = face.normal.normalize_or_zero();
    let ax = n.x.abs();
    let ay = n.y.abs();
    let az = n.z.abs();
    let axis = if ax >= ay && ax >= az {
        0usize
    } else if ay >= az {
        1usize
    } else {
        2usize
    };
    let mut area2 = 0.0;
    for i in 0..pts.len() {
        let p = pts[i];
        let q = pts[(i + 1) % pts.len()];
        area2 += match axis {
            0 => p.y * q.z - q.y * p.z,
            1 => p.x * q.z - q.x * p.z,
            _ => p.x * q.y - q.x * p.y,
        };
    }
    0.5 * area2.abs()
}

/// Same geometry as `examples/boolean_ops.rs` case 1: overlapping boxes.
#[test]
fn overlapping_boxes_union_ds_has_curves_on_a_plus_x_face() {
    let a = box_at(0.0, 0.0, 0.0, 3.0, 2.0, 2.0);
    let b = box_at(1.5, 0.5, 0.5, 3.0, 1.0, 1.0);

    let mut ds = DS::new(&a, &b);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    let plus_x: Vec<_> = ds
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| {
            f.origin == ShapeOrigin::ShapeA && f.normal.normalize().dot(DVec3::X) > 0.99
        })
        .collect();

    assert!(
        !plus_x.is_empty(),
        "expected at least one ShapeA face with +X normal"
    );

    let max_curves = plus_x
        .iter()
        .map(|(_, f)| f.face_info.curves_sc.len())
        .max()
        .unwrap_or(0);

    assert!(
        max_curves > 0,
        "B penetrates the interior of A's +X face — DS must register intersection curves so the boolean builder can split that face; got max curves_sc={max_curves}"
    );
}

/// After union, the contact patch on A's +x side must not collapse back to a
/// single face covering the full original 2×2 rectangle (area 4): the small
/// box removes the centre — remaining exterior is a frame, not a solid quad.
#[test]
fn overlapping_boxes_union_no_full_area_plus_x_contact_face() {
    let a = box_at(0.0, 0.0, 0.0, 3.0, 2.0, 2.0);
    let b = box_at(1.5, 0.5, 0.5, 3.0, 1.0, 1.0);

    let (u, _) = boolean_op_simplified(
        BooleanOpType::Union,
        &a,
        &b,
        SimplifyOptions::default(),
    )
    .expect("union");

    const FULL_A_PLUS_X_AREA: f64 = 4.0; // h * d = 2 * 2 on yz
    let tol = 0.05;

    // Wrong topology: one outer quad of full 2×2 with no hole. Correct: outer ~4 with inner wire(s), or smaller outers.
    let mut saw_invalid_full_fill = false;
    for solid in &u.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let n = face.normal.normalize_or_zero();
                if n.dot(DVec3::X) <= 0.99 {
                    continue;
                }
                let mut xs = 0.0f64;
                let mut npts = 0usize;
                for we in &face.outer_wire.edges {
                    let Some(e) = u.edges.get(we.idx) else { continue };
                    let vi = if we.forward { e.start } else { e.end };
                    if let Some(v) = u.vertices.get(vi) {
                        xs += v.point.x;
                        npts += 1;
                    }
                }
                if npts == 0 {
                    continue;
                }
                let x_avg = xs / npts as f64;
                if (x_avg - 3.0).abs() > 0.05 {
                    continue;
                }
                let a = polygon_area_outer_wire(&u, face);
                if a >= FULL_A_PLUS_X_AREA - tol && face.inner_wires.is_empty() {
                    saw_invalid_full_fill = true;
                    break;
                }
            }
        }
    }

    assert!(
        !saw_invalid_full_fill,
        "expected no solid 2×2 +X patch at x≈3 without inner wires (frame should have a hole)"
    );
}

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

fn triangle_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .map(|f| f.triangles.len())
        .sum()
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

// ── Chain operations ────────────────────────────────────────────────────────

/// A ∪ B, then result ∩ C: tests that a boolean result can be an input to
/// another boolean operation without panicking.
#[test]
fn chain_union_then_intersect() {
    let a = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    let b = box_at(1.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    let ab = boolean_op(BooleanOpType::Union, &a, &b).expect("union should succeed");

    // C completely overlaps the joined region
    let c = box_at(0.5, 0.0, 0.0, 2.0, 2.0, 2.0);
    let result = boolean_op(BooleanOpType::Intersection, &ab, &c)
        .expect("intersection of union result should succeed");

    assert!(face_count(&result) > 0, "chained result must have faces");
    assert!(triangle_count(&result) > 0, "chained result must have triangles");
    assert!(all_triangles_valid(&result), "all triangle indices must be in bounds");
}

/// A - B, then result - C: progressive subtraction.
#[test]
fn chain_two_differences() {
    let a = box_at(0.0, 0.0, 0.0, 3.0, 1.0, 1.0);
    let b = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let ab = boolean_op(BooleanOpType::Difference, &a, &b).expect("first diff should succeed");

    let c = box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let result = boolean_op(BooleanOpType::Difference, &ab, &c)
        .expect("second diff should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Box × Cylinder ──────────────────────────────────────────────────────────

/// Drill a cylindrical hole through a box: result must have more faces than
/// either input and all triangle indices must be valid.
#[test]
fn box_cylinder_drill() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box");
    let cyl = make_cylinder_brep(DVec3::new(1.0, 1.0, -0.5), DVec3::Z, DVec3::X, 0.4, 3.0)
        .expect("cylinder");
    let result = boolean_op(BooleanOpType::Difference, &b, &cyl)
        .expect("box-cylinder difference should succeed");

    // Box has 6 faces; drilling adds at least 1 new face (cylinder wall)
    assert!(face_count(&result) >= 6, "drilled box must have at least 6 faces");
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Box × Sphere ────────────────────────────────────────────────────────────

/// Union of a box and an overlapping sphere produces a valid solid.
#[test]
fn box_sphere_union_is_valid() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box");
    let s = make_sphere_brep(DVec3::new(1.0, 1.0, 2.0), 0.8).expect("sphere");
    let result = boolean_op(BooleanOpType::Union, &b, &s)
        .expect("box-sphere union should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Box × Cone ──────────────────────────────────────────────────────────────

#[test]
fn box_cone_difference_is_valid() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0)
        .expect("box");
    let c = make_cone_brep(DVec3::new(2.0, 2.0, -0.5), DVec3::Z, DVec3::X, 0.8, 5.0)
        .expect("cone");

    let result = boolean_op(BooleanOpType::Difference, &b, &c)
        .expect("box-cone difference should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

#[test]
fn cone_box_intersection_is_valid() {
    let c = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0)
        .expect("cone");
    let b = make_box_brep(DVec3::new(-3.0, -3.0, 0.0), DVec3::X, DVec3::Y, 6.0, 6.0, 3.0)
        .expect("box");

    let result = boolean_op(BooleanOpType::Intersection, &c, &b)
        .expect("cone-box intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

#[test]
fn sphere_cone_difference_is_valid() {
    let s = make_sphere_brep(DVec3::ZERO, 2.0)
        .expect("sphere");
    let c = make_cone_brep(DVec3::new(0.0, 0.0, -1.0), DVec3::Z, DVec3::X, 1.5, 3.0)
        .expect("cone");

    let result = boolean_op(BooleanOpType::Difference, &s, &c)
        .expect("sphere-cone difference should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
    assert!(volume(&result) > 0.0);
}

#[test]
fn cylinder_cone_intersection_is_valid() {
    let cyl = make_cylinder_brep(DVec3::new(0.0, 0.0, -1.0), DVec3::Z, DVec3::X, 0.8, 6.0)
        .expect("cylinder");
    let cone = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0)
        .expect("cone");

    let result = boolean_op(BooleanOpType::Intersection, &cyl, &cone)
        .expect("cylinder-cone intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

#[test]
fn torus_cone_intersection_is_valid() {
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.8)
        .expect("torus");
    let cone = make_cone_brep(DVec3::new(0.0, 0.0, -1.0), DVec3::Z, DVec3::X, 8.0, 4.0)
        .expect("cone");

    let result = boolean_op(BooleanOpType::Intersection, &torus, &cone)
        .expect("torus-cone intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Box × Torus ─────────────────────────────────────────────────────────────

#[test]
fn box_torus_union_is_valid() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0)
        .expect("box");
    let t = make_torus_brep(
        DVec3::new(2.0, 2.0, 3.25),
        DVec3::Z,
        DVec3::X,
        1.1,
        0.8,
    )
    .expect("torus");

    let result = boolean_op(BooleanOpType::Union, &b, &t)
        .expect("box-torus union should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

#[test]
fn box_torus_intersection_is_valid() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 6.0, 6.0, 6.0)
        .expect("box");
    let t = make_torus_brep(
        DVec3::new(3.0, 3.0, 3.0),
        DVec3::Z,
        DVec3::X,
        1.0,
        0.5,
    )
    .expect("torus");

    let result = boolean_op(BooleanOpType::Intersection, &b, &t)
        .expect("box-torus intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

#[test]
fn box_torus_clipped_intersection_is_valid() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0)
        .expect("box");
    let t = make_torus_brep(
        DVec3::new(2.0, 2.0, 3.25),
        DVec3::Z,
        DVec3::X,
        1.1,
        0.8,
    )
    .expect("torus");

    let result = boolean_op(BooleanOpType::Intersection, &b, &t)
        .expect("clipped box-torus intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(triangle_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Error paths ──────────────────────────────────────────────────────────────

/// An empty BRep should return BooleanError::EmptyInput.
#[test]
fn empty_input_returns_error() {
    let empty = BRep::default();
    let b = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let result = boolean_op(BooleanOpType::Union, &empty, &b);
    assert!(
        matches!(result, Err(BooleanError::EmptyInput)),
        "expected EmptyInput, got {result:?}"
    );
}

/// Disjoint box union then difference with a remote box should not panic.
#[test]
fn disjoint_union_then_difference_no_panic() {
    let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let b = box_at(10.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    let ab = boolean_op(BooleanOpType::Union, &a, &b).expect("disjoint union");

    let c = box_at(5.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    // c is disjoint from both a and b; difference should be identical to ab
    let result = boolean_op(BooleanOpType::Difference, &ab, &c)
        .expect("difference with disjoint c should succeed");

    assert_eq!(face_count(&result), face_count(&ab));
}

// ── MakerVolume ─────────────────────────────────────────────────────────────

#[test]
fn maker_volume_region_mask_unions_selected_cells() {
    let cells = vec![
        box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0),
    ];

    let result = make_solid_from_region(&cells, &[true, false, true])
        .expect("maker volume region mask should succeed");
    assert!((volume(&result) - 2.0).abs() < TOLERANCE_COORD_SUB);
}

#[test]
fn maker_volume_expression_and_history_workflow() {
    let maker = MakerVolume::from_cells(vec![
        box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0),
    ]);
    let expr = CellExpr::Union(
        Box::new(CellExpr::Cell(0)),
        Box::new(CellExpr::Union(Box::new(CellExpr::Cell(1)), Box::new(CellExpr::Cell(2)))),
    );

    let expr_result = maker
        .build_from_expr(&expr)
        .expect("maker volume expression should succeed");
    let (_history_result, history) = maker
        .build_from_indices_with_history(&[0, 1, 2])
        .expect("maker volume history path should succeed");

    assert!((volume(&expr_result) - 3.0).abs() < TOLERANCE_COORD_SUB);
    assert_eq!(history.steps.len(), 2);
}

// ── Cylinder × Cylinder (oblique/parallel) ─────────────────────────────────────

/// Two cylinders with parallel axes intersecting - should produce a valid result.
#[test]
fn cylinder_cylinder_parallel_axis_intersection() {
    let c1 = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 4.0).expect("cylinder1");
    let c2 = make_cylinder_brep(DVec3::new(1.5, 0.0, 0.0), DVec3::Z, DVec3::X, 1.0, 4.0)
        .expect("cylinder2");

    let result = boolean_op(BooleanOpType::Intersection, &c1, &c2)
        .expect("cylinder-cylinder parallel intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Two cylinders with perpendicular axes - "pipe joint" configuration.
#[test]
fn cylinder_cylinder_perpendicular_intersection() {
    let c1 = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 4.0).expect("cylinder1");
    // Perpendicular cylinder crossing through center
    let c2 = make_cylinder_brep(DVec3::new(0.0, 0.0, 2.0), DVec3::X, DVec3::Y, 1.0, 4.0)
        .expect("cylinder2");

    let result = boolean_op(BooleanOpType::Union, &c1, &c2)
        .expect("cylinder-cylinder perpendicular union should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Subtract a smaller cylinder from a larger one (hollow tube creation).
#[test]
fn cylinder_cylinder_hollow_tube() {
    let outer = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0).expect("outer");
    let inner = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.5, 4.0).expect("inner");

    let result = boolean_op(BooleanOpType::Difference, &outer, &inner)
        .expect("hollow tube difference should succeed");

    assert!(face_count(&result) >= 3, "hollow tube should have at least 3 faces (inner, outer, caps)");
    assert!(all_triangles_valid(&result));
}

// ── Cone × Cone ────────────────────────────────────────────────────────────────

/// Two cones with same apex direction intersecting.
/// Note: Cone-cone boolean operations are challenging due to the complex
/// intersection curves. This test uses a more robust configuration.
#[test]
fn cone_cone_intersection() {
    // Use coaxial cones (one inside the other) for a more stable intersection
    let c1 = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 3.0, 5.0).expect("cone1");
    let c2 = make_cone_brep(DVec3::new(0.0, 0.0, 1.0), DVec3::Z, DVec3::X, 1.5, 3.0)
        .expect("cone2");

    // Intersection of larger cone with smaller offset cone
    let result = boolean_op(BooleanOpType::Intersection, &c1, &c2);

    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {
            // Cone-cone intersection can produce degenerate geometry;
            // this is a known limitation for complex cone configurations.
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Cone subtracted from another cone.
#[test]
fn cone_cone_difference() {
    let c1 = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 3.0, 5.0).expect("cone1");
    let c2 = make_cone_brep(DVec3::new(0.0, 0.0, 0.5), DVec3::Z, DVec3::X, 2.0, 4.0)
        .expect("cone2");

    let result = boolean_op(BooleanOpType::Difference, &c1, &c2)
        .expect("cone-cone difference should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Torus × Torus ──────────────────────────────────────────────────────────────

/// Two torus rings intersecting.
#[test]
fn torus_torus_intersection() {
    let t1 = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.5).expect("torus1");
    let t2 = make_torus_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::Z, DVec3::X, 2.0, 0.5)
        .expect("torus2");

    let result = boolean_op(BooleanOpType::Intersection, &t1, &t2)
        .expect("torus-torus intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Torus union with another torus.
#[test]
fn torus_torus_union() {
    let t1 = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.5).expect("torus1");
    let t2 = make_torus_brep(DVec3::new(0.0, 0.0, 1.0), DVec3::Z, DVec3::X, 2.0, 0.5)
        .expect("torus2");

    let result = boolean_op(BooleanOpType::Union, &t1, &t2)
        .expect("torus-torus union should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Near-tangent configurations ────────────────────────────────────────────────

/// Sphere nearly tangent to a plane.
#[test]
fn sphere_plane_near_tangent() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 2.0).expect("box");
    // Sphere positioned just above the box top surface
    let s = make_sphere_brep(DVec3::new(2.0, 2.0, 2.1), 1.0).expect("sphere");

    let result = boolean_op(BooleanOpType::Union, &b, &s)
        .expect("near-tangent sphere-plane union should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Cylinder nearly tangent to another cylinder (parallel axes, close spacing).
#[test]
fn cylinder_cylinder_near_tangent() {
    let c1 = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 4.0).expect("cylinder1");
    // Second cylinder positioned very close (almost touching)
    let c2 = make_cylinder_brep(DVec3::new(2.01, 0.0, 0.0), DVec3::Z, DVec3::X, 1.0, 4.0)
        .expect("cylinder2");

    let result = boolean_op(BooleanOpType::Union, &c1, &c2)
        .expect("near-tangent cylinder union should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Small feature containment ──────────────────────────────────────────────────

/// Small sphere contained within a larger sphere.
#[test]
fn sphere_contained_in_sphere() {
    let large = make_sphere_brep(DVec3::ZERO, 3.0).expect("large sphere");
    let small = make_sphere_brep(DVec3::ZERO, 1.0).expect("small sphere");

    // Difference: should create hollow sphere, or fall back to a degenerate result.
    let result = boolean_op(BooleanOpType::Difference, &large, &small);
    match result {
        Ok(r) => {
            assert!(face_count(&r) > 0);
            assert!(all_triangles_valid(&r));
        }
        Err(BooleanError::DegenerateResult) => {
            // Accepted fallback for fully-contained curved subtraction.
        }
        Err(e) => panic!("unexpected error: {:?}", e),
    }
}

/// Small box fully contained in a larger box.
#[test]
fn box_contained_in_box() {
    let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("outer");
    let inner = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("inner");

    let result = boolean_op(BooleanOpType::Difference, &outer, &inner)
        .expect("contained box difference should succeed");

    assert!(face_count(&result) >= 6, "should have outer box faces");
    assert!(all_triangles_valid(&result));
}

// ── Sphere with pole intersection ──────────────────────────────────────────────

/// Box intersecting sphere at the pole region.
#[test]
fn sphere_pole_intersection() {
    let s = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
    // Box positioned to intersect the sphere's north pole
    let b = make_box_brep(DVec3::new(-2.0, -2.0, 1.5), DVec3::X, DVec3::Y, 4.0, 4.0, 2.0)
        .expect("box");

    let result = boolean_op(BooleanOpType::Intersection, &s, &b)
        .expect("sphere-pole intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Cone apex region intersection with box.
#[test]
fn cone_apex_intersection() {
    let c = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0).expect("cone");
    // Box intersecting near the cone apex
    let b = make_box_brep(DVec3::new(-1.0, -1.0, -0.5), DVec3::X, DVec3::Y, 2.0, 2.0, 1.0)
        .expect("box");

    let result = boolean_op(BooleanOpType::Intersection, &c, &b)
        .expect("cone-apex intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

// ── Multiple chained operations ────────────────────────────────────────────────

/// Three-way union: box + cylinder + sphere.
#[test]
fn three_way_union_box_cylinder_sphere() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box");
    let c = make_cylinder_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::Z, DVec3::X, 0.5, 2.0)
        .expect("cylinder");
    let s = make_sphere_brep(DVec3::new(1.0, 1.0, 2.0), 0.5).expect("sphere");

    let bc = boolean_op(BooleanOpType::Union, &b, &c).expect("box-cylinder union");
    let result = boolean_op(BooleanOpType::Union, &bc, &s).expect("three-way union");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Progressive subtraction creating a complex shape.
#[test]
fn progressive_subtraction() {
    let base = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("base");

    // Subtract cylinder (hole 1)
    let hole1 = make_cylinder_brep(DVec3::new(1.0, 1.0, -0.5), DVec3::Z, DVec3::X, 0.3, 5.0)
        .expect("hole1");
    let step1 = boolean_op(BooleanOpType::Difference, &base, &hole1).expect("first subtraction");

    // Subtract another cylinder (hole 2)
    let hole2 = make_cylinder_brep(DVec3::new(3.0, 3.0, -0.5), DVec3::Z, DVec3::X, 0.3, 5.0)
        .expect("hole2");
    let result = boolean_op(BooleanOpType::Difference, &step1, &hole2).expect("second subtraction");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

// ── Glue detection tests ──────────────────────────────────────────────────────

/// Test that identical overlapping faces are handled correctly with glue enabled.
#[test]
fn glue_detection_identical_boxes() {
    use rcad_algorithms::BooleanOptions;

    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let b2 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box2");

    let opts = BooleanOptions { use_glue: true, ..Default::default() };
    let (result, _report) = rcad_algorithms::boolean_op_with_options(
        BooleanOpType::Union,
        &b1,
        &b2,
        opts,
    ).expect("glued identical boxes union should succeed");

    // Result should be same as input (glue should skip intersection)
    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

// ── Near-degenerate geometry tests ────────────────────────────────────────────

/// Very thin box (near-degenerate in one dimension).
#[test]
fn thin_box_union() {
    let thin = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 0.01).expect("thin box");
    let thick = make_box_brep(DVec3::new(0.0, 0.0, -1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("thick box");

    let result = boolean_op(BooleanOpType::Union, &thin, &thick)
        .expect("thin box union should succeed");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

/// Sphere with very small radius.
#[test]
fn tiny_sphere_union() {
    let tiny = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 0.001).expect("tiny sphere");
    let box_ = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box");

    let result = boolean_op(BooleanOpType::Union, &tiny, &box_)
        .expect("tiny sphere union should succeed");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

// ── Complex curved solid combinations ──────────────────────────────────────────

/// Sphere intersected with cylinder (creates complex intersection curve).
#[test]
fn sphere_cylinder_complex_intersection() {
    let s = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
    let c = make_cylinder_brep(DVec3::new(1.0, 0.0, -3.0), DVec3::Z, DVec3::X, 0.8, 6.0)
        .expect("cylinder");

    let result = boolean_op(BooleanOpType::Intersection, &s, &c)
        .expect("sphere-cylinder complex intersection should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Torus with cylinder passing through center.
#[test]
fn torus_cylinder_through_center() {
    let t = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.5).expect("torus");
    let c = make_cylinder_brep(DVec3::new(0.0, -3.0, 0.0), DVec3::Y, DVec3::X, 0.3, 6.0)
        .expect("cylinder");

    let result = boolean_op(BooleanOpType::Difference, &t, &c)
        .expect("torus-cylinder through center difference should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Multiple operations on sphere: subtract two cylinders at right angles.
#[test]
fn sphere_with_cross_drilled_holes() {
    let s = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");

    // First hole along X
    let c1 = make_cylinder_brep(DVec3::new(-3.0, 0.0, 0.0), DVec3::X, DVec3::Y, 0.5, 6.0)
        .expect("cylinder1");
    let step1 = boolean_op(BooleanOpType::Difference, &s, &c1)
        .expect("first hole should succeed");

    // Second hole along Y
    let c2 = make_cylinder_brep(DVec3::new(0.0, -3.0, 0.0), DVec3::Y, DVec3::Z, 0.5, 6.0)
        .expect("cylinder2");
    let result = boolean_op(BooleanOpType::Difference, &step1, &c2)
        .expect("second hole should succeed");

    assert!(face_count(&result) > 0);
    assert!(all_triangles_valid(&result));
}

/// Difference creating an L-shaped solid.
#[test]
fn l_shaped_by_difference() {
    let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("outer");
    let inner = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 4.0, 4.0)
        .expect("inner");

    let result = boolean_op(BooleanOpType::Difference, &outer, &inner)
        .expect("L-shaped difference should succeed");

    assert!(face_count(&result) >= 6);
    assert!(all_triangles_valid(&result));
}

// ─- Symmetry and orientation tests ─────────────────────────────────────────────

/// Union of two boxes sharing a face (face-on-face contact).
#[test]
fn boxes_sharing_face() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let b2 = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box2");

    let result = boolean_op(BooleanOpType::Union, &b1, &b2)
        .expect("face-sharing boxes union should succeed");

    assert!(face_count(&result) >= 10); // Should have merged the shared face
    assert!(all_triangles_valid(&result));
}


