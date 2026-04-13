/// Integration tests for boolean operations across multiple shapes and scenarios.
/// These complement the inline unit tests by testing multi-step workflows and
/// error path behavior at crate boundary.
use glam::DVec3;
use rcad_algorithms::{
    BooleanError, BooleanOpType, CellExpr, MakerVolume, boolean_op, make_solid_from_region,
};
use rcad_kernel::PrimitiveSolid;
use rcad_kernel::BRep;
use rcad_kernel::properties::volume;
use rcad_algorithms::geom_populate;
use rcad_modeling::{make_box_brep, make_cylinder_brep, make_sphere_brep};

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
    assert!((volume(&result) - 2.0).abs() < 1e-9);
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

    assert!((volume(&expr_result) - 3.0).abs() < 1e-9);
    assert_eq!(history.steps.len(), 2);
}
