//! Example: Boolean operations (union, intersection, difference) between boxes.
//!
//! Run: cargo run --example boolean_ops

use glam::DVec3;
use rcad_algorithms::{BooleanOpType, boolean_op, geom_populate::populate_box_geom};
use rcad_kernel::BRep;
use rcad_modeling::*;
use rcad_step::writer::{ExportSelection, StepWriter};

fn main() {
    // ── 1. Union of two overlapping boxes ──────────────────────────────
    println!("1. Union of two overlapping boxes");
    let mut a = make_box_at(0.0, 0.0, 0.0, 3.0, 2.0, 2.0);
    let mut b = make_box_at(1.5, 0.5, 0.5, 3.0, 1.0, 1.0);
    populate_box_geom(&mut a);
    populate_box_geom(&mut b);

    let union = boolean_op(BooleanOpType::Union, &a, &b).expect("union");
    write_step(&union, "output_bool_union.step");

    // ── 2. Intersection of two overlapping boxes ───────────────────────
    println!("2. Intersection of two overlapping boxes");
    let intersection = boolean_op(BooleanOpType::Intersection, &a, &b).expect("intersection");
    write_step(&intersection, "output_bool_intersection.step");

    // ── 3. Difference A - B ────────────────────────────────────────────
    println!("3. Difference A - B");
    let difference = boolean_op(BooleanOpType::Difference, &a, &b).expect("difference");
    write_step(&difference, "output_bool_difference.step");

    // ── 4. Difference B - A (asymmetric) ──────────────────────────────
    println!("4. Difference B - A");
    let diff_ba = boolean_op(BooleanOpType::Difference, &b, &a).expect("difference B-A");
    write_step(&diff_ba, "output_bool_difference_ba.step");

    // ── 5. Box with a rectangular hole (contained subtraction) ────────
    println!("5. Box with rectangular slot");
    let mut outer = make_box_at(0.0, 0.0, 0.0, 6.0, 4.0, 4.0);
    let mut slot = make_box_at(2.0, 1.0, -0.5, 2.0, 2.0, 5.0);
    populate_box_geom(&mut outer);
    populate_box_geom(&mut slot);

    let slotted = boolean_op(BooleanOpType::Difference, &outer, &slot).expect("slot");
    write_step(&slotted, "output_bool_slot.step");

    // ── 6. Three-box union (chained) ──────────────────────────────────
    println!("6. Three-box cross (chained union)");
    let mut bx = make_box_at(-0.5, -2.0, -0.5, 1.0, 4.0, 1.0);
    let mut by = make_box_at(-2.0, -0.5, -0.5, 4.0, 1.0, 1.0);
    let mut bz = make_box_at(-0.5, -0.5, -2.0, 1.0, 1.0, 4.0);
    populate_box_geom(&mut bx);
    populate_box_geom(&mut by);
    populate_box_geom(&mut bz);

    let mut cross = boolean_op(BooleanOpType::Union, &bx, &by).expect("cross xy");
    // The result of union doesn't have GeomStore populated for further booleans,
    // so we export the two-arm cross and the third arm separately.
    // For full chaining, populate_box_geom would need to be generalized.
    // Instead, combine visually via append_brep:
    rcad_scene::append_brep(&mut cross, bz);
    write_step(&cross, "output_bool_cross.step");

    println!("Exported 6 boolean operation STEP files.");
}

/// Helper: create an axis-aligned box at (x, y, z) with given dimensions.
fn make_box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
    let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, w, h, d).expect("make_box_brep");
    for v in &mut brep.vertices {
        v.point += DVec3::new(x, y, z);
    }
    brep
}

fn write_step(brep: &BRep, path: &str) {
    let step = StepWriter::write_string(
        brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );
    std::fs::write(path, step).expect("write STEP file");
    println!("  -> {path}");
}
