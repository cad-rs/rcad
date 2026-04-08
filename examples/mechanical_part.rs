//! Example: Build a simple mechanical bracket using primitives + booleans.
//!
//! The part is an L-shaped bracket with mounting holes.
//!
//! Run: cargo run --example mechanical_part

use glam::DVec3;
use rcad_algorithms::{BooleanOpType, boolean_op, geom_populate::populate_box_geom};
use rcad_kernel::BRep;
use rcad_modeling::*;
use rcad_scene::append_brep;
use rcad_step::writer::{ExportSelection, StepWriter};

fn main() {
    // ── Step 1: Create the L-shaped bracket base ──────────────────────
    // Horizontal plate: 10 x 2 x 6
    println!("Building L-bracket...");
    let mut plate_h = make_box_at(0.0, 0.0, 0.0, 10.0, 2.0, 6.0);
    // Vertical plate: 2 x 8 x 6  (rises from the left end)
    let mut plate_v = make_box_at(0.0, 0.0, 0.0, 2.0, 8.0, 6.0);
    populate_box_geom(&mut plate_h);
    populate_box_geom(&mut plate_v);

    let bracket = boolean_op(BooleanOpType::Union, &plate_h, &plate_v).expect("L-bracket union");
    write_step(&bracket, "output_bracket.step");
    println!("  -> output_bracket.step (L-bracket base)");

    // ── Step 2: Cut mounting holes as rectangular slots ────────────────
    // Slot 1 in horizontal plate
    let mut slot1 = make_box_at(5.0, -0.5, 2.0, 1.5, 3.0, 2.0);
    // Slot 2 in horizontal plate
    let mut slot2 = make_box_at(8.0, -0.5, 2.0, 1.5, 3.0, 2.0);
    populate_box_geom(&mut slot1);
    populate_box_geom(&mut slot2);

    // Need to populate geom on the bracket result for further booleans.
    // Since bracket is a complex result, we build it step by step:
    // Cut slot1 from horizontal plate first
    let r1 = boolean_op(BooleanOpType::Difference, &plate_h, &slot1).expect("slot1 cut");
    write_step(&r1, "output_bracket_slot1.step");
    println!("  -> output_bracket_slot1.step (horizontal plate with slot 1)");

    // ── Step 3: Combine bracket + separate primitives for display ─────
    // Assemble: L-bracket + cylinder pins at slot locations
    let mut assembly = bracket;

    // Decorative sphere at the bracket corner
    let sphere = make_sphere_brep(DVec3::new(1.0, 9.0, 3.0), 0.8).expect("sphere");
    append_brep(&mut assembly, sphere);

    // Cylinder pin in slot 1 location
    let pin1 = make_cylinder_brep(DVec3::new(5.75, -0.5, 3.0), DVec3::Y, DVec3::X, 0.5, 3.0)
        .expect("pin1");
    append_brep(&mut assembly, pin1);

    // Cylinder pin in slot 2 location
    let pin2 = make_cylinder_brep(DVec3::new(8.75, -0.5, 3.0), DVec3::Y, DVec3::X, 0.5, 3.0)
        .expect("pin2");
    append_brep(&mut assembly, pin2);

    write_step(&assembly, "output_bracket_assembly.step");
    println!("  -> output_bracket_assembly.step (complete assembly)");

    // ── Step 4: Stacked boxes demo (multi-step boolean) ───────────────
    println!("\nBuilding stacked blocks demo...");
    let mut base = make_box_at(0.0, 0.0, 0.0, 8.0, 2.0, 8.0);
    let mut mid = make_box_at(1.0, 2.0, 1.0, 6.0, 2.0, 6.0);
    let mut top = make_box_at(2.0, 4.0, 2.0, 4.0, 2.0, 4.0);
    populate_box_geom(&mut base);
    populate_box_geom(&mut mid);
    populate_box_geom(&mut top);

    let step1 = boolean_op(BooleanOpType::Union, &base, &mid).expect("stack 1");
    // Append top (can't chain boolean without GeomStore on result)
    let mut pyramid = step1;
    append_brep(&mut pyramid, top);

    // Add a torus crown on top
    let crown =
        make_torus_brep(DVec3::new(4.0, 7.0, 4.0), DVec3::Y, DVec3::X, 1.5, 0.3).expect("crown");
    append_brep(&mut pyramid, crown);

    write_step(&pyramid, "output_pyramid.step");
    println!("  -> output_pyramid.step (stepped pyramid with torus crown)");

    // ── Step 5: Nested boxes difference (hollow box) ──────────────────
    println!("\nBuilding hollow box...");
    let mut outer = make_box_at(0.0, 0.0, 0.0, 6.0, 6.0, 6.0);
    let mut inner = make_box_at(0.5, 0.5, 0.5, 5.0, 5.0, 5.0);
    populate_box_geom(&mut outer);
    populate_box_geom(&mut inner);

    let hollow = boolean_op(BooleanOpType::Difference, &outer, &inner).expect("hollow box");
    write_step(&hollow, "output_hollow_box.step");
    println!("  -> output_hollow_box.step (hollow box)");

    println!("\nAll mechanical part examples exported.");
}

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
}
