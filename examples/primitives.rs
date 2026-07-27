//! Example: Create all 5 primitive solids and export each to STEP.
//!
//! Run: cargo run --example primitives

use glam::DVec3;
use rcad_modeling::*;
use rcad_scene::append_brep;
use rcad_step::writer::{ExportSelection, StepWriter};

fn main() {
    // ── Box ────────────────────────────────────────────────────────────
    let box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");
    write_step(&box_brep, "output_box.step");

    // ── Sphere ─────────────────────────────────────────────────────────
    let sphere = make_sphere_brep(DVec3::new(0.0, 0.0, 0.0), 2.0).expect("sphere");
    write_step(&sphere, "output_sphere.step");

    // ── Cylinder ───────────────────────────────────────────────────────
    let cylinder = make_cylinder_brep(DVec3::ZERO, DVec3::Y, DVec3::X, 1.5, 4.0).expect("cylinder");
    write_step(&cylinder, "output_cylinder.step");

    // ── Cone ───────────────────────────────────────────────────────────
    let cone = make_cone_brep(DVec3::ZERO, DVec3::Y, DVec3::X, 2.0, 0.0, 5.0).expect("cone");
    write_step(&cone, "output_cone.step");

    // ── Torus ──────────────────────────────────────────────────────────
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Y, DVec3::X, 3.0, 1.0).expect("torus");
    write_step(&torus, "output_torus.step");

    // ── All primitives combined in one STEP ────────────────────────────
    let mut all = make_box_brep(
        DVec3::new(-6.0, 0.0, 0.0),
        DVec3::X,
        DVec3::Y,
        2.0,
        3.0,
        4.0,
    )
    .expect("box");

    let sphere2 = make_sphere_brep(DVec3::new(0.0, 1.5, 0.0), 2.0).expect("sphere");
    append_brep(&mut all, sphere2);

    let cyl2 = make_cylinder_brep(DVec3::new(6.0, 0.0, 0.0), DVec3::Y, DVec3::X, 1.5, 4.0)
        .expect("cylinder");
    append_brep(&mut all, cyl2);

    let cone2 = make_cone_brep(
        DVec3::new(12.0, 0.0, 0.0),
        DVec3::Y,
        DVec3::X,
        2.0,
        0.0,
        5.0,
    )
    .expect("cone");
    append_brep(&mut all, cone2);

    let torus2 =
        make_torus_brep(DVec3::new(18.0, 2.0, 0.0), DVec3::Y, DVec3::X, 3.0, 1.0).expect("torus");
    append_brep(&mut all, torus2);

    write_step(&all, "output_all_primitives.step");
    println!("Exported 6 STEP files.");
}

fn write_step(brep: &rcad_kernel::BRep, path: &str) {
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
