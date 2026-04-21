//! Example: Union of two cubes and export to STEP.
//!
//! Run:
//!   cargo run --example boolean_union_two_boxes

use glam::DVec3;
use rcad_algorithms::{
    BooleanOpType,
    SimplifyOptions,
    boolean_op,
    boolean_op_simplified,
    geom_populate::populate_box_geom,
};
use rcad_kernel::BRep;
use rcad_modeling::make_box_brep;
use rcad_step::writer::{ExportSelection, StepWriter};

fn main() {
    let mut left = make_box_at(0.0, 0.0, 0.0, 20.0);
    let mut right = make_box_at(12.0, 0.0, 0.0, 12.0);

    // Boolean ops rely on geometric payload on the input boxes.
    populate_box_geom(&mut left);
    populate_box_geom(&mut right);

    let raw_union = boolean_op(BooleanOpType::Union, &left, &right).expect("raw boolean union failed");
    println!(
        "raw union: faces={}, edges={}, vertices={}",
        face_count(&raw_union),
        raw_union.edges.len(),
        raw_union.vertices.len()
    );

    let (union, _) = boolean_op_simplified(
        BooleanOpType::Union,
        &left,
        &right,
        SimplifyOptions::default(),
    )
    .expect("boolean union failed");
    println!(
        "simplified union: faces={}, edges={}, vertices={}",
        face_count(&union),
        union.edges.len(),
        union.vertices.len()
    );

    write_step(&union, "output_union_two_boxes.step");
    println!("done: output_union_two_boxes.step");
}

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum()
}

fn make_box_at(x: f64, y: f64, z: f64, size: f64) -> BRep {
    let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, size, size, size).expect("make box");
    for vertex in &mut brep.vertices {
        vertex.point += DVec3::new(x, y, z);
    }
    brep
}

fn write_step(brep: &BRep, path: &str) {
    let step_text = StepWriter::write_string(
        brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );
    std::fs::write(path, step_text).expect("write step file");
}
