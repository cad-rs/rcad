//! Smoke test: OCCT-style STEP interchange (surface-scoped edge curves + assembly scaffolding on full export).
//!
//! Run: `cargo run --example step_occt_export`

use glam::DVec3;
use rcad_modeling::*;
use rcad_step::writer::{ExportSelection, StepWriter};

fn main() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");

    let step = StepWriter::write_string(
        &brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );
    let path = "test_occt_interchange.step";
    std::fs::write(path, &step).expect("write STEP");
    println!(
        "Exported {} (ADVANCED_FACE: {}, SURFACE_CURVE: {})",
        path,
        step.matches("ADVANCED_FACE").count(),
        step.matches("SURFACE_CURVE").count()
    );
}
