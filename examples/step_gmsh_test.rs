//! Test GMSH-compatible STEP export

use glam::DVec3;
use rcad_modeling::*;
use rcad_step::writer::{ExportSelection, StepWriteOptions, StepWriter, StepProtocol};

fn main() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");

    // Export with default options
    let step_default = StepWriter::write_string(&brep, ExportSelection { selected_faces: &[], selected_edges: &[] });
    std::fs::write("test_default.step", &step_default).expect("write");
    println!("Exported test_default.step (default)");

    // Export with gmsh_strict: true
    let step_gmsh = StepWriter::write_string_with_options(
        &brep,
        ExportSelection { selected_faces: &[], selected_edges: &[] },
        &StepWriteOptions {
            protocol: StepProtocol::Ap214,
            gmsh_strict: true,
            ..Default::default()
        },
    );
    std::fs::write("test_gmsh_strict.step", &step_gmsh).expect("write");
    println!("Exported test_gmsh_strict.step (gmsh_strict)");

    // Compare
    let default_faces = step_default.matches("ADVANCED_FACE").count();
    let gmsh_faces = step_gmsh.matches("ADVANCED_FACE").count();
    println!("Default: {} faces, GMSH strict: {} faces", default_faces, gmsh_faces);

    println!("Default has GEOMETRIC_CURVE_SET: {}", step_default.contains("GEOMETRIC_CURVE_SET"));
    println!("GMSH strict has GEOMETRIC_CURVE_SET: {}", step_gmsh.contains("GEOMETRIC_CURVE_SET"));

    println!("Default has SHAPE_REPRESENTATION_RELATIONSHIP: {}", step_default.contains("SHAPE_REPRESENTATION_RELATIONSHIP"));
    println!("GMSH strict has SHAPE_REPRESENTATION_RELATIONSHIP: {}", step_gmsh.contains("SHAPE_REPRESENTATION_RELATIONSHIP"));
}
