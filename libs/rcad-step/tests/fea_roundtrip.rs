/// STEP write/parse round-trip checks **specific to downstream FEA meshing**:
/// minimum vertex budgets after translation (STEP tessellation may add sample vertices)
/// and loose bounding checks on analytic primitives.
///
/// Generic topology, coordinates, headers, and error paths are covered in `occt_alignment.rs`.
use glam::DVec3;
use rcad_modeling::{make_box_brep, make_cylinder_brep, make_sphere_brep};
use rcad_step::{ExportSelection, StepReader, StepWriter};

fn all_faces_selection() -> ExportSelection<'static> {
    ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    }
}

fn face_count(brep: &rcad_kernel::BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn vertex_count(brep: &rcad_kernel::BRep) -> usize {
    brep.vertices.len()
}

#[test]
fn box_fea_vertex_budget_after_round_trip() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 20.0, 5.0).expect("box");
    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    assert_eq!(face_count(&parsed), 6, "box should still present six faces");
    // Triangulation sample nodes may appear as extra `BRep.vertices`; mesh pipelines need ≥ 8 corners.
    assert!(
        vertex_count(&parsed) >= 8,
        "box should have at least eight vertices after STEP round-trip, got {}",
        vertex_count(&parsed)
    );
}

#[test]
fn cylinder_fea_vertex_budget_after_round_trip() {
    let brep =
        make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 5.0, 15.0).expect("cylinder");
    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    assert_eq!(face_count(&parsed), 3, "cylinder caps + lateral");
    assert!(
        vertex_count(&parsed) >= 2,
        "cylinder should retain enough vertices for seam/cap meshing, got {}",
        vertex_count(&parsed)
    );
}

#[test]
fn sphere_fea_tessellated_bbox_after_round_trip() {
    let radius = 7.5;
    let brep = make_sphere_brep(DVec3::ZERO, radius).expect("sphere");
    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    assert!(face_count(&parsed) >= 1, "sphere should have at least one face");

    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    let diameter = 2.0 * radius;
    let max_dim = (max.x - min.x).max(max.y - min.y).max(max.z - min.z);
    assert!(
        max_dim >= diameter * 0.8,
        "tessellated sphere bbox should reach most of analytic diameter (mesh/STEP tolerance)"
    );
}
