//! Integration test: OCCT documentation appendix ASCII BREP (topology V1 box).

use rcad_step::OcctBrepReader;

fn face_count(brep: &rcad_kernel::BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

#[test]
fn parses_occt_manual_appendix_box() {
    let src = include_str!("fixtures/occt_doc_box.brep");
    let brep = OcctBrepReader::parse_string(src).expect("parse OCCT doc appendix box BREP");

    assert_eq!(brep.solids.len(), 1, "single solid");
    assert_eq!(brep.solids[0].shells.len(), 1, "single shell");
    assert_eq!(face_count(&brep), 6, "box has six faces");

    assert!(
        brep.vertices.len() >= 8,
        "at least eight corner vertices; triangulation may add more"
    );
    assert!(brep.edges.len() >= 12, "at least twelve edges");

    let tri_total: usize = brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .map(|f| f.triangles.len())
        .sum();
    assert_eq!(tri_total, 12, "six faces × two triangles each from fixture");
}
