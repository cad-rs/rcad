use rcad_kernel::Surface3;
use rcad_kernel::topods;
use rcad_step::StepReader;
use std::path::Path;

fn count_faces(t: &topods::BRep) -> (usize, usize, usize) {
    let mut nf = 0usize;
    let mut np = 0usize;
    let mut nb = 0usize;
    for ts in &t.tshapes {
        if let topods::TShape::Face(fd) = &**ts {
            nf += 1;
            match fd.surface.as_ref() {
                Some(Surface3::Plane(_)) => np += 1,
                Some(Surface3::BSpline(_)) => nb += 1,
                _ => {}
            }
        }
    }
    (nf, np, nb)
}

#[test]
fn check_all_references() {
    for i in 1..=9 {
        let path = format!(
            "../../../tests/occt/step_output/occt_boolean_bfuse_simple_b{}.step",
            i
        );
        let t = StepReader::read_file(Path::new(&path)).expect("load step");
        let (nf, np, nb) = count_faces(&t);
        println!("[OCCT REF B{}] {}f {}PLANE+{}BS", i, nf, np, nb);
    }
}
