use rcad_kernel::{Surface3, BRep};
use std::path::Path;

fn count_faces(r: &BRep) -> (usize, usize, usize) {
    let nf = r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).count();
    let mut np = 0usize; let mut nb = 0usize;
    for (fi,_) in r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        match r.geom.face_surface.get(fi).copied().flatten().and_then(|si| r.geom.surfaces.get(si)) {
            Some(Surface3::Plane(_)) => np += 1,
            Some(Surface3::BSpline(_)) => nb += 1,
            _ => {}
        }
    }
    (nf, np, nb)
}

#[test]
fn check_all_references() {
    for i in 1..=9 {
        let path = format!("../../../tests/occt/step_output/occt_boolean_bfuse_simple_b{}.step", i);
        let r = rcad_step::load_step(Path::new(&path)).expect("load step");
        let (nf, np, nb) = count_faces(&r);
        println!("[OCCT REF B{}] {}f {}PLANE+{}BS", i, nf, np, nb);
    }
}
