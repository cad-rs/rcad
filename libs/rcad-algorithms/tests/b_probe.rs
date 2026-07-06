use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::pave_filler::PaveFiller;
use rcad_algorithms::builder::BooleanBuilder;
use rcad_algorithms::BooleanOpType;
use rcad_algorithms::geom_convert::surface_to_bspline;
use rcad_kernel::{topods, Surface3, BRep};
use rcad_modeling::make_box_brep;

fn nurbsconvert_brep(mut brep: BRep) -> BRep {
    let params = rcad_algorithms::geom_convert::ConvertParams::default();
    brep.geom.surfaces = brep.geom.surfaces.into_iter()
        .map(|s| Surface3::BSpline(surface_to_bspline(&s, &params)))
        .collect();
    brep
}

#[test]
fn probe_bspline_boolean() {
    use std::time::Instant;
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert_brep(ba);
    let bb = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).unwrap();

    let mut ds = DS::new(&ba, &bb);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();
    eprintln!("PaveFiller done");

    let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
    eprintln!("Starting build()...");
    let r = builder.build().expect("build");
    let nf = r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).count();
    let mut np = 0usize; let mut nb = 0usize;
    for (fi,_) in r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        match r.geom.face_surface.get(fi).copied().flatten().and_then(|si| r.geom.surfaces.get(si)) {
            Some(Surface3::Plane(_)) => np += 1,
            Some(Surface3::BSpline(_)) => nb += 1,
            _ => {}
        }
    }
    println!("RCAD B1: {}f {}PLANE+{}BS", nf, np, nb);
}
