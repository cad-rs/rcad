use glam::DVec3;
use rcad_algorithms::boolean_op;
use rcad_algorithms::BooleanOpType;
use rcad_algorithms::geom_convert::surface_to_bspline;
use rcad_kernel::{topods, Surface3, BRep};
use rcad_modeling::make_box_brep;

fn nurbsconvert(mut brep: BRep) -> BRep {
    let params = rcad_algorithms::geom_convert::ConvertParams::default();
    brep.geom.surfaces = brep.geom.surfaces.into_iter()
        .map(|s| Surface3::BSpline(surface_to_bspline(&s, &params))).collect();
    brep
}

#[test]
fn b1_fresh() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).unwrap();
    let r = boolean_op(BooleanOpType::Union, &ba, &bb).expect("bfuse");
    let nf = r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).count();
    let mut np = 0; let mut nb = 0;
    for (fi, _) in r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        match r.geom.face_surface.get(fi).copied().flatten().and_then(|si| r.geom.surfaces.get(si)) {
            Some(Surface3::Plane(_)) => np += 1,
            Some(Surface3::BSpline(_)) => nb += 1,
            _ => {}
        }
    }
    println!("B1: {}f {}PLANE+{}BS", nf, np, nb);
}
