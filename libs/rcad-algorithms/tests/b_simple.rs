use glam::DVec3;
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
fn bspline_plane_boolean() {
    use rcad_algorithms::boolean_op;
    use rcad_algorithms::BooleanOpType;
    use std::time::Instant;
    let t0 = Instant::now();
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert_brep(ba);
    let bb = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).unwrap();
    let t1 = Instant::now();
    let r = boolean_op(BooleanOpType::Union, &ba, &bb).expect("bfuse");
    let t2 = Instant::now();
    println!("bspline_plane_boolean: {} faces, prep={:.3}s boolean_op={:.3}s",
        r.solids[0].shells[0].faces.len(),
        (t1-t0).as_secs_f64(), (t2-t1).as_secs_f64());
}
