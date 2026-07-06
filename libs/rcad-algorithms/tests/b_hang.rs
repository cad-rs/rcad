use glam::DVec3;
use rcad_algorithms::boolean_op;
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

fn count_faces(r: &BRep) -> (usize, usize, usize) {
    let nf = r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).count();
    let mut np = 0; let mut nb = 0;
    for (fi,_) in r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        match r.geom.face_surface.get(fi).copied().flatten().and_then(|si| r.geom.surfaces.get(si)) {
            Some(Surface3::Plane(_)) => np += 1,
            Some(Surface3::BSpline(_)) => nb += 1, _ => {}
        }
    }
    (nf, np, nb)
}

macro_rules! btest {
    ($name:ident, $label:expr, $b2:expr) => {
        #[test]
        fn $name() {
            let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
            let ba = nurbsconvert_brep(ba);
            let (x,y,z,dx,dy,dz) = $b2;
            let bb = make_box_brep(DVec3::new(x,y,z), DVec3::X, DVec3::Y, dx, dy, dz).unwrap();
            let r = boolean_op(BooleanOpType::Union, &ba, &bb).expect("bfuse");
            let (nf, np, nb) = count_faces(&r);
            println!("[RCAD {}] {}f {}PLANE+{}BS", $label, nf, np, nb);
        }
    };
}

btest!(b1, "B1", (0.0, 0.0, 0.0, 0.5, 1.0, 0.5));
btest!(b2, "B2", (0.0, -0.5, 0.0, 0.5, 0.5, 1.0));
btest!(b3, "B3", (0.0, -0.5, 0.0, 0.5, 1.5, 1.0));
btest!(b4, "B4", (0.0, 0.5, 0.0, 1.0, 1.0, 1.0));
btest!(b5, "B5", (0.0, 0.25, 0.0, 1.0, 0.5, 1.0));
btest!(b6, "B6", (0.0, 0.0, 0.0, 0.5, 0.5, 0.5));
btest!(b7, "B7", (0.0, -0.5, 0.0, 0.5, 0.5, 0.5));
btest!(b8, "B8", (0.0, -0.5, -0.5, 0.5, 0.5, 0.5));
btest!(b9, "B9", (-0.5, -0.5, -0.5, 0.5, 0.5, 0.5));
