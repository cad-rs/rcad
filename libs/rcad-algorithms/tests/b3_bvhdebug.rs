use glam::DVec3;
use rcad_algorithms::bvh::Bvh;
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
fn b3_bvh_pairs() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(0.0, -0.5, 0.0), DVec3::X, DVec3::Y, 0.5, 1.5, 1.0).unwrap();

    let bvh_a = Bvh::build(&ba);
    let bvh_b = Bvh::build(&bb);

    let pairs = Bvh::candidate_pairs(&bvh_a, &bvh_b);

    // Check for specific pairs of interest
    eprintln!("Total candidate pairs: {}", pairs.len());
    eprintln!("bvh_a faces={}, bvh_b faces={}", bvh_a.face_count(), bvh_b.face_count());

    // Check for (0, 3) = b1 face[0] (z=0) with b2 face[3] (y=1.0)
    eprintln!("(0,3) present: {}", pairs.iter().any(|(a,b)| *a==0 && *b==3));
    eprintln!("(0,5) present: {}", pairs.iter().any(|(a,b)| *a==0 && *b==5)); // z=0 vs x=0.5
    eprintln!("(5,3) present: {}", pairs.iter().any(|(a,b)| *a==5 && *b==3)); // x=1 vs y=1.0
    
    // List all pairs for b1 face[0]
    let p0: Vec<usize> = pairs.iter().filter(|(a,_)| *a==0).map(|(_,b)| *b).collect();
    eprintln!("b1[0] pairs: {:?}", p0);
    
    // List all pairs for b1 face[5]
    let p5: Vec<usize> = pairs.iter().filter(|(a,_)| *a==5).map(|(_,b)| *b).collect();
    eprintln!("b1[5] pairs: {:?}", p5);
}
