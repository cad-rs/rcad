use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::bopds::ds::Interference;
use rcad_algorithms::pave_filler::PaveFiller;
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
fn b5_bvh_check() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(0.0, 0.25, 0.0), DVec3::X, DVec3::Y, 1.0, 0.5, 1.0).unwrap();
    
    // Build BVH for both BReps and check candidate pairs
    let bvh_a = Bvh::build(&ba);
    let bvh_b = Bvh::build(&bb);
    
    // Get all pair candidates
    let pairs = Bvh::candidate_pairs(&bvh_a, &bvh_b);
    eprintln!("BVH candidates: {} total", pairs.len());
    
    let has_4_2 = pairs.iter().any(|(a, b)| *a == 4 && *b == 2);
    let has_5_2 = pairs.iter().any(|(a, b)| *a == 5 && *b == 2);
    let has_4_0 = pairs.iter().any(|(a, b)| *a == 4 && *b == 0);
    eprintln!("  (4,2)={} (5,2)={} (4,0)={}", has_4_2, has_5_2, has_4_0);
    
    // Also check brute-force: manually compute
    let bvh_a_faces = bvh_a.face_count();
    let bvh_b_faces = bvh_b.face_count();
    eprintln!("  bvh_a faces={}, bvh_b faces={}", bvh_a_faces, bvh_b_faces);
}
