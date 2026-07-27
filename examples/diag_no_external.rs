use glam::DVec3;
use rcad_algorithms::{
    BooleanOpType, BooleanOptions, boolean_op_with_options, boolean_op_with_retry,
    extrude_polygon_solid, total_surface_area,
};

fn main() {
    // Match I3 test case
    let pr1 = extrude_polygon_solid(
        &[
            DVec3::new(0.0, 50.0, 40.0),
            DVec3::new(0.0, 150.0, 40.0),
            DVec3::new(150.0, 150.0, 40.0),
            DVec3::new(150.0, 50.0, 40.0),
        ],
        DVec3::new(0.0, 0.0, -1.0),
        40.0,
    )
    .expect("pr1");
    eprintln!("pr1 SA: {}", total_surface_area(&pr1));

    let pr2 = extrude_polygon_solid(
        &[
            DVec3::new(25.0, 25.0, 50.0),
            DVec3::new(100.0, 25.0, 50.0),
            DVec3::new(100.0, 125.0, 50.0),
            DVec3::new(25.0, 125.0, 50.0),
        ],
        DVec3::new(0.0, 0.0, -1.0),
        30.0,
    )
    .expect("pr2");
    eprintln!("pr2 SA: {}", total_surface_area(&pr2));

    let po1 =
        boolean_op_with_retry(BooleanOpType::Difference, &pr1, &pr2).expect("po1 = pr1 - pr2");
    let po1_sa = total_surface_area(&po1);
    eprintln!("po1 SA: {} (expected ~53000)", po1_sa);

    // Print face surfaces of po1
    let shell = &po1.solids[0].shells[0];
    eprintln!("po1 has {} faces", shell.faces.len());
    for (i, face) in shell.faces.iter().enumerate() {
        if let Some(surf_idx) = po1.geom.face_surface[i] {
            if let Some(surf) = po1.geom.surfaces.get(surf_idx) {
                eprintln!(
                    "  face {}: {:?} at surface idx {}",
                    i,
                    std::mem::discriminant(surf),
                    surf_idx
                );
            }
        } else {
            eprintln!("  face {}: no surface", i);
        }
    }

    let pr3 = extrude_polygon_solid(
        &[
            DVec3::new(50.0, 75.0, 50.0),
            DVec3::new(125.0, 75.0, 50.0),
            DVec3::new(125.0, 175.0, 50.0),
            DVec3::new(50.0, 175.0, 50.0),
        ],
        DVec3::new(0.0, 0.0, -1.0),
        30.0,
    )
    .expect("pr3");
    eprintln!("pr3 SA: {}", total_surface_area(&pr3));

    let result =
        boolean_op_with_retry(BooleanOpType::Difference, &po1, &pr3).expect("result = po1 - pr3");

    // Also try without BVH
    let opts_no_bvh = BooleanOptions {
        use_bvh: false,
        ..Default::default()
    };
    let result_no_bvh = boolean_op_with_options(BooleanOpType::Difference, &po1, &pr3, opts_no_bvh)
        .expect("result = po1 - pr3 (no bvh)");
    let sa_nb = total_surface_area(&result_no_bvh.0);
    eprintln!("result SA (no BVH): {} (expected 52000)", sa_nb);
    eprintln!(
        "Result faces (no BVH): {}",
        result_no_bvh.0.solids[0].shells[0].faces.len()
    );

    let sa = total_surface_area(&result);
    eprintln!("result SA: {} (expected 52000)", sa);
    eprintln!("Result faces: {}", result.solids[0].shells[0].faces.len());
}
