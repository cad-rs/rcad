/// Debug test: box-sphere union matching the GMSH example.
/// Box (0,0,0)-(1.6,1.2,1.0) + Sphere (1.0,0.6,0.5) r=0.7 Fuse.
use glam::DVec3;
use rcad_kernel::{face_surface_area, surface_area, BRep, Surface3};
use rcad_algorithms::{boolean_op, boolean_op_simplified, BooleanOpType, SimplifyOptions};

fn print_stats(label: &str, brep: &BRep) {
    let sa = surface_area(brep);
    let n_faces: usize = brep.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    let n_edges: usize = brep.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .map(|f| {
            f.outer_wire.edges.len()
            + f.inner_wires.iter().map(|w| w.edges.len()).sum::<usize>()
        })
        .sum();
    let n_verts = brep.vertices.len();

    println!("[{label}] faces={n_faces} edges={n_edges} verts={n_verts} SA={sa:.10}");

    let mut fi = 0usize;
    for s in &brep.solids {
        for sh in &s.shells {
            for f in &sh.faces {
                let st = brep.geom.face_surface.get(fi)
                    .and_then(|&si| si)
                    .and_then(|si| brep.geom.surfaces.get(si))
                    .map(|s| match s {
                        Surface3::Plane(..) => "Plane",
                        Surface3::Sphere(..) => "Sphere",
                        Surface3::Cylinder(..) => "Cylinder",
                        Surface3::Cone(..) => "Cone",
                        Surface3::Torus(..) => "Torus",
                        _ => "Other",
                    })
                    .unwrap_or("?");
                let fsa = face_surface_area(brep, f, fi);
                let ne = f.outer_wire.edges.len();
                let niw = f.inner_wires.len();
                println!("  F[{fi}] {st} outer_edges={ne} inner_wires={niw} SA={fsa:.10}");
                fi += 1;
            }
        }
    }
}

#[test]
fn debug_box_sphere_union_exact() {
    let a = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.6, 1.2, 1.0).unwrap();
    let b = rcad_modeling::make_sphere_brep(DVec3::new(1.0, 0.6, 0.5), 0.7).unwrap();

    // Raw boolean_op — no post-processing
    let raw = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
    println!("=== Raw boolean_op result ===");
    print_stats("raw", &raw);

    // Every face MUST have an analytic surface (not mesh-only!)
    let all_have_surface = raw.geom.face_surface.iter().all(|s| s.is_some());
    assert!(all_have_surface,
        "All faces must have analytic surfaces — found {} without",
        raw.geom.face_surface.iter().filter(|s| s.is_none()).count());

    // Must have edges (not mesh-only)
    assert!(raw.edges.len() > 0, "Result must have edges");

    // Must have box faces + spherical caps = 6 + N intersecting planes
    // For this case: x=0 (no intersect), 5 others → 11 total
    assert_eq!(raw.solids[0].shells[0].faces.len(), 11,
        "Expected 11 faces (6 box + 5 sphere caps), got {}",
        raw.solids[0].shells[0].faces.len());

    // Verify at least one sphere surface and one plane surface
    let has_sphere = raw.geom.surfaces.iter().any(|s| matches!(s, Surface3::Sphere(_)));
    let has_plane = raw.geom.surfaces.iter().any(|s| matches!(s, Surface3::Plane(_)));
    assert!(has_sphere, "Must have at least one sphere surface");
    assert!(has_plane, "Must have at least one plane surface");

    // Simplified
    let sopts = SimplifyOptions::default();
    let (simplified, report) = boolean_op_simplified(BooleanOpType::Union, &a, &b, sopts).unwrap();
    println!("=== Simplified result ===");
    print_stats("simplified", &simplified);

    // Simplified must have valid result
    let simp_faces: usize = simplified.solids.iter()
        .flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
    assert!(simp_faces > 0, "Simplified must have at least one face");
    let sa = surface_area(&simplified);
    assert!(sa > 0.0, "Simplified SA must be positive, got {sa}");
    assert!(sa > 9.0, "Union SA should be > box SA (9.44), got {sa}");

    println!("=== Report ===");
    println!("  issues_before={} issues_after={}", report.issues_before, report.issues_after);
    println!("  vertices_merged={}", report.vertices_merged);
}
