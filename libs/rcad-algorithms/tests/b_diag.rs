use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::pave_filler::PaveFiller;
use rcad_algorithms::builder::BooleanBuilder;
use rcad_algorithms::BooleanOpType;
use rcad_algorithms::geom_convert::surface_to_bspline;
use rcad_kernel::{Surface3, BRep};
use rcad_modeling::make_box_brep;

fn nurbsconvert(mut brep: BRep) -> BRep {
    let params = rcad_algorithms::geom_convert::ConvertParams::default();
    brep.geom.surfaces = brep.geom.surfaces.into_iter()
        .map(|s| Surface3::BSpline(surface_to_bspline(&s, &params))).collect();
    brep
}

#[test]
fn diag_subface_boundaries() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).unwrap();
    let mut ds = DS::new(&ba, &bb);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    // Manually check split_face for B-face[6] (Plane z=0, curves_in=[1])
    let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
    let fi = 6; // B-face, Plane z=0

    // Create a test that calls split_face...
    // We can't access split_face directly, but build() calls it internally
    let r = builder.build().expect("build");
    
    // Print result faces
    let nf = r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).count();
    println!("Result: {} faces", nf);
    for (fi, face) in r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        let surf = r.geom.face_surface.get(fi).copied().flatten()
            .and_then(|si| r.geom.surfaces.get(si))
            .map(|s| match s { Surface3::Plane(_) => "Plane", Surface3::BSpline(_) => "BSpline", _ => "?" })
            .unwrap_or("?");
        println!("  face[{}]: {} nverts={} nedges={}", fi, surf, face.outer_wire.edges.len(),
            face.outer_wire.edges.len());
        for we in &face.outer_wire.edges {
            if let Some(e) = r.edges.get(we.idx) {
                let p1 = r.vertices[e.start].point;
                let p2 = r.vertices[e.end].point;
                println!("    edge[{}]: ({:.6},{:.6},{:.6}) - ({:.6},{:.6},{:.6})", we.idx,
                    p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
            }
        }
    }
}
