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
fn check_merge() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).unwrap();
    let mut ds = DS::new(&ba, &bb);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    // B-face[6]: Plane z=0, curves_in=[1]
    // Use RCAD_DEBUG_BUILDER output to see the sub-face boundaries
    
    let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
    let r = builder.build().expect("build");
    
    // Print all result faces with their boundary vertices
    let nf = r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).count();
    println!("Total faces: {}", nf);
    for (fi, face) in r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        let surf = r.geom.face_surface.get(fi).copied().flatten()
            .and_then(|si| r.geom.surfaces.get(si))
            .map(|s| match s { Surface3::Plane(p) => format!("Plane z={}", p.origin.z),
                               Surface3::BSpline(_) => "BSpline".into(), _ => "?".into() })
            .unwrap_or("?".into());
        println!("Face[{}]: {} nverts={}", fi, surf, face.outer_wire.edges.len());
        for we in &face.outer_wire.edges {
            if let Some(e) = r.edges.get(we.idx) {
                let p = r.vertices[e.start].point;
                println!("  V: ({:.6},{:.6},{:.6})", p.x, p.y, p.z);
            }
        }
        // Last vertex from edge endpoint
        if let Some(we) = face.outer_wire.edges.last() {
            if let Some(e) = r.edges.get(we.idx) {
                let p = r.vertices[e.end].point;
                println!("  V: ({:.6},{:.6},{:.6})", p.x, p.y, p.z);
            }
        }
    }
}
