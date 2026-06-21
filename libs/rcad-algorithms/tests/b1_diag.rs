use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::bopds::ds::Interference;
use rcad_algorithms::pave_filler::PaveFiller;
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
fn b1_diag() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).unwrap();
    let mut ds = DS::new(&ba, &bb);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    eprintln!("B1 A-faces:");
    for fi in 0..ds.a_face_count {
        let f = &ds.faces[fi];
        let pts: Vec<DVec3> = f.boundary_verts.iter().map(|&vi| ds.vertices[vi].point).collect();
        let normal = f.normal;
        let cis: Vec<String> = f.face_info.curves_sc.iter().map(|c| format!("{}", c)).collect();
        eprintln!("  A[{}]: n=({:.0},{:.0},{:.0}) curves=[{}] boundary={:?}",
            fi, normal.x, normal.y, normal.z, cis.join(","), 
            pts.iter().map(|p| format!("({:.1},{:.1},{:.1})", p.x, p.y, p.z)).collect::<Vec<_>>());
    }
    
    // Check all FF curves
    eprintln!("\nFF with curves:");
    for inf in &ds.interferences {
        if let Interference::FaceFace { f1, f2, curves, .. } = inf {
            if !curves.is_empty() {
                eprintln!("  FF({},{}): {}", f1, f2, curves.len());
            }
        }
    }
}
