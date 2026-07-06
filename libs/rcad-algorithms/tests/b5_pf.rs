use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::bopds::ds::{Interference, ShapeOrigin};
use rcad_algorithms::pave_filler::PaveFiller;
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
fn b5_pf_check() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(0.0, 0.25, 0.0), DVec3::X, DVec3::Y, 1.0, 0.5, 1.0).unwrap();
    let mut ds = DS::new(&ba, &bb);
    
    // Check a_faces and b_faces order
    let a_faces: Vec<usize> = (0..ds.a_face_count).collect();
    let b_faces: Vec<usize> = (ds.a_face_count..ds.faces.len()).collect();
    eprintln!("B5: {} A-faces, {} B-faces", a_faces.len(), b_faces.len());
    
    // List all faces with their surfaces and source_face_idx
    for fi in 0..ds.faces.len() {
        let f = &ds.faces[fi];
        let surf_name = match &f.surface {
            Surface3::Plane(p) => format!("Plane n=({:.1},{:.1},{:.1})", p.normal.x, p.normal.y, p.normal.z),
            Surface3::BSpline(bsp) => "BSpline".to_string(),
            _ => "?".to_string(),
        };
        eprintln!("  DS face[{}]: {} src={}", fi, surf_name, f.source_face_idx);
    }
    
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();
    
    // Check ALL FaceFace pairs that were processed
    eprintln!("\nFaceFace interferences after PaveFiller:");
    for inf in &ds.interferences {
        match inf {
            Interference::FaceFace { f1, f2, curves, .. } => {
                eprintln!("  FF({},{}): ncurves={}", f1, f2, curves.len());
            }
            _ => {}
        }
    }
    
    // Check curves_sc for each A-face
    for fi in 0..ds.a_face_count {
        eprintln!("  A-face[{}] curves_sc={:?}", fi, ds.faces[fi].face_info.curves_sc.iter().collect::<Vec<_>>());
    }
}
