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
fn b5_ff_check() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(0.0, 0.25, 0.0), DVec3::X, DVec3::Y, 1.0, 0.5, 1.0).unwrap();
    let mut ds = DS::new(&ba, &bb);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    eprintln!("A-face[0] curves_in: {:?}", ds.faces[0].face_info.curves_in.iter().collect::<Vec<_>>());
    eprintln!("A-face[4] curves_in: {:?}", ds.faces[4].face_info.curves_in.iter().collect::<Vec<_>>());

    for inf in &ds.interferences {
        match inf {
            Interference::FaceFace { f1, f2, curves, .. } => {
                eprintln!("FF({},{}): ncurves={}", f1, f2, curves.len());
            }
            _ => {}
        }
    }
}
