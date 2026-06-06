use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::bopds::ds::{Interference, ShapeOrigin};
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

fn check_grid(label: &str, b2: (f64,f64,f64,f64,f64,f64)) {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(b2.0,b2.1,b2.2), DVec3::X, DVec3::Y, b2.3,b2.4,b2.5).unwrap();
    let mut ds = DS::new(&ba, &bb);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();
    
    eprintln!("=== {} ===", label);
    for fi in 0..ds.a_face_count {
        let f = &ds.faces[fi];
        let nv = f.boundary_verts.len();
        let cis: Vec<String> = f.face_info.curves_in.iter().map(|c| format!("{}", c)).collect();
        eprintln!("  A-face[{}] src={} nv={} curves=[{}]", fi, f.source_face_idx, nv, cis.join(","));
    }
}

#[test]
fn b3b4_ff_check() {
    // B3: b2=(0, -0.5, 0, 0.5, 1.5, 1)
    // B4: b2=(0, 0.5, 0, 1, 1, 1)
    check_grid("B3", (0.0, -0.5, 0.0, 0.5, 1.5, 1.0));
    check_grid("B4", (0.0, 0.5, 0.0, 1.0, 1.0, 1.0));
}
