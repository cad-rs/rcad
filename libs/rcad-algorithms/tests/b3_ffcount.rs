use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::bopds::ds::Interference;
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
fn b3_ff_list() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(0.0, -0.5, 0.0), DVec3::X, DVec3::Y, 0.5, 1.5, 1.0).unwrap();
    let mut ds = DS::new(&ba, &bb);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    eprintln!("B3 A-face counts:");
    for fi in 0..ds.a_face_count {
        let f = &ds.faces[fi];
        let bnd = &f.boundary_verts;
        let pts: Vec<DVec3> = bnd.iter().map(|&vi| ds.vertices[vi].point).collect();
        let xr = (pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min), pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max));
        let yr = (pts.iter().map(|p| p.y).fold(f64::INFINITY, f64::min), pts.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max));
        let zr = (pts.iter().map(|p| p.z).fold(f64::INFINITY, f64::min), pts.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max));
        let normal = f.normal;
        eprintln!("  A[{}]: x=[{:.1},{:.1}] y=[{:.1},{:.1}] z=[{:.1},{:.1}] n=({:.0},{:.0},{:.0}) curves={:?}",
            fi, xr.0, xr.1, yr.0, yr.1, zr.0, zr.1,
            normal.x, normal.y, normal.z,
            f.face_info.curves_sc.iter().collect::<Vec<_>>());
    }

    eprintln!("\nB3 FF pairs with curves:");
    for inf in &ds.interferences {
        if let Interference::FaceFace { f1, f2, curves, .. } = inf {
            if !curves.is_empty() {
                eprintln!("  FF({},{}): ncurves={}", f1, f2, curves.len());
            }
        }
    }

    // Also check A-face surfaces
    for fi in 0..ds.a_face_count {
        if let Surface3::BSpline(bsp) = &ds.faces[fi].surface {
            eprintln!("  A[{}] BSpline deg=({},{}) cpts={}x{}",
                fi, bsp.degree_u, bsp.degree_v,
                bsp.control_points.len(),
                if !bsp.control_points.is_empty() { bsp.control_points[0].len() } else { 0 });
        }
    }
}
