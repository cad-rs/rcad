use glam::DVec3;
use rcad_algorithms::boolean_op;
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
fn b5_detail() {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(0.0, 0.25, 0.0), DVec3::X, DVec3::Y, 1.0, 0.5, 1.0).unwrap();
    let r = boolean_op(BooleanOpType::Union, &ba, &bb).expect("bfuse");

    let nf = r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).count();
    let mut np = 0; let mut nb = 0;
    for (fi, face) in r.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        let sname = r.geom.face_surface.get(fi).copied().flatten()
            .and_then(|si| r.geom.surfaces.get(si))
            .map(|s| match s { Surface3::Plane(_) => "Plane", Surface3::BSpline(_) => "BSpline", _ => "?" })
            .unwrap_or("?");
        if sname == "Plane" { np += 1; } else { nb += 1; }
        let nv = face.outer_wire.edges.len();
        let pts: Vec<DVec3> = face.outer_wire.edges.iter().filter_map(|we| r.edges.get(we.idx)).map(|e| r.vertices[e.start].point).collect();
        let xmin = pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let xmax = pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let ymin = pts.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let ymax = pts.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        let zmin = pts.iter().map(|p| p.z).fold(f64::INFINITY, f64::min);
        let zmax = pts.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max);
        println!("[B5] face[{}]: {} nv={} x=[{:.2},{:.2}] y=[{:.2},{:.2}] z=[{:.2},{:.2}]",
            fi, sname, nv, xmin, xmax, ymin, ymax, zmin, zmax);
    }
    println!("[B5] TOTAL: {}f {}PLANE+{}BS", nf, np, nb);
}
