use glam::DVec3;
use rcad_algorithms::boolean_op;
use rcad_algorithms::BooleanOpType;
use rcad_algorithms::geom_convert::surface_to_bspline;
use rcad_kernel::{topods, Surface3, BRep};
use rcad_modeling::make_box_brep;

fn nurbsconvert(mut brep: BRep) -> BRep {
    let params = rcad_algorithms::geom_convert::ConvertParams::default();
    brep.geom.surfaces = brep.geom.surfaces.into_iter()
        .map(|s| Surface3::BSpline(surface_to_bspline(&s, &params))).collect();
    brep
}

fn run_case(label: &str, b2: (f64,f64,f64,f64,f64,f64)) {
    let ba = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let ba = nurbsconvert(ba);
    let bb = make_box_brep(DVec3::new(b2.0,b2.1,b2.2), DVec3::X, DVec3::Y, b2.3,b2.4,b2.5).unwrap();
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
        println!("[{}] face[{}]: {} nv={} x=[{:.1},{:.1}] y=[{:.1},{:.1}] z=[{:.1},{:.1}]",
            label, fi, sname, nv, xmin, xmax, ymin, ymax, zmin, zmax);
    }
    println!("[{}] TOTAL: {}f {}PLANE+{}BS (OCCT={})", label, nf, np, nb, 
        match label { "B3" => "14f", "B4" => "14f", "B5" => "14f", _ => "?" });
}

#[test]
fn b3_detail() { run_case("B3", (0.0, -0.5, 0.0, 0.5, 1.5, 1.0)); }
#[test]
fn b4_detail() { run_case("B4", (0.0, 0.5, 0.0, 1.0, 1.0, 1.0)); }
