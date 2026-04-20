use glam::DVec3;
use rcad_algorithms::{BooleanOpType, SimplifyOptions, boolean_op_simplified, geom_populate::populate_box_geom};
use rcad_kernel::BRep;
use rcad_modeling::*;

fn make_box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
    let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, w, h, d).expect("make_box_brep");
    for v in &mut brep.vertices { v.point += DVec3::new(x, y, z); }
    brep
}

fn main() {
    let mut a = make_box_at(0.0, 0.0, 0.0, 3.0, 2.0, 2.0);
    let mut b = make_box_at(1.5, 0.5, 0.5, 3.0, 1.0, 1.0);
    populate_box_geom(&mut a);
    populate_box_geom(&mut b);

    println!("A: box [0..3] x [0..2] x [0..2]");
    println!("B: box [1.5..4.5] x [0.5..1.5] x [0.5..1.5]");
    println!("Expected A∩B: box [1.5..3] x [0.5..1.5] x [0.5..1.5]  (6 faces, each 4 edges)");
    println!();

    let (result, report) = boolean_op_simplified(
        BooleanOpType::Intersection, &a, &b, SimplifyOptions::default()
    ).expect("intersection");

    println!("Simplify report: merges={} internal_removed={} wires_fixed={} vertices_merged={}",
        report.same_domain_face_merges, report.internal_faces_removed,
        report.wires_fixed, report.vertices_merged);
    println!();
    println!("vertices: {}", result.vertices.len());
    println!("edges: {}", result.edges.len());

    for (si, solid) in result.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            println!("solid[{}] shell[{}]: {} faces", si, shi, shell.faces.len());
            for (fi, face) in shell.faces.iter().enumerate() {
                let n = face.normal;
                let outer_edges = face.outer_wire.edges.len();
                let inner_wires = face.inner_wires.len();

                let mut pts: Vec<DVec3> = Vec::new();
                for we in &face.outer_wire.edges {
                    if let Some(e) = result.edges.get(we.idx) {
                        let vi = if we.forward { e.start } else { e.end };
                        if let Some(v) = result.vertices.get(vi) { pts.push(v.point); }
                    }
                }
                let x_min = pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                let x_max = pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                let y_min = pts.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                let y_max = pts.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                let z_min = pts.iter().map(|p| p.z).fold(f64::INFINITY, f64::min);
                let z_max = pts.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max);

                println!("  face[{:2}] n=({:+.2},{:+.2},{:+.2}) edges={} inner={} x=[{:.3}..{:.3}] y=[{:.3}..{:.3}] z=[{:.3}..{:.3}]",
                    fi, n.x, n.y, n.z, outer_edges, inner_wires,
                    x_min, x_max, y_min, y_max, z_min, z_max);

                if outer_edges != 4 || inner_wires > 0 {
                    println!("         [SUSPICIOUS] outer wire sequence:");
                    for (ei, we) in face.outer_wire.edges.iter().enumerate() {
                        if let Some(e) = result.edges.get(we.idx) {
                            let (sv, ev) = if we.forward { (e.start, e.end) } else { (e.end, e.start) };
                            let sp = result.vertices.get(sv).map(|v| v.point).unwrap_or_default();
                            let ep = result.vertices.get(ev).map(|v| v.point).unwrap_or_default();
                            println!("           [{:2}] ({:.3},{:.3},{:.3}) -> ({:.3},{:.3},{:.3})",
                                ei, sp.x, sp.y, sp.z, ep.x, ep.y, ep.z);
                        }
                    }
                }
            }
        }
    }
}
