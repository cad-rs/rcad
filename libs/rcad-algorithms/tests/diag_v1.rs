// Diagnostic for v1 offset test
use rcad_algorithms::{
    offset_shape, OffsetOptions, JoinType,
    brep_algo::total_volume,
};
use rcad_algorithms::extrude_polygon_solid;
use glam::DVec3;

#[test]
fn diag_v1_offset() {
    let s = extrude_polygon_solid(&[
        DVec3::new(4.0, 0.0, 0.0),
        DVec3::new(10.0, 0.0, 0.0),
        DVec3::new(10.0, 0.0, 10.0),
        DVec3::new(8.0, 0.0, 10.0),
        DVec3::new(8.0, 0.0, 4.0),
        DVec3::new(6.0, 0.0, 4.0),
        DVec3::new(6.0, 0.0, 6.0),
        DVec3::new(4.0, 0.0, 6.0),
        DVec3::new(4.0, 0.0, 0.0),
    ], DVec3::new(0.0, 1.0, 0.0), 5.0).expect("extrude");

    println!("=== SOURCE ===");
    println!("Source volume: {:.6}", total_volume(&s));
    println!("Source vertices: {}", s.vertices.len());

    let shell = &s.solids[0].shells[0];

    // For each vertex, list incident face indices
    println!("\n=== VERTEX INCIDENT FACES ===");
    for vi in 0..s.vertices.len() {
        let pt = s.vertices[vi].point;
        let mut incident: Vec<(usize, DVec3)> = Vec::new();
        for (fi, face) in shell.faces.iter().enumerate() {
            if face.normal.is_nan() { continue; }
            let uses_vert = face.outer_wire.edges.iter().any(|we| {
                let e = &s.edges[we.idx];
                e.start == vi || e.end == vi
            }) || face.inner_wires.iter().any(|wire| {
                wire.edges.iter().any(|we| {
                    let e = &s.edges[we.idx];
                    e.start == vi || e.end == vi
                })
            });
            if uses_vert {
                incident.push((fi, face.normal));
            }
        }
        let incident_info: Vec<String> = incident.iter()
            .map(|(fi, n)| format!("F{} ({:.2},{:.2},{:.2})", fi, n.x, n.y, n.z))
            .collect();
        println!("  V{vi} ({:.4},{:.4},{:.4}) → {}", pt.x, pt.y, pt.z, incident_info.join(", "));
    }

    // Also show which face corresponds to which polygon edge
    println!("\n=== FACE EDGE MAP ===");
    for (fi, face) in shell.faces.iter().enumerate() {
        if face.normal.is_nan() { continue; }
        let edges_info: Vec<String> = face.outer_wire.edges.iter().map(|we| {
            let e = &s.edges[we.idx];
            format!("E{} V{}→V{}", we.idx, e.start, e.end)
        }).collect();
        println!("  F{fi}: normal=({:.2},{:.2},{:.2}) edges=[{}]",
            face.normal.x, face.normal.y, face.normal.z,
            edges_info.join(", "));
    }

    // Run offset
    let opts = OffsetOptions::new(2.0)
        .with_tolerance(1e-7)
        .with_join_type(JoinType::Intersection);
    let result_raw = offset_shape(&s, opts).expect("offset");
    let result = result_raw.brep;

    println!("\n=== OFFSET RESULT ===");
    println!("Result volume: {:.6} (expected: 1116.0)", total_volume(&result));
    println!("Result vertices: {}", result.vertices.len());
    for (vi, v) in result.vertices.iter().enumerate() {
        println!("  V{vi}: ({:.4},{:.4},{:.4})", v.point.x, v.point.y, v.point.z);
    }
}
