//! Debug: sphere-box union — examine built BRep faces and compute per-face geometry
use glam::DVec3;
use rcad_algorithms::geom_populate::populate_box_geom;
use rcad_algorithms::{boolean_op, boolean_op_with_history, BooleanOpType};
use rcad_kernel::BRep;
use rcad_modeling::{make_box_brep, make_sphere_brep};

fn face_count(brep: &BRep) -> usize {
    brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count()
}

fn vertex(brep: &BRep, idx: usize) -> DVec3 {
    brep.vertices[idx].point
}

fn centroid_of_face_triangles(brep: &BRep, face_idx: usize) -> Option<DVec3> {
    let face = &brep.solids[0].shells[0].faces[face_idx];
    if face.triangles.is_empty() {
        return None;
    }
    let mut sum = DVec3::ZERO;
    let mut count = 0usize;
    for &[i0, i1, i2] in &face.triangles {
        sum += vertex(brep, i0) + vertex(brep, i1) + vertex(brep, i2);
        count += 3;
    }
    Some(sum / count as f64)
}

fn face_volume_contribution(brep: &BRep, face_idx: usize) -> f64 {
    let face = &brep.solids[0].shells[0].faces[face_idx];
    if face.triangles.is_empty() {
        return 0.0;
    }
    let mut vol = 0.0;
    for &[i0, i1, i2] in &face.triangles {
        let v0 = vertex(brep, i0);
        let v1 = vertex(brep, i1);
        let v2 = vertex(brep, i2);
        vol += v0.x * (v1.y * v2.z - v2.y * v1.z)
             + v1.x * (v2.y * v0.z - v0.y * v2.z)
             + v2.x * (v0.y * v1.z - v1.y * v0.z);
    }
    vol / 6.0
}

fn face_surface_area(brep: &BRep, face_idx: usize) -> f64 {
    let face = &brep.solids[0].shells[0].faces[face_idx];
    if face.triangles.is_empty() {
        return 0.0;
    }
    let mut area = 0.0;
    for &[i0, i1, i2] in &face.triangles {
        let e1 = vertex(brep, i1) - vertex(brep, i0);
        let e2 = vertex(brep, i2) - vertex(brep, i0);
        area += e1.cross(e2).length();
    }
    area * 0.5
}

#[test]
fn debug_sphere_classification() {
    let s = make_sphere_brep(DVec3::ZERO, 1.0).expect("sphere");
    let mut b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    populate_box_geom(&mut b);

    let sphere_sa = rcad_algorithms::total_surface_area(&s);
    let sphere_vol = rcad_algorithms::total_volume(&s);
    let box_sa = rcad_algorithms::total_surface_area(&b);
    let box_vol = rcad_algorithms::total_volume(&b);
    eprintln!("sphere: SA={sphere_sa:.6} vol={sphere_vol:.6} faces=1");
    eprintln!("box: SA={box_sa:.6} vol={box_vol:.6} faces=6");
    eprintln!();

    // Compare raw builder vs fuse
    let (hist_result, _history) = boolean_op_with_history(BooleanOpType::Union, &s, &b)
        .expect("boolean_op_with_history");
    let nf = face_count(&hist_result);
    let hist_sa = rcad_algorithms::total_surface_area(&hist_result);
    let hist_vol = rcad_algorithms::total_volume(&hist_result);
    eprintln!("boolean_op_with_history: SA={hist_sa:.6} vol={hist_vol:.6} faces={nf}");

    // Per-face analysis
    for (fi, face) in hist_result.solids[0].shells[0].faces.iter().enumerate() {
        let n_tris = face.triangles.len();
        let face_area = face_surface_area(&hist_result, fi);
        let face_vol = face_volume_contribution(&hist_result, fi);
        let centroid = centroid_of_face_triangles(&hist_result, fi);
        let c_str = centroid.map(|c| format!("({:.4}, {:.4}, {:.4})", c.x, c.y, c.z))
            .unwrap_or_else(|| "N/A".into());
        let dist_origin = centroid.map(|c| c.length()).unwrap_or(0.0);
        let normal_len = face.normal.length();
        eprintln!("  face[{fi}]: tris={n_tris} area={face_area:.6} vol_contrib={face_vol:.6} centroid={c_str} dist_from_origin={dist_origin:.4} normal_len={normal_len:.6}");
    }

    let fused = boolean_op(BooleanOpType::Union, &s, &b).expect("fuse");
    let nf_fused = face_count(&fused);
    let fused_sa = rcad_algorithms::total_surface_area(&fused);
    let fused_vol = rcad_algorithms::total_volume(&fused);
    eprintln!("\nboolean_op (fuse): SA={fused_sa:.6} vol={fused_vol:.6} faces={nf_fused}");

    for (fi, face) in fused.solids[0].shells[0].faces.iter().enumerate() {
        let n_tris = face.triangles.len();
        let face_area = face_surface_area(&fused, fi);
        let face_vol = face_volume_contribution(&fused, fi);
        let centroid = centroid_of_face_triangles(&fused, fi);
        let c_str = centroid.map(|c| format!("({:.4}, {:.4}, {:.4})", c.x, c.y, c.z))
            .unwrap_or_else(|| "N/A".into());
        let dist_origin = centroid.map(|c| c.length()).unwrap_or(0.0);
        eprintln!("  face[{fi}]: tris={n_tris} area={face_area:.6} vol_contrib={face_vol:.6} centroid={c_str} dist_origin={dist_origin:.4}");
    }

    eprintln!();
    eprintln!("Expected (A1) SA ~14.6394, vol ~4.665");
}
