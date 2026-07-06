//! Debug: sphere-box union — examine built BRep faces and compute per-face geometry
use glam::{DVec2, DVec3};
use rcad_algorithms::geom_populate::populate_box_geom;
use rcad_algorithms::{boolean_op, boolean_op_with_history, BooleanOpType};
use rcad_kernel::geom::Surface3;
use rcad_kernel::properties::face_triangles_pub;
use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_modeling::{make_box_brep, make_sphere_brep};
use std::f64::consts::PI;

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
    eprintln!("\n--- Wire / surface diagnostics ---");
    for (fi, face) in hist_result.solids[0].shells[0].faces.iter().enumerate() {
        let n_we = face.outer_wire.edges.len();
        let sidx = hist_result.geom.face_surface.get(fi).and_then(|o| *o);
        let surf_name = sidx.and_then(|i| hist_result.geom.surfaces.get(i))
            .map(|s| match s {
                rcad_kernel::geom::Surface3::Sphere(_) => "Sphere",
                rcad_kernel::geom::Surface3::Plane(_) => "Plane",
                _ => "Other",
            }).unwrap_or("None");
        let mut sample_pts = {
            use rcad_kernel::topology::Wire;
            // Sample wire polyline manually (same logic as sample_wire_polyline_3d)
            let mut pts = Vec::new();
            for we in &face.outer_wire.edges {
                if let Some(edge) = hist_result.edges.get(we.idx) {
                    let vidx = if we.forward { edge.start } else { edge.end };
                    if let Some(v) = hist_result.vertices.get(vidx) {
                        pts.push(v.point);
                    }
                }
            }
            pts
        };
        // dedup consecutive
        sample_pts.dedup_by(|a, b| (*a - *b).length_squared() < 1e-10);
        eprintln!("  face[{fi}]: surf={surf_name} wire_edges={n_we} sample_pts={} stored_tris={}", sample_pts.len(), face.triangles.len());
    }
    eprintln!();

    // UV polygon area diagnostic for sphere faces
    {
        use std::f64::consts::PI;
        let sphere_s = match hist_result.geom.surfaces.get(
            hist_result.geom.face_surface.get(0).and_then(|o| *o).unwrap_or(usize::MAX)
        ) {
            Some(rcad_kernel::geom::Surface3::Sphere(s)) => *s,
            _ => panic!("first face not sphere"),
        };
        eprintln!("--- UV polygon area diagnostic ---");
        let mut total_uv_area = 0.0_f64;
        for (fi, face) in hist_result.solids[0].shells[0].faces.iter().enumerate() {
            let sidx = hist_result.geom.face_surface.get(fi).and_then(|o| *o);
            let is_sphere = sidx.and_then(|i| hist_result.geom.surfaces.get(i))
                .map(|s| matches!(s, rcad_kernel::geom::Surface3::Sphere(_)))
                .unwrap_or(false);
            if !is_sphere { continue; }
            let mut pts: Vec<DVec3> = Vec::new();
            for we in &face.outer_wire.edges {
                if let Some(edge) = hist_result.edges.get(we.idx) {
                    let vidx = if we.forward { edge.start } else { edge.end };
                    if let Some(v) = hist_result.vertices.get(vidx) {
                        pts.push(v.point);
                    }
                }
            }
            pts.dedup_by(|a, b| (*a - *b).length_squared() < 1e-10);
            let uv_pts: Vec<DVec2> = pts.iter().map(|p| sphere_s.world_to_uv(*p)).collect();
            let mut uv_area = 0.0_f64;
            let n = uv_pts.len();
            for i in 0..n {
                let j = (i + 1) % n;
                uv_area += uv_pts[i].x * uv_pts[j].y - uv_pts[j].x * uv_pts[i].y;
            }
            uv_area = uv_area.abs() / 2.0;
            total_uv_area += uv_area;
            let u_min = uv_pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
            let u_max = uv_pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
            let v_min = uv_pts.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
            let v_max = uv_pts.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
            eprintln!("  face[{fi}]: UV_area={uv_area:.4} n_uv={} u=[{u_min:.4},{u_max:.4}] v=[{v_min:.4},{v_max:.4}]",
                uv_pts.len(), fi=fi, uv_area=uv_area);
        }
        let full_uv_area = 2.0 * PI * PI;
        eprintln!("  TOTAL sphere face UV area: {total_uv_area:.6} (full={full_uv_area:.6}, 7/8={:.6})",
            full_uv_area * 7.0 / 8.0);
        eprintln!();
    }

    eprintln!("--- Per-face from face_triangles (orient_tri + UV ear-cut) ---");
    for (fi, face) in hist_result.solids[0].shells[0].faces.iter().enumerate() {
        let tris = face_triangles_pub(&hist_result, face, fi);
        let n_tris = tris.len();
        let vol: f64 = tris.iter().map(|&[a, b, c]| {
            a.x * (b.y * c.z - c.y * b.z)
            + b.x * (c.y * a.z - a.y * c.z)
            + c.x * (a.y * b.z - b.y * a.z)
        }).sum::<f64>() / 6.0;
        eprintln!("  face[{fi}]: tris={n_tris} vol={vol:.6}");
    }

    // Deep diagnostic: for face[0], trace UV earcut
    if let Some(face0) = hist_result.solids[0].shells[0].faces.get(0) {
        let sidx = hist_result.geom.face_surface.get(0).and_then(|o| *o);
        if let Some(Surface3::Sphere(s)) = sidx.and_then(|i| hist_result.geom.surfaces.get(i)) {
            let mut outer: Vec<DVec3> = Vec::new();
            for we in &face0.outer_wire.edges {
                if let Some(edge) = hist_result.edges.get(we.idx) {
                    let vidx = if we.forward { edge.start } else { edge.end };
                    if let Some(v) = hist_result.vertices.get(vidx) {
                        outer.push(v.point);
                    }
                }
            }
            outer.dedup_by(|a, b| (*a - *b).length_squared() < 1e-10);
            eprintln!("\n--- UV deep diagnostic for face[0]: {} boundary pts ---", outer.len());
            let outer_uv: Vec<DVec2> = outer.iter().map(|p| s.world_to_uv(*p)).collect();
            let u_min = outer_uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
            let u_max = outer_uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
            let v_min = outer_uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
            let v_max = outer_uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
            eprintln!("  UV bbox: u=[{:.6}, {:.6}] v=[{:.6}, {:.6}]", u_min, u_max, v_min, v_max);
            eprintln!("  First 5 UV: {:?}", &outer_uv[..5.min(outer_uv.len())]);
            eprintln!("  Last 5 UV: {:?}", &outer_uv[outer_uv.len().saturating_sub(5)..]);
            // Check UV area via ear-cut
            let mut flat: Vec<f64> = Vec::new();
            for uv in &outer_uv {
                flat.push(uv.x);
                flat.push(uv.y);
            }
            let uv_coords: Vec<[f64; 2]> = flat.chunks_exact(2).map(|c| [c[0], c[1]]).collect();
            let mut indices = Vec::new();
            let mut ear_inst = earcut::Earcut::new();
            ear_inst.earcut(uv_coords, &[] as &[usize], &mut indices);
            eprintln!("  Earcut indices: {} ({} tris)", indices.len(), indices.len() / 3);
            if !indices.is_empty() {
                let all_3d = outer.clone();
                let mut total_vol = 0.0;
                for tri in indices.chunks_exact(3) {
                    let a = all_3d[tri[0]];
                    let b = all_3d[tri[1]];
                    let c = all_3d[tri[2]];
                    total_vol += (a.x * (b.y * c.z - c.y * b.z)
                        + b.x * (c.y * a.z - a.y * c.z)
                        + c.x * (a.y * b.z - b.y * a.z)) / 6.0;
                }
                eprintln!("  Volume from UV earcut: {:.6}", total_vol);
            }
        }
    }
    eprintln!();

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
