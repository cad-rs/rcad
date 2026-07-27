use glam::{DVec2, DVec3};
use rcad_algorithms::{TessellationParams, mesh_brep};
use rcad_kernel::{BRep, Curve2dEval, Surface3, Wire};
use rcad_render::Tessellator;
use rcad_step::StepReader;
use std::collections::{BTreeMap, BTreeSet};

fn round_key3(v: [f64; 3], scale: f64) -> [i64; 3] {
    [
        (v[0] * scale).round() as i64,
        (v[1] * scale).round() as i64,
        (v[2] * scale).round() as i64,
    ]
}

fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
}

fn triangle_unit_normal(a: DVec3, b: DVec3, c: DVec3) -> Option<DVec3> {
    let normal = (b - a).cross(c - a).normalize_or_zero();
    if normal.length_squared() <= 1e-20 {
        None
    } else {
        Some(normal)
    }
}

fn quantize_normal(normal: DVec3) -> [i32; 3] {
    [
        (normal.x * 1000.0).round() as i32,
        (normal.y * 1000.0).round() as i32,
        (normal.z * 1000.0).round() as i32,
    ]
}

fn curve2d_default_range(curve: &rcad_kernel::Curve2d) -> [f64; 2] {
    match curve {
        rcad_kernel::Curve2d::Line(_) => [-10.0, 10.0],
        rcad_kernel::Curve2d::Circle(_) | rcad_kernel::Curve2d::Ellipse(_) => {
            [0.0, 2.0 * std::f64::consts::PI]
        }
        rcad_kernel::Curve2d::CircleInvolute(_) => [0.0, 4.0 * std::f64::consts::PI],
        rcad_kernel::Curve2d::ArchimedeanSpiral(_) | rcad_kernel::Curve2d::LogarithmicSpiral(_) => {
            [0.0, 4.0 * std::f64::consts::PI]
        }
        rcad_kernel::Curve2d::SineWave(_) => [0.0, 1.0],
        rcad_kernel::Curve2d::BSpline(_) | rcad_kernel::Curve2d::Bezier(_) => [0.0, 1.0],
    }
}

fn sample_wire_uv_points(brep: &BRep, wire: &Wire, surface_idx: usize) -> Vec<DVec2> {
    let mut points = Vec::new();

    for we in &wire.edges {
        let Some(edge) = brep.edges.get(we.idx) else {
            continue;
        };
        let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) else {
            continue;
        };
        let Some(pcurve) = pcurves.iter().find(|pc| pc.surface_idx == surface_idx) else {
            continue;
        };
        let Some(curve2d) = brep.geom.curve2ds.get(pcurve.curve2d_idx) else {
            continue;
        };

        let mut range = brep
            .geom
            .curve2d_range
            .get(pcurve.curve2d_idx)
            .and_then(|value| *value)
            .or_else(|| {
                brep.geom
                    .edge_curve_range
                    .get(we.idx)
                    .and_then(|value| *value)
            })
            .unwrap_or_else(|| curve2d_default_range(curve2d));
        if !we.forward {
            range.swap(0, 1);
        }

        let start = if we.forward { edge.start } else { edge.end };
        let end = if we.forward { edge.end } else { edge.start };
        let span = (range[1] - range[0]).abs();
        let segments = match curve2d {
            rcad_kernel::Curve2d::Line(_) => 1,
            rcad_kernel::Curve2d::Circle(_) | rcad_kernel::Curve2d::Ellipse(_) => {
                ((span / (2.0 * std::f64::consts::PI) * 48.0).ceil() as usize).clamp(8, 64)
            }
            _ => 24,
        };

        for step in 0..=segments {
            if !points.is_empty() && step == 0 {
                continue;
            }
            let t = range[0] + (range[1] - range[0]) * (step as f64 / segments as f64);
            points.push(curve2d.point_at(t));
        }

        if segments == 1 {
            if points.is_empty() {
                points.push(curve2d.point_at(range[0]));
            }
            points.push(curve2d.point_at(range[1]));
        }

        let _ = (start, end);
    }

    if points.len() >= 2 && (points[0] - points[points.len() - 1]).length() < 1e-9 {
        points.pop();
    }

    points
}

fn point_in_polygon_2d(point: DVec2, polygon: &[DVec2]) -> bool {
    if polygon.len() < 3 {
        return false;
    }

    let mut inside = false;
    let mut prev = polygon.len() - 1;
    for current in 0..polygon.len() {
        let a = polygon[current];
        let b = polygon[prev];
        let intersects = ((a.y > point.y) != (b.y > point.y))
            && (point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y).abs().max(1e-20) + a.x);
        if intersects {
            inside = !inside;
        }
        prev = current;
    }
    inside
}

fn point_is_inside_face_uv(point: DVec2, outer: &[DVec2], holes: &[Vec<DVec2>]) -> bool {
    point_in_polygon_2d(point, outer) && !holes.iter().any(|hole| point_in_polygon_2d(point, hole))
}

fn project_to_uv(surface: &Surface3, point: DVec3) -> Option<DVec2> {
    match surface {
        Surface3::Plane(plane) => {
            let u = rcad_kernel::any_perpendicular(plane.normal);
            let v = plane.normal.cross(u).normalize_or_zero();
            let delta = point - plane.origin;
            Some(DVec2::new(delta.dot(u), delta.dot(v)))
        }
        Surface3::Cylinder(cylinder) => {
            let axis = cylinder.axis.normalize_or_zero();
            let x_axis = rcad_kernel::any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let radial = point - cylinder.origin - axis * (point - cylinder.origin).dot(axis);
            Some(DVec2::new(
                radial.dot(y_axis).atan2(radial.dot(x_axis)),
                (point - cylinder.origin).dot(axis),
            ))
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis.normalize_or_zero();
            let x_axis = rcad_kernel::any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let radial = point - cone.apex - axis * (point - cone.apex).dot(axis);
            Some(DVec2::new(
                radial.dot(y_axis).atan2(radial.dot(x_axis)),
                (point - cone.apex).dot(axis),
            ))
        }
        _ => None,
    }
}

fn parse_step_arg(args: &[String]) -> Result<String, String> {
    let mut path: Option<&str> = None;
    let mut iter = args.iter().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--step" | "-s" => {
                let Some(value) = iter.next() else {
                    return Err("missing file path after --step".to_string());
                };
                path = Some(value);
            }
            value if value.starts_with('-') => return Err(format!("unknown option: {value}")),
            value => path = Some(value),
        }
    }
    path.map(ToOwned::to_owned).ok_or_else(|| {
        "usage: cargo run -p rcad-examples --example diagnose_step_mesh -- --step <file.step>"
            .to_string()
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = match parse_step_arg(&args) {
        Ok(path) => path,
        Err(message) => {
            eprintln!("{message}");
            std::process::exit(2);
        }
    };

    let step = std::fs::read_to_string(&path).expect("failed to read STEP file");
    let mut brep: BRep = StepReader::parse_string(&step).expect("failed to parse STEP file");
    mesh_brep(&mut brep, &TessellationParams::default());

    let vertex_scale = 1_000_000.0;
    let mut coincident_vertices: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
    for (index, vertex) in brep.vertices.iter().enumerate() {
        coincident_vertices
            .entry(round_key3(
                [vertex.point.x, vertex.point.y, vertex.point.z],
                vertex_scale,
            ))
            .or_default()
            .push(index);
    }
    let repeated_vertex_groups: Vec<_> = coincident_vertices
        .iter()
        .filter(|(_, indices)| indices.len() > 1)
        .collect();

    let mut triangle_count = 0usize;
    let mut degenerate_triangles = 0usize;
    let mut duplicate_triangle_groups: BTreeMap<[[i64; 3]; 3], usize> = BTreeMap::new();
    let mut segment_groups: BTreeMap<[[i64; 3]; 2], usize> = BTreeMap::new();
    let mut face_index = 0usize;
    let mut outside_trim_triangles = 0usize;
    let mut outside_faces = Vec::new();
    let mut vertex_normal_groups: BTreeMap<[i64; 3], BTreeSet<[i32; 3]>> = BTreeMap::new();
    let mut hanging_edges_by_face = Vec::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let face_surface = brep
                    .geom
                    .face_surface
                    .get(face_index)
                    .and_then(|value| *value)
                    .and_then(|surface_idx| brep.geom.surfaces.get(surface_idx).cloned())
                    .zip(
                        brep.geom
                            .face_surface
                            .get(face_index)
                            .and_then(|value| *value),
                    );

                let trim_data = face_surface.and_then(|(surface, surface_idx)| {
                    let outer = sample_wire_uv_points(&brep, &face.outer_wire, surface_idx);
                    if outer.len() < 3 {
                        return None;
                    }
                    let holes: Vec<Vec<DVec2>> = face
                        .inner_wires
                        .iter()
                        .map(|wire| sample_wire_uv_points(&brep, wire, surface_idx))
                        .filter(|poly| poly.len() >= 3)
                        .collect();
                    Some((surface, outer, holes))
                });

                let mut face_outside = 0usize;
                let mut face_edges: BTreeMap<[usize; 2], usize> = BTreeMap::new();
                for &[a_idx, b_idx, c_idx] in &face.triangles {
                    triangle_count += 1;
                    let a = brep.vertices[a_idx].point;
                    let b = brep.vertices[b_idx].point;
                    let c = brep.vertices[c_idx].point;
                    let pa = [a.x, a.y, a.z];
                    let pb = [b.x, b.y, b.z];
                    let pc = [c.x, c.y, c.z];
                    if triangle_area(pa, pb, pc) < 1e-12 {
                        degenerate_triangles += 1;
                    }
                    let mut tri = [
                        round_key3(pa, vertex_scale),
                        round_key3(pb, vertex_scale),
                        round_key3(pc, vertex_scale),
                    ];
                    tri.sort();
                    *duplicate_triangle_groups.entry(tri).or_default() += 1;

                    for (p0, p1) in [(pa, pb), (pb, pc), (pc, pa)] {
                        let mut seg = [round_key3(p0, vertex_scale), round_key3(p1, vertex_scale)];
                        seg.sort();
                        *segment_groups.entry(seg).or_default() += 1;
                    }

                    for mut edge in [[a_idx, b_idx], [b_idx, c_idx], [c_idx, a_idx]] {
                        edge.sort();
                        *face_edges.entry(edge).or_default() += 1;
                    }

                    if let Some(unit_normal) = triangle_unit_normal(a, b, c) {
                        let normal_key = quantize_normal(unit_normal);
                        for (vertex_idx, point) in [(a_idx, a), (b_idx, b), (c_idx, c)] {
                            let _ = vertex_idx;
                            vertex_normal_groups
                                .entry(round_key3([point.x, point.y, point.z], vertex_scale))
                                .or_default()
                                .insert(normal_key);
                        }
                    }

                    if let Some((surface, outer, holes)) = &trim_data {
                        let centroid = (a + b + c) / 3.0;
                        if let Some(uv) = project_to_uv(surface, centroid)
                            && !point_is_inside_face_uv(uv, outer, holes)
                        {
                            outside_trim_triangles += 1;
                            face_outside += 1;
                        }
                    }
                }
                if face_outside > 0 {
                    outside_faces.push((face_index, face_outside, face.triangles.len()));
                }
                let hanging_edges = face_edges.values().filter(|count| **count == 1).count();
                if hanging_edges > face.outer_wire.edges.len() {
                    hanging_edges_by_face.push((face_index, hanging_edges, face.triangles.len()));
                }
                face_index += 1;
            }
        }
    }

    let duplicate_triangles: Vec<_> = duplicate_triangle_groups
        .iter()
        .filter(|(_, count)| **count > 1)
        .collect();
    let suspicious_segments: Vec<_> = segment_groups
        .iter()
        .filter(|(_, count)| **count > 2)
        .collect();

    println!("STEP: {path}");
    println!(
        "vertices={} edges={} solids={}",
        brep.vertices.len(),
        brep.edges.len(),
        brep.solids.len()
    );
    println!("triangles={triangle_count} degenerate_triangles={degenerate_triangles}");
    println!("coincident_vertex_groups={}", repeated_vertex_groups.len());
    for (key, indices) in repeated_vertex_groups.iter().take(10) {
        println!("  coincident vertex {:?}: {:?}", key, indices);
    }
    println!("duplicate_triangle_groups={}", duplicate_triangles.len());
    for (tri, count) in duplicate_triangles.iter().take(10) {
        println!("  duplicate triangle x{count}: {:?}", tri);
    }
    println!(
        "segments_shared_by_more_than_two_triangles={}",
        suspicious_segments.len()
    );
    for (seg, count) in suspicious_segments.iter().take(20) {
        println!("  segment x{count}: {:?}", seg);
    }
    println!("triangles_with_centroid_outside_trim={outside_trim_triangles}");
    for (face_idx, outside, total) in outside_faces.iter().take(20) {
        println!("  face {face_idx}: outside_trim={outside}/{total}");
    }
    let split_normal_vertices: Vec<_> = vertex_normal_groups
        .iter()
        .filter(|(_, normals)| normals.len() > 2)
        .collect();
    println!(
        "coincident_vertices_with_multiple_triangle_normals={}",
        split_normal_vertices.len()
    );
    for (key, normals) in split_normal_vertices.iter().take(20) {
        println!(
            "  vertex {:?}: distinct_triangle_normals={}",
            key,
            normals.len()
        );
    }
    println!(
        "faces_with_extra_hanging_edges={}",
        hanging_edges_by_face.len()
    );
    for (face_idx, hanging_edges, total) in hanging_edges_by_face.iter().take(20) {
        println!("  face {face_idx}: hanging_edges={hanging_edges}, triangles={total}");
    }

    let mut face_triangle_counts = Vec::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                face_triangle_counts.push(face.triangles.len());
            }
        }
    }
    face_triangle_counts.sort_unstable();
    println!("face_triangle_count_histogram:");
    let mut hist = BTreeMap::<usize, usize>::new();
    for count in face_triangle_counts {
        *hist.entry(count).or_default() += 1;
    }
    for (triangles, faces) in hist {
        println!("  {triangles} -> {faces}");
    }

    let unique_edge_keys: BTreeSet<_> = brep
        .edges
        .iter()
        .map(|edge| {
            let mut seg = [
                round_key3(
                    [
                        brep.vertices[edge.start].point.x,
                        brep.vertices[edge.start].point.y,
                        brep.vertices[edge.start].point.z,
                    ],
                    vertex_scale,
                ),
                round_key3(
                    [
                        brep.vertices[edge.end].point.x,
                        brep.vertices[edge.end].point.y,
                        brep.vertices[edge.end].point.z,
                    ],
                    vertex_scale,
                ),
            ];
            seg.sort();
            seg
        })
        .collect();
    println!(
        "topology_edges={} unique_topology_segments={}",
        brep.edges.len(),
        unique_edge_keys.len()
    );

    let render_mesh = Tessellator::tessellate(&brep);
    let zero_normals = render_mesh
        .normals
        .iter()
        .filter(|normal| normal.iter().all(|component| component.abs() < 1e-6))
        .count();
    let nonzero_normals = render_mesh.normals.len().saturating_sub(zero_normals);
    println!(
        "render_mesh: nodes={} triangles={} line_indices={} normals={} nonzero_normals={} zero_normals={}",
        render_mesh.nodes.len(),
        render_mesh.indices.len() / 3,
        render_mesh.line_indices.len() / 2,
        render_mesh.normals.len(),
        nonzero_normals,
        zero_normals
    );
}
