use rcad_kernel::{BRep, Curve3, GeomStore, Surface3};
use rcad_kernel::{Edge, Face, Shell, Solid, Vertex, Wire};
use std::collections::{BTreeSet, HashMap};
use std::path::Path;

#[derive(Debug, Clone)]
struct AdvancedFaceRecord {
    bounds: Vec<u64>,
    surface: Option<u64>,
}

#[derive(Debug, Clone)]
struct ParsedStep {
    cartesian_points: HashMap<u64, [f64; 3]>,
    directions: HashMap<u64, [f64; 3]>,
    vectors: HashMap<u64, u64>,
    axis2_placements: HashMap<u64, (u64, u64)>,
    lines: HashMap<u64, (u64, u64)>,
    circles: HashMap<u64, (u64, f64)>,
    planes: HashMap<u64, u64>,
    cylindrical_surfaces: HashMap<u64, (u64, f64)>,
    spherical_surfaces: HashMap<u64, (u64, f64)>,
    conical_surfaces: HashMap<u64, (u64, f64)>,
    toroidal_surfaces: HashMap<u64, (u64, f64, f64)>,
    vertex_points: HashMap<u64, u64>,
    edge_curves: HashMap<u64, (u64, u64, Option<u64>)>,
    oriented_edges: HashMap<u64, (u64, bool)>,
    edge_loops: HashMap<u64, Vec<u64>>,
    face_bounds: HashMap<u64, u64>,
    advanced_faces: HashMap<u64, AdvancedFaceRecord>,
    closed_shells: HashMap<u64, Vec<u64>>,
    manifold_solids: Vec<u64>,
}

impl ParsedStep {
    fn new() -> Self {
        Self {
            cartesian_points: HashMap::new(),
            directions: HashMap::new(),
            vectors: HashMap::new(),
            axis2_placements: HashMap::new(),
            lines: HashMap::new(),
            circles: HashMap::new(),
            planes: HashMap::new(),
            cylindrical_surfaces: HashMap::new(),
            spherical_surfaces: HashMap::new(),
            conical_surfaces: HashMap::new(),
            toroidal_surfaces: HashMap::new(),
            vertex_points: HashMap::new(),
            edge_curves: HashMap::new(),
            oriented_edges: HashMap::new(),
            edge_loops: HashMap::new(),
            face_bounds: HashMap::new(),
            advanced_faces: HashMap::new(),
            closed_shells: HashMap::new(),
            manifold_solids: Vec::new(),
        }
    }
}

pub struct StepReader;

impl StepReader {
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, String> {
        let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        Self::parse_string(&content)
    }

    pub fn parse_string(content: &str) -> Result<BRep, String> {
        if !content.contains("ISO-10303-21") {
            return Err("Invalid STEP file format".to_string());
        }

        let entities = parse_entities(content)?;
        build_brep_from_parsed(&entities)
    }
}

fn parse_entities(content: &str) -> Result<ParsedStep, String> {
    let data = extract_data_section(content)?;
    let records = split_records(data);
    let mut parsed = ParsedStep::new();

    for record in records {
        let Some((id, body)) = parse_entity_record(&record)? else {
            continue;
        };
        if let Some((entity, args)) = parse_entity_body(body) {
            match entity {
                "CARTESIAN_POINT" => {
                    if let Some(coords) = parse_cartesian_point(args) {
                        parsed.cartesian_points.insert(id, coords);
                    }
                }
                "DIRECTION" => {
                    if let Some(coords) = parse_cartesian_point(args) {
                        parsed.directions.insert(id, coords);
                    }
                }
                "VECTOR" => {
                    if let Some(dir_ref) = parse_single_ref_after_name(args) {
                        parsed.vectors.insert(id, dir_ref);
                    }
                }
                "AXIS2_PLACEMENT_3D" => {
                    if let Some((origin, axis)) = parse_axis2_placement(args) {
                        parsed.axis2_placements.insert(id, (origin, axis));
                    }
                }
                "LINE" => {
                    if let Some((origin, vector_ref)) = parse_curve_basis(args) {
                        parsed.lines.insert(id, (origin, vector_ref));
                    }
                }
                "CIRCLE" => {
                    if let Some((placement, radius)) = parse_placement_radius(args) {
                        parsed.circles.insert(id, (placement, radius));
                    }
                }
                "PLANE" => {
                    if let Some(placement) = parse_single_ref_after_name(args) {
                        parsed.planes.insert(id, placement);
                    }
                }
                "CYLINDRICAL_SURFACE" => {
                    if let Some((placement, radius)) = parse_placement_radius(args) {
                        parsed.cylindrical_surfaces.insert(id, (placement, radius));
                    }
                }
                "SPHERICAL_SURFACE" => {
                    if let Some((placement, radius)) = parse_placement_radius(args) {
                        parsed.spherical_surfaces.insert(id, (placement, radius));
                    }
                }
                "CONICAL_SURFACE" => {
                    if let Some((placement, radius, half_angle_rad)) = parse_conical_surface(args) {
                        parsed.conical_surfaces.insert(id, (placement, half_angle_rad));
                        // Keep radius parsed for forward compatibility even if not used yet.
                        let _ = radius;
                    }
                }
                "TOROIDAL_SURFACE" => {
                    if let Some((placement, major, minor)) = parse_toroidal_surface(args) {
                        parsed.toroidal_surfaces.insert(id, (placement, major, minor));
                    }
                }
                "VERTEX_POINT" => {
                    if let Some(point_ref) = parse_single_ref_after_name(args) {
                        parsed.vertex_points.insert(id, point_ref);
                    }
                }
                "EDGE_CURVE" => {
                    if let Some((start, end, curve_ref)) = parse_edge_curve_vertices(args) {
                        parsed.edge_curves.insert(id, (start, end, curve_ref));
                    }
                }
                "ORIENTED_EDGE" => {
                    if let Some((edge_ref, orientation)) = parse_oriented_edge(args) {
                        parsed.oriented_edges.insert(id, (edge_ref, orientation));
                    }
                }
                "EDGE_LOOP" => {
                    if let Some(items) = parse_ref_list_after_name(args) {
                        parsed.edge_loops.insert(id, items);
                    }
                }
                "FACE_BOUND" | "FACE_OUTER_BOUND" => {
                    if let Some(loop_ref) = parse_single_ref_after_name(args) {
                        parsed.face_bounds.insert(id, loop_ref);
                    }
                }
                "ADVANCED_FACE" => {
                    if let Some((bounds, surface)) = parse_advanced_face(args) {
                        parsed.advanced_faces.insert(id, AdvancedFaceRecord { bounds, surface });
                    }
                }
                "CLOSED_SHELL" => {
                    if let Some(face_refs) = parse_ref_list_after_name(args) {
                        parsed.closed_shells.insert(id, face_refs);
                    }
                }
                "MANIFOLD_SOLID_BREP" => {
                    if let Some(shell_ref) = parse_single_ref_after_name(args) {
                        parsed.manifold_solids.push(shell_ref);
                    }
                }
                _ => {}
            }
        }
    }

    Ok(parsed)
}

fn build_brep_from_parsed(parsed: &ParsedStep) -> Result<BRep, String> {
    let shell_face_sets = collect_shell_faces(parsed);
    if shell_face_sets.is_empty() {
        if let Some(brep) = brep_from_points_bbox(parsed) {
            return Ok(brep);
        }
        return Err("STEP parse produced no shell faces".to_string());
    }

    let used_vertex_ids = collect_used_vertices(parsed, &shell_face_sets)?;
    if used_vertex_ids.is_empty() {
        return Err("STEP parse produced no vertices".to_string());
    }

    let mut vertex_ids: Vec<u64> = used_vertex_ids.into_iter().collect();
    vertex_ids.sort_unstable();

    let mut vertex_index_by_id: HashMap<u64, usize> = HashMap::new();
    let mut vertices = Vec::with_capacity(vertex_ids.len());
    for (idx, vertex_id) in vertex_ids.iter().enumerate() {
        let point_id = *parsed
            .vertex_points
            .get(vertex_id)
            .ok_or_else(|| format!("Missing VERTEX_POINT reference for #{}", vertex_id))?;
        let point = *parsed
            .cartesian_points
            .get(&point_id)
            .ok_or_else(|| format!("Missing CARTESIAN_POINT #{}", point_id))?;
        vertices.push(Vertex {
            point: glam::DVec3::new(point[0], point[1], point[2]),
        });
        vertex_index_by_id.insert(*vertex_id, idx);
    }

    let mut edges: Vec<Edge> = Vec::new();
    let mut edge_index_by_curve: HashMap<u64, usize> = HashMap::new();
    let mut curve_store_index_by_step: HashMap<u64, usize> = HashMap::new();
    let mut surface_store_index_by_step: HashMap<u64, usize> = HashMap::new();
    let mut solids: Vec<Solid> = Vec::new();
    let mut geom = GeomStore::default();

    for shell_faces in shell_face_sets {
        let mut faces: Vec<Face> = Vec::new();
        for face_id in shell_faces {
            if let Some((face, surface_ref)) = build_face(
                parsed,
                face_id,
                &vertex_index_by_id,
                &mut edges,
                &mut edge_index_by_curve,
                &mut geom,
                &mut curve_store_index_by_step,
            ) {
                let surface_binding = surface_ref.and_then(|step_surface| {
                    if let Some(idx) = surface_store_index_by_step.get(&step_surface) {
                        return Some(*idx);
                    }
                    let surface = resolve_surface(parsed, step_surface)?;
                    let idx = geom.surfaces.len();
                    geom.surfaces.push(surface);
                    surface_store_index_by_step.insert(step_surface, idx);
                    Some(idx)
                });
                geom.face_surface.push(surface_binding);
                faces.push(face);
            }
        }

        if !faces.is_empty() {
            solids.push(Solid {
                shells: vec![Shell { faces }],
            });
        }
    }

    if solids.is_empty() {
        if let Some(brep) = brep_from_points_bbox(parsed) {
            return Ok(brep);
        }
        return Err("STEP parse produced no triangulated faces".to_string());
    }

    Ok(BRep {
        vertices,
        edges,
        solids,
        geom,
    })
}

fn brep_from_points_bbox(parsed: &ParsedStep) -> Option<BRep> {
    if parsed.cartesian_points.is_empty() {
        return None;
    }

    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for point in parsed.cartesian_points.values() {
        for i in 0..3 {
            if point[i] < min[i] {
                min[i] = point[i];
            }
            if point[i] > max[i] {
                max[i] = point[i];
            }
        }
    }

    if min.iter().any(|v| !v.is_finite()) || max.iter().any(|v| !v.is_finite()) {
        return None;
    }

    let mut width = max[0] - min[0];
    let mut height = max[1] - min[1];
    let mut depth = max[2] - min[2];

    if width.abs() < 1e-9 {
        width = 1.0;
    }
    if height.abs() < 1e-9 {
        height = 1.0;
    }
    if depth.abs() < 1e-9 {
        depth = 1.0;
    }

    let mut brep = BRep::create_box(width, height, depth);
    for vertex in &mut brep.vertices {
        vertex.point.x += min[0];
        vertex.point.y += min[1];
        vertex.point.z += min[2];
    }
    Some(brep)
}

fn collect_shell_faces(parsed: &ParsedStep) -> Vec<Vec<u64>> {
    let mut shells = Vec::new();

    if !parsed.manifold_solids.is_empty() {
        for shell_id in &parsed.manifold_solids {
            if let Some(face_ids) = parsed.closed_shells.get(shell_id) {
                shells.push(face_ids.clone());
            }
        }
    }

    if shells.is_empty() && !parsed.closed_shells.is_empty() {
        for face_ids in parsed.closed_shells.values() {
            shells.push(face_ids.clone());
        }
    }

    if shells.is_empty() && !parsed.advanced_faces.is_empty() {
        let mut all_faces: Vec<u64> = parsed.advanced_faces.keys().copied().collect();
        all_faces.sort_unstable();
        shells.push(all_faces);
    }

    shells
}

fn collect_used_vertices(parsed: &ParsedStep, shell_face_sets: &[Vec<u64>]) -> Result<BTreeSet<u64>, String> {
    let mut used = BTreeSet::new();

    for shell in shell_face_sets {
        for face_id in shell {
            let bound_ids = parsed
                .advanced_faces
                .get(face_id)
                .ok_or_else(|| format!("Missing ADVANCED_FACE #{}", face_id))?;
            for bound_id in &bound_ids.bounds {
                let loop_id = parsed
                    .face_bounds
                    .get(bound_id)
                    .ok_or_else(|| format!("Missing FACE_BOUND #{}", bound_id))?;
                let oriented_ids = parsed
                    .edge_loops
                    .get(loop_id)
                    .ok_or_else(|| format!("Missing EDGE_LOOP #{}", loop_id))?;
                for oriented_id in oriented_ids {
                    let (edge_curve_id, _) = parsed
                        .oriented_edges
                        .get(oriented_id)
                        .ok_or_else(|| format!("Missing ORIENTED_EDGE #{}", oriented_id))?;
                    let (start, end, _) = parsed
                        .edge_curves
                        .get(edge_curve_id)
                        .ok_or_else(|| format!("Missing EDGE_CURVE #{}", edge_curve_id))?;
                    used.insert(*start);
                    used.insert(*end);
                }
            }
        }
    }

    Ok(used)
}

fn build_face(
    parsed: &ParsedStep,
    face_id: u64,
    vertex_index_by_id: &HashMap<u64, usize>,
    edges: &mut Vec<Edge>,
    edge_index_by_curve: &mut HashMap<u64, usize>,
    geom: &mut GeomStore,
    curve_store_index_by_step: &mut HashMap<u64, usize>,
) -> Option<(Face, Option<u64>)> {
    let bound_ids = parsed.advanced_faces.get(&face_id)?;
    let outer_bound = *bound_ids.bounds.first()?;
    let loop_id = *parsed.face_bounds.get(&outer_bound)?;
    let oriented_ids = parsed.edge_loops.get(&loop_id)?;

    let mut polygon: Vec<usize> = Vec::new();
    let mut wire_edge_indices: Vec<usize> = Vec::new();

    for oriented_id in oriented_ids {
        let (edge_curve_id, orientation) = *parsed.oriented_edges.get(oriented_id)?;
        let (start_id, end_id, curve_ref) = *parsed.edge_curves.get(&edge_curve_id)?;

        let (from_vertex_id, to_vertex_id) = if orientation {
            (start_id, end_id)
        } else {
            (end_id, start_id)
        };

        let from = *vertex_index_by_id.get(&from_vertex_id)?;
        let to = *vertex_index_by_id.get(&to_vertex_id)?;

        if polygon.is_empty() {
            polygon.push(from);
            polygon.push(to);
        } else {
            let last = *polygon.last()?;
            if last == from {
                polygon.push(to);
            } else if last == to {
                polygon.push(from);
            } else {
                polygon.push(from);
                polygon.push(to);
            }
        }

        let edge_index = if let Some(idx) = edge_index_by_curve.get(&edge_curve_id) {
            *idx
        } else {
            let idx = edges.len();
            edges.push(Edge {
                start: *vertex_index_by_id.get(&start_id)?,
                end: *vertex_index_by_id.get(&end_id)?,
            });
            edge_index_by_curve.insert(edge_curve_id, idx);

            if geom.edge_curve.len() <= idx {
                geom.edge_curve.resize(idx + 1, None);
            }
            geom.edge_curve[idx] = curve_ref.and_then(|step_curve| {
                if let Some(existing) = curve_store_index_by_step.get(&step_curve) {
                    return Some(*existing);
                }
                let curve = resolve_curve(parsed, step_curve)?;
                let cidx = geom.curves.len();
                geom.curves.push(curve);
                curve_store_index_by_step.insert(step_curve, cidx);
                Some(cidx)
            });
            idx
        };
        wire_edge_indices.push(edge_index);
    }

    while polygon.len() > 1 && polygon.first() == polygon.last() {
        polygon.pop();
    }
    dedup_consecutive(&mut polygon);

    if polygon.len() < 3 {
        return None;
    }

    let triangles = triangulate_fan(&polygon);
    if triangles.is_empty() {
        return None;
    }

    let normal = compute_normal(&triangles[0], vertex_index_by_id, parsed);

    Some((Face {
        outer_wire: Wire {
            edges: wire_edge_indices,
        },
        inner_wires: Vec::new(),
        normal,
        triangles,
    }, bound_ids.surface))
}

fn compute_normal(
    _triangle: &[usize; 3],
    _vertex_index_by_id: &HashMap<u64, usize>,
    _parsed: &ParsedStep,
) -> glam::DVec3 {
    // Normal is optional for current renderer, keep a stable default if degenerate.
    glam::DVec3::new(0.0, 0.0, 1.0)
}

fn triangulate_fan(polygon: &[usize]) -> Vec<[usize; 3]> {
    if polygon.len() < 3 {
        return Vec::new();
    }

    let mut triangles = Vec::with_capacity(polygon.len().saturating_sub(2));
    for i in 1..(polygon.len() - 1) {
        let tri = [polygon[0], polygon[i], polygon[i + 1]];
        if tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2] {
            triangles.push(tri);
        }
    }
    triangles
}

fn dedup_consecutive(polygon: &mut Vec<usize>) {
    if polygon.len() < 2 {
        return;
    }
    let mut deduped = Vec::with_capacity(polygon.len());
    let mut last = None;
    for &idx in polygon.iter() {
        if Some(idx) != last {
            deduped.push(idx);
            last = Some(idx);
        }
    }
    *polygon = deduped;
}

fn extract_data_section(content: &str) -> Result<&str, String> {
    let start = content
        .find("DATA;")
        .ok_or_else(|| "STEP file missing DATA section".to_string())?;
    let after_start = &content[start + "DATA;".len()..];
    let end = after_start
        .find("ENDSEC;")
        .ok_or_else(|| "STEP file missing ENDSEC after DATA".to_string())?;
    Ok(&after_start[..end])
}

fn split_records(data: &str) -> Vec<String> {
    let mut records = Vec::new();
    let mut current = String::new();
    let mut in_comment = false;
    let mut chars = data.chars().peekable();

    while let Some(ch) = chars.next() {
        if in_comment {
            if ch == '*' && chars.peek() == Some(&'/') {
                let _ = chars.next();
                in_comment = false;
            }
            continue;
        }

        if ch == '/' && chars.peek() == Some(&'*') {
            let _ = chars.next();
            in_comment = true;
            continue;
        }

        if ch == ';' {
            let trimmed = current.trim();
            if !trimmed.is_empty() {
                records.push(trimmed.to_string());
            }
            current.clear();
            continue;
        }

        current.push(ch);
    }

    let trailing = current.trim();
    if !trailing.is_empty() {
        records.push(trailing.to_string());
    }

    records
}

fn parse_entity_record(record: &str) -> Result<Option<(u64, &str)>, String> {
    let line = record.trim();
    if !line.starts_with('#') {
        return Ok(None);
    }

    let eq = line
        .find('=')
        .ok_or_else(|| format!("Invalid STEP entity record: {}", line))?;
    let id_str = line[1..eq].trim();
    let id = id_str
        .parse::<u64>()
        .map_err(|e| format!("Invalid STEP entity id {}: {}", id_str, e))?;
    Ok(Some((id, line[eq + 1..].trim())))
}

fn parse_entity_body(body: &str) -> Option<(&str, &str)> {
    let mut payload = body.trim();
    if payload.starts_with('(') {
        payload = payload.strip_prefix('(')?.strip_suffix(')')?.trim();
    }

    let open = payload.find('(')?;
    let close = payload.rfind(')')?;
    let entity = payload[..open].trim();
    let args = &payload[open + 1..close];
    Some((entity, args))
}

fn parse_cartesian_point(args: &str) -> Option<[f64; 3]> {
    let list = parse_coord_list(args)?;
    if list.len() < 3 {
        return None;
    }
    Some([list[0], list[1], list[2]])
}

fn parse_coord_list(args: &str) -> Option<Vec<f64>> {
    let open = args.rfind('(')?;
    let close = args.rfind(')')?;
    if close <= open {
        return None;
    }
    let raw = &args[open + 1..close];
    let mut coords = Vec::new();
    for item in raw.split(',') {
        let v = item.trim().parse::<f64>().ok()?;
        coords.push(v);
    }
    Some(coords)
}

fn parse_single_ref_after_name(args: &str) -> Option<u64> {
    let parts = split_top_level(args, ',');
    for part in parts.into_iter().skip(1) {
        if let Some(reference) = parse_ref(part) {
            return Some(reference);
        }
    }
    None
}

fn parse_edge_curve_vertices(args: &str) -> Option<(u64, u64, Option<u64>)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    let start = parse_ref(parts[1])?;
    let end = parse_ref(parts[2])?;
    let curve_ref = parse_ref(parts[3]);
    Some((start, end, curve_ref))
}

fn parse_axis2_placement(args: &str) -> Option<(u64, u64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    Some((parse_ref(parts[1])?, parse_ref(parts[2])?))
}

fn parse_curve_basis(args: &str) -> Option<(u64, u64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    Some((parse_ref(parts[1])?, parse_ref(parts[2])?))
}

fn parse_placement_radius(args: &str) -> Option<(u64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    Some((parse_ref(parts[1])?, parts[2].trim().parse::<f64>().ok()?))
}

fn parse_conical_surface(args: &str) -> Option<(u64, f64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    Some((
        parse_ref(parts[1])?,
        parts[2].trim().parse::<f64>().ok()?,
        parts[3].trim().parse::<f64>().ok()?.to_radians(),
    ))
}

fn parse_toroidal_surface(args: &str) -> Option<(u64, f64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    Some((
        parse_ref(parts[1])?,
        parts[2].trim().parse::<f64>().ok()?,
        parts[3].trim().parse::<f64>().ok()?,
    ))
}

fn parse_advanced_face(args: &str) -> Option<(Vec<u64>, Option<u64>)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    let bounds = parse_ref_list(parts[1]);
    let surface = parse_ref(parts[2]);
    if bounds.is_empty() {
        return None;
    }
    Some((bounds, surface))
}

fn resolve_curve(parsed: &ParsedStep, curve_ref: u64) -> Option<Curve3> {
    if let Some((origin_point, vector_ref)) = parsed.lines.get(&curve_ref) {
        let origin = point_from_ref(parsed, *origin_point)?;
        let direction_ref = *parsed.vectors.get(vector_ref)?;
        let direction = direction_from_ref(parsed, direction_ref)?;
        return Some(Curve3::Line(rcad_kernel::geom::Line3 { origin, direction }));
    }

    if let Some((placement_ref, radius)) = parsed.circles.get(&curve_ref) {
        let (center, normal) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Curve3::Circle(rcad_kernel::geom::Circle3 {
            center,
            normal,
            radius: *radius,
        }));
    }

    None
}

fn resolve_surface(parsed: &ParsedStep, surface_ref: u64) -> Option<Surface3> {
    if let Some(placement_ref) = parsed.planes.get(&surface_ref) {
        let (origin, normal) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Plane(rcad_kernel::geom::Plane { origin, normal }));
    }

    if let Some((placement_ref, radius)) = parsed.cylindrical_surfaces.get(&surface_ref) {
        let (origin, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin,
            axis,
            radius: *radius,
        }));
    }

    if let Some((placement_ref, radius)) = parsed.spherical_surfaces.get(&surface_ref) {
        let (center, _) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center,
            radius: *radius,
        }));
    }

    if let Some((placement_ref, half_angle_rad)) = parsed.conical_surfaces.get(&surface_ref) {
        let (apex, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Cone(rcad_kernel::geom::ConicalSurface {
            apex,
            axis,
            half_angle_rad: *half_angle_rad,
        }));
    }

    if let Some((placement_ref, major_radius, minor_radius)) = parsed.toroidal_surfaces.get(&surface_ref) {
        let (center, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center,
            axis,
            major_radius: *major_radius,
            minor_radius: *minor_radius,
        }));
    }

    None
}

fn point_from_ref(parsed: &ParsedStep, point_ref: u64) -> Option<glam::DVec3> {
    let p = parsed.cartesian_points.get(&point_ref)?;
    Some(glam::DVec3::new(p[0], p[1], p[2]))
}

fn direction_from_ref(parsed: &ParsedStep, direction_ref: u64) -> Option<glam::DVec3> {
    let d = parsed.directions.get(&direction_ref)?;
    Some(glam::DVec3::new(d[0], d[1], d[2]).normalize_or_zero())
}

fn placement_from_ref(parsed: &ParsedStep, placement_ref: u64) -> Option<(glam::DVec3, glam::DVec3)> {
    let (origin_ref, axis_ref) = *parsed.axis2_placements.get(&placement_ref)?;
    Some((point_from_ref(parsed, origin_ref)?, direction_from_ref(parsed, axis_ref)?))
}

fn parse_oriented_edge(args: &str) -> Option<(u64, bool)> {
    let parts = split_top_level(args, ',');
    let mut edge_ref = None;
    let mut orientation = None;

    for part in &parts {
        if edge_ref.is_none() {
            edge_ref = parse_ref(part);
            if edge_ref.is_some() {
                continue;
            }
        }

        if orientation.is_none() {
            let v = part.trim();
            if v == ".T." {
                orientation = Some(true);
            } else if v == ".F." {
                orientation = Some(false);
            }
        }
    }

    Some((edge_ref?, orientation?))
}

fn parse_ref_list_after_name(args: &str) -> Option<Vec<u64>> {
    let open = args.find('(')?;
    let close = args.rfind(')')?;
    if close <= open {
        return None;
    }
    let inside = args[open + 1..close].trim();
    if !inside.starts_with('#') && !inside.contains('#') {
        return None;
    }
    Some(parse_ref_list(inside))
}

fn parse_ref_list(input: &str) -> Vec<u64> {
    let mut refs = Vec::new();
    let mut i = 0usize;
    let bytes = input.as_bytes();

    while i < bytes.len() {
        if bytes[i] == b'#' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            if start < i && let Ok(v) = input[start..i].parse::<u64>() {
                refs.push(v);
            }
        } else {
            i += 1;
        }
    }

    refs
}

fn parse_ref(input: &str) -> Option<u64> {
    let trimmed = input.trim();
    let hash = trimmed.find('#')?;
    let digits = &trimmed[hash + 1..];
    let mut end = 0usize;
    for ch in digits.chars() {
        if ch.is_ascii_digit() {
            end += ch.len_utf8();
        } else {
            break;
        }
    }
    if end == 0 {
        return None;
    }
    digits[..end].parse::<u64>().ok()
}

fn split_top_level(input: &str, delimiter: char) -> Vec<&str> {
    let mut result = Vec::new();
    let mut depth = 0i32;
    let mut in_string = false;
    let mut start = 0usize;

    for (idx, ch) in input.char_indices() {
        match ch {
            '\'' => in_string = !in_string,
            '(' if !in_string => depth += 1,
            ')' if !in_string => depth -= 1,
            _ => {}
        }

        if ch == delimiter && depth == 0 && !in_string {
            result.push(input[start..idx].trim());
            start = idx + ch.len_utf8();
        }
    }

    if start <= input.len() {
        result.push(input[start..].trim());
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    const HFSS_STEP: &str = include_str!("../../../assets/hfss.step");
    const BOX_STEP: &str = include_str!("../../../assets/box.step");

    #[test]
    fn parses_hfss_into_non_trivial_brep() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        assert!(brep.vertices.len() > 8, "hfss should have more than box vertices");
        assert!(!brep.edges.is_empty(), "hfss should produce edges");

        let triangle_count: usize = brep
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .flat_map(|sh| sh.faces.iter())
            .map(|f| f.triangles.len())
            .sum();

        assert!(triangle_count > 0, "hfss should produce triangulated faces");
    }

    #[test]
    fn parses_box_example() {
        let brep = StepReader::parse_string(BOX_STEP).expect("box.step should parse");
        assert!(!brep.vertices.is_empty());
    }
}
