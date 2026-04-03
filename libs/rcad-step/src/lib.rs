use rcad_kernel::{BRep, Curve3, GeomStore, Surface3};
use rcad_kernel::{Edge, Face, Shell, Solid, Vertex, Wire};
use rcad_modeling::make_box_brep;
use std::collections::{BTreeSet, HashMap};
use std::path::Path;

pub mod writer;

#[derive(Debug, Clone)]
struct AdvancedFaceRecord {
    bounds: Vec<u64>,
    surface: Option<u64>,
}

#[derive(Debug, Clone)]
struct ParsedStep {
    cartesian_points: HashMap<u64, [f64; 3]>,
    directions: HashMap<u64, [f64; 3]>,
    vectors: HashMap<u64, (u64, f64)>,
    axis2_placements: HashMap<u64, (u64, u64, Option<u64>)>,
    lines: HashMap<u64, (u64, u64)>,
    circles: HashMap<u64, (u64, f64)>,
    ellipses: HashMap<u64, (u64, f64, f64)>,
    b_spline_curves: HashMap<u64, Vec<u64>>,
    planes: HashMap<u64, u64>,
    cylindrical_surfaces: HashMap<u64, (u64, f64)>,
    spherical_surfaces: HashMap<u64, (u64, f64)>,
    conical_surfaces: HashMap<u64, (u64, f64, f64)>,
    toroidal_surfaces: HashMap<u64, (u64, f64, f64)>,
    vertex_points: HashMap<u64, u64>,
    edge_curves: HashMap<u64, (u64, u64, Option<u64>)>,
    oriented_edges: HashMap<u64, (u64, bool)>,
    edge_loops: HashMap<u64, Vec<u64>>,
    face_bounds: HashMap<u64, u64>,
    advanced_faces: HashMap<u64, AdvancedFaceRecord>,
    closed_shells: HashMap<u64, Vec<u64>>,
    open_shells: HashMap<u64, Vec<u64>>,
    manifold_solids: Vec<u64>,
    shell_based_surface_models: Vec<Vec<u64>>,
    trimmed_curves: HashMap<u64, (u64, f64, f64)>,
    geometric_curve_sets: Vec<Vec<u64>>,
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
            ellipses: HashMap::new(),
            b_spline_curves: HashMap::new(),
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
            open_shells: HashMap::new(),
            manifold_solids: Vec::new(),
            shell_based_surface_models: Vec::new(),
            trimmed_curves: HashMap::new(),
            geometric_curve_sets: Vec::new(),
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
                    if let Some((dir_ref, magnitude)) = parse_vector(args) {
                        parsed.vectors.insert(id, (dir_ref, magnitude));
                    }
                }
                "AXIS2_PLACEMENT_3D" => {
                    if let Some((origin, axis, ref_dir)) = parse_axis2_placement(args) {
                        parsed.axis2_placements.insert(id, (origin, axis, ref_dir));
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
                "ELLIPSE" => {
                    if let Some((placement, major, minor)) = parse_placement_two_radii(args) {
                        parsed.ellipses.insert(id, (placement, major, minor));
                    }
                }
                "B_SPLINE_CURVE_WITH_KNOTS" => {
                    if let Some(points) = parse_bspline_control_points(args)
                        && !points.is_empty()
                    {
                        parsed.b_spline_curves.insert(id, points);
                    }
                }
                "TRIMMED_CURVE" => {
                    if let Some((curve_ref, t0, t1)) = parse_trimmed_curve(args) {
                        parsed.trimmed_curves.insert(id, (curve_ref, t0, t1));
                    }
                }
                "GEOMETRIC_CURVE_SET" => {
                    if let Some(curve_refs) = parse_ref_list_after_name(args)
                        && !curve_refs.is_empty()
                    {
                        parsed.geometric_curve_sets.push(curve_refs);
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
                        parsed
                            .conical_surfaces
                            .insert(id, (placement, radius, half_angle_rad));
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
                "OPEN_SHELL" => {
                    if let Some(face_refs) = parse_ref_list_after_name(args) {
                        parsed.open_shells.insert(id, face_refs);
                    }
                }
                "MANIFOLD_SOLID_BREP" => {
                    if let Some(shell_ref) = parse_single_ref_after_name(args) {
                        parsed.manifold_solids.push(shell_ref);
                    }
                }
                "SHELL_BASED_SURFACE_MODEL" => {
                    if let Some(shell_refs) = parse_ref_list_after_name(args)
                        && !shell_refs.is_empty()
                    {
                        parsed.shell_based_surface_models.push(shell_refs);
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
    let used_vertex_ids = if shell_face_sets.is_empty() {
        collect_edge_vertices(parsed)
    } else {
        let mut used = collect_used_vertices(parsed, &shell_face_sets)?;
        used.extend(collect_edge_vertices(parsed));
        used
    };
    if used_vertex_ids.is_empty() && parsed.geometric_curve_sets.is_empty() {
        if let Some(brep) = brep_from_points_bbox(parsed) {
            return Ok(brep);
        }
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
                &mut vertices,
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

    // Preserve standalone edges that are not part of any face loop.
    for (edge_curve_id, (start_id, end_id, curve_ref)) in &parsed.edge_curves {
        if edge_index_by_curve.contains_key(edge_curve_id) {
            continue;
        }
        let Some(&start) = vertex_index_by_id.get(start_id) else {
            continue;
        };
        let Some(&end) = vertex_index_by_id.get(end_id) else {
            continue;
        };

        let idx = edges.len();
        edges.push(Edge { start, end });
        edge_index_by_curve.insert(*edge_curve_id, idx);

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
    }

    // Sample standalone 1D curves from GEOMETRIC_CURVE_SET (Polyline1/2/3 etc.)
    for curve_set in &parsed.geometric_curve_sets {
        for &curve_ref in curve_set {
            let points = sample_standalone_curve(parsed, curve_ref);
            if points.len() < 2 {
                continue;
            }
            let base = vertices.len();
            for p in &points {
                vertices.push(Vertex { point: *p });
            }
            for i in 0..points.len() - 1 {
                let idx = edges.len();
                edges.push(Edge {
                    start: base + i,
                    end: base + i + 1,
                });
                if geom.edge_curve.len() <= idx {
                    geom.edge_curve.resize(idx + 1, None);
                }
                // No curve binding needed — already tessellated into polyline
            }
        }
    }

    if solids.is_empty() && edges.is_empty() {
        if let Some(brep) = brep_from_points_bbox(parsed) {
            return Ok(brep);
        }
        return Err("STEP parse produced no triangulated faces or edges".to_string());
    }

    Ok(BRep {
        vertices,
        edges,
        solids,
        geom,
    })
}

fn collect_edge_vertices(parsed: &ParsedStep) -> BTreeSet<u64> {
    let mut used = BTreeSet::new();
    for (start, end, _) in parsed.edge_curves.values() {
        used.insert(*start);
        used.insert(*end);
    }
    used
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

    make_box_brep(
        glam::DVec3::new(min[0], min[1], min[2]),
        glam::DVec3::X,
        glam::DVec3::Y,
        width,
        height,
        depth,
    )
    .ok()
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

    if !parsed.shell_based_surface_models.is_empty() {
        for shell_refs in &parsed.shell_based_surface_models {
            for shell_id in shell_refs {
                if let Some(face_ids) = parsed.open_shells.get(shell_id) {
                    shells.push(face_ids.clone());
                } else if let Some(face_ids) = parsed.closed_shells.get(shell_id) {
                    shells.push(face_ids.clone());
                }
            }
        }
    }

    if shells.is_empty() && !parsed.closed_shells.is_empty() {
        for face_ids in parsed.closed_shells.values() {
            shells.push(face_ids.clone());
        }
    }

    if shells.is_empty() && !parsed.open_shells.is_empty() {
        for face_ids in parsed.open_shells.values() {
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
    vertices: &mut Vec<Vertex>,
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
    let mut face_vertex_indices: Vec<usize> = Vec::new();
    let mut sampled_loop_points: Vec<glam::DVec3> = Vec::new();

    // Detect seam edges: an edge_curve that appears twice in the same face boundary
    let mut edge_curve_count: HashMap<u64, usize> = HashMap::new();
    for oriented_id in oriented_ids {
        let (edge_curve_id, _) = *parsed.oriented_edges.get(oriented_id)?;
        *edge_curve_count.entry(edge_curve_id).or_insert(0) += 1;
    }
    let has_seam = edge_curve_count.values().any(|&c| c >= 2);

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
        face_vertex_indices.push(from);
        face_vertex_indices.push(to);

        if let Some(points) = sample_oriented_edge_points(parsed, edge_curve_id, orientation) {
            append_edge_points(&mut sampled_loop_points, &points);
        }

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

    let triangles = if polygon.len() >= 3 {
        triangulate_fan(&polygon)
    } else if sampled_loop_points.len() >= 3
        && bound_ids
            .surface
            .map(|sid| parsed.planes.contains_key(&sid))
            .unwrap_or(true)
    {
        triangulate_point_loop(vertices, &sampled_loop_points)
    } else if let Some(surface_ref) = bound_ids.surface {
        triangulate_surface_fallback(parsed, surface_ref, vertices, &face_vertex_indices, has_seam)
    } else {
        Vec::new()
    };
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

fn triangulate_point_loop(vertices: &mut Vec<Vertex>, points: &[glam::DVec3]) -> Vec<[usize; 3]> {
    let mut loop_points: Vec<glam::DVec3> = points.to_vec();
    dedup_consecutive_points(&mut loop_points);
    if loop_points.len() > 1 {
        let first = loop_points[0];
        let last = *loop_points.last().unwrap_or(&first);
        if first.distance(last) < 1e-6 {
            loop_points.pop();
        }
    }

    if loop_points.len() < 3 {
        return Vec::new();
    }

    let base = vertices.len();
    for p in &loop_points {
        vertices.push(Vertex { point: *p });
    }

    let mut triangles = Vec::with_capacity(loop_points.len().saturating_sub(2));
    for i in 1..(loop_points.len() - 1) {
        triangles.push([base, base + i, base + i + 1]);
    }
    triangles
}

fn append_edge_points(dst: &mut Vec<glam::DVec3>, src: &[glam::DVec3]) {
    for p in src {
        if let Some(last) = dst.last()
            && last.distance(*p) < 1e-7
        {
            continue;
        }
        dst.push(*p);
    }
}

fn dedup_consecutive_points(points: &mut Vec<glam::DVec3>) {
    if points.len() < 2 {
        return;
    }
    let mut deduped = Vec::with_capacity(points.len());
    for p in points.iter().copied() {
        if deduped
            .last()
            .map(|last: &glam::DVec3| last.distance(p) < 1e-7)
            .unwrap_or(false)
        {
            continue;
        }
        deduped.push(p);
    }
    *points = deduped;
}

fn sample_oriented_edge_points(
    parsed: &ParsedStep,
    edge_curve_id: u64,
    orientation: bool,
) -> Option<Vec<glam::DVec3>> {
    let (start_id, end_id, curve_ref) = *parsed.edge_curves.get(&edge_curve_id)?;
    let start = vertex_point_from_ref(parsed, start_id)?;
    let end = vertex_point_from_ref(parsed, end_id)?;

    let mut points = if let Some(curve_id) = curve_ref {
        if let Some((placement_ref, radius)) = parsed.circles.get(&curve_id) {
            sample_circle_edge(parsed, *placement_ref, *radius, start, end)?
        } else if let Some((placement_ref, major, minor)) = parsed.ellipses.get(&curve_id) {
            sample_ellipse_edge(parsed, *placement_ref, *major, *minor, start, end)?
        } else if let Some(control_refs) = parsed.b_spline_curves.get(&curve_id) {
            sample_bspline_polyline(parsed, control_refs)
        } else {
            vec![start, end]
        }
    } else {
        vec![start, end]
    };

    if !orientation {
        points.reverse();
    }
    Some(points)
}

fn sample_circle_edge(
    parsed: &ParsedStep,
    placement_ref: u64,
    radius: f64,
    start: glam::DVec3,
    end: glam::DVec3,
) -> Option<Vec<glam::DVec3>> {
    if !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    let (center, axis, u) = placement_frame_from_ref(parsed, placement_ref)?;
    let v = axis.cross(u).normalize_or_zero();
    let a0 = angle_on_basis(start - center, u, v);
    let a1 = angle_on_basis(end - center, u, v);
    let mut sweep = a1 - a0;
    if start.distance(end) < 1e-6 {
        sweep = std::f64::consts::TAU;
    } else if sweep <= 0.0 {
        sweep += std::f64::consts::TAU;
    }
    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0).ceil().max(8.0) as usize;
    let mut points = Vec::with_capacity(seg + 1);
    for i in 0..=seg {
        let t = a0 + sweep * (i as f64 / seg as f64);
        points.push(center + u * (radius * t.cos()) + v * (radius * t.sin()));
    }
    Some(points)
}

fn sample_ellipse_edge(
    parsed: &ParsedStep,
    placement_ref: u64,
    major: f64,
    minor: f64,
    start: glam::DVec3,
    end: glam::DVec3,
) -> Option<Vec<glam::DVec3>> {
    if !major.is_finite() || !minor.is_finite() || major <= 0.0 || minor <= 0.0 {
        return None;
    }
    let (center, axis, u) = placement_frame_from_ref(parsed, placement_ref)?;
    let v = axis.cross(u).normalize_or_zero();

    let param = |p: glam::DVec3| {
        let d = p - center;
        let x = d.dot(u) / major;
        let y = d.dot(v) / minor;
        y.atan2(x)
    };

    let t0 = param(start);
    let t1 = param(end);
    let mut sweep = t1 - t0;
    if start.distance(end) < 1e-6 {
        sweep = std::f64::consts::TAU;
    } else if sweep <= 0.0 {
        sweep += std::f64::consts::TAU;
    }
    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0).ceil().max(8.0) as usize;
    let mut points = Vec::with_capacity(seg + 1);
    for i in 0..=seg {
        let t = t0 + sweep * (i as f64 / seg as f64);
        points.push(center + u * (major * t.cos()) + v * (minor * t.sin()));
    }
    Some(points)
}

fn sample_bspline_polyline(parsed: &ParsedStep, control_refs: &[u64]) -> Vec<glam::DVec3> {
    let mut points = Vec::new();
    for cref in control_refs {
        if let Some(p) = point_from_ref(parsed, *cref) {
            points.push(p);
        }
    }
    points
}

fn angle_on_basis(v: glam::DVec3, u: glam::DVec3, w: glam::DVec3) -> f64 {
    v.dot(w).atan2(v.dot(u))
}

fn triangulate_spherical_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, radius)) = parsed.spherical_surfaces.get(&surface_ref) else {
        return Vec::new();
    };
    let Some((center, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !radius.is_finite() || *radius <= 0.0 {
        return Vec::new();
    }

    let u = ref_dir;
    let v = axis.cross(u).normalize_or_zero();
    let w = axis.normalize_or_zero();

    let (theta_min, theta_max) = infer_angular_range(vertices, face_vertex_indices, center, w, u, v, has_seam);
    let (phi_min, phi_max) = infer_polar_range(vertices, face_vertex_indices, center, w, *radius, has_seam);

    let u_segments = 32usize;
    let v_segments = 16usize;
    let stride = u_segments + 1;
    let base = vertices.len();

    let at_north_pole = phi_min.abs() < 1e-6;
    let at_south_pole = (phi_max - std::f64::consts::PI).abs() < 1e-6;

    for vi in 0..=v_segments {
        let phi = phi_min + (phi_max - phi_min) * vi as f64 / v_segments as f64;
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        for ui in 0..=u_segments {
            let theta = theta_min + (theta_max - theta_min) * ui as f64 / u_segments as f64;
            let point = center
                + w * (radius * cos_phi)
                + u * (radius * theta.cos() * sin_phi)
                + v * (radius * theta.sin() * sin_phi);
            vertices.push(Vertex { point });
        }
    }

    let mut triangles = Vec::with_capacity(u_segments * v_segments * 2);
    for vi in 0..v_segments {
        for ui in 0..u_segments {
            let i0 = base + vi * stride + ui;
            let i1 = i0 + 1;
            let i2 = i0 + stride;
            let i3 = i2 + 1;

            // Skip degenerate triangles at poles
            if !(at_north_pole && vi == 0) {
                triangles.push([i0, i2, i1]);
            }
            if !(at_south_pole && vi == v_segments - 1) {
                triangles.push([i1, i2, i3]);
            }
        }
    }

    triangles
}

fn triangulate_surface_fallback(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    if parsed.spherical_surfaces.contains_key(&surface_ref) {
        return triangulate_spherical_surface(parsed, surface_ref, vertices, face_vertex_indices, has_seam);
    }
    if parsed.cylindrical_surfaces.contains_key(&surface_ref) {
        return triangulate_cylindrical_surface(parsed, surface_ref, vertices, face_vertex_indices, has_seam);
    }
    if parsed.conical_surfaces.contains_key(&surface_ref) {
        return triangulate_conical_surface(parsed, surface_ref, vertices, face_vertex_indices, has_seam);
    }
    if parsed.toroidal_surfaces.contains_key(&surface_ref) {
        return triangulate_toroidal_surface(parsed, surface_ref, vertices, face_vertex_indices, has_seam);
    }
    Vec::new()
}

fn triangulate_cylindrical_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, radius)) = parsed.cylindrical_surfaces.get(&surface_ref) else {
        return Vec::new();
    };
    let Some((origin, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !radius.is_finite() || *radius <= 0.0 {
        return Vec::new();
    }

    let u = ref_dir;
    let v = axis.cross(u).normalize_or_zero();
    let w = axis.normalize_or_zero();
    let (t_min, t_max) = infer_axis_range(vertices, face_vertex_indices, origin, w, -*radius, *radius);
    let (theta_min, theta_max) = infer_angular_range(vertices, face_vertex_indices, origin, w, u, v, has_seam);

    let radial_segments = 40usize;
    let height_segments = 12usize;
    let stride = radial_segments + 1;
    let base = vertices.len();

    for j in 0..=height_segments {
        let tj = t_min + (t_max - t_min) * (j as f64 / height_segments as f64);
        let center = origin + w * tj;
        for i in 0..=radial_segments {
            let theta = theta_min + (theta_max - theta_min) * (i as f64 / radial_segments as f64);
            let ring_dir = u * theta.cos() + v * theta.sin();
            vertices.push(Vertex {
                point: center + ring_dir * *radius,
            });
        }
    }

    triangulate_grid(base, height_segments, radial_segments, stride)
}

fn triangulate_conical_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, ref_radius, half_angle_rad)) = parsed.conical_surfaces.get(&surface_ref)
    else {
        return Vec::new();
    };
    let Some((origin, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !ref_radius.is_finite() || !half_angle_rad.is_finite() {
        return Vec::new();
    }

    let tan_a = half_angle_rad.tan();
    if !tan_a.is_finite() || tan_a.abs() < 1e-8 {
        return Vec::new();
    }

    let u = ref_dir;
    let v = axis.cross(u).normalize_or_zero();
    let w = axis.normalize_or_zero();
    let default_h = ref_radius.abs().max(0.1) / tan_a.abs();
    let (t_min, t_max) = infer_axis_range(vertices, face_vertex_indices, origin, w, 0.0, default_h);
    let (theta_min, theta_max) = infer_angular_range(vertices, face_vertex_indices, origin, w, u, v, has_seam);

    let radial_segments = 40usize;
    let height_segments = 12usize;
    let stride = radial_segments + 1;
    let base = vertices.len();

    for j in 0..=height_segments {
        let tj = t_min + (t_max - t_min) * (j as f64 / height_segments as f64);
        let center = origin + w * tj;
        let rj = (ref_radius + tj * tan_a).abs().max(1e-5);
        for i in 0..=radial_segments {
            let theta = theta_min + (theta_max - theta_min) * (i as f64 / radial_segments as f64);
            let ring_dir = u * theta.cos() + v * theta.sin();
            vertices.push(Vertex {
                point: center + ring_dir * rj,
            });
        }
    }

    triangulate_grid(base, height_segments, radial_segments, stride)
}

fn triangulate_toroidal_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, major_radius, minor_radius)) = parsed.toroidal_surfaces.get(&surface_ref)
    else {
        return Vec::new();
    };
    let Some((center, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !major_radius.is_finite()
        || !minor_radius.is_finite()
        || *major_radius <= 0.0
        || *minor_radius <= 0.0
    {
        return Vec::new();
    }

    let u_dir = ref_dir;
    let v_dir = axis.cross(u_dir).normalize_or_zero();
    let w = axis.normalize_or_zero();

    // Infer major angle (theta) range from boundary vertices
    let (theta_min, theta_max) = infer_angular_range(vertices, face_vertex_indices, center, w, u_dir, v_dir, has_seam);

    // Infer minor angle (phi) range: project boundary vertices into the minor circle plane
    let (phi_min, phi_max) = infer_torus_minor_range(
        vertices, face_vertex_indices, center, w, u_dir, v_dir,
        *major_radius, has_seam,
    );

    let major_segments = 48usize;
    let minor_segments = 24usize;
    let stride = minor_segments + 1;
    let base = vertices.len();

    for i in 0..=major_segments {
        let theta = theta_min + (theta_max - theta_min) * i as f64 / major_segments as f64;
        let ring_dir = u_dir * theta.cos() + v_dir * theta.sin();
        let ring_center = center + ring_dir * *major_radius;
        for j in 0..=minor_segments {
            let phi = phi_min + (phi_max - phi_min) * j as f64 / minor_segments as f64;
            let minor_dir = ring_dir * phi.cos() + w * phi.sin();
            vertices.push(Vertex {
                point: ring_center + minor_dir * *minor_radius,
            });
        }
    }

    triangulate_grid(base, major_segments, minor_segments, stride)
}

/// Infer the angular (theta) range of boundary vertices projected onto a
/// surface local frame (origin, u_dir, v_dir, axis).
/// Returns (theta_min, theta_max) in radians.
/// If `has_seam` is true, the face wraps the full period → returns (0, 2π).
fn infer_angular_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    origin: glam::DVec3,
    axis: glam::DVec3,
    u_dir: glam::DVec3,
    v_dir: glam::DVec3,
    has_seam: bool,
) -> (f64, f64) {
    use std::f64::consts::TAU;
    if has_seam {
        return (0.0, TAU);
    }

    let mut angles: Vec<f64> = Vec::new();
    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let d = v.point - origin;
            let proj = d - axis * d.dot(axis);
            let a = proj.dot(v_dir).atan2(proj.dot(u_dir));
            angles.push(a);
        }
    }

    if angles.is_empty() {
        return (0.0, TAU);
    }

    // Sort angles and find the largest gap; the face covers the complement.
    angles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    angles.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

    if angles.len() < 2 {
        return (0.0, TAU);
    }

    let mut max_gap = 0.0f64;
    let mut gap_end = 0usize;
    for i in 1..angles.len() {
        let gap = angles[i] - angles[i - 1];
        if gap > max_gap {
            max_gap = gap;
            gap_end = i;
        }
    }
    // Wrap-around gap
    let wrap_gap = (angles[0] + TAU) - angles[angles.len() - 1];
    if wrap_gap > max_gap {
        // Face goes from angles[0] to angles[last]
        return (angles[0], angles[angles.len() - 1]);
    }

    // Face goes from angles[gap_end] around to angles[gap_end-1]
    let theta_min = angles[gap_end];
    let theta_max = angles[gap_end - 1] + TAU;
    (theta_min, theta_max)
}

/// Infer polar (phi) range for spherical surfaces from boundary vertices.
/// Returns (phi_min, phi_max) in [0, π] (north-pole to south-pole).
fn infer_polar_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    center: glam::DVec3,
    axis: glam::DVec3,
    radius: f64,
    has_seam: bool,
) -> (f64, f64) {
    let _ = has_seam;
    let mut phi_min = std::f64::consts::PI;
    let mut phi_max = 0.0f64;

    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let d = (v.point - center).normalize_or_zero();
            let cos_phi = d.dot(axis).clamp(-1.0, 1.0);
            let phi = cos_phi.acos();
            phi_min = phi_min.min(phi);
            phi_max = phi_max.max(phi);
        }
    }

    let _ = radius;
    if !phi_min.is_finite() || !phi_max.is_finite() || (phi_max - phi_min).abs() < 1e-8 {
        (0.0, std::f64::consts::PI)
    } else {
        (phi_min, phi_max)
    }
}

/// Infer the minor (phi) angular range for a torus face.
/// Projects each boundary vertex into the local minor circle plane at the
/// closest major angle, then computes the minor angle.
fn infer_torus_minor_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    center: glam::DVec3,
    axis: glam::DVec3,
    u_dir: glam::DVec3,
    v_dir: glam::DVec3,
    major_radius: f64,
    has_seam: bool,
) -> (f64, f64) {
    use std::f64::consts::TAU;
    if has_seam {
        return (0.0, TAU);
    }

    let mut angles: Vec<f64> = Vec::new();
    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let d = v.point - center;
            // Project onto the equatorial plane to find the major angle
            let proj_equat = d - axis * d.dot(axis);
            let theta = proj_equat.dot(v_dir).atan2(proj_equat.dot(u_dir));
            // Reconstruct ring_dir and ring_center for this major angle
            let ring_dir = u_dir * theta.cos() + v_dir * theta.sin();
            let ring_center = center + ring_dir * major_radius;
            // Vector from ring center to the vertex, projected onto the minor circle plane
            let to_v = v.point - ring_center;
            let phi = to_v.dot(axis).atan2(to_v.dot(ring_dir));
            angles.push(phi);
        }
    }

    if angles.len() < 2 {
        return (0.0, TAU);
    }

    angles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    angles.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

    if angles.len() < 2 {
        return (0.0, TAU);
    }

    let mut max_gap = 0.0f64;
    let mut gap_end = 0usize;
    for i in 1..angles.len() {
        let gap = angles[i] - angles[i - 1];
        if gap > max_gap {
            max_gap = gap;
            gap_end = i;
        }
    }
    let wrap_gap = (angles[0] + TAU) - angles[angles.len() - 1];
    if wrap_gap > max_gap {
        return (angles[0], angles[angles.len() - 1]);
    }
    let phi_min = angles[gap_end];
    let phi_max = angles[gap_end - 1] + TAU;
    (phi_min, phi_max)
}

fn triangulate_grid(
    base: usize,
    rows: usize,
    cols: usize,
    stride: usize,
) -> Vec<[usize; 3]> {
    let mut triangles = Vec::with_capacity(rows * cols * 2);
    for r in 0..rows {
        for c in 0..cols {
            let i0 = base + r * stride + c;
            let i1 = i0 + 1;
            let i2 = i0 + stride;
            let i3 = i2 + 1;
            triangles.push([i0, i2, i1]);
            triangles.push([i1, i2, i3]);
        }
    }
    triangles
}

fn infer_axis_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    origin: glam::DVec3,
    axis: glam::DVec3,
    default_min: f64,
    default_max: f64,
) -> (f64, f64) {
    let mut min_t = f64::INFINITY;
    let mut max_t = f64::NEG_INFINITY;

    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let t = (v.point - origin).dot(axis);
            min_t = min_t.min(t);
            max_t = max_t.max(t);
        }
    }

    if !min_t.is_finite() || !max_t.is_finite() || (max_t - min_t).abs() < 1e-8 {
        (default_min, default_max)
    } else {
        (min_t, max_t)
    }
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

fn parse_axis2_placement(args: &str) -> Option<(u64, u64, Option<u64>)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    let ref_dir = parts.get(3).and_then(|s| parse_ref(s));
    Some((parse_ref(parts[1])?, parse_ref(parts[2])?, ref_dir))
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

fn parse_placement_two_radii(args: &str) -> Option<(u64, f64, f64)> {
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

fn parse_bspline_control_points(args: &str) -> Option<Vec<u64>> {
    let parts = split_top_level(args, ',');
    let refs = parts.get(2).map(|s| parse_ref_list(s))?;
    if refs.is_empty() {
        return None;
    }
    Some(refs)
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
        let (direction_ref, _magnitude) = *parsed.vectors.get(vector_ref)?;
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

    if let Some((placement_ref, major_radius, minor_radius)) = parsed.ellipses.get(&curve_ref) {
        let (center, normal, major_dir) = placement_frame_from_ref(parsed, *placement_ref)?;
        return Some(Curve3::Ellipse(rcad_kernel::geom::Ellipse3 {
            center,
            normal,
            major_dir,
            major_radius: *major_radius,
            minor_radius: *minor_radius,
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
        let (center, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center,
            axis,
            radius: *radius,
        }));
    }

    if let Some((placement_ref, ref_radius, half_angle_rad)) = parsed.conical_surfaces.get(&surface_ref) {
        let (apex, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Cone(rcad_kernel::geom::ConicalSurface {
            apex,
            axis,
            radius: *ref_radius,
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

fn vertex_point_from_ref(parsed: &ParsedStep, vertex_ref: u64) -> Option<glam::DVec3> {
    let point_ref = *parsed.vertex_points.get(&vertex_ref)?;
    point_from_ref(parsed, point_ref)
}

fn direction_from_ref(parsed: &ParsedStep, direction_ref: u64) -> Option<glam::DVec3> {
    let d = parsed.directions.get(&direction_ref)?;
    Some(glam::DVec3::new(d[0], d[1], d[2]).normalize_or_zero())
}

fn placement_from_ref(parsed: &ParsedStep, placement_ref: u64) -> Option<(glam::DVec3, glam::DVec3)> {
    let (origin_ref, axis_ref, _) = *parsed.axis2_placements.get(&placement_ref)?;
    Some((point_from_ref(parsed, origin_ref)?, direction_from_ref(parsed, axis_ref)?))
}

fn placement_frame_from_ref(
    parsed: &ParsedStep,
    placement_ref: u64,
) -> Option<(glam::DVec3, glam::DVec3, glam::DVec3)> {
    let (origin_ref, axis_ref, ref_dir_ref) = *parsed.axis2_placements.get(&placement_ref)?;
    let origin = point_from_ref(parsed, origin_ref)?;
    let axis = direction_from_ref(parsed, axis_ref)?;

    let major_dir = if let Some(dref) = ref_dir_ref {
        let d = direction_from_ref(parsed, dref)?;
        let d_proj = (d - axis * d.dot(axis)).normalize_or_zero();
        if d_proj.length_squared() > 1e-8 {
            d_proj
        } else {
            any_perpendicular(axis)
        }
    } else {
        any_perpendicular(axis)
    };

    Some((origin, axis, major_dir))
}

fn any_perpendicular(axis: glam::DVec3) -> glam::DVec3 {
    let helper = if axis.dot(glam::DVec3::Y).abs() < 0.9 {
        glam::DVec3::Y
    } else {
        glam::DVec3::X
    };
    axis.cross(helper).normalize_or_zero()
}

fn parse_vector(args: &str) -> Option<(u64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    let dir_ref = parse_ref(parts[1])?;
    let magnitude = parts[2].trim().parse::<f64>().ok()?;
    Some((dir_ref, magnitude))
}

fn parse_trimmed_curve(args: &str) -> Option<(u64, f64, f64)> {
    // TRIMMED_CURVE('name', curve_ref, (PARAMETER_VALUE(t0)), (PARAMETER_VALUE(t1)), ...)
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    let curve_ref = parse_ref(parts[1])?;
    let t0 = parse_parameter_value(parts[2])?;
    let t1 = parse_parameter_value(parts[3])?;
    Some((curve_ref, t0, t1))
}

fn parse_parameter_value(s: &str) -> Option<f64> {
    // s looks like "(PARAMETER_VALUE(0.))" — find the float inside PARAMETER_VALUE(...)
    let cursor = s.to_uppercase();
    let pv_pos = cursor.find("PARAMETER_VALUE(")?;
    let after = &s[pv_pos + "PARAMETER_VALUE(".len()..];
    let end = after.find(')')?;
    after[..end].trim().parse::<f64>().ok()
}

/// Evaluate a STEP LINE curve at parameter t: p(t) = origin + dir * magnitude * t
fn eval_line_at(parsed: &ParsedStep, line_ref: u64, t: f64) -> Option<glam::DVec3> {
    let &(origin_ref, vec_ref) = parsed.lines.get(&line_ref)?;
    let origin = point_from_ref(parsed, origin_ref)?;
    let &(dir_ref, magnitude) = parsed.vectors.get(&vec_ref)?;
    let dir = direction_from_ref(parsed, dir_ref)?;
    Some(origin + dir * (magnitude * t))
}

/// Sample a standalone curve referenced from a GEOMETRIC_CURVE_SET into polyline points.
fn sample_standalone_curve(parsed: &ParsedStep, curve_ref: u64) -> Vec<glam::DVec3> {
    // Handle TRIMMED_CURVE wrapper
    if let Some(&(underlying_ref, t0, t1)) = parsed.trimmed_curves.get(&curve_ref) {
        return sample_trimmed_curve_geom(parsed, underlying_ref, t0, t1);
    }
    // Handle bare B_SPLINE_CURVE
    if let Some(control_refs) = parsed.b_spline_curves.get(&curve_ref) {
        return sample_bspline_polyline(parsed, control_refs);
    }
    // Handle bare LINE (t 0..1)
    if parsed.lines.contains_key(&curve_ref) {
        let p0 = eval_line_at(parsed, curve_ref, 0.0);
        let p1 = eval_line_at(parsed, curve_ref, 1.0);
        return match (p0, p1) {
            (Some(a), Some(b)) => vec![a, b],
            _ => Vec::new(),
        };
    }
    Vec::new()
}

/// Sample the underlying geometry of a TRIMMED_CURVE at [t0, t1].
fn sample_trimmed_curve_geom(
    parsed: &ParsedStep,
    curve_ref: u64,
    t0: f64,
    t1: f64,
) -> Vec<glam::DVec3> {
    if parsed.lines.contains_key(&curve_ref) {
        // LINE: p(t) = origin + dir * magnitude * t
        let p0 = eval_line_at(parsed, curve_ref, t0);
        let p1 = eval_line_at(parsed, curve_ref, t1);
        return match (p0, p1) {
            (Some(a), Some(b)) => vec![a, b],
            _ => Vec::new(),
        };
    }

    if let Some(&(placement_ref, radius)) = parsed.circles.get(&curve_ref) {
        // CIRCLE: HFSS exports parameters in degrees (0..360 is full circle)
        return sample_standalone_circle(parsed, placement_ref, radius, t0, t1);
    }

    Vec::new()
}

/// Sample a CIRCLE arc from t0_deg to t1_deg (degrees, HFSS convention).
fn sample_standalone_circle(
    parsed: &ParsedStep,
    placement_ref: u64,
    radius: f64,
    t0_deg: f64,
    t1_deg: f64,
) -> Vec<glam::DVec3> {
    if !radius.is_finite() || radius <= 0.0 {
        return Vec::new();
    }
    let Some((center, axis, u)) = placement_frame_from_ref(parsed, placement_ref) else {
        return Vec::new();
    };
    let v = axis.cross(u).normalize_or_zero();

    let t0 = t0_deg.to_radians();
    let mut sweep = (t1_deg - t0_deg).to_radians();
    if sweep.abs() < 1e-9 {
        sweep = std::f64::consts::TAU; // full circle
    } else if sweep < 0.0 {
        sweep += std::f64::consts::TAU;
    }

    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0).ceil().max(8.0) as usize;
    let mut points = Vec::with_capacity(seg + 1);
    for i in 0..=seg {
        let t = t0 + sweep * (i as f64 / seg as f64);
        points.push(center + u * (radius * t.cos()) + v * (radius * t.sin()));
    }
    points
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
    const EDGE_ONLY_STEP: &str = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=CARTESIAN_POINT('',(1.,0.,0.));\n#3=VERTEX_POINT('',#1);\n#4=VERTEX_POINT('',#2);\n#5=EDGE_CURVE('',#3,#4,$,.T.);\nENDSEC;\nEND-ISO-10303-21;\n";

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
    fn triangulates_spherical_face_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        let mut face_triangles = Vec::new();
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    face_triangles.push(face.triangles.len());
                }
            }
        }

        let max_face_triangles = face_triangles.into_iter().max().unwrap_or(0);
        assert!(
            max_face_triangles >= 200,
            "expected a tessellated spherical face, got max triangles={max_face_triangles}"
        );
    }

    #[test]
    fn triangulates_toroidal_face_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        let mut face_idx = 0usize;
        let mut found_torus = false;
        let mut torus_has_triangles = false;

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let is_torus = brep
                        .geom
                        .face_surface
                        .get(face_idx)
                        .and_then(|binding| *binding)
                        .and_then(|sid| brep.geom.surfaces.get(sid))
                        .map(|s| matches!(s, Surface3::Torus(_)))
                        .unwrap_or(false);

                    if is_torus {
                        found_torus = true;
                        if !face.triangles.is_empty() {
                            torus_has_triangles = true;
                        }
                    }

                    face_idx += 1;
                }
            }
        }

        assert!(found_torus, "expected hfss.step to contain a toroidal face");
        assert!(
            torus_has_triangles,
            "expected toroidal faces to be triangulated"
        );
    }

    #[test]
    fn triangulates_single_edge_planar_faces_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let mut found = false;

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    if face.outer_wire.edges.len() == 1 && face.triangles.len() >= 16 {
                        found = true;
                        break;
                    }
                }
            }
        }

        assert!(
            found,
            "expected at least one circular/elliptic single-edge planar face with triangulation"
        );
    }

    #[test]
    fn parses_box_example() {
        let brep = StepReader::parse_string(BOX_STEP).expect("box.step should parse");
        assert!(!brep.vertices.is_empty());
    }

    #[test]
    fn preserves_standalone_edge_only_geometry() {
        let brep = StepReader::parse_string(EDGE_ONLY_STEP).expect("edge-only STEP should parse");
        assert_eq!(brep.vertices.len(), 2);
        assert_eq!(brep.edges.len(), 1);
        assert!(brep.solids.is_empty(), "edge-only data should not fabricate solids");
    }

    #[test]
    fn samples_geometric_curve_sets_from_hfss() {
        // hfss.step has GEOMETRIC_CURVE_SET with Polyline1 (2 trimmed lines),
        // Polyline2 (b-spline), and Polyline3 (trimmed circle arc 0..135 deg).
        // All should produce additional edges in the BRep.
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        // The b-spline alone has 7 control points → at least 6 edges
        // The trimmed lines each contribute 1 edge (2 total)
        // The circle arc contributes 8+ edges
        // Total from curve sets: at least 16 edges beyond the face topology edges
        let total_edges = brep.edges.len();
        assert!(
            total_edges >= 20,
            "expected geometric curve set edges, got total edge count = {total_edges}"
        );
    }
}
