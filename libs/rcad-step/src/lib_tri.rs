
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

    // BREP_WITH_VOIDS: insert all shells (outer + voids) in order so that
    // the caller can detect multi-shell groupings from parsed.brep_with_voids.
    if !parsed.brep_with_voids.is_empty() {
        for (_, (outer, voids)) in &parsed.brep_with_voids {
            if let Some(face_ids) = parsed.closed_shells.get(outer) {
                shells.push(face_ids.clone());
            }
            for void_ref in voids {
                if let Some(face_ids) = parsed.closed_shells.get(void_ref) {
                    shells.push(face_ids.clone());
                }
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

fn collect_used_vertices(
    parsed: &ParsedStep,
    shell_face_sets: &[Vec<u64>],
) -> Result<BTreeSet<u64>, StepError> {
    let mut used = BTreeSet::new();

    for shell in shell_face_sets {
        for face_id in shell {
            let bound_ids = parsed
                .advanced_faces
                .get(face_id)
                .ok_or(StepError::MissingEntity {
                    entity_type: "ADVANCED_FACE",
                    id: Some(*face_id),
                })?;
            for bound_id in &bound_ids.bounds {
                let (loop_id, _) = parsed
                    .face_bounds
                    .get(bound_id)
                    .ok_or(StepError::MissingEntity {
                        entity_type: "FACE_BOUND",
                        id: Some(*bound_id),
                    })?;
                if let Some(oriented_ids) = parsed.edge_loops.get(loop_id) {
                    for oriented_id in oriented_ids {
                        let (edge_curve_id, _) =
                            parsed
                                .oriented_edges
                                .get(oriented_id)
                                .ok_or(StepError::MissingEntity {
                                    entity_type: "ORIENTED_EDGE",
                                    id: Some(*oriented_id),
                                })?;
                        let (start, end, _, _) =
                            parsed
                                .edge_curves
                                .get(edge_curve_id)
                                .ok_or(StepError::MissingEntity {
                                    entity_type: "EDGE_CURVE",
                                    id: Some(*edge_curve_id),
                                })?;
                        used.insert(*start);
                        used.insert(*end);
                    }
                } else if let Some(vp_id) = parsed.vertex_loops.get(loop_id) {
                    used.insert(*vp_id);
                } else {
                    return Err(StepError::MissingEntity {
                        entity_type: "EDGE_LOOP",
                        id: Some(*loop_id),
                    });
                }
            }
        }
    }

    Ok(used)
}

#[allow(clippy::too_many_arguments)]
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
    // Determine outer loop from FACE_OUTER_BOUND marker when available.
    // Some STEP writers do not guarantee bound list ordering.
    let outer_bound = bound_ids
        .bounds
        .iter()
        .copied()
        .find(|bid| {
            parsed
                .face_bounds
                .get(bid)
                .map(|(_, is_outer)| *is_outer)
                .unwrap_or(false)
        })
        .unwrap_or(*bound_ids.bounds.first()?);
    let (loop_id, _) = *parsed.face_bounds.get(&outer_bound)?;

    // Single-vertex outer bound (e.g. spherical face in interchange STEP writer).
    if let Some(&_vp_id) = parsed.vertex_loops.get(&loop_id) {
        let face_surface = bound_ids.surface;
        let triangles = if let Some(surface_ref) = face_surface {
            triangulate_surface_fallback(
                parsed,
                surface_ref,
                vertices,
                &[],
                false,
            )
        } else {
            Vec::new()
        };
        return Some((
            Face {
                outer_wire: Wire { edges: vec![] },
                inner_wires: Vec::new(),
                normal: glam::DVec3::new(0.0, 0.0, 1.0),
                triangles,
                sample_point: None,
                mesh_dirty: false,
                surface_idx: None,
            },
            face_surface,
        ));
    }

    let oriented_ids = parsed.edge_loops.get(&loop_id)?;

    let mut polygon: Vec<usize> = Vec::new();
    let mut wire_edge_indices: Vec<WireEdge> = Vec::new();
    let mut face_vertex_indices: Vec<usize> = Vec::new();
    let mut sampled_loop_points: Vec<glam::DVec3> = Vec::new();
    let mut sampled_loop_uv_points: Vec<glam::DVec2> = Vec::new();

    // Detect seam edges: an edge_curve that appears twice in the same face boundary
    let mut edge_curve_count: HashMap<u64, usize> = HashMap::new();
    for oriented_id in oriented_ids {
        let (edge_curve_id, _) = *parsed.oriented_edges.get(oriented_id)?;
        *edge_curve_count.entry(edge_curve_id).or_insert(0) += 1;
    }
    let has_seam = edge_curve_count.values().any(|&c| c >= 2);

    for oriented_id in oriented_ids {
        let (edge_curve_id, orientation) = *parsed.oriented_edges.get(oriented_id)?;
        let (start_id, end_id, curve_ref, same_sense) = *parsed.edge_curves.get(&edge_curve_id)?;

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
        if let Some(surface_ref) = bound_ids.surface
            && let Some(points) = sample_oriented_edge_uv_points(
                parsed,
                edge_curve_id,
                orientation,
                surface_ref,
            )
        {
            append_edge_uv_points(&mut sampled_loop_uv_points, &points);
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

        let edge_index = ensure_edge_index(
            parsed,
            start_id,
            end_id,
            curve_ref,
            same_sense,
            vertices,
            vertex_index_by_id,
            edges,
            edge_index_by_curve,
            geom,
            curve_store_index_by_step,
            edge_curve_id,
        )?;
        wire_edge_indices.push(WireEdge {
            idx: edge_index,
            forward: orientation,
        });
    }

    // Import all non-outer bounds as inner wires (holes).
    let mut inner_wires: Vec<Wire> = Vec::new();
    for inner_bound in bound_ids.bounds.iter().copied().filter(|bid| *bid != outer_bound) {
        let Some((inner_loop_id, _)) = parsed.face_bounds.get(&inner_bound).copied() else {
            continue;
        };
        let Some(inner_oriented_ids) = parsed.edge_loops.get(&inner_loop_id) else {
            continue;
        };
        let mut inner_edges = Vec::new();
        for oriented_id in inner_oriented_ids {
            let (edge_curve_id, orientation) = *parsed.oriented_edges.get(oriented_id)?;
            let (start_id, end_id, curve_ref, same_sense) = *parsed.edge_curves.get(&edge_curve_id)?;
            let edge_index = ensure_edge_index(
                parsed,
                start_id,
                end_id,
                curve_ref,
                same_sense,
                vertices,
                vertex_index_by_id,
                edges,
                edge_index_by_curve,
                geom,
                curve_store_index_by_step,
                edge_curve_id,
            )?;
            inner_edges.push(WireEdge {
                idx: edge_index,
                forward: orientation,
            });
        }
        if !inner_edges.is_empty() {
            inner_wires.push(Wire { edges: inner_edges });
        }
    }

    while polygon.len() > 1 && polygon.first() == polygon.last() {
        polygon.pop();
    }
    dedup_consecutive(&mut polygon);

    let face_surface = bound_ids.surface;
    let is_planar_face = face_surface
        .map(|sid| parsed.planes.contains_key(&sid))
        .unwrap_or(true);

    if let Some(surface_ref) = face_surface
        && sampled_loop_uv_points.len() < 3
        && sampled_loop_points.len() >= 3
        && let Some(projected) = project_boundary_points_to_surface_uv(
            parsed,
            surface_ref,
            &sampled_loop_points,
        ) {
            sampled_loop_uv_points = projected;
        }

    let triangles = if sampled_loop_points.len() >= 3 && is_planar_face {
        triangulate_point_loop(vertices, &sampled_loop_points)
    } else if let Some(surface_ref) = face_surface
        && sampled_loop_uv_points.len() >= 3
    {
        let tris = triangulate_surface_trim_loop(parsed, surface_ref, vertices, &sampled_loop_uv_points);
        eprintln!(
            "[rcad-step][diag] face #{face_id} curved uv-loop surface=#{surface_ref} uv_pts={} tris={}",
            sampled_loop_uv_points.len(),
            tris.len()
        );
        tris
    } else if let Some(surface_ref) = face_surface
        && !parsed.planes.contains_key(&surface_ref)
    {
        let tris = triangulate_surface_fallback(
            parsed,
            surface_ref,
            vertices,
            &face_vertex_indices,
            has_seam,
        );
        eprintln!(
            "[rcad-step][diag] face #{face_id} curved fallback surface=#{surface_ref} uv_pts={} sampled_pts={} tris={}",
            sampled_loop_uv_points.len(),
            sampled_loop_points.len(),
            tris.len()
        );
        tris
    } else if polygon.len() >= 3 {
        triangulate_fan(&polygon)
    } else if let Some(surface_ref) = face_surface {
        let tris = triangulate_surface_fallback(
            parsed,
            surface_ref,
            vertices,
            &face_vertex_indices,
            has_seam,
        );
        eprintln!(
            "[rcad-step][diag] face #{face_id} fallback-no-polygon surface=#{surface_ref} uv_pts={} sampled_pts={} tris={}",
            sampled_loop_uv_points.len(),
            sampled_loop_points.len(),
            tris.len()
        );
        tris
    } else {
        Vec::new()
    };
    let force_rebuild = triangles.is_empty() || !inner_wires.is_empty() || is_planar_face;
    let normal = if force_rebuild {
        estimate_loop_normal(&sampled_loop_points).unwrap_or(glam::DVec3::new(0.0, 0.0, 1.0))
    } else {
        triangles
            .first()
            .map(|tri| compute_normal(tri, vertex_index_by_id, parsed))
            .unwrap_or(glam::DVec3::new(0.0, 0.0, 1.0))
    };
    let parser_triangles = if force_rebuild { Vec::new() } else { triangles };

    Some((
        Face {
            outer_wire: Wire {
                edges: wire_edge_indices,
            },
            inner_wires: inner_wires.clone(),
            normal,
            triangles: parser_triangles,
            // Stable transitional mode:
            // - keep parser triangles for curved faces that already tessellate correctly
            // - force rebuild for planar faces and hole faces
            sample_point: None,
            mesh_dirty: force_rebuild,
            surface_idx: None,
        },
        bound_ids.surface,
    ))
}

fn estimate_loop_normal(points: &[glam::DVec3]) -> Option<glam::DVec3> {
    if points.len() < 3 {
        return None;
    }
    // Newell normal estimate for a 3D boundary loop.
    let mut n = glam::DVec3::ZERO;
    for i in 0..points.len() {
        let a = points[i];
        let b = points[(i + 1) % points.len()];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    let len2 = n.length_squared();
    if len2 <= 1e-20 {
        None
    } else {
        Some(n / len2.sqrt())
    }
}

#[allow(clippy::too_many_arguments)]
fn ensure_edge_index(
    parsed: &ParsedStep,
    start_id: u64,
    end_id: u64,
    curve_ref: Option<u64>,
    same_sense: bool,
    vertices: &[Vertex],
    vertex_index_by_id: &HashMap<u64, usize>,
    edges: &mut Vec<Edge>,
    edge_index_by_curve: &mut HashMap<u64, usize>,
    geom: &mut GeomStore,
    curve_store_index_by_step: &mut HashMap<u64, usize>,
    edge_curve_id: u64,
) -> Option<usize> {
    if let Some(idx) = edge_index_by_curve.get(&edge_curve_id) {
        return Some(*idx);
    }
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

    if geom.edge_curve_range.len() <= idx {
        geom.edge_curve_range.resize(idx + 1, None);
    }
    if geom.edge_degenerated.len() <= idx {
        geom.edge_degenerated.resize(idx + 1, false);
    }
    let explicit_trim_range = curve_ref.and_then(|step_curve| {
        if let Some((_, t0, t1)) = parsed.trimmed_curves.get(&step_curve) {
            return Some([*t0, *t1]);
        }
        parsed
            .surface_curves
            .get(&step_curve)
            .and_then(|(inner_ref, _, _)| {
                parsed
                    .trimmed_curves
                    .get(inner_ref)
                    .map(|(_, t0, t1)| [*t0, *t1])
            })
    });
    if let Some(Some(cidx)) = geom.edge_curve.get(idx)
        && let Some(curve) = geom.curves.get(*cidx)
    {
        let p0 = vertices
            .get(*vertex_index_by_id.get(&start_id)?)
            .map(|v| v.point);
        let p1 = vertices
            .get(*vertex_index_by_id.get(&end_id)?)
            .map(|v| v.point);
        if let (Some(p0), Some(p1)) = (p0, p1) {
            let mut t_range = explicit_trim_range.unwrap_or_else(|| match curve {
                Curve3::Line(line) => {
                    let t0 = (p0 - line.origin).dot(line.direction);
                    let t1 = (p1 - line.origin).dot(line.direction);
                    [t0, t1]
                }
                _ => curve.default_domain(),
            });
            if !same_sense {
                t_range.swap(0, 1);
            }
            geom.edge_curve_range[idx] = Some(t_range);
            let len = (p1 - p0).length();
            geom.edge_degenerated[idx] = len <= 1e-12;
        }
    }
    Some(idx)
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

fn append_edge_uv_points(dst: &mut Vec<glam::DVec2>, src: &[glam::DVec2]) {
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
    let (start_id, end_id, curve_ref, same_sense) = *parsed.edge_curves.get(&edge_curve_id)?;
    let start = vertex_point_from_ref(parsed, start_id)?;
    let end = vertex_point_from_ref(parsed, end_id)?;

    let mut points = if let Some(curve_id) = curve_ref {
        let actual_curve_id = if let Some(&(inner_ref, _, _)) = parsed.surface_curves.get(&curve_id) {
            inner_ref
        } else {
            curve_id
        };

        if let Some(&(underlying_ref, t0, t1)) = parsed.trimmed_curves.get(&actual_curve_id) {
            sample_trimmed_curve_geom(parsed, underlying_ref, t0, t1)
        } else if let Some((placement_ref, radius)) = parsed.circles.get(&actual_curve_id) {
            sample_circle_edge(parsed, *placement_ref, *radius, start, end)?
        } else if let Some((placement_ref, major, minor)) = parsed.ellipses.get(&actual_curve_id) {
            sample_ellipse_edge(parsed, *placement_ref, *major, *minor, start, end)?
        } else if let Some(points) = sample_spline_curve(parsed, actual_curve_id, None) {
            points
        } else if let Some(control_refs) = parsed.b_spline_curves.get(&actual_curve_id) {
            sample_bspline_polyline(parsed, control_refs)
        } else {
            vec![start, end]
        }
    } else {
        vec![start, end]
    };

    // Keep sampled geometry anchored to topological edge endpoints.
    if !points.is_empty() {
        points[0] = start;
        let last = points.len() - 1;
        points[last] = end;
    }

    if orientation != same_sense {
        points.reverse();
    }
    Some(points)
}

fn sample_oriented_edge_uv_points(
    parsed: &ParsedStep,
    edge_curve_id: u64,
    orientation: bool,
    surface_ref: u64,
) -> Option<Vec<glam::DVec2>> {
    let (_, _, curve_ref, same_sense) = *parsed.edge_curves.get(&edge_curve_id)?;
    let step_curve_id = curve_ref?;
    let (_, pcurve_refs, _) = parsed
        .surface_curves
        .get(&step_curve_id)
        .or_else(|| {
            parsed.surface_curves.iter().find_map(|(_, value)| {
                if value.0 == step_curve_id {
                    Some(value)
                } else {
                    None
                }
            })
        })?;

    let pcurve_ref = pcurve_refs.iter().find(|pcurve_ref| {
        parsed
            .pcurves
            .get(pcurve_ref)
            .map(|(candidate_surface, _)| *candidate_surface == surface_ref)
            .unwrap_or(false)
    })?;
    let (_, def_ref) = *parsed.pcurves.get(pcurve_ref)?;
    let curve2d_ref = *parsed.definitional_reps.get(&def_ref)?;

    let mut points = sample_curve2d_points(parsed, curve2d_ref)?;
    if orientation != same_sense {
        points.reverse();
    }
    Some(points)
}

fn curve2d_default_range(curve: &Curve2d) -> [f64; 2] {
    match curve {
        Curve2d::Trimmed(tc) => {
            // For trimmed curves, the effective range is [t_min, t_max].
            [tc.t_min, tc.t_max]
        }
        Curve2d::Line(_) => [0.0, 1.0],
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => [0.0, std::f64::consts::TAU],
        Curve2d::CircleInvolute(_) => [0.0, 4.0 * std::f64::consts::PI],
        Curve2d::ArchimedeanSpiral(_) | Curve2d::LogarithmicSpiral(_) => {
            [0.0, 4.0 * std::f64::consts::PI]
        }
        Curve2d::SineWave(_) | Curve2d::BSpline(_) | Curve2d::Bezier(_)
        | Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => [0.0, 1.0],
    }
}

fn sample_curve2d_points(parsed: &ParsedStep, curve_ref: u64) -> Option<Vec<glam::DVec2>> {
    let (base_curve_ref, range) = if let Some(&(underlying_ref, t0, t1)) = parsed.trimmed_curves.get(&curve_ref) {
        (underlying_ref, Some([t0, t1]))
    } else {
        (curve_ref, None)
    };
    let curve = resolve_curve2d(parsed, base_curve_ref)?;
    let [t0, t1] = range.unwrap_or_else(|| curve2d_default_range(&curve));
    let segments = match curve {
        Curve2d::Line(_) => 1usize,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => 48usize,
        Curve2d::BSpline(_) | Curve2d::Bezier(_) => 64usize,
        _ => 32usize,
    };

    let mut points = Vec::with_capacity(segments + 1);
    if (t1 - t0).abs() < 1e-12 {
        points.push(curve.point_at(t0));
        return Some(points);
    }

    for index in 0..=segments {
        let t = t0 + (t1 - t0) * (index as f64 / segments as f64);
        points.push(curve.point_at(t));
    }
    dedup_consecutive_uv_points(&mut points);
    Some(points)
}

fn dedup_consecutive_uv_points(points: &mut Vec<glam::DVec2>) {
    if points.len() < 2 {
        return;
    }
    let mut deduped = Vec::with_capacity(points.len());
    for point in points.iter().copied() {
        if deduped
            .last()
            .map(|last: &glam::DVec2| last.distance(point) < 1e-7)
            .unwrap_or(false)
        {
            continue;
        }
        deduped.push(point);
    }
    *points = deduped;
}

fn triangulate_surface_trim_loop(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    uv_points: &[glam::DVec2],
) -> Vec<[usize; 3]> {
    let Some(surface) = resolve_surface_for_trim_ops(parsed, surface_ref) else {
        return Vec::new();
    };

    let mut loop_uv = uv_points.to_vec();
    dedup_consecutive_uv_points(&mut loop_uv);
    if loop_uv.len() > 1 {
        let first = loop_uv[0];
        let last = *loop_uv.last().unwrap_or(&first);
        if first.distance(last) < 1e-6 {
            loop_uv.pop();
        }
    }
    unwrap_periodic_uv_loop(&surface, &mut loop_uv);
    if loop_uv.len() < 3 {
        return Vec::new();
    }

    if !matches!(surface, Surface3::Plane(_))
        && let Some(triangles) = triangulate_surface_trim_grid(&surface, vertices, &loop_uv)
            && !triangles.is_empty() {
                return triangles;
            }

    let uv_poly_3d: Vec<glam::DVec3> = loop_uv
        .iter()
        .map(|uv| glam::DVec3::new(uv.x, uv.y, 0.0))
        .collect();
    let local_tris = rcad_algorithms::triangulate::triangulate_polygon(&uv_poly_3d, glam::DVec3::Z);
    if local_tris.is_empty() {
        return Vec::new();
    }

    let base = vertices.len();
    for uv in &loop_uv {
        vertices.push(Vertex {
            point: surface.point_at(uv.x, uv.y),
        });
    }

    local_tris
        .into_iter()
        .filter(|[a, b, c]| a != b && b != c && a != c)
        .map(|[a, b, c]| [base + a, base + b, base + c])
        .collect()
}

fn triangulate_surface_trim_grid(
    surface: &Surface3,
    vertices: &mut Vec<Vertex>,
    uv_loop: &[glam::DVec2],
) -> Option<Vec<[usize; 3]>> {
    if uv_loop.len() < 3 {
        return None;
    }

    let (segments_u, segments_v) = match surface {
        Surface3::Sphere(_) => (32usize, 16usize),
        Surface3::Cylinder(_) | Surface3::Cone(_) => (40usize, 12usize),
        Surface3::Torus(_) => (48usize, 24usize),
        Surface3::BSpline(_) | Surface3::Bezier(_) | Surface3::Trimmed(_) => (24usize, 24usize),
        _ => (24usize, 16usize),
    };

    let mut min_u = f64::INFINITY;
    let mut max_u = f64::NEG_INFINITY;
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    for uv in uv_loop {
        min_u = min_u.min(uv.x);
        max_u = max_u.max(uv.x);
        min_v = min_v.min(uv.y);
        max_v = max_v.max(uv.y);
    }
    if !min_u.is_finite() || !max_u.is_finite() || !min_v.is_finite() || !max_v.is_finite() {
        return None;
    }
    if (max_u - min_u).abs() < 1e-10 || (max_v - min_v).abs() < 1e-10 {
        return None;
    }

    let mut grid_indices = vec![None; (segments_u + 1) * (segments_v + 1)];
    for j in 0..=segments_v {
        let v = min_v + (max_v - min_v) * (j as f64 / segments_v as f64);
        for i in 0..=segments_u {
            let u = min_u + (max_u - min_u) * (i as f64 / segments_u as f64);
            let uv = glam::DVec2::new(u, v);
            if !point_in_polygon_2d(uv, uv_loop) && !point_near_polygon_2d(uv, uv_loop, 1e-6) {
                continue;
            }
            let idx = vertices.len();
            vertices.push(Vertex {
                point: surface.point_at(u, v),
            });
            grid_indices[j * (segments_u + 1) + i] = Some(idx);
        }
    }

    let mut triangles = Vec::new();
    for j in 0..segments_v {
        for i in 0..segments_u {
            let idx = |ii: usize, jj: usize| grid_indices[jj * (segments_u + 1) + ii];
            let uv0 = glam::DVec2::new(
                min_u + (max_u - min_u) * (i as f64 / segments_u as f64),
                min_v + (max_v - min_v) * (j as f64 / segments_v as f64),
            );
            let uv1 = glam::DVec2::new(
                min_u + (max_u - min_u) * ((i + 1) as f64 / segments_u as f64),
                min_v + (max_v - min_v) * (j as f64 / segments_v as f64),
            );
            let uv2 = glam::DVec2::new(
                min_u + (max_u - min_u) * (i as f64 / segments_u as f64),
                min_v + (max_v - min_v) * ((j + 1) as f64 / segments_v as f64),
            );
            let uv3 = glam::DVec2::new(
                min_u + (max_u - min_u) * ((i + 1) as f64 / segments_u as f64),
                min_v + (max_v - min_v) * ((j + 1) as f64 / segments_v as f64),
            );

            let Some(i0) = idx(i, j) else { continue };
            let Some(i1) = idx(i + 1, j) else { continue };
            let Some(i2) = idx(i, j + 1) else { continue };
            let Some(i3) = idx(i + 1, j + 1) else { continue };

            let p0 = vertices[i0].point;
            let p1 = vertices[i1].point;
            let p2 = vertices[i2].point;
            let p3 = vertices[i3].point;
            let diag_03 = (p3 - p0).length_squared();
            let diag_12 = (p2 - p1).length_squared();
            let use_diag_03 = if (diag_03 - diag_12).abs() <= 1e-12 {
                ((i + j) & 1) == 0
            } else {
                diag_03 < diag_12
            };

            if use_diag_03 {
                let centroid_a = (uv0 + uv2 + uv3) / 3.0;
                if point_in_polygon_2d(centroid_a, uv_loop)
                    || point_near_polygon_2d(centroid_a, uv_loop, 1e-6)
                {
                    triangles.push([i0, i2, i3]);
                }

                let centroid_b = (uv0 + uv3 + uv1) / 3.0;
                if point_in_polygon_2d(centroid_b, uv_loop)
                    || point_near_polygon_2d(centroid_b, uv_loop, 1e-6)
                {
                    triangles.push([i0, i3, i1]);
                }
            } else {
                let centroid_a = (uv0 + uv2 + uv1) / 3.0;
                if point_in_polygon_2d(centroid_a, uv_loop)
                    || point_near_polygon_2d(centroid_a, uv_loop, 1e-6)
                {
                    triangles.push([i0, i2, i1]);
                }

                let centroid_b = (uv1 + uv2 + uv3) / 3.0;
                if point_in_polygon_2d(centroid_b, uv_loop)
                    || point_near_polygon_2d(centroid_b, uv_loop, 1e-6)
                {
                    triangles.push([i1, i2, i3]);
                }
            }
        }
    }

    Some(triangles)
}

fn point_in_polygon_2d(point: glam::DVec2, polygon: &[glam::DVec2]) -> bool {
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

fn point_near_polygon_2d(point: glam::DVec2, polygon: &[glam::DVec2], tolerance: f64) -> bool {
    if polygon.len() < 2 {
        return false;
    }
    for index in 0..polygon.len() {
        let a = polygon[index];
        let b = polygon[(index + 1) % polygon.len()];
        if point_segment_distance_2d(point, a, b) <= tolerance {
            return true;
        }
    }
    false
}

fn point_segment_distance_2d(point: glam::DVec2, a: glam::DVec2, b: glam::DVec2) -> f64 {
    let ab = b - a;
    let denom = ab.length_squared();
    if denom <= 1e-20 {
        return point.distance(a);
    }
    let t = ((point - a).dot(ab) / denom).clamp(0.0, 1.0);
    point.distance(a + ab * t)
}

fn project_boundary_points_to_surface_uv(
    parsed: &ParsedStep,
    surface_ref: u64,
    boundary_points: &[glam::DVec3],
) -> Option<Vec<glam::DVec2>> {
    let surface = resolve_surface_for_trim_ops(parsed, surface_ref)?;
    let mut uv_points = Vec::with_capacity(boundary_points.len());
    for &point in boundary_points {
        uv_points.push(project_point_to_surface_uv(&surface, point)?);
    }
    unwrap_periodic_uv_loop(&surface, &mut uv_points);
    dedup_consecutive_uv_points(&mut uv_points);
    Some(uv_points)
}

fn project_point_to_surface_uv(surface: &Surface3, point: glam::DVec3) -> Option<glam::DVec2> {
    match surface {
        Surface3::Plane(plane) => {
            let u_axis = any_perpendicular(plane.normal);
            let v_axis = plane.normal.cross(u_axis).normalize_or_zero();
            let delta = point - plane.origin;
            Some(glam::DVec2::new(delta.dot(u_axis), delta.dot(v_axis)))
        }
        Surface3::Cylinder(cylinder) => {
            let axis = cylinder.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let delta = point - cylinder.origin;
            let radial = delta - axis * delta.dot(axis);
            Some(glam::DVec2::new(
                radial.dot(v_axis).atan2(radial.dot(u_axis)),
                delta.dot(axis),
            ))
        }
        Surface3::Sphere(sphere) => {
            let axis = sphere.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let radial = (point - sphere.center).normalize_or_zero();
            if radial.length_squared() <= 1e-20 {
                return None;
            }
            Some(glam::DVec2::new(
                radial.dot(v_axis).atan2(radial.dot(u_axis)),
                radial.dot(axis).clamp(-1.0, 1.0).acos(),
            ))
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis_dir();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let delta = point - cone.apex;
            let height = delta.dot(axis);
            let radial = delta - axis * height;
            Some(glam::DVec2::new(
                radial.dot(v_axis).atan2(radial.dot(u_axis)),
                cone.slant_from_axial(height),
            ))
        }
        Surface3::Torus(torus) => {
            let axis = torus.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let delta = point - torus.center;
            let planar = delta - axis * delta.dot(axis);
            if planar.length_squared() <= 1e-20 {
                return None;
            }
            let theta = planar.dot(v_axis).atan2(planar.dot(u_axis));
            let ring_dir = (u_axis * theta.cos() + v_axis * theta.sin()).normalize_or_zero();
            let ring_center = torus.center + torus.major_radius * ring_dir;
            let tube_vec = (point - ring_center).normalize_or_zero();
            Some(glam::DVec2::new(
                theta,
                tube_vec.dot(axis).atan2(tube_vec.dot(ring_dir)),
            ))
        }
        Surface3::Trimmed(trimmed) => project_point_to_surface_uv(trimmed.basis.as_ref(), point),
        _ => None,
    }
}

fn unwrap_periodic_uv_loop(surface: &Surface3, points: &mut [glam::DVec2]) {
    let period = match surface {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
            Some(std::f64::consts::TAU)
        }
        Surface3::Trimmed(trimmed) => match trimmed.basis.as_ref() {
            Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
                Some(std::f64::consts::TAU)
            }
            _ => None,
        },
        _ => None,
    };

    let Some(period) = period else {
        return;
    };
    if points.len() < 2 {
        return;
    }

    let mut offset = 0.0;
    let mut previous = points[0].x;
    for point in points.iter_mut().skip(1) {
        let raw = point.x + offset;
        let delta = raw - previous;
        if delta > period * 0.5 {
            offset -= period;
        } else if delta < -period * 0.5 {
            offset += period;
        }
        point.x += offset;
        previous = point.x;
    }
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
    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0)
        .ceil()
        .max(8.0) as usize;
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
    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0)
        .ceil()
        .max(8.0) as usize;
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

fn sample_spline_curve(
    parsed: &ParsedStep,
    curve_ref: u64,
    range: Option<[f64; 2]>,
) -> Option<Vec<glam::DVec3>> {
    let curve = resolve_curve(parsed, curve_ref)?;
    if !matches!(curve, Curve3::BSpline(_) | Curve3::Bezier(_)) {
        return None;
    }

    let [t0, t1] = range.unwrap_or_else(|| curve.default_domain());
    if !t0.is_finite() || !t1.is_finite() {
        return None;
    }

    let seg = 64usize;
    let mut points = Vec::with_capacity(seg + 1);
    if (t1 - t0).abs() < 1e-12 {
        points.push(curve.point_at(t0));
        return Some(points);
    }

    for i in 0..=seg {
        let t = t0 + (t1 - t0) * (i as f64 / seg as f64);
        points.push(curve.point_at(t));
    }
    dedup_consecutive_points(&mut points);
    Some(points)
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

    let (theta_min, theta_max) =
        infer_angular_range(vertices, face_vertex_indices, center, w, u, v, has_seam);
    let (phi_min, phi_max) =
        infer_polar_range(vertices, face_vertex_indices, center, w, *radius, has_seam);

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
    if parsed.b_spline_surfaces.contains_key(&surface_ref) {
        return triangulate_bspline_surface(parsed, surface_ref, vertices);
    }
    if parsed.spherical_surfaces.contains_key(&surface_ref) {
        return triangulate_spherical_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
    }
    if parsed.cylindrical_surfaces.contains_key(&surface_ref) {
        return triangulate_cylindrical_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
    }
    if parsed.conical_surfaces.contains_key(&surface_ref) {
        return triangulate_conical_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
    }
    if parsed.toroidal_surfaces.contains_key(&surface_ref) {
        return triangulate_toroidal_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
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
    let (t_min, t_max) =
        infer_axis_range(vertices, face_vertex_indices, origin, w, -*radius, *radius);
    let (theta_min, theta_max) =
        infer_angular_range(vertices, face_vertex_indices, origin, w, u, v, has_seam);

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

    triangulate_grid(vertices, base, height_segments, radial_segments, stride)
}

fn triangulate_conical_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, ref_radius, half_angle_rad)) =
        parsed.conical_surfaces.get(&surface_ref)
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
    let (theta_min, theta_max) =
        infer_angular_range(vertices, face_vertex_indices, origin, w, u, v, has_seam);

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

    triangulate_grid(vertices, base, height_segments, radial_segments, stride)
}

fn triangulate_toroidal_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, major_radius, minor_radius)) =
        parsed.toroidal_surfaces.get(&surface_ref)
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
    let (theta_min, theta_max) = infer_angular_range(
        vertices,
        face_vertex_indices,
        center,
        w,
        u_dir,
        v_dir,
        has_seam,
    );

    // Infer minor angle (phi) range: project boundary vertices into the minor circle plane
    let (phi_min, phi_max) = infer_torus_minor_range(
        vertices,
        face_vertex_indices,
        center,
        w,
        u_dir,
        v_dir,
        *major_radius,
        has_seam,
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

    triangulate_grid(vertices, base, major_segments, minor_segments, stride)
}

/// Infer the angular (theta) range of boundary vertices projected onto a
/// surface local frame (origin, u_dir, v_dir, axis).
/// Returns (theta_min, theta_max) in radians.
/// If `has_seam` is true, the face wraps the full period ->returns (0, 2??).
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
/// Returns (phi_min, phi_max) in [0, ??] (north-pole to south-pole).
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
#[allow(clippy::too_many_arguments)]
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
    vertices: &[Vertex],
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
            let p0 = vertices[i0].point;
            let p1 = vertices[i1].point;
            let p2 = vertices[i2].point;
            let p3 = vertices[i3].point;
            let diag_03 = (p3 - p0).length_squared();
            let diag_12 = (p2 - p1).length_squared();

            if diag_03 <= diag_12 {
                triangles.push([i0, i2, i3]);
                triangles.push([i0, i3, i1]);
            } else {
                triangles.push([i0, i2, i1]);
                triangles.push([i1, i2, i3]);
            }
        }
    }
    triangles
}

/// Triangulate a B-Spline surface by uniform (u,v) grid sampling.
fn triangulate_bspline_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
) -> Vec<[usize; 3]> {
    use rcad_kernel::geom::SurfaceEval;
    let Some(surface) = resolve_surface(parsed, surface_ref) else {
        return Vec::new();
    };
    let [u0, u1, v0, v1] = surface.default_domain();
    if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
        return Vec::new();
    }
    let nu = 20usize;
    let nv = 20usize;
    let base = vertices.len();
    for j in 0..=nv {
        let v = v0 + (v1 - v0) * (j as f64 / nv as f64);
        for i in 0..=nu {
            let u = u0 + (u1 - u0) * (i as f64 / nu as f64);
            vertices.push(Vertex {
                point: surface.point_at(u, v),
            });
        }
    }
    triangulate_grid(vertices, base, nv, nu, nu + 1)
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
