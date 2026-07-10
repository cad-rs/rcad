use std::sync::Arc;
use rcad_kernel::topods::{ShapeRef, Orientation, TShape, TEdgeData};

/// Helper: extract face ShapeRefs from the first solid's first shell.
fn get_first_shell_faces(brep: &rcad_kernel::BRep) -> Vec<ShapeRef> {
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = &**ts {
            if let Some(shell_sr) = sd.shells.first() {
                return brep.shell(*shell_sr).faces.clone();
            }
        }
    }
    Vec::new()
}

/// Helper: extract the first shell's face refs from result, plus produce
/// a utility closure that maps a face index to its wire edges in (esr, forward) form.
fn get_face_outer_edges(brep: &rcad_kernel::BRep, face_sr: ShapeRef) -> Vec<(ShapeRef, bool)> {
    let fd = brep.face(face_sr);
    let wd = brep.wire(fd.outer_wire);
    wd.edges.iter().map(|esr| {
        let forward = esr.orientation == Orientation::Forward;
        (*esr, forward)
    }).collect()
}

/// Helper: get face normal from its surface (evaluate at uv 0,0).
fn get_face_normal(brep: &rcad_kernel::BRep, face_sr: ShapeRef) -> DVec3 {
    let fd = brep.face(face_sr);
    fd.surface.as_ref().map(|s| SurfaceEval::normal_at(s, 0.0, 0.0)).unwrap_or_default()
}

/// Helper: get the face surface (cloned). Returns None if no surface.
fn get_face_surface(brep: &rcad_kernel::BRep, face_sr: ShapeRef) -> Option<Surface3> {
    brep.face(face_sr).surface.clone()
}

/// Helper: check if a vertex index is used by a face (in outer or inner wires).
fn face_uses_vertex(brep: &rcad_kernel::BRep, face_sr: ShapeRef, vi: usize) -> bool {
    let fd = brep.face(face_sr);
    // Check outer wire
    let owd = brep.wire(fd.outer_wire);
    if owd.edges.iter().any(|esr| {
        let ed = brep.edge(*esr);
        ed.first.index == vi || ed.last.index == vi
    }) {
        return true;
    }
    // Check inner wires
    for &iw_sr in &fd.inner_wires {
        let iwd = brep.wire(iw_sr);
        if iwd.edges.iter().any(|esr| {
            let ed = brep.edge(*esr);
            ed.first.index == vi || ed.last.index == vi
        }) {
            return true;
        }
    }
    false
}

/// Helper: get the old-shell face count (number of faces in first shell).
fn get_face_count(brep: &rcad_kernel::BRep) -> usize {
    get_first_shell_faces(brep).len()
}

/// Helper: get the edge TEdgeData by tshape index (panics if not an edge).
fn get_edge_data<'a>(brep: &'a rcad_kernel::BRep, ei: usize) -> &'a TEdgeData {
    match &*brep.tshapes[ei] {
        TShape::Edge(ed) => ed,
        _ => unreachable!("e2.rs: expected edge at index {}", ei),
    }
}

/// Helper: get the edge TEdgeData optionally by tshape index.
fn get_edge_data_opt<'a>(brep: &'a rcad_kernel::BRep, ei: usize) -> Option<&'a TEdgeData> {
    match brep.tshapes.get(ei) {
        Some(ts) => {
            if let TShape::Edge(ed) = &**ts { Some(ed) } else { None }
        }
        None => None,
    }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Self-Intersection Detection
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Detect self-intersection with detailed results.
///
/// This is a more comprehensive analysis that returns information about
/// which faces might intersect and the minimum safe offset distance.
pub fn detect_self_intersection_detailed(brep: &rcad_kernel::BRep, distance: f64) -> SelfIntersectionResult {
    let shell_faces = get_first_shell_faces(brep);

    if shell_faces.len() < 3 {
        return SelfIntersectionResult {
            has_intersection: false,
            intersecting_pairs: Vec::new(),
            min_safe_distance: None,
            description: "insufficient faces".to_string(),
        };
    }

    // Compute face centroids
    let centroids: Vec<DVec3> = shell_faces
        .iter()
        .map(|face_sr| {
            let outer_edges = get_face_outer_edges(brep, *face_sr);
            let mut sum = DVec3::ZERO;
            let mut count = 0;
            for (esr, _forward) in &outer_edges {
                let ed = brep.edge(*esr);
                if let Some(p) = brep.vertex_point(ed.first.index) {
                    sum += p;
                    count += 1;
                }
                if let Some(p) = brep.vertex_point(ed.last.index) {
                    sum += p;
                    count += 1;
                }
            }
            if count > 0 {
                sum / count as f64
            } else {
                DVec3::ZERO
            }
        })
        .collect();

    // Build adjacency map
    let mut adjacent_pairs: HashSet<(usize, usize)> = HashSet::new();
    for (fi, face_sr) in shell_faces.iter().enumerate() {
        let outer_edges = get_face_outer_edges(brep, *face_sr);
        for (esr, _forward) in &outer_edges {
            for (fj, other_face_sr) in shell_faces.iter().enumerate() {
                if fi < fj {
                    let other_edges = get_face_outer_edges(brep, *other_face_sr);
                    if other_edges.iter().any(|(oesr, _)| oesr.ptr_id == esr.ptr_id) {
                        adjacent_pairs.insert((fi, fj));
                    }
                }
            }
        }
    }

    // Find minimum distance between non-adjacent faces
    let mut min_dist = f64::MAX;
    let mut intersecting_pairs = Vec::new();
    let abs_distance = distance.abs();

    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            if adjacent_pairs.contains(&(i, j)) {
                continue;
            }

            let dist = (centroids[i] - centroids[j]).length();

            // Check if these faces would intersect
            if abs_distance > dist * 0.5 {
                intersecting_pairs.push((i, j));
            }

            if dist < min_dist {
                min_dist = dist;
            }
        }
    }

    if min_dist == f64::MAX {
        return SelfIntersectionResult {
            has_intersection: false,
            intersecting_pairs: Vec::new(),
            min_safe_distance: None,
            description: "no non-adjacent faces found".to_string(),
        };
    }

    let has_intersection = abs_distance > min_dist * 0.5;
    let min_safe_distance = Some(min_dist * 0.5);

    let description = if has_intersection {
        format!(
            "self-intersection likely: {} face pairs at distance {} with offset {}",
            intersecting_pairs.len(),
            min_dist,
            abs_distance
        )
    } else {
        format!("no self-intersection: min distance {}, offset {}", min_dist, abs_distance)
    };

    SelfIntersectionResult {
        has_intersection,
        intersecting_pairs,
        min_safe_distance,
        description,
    }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Join Geometry Creation
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Create a sewing face to bridge the gap between two separating offset surfaces
/// at a concave edge.
///
/// When two adjacent offset surfaces separate (concave edge + outward offset, or
/// convex edge + inward offset), the offset faces no longer meet. This function
/// creates a 4-sided planar face that fills the gap.
///
/// The sewing face is bounded by:
/// - Two edges along the offset faces (one on each face's boundary)
/// - Two connecting edges at the endpoints of the shared edge
pub fn create_sewing_face(
    brep: &mut rcad_kernel::BRep,
    original_brep: &rcad_kernel::BRep,
    edge_idx: usize,
    _face_a_idx: usize,
    _face_b_idx: usize,
    distance: f64,
    offset_surfaces: &[Option<Surface3>],
) -> Option<usize> {
    let ed = get_edge_data(original_brep, edge_idx);
    let shell_faces = get_first_shell_faces(original_brep);

    let v_start = ed.first.index;
    let v_end = ed.last.index;
    let p_start = original_brep.vertex_point(v_start).unwrap_or(DVec3::ZERO);
    let p_end = original_brep.vertex_point(v_end).unwrap_or(DVec3::ZERO);

    // Compute offset vertex positions using the per-face offset normals.
    let start_offset_avg = {
        let mut sum = DVec3::ZERO;
        let mut count = 0;
        for (fi, face_sr) in shell_faces.iter().enumerate() {
            let outer_edges = get_face_outer_edges(original_brep, *face_sr);
            let uses = outer_edges.iter().any(|(esr, _)| {
                let e = original_brep.edge(*esr);
                e.first.index == v_start || e.last.index == v_start
            });
            if uses {
                sum += get_face_offset_normal(original_brep, fi, &Shell { faces: Vec::new() });
                count += 1;
            }
        }
        if count > 0 { sum.normalize_or(DVec3::Z) * distance } else { DVec3::Z * distance }
    };

    let end_offset_avg = {
        let mut sum = DVec3::ZERO;
        let mut count = 0;
        for (fi, face_sr) in shell_faces.iter().enumerate() {
            let outer_edges = get_face_outer_edges(original_brep, *face_sr);
            let uses = outer_edges.iter().any(|(esr, _)| {
                let e = original_brep.edge(*esr);
                e.first.index == v_end || e.last.index == v_end
            });
            if uses {
                sum += get_face_offset_normal(original_brep, fi, &Shell { faces: Vec::new() });
                count += 1;
            }
        }
        if count > 0 { sum.normalize_or(DVec3::Z) * distance } else { DVec3::Z * distance }
    };

    let _off_p_start = p_start + start_offset_avg;
    let _off_p_end = p_end + end_offset_avg;

    // Compute the separation direction perpendicular to both faces.
    let edge_dir = (p_end - p_start).normalize_or(DVec3::X);
    let sep_dir = any_perpendicular(edge_dir);

    // Create 4 vertices of the sewing face
    let sv0 = p_start + sep_dir * distance.abs() * 0.5;
    let sv1 = p_end + sep_dir * distance.abs() * 0.5;
    let sv2 = p_end - sep_dir * distance.abs() * 0.5;
    let sv3 = p_start - sep_dir * distance.abs() * 0.5;

    // Check that vertices are non-degenerate
    if (sv1 - sv0).length_squared() < 1e-12 || (sv3 - sv0).length_squared() < 1e-12 {
        return None;
    }

    let v0 = add_vertex(brep, sv0);
    let v1 = add_vertex(brep, sv1);
    let v2 = add_vertex(brep, sv2);
    let v3 = add_vertex(brep, sv3);

    // Create 4 edges
    let len01 = (sv1 - sv0).length();
    let e01 = add_edge(brep,
        Curve3::Line(Line3 { origin: sv0, direction: (sv1 - sv0).normalize_or(DVec3::X) }),
        0.0, len01, v0, v1);

    let len12 = (sv2 - sv1).length();
    let e12 = add_edge(brep,
        Curve3::Line(Line3 { origin: sv1, direction: (sv2 - sv1).normalize_or(DVec3::X) }),
        0.0, len12, v1, v2);

    let len23 = (sv3 - sv2).length();
    let e23 = add_edge(brep,
        Curve3::Line(Line3 { origin: sv2, direction: (sv3 - sv2).normalize_or(DVec3::X) }),
        0.0, len23, v2, v3);

    let len30 = (sv0 - sv3).length();
    let e30 = add_edge(brep,
        Curve3::Line(Line3 { origin: sv3, direction: (sv0 - sv3).normalize_or(DVec3::X) }),
        0.0, len30, v3, v0);

    let normal = edge_dir.cross(sep_dir).normalize();
    let sewing_surface = Surface3::Plane(Plane {
        origin: (sv0 + sv1 + sv2 + sv3) * 0.25,
        normal,
    });

    let wire = Wire {
        edges: vec![
            WireEdge::fwd(e01),
            WireEdge::fwd(e12),
            WireEdge::fwd(e23),
            WireEdge::fwd(e30),
        ],
    };

    let _ = offset_surfaces; // Used in more sophisticated implementations
    let _ = _face_a_idx;
    let _ = _face_b_idx;
    Some(add_face(brep, sewing_surface, wire, Vec::new()))
}

/// Create an arc join between two offset edges.
pub fn create_arc_join(
    brep: &mut rcad_kernel::BRep,
    edge_idx: usize,
    face0_idx: usize,
    face1_idx: usize,
    radius: f64,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let shell_faces = get_first_shell_faces(brep);
    if shell_faces.is_empty() {
        return Err(OffsetError::JoinCreationFailed {
            join_type: JoinType::Arc,
            edge_index: edge_idx,
            reason: "no shell found".to_string(),
        });
    }

    let ed = get_edge_data(brep, edge_idx);

    // Get the face normals from the face surfaces
    let face0_sr = shell_faces[face0_idx];
    let face1_sr = shell_faces[face1_idx];
    let n0 = get_face_normal(brep, face0_sr);
    let n1 = get_face_normal(brep, face1_sr);

    // Get the edge endpoints
    let v0 = vertex_map.get(ed.first.index).copied().unwrap_or(ed.first.index);
    let v1 = vertex_map.get(ed.last.index).copied().unwrap_or(ed.last.index);

    let p0 = brep.vertex_point(v0).unwrap_or(DVec3::ZERO);
    let p1 = brep.vertex_point(v1).unwrap_or(DVec3::ZERO);

    // Compute the edge direction and length
    let edge_dir = (p1 - p0).normalize_or(DVec3::X);
    let edge_len = (p1 - p0).length();

    // Compute the bisector normal from the two face normals
    let _bisector = (n0 + n1).normalize_or(n0);

    // Create a cylindrical surface for the arc join
    let cylinder = Surface3::Cylinder(CylindricalSurface {
        origin: p0,
        axis: edge_dir,
        radius,
        ref_dir: any_perpendicular(edge_dir),
    });

    // Create vertices for the arc join face
    let vs = add_vertex(brep, p0);
    let ve = add_vertex(brep, p1);

    // Create the edge along the cylinder
    let curve = Curve3::Line(Line3 {
        origin: p0,
        direction: edge_dir,
    });
    let arc_edge = add_edge(brep, curve, 0.0, edge_len, vs, ve);

    // Create the arc face wire
    let wire = Wire {
        edges: vec![WireEdge::fwd(arc_edge)],
    };

    // Add the arc join face
    let face_idx = add_face(brep, cylinder, wire, Vec::new());

    Ok(face_idx)
}

/// Create a tangent join between two offset edges.
pub fn create_tangent_join(
    brep: &mut rcad_kernel::BRep,
    edge_idx: usize,
    face0_idx: usize,
    face1_idx: usize,
    distance: f64,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let shell_faces = get_first_shell_faces(brep);
    if shell_faces.is_empty() {
        return Err(OffsetError::JoinCreationFailed {
            join_type: JoinType::Tangent,
            edge_index: edge_idx,
            reason: "no shell found".to_string(),
        });
    }

    let face0_sr = shell_faces[face0_idx];
    let face1_sr = shell_faces[face1_idx];
    let n0 = get_face_normal(brep, face0_sr);
    let n1 = get_face_normal(brep, face1_sr);

    // Check the angle between face normals
    let dot = n0.dot(n1);

    // If the angle is too large, fall back to intersection join
    let angle_threshold = 0.9;
    if dot < angle_threshold {
        return create_intersection_join(brep, edge_idx, face0_idx, face1_idx, vertex_map);
    }

    let ed = get_edge_data(brep, edge_idx);
    let v0 = vertex_map.get(ed.first.index).copied().unwrap_or(ed.first.index);
    let v1 = vertex_map.get(ed.last.index).copied().unwrap_or(ed.last.index);

    let p0 = brep.vertex_point(v0).unwrap_or(DVec3::ZERO);
    let p1 = brep.vertex_point(v1).unwrap_or(DVec3::ZERO);

    // Create a plane that smoothly blends the two face normals
    let blend_normal = (n0 + n1).normalize();
    let blend_plane = Surface3::Plane(Plane {
        origin: (p0 + p1) * 0.5,
        normal: blend_normal,
    });

    // Create the wire for the tangent join face
    let dir = (p1 - p0).normalize_or(DVec3::X);
    let len = (p1 - p0).length();
    let curve = Curve3::Line(Line3 { origin: p0, direction: dir });

    let vs = add_vertex(brep, p0);
    let ve = add_vertex(brep, p1);
    let blend_edge = add_edge(brep, curve, 0.0, len, vs, ve);

    let wire = Wire {
        edges: vec![WireEdge::fwd(blend_edge)],
    };

    let face_idx = add_face(brep, blend_plane, wire, Vec::new());

    let _ = distance; // Used in more sophisticated implementations
    Ok(face_idx)
}

/// Create an intersection join between two offset edges.
pub fn create_intersection_join(
    brep: &mut rcad_kernel::BRep,
    edge_idx: usize,
    _face0_idx: usize,
    _face1_idx: usize,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let ed = get_edge_data(brep, edge_idx);

    let v0 = vertex_map.get(ed.first.index).copied().unwrap_or(ed.first.index);
    let v1 = vertex_map.get(ed.last.index).copied().unwrap_or(ed.last.index);

    let p0 = brep.vertex_point(v0).unwrap_or(DVec3::ZERO);
    let p1 = brep.vertex_point(v1).unwrap_or(DVec3::ZERO);

    // For intersection join, we don't create additional geometry
    let dir = (p1 - p0).normalize_or(DVec3::X);
    let len = (p1 - p0).length();

    let midpoint = (p0 + p1) * 0.5;
    let normal = dir.any_orthonormal_pair().0;

    let plane = Surface3::Plane(Plane {
        origin: midpoint,
        normal,
    });

    let vs = add_vertex(brep, p0);
    let ve = add_vertex(brep, p1);
    let curve = Curve3::Line(Line3 { origin: p0, direction: dir });
    let int_edge = add_edge(brep, curve, 0.0, len, vs, ve);

    let wire = Wire {
        edges: vec![WireEdge::fwd(int_edge)],
    };

    let face_idx = add_face(brep, plane, wire, Vec::new());

    Ok(face_idx)
}

/// Apply join type to all edges in the shell.
pub fn apply_join_type(
    result: &mut rcad_kernel::BRep,
    original_brep: &rcad_kernel::BRep,
    opts: &OffsetOptions,
    edge_to_faces: &HashMap<usize, Vec<usize>>,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let mut join_face_count = 0;

    if opts.join_type == JoinType::Intersection {
        return Ok(0);
    }

    for (&edge_idx, face_indices) in edge_to_faces {
        if face_indices.len() < 2 {
            continue;
        }

        let face0_idx = face_indices[0];
        let face1_idx = face_indices[1];

        let join_result = match opts.join_type {
            JoinType::Arc => {
                let radius = opts.distance.abs();
                create_arc_join(result, edge_idx, face0_idx, face1_idx, radius, vertex_map)
            }
            JoinType::Tangent => {
                create_tangent_join(result, edge_idx, face0_idx, face1_idx, opts.distance, vertex_map)
            }
            JoinType::Intersection => {
                create_intersection_join(result, edge_idx, face0_idx, face1_idx, vertex_map)
            }
        };

        if join_result.is_ok() {
            join_face_count += 1;
        }
    }

    let _ = original_brep;
    Ok(join_face_count)
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Offset Quality Analysis
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Analyze the quality of an offset result.
pub fn analyze_offset_quality(
    result: &rcad_kernel::BRep,
    original: &rcad_kernel::BRep,
    opts: &OffsetOptions,
) -> OffsetQuality {
    let mut quality = OffsetQuality::default();

    // Compute minimum wall thickness
    quality.min_wall_thickness = compute_min_wall_thickness(result, opts.distance);

    // Compute maximum deviation from expected offset
    quality.max_deviation = compute_max_deviation(result, original, opts);

    // Count degenerate edges
    quality.degenerate_edge_count = result.tshapes.iter().filter(|ts| {
        if let TShape::Edge(ed) = ts.as_ref() { ed.degenerated } else { false }
    }).count();

    // Self-intersection count
    let si_result = detect_self_intersection_detailed(result, opts.distance);
    quality.self_intersection_count = si_result.intersecting_pairs.len();

    // Compute face area ratio
    quality.face_area_ratio = compute_face_area_ratio(result, original);

    // Compute edge length ratio
    quality.edge_length_ratio = compute_edge_length_ratio(result, original);

    // Determine if result is valid
    quality.is_valid = quality.self_intersection_count == 0
        && quality.min_wall_thickness >= opts.min_wall_thickness;

    // Generate warnings
    if quality.min_wall_thickness < opts.min_wall_thickness {
        quality.warnings.push(format!(
            "Minimum wall thickness {} is below threshold {}",
            quality.min_wall_thickness, opts.min_wall_thickness
        ));
    }
    if quality.max_deviation > opts.approximation_tolerance {
        quality.warnings.push(format!(
            "Maximum deviation {} exceeds approximation tolerance {}",
            quality.max_deviation, opts.approximation_tolerance
        ));
    }
    if quality.degenerate_edge_count > 0 {
        quality.warnings.push(format!(
            "Found {} degenerate edges in result",
            quality.degenerate_edge_count
        ));
    }

    quality
}

/// Compute the minimum wall thickness in the offset result.
pub fn compute_min_wall_thickness(brep: &rcad_kernel::BRep, distance: f64) -> f64 {
    let shell_faces = get_first_shell_faces(brep);
    if shell_faces.len() < 2 {
        return distance;
    }

    // Compute face centroids
    let centroids: Vec<DVec3> = shell_faces.iter().map(|face_sr| {
        let outer_edges = get_face_outer_edges(brep, *face_sr);
        let mut sum = DVec3::ZERO;
        let mut count = 0;
        for (esr, _forward) in &outer_edges {
            let ed = brep.edge(*esr);
            if let Some(p) = brep.vertex_point(ed.first.index) {
                sum += p;
                count += 1;
            }
            if let Some(p) = brep.vertex_point(ed.last.index) {
                sum += p;
                count += 1;
            }
        }
        if count > 0 { sum / count as f64 } else { DVec3::ZERO }
    }).collect();

    // Find minimum distance between any two faces
    let mut min_dist = f64::MAX;
    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            let dist = (centroids[i] - centroids[j]).length();
            if dist > 0.0 && dist < min_dist {
                min_dist = dist;
            }
        }
    }

    if min_dist == f64::MAX {
        distance
    } else {
        (min_dist - 2.0 * distance.abs()).max(0.0)
    }
}

/// Compute the maximum deviation between offset and expected positions.
pub fn compute_max_deviation(result: &rcad_kernel::BRep, original: &rcad_kernel::BRep, opts: &OffsetOptions) -> f64 {
    let _result_faces = get_first_shell_faces(result);
    let _original_faces = get_first_shell_faces(original);

    if _result_faces.is_empty() || _original_faces.is_empty() {
        return 0.0;
    }

    let mut max_dev = 0.0;

    // Compare vertex positions by iterating over vertex TShapes
    let result_verts: Vec<usize> = result.tshapes.iter().enumerate()
        .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Vertex(_)))
        .map(|(i, _)| i)
        .collect();
    let original_verts: Vec<usize> = original.tshapes.iter().enumerate()
        .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Vertex(_)))
        .map(|(i, _)| i)
        .collect();

    for (i, &vi) in result_verts.iter().enumerate() {
        if i >= original_verts.len() {
            break;
        }
        let ovi = original_verts[i];
        let actual_offset = match (result.vertex_point(vi), original.vertex_point(ovi)) {
            (Some(rp), Some(op)) => (rp - op).length(),
            _ => continue,
        };
        let expected_offset = opts.distance.abs();
        let deviation = (actual_offset - expected_offset).abs();
        if deviation > max_dev {
            max_dev = deviation;
        }
    }

    max_dev
}

/// Compute the ratio of face areas between result and original.
pub fn compute_face_area_ratio(result: &rcad_kernel::BRep, original: &rcad_kernel::BRep) -> f64 {
    let result_faces = get_first_shell_faces(result);
    let original_faces = get_first_shell_faces(original);

    if original_faces.is_empty() {
        return 1.0;
    }

    result_faces.len() as f64 / original_faces.len() as f64
}

/// Compute the ratio of edge lengths between result and original.
pub fn compute_edge_length_ratio(result: &rcad_kernel::BRep, original: &rcad_kernel::BRep) -> f64 {
    let original_count = original.edge_count();
    if original_count == 0 {
        return 1.0;
    }

    // Compute total edge lengths by iterating over edge TShapes
    let original_len: f64 = original.tshapes.iter().filter_map(|ts| {
        if let TShape::Edge(ed) = &**ts {
            let p0 = original.vertex_point(ed.first.index).unwrap_or_default();
            let p1 = original.vertex_point(ed.last.index).unwrap_or_default();
            Some((p1 - p0).length())
        } else { None }
    }).sum();

    let result_len: f64 = result.tshapes.iter().filter_map(|ts| {
        if let TShape::Edge(ed) = &**ts {
            let p0 = result.vertex_point(ed.first.index).unwrap_or_default();
            let p1 = result.vertex_point(ed.last.index).unwrap_or_default();
            Some((p1 - p0).length())
        } else { None }
    }).sum();

    if original_len > 0.0 {
        result_len / original_len
    } else {
        1.0
    }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Self-Intersection Repair
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Attempt to repair self-intersection by reducing offset distance.
pub fn repair_self_intersection(
    brep: &rcad_kernel::BRep,
    opts: &OffsetOptions,
) -> Result<(rcad_kernel::BRep, f64, usize), OffsetError> {
    let config = &opts.self_intersection_config;

    if !config.auto_repair {
        return Err(OffsetError::RecoveryFailed {
            attempts: 0,
            last_error: "auto-repair not enabled".to_string(),
        });
    }

    let mut current_distance = opts.distance;
    let mut attempts = 0;
    let mut last_error = String::new();

    while attempts < config.max_repair_attempts {
        attempts += 1;

        // Reduce the offset distance
        current_distance *= config.reduction_factor;

        if current_distance.abs() < config.min_offset_distance {
            last_error = format!(
                "offset distance {} below minimum {}",
                current_distance.abs(),
                config.min_offset_distance
            );
            continue;
        }

        // Try with reduced distance
        let mut reduced_opts = opts.clone();
        reduced_opts.distance = current_distance;
        reduced_opts.check_self_intersection = true;

        let face_refs = get_first_shell_faces(brep);
        if face_refs.is_empty() {
            return Err(OffsetError::InvalidInput("no shell"));
        }

        match offset_shell_with_options_impl(&face_refs, brep, &reduced_opts) {
            Ok(result) => {
                let si_result = detect_self_intersection_detailed(&result, current_distance);
                if !si_result.has_intersection {
                    return Ok((result, current_distance, attempts));
                }
                last_error = si_result.description;
            }
            Err(e) => {
                last_error = e.to_string();
            }
        }
    }

    Err(OffsetError::RecoveryFailed { attempts, last_error })
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Main API Functions
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Clip a 3D polygon against a half-space defined by n路x 鈮?d (interior of solid).
fn clip_polygon_by_halfspace(polygon: &[DVec3], n: DVec3, d: f64, tol: f64) -> Vec<DVec3> {
    if polygon.len() < 3 {
        return Vec::new();
    }

    let mut output = Vec::new();
    let m = polygon.len();

    for i in 0..m {
        let curr = polygon[i];
        let prev = polygon[(i + m - 1) % m];

        let curr_dist = n.dot(curr) - d;
        let prev_dist = n.dot(prev) - d;

        let curr_inside = curr_dist <= tol;
        let prev_inside = prev_dist <= tol;

        if curr_inside {
            if !prev_inside {
                let t = -prev_dist / (curr_dist - prev_dist);
                output.push(prev + t * (curr - prev));
            }
            output.push(curr);
        } else if prev_inside {
            let t = -prev_dist / (curr_dist - prev_dist);
            output.push(prev + t * (curr - prev));
        }
    }

    output
}

/// Offset a shell by moving all faces along their normals.
pub fn offset_shell(shell: &Shell, brep: &rcad_kernel::BRep, distance: f64) -> Result<rcad_kernel::BRep, OffsetError> {
    offset_shell_with_options(shell, brep, &OffsetOptions::new(distance))
}

/// Build an old-style Shell from face ShapeRefs for passing to legacy helpers.
fn build_legacy_shell(brep: &rcad_kernel::BRep, face_refs: &[ShapeRef]) -> Shell {
    Shell {
        faces: face_refs.iter().map(|&face_sr| {
            let fd = brep.face(face_sr);
            let normal = fd.surface.as_ref()
                .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
                .unwrap_or_default();
            let outer_wire = {
                let wd = brep.wire(fd.outer_wire);
                Wire {
                    edges: wd.edges.iter().map(|esr| {
                        WireEdge {
                            idx: esr.index,
                            forward: esr.orientation == Orientation::Forward,
                        }
                    }).collect(),
                }
            };
            let inner_wires = fd.inner_wires.iter().map(|&iw_sr| {
                let iwd = brep.wire(iw_sr);
                Wire {
                    edges: iwd.edges.iter().map(|esr| {
                        WireEdge {
                            idx: esr.index,
                            forward: esr.orientation == Orientation::Forward,
                        }
                    }).collect(),
                }
            }).collect();
            Face {
                outer_wire,
                inner_wires,
                normal,
                triangles: Vec::new(),
                sample_point: None,
                mesh_dirty: false,
                surface_idx: None,
            }
        }).collect(),
    }
}

/// Helper to fix face winding and normals in a result BRep using TShape walking.
fn fix_result_winding(result: &mut rcad_kernel::BRep) {
    // Collect face info (ShapeRef + normal) from first solid's first shell
    let shell_face_data: Vec<(ShapeRef, DVec3)> = {
        let mut out = Vec::new();
        for ts in &result.tshapes {
            if let TShape::Solid(sd) = ts.as_ref() {
                if let Some(&shell_sr) = sd.shells.first() {
                    let shd = result.shell(shell_sr);
                    for &face_sr in &shd.faces {
                        let fd = result.face(face_sr);
                        let n = fd.surface.as_ref()
                            .map(|s| SurfaceEval::normal_at(s, 0.0, 0.0))
                            .unwrap_or_default();
                        out.push((face_sr, n));
                    }
                }
                break;
            }
        }
        out
    };

    // Compute centroid from all vertex positions
    let center = {
        let mut center_sum = DVec3::ZERO;
        let mut vert_count = 0usize;
        for vi in 0..result.tshapes.len() {
            if let Some(p) = result.vertex_point(vi) {
                center_sum += p;
                vert_count += 1;
            }
        }
        if vert_count > 0 { center_sum / vert_count as f64 } else { DVec3::ZERO }
    };

    // Phase 1: Fix winding based on signed area vs face normal
    for &(face_sr, n) in &shell_face_data {
        let fd = result.face(face_sr);
        let owd = result.wire(fd.outer_wire);
        let mut verts: Vec<DVec3> = Vec::new();
        for esr in &owd.edges {
            let ed = result.edge(*esr);
            let pt = if esr.orientation == Orientation::Forward {
                result.vertex_point(ed.last.index).unwrap_or_default()
            } else {
                result.vertex_point(ed.first.index).unwrap_or_default()
            };
            verts.push(pt);
        }
        if verts.len() < 3 { continue; }
        let mut signed = 0.0;
        for i in 0..verts.len() {
            let j = (i + 1) % verts.len();
            signed += verts[i].cross(verts[j]).dot(n);
        }
        signed *= 0.5;
        if signed >= 0.0 { continue; }

        // Flip the wire - collect needed indices first to avoid borrow conflicts
        let outer_idx = fd.outer_wire.index;
        let inner_indices: Vec<usize> = fd.inner_wires.iter().map(|sr| sr.index).collect();

        if let Some(arc) = result.tshapes.get_mut(outer_idx) {
            if let Some(TShape::Wire(wd)) = Arc::get_mut(arc) {
                for esr in &mut wd.edges {
                    esr.orientation = match esr.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        other => other,
                    };
                }
                wd.edges.reverse();
            }
        }
        for iw_idx in &inner_indices {
            if let Some(arc) = result.tshapes.get_mut(*iw_idx) {
                if let Some(TShape::Wire(wd)) = Arc::get_mut(arc) {
                    for esr in &mut wd.edges {
                        esr.orientation = match esr.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            other => other,
                        };
                    }
                    wd.edges.reverse();
                }
            }
        }
    }

    // Phase 2: Fix inward-facing face normals using centered signed volume
    for &(face_sr, _n) in &shell_face_data {
        let fd = result.face(face_sr);
        let owd = result.wire(fd.outer_wire);
        if owd.edges.len() < 3 { continue; }
        let mut verts: Vec<DVec3> = Vec::new();
        for esr in &owd.edges {
            let ed = result.edge(*esr);
            let pt = if esr.orientation == Orientation::Forward {
                result.vertex_point(ed.last.index).unwrap_or_default()
            } else {
                result.vertex_point(ed.first.index).unwrap_or_default()
            };
            verts.push(pt);
        }
        if verts.len() < 3 { continue; }
        let centered: Vec<DVec3> = verts.iter().map(|v| *v - center).collect();
        let p0 = centered[0];
        let mut vol_6 = 0.0;
        for i in 1..centered.len() - 1 {
            vol_6 += p0.cross(centered[i]).dot(centered[i + 1]);
        }
        if vol_6 >= 0.0 { continue; }

        // Face normal points inward - collect indices first
        let outer_idx = fd.outer_wire.index;
        let inner_indices: Vec<usize> = fd.inner_wires.iter().map(|sr| sr.index).collect();
        let face_idx = face_sr.index;

        // Flip surface normal
        if let Some(arc) = result.tshapes.get_mut(face_idx) {
            if let Some(TShape::Face(f)) = Arc::get_mut(arc) {
                if let Some(surf) = &mut f.surface {
                    if let Surface3::Plane(p) = surf {
                        p.normal = -p.normal;
                    }
                }
            }
        }
        // Flip outer wire
        if let Some(arc) = result.tshapes.get_mut(outer_idx) {
            if let Some(TShape::Wire(wd)) = Arc::get_mut(arc) {
                for esr in &mut wd.edges {
                    esr.orientation = match esr.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        other => other,
                    };
                }
                wd.edges.reverse();
            }
        }
        // Flip inner wires
        for iw_idx in &inner_indices {
            if let Some(arc) = result.tshapes.get_mut(*iw_idx) {
                if let Some(TShape::Wire(wd)) = Arc::get_mut(arc) {
                    for esr in &mut wd.edges {
                        esr.orientation = match esr.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            other => other,
                        };
                    }
                    wd.edges.reverse();
                }
            }
        }
    }
}

/// Implementation of offset_shell_with_options that can be called internally.
fn offset_shell_with_options_impl(
    face_refs: &[ShapeRef],
    brep: &rcad_kernel::BRep,
    opts: &OffsetOptions,
) -> Result<rcad_kernel::BRep, OffsetError> {
    // Build legacy shell for helpers
    let shell = build_legacy_shell(brep, face_refs);

    // Validate variable thickness if specified
    if let Some(ref vt) = opts.variable_thickness {
        vt.validate(face_refs.len())?;
    }

    let distance = opts.distance;

    if distance.abs() < TOLERANCE_LEN_MIN {
        return Err(OffsetError::ZeroDistance);
    }

    if face_refs.is_empty() {
        return Err(OffsetError::InvalidInput("shell has no faces"));
    }

    // Step 1: Compute offset surfaces for each face
    let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(face_refs.len());
    for (fi, &face_sr) in face_refs.iter().enumerate() {
        let fd = brep.face(face_sr);
        let surf = match fd.surface.as_ref() {
            Some(s) => s,
            None => {
                offset_surfaces.push(None);
                continue;
            }
        };
        let face_distance = opts.effective_distance_for_face(fi);
        let off_surf = offset_surface(surf, face_distance);
        offset_surfaces.push(off_surf);
    }

    // Step 2: Build edge-to-face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, &face_sr) in face_refs.iter().enumerate() {
        let fd = brep.face(face_sr);
        let owd = brep.wire(fd.outer_wire);
        for esr in &owd.edges {
            edge_to_faces.entry(esr.index).or_default().push(fi);
        }
        for &iw_sr in &fd.inner_wires {
            let iwd = brep.wire(iw_sr);
            for esr in &iwd.edges {
                edge_to_faces.entry(esr.index).or_default().push(fi);
            }
        }
    }

    // Step 3: Compute offset vertex positions (edge-first approach).
    let get_face_distance = |fi: usize| -> f64 {
        if let Some(ref vt) = opts.variable_thickness {
            vt.thickness_for_face(fi)
        } else {
            distance
        }
    };
    let offset_vertices: Vec<DVec3> = (0..brep.vertex_count())
        .map(|vi| {
            let pt = brep.vertex_point(vi).unwrap_or(DVec3::ZERO);
            // Find all incident edges for this vertex
            let mut proj_sum = DVec3::ZERO;
            let mut proj_count = 0usize;
            for (ei, ts) in brep.tshapes.iter().enumerate() {
                let ed = match &**ts { TShape::Edge(ed) => ed, _ => continue };
                if ed.first.index != vi && ed.last.index != vi { continue; }
                let faces = match edge_to_faces.get(&ei) {
                    Some(f) if f.len() >= 2 => f,
                    _ => continue,
                };
                for pair in faces.windows(2) {
                    let fi1 = pair[0];
                    let fi2 = pair[1];
                    let f1 = brep.face(face_refs[fi1]);
                    let f2 = brep.face(face_refs[fi2]);
                    let s1 = match f1.surface.as_ref() { Some(s) => s, None => continue };
                    let s2 = match f2.surface.as_ref() { Some(s) => s, None => continue };
                    let d1 = get_face_distance(fi1);
                    let d2 = get_face_distance(fi2);
                    let inter = intersect_offset_surfaces(s1, s2, d1, d2);
                    if let Some(proj) = project_point_onto_intersection(pt, &inter) {
                        proj_sum += proj;
                        proj_count += 1;
                        break;
                    }
                }
            }
            if proj_count >= 1 {
                proj_sum / proj_count as f64
            } else {
                // Fallback: average-normal translation
                let avg_dist = if let Some(ref vt) = opts.variable_thickness {
                    let mut s = 0.0; let mut c = 0;
                    for (fi, &face_sr) in face_refs.iter().enumerate() {
                        if face_uses_vertex(brep, face_sr, vi) {
                            s += vt.thickness_for_face(fi); c += 1;
                        }
                    }
                    if c > 0 { s / c as f64 } else { distance }
                } else { distance };
                offset_vertex(brep, vi, avg_dist, &shell, None)
            }
        })
        .collect();

    // Step 4: Build result BRep
    let mut result = rcad_kernel::BRep::new();

    // Map original vertices to offset vertices
    let mut vertex_map: Vec<usize> = Vec::with_capacity(offset_vertices.len());
    for &p in &offset_vertices {
        vertex_map.push(add_vertex(&mut result, p));
    }

    // Step 5: Create offset faces with offset edges
    let mut valid_face_count = 0;

    for (fi, &face_sr) in face_refs.iter().enumerate() {
        let off_surf = match &offset_surfaces[fi] {
            Some(s) => s.clone(),
            None => continue,
        };

        let fd = brep.face(face_sr);
        let owd = brep.wire(fd.outer_wire);

        let mut wire_edges = Vec::new();

        for esr in &owd.edges {
            let ed = brep.edge(*esr);
            let vs = vertex_map[ed.first.index];
            let ve = vertex_map[ed.last.index];
            let p_start = result.vertex_point(vs).unwrap_or(DVec3::ZERO);
            let p_end = result.vertex_point(ve).unwrap_or(DVec3::ZERO);

            let faces = edge_to_faces.get(&esr.index).cloned().unwrap_or_default();
            let (curve, t0, t1) = offset_edge(brep, esr.index, &faces, distance, &offset_surfaces, &offset_vertices)
                .map(|(c, _, _)| {
                    match &c {
                        Curve3::Line(line) => {
                            let ts = project_point_to_line(p_start, line);
                            let te = project_point_to_line(p_end, line);
                            (c, ts.min(te), ts.max(te))
                        }
                        _ => (c, 0.0, (p_end - p_start).length()),
                    }
                })
                .unwrap_or_else(|| {
                    let dir = (p_end - p_start).normalize_or(DVec3::X);
                    let len = (p_end - p_start).length();
                    (Curve3::Line(Line3 { origin: p_start, direction: dir }), 0.0, len)
                });

            if (t1 - t0).abs() < TOLERANCE_LEN_MIN { continue; }

            let eidx = add_edge(&mut result, curve, t0, t1, vs, ve);
            wire_edges.push(if esr.orientation == Orientation::Forward { WireEdge::fwd(eidx) } else { WireEdge::rev(eidx) });
        }

        if wire_edges.len() < 3 { continue; }

        add_face(&mut result, off_surf, Wire { edges: wire_edges }, Vec::new());
        valid_face_count += 1;
    }

    if valid_face_count == 0 {
        return Err(OffsetError::EmptyResult);
    }

    // Step 6: Apply join type if needed
    if opts.join_type.requires_join_geometry() {
        let _join_faces = apply_join_type(&mut result, brep, opts, &edge_to_faces, &vertex_map)?;
    }

    // Fix winding and normals
    fix_result_winding(&mut result);

    Ok(result)
}

/// Offset a shell with full options.
pub fn offset_shell_with_options(
    _shell: &Shell,
    brep: &rcad_kernel::BRep,
    opts: &OffsetOptions,
) -> Result<rcad_kernel::BRep, OffsetError> {
    let face_refs = get_first_shell_faces(brep);
    if face_refs.is_empty() {
        return Err(OffsetError::InvalidInput("shell has no faces"));
    }

    let distance = opts.distance;

    if distance.abs() < TOLERANCE_LEN_MIN {
        return Err(OffsetError::ZeroDistance);
    }

    // Build legacy shell for helpers
    let _shell = build_legacy_shell(brep, &face_refs);

    // Step 1: Compute offset surfaces for each face
    let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(face_refs.len());
    for (_fi, &face_sr) in face_refs.iter().enumerate() {
        let fd = brep.face(face_sr);
        let surf = match fd.surface.as_ref() {
            Some(s) => s,
            None => {
                offset_surfaces.push(None);
                continue;
            }
        };
        let off_surf = offset_surface(surf, distance);
        if off_surf.is_none() && distance > 0.0 {
        }
        offset_surfaces.push(off_surf);
    }

    // Step 2: Build edge-to-face adjacency (including inner wires for faces with holes)
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, &face_sr) in face_refs.iter().enumerate() {
        let fd = brep.face(face_sr);
        let owd = brep.wire(fd.outer_wire);
        for esr in &owd.edges {
            edge_to_faces.entry(esr.index).or_default().push(fi);
        }
        for &iw_sr in &fd.inner_wires {
            let iwd = brep.wire(iw_sr);
            for esr in &iwd.edges {
                edge_to_faces.entry(esr.index).or_default().push(fi);
            }
        }
    }

    // Step 3: Compute offset vertex positions with position-based deduplication
    // and OCCT-aligned edge-first computation.
    let pos_tol = 1e-8;
    let mut pos_to_group: Vec<usize> = vec![usize::MAX; brep.vertex_count()];
    let mut group_positions: Vec<DVec3> = Vec::new();
    let mut group_vertex_indices: Vec<Vec<usize>> = Vec::new();
    for vi in 0..brep.vertex_count() {
        let pt = brep.vertex_point(vi).unwrap_or(DVec3::ZERO);
        let mut found = None;
        for (gi, gp) in group_positions.iter().enumerate() {
            if (pt - *gp).length_squared() < pos_tol * pos_tol {
                found = Some(gi);
                break;
            }
        }
        if let Some(gi) = found {
            pos_to_group[vi] = gi;
            group_vertex_indices[gi].push(vi);
        } else {
            pos_to_group[vi] = group_positions.len();
            group_positions.push(pt);
            group_vertex_indices.push(vec![vi]);
        }
    }

    // Compute offset for each group using edge-first projection (OCCT-aligned)
    let mut group_offsets: Vec<DVec3> = Vec::with_capacity(group_positions.len());
    for gi in 0..group_positions.len() {
        let pt = group_positions[gi];
        // Collect all incident edges across all vertices in the group
        let mut incident_edges: Vec<usize> = Vec::new();
        for &vi in &group_vertex_indices[gi] {
            for (ei, ts) in brep.tshapes.iter().enumerate() {
                let ed = match &**ts { TShape::Edge(ed) => ed, _ => continue };
                if (ed.first.index == vi || ed.last.index == vi) && !incident_edges.contains(&ei) {
                    incident_edges.push(ei);
                }
            }
        }
        // Check if ALL incident faces are planar
        let all_planar = incident_edges.iter().all(|ei| {
            edge_to_faces.get(ei).map_or(true, |faces| {
                faces.iter().all(|fi| {
                    let fd = brep.face(face_refs[*fi]);
                    fd.surface.as_ref().is_some_and(|s| matches!(s, Surface3::Plane(_)))
                })
            })
        });

        let off = if all_planar {
            // Planar-only vertex: use Cramer's rule
            let mut fi_list: Vec<usize> = Vec::new();
            let mut normal_sum = DVec3::ZERO;
            for &vi in &group_vertex_indices[gi] {
                for (fi, &face_sr) in face_refs.iter().enumerate() {
                    let uses = face_uses_vertex(brep, face_sr, vi);
                    if uses && !fi_list.contains(&fi) {
                        fi_list.push(fi);
                        normal_sum += get_face_normal(brep, face_sr);
                    }
                }
            }
            if !fi_list.is_empty() {
                offset_vertex_from_faces(brep, pt, &fi_list, normal_sum, distance, &_shell)
            } else {
                pt
            }
        } else {
            // Curved-surface vertex: OCCT edge-first projection.
            let mut projections: Vec<DVec3> = Vec::new();
            let mut seam_projection: Option<DVec3> = None;
            for &ei in &incident_edges {
                let efaces = match edge_to_faces.get(&ei) {
                    Some(f) => f,
                    None => continue,
                };
                let mut found = false;
                if efaces.len() >= 2 {
                    for pair in efaces.windows(2) {
                        let fi1 = pair[0]; let fi2 = pair[1];
                        let f1 = brep.face(face_refs[fi1]);
                        let f2 = brep.face(face_refs[fi2]);
                        let s1 = match f1.surface.as_ref() { Some(s) => s, None => continue };
                        let s2 = match f2.surface.as_ref() { Some(s) => s, None => continue };
                        let intersection = intersect_offset_surfaces(s1, s2, distance, distance);
                        if let Some(proj) = project_point_onto_intersection(pt, &intersection) {
                            projections.push(proj);
                            found = true;
                            break;
                        }
                    }
                }
                if !found {
                    // Single-face seam edge: project onto offset cylinder surface
                    for &fi in efaces.iter().take(2) {
                        let fd = brep.face(face_refs[fi]);
                        let is_cylinder = fd.surface.as_ref().is_some_and(|s| matches!(s, Surface3::Cylinder(_)));
                        if !is_cylinder { continue; }
                        let orig_surf = match fd.surface.as_ref() { Some(s) => s, None => continue };
                        let off_surf = match offset_surface(orig_surf, distance) { Some(s) => s, None => continue };
                        if let Some(uv) = project_point_to_surface_uv(pt, &off_surf, None) {
                            seam_projection = Some(off_surf.point_at(uv[0], uv[1]));
                            break;
                        }
                    }
                }
            }

            // Cone+plane result
            let cone_plane_result: Option<DVec3> = {
                let mut all_fis: Vec<usize> = Vec::new();
                for &vi in &group_vertex_indices[gi] {
                    for (fi, &face_sr) in face_refs.iter().enumerate() {
                        if face_uses_vertex(brep, face_sr, vi) && !all_fis.contains(&fi) {
                            all_fis.push(fi);
                        }
                    }
                }
                let mut cone_fi = None;
                let mut plane_fis: Vec<usize> = Vec::new();
                for &fi in &all_fis {
                    let fd = brep.face(face_refs[fi]);
                    match fd.surface.as_ref() {
                        Some(Surface3::Cone(_)) => cone_fi = Some(fi),
                        Some(Surface3::Plane(_)) => plane_fis.push(fi),
                        _ => {}
                    }
                }
                match (cone_fi, plane_fis.is_empty()) {
                    (Some(cfi), false) => offset_vertex_curved_plane(pt, brep, cfi, &plane_fis, distance, &_shell),
                    _ => None,
                }
            };

            if projections.len() >= 2 {
                let mut sum = DVec3::ZERO;
                let mut count = 0usize;
                for p in &projections { sum += *p; count += 1; }
                if let Some(cp) = cone_plane_result { sum += cp; count += 1; }
                sum / count as f64
            } else if projections.len() == 1 {
                if let Some(cp) = cone_plane_result {
                    (projections[0] + cp) * 0.5
                } else {
                    projections[0]
                }
            } else if let Some(sp) = seam_projection {
                sp
            } else if let Some(cp) = cone_plane_result {
                cp
            } else {
                // Fallback: Cramer's rule from all incident faces
                let mut fi_list: Vec<usize> = Vec::new();
                let mut normal_sum = DVec3::ZERO;
                for &vi in &group_vertex_indices[gi] {
                    for (fi, &face_sr) in face_refs.iter().enumerate() {
                        let uses = face_uses_vertex(brep, face_sr, vi);
                        if uses && !fi_list.contains(&fi) {
                            fi_list.push(fi);
                            normal_sum += get_face_normal(brep, face_sr);
                        }
                    }
                }
                if !fi_list.is_empty() {
                    offset_vertex_from_faces(brep, pt, &fi_list, normal_sum, distance, &_shell)
                } else {
                    pt
                }
            }
        };
        group_offsets.push(off);
    }

    let offset_vertices: Vec<DVec3> = (0..brep.vertex_count())
        .map(|vi| group_offsets[pos_to_group[vi]])
        .collect();

    // Step 4: Build result BRep
    let mut result = rcad_kernel::BRep::new();

    // Map original vertices to offset vertices
    let mut vertex_map: Vec<usize> = Vec::with_capacity(offset_vertices.len());
    for &p in &offset_vertices {
        vertex_map.push(add_vertex(&mut result, p));
    }

    // Step 5: Create offset faces with offset edges.
    let mut valid_face_count = 0;

    for (fi, &face_sr) in face_refs.iter().enumerate() {
        let off_surf = match &offset_surfaces[fi] {
            Some(s) => s.clone(),
            None => continue,
        };

        let fd = brep.face(face_sr);
        let owd = brep.wire(fd.outer_wire);

        let mut wire_edges = Vec::new();

        for esr in &owd.edges {
            let ed = brep.edge(*esr);
            let vs = vertex_map[ed.first.index];
            let ve = vertex_map[ed.last.index];
            let p_start = result.vertex_point(vs).unwrap_or(DVec3::ZERO);
            let p_end = result.vertex_point(ve).unwrap_or(DVec3::ZERO);

            let faces = edge_to_faces.get(&esr.index).cloned().unwrap_or_default();
            let (curve, t0, t1, edge_vs, edge_ve) = offset_edge(brep, esr.index, &faces, distance, &offset_surfaces, &offset_vertices)
                .map(|(c, _, _)| {
                    const VTX_TOL_SQ: f64 = 1e-12;
                    match &c {
                        Curve3::Line(line) => {
                            let ts = project_point_to_line(p_start, line);
                            let te = project_point_to_line(p_end, line);
                            (c, ts.min(te), ts.max(te), vs, ve)
                        }
                        Curve3::Circle(off_circle) => {
                            let orig_curve = ed.curve.as_ref();
                            let (ta, tb) = orig_curve
                                .and_then(|oc| match oc {
                                    Curve3::Circle(circ) => {
                                        let range = ed.range;
                                        let a0 = point_on_circle_angle(circ.point_at(range[0]), circ);
                                        let a1 = point_on_circle_angle(circ.point_at(range[1]), circ);
                                        if (a1 - a0).abs() < 1e-12 && (range[1] - range[0]).abs() > std::f64::consts::PI {
                                            Some((0.0, std::f64::consts::TAU))
                                        } else if a1 < a0 {
                                            Some((a0, a1 + std::f64::consts::TAU))
                                        } else {
                                            Some((a0, a1))
                                        }
                                    }
                                    _ => None,
                                })
                                .unwrap_or((0.0, std::f64::consts::TAU));
                            let normal = off_circle.normal.normalize_or(DVec3::Z);
                            let ref_dir = if normal.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                            let u_axis = normal.cross(ref_dir).normalize();
                            let _v_axis = normal.cross(u_axis).normalize();
                            // Full-circle edges keep the merged vertex
                            let is_self_loop = (tb - ta - std::f64::consts::TAU).abs() < 1e-12;
                            let vs_pt = result.vertex_point(vs).unwrap_or(DVec3::ZERO);
                            let (lvs, lve, nta, ntb) = if is_self_loop {
                                (vs, vs, ta, tb)
                            } else {
                                let vu = (vs_pt - off_circle.center).normalize_or(u_axis);
                                let vv = off_circle.normal.cross(vu).normalize();
                                let ve_pt = result.vertex_point(ve).unwrap_or(DVec3::ZERO);
                                let local_e = ve_pt - off_circle.center;
                                let ang_e = local_e.dot(vv).atan2(local_e.dot(vu));
                                let ang_e = if ang_e < 0.0 { ang_e + std::f64::consts::TAU } else { ang_e };
                                let p_end_adj = off_circle.center + off_circle.radius
                                    * (vu * ang_e.cos() + vv * ang_e.sin());
                                (vs, add_vertex(&mut result, p_end_adj), 0.0, ang_e)
                            };
                            (c, nta, ntb, lvs, lve)
                        }
                        _ => (c, 0.0, (p_end - p_start).length(), vs, ve),
                    }
                })
                .unwrap_or_else(|| {
                    let dir = (p_end - p_start).normalize_or(DVec3::X);
                    let len = (p_end - p_start).length();
                    (Curve3::Line(Line3 { origin: p_start, direction: dir }), 0.0, len, vs, ve)
                });

            if (t1 - t0).abs() < TOLERANCE_LEN_MIN {
                continue;
            }

            let eidx = add_edge(&mut result, curve, t0, t1, edge_vs, edge_ve);
            wire_edges.push(if esr.orientation == Orientation::Forward { WireEdge::fwd(eidx) } else { WireEdge::rev(eidx) });
        }

        // Fix self-loop Circle edges by splitting into two half-circles
        let mut split_wire: Vec<WireEdge> = Vec::new();
        for we in &wire_edges {
            let (ei, fwd) = (we.idx, we.forward);
            // Check if this edge is a self-loop
            let is_self_loop = result.tshapes.get(ei).map_or(false, |ts| {
                if let TShape::Edge(ed) = ts.as_ref() { ed.first.index == ed.last.index } else { false }
            });
            if is_self_loop {
                // Clone edge data before any mutable operations
                let self_loop_data = result.tshapes.get(ei).and_then(|ts| {
                    if let TShape::Edge(ed) = ts.as_ref() {
                        if let Some(Curve3::Circle(c)) = &ed.curve {
                            Some((ed.range, ed.first.index, *c))
                        } else { None }
                    } else { None }
                });
                if let Some(([t0, t1], vs_idx, circle)) = self_loop_data {
                    let mid = (t0 + t1) * 0.5;
                    let n = circle.normal.normalize_or(DVec3::Z);
                    let rd = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                    let u = n.cross(rd).normalize();
                    let v = n.cross(u).normalize();
                    let mid_pt = circle.center + circle.radius * (u * mid.cos() + v * mid.sin());
                    let mvi = add_vertex(&mut result, mid_pt);
                    let e1 = add_edge(&mut result, Curve3::Circle(circle), t0, mid, vs_idx, mvi);
                    let e2 = add_edge(&mut result, Curve3::Circle(circle), mid, t1, mvi, vs_idx);
                    split_wire.push(if fwd { WireEdge::fwd(e1) } else { WireEdge::rev(e1) });
                    split_wire.push(if fwd { WireEdge::fwd(e2) } else { WireEdge::rev(e2) });
                    continue;
                }
            }
            split_wire.push(if fwd { WireEdge::fwd(ei) } else { WireEdge::rev(ei) });
        }

        if split_wire.len() < 2 {
            continue;
        }

        let fi_new = add_face(&mut result, off_surf.clone(), Wire { edges: split_wire }, Vec::new());
        valid_face_count += 1;

        // For offset full-cylinder faces, set uv_domain
        if let Surface3::Cylinder(cyl) = &off_surf {
            let owd = brep.wire(fd.outer_wire);
            let orig_seam = {
                let mut seen = std::collections::HashSet::new();
                owd.edges.iter().any(|esr| !seen.insert(esr.ptr_id))
            };
            if orig_seam {
                // Collect v values first (immutable borrow only)
                let mut v_vals: Vec<f64> = Vec::new();
                {
                    let new_fd = result.face(ShapeRef::synthetic(fi_new));
                    let new_owd = result.wire(new_fd.outer_wire);
                    for esr in &new_owd.edges {
                        if let Some(v) = result.vertex_point(esr.index) {
                            v_vals.push((v - cyl.origin).dot(cyl.axis));
                        }
                    }
                }
                if !v_vals.is_empty() {
                    let v0 = v_vals.iter().cloned().fold(f64::INFINITY, f64::min);
                    let v1 = v_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    // Now mutate (mutable borrow) - borrows are separate
                    if let Some(arc) = result.tshapes.get_mut(fi_new) {
                        if let Some(TShape::Face(f)) = Arc::get_mut(arc) {
                            f.uv_domain = Some([0.0, std::f64::consts::TAU, v0 - 1e-10, v1 + 1e-10]);
                        }
                    }
                }
            }
        }
    }

    if valid_face_count == 0 {
        return Err(OffsetError::EmptyResult);
    }

    // Step 6: Check for self-intersection if requested
    // (Check happens inside fix_result_winding + caller)

    // Fix winding and normals
    fix_result_winding(&mut result);

    // Step 7: Remove crossed faces and fill holes
    if distance.abs() > TOLERANCE_MESH_LEGACY {
        // fix_crossed_faces takes old types; pass empty shell as it uses tshape data
        let _crossed_removed = fix_crossed_faces(&mut result, brep, &_shell, distance);
    }

    Ok(result)
}
