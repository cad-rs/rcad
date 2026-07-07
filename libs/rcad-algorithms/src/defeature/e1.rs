
// =============================================================================
// ENHANCED POCKET DETECTION
// =============================================================================

/// Detect all pocket features in `solids[0].shells[0]` with enhanced classification.
///
/// This function detects both circular and rectangular pockets, classifying them
/// as through-pockets or blind-pockets based on topology analysis.
///
/// # Arguments
/// * `brep` - The B-Rep to analyze.
/// * `config` - Pocket detection configuration.
///
/// # Returns
/// A list of detected pocket features with through/blind classification.
pub fn detect_pockets(brep: &BRep, config: &PocketDetectionConfig) -> Vec<PocketFeature> {
    if config.max_diameter <= 0.0 || config.max_depth <= 0.0 {
        return Vec::new();
    }

    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    // Strategy: Find groups of connected faces that form pocket-like shapes
    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face could be part of a pocket
        let is_pocket_candidate = (config.detect_rectangular && face_plane(brep, si, shi, start).is_some())
            || (config.detect_circular
                && face_cylinder(brep, si, shi, start)
                    .map(|c| c.radius <= config.max_diameter)
                    .unwrap_or(false));

        if !is_pocket_candidate {
            continue;
        }

        // BFS to find connected pocket-like faces
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        let mut cylindrical_walls: Vec<(CylindricalSurface, usize)> = Vec::new();
        let mut planar_faces: Vec<(Plane, usize)> = Vec::new();

        while let Some(fi) = queue.pop_front() {
            group.push(fi);

            if let Some(plane) = face_plane(brep, si, shi, fi) {
                planar_faces.push((plane, fi));
            }
            if let Some(cyl) = face_cylinder(brep, si, shi, fi)
                && cyl.radius <= config.max_diameter {
                    cylindrical_walls.push((cyl, fi));
                }

            let face_edges: Vec<usize> = {
                let f = &shell.faces[fi];
                let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                for iw in &f.inner_wires {
                    es.extend(iw.edges.iter().map(|we| we.idx));
                }
                es
            };

            for ei in face_edges {
                let Some(neighbours) = edge_to_faces.get(&ei) else {
                    continue;
                };
                for &nfi in neighbours {
                    if visited[nfi] {
                        continue;
                    }
                    let is_neighbor_candidate =
                        (config.detect_rectangular && face_plane(brep, si, shi, nfi).is_some())
                            || (config.detect_circular
                                && face_cylinder(brep, si, shi, nfi)
                                    .map(|c| c.radius <= config.max_diameter)
                                    .unwrap_or(false));

                    if is_neighbor_candidate {
                        visited[nfi] = true;
                        queue.push_back(nfi);
                    }
                }
            }
        }

        // Analyze the group for pocket characteristics
        if let Some(pocket) = analyze_pocket_enhanced(
            brep,
            si,
            shi,
            &group,
            &cylindrical_walls,
            &planar_faces,
            config,
        ) {
            features.push(pocket);
        }
    }

    features
}

/// Analyze a group of faces to determine if they form an enhanced pocket.
fn analyze_pocket_enhanced(
    brep: &BRep,
    si: usize,
    shi: usize,
    group: &[usize],
    cylindrical_walls: &[(CylindricalSurface, usize)],
    planar_faces: &[(Plane, usize)],
    config: &PocketDetectionConfig,
) -> Option<PocketFeature> {
    if group.len() < 2 {
        return None;
    }

    // Collect all vertices
    let mut vertices: Vec<DVec3> = Vec::new();
    for &fi in group {
        let face = &brep.solids[si].shells[shi].faces[fi];
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                if let Some(v) = brep.vertices.get(edge.start) {
                    vertices.push(v.point);
                }
                if let Some(v) = brep.vertices.get(edge.end) {
                    vertices.push(v.point);
                }
            }
        }
    }

    if vertices.len() < 4 {
        return None;
    }

    // Compute bounding box
    let mut min_pt = vertices[0];
    let mut max_pt = vertices[0];
    for pt in &vertices[1..] {
        min_pt = min_pt.min(*pt);
        max_pt = max_pt.max(*pt);
    }

    let dimensions = max_pt - min_pt;
    let mut dims = [dimensions.x, dimensions.y, dimensions.z];
    dims.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

    let depth = dims[2]; // Smallest dimension is likely depth

    if depth < config.min_depth || depth > config.max_depth {
        return None;
    }

    let center = (min_pt + max_pt) * 0.5;

    // Determine if circular or rectangular
    let is_circular = !cylindrical_walls.is_empty();

    let (diameter, width, length) = if is_circular {
        // Use average radius from cylindrical walls
        let avg_radius: f64 = cylindrical_walls.iter().map(|(c, _)| c.radius).sum::<f64>()
            / cylindrical_walls.len() as f64;
        (avg_radius * 2.0, 0.0, 0.0)
    } else {
        (dims[0], dims[1], dims[0])
    };

    if diameter > config.max_diameter {
        return None;
    }

    // Compute normal from planar faces (likely bottom face)
    let normal = planar_faces
        .iter()
        .filter_map(|(p, _)| {
            let n = p.normal.normalize_or_zero();
            if n.length_squared() > 0.5 {
                Some(n)
            } else {
                None
            }
        })
        .next()
        .unwrap_or(DVec3::Z);

    // Determine through vs blind pocket
    let (is_through, bottom_face_index, wall_face_indices) =
        classify_pocket_type(brep, si, shi, group, cylindrical_walls, planar_faces, config);

    Some(PocketFeature {
        face_indices: group.to_vec(),
        is_recess: true,
        diameter,
        depth,
        center,
        normal,
        is_circular,
        width,
        length,
        is_through,
        bottom_face_index,
        wall_face_indices,
    })
}

/// Classify a pocket as through or blind, and identify wall/bottom faces.
fn classify_pocket_type(
    brep: &BRep,
    si: usize,
    shi: usize,
    group: &[usize],
    cylindrical_walls: &[(CylindricalSurface, usize)],
    planar_faces: &[(Plane, usize)],
    config: &PocketDetectionConfig,
) -> (bool, Option<usize>, Vec<usize>) {
    let _group_set: HashSet<usize> = group.iter().copied().collect();

    // Collect wall face indices (cylindrical faces are typically walls)
    let wall_face_indices: Vec<usize> = cylindrical_walls
        .iter()
        .map(|(_, fi)| *fi)
        .collect();

    // Find potential bottom face (planar face with normal perpendicular to walls)
    let mut bottom_face_index: Option<usize> = None;

    // For circular pockets, check if there's a planar bottom
    if !cylindrical_walls.is_empty() {
        // Get cylinder axis direction
        let cylinder_axis = cylindrical_walls[0].0.axis.normalize_or_zero();

        // Look for planar face perpendicular to cylinder axis
        for (plane, fi) in planar_faces {
            let plane_normal = plane.normal.normalize_or_zero();
            // Bottom face normal should be opposite to cylinder axis for a blind hole
            if plane_normal.dot(cylinder_axis).abs() > 0.9 {
                bottom_face_index = Some(*fi);
                break;
            }
        }
    }

    // Determine if through-pocket
    // A through-pocket has openings on both sides of the solid
    let is_through = if bottom_face_index.is_none() {
        // No bottom face found - check if pocket opens on both sides
        // This is a heuristic based on topology
        if !cylindrical_walls.is_empty() {
            let cyl = &cylindrical_walls[0].0;
            let axis = cyl.axis.normalize_or_zero();

            // Check if cylindrical walls extend across the solid
            let mut t_values: Vec<f64> = Vec::new();
            for (_, fi) in cylindrical_walls {
                let face = &brep.solids[si].shells[shi].faces[*fi];
                for we in &face.outer_wire.edges {
                    if let Some(edge) = brep.edges.get(we.idx) {
                        for &vi in &[edge.start, edge.end] {
                            if let Some(v) = brep.vertices.get(vi) {
                                let t = (v.point - cyl.origin).dot(axis);
                                t_values.push(t);
                            }
                        }
                    }
                }
            }

            if !t_values.is_empty() {
                let t_min = t_values.iter().cloned().fold(f64::INFINITY, f64::min);
                let t_max = t_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let height = t_max - t_min;

                // If height is significant and no bottom face, likely through
                height > config.max_depth * 0.5
            } else {
                false
            }
        } else {
            false
        }
    } else {
        false
    };

    (is_through, bottom_face_index, wall_face_indices)
}

// =============================================================================
// BOSS DETECTION
// =============================================================================

/// Detect all boss features in `solids[0].shells[0]`.
///
/// Bosses are protruding cylindrical or rectangular features. This function
/// identifies bosses based on geometry and normal orientation.
///
/// # Arguments
/// * `brep` - The B-Rep to analyze.
/// * `max_diameter` - Maximum boss diameter to detect.
/// * `max_height` - Maximum boss height to detect.
///
/// # Returns
/// A list of detected boss features with height analysis.
pub fn detect_bosses(brep: &BRep, max_diameter: f64, max_height: f64) -> Vec<BossFeature> {
    if max_diameter <= 0.0 || max_height <= 0.0 {
        return Vec::new();
    }

    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    // Find cylindrical bosses first
    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face is a cylinder that could be a boss
        let Some(cyl) = face_cylinder(brep, si, shi, start) else {
            continue;
        };

        if cyl.radius > max_diameter {
            continue;
        }

        // Determine if this is a boss by checking normal direction
        let face = &shell.faces[start];
        let is_boss = !is_hole_face(face, brep, &cyl);

        if !is_boss {
            continue;
        }

        // BFS to find connected cylindrical faces on the same axis
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(fi) = queue.pop_front() {
            group.push(fi);

            let face_edges: Vec<usize> = {
                let f = &shell.faces[fi];
                let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                for iw in &f.inner_wires {
                    es.extend(iw.edges.iter().map(|we| we.idx));
                }
                es
            };

            for ei in face_edges {
                let Some(neighbours) = edge_to_faces.get(&ei) else {
                    continue;
                };
                for &nfi in neighbours {
                    if visited[nfi] {
                        continue;
                    }
                    let Some(ncyl) = face_cylinder(brep, si, shi, nfi) else {
                        continue;
                    };
                    if (ncyl.radius - cyl.radius).abs() > RADIUS_TOL {
                        continue;
                    }
                    if !axes_same_line(cyl.origin, cyl.axis, ncyl.origin, ncyl.axis) {
                        continue;
                    }
                    visited[nfi] = true;
                    queue.push_back(nfi);
                }
            }
        }

        // Analyze boss geometry
        if let Some(boss) = analyze_boss_group(brep, si, shi, &group, &cyl, max_height) {
            features.push(boss);
        }
    }

    // Also detect rectangular bosses (pads)
    // These are groups of planar faces forming a protrusion
    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        let Some(_plane) = face_plane(brep, si, shi, start) else {
            continue;
        };

        // Check if this could be part of a rectangular boss
        // Look for groups of planar faces that form a protruding shape
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        let mut planar_faces: Vec<(Plane, usize)> = Vec::new();

        while let Some(fi) = queue.pop_front() {
            group.push(fi);

            if let Some(plane) = face_plane(brep, si, shi, fi) {
                planar_faces.push((plane, fi));
            }

            let face_edges: Vec<usize> = {
                let f = &shell.faces[fi];
                let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                for iw in &f.inner_wires {
                    es.extend(iw.edges.iter().map(|we| we.idx));
                }
                es
            };

            for ei in face_edges {
                let Some(neighbours) = edge_to_faces.get(&ei) else {
                    continue;
                };
                for &nfi in neighbours {
                    if visited[nfi] {
                        continue;
                    }
                    if face_plane(brep, si, shi, nfi).is_some() {
                        visited[nfi] = true;
                        queue.push_back(nfi);
                    }
                }
            }
        }

        // Analyze as potential rectangular boss
        if group.len() >= 5 {
            // Need at least top + 4 sides
            if let Some(boss) = analyze_rectangular_boss(brep, si, shi, &group, &planar_faces, max_diameter, max_height) {
                features.push(boss);
            }
        }
    }

    features
}

/// Analyze a group of cylindrical faces to determine boss properties.
fn analyze_boss_group(
    brep: &BRep,
    si: usize,
    shi: usize,
    group: &[usize],
    cyl: &CylindricalSurface,
    max_height: f64,
) -> Option<BossFeature> {
    if group.is_empty() {
        return None;
    }

    // Compute height from vertex extents
    let ax = cyl.axis.normalize_or_zero();
    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;

    for &fi in group {
        let face = &brep.solids[si].shells[shi].faces[fi];
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                for &vi in &[edge.start, edge.end] {
                    if let Some(v) = brep.vertices.get(vi) {
                        let t = (v.point - cyl.origin).dot(ax);
                        t_min = t_min.min(t);
                        t_max = t_max.max(t);
                    }
                }
            }
        }
    }

    let height = t_max - t_min;
    if height <= 0.0 || height > max_height {
        return None;
    }

    let base_center = cyl.origin + ax * t_min;
    let diameter = cyl.radius * 2.0;

    // Find top face (planar face at t_max)
    let top_face_index = find_top_face(brep, si, shi, group, ax, t_max);

    Some(BossFeature {
        face_indices: group.to_vec(),
        diameter,
        height,
        base_center,
        normal: ax,
        is_circular: true,
        width: 0.0,
        length: 0.0,
        wall_face_indices: group.to_vec(),
        top_face_index,
    })
}

/// Find the top face of a boss (planar face at the maximum extent).
fn find_top_face(
    brep: &BRep,
    _si: usize,
    _shi: usize,
    wall_faces: &[usize],
    axis: DVec3,
    t_max: f64,
) -> Option<usize> {
    let _group_set: HashSet<usize> = wall_faces.iter().copied().collect();

    // Find faces adjacent to the top edge of the cylindrical wall
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return None;
    };

    for &fi in wall_faces {
        let face = &shell.faces[fi];
        for we in &face.outer_wire.edges {
            // Check if this edge is at the top of the cylinder
            if let Some(edge) = brep.edges.get(we.idx) {
                let mid_point = if let (Some(v1), Some(v2)) =
                    (brep.vertices.get(edge.start), brep.vertices.get(edge.end))
                {
                    (v1.point + v2.point) * 0.5
                } else {
                    continue;
                };

                let t = mid_point.dot(axis);
                if (t - t_max).abs() < TOLERANCE_ABS * 10.0 {
                    // This edge is at the top - find adjacent planar faces
                    // (This would require edge-face adjacency which we'd need to compute)
                }
            }
        }
    }

    None
}

/// Analyze a group of planar faces to determine if they form a rectangular boss.
fn analyze_rectangular_boss(
    brep: &BRep,
    si: usize,
    shi: usize,
    group: &[usize],
    planar_faces: &[(Plane, usize)],
    max_diameter: f64,
    max_height: f64,
) -> Option<BossFeature> {
    if planar_faces.len() < 5 {
        return None;
    }

    // Collect all vertices
    let mut vertices: Vec<DVec3> = Vec::new();
    for &fi in group {
        let face = &brep.solids[si].shells[shi].faces[fi];
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                if let Some(v) = brep.vertices.get(edge.start) {
                    vertices.push(v.point);
                }
                if let Some(v) = brep.vertices.get(edge.end) {
                    vertices.push(v.point);
                }
            }
        }
    }

    if vertices.len() < 8 {
        return None;
    }

    // Compute bounding box
    let mut min_pt = vertices[0];
    let mut max_pt = vertices[0];
    for pt in &vertices[1..] {
        min_pt = min_pt.min(*pt);
        max_pt = max_pt.max(*pt);
    }

    let dimensions = max_pt - min_pt;
    let mut dims = [dimensions.x, dimensions.y, dimensions.z];
    dims.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

    let height = dims[2]; // Smallest dimension is height
    let length = dims[0];
    let width = dims[1];

    if height > max_height || length > max_diameter || width > max_diameter {
        return None;
    }

    let base_center = (min_pt + max_pt) * 0.5 - DVec3::Z * height * 0.5;
    let normal = DVec3::Z; // Simplified - should compute from face normals

    Some(BossFeature {
        face_indices: group.to_vec(),
        diameter: length,
        height,
        base_center,
        normal,
        is_circular: false,
        width,
        length,
        wall_face_indices: Vec::new(), // Would need more analysis to identify
        top_face_index: None,
    })
}

// =============================================================================
// FILLET AND CHAMFER DETECTION
// =============================================================================

/// Detect all fillet features in `solids[0].shells[0]`.
///
/// Fillets are identified by toroidal, spherical, or cylindrical faces with
/// small radii that connect adjacent faces smoothly.
///
/// # Arguments
/// * `brep` - The B-Rep to analyze.
/// * `max_radius` - Maximum fillet radius to detect.
///
/// # Returns
/// A list of detected fillet features.
pub fn detect_fillets(brep: &BRep, max_radius: f64) -> Vec<FilletFeature> {
    if max_radius <= 0.0 {
        return Vec::new();
    }

    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face is a fillet (torus, sphere, or small cylinder)
        let fillet_info = detect_fillet_face(brep, si, shi, start, max_radius);

        if let Some((radius, axis, sample_point)) = fillet_info {
            // BFS to find connected fillet faces with similar radius
            let mut group: Vec<usize> = Vec::new();
            let mut queue: VecDeque<usize> = VecDeque::new();
            queue.push_back(start);
            visited[start] = true;

            let mut total_radius = radius;
            let mut min_radius = radius;
            let mut max_radius_found = radius;
            let mut count = 1usize;

            while let Some(fi) = queue.pop_front() {
                group.push(fi);

                let face_edges: Vec<usize> = {
                    let f = &shell.faces[fi];
                    let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                    for iw in &f.inner_wires {
                        es.extend(iw.edges.iter().map(|we| we.idx));
                    }
                    es
                };

                for ei in face_edges {
                    let Some(neighbours) = edge_to_faces.get(&ei) else {
                        continue;
                    };
                    for &nfi in neighbours {
                        if visited[nfi] {
                            continue;
                        }
                        if let Some((nr, _, _)) = detect_fillet_face(brep, si, shi, nfi, max_radius) {
                            // Check if similar radius
                            if (nr - radius).abs() < TOLERANCE_RETRY_LADDER_COARSE {
                                visited[nfi] = true;
                                queue.push_back(nfi);
                                total_radius += nr;
                                min_radius = min_radius.min(nr);
                                max_radius_found = max_radius_found.max(nr);
                                count += 1;
                            }
                        }
                    }
                }
            }

            let avg_radius = total_radius / count as f64;
            let is_variable = (max_radius_found - min_radius) > TOLERANCE_RETRY_LADDER_COARSE;

            // Find adjacent faces
            let adjacent_faces = find_adjacent_faces(&edge_to_faces, &group);

            features.push(FilletFeature {
                face_indices: group,
                radius: avg_radius,
                sample_point,
                axis,
                is_variable,
                min_radius,
                max_radius: max_radius_found,
                adjacent_faces,
            });
        }
    }

    features
}

/// Check if a face is a fillet face and extract its properties.
fn detect_fillet_face(
    brep: &BRep,
    si: usize,
    shi: usize,
    fi: usize,
    max_radius: f64,
) -> Option<(f64, DVec3, DVec3)> {
    // Check for torus (typical fillet)
    if let Some(torus) = face_torus(brep, si, shi, fi)
        && torus.minor_radius > 0.0 && torus.minor_radius <= max_radius {
            let sample_point = torus.center;
            return Some((torus.minor_radius, torus.axis.normalize_or_zero(), sample_point));
        }

    // Check for sphere (ball-end fillet)
    if let Some(sphere) = face_sphere(brep, si, shi, fi)
        && sphere.radius > 0.0 && sphere.radius <= max_radius {
            return Some((sphere.radius, DVec3::Z, sphere.center));
        }

    // Check for cylinder with small radius (edge fillet)
    if let Some(cyl) = face_cylinder(brep, si, shi, fi)
        && cyl.radius > 0.0 && cyl.radius <= max_radius {
            // Verify this is actually a fillet and not a hole/boss
            // by checking the normal orientation
            let sample_point = cyl.origin;
            return Some((cyl.radius, cyl.axis.normalize_or_zero(), sample_point));
        }

    None
}

/// Detect all chamfer features in `solids[0].shells[0]`.
///
/// Chamfers are identified by small planar faces that connect two
/// non-parallel faces at an angle.
///
/// # Arguments
/// * `brep` - The B-Rep to analyze.
/// * `max_distance` - Maximum chamfer distance to detect.
///
/// # Returns
/// A list of detected chamfer features.
pub fn detect_chamfers(brep: &BRep, max_distance: f64) -> Vec<ChamferFeature> {
    if max_distance <= 0.0 {
        return Vec::new();
    }

    let si = 0;
    let shi: usize = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face is a chamfer (small planar face connecting two other faces)
        let chamfer_info = detect_chamfer_face(brep, si, shi, start, max_distance, &edge_to_faces);

        if let Some((distance, angle, sample_point, normal)) = chamfer_info {
            // BFS to find connected chamfer faces
            let mut group: Vec<usize> = Vec::new();
            let mut queue: VecDeque<usize> = VecDeque::new();
            queue.push_back(start);
            visited[start] = true;

            while let Some(fi) = queue.pop_front() {
                group.push(fi);

                let face_edges: Vec<usize> = {
                    let f = &shell.faces[fi];
                    let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                    for iw in &f.inner_wires {
                        es.extend(iw.edges.iter().map(|we| we.idx));
                    }
                    es
                };

                for ei in face_edges {
                    let Some(neighbours) = edge_to_faces.get(&ei) else {
                        continue;
                    };
                    for &nfi in neighbours {
                        if visited[nfi] {
                            continue;
                        }
                        if let Some((nd, na, _, _)) =
                            detect_chamfer_face(brep, si, shi, nfi, max_distance, &edge_to_faces)
                        {
                            // Check if similar chamfer
                            if (nd - distance).abs() < TOLERANCE_RETRY_LADDER_COARSE && (na - angle).abs() < 0.1 {
                                visited[nfi] = true;
                                queue.push_back(nfi);
                            }
                        }
                    }
                }
            }

            // Find adjacent faces
            let adjacent_faces = find_adjacent_faces(&edge_to_faces, &group);

            features.push(ChamferFeature {
                face_indices: group,
                distance,
                distance2: distance, // Equal for symmetric chamfer
                angle,
                sample_point,
                normal,
                adjacent_faces,
            });
        }
    }

    features
}

/// Check if a face is a chamfer face and extract its properties.
fn detect_chamfer_face(
    brep: &BRep,
    si: usize,
    shi: usize,
    fi: usize,
    max_distance: f64,
    edge_to_faces: &HashMap<usize, Vec<usize>>,
) -> Option<(f64, f64, DVec3, DVec3)> {
    // Chamfers are planar faces
    let plane = face_plane(brep, si, shi, fi)?;

    let shell = brep.solids.first()?.shells.first()?;
    let face = &shell.faces[fi];

    // Estimate chamfer size from face dimensions
    let face_area = estimate_face_area(brep, si, shi, fi);
    let chamfer_estimate = (face_area / 2.0).sqrt(); // Approximate for typical chamfer

    if chamfer_estimate <= 0.0 || chamfer_estimate > max_distance {
        return None;
    }

    // Compute chamfer angle by analyzing adjacent faces
    let mut adjacent_normals: Vec<DVec3> = Vec::new();

    for we in &face.outer_wire.edges {
        if let Some(neighbours) = edge_to_faces.get(&we.idx) {
            for &nfi in neighbours {
                if nfi == fi {
                    continue;
                }
                if let Some(nplane) = face_plane(brep, si, shi, nfi) {
                    adjacent_normals.push(nplane.normal.normalize_or_zero());
                }
            }
        }
    }

    // Estimate angle (typically 45 degrees for standard chamfers)
    let angle = if adjacent_normals.len() >= 2 {
        // Compute angle between adjacent face normals
        let dot = adjacent_normals[0].dot(adjacent_normals[1]);
        (1.0 - dot.abs()).acos() / 2.0
    } else {
        std::f64::consts::FRAC_PI_4 // Default to 45 degrees
    };

    let sample_point = get_face_sample_point(brep, si, shi, fi).unwrap_or_default();
    let normal = plane.normal.normalize_or_zero();

    Some((chamfer_estimate, angle, sample_point, normal))
}

/// Find faces adjacent to a group of faces.
fn find_adjacent_faces(
    edge_to_faces: &HashMap<usize, Vec<usize>>,
    group: &[usize],
) -> Vec<usize> {
    let group_set: HashSet<usize> = group.iter().copied().collect();
    let mut adjacent: HashSet<usize> = HashSet::new();

    for &fi in group {
        // Find all edges for this face (would need to access the face)
        // For now, iterate through all edges
        for faces in edge_to_faces.values() {
            if faces.contains(&fi) {
                for &nfi in faces {
                    if !group_set.contains(&nfi) {
                        adjacent.insert(nfi);
                    }
                }
            }
        }
    }

    adjacent.into_iter().collect()
}

// =============================================================================
// FEATURE REMOVAL WITH HEALING
// =============================================================================

/// Remove a feature from a B-Rep with automatic topology healing.
///
/// This function removes a detected feature by index, fills the resulting
/// void, and heals the surrounding geometry.
///
/// # Arguments
/// * `brep` - The B-Rep to modify.
/// * `feature_idx` - Index of the feature to remove.
/// * `feature_type` - Type of the feature to remove.
/// * `features` - The collection of detected features.
/// * `healing_tolerance` - Tolerance for post-removal healing.
///
/// # Returns
/// The modified B-Rep with the feature removed and geometry healed.
pub fn remove_feature_with_healing<F>(
    brep: &BRep,
    feature_idx: usize,
    _feature_type: FeatureType,
    features: &[F],
    healing_tolerance: f64,
) -> BRep
where
    F: FeatureToBRep,
{
    let Some(feature) = features.get(feature_idx) else {
        return brep.clone();
    };

    // Build the fill solid
    let fill_brep = feature.to_fill_brep();
    let fill_old = rcad_kernel::BRep::from_topods(&fill_brep);

    // Perform the boolean operation
    let result = if feature.is_removal_by_union() {
        boolean_op(BooleanOpType::Union, brep, &fill_old)
    } else {
        boolean_op(BooleanOpType::Difference, brep, &fill_old)
    };

    let mut result_brep = match result {
        Ok(b) => rcad_kernel::BRep::from_topods(&b),
        Err(_) => return brep.clone(),
    };

    // Apply healing
    let healing_opts = PostSuppressionHealingOptions {
        gap_tolerance: healing_tolerance,
        merge_tolerance: healing_tolerance,
        ..PostSuppressionHealingOptions::default()
    };

    let (healed, _report) = heal_after_suppression(&result_brep, &healing_opts);
    result_brep = healed;

    result_brep
}

/// Trait for converting a feature to a B-Rep for removal operations.
pub trait FeatureToBRep {
    /// Convert the feature to a B-Rep solid for filling/removal.
    fn to_fill_brep(&self) -> topods::BRep;

    /// Whether the feature is removed by union (true) or difference (false).
    fn is_removal_by_union(&self) -> bool;
}

impl FeatureToBRep for CylindricalFeature {
    fn to_fill_brep(&self) -> topods::BRep {
        make_fill_cylinder(self, DEFAULT_FILL_MARGIN).unwrap_or_default()
    }

    fn is_removal_by_union(&self) -> bool {
        self.is_hole
    }
}

impl FeatureToBRep for PocketFeature {
    fn to_fill_brep(&self) -> topods::BRep {
        // Build a fill solid for the pocket
        // For circular pockets, use a cylinder
        // For rectangular pockets, use a box
        if self.is_circular {
            let radius = self.diameter / 2.0;
            let height = self.depth + DEFAULT_FILL_MARGIN * 2.0;
            let base_pt = self.center - self.normal * (self.depth + DEFAULT_FILL_MARGIN);
            make_cylinder_brep(
                base_pt,
                self.normal,
                any_perpendicular(self.normal),
                radius + TOLERANCE_ABS * 10.0,
                height,
            )
            .unwrap_or_default()
        } else {
            // Rectangular pocket - use a box
            let height = self.depth + DEFAULT_FILL_MARGIN * 2.0;
            rcad_modeling::make_box_brep(
                self.center - DVec3::new(self.length / 2.0, self.width / 2.0, 0.0)
                    - self.normal * DEFAULT_FILL_MARGIN,
                DVec3::X,
                DVec3::Y,
                self.length,
                self.width,
                height,
            )
            .unwrap_or_default()
        }
    }

    fn is_removal_by_union(&self) -> bool {
        self.is_recess
    }
}

impl FeatureToBRep for BossFeature {
    fn to_fill_brep(&self) -> topods::BRep {
        // Build a solid representing the boss for removal
        if self.is_circular {
            let radius = self.diameter / 2.0;
            let height = self.height + DEFAULT_FILL_MARGIN * 2.0;
            let base_pt = self.base_center - self.normal * DEFAULT_FILL_MARGIN;
            make_cylinder_brep(
                base_pt,
                self.normal,
                any_perpendicular(self.normal),
                radius + TOLERANCE_ABS * 10.0,
                height,
            )
            .unwrap_or_default()
        } else {
            // Rectangular boss
            let height = self.height + DEFAULT_FILL_MARGIN * 2.0;
            rcad_modeling::make_box_brep(
                self.base_center
                    - DVec3::new(self.length / 2.0, self.width / 2.0, 0.0)
                    - self.normal * DEFAULT_FILL_MARGIN,
                DVec3::X,
                DVec3::Y,
                self.length,
                self.width,
                height,
            )
            .unwrap_or_default()
        }
    }

    fn is_removal_by_union(&self) -> bool {
        false // Bosses are removed by difference
    }
}

impl FeatureToBRep for FilletFeature {
    fn to_fill_brep(&self) -> topods::BRep {
        // Fillet removal is complex - would need to reconstruct the sharp edge
        // For now, return an empty BRep
        topods::BRep::new()
    }

    fn is_removal_by_union(&self) -> bool {
        true // Simplified - actual operation depends on geometry
    }
}

impl FeatureToBRep for ChamferFeature {
    fn to_fill_brep(&self) -> topods::BRep {
        // Chamfer removal is complex - would need to reconstruct the sharp edge
        // For now, return an empty BRep
        topods::BRep::new()
    }

    fn is_removal_by_union(&self) -> bool {
        true // Simplified - actual operation depends on geometry
    }
}

impl FeatureToBRep for SlotFeature {
    fn to_fill_brep(&self) -> topods::BRep {
        // Build a fill solid for the slot
        let height = self.depth + DEFAULT_FILL_MARGIN * 2.0;
        rcad_modeling::make_box_brep(
            self.origin - self.depth_dir * DEFAULT_FILL_MARGIN,
            self.length_dir,
            self.width_dir,
            self.length,
            self.width,
            height,
        )
        .unwrap_or_default()
    }

    fn is_removal_by_union(&self) -> bool {
        self.is_recess
    }
}

impl FeatureToBRep for BlendFeature {
    fn to_fill_brep(&self) -> topods::BRep {
        topods::BRep::new()
    }

    fn is_removal_by_union(&self) -> bool {
        true
    }
}

/// Detect all blend (fillet/chamfer) features in `solids[0].shells[0]`.
///
/// Parameters:
/// - `max_blend_radius`: Maximum fillet radius to detect
/// - `max_chamfer_distance`: Maximum chamfer distance to detect
///
/// Returns a list of [`BlendFeature`] objects.
pub fn detect_blend_features(
    brep: &BRep,
    max_blend_radius: f64,
    max_chamfer_distance: f64,
) -> Vec<BlendFeature> {
    if max_blend_radius <= 0.0 && max_chamfer_distance <= 0.0 {
        return Vec::new();
    }

    let si = 0;
    let shi = 0;
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };
    let n_faces = shell.faces.len();

    // Build edge -> face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    let mut visited = vec![false; n_faces];
    let mut features = Vec::new();

    for start in 0..n_faces {
        if visited[start] {
            continue;
        }

        // Check if this face is a blend
        if let Some(blend) = detect_blend_face(brep, si, shi, start, max_blend_radius, max_chamfer_distance) {
            // BFS to find connected blend faces with similar characteristics
            let mut group: Vec<usize> = Vec::new();
            let mut queue: VecDeque<usize> = VecDeque::new();
            queue.push_back(start);
            visited[start] = true;

            let mut total_radius = blend.radius;
            let mut count = 1usize;

            while let Some(fi) = queue.pop_front() {
                group.push(fi);

                let face_edges: Vec<usize> = {
                    let f = &shell.faces[fi];
                    let mut es: Vec<usize> = f.outer_wire.edges.iter().map(|we| we.idx).collect();
                    for iw in &f.inner_wires {
                        es.extend(iw.edges.iter().map(|we| we.idx));
                    }
                    es
                };

                for ei in face_edges {
                    let Some(neighbours) = edge_to_faces.get(&ei) else {
                        continue;
                    };
                    for &nfi in neighbours {
                        if visited[nfi] {
                            continue;
                        }
                        if let Some(nblend) = detect_blend_face(brep, si, shi, nfi, max_blend_radius, max_chamfer_distance) {
                            // Check if similar blend
                            if (nblend.is_fillet == blend.is_fillet)
                                && (nblend.radius - blend.radius).abs() < TOLERANCE_RETRY_LADDER_COARSE
                            {
                                visited[nfi] = true;
                                queue.push_back(nfi);
                                total_radius += nblend.radius;
                                count += 1;
                            }
                        }
                    }
                }
            }

            let avg_radius = if count > 0 { total_radius / count as f64 } else { blend.radius };

            features.push(BlendFeature {
                face_indices: group,
                is_fillet: blend.is_fillet,
                radius: avg_radius,
                chamfer_distance: blend.chamfer_distance,
                sample_point: blend.sample_point,
                normal: blend.normal,
            });
        }
    }

    features
}

/// Detect connected groups of features that should be processed together.
///
/// This function analyzes spatial relationships between features and groups
/// those that share edges or vertices.
///
/// Returns a map from face index to group ID.
pub fn detect_connected_feature_groups(
    brep: &BRep,
    cylindrical_features: &[CylindricalFeature],
    conical_features: &[ConicalFeature],
    slot_features: &[SlotFeature],
    pocket_features: &[PocketFeature],
    blend_features: &[BlendFeature],
) -> (Vec<FeatureGroup>, HashMap<usize, usize>) {
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return (Vec::new(), HashMap::new());
    };

    // Build face -> feature indices mapping
    let mut face_to_features: HashMap<usize, Vec<(usize, FeatureType)>> = HashMap::new();

    for (i, f) in cylindrical_features.iter().enumerate() {
        for &fi in &f.face_indices {
            face_to_features.entry(fi).or_default().push((i, FeatureType::Cylindrical));
        }
    }
    for (i, f) in conical_features.iter().enumerate() {
        for &fi in &f.face_indices {
            face_to_features.entry(fi).or_default().push((i, FeatureType::Conical));
        }
    }
    for (i, f) in slot_features.iter().enumerate() {
        for &fi in &f.face_indices {
            face_to_features.entry(fi).or_default().push((i, FeatureType::Slot));
        }
    }
    for (i, f) in pocket_features.iter().enumerate() {
        for &fi in &f.face_indices {
            face_to_features.entry(fi).or_default().push((i, FeatureType::Pocket));
        }
    }
    for (i, f) in blend_features.iter().enumerate() {
        for &fi in &f.face_indices {
            face_to_features.entry(fi).or_default().push((i, FeatureType::Blend));
        }
    }

    // Build edge adjacency between faces
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                edge_to_faces.entry(we.idx).or_default().push(fi);
            }
        }
    }

    // Build feature adjacency graph through shared edges
    let mut feature_adjacency: HashMap<(usize, FeatureType), HashSet<(usize, FeatureType)>> = HashMap::new();

    for feature_list in face_to_features.values() {
        // All features sharing a face are connected
        for i in 0..feature_list.len() {
            for j in (i + 1)..feature_list.len() {
                feature_adjacency
                    .entry(feature_list[i])
                    .or_default()
                    .insert(feature_list[j]);
                feature_adjacency
                    .entry(feature_list[j])
                    .or_default()
                    .insert(feature_list[i]);
            }
        }
    }

    // Also check edge-sharing between features
    for (fi, features1) in &face_to_features {
        // Find neighboring faces through edges
        let face = &shell.faces[*fi];
        for we in &face.outer_wire.edges {
            if let Some(neighbors) = edge_to_faces.get(&we.idx) {
                for &nfi in neighbors {
                    if nfi == *fi {
                        continue;
                    }
                    if let Some(features2) = face_to_features.get(&nfi) {
                        for f1 in features1 {
                            for f2 in features2 {
                                if f1 != f2 {
                                    feature_adjacency.entry(*f1).or_default().insert(*f2);
                                    feature_adjacency.entry(*f2).or_default().insert(*f1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Find connected components using BFS
    let mut visited: HashSet<(usize, FeatureType)> = HashSet::new();
    let mut groups: Vec<FeatureGroup> = Vec::new();
    let mut face_to_group: HashMap<usize, usize> = HashMap::new();

    let all_features: Vec<(usize, FeatureType)> = face_to_features
        .values()
        .flat_map(|v| v.iter().copied())
        .collect();

    for start in all_features {
        if visited.contains(&start) {
            continue;
        }

        let mut group = FeatureGroup {
            id: groups.len(),
            cylindrical_indices: Vec::new(),
            conical_indices: Vec::new(),
            slot_indices: Vec::new(),
            pocket_indices: Vec::new(),
            blend_indices: Vec::new(),
            total_faces: 0,
        };

        let mut queue: VecDeque<(usize, FeatureType)> = VecDeque::new();
        queue.push_back(start);
        visited.insert(start);

        while let Some((idx, ftype)) = queue.pop_front() {
            match ftype {
                FeatureType::Cylindrical => group.cylindrical_indices.push(idx),
                FeatureType::Conical => group.conical_indices.push(idx),
                FeatureType::Slot => group.slot_indices.push(idx),
                FeatureType::Pocket => group.pocket_indices.push(idx),
                FeatureType::Blend => group.blend_indices.push(idx),
                FeatureType::Boss | FeatureType::Fillet | FeatureType::Chamfer => {
                    // These feature types don't have dedicated indices in the group
                }
            }

            // Add faces to group map
            let face_indices: &Vec<usize> = match ftype {
                FeatureType::Cylindrical => &cylindrical_features[idx].face_indices,
                FeatureType::Conical => &conical_features[idx].face_indices,
                FeatureType::Slot => &slot_features[idx].face_indices,
                FeatureType::Pocket => &pocket_features[idx].face_indices,
                FeatureType::Blend => &blend_features[idx].face_indices,
                FeatureType::Boss | FeatureType::Fillet | FeatureType::Chamfer => &Vec::new(),
            };
            for &fi in face_indices {
                face_to_group.insert(fi, group.id);
                group.total_faces += 1;
            }

            // Explore neighbors
            if let Some(neighbors) = feature_adjacency.get(&(idx, ftype)) {
                for &neighbor in neighbors {
                    if !visited.contains(&neighbor) {
                        visited.insert(neighbor);
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        groups.push(group);
    }

    (groups, face_to_group)
}

/// Detect hole patterns (arrays of similar cylindrical holes) from a list of
/// cylindrical features.
///
/// This function groups cylindrical features that share similar radii and are
/// arranged in regular patterns (linear, circular, rectangular grid).
///
/// # Arguments
/// * `features` - List of cylindrical features to analyze.
/// * `radius_tolerance` - Maximum difference in radii for holes to be grouped.
/// * `spacing_tolerance` - Maximum deviation from expected pattern spacing (as fraction).
///
/// # Returns
/// A list of detected hole patterns.
pub fn detect_hole_patterns(
    features: &[CylindricalFeature],
    radius_tolerance: f64,
    spacing_tolerance: f64,
) -> Vec<HolePattern> {
    if features.len() < 2 {
        return Vec::new();
    }

    let radius_tol = radius_tolerance.max(TOLERANCE_MESH_LEGACY);
    let _spacing_tol = spacing_tolerance.max(0.01).min(0.5); // Clamp to reasonable range

    // Group features by similar radius and axis direction
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut assigned = vec![false; features.len()];

    for i in 0..features.len() {
        if assigned[i] {
            continue;
        }
        if !features[i].is_hole {
            continue;
        }

        let mut group = vec![i];
        assigned[i] = true;

        for j in (i + 1)..features.len() {
            if assigned[j] || !features[j].is_hole {
                continue;
            }

            // Check similar radius
            if (features[i].radius - features[j].radius).abs() > radius_tol {
                continue;
            }

            // Check parallel axes
            if !axes_parallel(features[i].axis, features[j].axis) {
                continue;
            }

            group.push(j);
            assigned[j] = true;
        }

        if group.len() >= 2 {
            groups.push(group);
        }
    }

    // Analyze each group for pattern type
    let mut patterns: Vec<HolePattern> = Vec::new();

    for group in groups {
        if group.len() < 2 {
            continue;
        }

        // Get centers of all holes in the group
        let centers: Vec<DVec3> = group
            .iter()
            .map(|&idx| {
                let f = &features[idx];
                f.origin + f.axis * (f.t_min + f.t_max) * 0.5
            })
            .collect();

        // Try to detect pattern type
        let pattern_type = classify_pattern_type(&centers, spacing_tolerance);

        let common_radius = features[group[0]].radius;
        let common_depth = features[group[0]].height();

        // Compute pattern properties
        let (origin, direction, spacing) = compute_pattern_properties(&centers, &pattern_type);

        patterns.push(HolePattern {
            feature_indices: group.clone(),
            pattern_type,
            count: group.len(),
            spacing,
            origin,
            direction,
            common_radius,
            common_depth,
        });
    }

    patterns
}

/// Classify the pattern type from a set of hole centers.
fn classify_pattern_type(centers: &[DVec3], spacing_tolerance: f64) -> HolePatternType {
    let n = centers.len();
    if n < 2 {
        return HolePatternType::Irregular;
    }

    // Compute centroid
    let centroid = centers.iter().fold(DVec3::ZERO, |acc, &p| acc + p) / n as f64;

    // Check for circular pattern: all points equidistant from centroid
    if n >= 3 {
        let distances: Vec<f64> = centers.iter().map(|p| (*p - centroid).length()).collect();
        let avg_dist = distances.iter().sum::<f64>() / n as f64;
        let max_deviation = distances
            .iter()
            .map(|&d| (d - avg_dist).abs())
            .fold(0.0, f64::max);

        if avg_dist > TOLERANCE_MESH_LEGACY && max_deviation / avg_dist < spacing_tolerance {
            return HolePatternType::Circular;
        }
    }

    // Check for linear pattern: points lie along a line
    if n >= 2 {
        // Compute best-fit line through points
        let direction = (centers[n - 1] - centers[0]).normalize_or_zero();
        if direction.length_squared() > 0.5 {
            let mut max_dist_from_line = 0.0f64;
            for p in centers {
                let to_point = *p - centers[0];
                let proj = to_point.dot(direction);
                let perp = to_point - proj * direction;
                max_dist_from_line = max_dist_from_line.max(perp.length());
            }

            // If all points are close to the line, it's a linear pattern
            let line_length = (centers[n - 1] - centers[0]).length();
            if line_length > TOLERANCE_MESH_LEGACY && max_dist_from_line / line_length < spacing_tolerance {
                return HolePatternType::Linear;
            }
        }
    }

    // Check for rectangular grid
    if n >= 4 {
        // Compute bounding box and check regular spacing
        let mut min_pt = centers[0];
        let mut max_pt = centers[0];
        for p in &centers[1..] {
            min_pt = min_pt.min(*p);
            max_pt = max_pt.max(*p);
        }

        let dims = max_pt - min_pt;
        let mut dim_count = 0;
        for &d in &[dims.x, dims.y, dims.z] {
            if d > TOLERANCE_MESH_LEGACY {
                dim_count += 1;
            }
        }

        if dim_count >= 2 {
            // Check if points form a regular grid
            let spacing_x = if dims.x > TOLERANCE_MESH_LEGACY {
                let unique_x: std::collections::BTreeSet<i64> = centers
                    .iter()
                    .map(|p| (p.x / dims.x * 100.0).round() as i64)
                    .collect();
                if unique_x.len() > 1 {
                    dims.x / (unique_x.len() - 1) as f64
                } else {
                    0.0
                }
            } else {
                0.0
            };

            let spacing_y = if dims.y > TOLERANCE_MESH_LEGACY {
                let unique_y: std::collections::BTreeSet<i64> = centers
                    .iter()
                    .map(|p| (p.y / dims.y * 100.0).round() as i64)
                    .collect();
                if unique_y.len() > 1 {
                    dims.y / (unique_y.len() - 1) as f64
                } else {
                    0.0
                }
            } else {
                0.0
            };

            if spacing_x > TOLERANCE_MESH_LEGACY || spacing_y > TOLERANCE_MESH_LEGACY {
                return HolePatternType::RectangularGrid;
            }
        }
    }

    HolePatternType::Irregular
}

/// Compute pattern origin, direction, and spacing from centers and pattern type.
fn compute_pattern_properties(
    centers: &[DVec3],
    pattern_type: &HolePatternType,
) -> (DVec3, DVec3, f64) {
    if centers.is_empty() {
        return (DVec3::ZERO, DVec3::Z, 0.0);
    }

    let origin = centers[0];

    match pattern_type {
        HolePatternType::Linear => {
            let direction = (centers[centers.len() - 1] - centers[0]).normalize_or_zero();
            let spacing = if centers.len() > 1 {
                (centers[centers.len() - 1] - centers[0]).length() / (centers.len() - 1) as f64
            } else {
                0.0
            };
            (origin, direction, spacing)
        }
        HolePatternType::Circular => {
            let centroid = centers.iter().fold(DVec3::ZERO, |acc, &p| acc + p) / centers.len() as f64;
            // Compute normal to the plane containing all centers
            let mut normal = DVec3::Z;
            if centers.len() >= 3 {
                let v1 = centers[1] - centers[0];
                let v2 = centers[2] - centers[0];
                normal = v1.cross(v2).normalize_or_zero();
                if normal.length_squared() < 0.5 {
                    normal = DVec3::Z;
                }
            }
            let spacing = if !centers.is_empty() {
                // Approximate angular spacing
                std::f64::consts::TAU / centers.len() as f64
            } else {
                0.0
            };
            (centroid, normal, spacing)
        }
        HolePatternType::RectangularGrid => {
            let mut min_pt = centers[0];
            let mut max_pt = centers[0];
            for p in &centers[1..] {
                min_pt = min_pt.min(*p);
                max_pt = max_pt.max(*p);
            }
            let center = (min_pt + max_pt) * 0.5;
            let dims = max_pt - min_pt;
            let spacing = dims.x.max(dims.y).max(dims.z)
                / (centers.len() as f64).sqrt().max(1.0);
            (center, DVec3::Z, spacing)
        }
        HolePatternType::Irregular => {
            let centroid = centers.iter().fold(DVec3::ZERO, |acc, &p| acc + p) / centers.len() as f64;
            (centroid, DVec3::Z, 0.0)
        }
    }
}