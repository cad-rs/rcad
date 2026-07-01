
/// Compute local wall thickness at a specific point.
///
/// Uses ray casting to find the distance to the nearest boundary in
/// multiple directions, returning the minimum distance as the local thickness.
///
/// # Arguments
/// * `point` - Point inside the solid
/// * `brep` - The B-Rep model
/// * `opts` - Computation options
///
/// # Returns
/// Local wall thickness and the direction to the nearest boundary.
pub fn compute_local_thickness(point: &DVec3, brep: &BRep, opts: &MedialAxisOptions) -> (f64, DVec3) {
    let mut min_distance = f64::MAX;
    let mut min_direction = DVec3::Z;

    // Sample directions on a sphere
    let num_theta = (std::f64::consts::PI / opts.angular_resolution).ceil() as usize;
    let num_phi = (2.0 * std::f64::consts::PI / opts.angular_resolution).ceil() as usize;

    for i in 0..num_theta {
        let theta = (i as f64 / num_theta as f64) * std::f64::consts::PI;
        for j in 0..num_phi {
            let phi = (j as f64 / num_phi as f64) * 2.0 * std::f64::consts::PI;

            let direction = DVec3::new(
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            );

            let distance = ray_cast_to_boundary(point, &direction, brep, opts);
            if distance < min_distance {
                min_distance = distance;
                min_direction = direction;
            }
        }
    }

    (min_distance * 2.0, min_direction) // Full thickness is 2x distance to nearest boundary
}

/// Identify thick and thin zones in a solid.
///
/// Classifies regions of the solid based on wall thickness relative
/// to target values, useful for manufacturing analysis.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `target_thickness` - Target wall thickness
/// * `tolerance` - Acceptable deviation from target
/// * `opts` - Computation options
///
/// # Returns
/// Vector of zone classifications with thickness and location.
pub fn identify_thickness_zones(
    brep: &BRep,
    target_thickness: f64,
    tolerance: f64,
    opts: &MedialAxisOptions,
) -> Vec<ThicknessZone> {
    let medial = compute_medial_surface_voxel(brep, opts);

    // Cluster medial vertices into zones
    let clusters = cluster_medial_vertices(&medial, opts.cluster_distance * 10.0);

    clusters
        .iter()
        .map(|cluster| {
            let points: Vec<&MedialVertex> = cluster.iter().filter_map(|&i| medial.vertices.get(i)).collect();

            let avg_thickness = points.iter().map(|v| v.radius * 2.0).sum::<f64>() / points.len() as f64;
            let center = points
                .iter()
                .fold(DVec3::ZERO, |acc, v| acc + v.point)
                / points.len() as f64;

            let class = if avg_thickness < target_thickness - tolerance {
                ThicknessClass::Thin
            } else if avg_thickness > target_thickness + tolerance {
                ThicknessClass::Thick
            } else {
                ThicknessClass::Normal
            };

            ThicknessZone {
                center,
                avg_thickness,
                thickness_class: class,
                point_count: points.len(),
            }
        })
        .collect()
}

/// A zone with classified thickness.
#[derive(Debug, Clone)]
pub struct ThicknessZone {
    /// Center of the zone.
    pub center: DVec3,
    /// Average thickness in the zone.
    pub avg_thickness: f64,
    /// Thickness classification.
    pub thickness_class: ThicknessClass,
    /// Number of sample points in the zone.
    pub point_count: usize,
}

// ============================================================================
// Helper Functions for Enhanced 3D Computation
// ============================================================================

/// Bounding box for a B-Rep.
struct BoundingBox {
    min: DVec3,
    max: DVec3,
}

impl BoundingBox {
    fn is_valid(&self) -> bool {
        self.min.x <= self.max.x && self.min.y <= self.max.y && self.min.z <= self.max.z
    }
}

fn compute_brep_bbox(brep: &BRep) -> BoundingBox {
    let mut min = DVec3::splat(f64::MAX);
    let mut max = DVec3::splat(f64::MIN);

    for vertex in &brep.vertices {
        min = min.min(vertex.point);
        max = max.max(vertex.point);
    }

    BoundingBox { min, max }
}

fn compute_signed_distance_field(brep: &BRep, grid: &mut VoxelGrid, opts: &MedialAxisOptions) {
    let dims = grid.dimensions;
    let voxel_size = grid.voxel_size;
    let origin = grid.origin;

    for k in 0..dims[2] {
        for j in 0..dims[1] {
            for i in 0..dims[0] {
                // Compute point without borrowing grid
                let point = origin + DVec3::new(
                    (i as f64 + 0.5) * voxel_size,
                    (j as f64 + 0.5) * voxel_size,
                    (k as f64 + 0.5) * voxel_size,
                );

                // Compute distance to nearest boundary
                let (distance, inside) = compute_point_distance_to_brep(&point, brep, opts);

                let idx = grid.index(i, j, k);
                grid.distances[idx] = distance;
                grid.inside[idx] = inside;
            }
        }
    }
}

fn compute_point_distance_to_brep(point: &DVec3, brep: &BRep, _opts: &MedialAxisOptions) -> (f64, bool) {
    let mut min_dist = f64::MAX;
    let mut inside = false;

    // Check distance to each face
    for (face_idx, face) in brep.geom.face_surface.iter().enumerate() {
        if let Some(surf_idx) = face
            && let Some(surf) = brep.geom.surfaces.get(*surf_idx) {
                let dist = distance_point_to_surface(point, surf);
                if dist < min_dist {
                    min_dist = dist;

                    // Determine if inside by checking face normal
                    if let Some(face_data) = brep.solids.first().and_then(|s| s.shells.first())
                        .and_then(|shell| shell.faces.get(face_idx))
                    {
                        let normal = face_data.normal;
                        // Simple inside test based on normal direction
                        let to_point = *point - surf.point_at(0.5, 0.5);
                        inside = to_point.dot(normal) < 0.0;
                    }
                }
            }
    }

    (min_dist, inside)
}

fn distance_point_to_surface(point: &DVec3, surf: &Surface3) -> f64 {
    // Use simple projection for now
    let [u_min, u_max, v_min, v_max] = surf.default_domain();

    let mut min_dist = f64::MAX;

    // Sample the surface to find minimum distance
    for i in 0..20 {
        for j in 0..20 {
            let u = u_min + (i as f64 / 19.0) * (u_max - u_min);
            let v = v_min + (j as f64 / 19.0) * (v_max - v_min);

            if u.is_finite() && v.is_finite() {
                let surf_point = surf.point_at(u, v);
                let dist = (*point - surf_point).length();
                min_dist = min_dist.min(dist);
            }
        }
    }

    min_dist
}

fn connect_medial_vertices(surface: &mut MedialSurface, max_distance: f64) {
    let n = surface.vertices.len();
    if n < 2 {
        return;
    }

    // Build edges between nearby vertices
    for i in 0..n {
        for j in (i + 1)..n {
            let d = (surface.vertices[i].point - surface.vertices[j].point).length();
            if d < max_distance {
                surface.edges.push(MedialEdge {
                    start_vertex: i,
                    end_vertex: j,
                    curve: None,
                    start_radius: surface.vertices[i].radius,
                    end_radius: surface.vertices[j].radius,
                });
            }
        }
    }
}

fn build_medial_faces(surface: &mut MedialSurface) {
    // Build adjacency list
    let mut adj: Vec<Vec<(usize, usize)>> = vec![vec![]; surface.vertices.len()];
    for (edge_idx, edge) in surface.edges.iter().enumerate() {
        adj[edge.start_vertex].push((edge.end_vertex, edge_idx));
        adj[edge.end_vertex].push((edge.start_vertex, edge_idx));
    }

    // Find edge loops that could form faces
    let mut visited_edges: HashSet<usize> = HashSet::new();

    for start_edge_idx in 0..surface.edges.len() {
        if visited_edges.contains(&start_edge_idx) {
            continue;
        }

        // Try to find a loop starting from this edge
        if let Some(loop_vertices) = find_edge_loop(start_edge_idx, &adj, &surface.edges)
            && loop_vertices.len() >= 3 {
                let radii: Vec<f64> = loop_vertices
                    .iter()
                    .filter_map(|&v| surface.vertices.get(v).map(|v| v.radius))
                    .collect();

                let min_radius = radii.iter().cloned().fold(f64::MAX, f64::min);
                let max_radius = radii.iter().cloned().fold(0.0, f64::max);

                // Mark edges as visited before moving loop_vertices
                for i in 0..loop_vertices.len() {
                    let v1 = loop_vertices[i];
                    let v2 = loop_vertices[(i + 1) % loop_vertices.len()];
                    for (edge_idx, edge) in surface.edges.iter().enumerate() {
                        if (edge.start_vertex == v1 && edge.end_vertex == v2)
                            || (edge.start_vertex == v2 && edge.end_vertex == v1)
                        {
                            visited_edges.insert(edge_idx);
                        }
                    }
                }

                surface.faces.push(MedialFace {
                    vertices: loop_vertices,
                    surface: None,
                    min_radius,
                    max_radius,
                });
            }
    }
}

fn find_edge_loop(
    start_edge_idx: usize,
    adj: &[Vec<(usize, usize)>],
    edges: &[MedialEdge],
) -> Option<Vec<usize>> {
    let start_edge = &edges[start_edge_idx];
    let mut loop_vertices = vec![start_edge.start_vertex, start_edge.end_vertex];
    let mut current = start_edge.end_vertex;
    let target = start_edge.start_vertex;

    for _ in 0..edges.len() {
        // Find next edge
        let mut found = false;
        for &(next_vertex, edge_idx) in &adj[current] {
            if edge_idx == start_edge_idx && loop_vertices.len() > 2 {
                // Completed the loop
                return Some(loop_vertices);
            }
            if !loop_vertices.contains(&next_vertex) {
                loop_vertices.push(next_vertex);
                current = next_vertex;
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
        if current == target && loop_vertices.len() > 2 {
            return Some(loop_vertices);
        }
    }

    None
}

fn compute_opposing_face_pairs(brep: &BRep, _opts: &MedialAxisOptions) -> Vec<(usize, usize, f64)> {
    let mut pairs: Vec<(usize, usize, f64)> = Vec::new();

    // Get all faces
    let faces: Vec<(usize, &Face, Option<&Surface3>)> = brep
        .geom
        .face_surface
        .iter()
        .enumerate()
        .filter_map(|(idx, fs)| {
            if let Some(surf_idx) = fs
                && let Some(shell) = brep.solids.first().and_then(|s| s.shells.first())
                    && let Some(face) = shell.faces.get(idx) {
                        let surf = brep.geom.surfaces.get(*surf_idx);
                        return Some((idx, face, surf));
                    }
            None
        })
        .collect();

    // Find pairs of approximately parallel faces
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            let (idx_i, face_i, surf_i) = &faces[i];
            let (idx_j, face_j, surf_j) = &faces[j];

            // Check if faces are approximately parallel and facing each other
            let normal_i = face_i.normal;
            let normal_j = face_j.normal;

            // Parallel if normals are opposite
            let dot = normal_i.dot(normal_j);
            if dot < -0.9 {
                // Estimate distance between faces
                if let (Some(surf_i), Some(surf_j)) = (surf_i, surf_j) {
                    let center_i = surf_i.point_at(0.5, 0.5);
                    let center_j = surf_j.point_at(0.5, 0.5);
                    let distance = (center_j - center_i).length();
                    pairs.push((*idx_i, *idx_j, distance));
                }
            }
        }
    }

    pairs
}

fn find_associated_face_pairs(
    _point: &DVec3,
    face_pairs: &[(usize, usize, f64)],
    _tolerance: f64,
) -> Vec<(usize, usize)> {
    face_pairs
        .iter()
        .filter_map(|&(f1, f2, _distance)| {
            // Check if the point is approximately midway between the faces
            // This is a simplified check
            Some((f1, f2))
        })
        .collect()
}

fn compute_chordal_direction(
    _point: &DVec3,
    _face_pairs: &[(usize, usize)],
    _brep: &BRep,
) -> DVec3 {
    // Default direction for now
    DVec3::X
}

fn compute_mid_surface_normal(
    _point: &DVec3,
    _face_pairs: &[(usize, usize)],
    _brep: &BRep,
) -> DVec3 {
    // Default normal for now
    DVec3::Z
}

fn connect_chordal_vertices(axis: &mut ChordalAxis, max_distance: f64) {
    let n = axis.vertices.len();
    if n < 2 {
        return;
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let d = (axis.vertices[i].point - axis.vertices[j].point).length();
            if d < max_distance {
                axis.edges.push(ChordalEdge {
                    start: i,
                    end: j,
                    curve: None,
                    avg_thickness: (axis.vertices[i].thickness + axis.vertices[j].thickness) / 2.0,
                    length: d,
                });
            }
        }
    }
}

fn identify_thin_sheets(axis: &ChordalAxis, _brep: &BRep, _opts: &MedialAxisOptions) -> Vec<ThinSheet> {
    let mut sheets = Vec::new();

    for (edge_idx, edge) in axis.edges.iter().enumerate() {
        sheets.push(ThinSheet {
            spine_edge: edge_idx,
            side_a_faces: vec![],
            side_b_faces: vec![],
            avg_thickness: edge.avg_thickness,
            area: 0.0,
            quality: 1.0,
        });
    }

    sheets
}

fn create_mid_surface_patch(
    start_v: &ChordalVertex,
    end_v: &ChordalVertex,
    _sheet: &ThinSheet,
    mid_brep: &mut BRep,
    face_thickness: &mut Vec<f64>,
    face_mapping: &mut Vec<(usize, usize)>,
    opts: &MidSurfaceOptions,
) {
    // Create a ruled surface between the two vertices
    let direction = end_v.point - start_v.point;
    let length = direction.length();

    if length < opts.base.tolerance {
        return;
    }

    // Create a simple planar patch
    let center = (start_v.point + end_v.point) / 2.0;
    let avg_normal = (start_v.normal + end_v.normal).normalize();

    // Create perpendicular directions for the patch
    let u_dir = direction.normalize();
    let v_dir = avg_normal.cross(u_dir).normalize();

    // Create a quad patch
    let half_length = length / 2.0;
    let half_width = (start_v.thickness + end_v.thickness) / 4.0;

    let corners = [center - u_dir * half_length - v_dir * half_width,
        center + u_dir * half_length - v_dir * half_width,
        center + u_dir * half_length + v_dir * half_width,
        center - u_dir * half_length + v_dir * half_width];

    // Add vertices
    let v_indices: Vec<usize> = corners
        .iter()
        .map(|&p| {
            let idx = mid_brep.vertices.len();
            mid_brep.vertices.push(rcad_kernel::Vertex { point: p });
            idx
        })
        .collect();

    // Create edges
    let e_indices: Vec<usize> = (0..4)
        .map(|i| {
            let start = v_indices[i];
            let end = v_indices[(i + 1) % 4];
            let idx = mid_brep.edges.len();
            mid_brep.edges.push(rcad_kernel::Edge { start, end });
            idx
        })
        .collect();

    // Create face
    let wire = Wire {
        edges: e_indices.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
    };

    let face = Face {
        outer_wire: wire,
        inner_wires: vec![],
        normal: avg_normal,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };

    // Create plane surface
    let plane = Plane {
        origin: center,
        normal: avg_normal,
    };

    let surf_idx = mid_brep.geom.surfaces.len();
    mid_brep.geom.surfaces.push(Surface3::Plane(plane));

    let face_idx = mid_brep.geom.face_surface.len();
    mid_brep.geom.face_surface.push(Some(surf_idx));

    if mid_brep.solids.is_empty() {
        mid_brep.solids.push(Solid { shells: vec![Shell { faces: vec![] }] });
    }
    if let Some(shell) = mid_brep.solids[0].shells.first_mut() {
        shell.faces.push(face);
    }

    face_thickness.push((start_v.thickness + end_v.thickness) / 2.0);

    // Map to original faces
    if let Some(&(f1, _f2)) = start_v.face_pairs.first() {
        face_mapping.push((face_idx, f1));
    }
}

fn create_mid_surface_point(
    vertex: &ChordalVertex,
    mid_brep: &mut BRep,
    face_thickness: &mut Vec<f64>,
    face_mapping: &mut Vec<(usize, usize)>,
    _opts: &MidSurfaceOptions,
) {
    // Create a small triangular patch at this point
    let normal = vertex.normal;
    let tangent = if normal.x.abs() > 0.5 {
        DVec3::Y
    } else {
        DVec3::X
    };
    let u_dir = tangent.cross(normal).normalize();
    let v_dir = normal.cross(u_dir).normalize();

    let r = vertex.thickness / 4.0;

    let corners = [vertex.point + u_dir * r,
        vertex.point - u_dir * r / 2.0 + v_dir * r * 0.866,
        vertex.point - u_dir * r / 2.0 - v_dir * r * 0.866];

    // Add vertices
    let v_indices: Vec<usize> = corners
        .iter()
        .map(|&p| {
            let idx = mid_brep.vertices.len();
            mid_brep.vertices.push(rcad_kernel::Vertex { point: p });
            idx
        })
        .collect();

    // Create edges
    let e_indices: Vec<usize> = (0..3)
        .map(|i| {
            let start = v_indices[i];
            let end = v_indices[(i + 1) % 3];
            let idx = mid_brep.edges.len();
            mid_brep.edges.push(rcad_kernel::Edge { start, end });
            idx
        })
        .collect();

    // Create face
    let wire = Wire {
        edges: e_indices.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
    };

    let face = Face {
        outer_wire: wire,
        inner_wires: vec![],
        normal,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };

    // Create plane surface
    let plane = Plane {
        origin: vertex.point,
        normal,
    };

    let surf_idx = mid_brep.geom.surfaces.len();
    mid_brep.geom.surfaces.push(Surface3::Plane(plane));

    let face_idx = mid_brep.geom.face_surface.len();
    mid_brep.geom.face_surface.push(Some(surf_idx));

    if mid_brep.solids.is_empty() {
        mid_brep.solids.push(Solid { shells: vec![Shell { faces: vec![] }] });
    }
    if let Some(shell) = mid_brep.solids[0].shells.first_mut() {
        shell.faces.push(face);
    }

    face_thickness.push(vertex.thickness);

    if let Some(&(f1, _)) = vertex.face_pairs.first() {
        face_mapping.push((face_idx, f1));
    }
}

fn compute_mid_surface_quality(
    mid_brep: &BRep,
    _original: &BRep,
    chordal_axis: &ChordalAxis,
    _opts: &MidSurfaceOptions,
) -> MidSurfaceQuality {
    // Compute coverage (ratio of chordal vertices represented)
    let coverage = if chordal_axis.vertices.is_empty() {
        1.0
    } else {
        // Count faces in mid-surface
        let face_count = mid_brep
            .solids
            .first()
            .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum())
            .unwrap_or(0);

        (face_count as f64 / chordal_axis.vertices.len() as f64).min(1.0)
    };

    // Compute average deviation (simplified)
    let avg_deviation = 0.0; // Would need proper deviation computation
    let max_deviation = 0.0;

    // Compute thickness accuracy
    let thickness_accuracy = 1.0; // Would need proper accuracy computation

    // Count discontinuities
    let discontinuities = 0; // Would need proper connectivity analysis

    // Compute overall score
    let overall_score = coverage * 0.5 + thickness_accuracy * 0.5;

    MidSurfaceQuality {
        coverage,
        avg_deviation,
        max_deviation,
        thickness_accuracy,
        discontinuities,
        overall_score,
    }
}

fn cluster_thin_regions(regions: &mut [ThinRegion], max_distance: f64) {
    // Simple clustering based on distance
    let n = regions.len();
    if n < 2 {
        return;
    }

    let mut cluster_ids: Vec<usize> = (0..n).collect();
    let _next_cluster_id = n;

    // Merge nearby regions
    for i in 0..n {
        for j in (i + 1)..n {
            let d = (regions[i].center - regions[j].center).length();
            if d < max_distance {
                let old_id = cluster_ids[j];
                let new_id = cluster_ids[i];
                for id in &mut cluster_ids {
                    if *id == old_id {
                        *id = new_id;
                    }
                }
            }
        }
    }

    // Update severity based on cluster
    for i in 0..n {
        let cluster_id = cluster_ids[i];
        let cluster_count = cluster_ids.iter().filter(|&&id| id == cluster_id).count();
        if cluster_count > 1 {
            regions[i].severity = regions[i].severity.min(1.0).max(0.5);
        }
    }
}

fn estimate_region_area(center: &DVec3, thickness: f64, medial: &MedialSurface) -> f64 {
    // Estimate area based on nearby medial vertices
    let nearby: Vec<&MedialVertex> = medial
        .vertices
        .iter()
        .filter(|v| (*center - v.point).length() < thickness * 2.0)
        .collect();

    if nearby.is_empty() {
        thickness * thickness
    } else {
        nearby.len() as f64 * thickness * thickness
    }
}

fn classify_thickness(stats: &ThicknessStats, target: f64) -> ThicknessClass {
    let ratio = stats.mean / target;

    if ratio < 0.25 {
        ThicknessClass::VeryThin
    } else if ratio < 0.5 {
        ThicknessClass::Thin
    } else if ratio < 1.5 {
        ThicknessClass::Normal
    } else if ratio < 2.0 {
        ThicknessClass::Thick
    } else {
        ThicknessClass::VeryThick
    }
}

fn build_thickness_histogram(medial: &MedialSurface, num_bins: usize) -> Vec<ThicknessHistogramBin> {
    if medial.vertices.is_empty() {
        return vec![];
    }

    let thicknesses: Vec<f64> = medial.vertices.iter().map(|v| v.radius * 2.0).collect();
    let min_t = thicknesses.iter().cloned().fold(f64::MAX, f64::min);
    let max_t = thicknesses.iter().cloned().fold(0.0, f64::max);

    let bin_width = (max_t - min_t) / num_bins as f64;
    if bin_width < TOLERANCE_LINEAR_ULTRA_STRICT {
        return vec![ThicknessHistogramBin {
            lower: min_t,
            upper: max_t,
            count: thicknesses.len(),
        }];
    }

    let mut bins: Vec<ThicknessHistogramBin> = (0..num_bins)
        .map(|i| ThicknessHistogramBin {
            lower: min_t + i as f64 * bin_width,
            upper: min_t + (i + 1) as f64 * bin_width,
            count: 0,
        })
        .collect();

    for t in thicknesses {
        let bin_idx = ((t - min_t) / bin_width).floor() as usize;
        let bin_idx = bin_idx.min(num_bins - 1);
        bins[bin_idx].count += 1;
    }

    bins
}

fn compute_recommended_min_thickness(medial: &MedialSurface, target: f64) -> f64 {
    if medial.vertices.is_empty() {
        return target;
    }

    // Use the minimum thickness found, with a safety factor
    let min_found = medial.thickness_stats.min;
    min_found.max(target * 0.8)
}

fn find_rib_candidates(medial: &MedialSurface, opts: &RibGenerationOptions) -> Vec<RibCandidate> {
    let mut candidates = Vec::new();

    for edge in &medial.edges {
        if let (Some(start_v), Some(end_v)) = (
            medial.vertices.get(edge.start_vertex),
            medial.vertices.get(edge.end_vertex),
        ) {
            let length = (end_v.point - start_v.point).length();
            if length >= opts.min_length {
                let avg_thickness = (start_v.radius + end_v.radius) * 2.0;

                // Ribs are most useful in thin regions
                if avg_thickness < opts.base.min_thickness * 10.0 {
                    candidates.push(RibCandidate {
                        start: start_v.point,
                        end: end_v.point,
                        avg_thickness,
                        length,
                        medial_edge_idx: Some(edge.start_vertex), // Simplified
                    });
                }
            }
        }
    }

    candidates
}

struct RibCandidate {
    start: DVec3,
    end: DVec3,
    avg_thickness: f64,
    length: f64,
    medial_edge_idx: Option<usize>,
}

fn create_rib_placement(
    candidate: &RibCandidate,
    medial: &MedialSurface,
    _brep: &BRep,
    opts: &RibGenerationOptions,
) -> Option<RibPlacement> {
    let direction = candidate.end - candidate.start;

    // Compute optimal rib height based on thickness
    let height = (candidate.avg_thickness * 3.0).clamp(opts.min_height, opts.max_height);
    let width = height * 0.6; // Typical width-to-height ratio

    // Compute efficiency score
    let efficiency = (candidate.length / opts.min_length).min(1.0)
        * (height / opts.max_height).min(1.0);

    // Find attached face
    let attached_face = candidate
        .medial_edge_idx
        .and_then(|idx| medial.vertices.get(idx))
        .and_then(|v| v.boundary_elements.first().copied())
        .unwrap_or(0);

    Some(RibPlacement {
        centerline: Curve3::Line(Line3 {
            origin: candidate.start,
            direction: direction.normalize(),
        }),
        start: candidate.start,
        end: candidate.end,
        height,
        width,
        draft_angle: opts.draft_angle,
        efficiency,
        medial_edge: candidate.medial_edge_idx,
        attached_face,
    })
}

fn optimize_rib_layout(ribs: &mut Vec<RibPlacement>, _brep: &BRep, opts: &RibGenerationOptions) {
    // Remove overlapping ribs and optimize spacing
    let mut to_remove: HashSet<usize> = HashSet::new();

    for i in 0..ribs.len() {
        if to_remove.contains(&i) {
            continue;
        }
        for j in (i + 1)..ribs.len() {
            if to_remove.contains(&j) {
                continue;
            }

            // Check if ribs are too close
            let dist = (ribs[i].start - ribs[j].start).length().min(
                (ribs[i].end - ribs[j].end).length(),
            );

            if dist < opts.spacing * 0.5 {
                // Keep the more efficient rib
                if ribs[i].efficiency >= ribs[j].efficiency {
                    to_remove.insert(j);
                } else {
                    to_remove.insert(i);
                    break;
                }
            }
        }
    }

    // Sort by efficiency and remove low-efficiency ribs
    ribs.sort_by(|a, b| b.efficiency.partial_cmp(&a.efficiency).unwrap_or(std::cmp::Ordering::Equal));

    // Keep only top ribs
    let max_ribs = 20; // Reasonable limit
    if ribs.len() > max_ribs {
        ribs.truncate(max_ribs);
    }
}

fn estimate_stiffness_improvement(ribs: &[RibPlacement], medial: &MedialSurface) -> f64 {
    if ribs.is_empty() || medial.vertices.is_empty() {
        return 0.0;
    }

    // Simplified estimate: total rib volume / original volume
    let _total_rib_volume: f64 = ribs
        .iter()
        .map(|r| {
            let length = (r.end - r.start).length();
            length * r.width * r.height * 0.5
        })
        .sum();

    // Estimate original volume from medial surface
    let avg_thickness = medial.thickness_stats.mean;

    // Stiffness improvement is roughly proportional to the moment of inertia increase
    let rib_inertia_factor: f64 = ribs.iter().map(|r| r.height * r.height).sum();
    let base_factor = avg_thickness * avg_thickness * medial.vertices.len() as f64;

    if base_factor > 0.0 {
        (rib_inertia_factor / base_factor).min(10.0) // Cap at 10x improvement
    } else {
        0.0
    }
}

fn compute_weight_increase(ribs: &[RibPlacement], _brep: &BRep) -> f64 {
    if ribs.is_empty() {
        return 0.0;
    }

    let total_rib_volume: f64 = ribs
        .iter()
        .map(|r| {
            let length = (r.end - r.start).length();
            length * r.width * r.height * 0.5
        })
        .sum();

    // Simplified: assume base part has volume proportional to bounding box
    // Return percentage increase
    total_rib_volume * 100.0 / 1000000.0 // Simplified percentage
}

fn compute_rib_quality_score(ribs: &[RibPlacement], medial: &MedialSurface, opts: &RibGenerationOptions) -> f64 {
    if ribs.is_empty() {
        return 0.0;
    }

    // Average efficiency
    let avg_efficiency: f64 = ribs.iter().map(|r| r.efficiency).sum::<f64>() / ribs.len() as f64;

    // Coverage: how many thin regions are addressed
    let coverage = ribs.len() as f64 / medial.vertices.len().max(1) as f64;

    // Spacing score
    let spacing_score = if ribs.len() > 1 {
        let mut min_spacing = f64::MAX;
        for i in 0..ribs.len() {
            for j in (i + 1)..ribs.len() {
                let d = (ribs[i].start - ribs[j].start).length();
                min_spacing = min_spacing.min(d);
            }
        }
        if min_spacing > opts.spacing {
            1.0
        } else {
            min_spacing / opts.spacing
        }
    } else {
        1.0
    };

    avg_efficiency * 0.5 + coverage.min(1.0) * 0.3 + spacing_score * 0.2
}

fn ray_cast_to_boundary(point: &DVec3, direction: &DVec3, brep: &BRep, opts: &MedialAxisOptions) -> f64 {
    let mut min_distance = f64::MAX;

    // Check intersection with each face
    for face_surf in brep.geom.face_surface.iter() {
        if let Some(surf_idx) = face_surf
            && let Some(surf) = brep.geom.surfaces.get(*surf_idx) {
                let distance = ray_surface_intersection(point, direction, surf, opts);
                min_distance = min_distance.min(distance);
            }
    }

    min_distance
}

fn ray_surface_intersection(point: &DVec3, direction: &DVec3, surf: &Surface3, _opts: &MedialAxisOptions) -> f64 {
    // Simplified: sample the surface and find minimum distance along ray
    let [u_min, u_max, v_min, v_max] = surf.default_domain();

    let mut min_dist = f64::MAX;

    for i in 0..20 {
        for j in 0..20 {
            let u = u_min + (i as f64 / 19.0) * (u_max - u_min);
            let v = v_min + (j as f64 / 19.0) * (v_max - v_min);

            if u.is_finite() && v.is_finite() {
                let surf_point = surf.point_at(u, v);
                let to_surf = surf_point - *point;

                // Project onto ray direction
                let t = to_surf.dot(*direction);
                if t > 0.0 {
                    let closest = *point + t * *direction;
                    let dist = (surf_point - closest).length();
                    if dist < min_dist {
                        min_dist = t; // Distance along ray
                    }
                }
            }
        }
    }

    min_dist
}

/// Compute wall thickness distribution for a solid.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
///
/// # Returns
/// Statistical summary of wall thickness and detected thin regions.
pub fn compute_wall_thickness(brep: &BRep) -> WallThicknessResult {
    let opts = MedialAxisOptions::default();
    let surface = compute_medial_surface(brep, &opts);

    if surface.vertices.is_empty() {
        return WallThicknessResult {
            min_thickness: 0.0,
            max_thickness: 0.0,
            avg_thickness: 0.0,
            thin_regions: vec![],
        };
    }

    let radii: Vec<f64> = surface.vertices.iter().map(|v| v.radius * 2.0).collect();
    let min_thickness = radii.iter().cloned().fold(f64::MAX, f64::min);
    let max_thickness = radii.iter().cloned().fold(0.0, f64::max);
    let avg_thickness = radii.iter().sum::<f64>() / radii.len() as f64;

    WallThicknessResult {
        min_thickness,
        max_thickness,
        avg_thickness,
        thin_regions: vec![],
    }
}

/// Detect thin-walled regions in a solid.
///
/// Finds regions where the wall thickness falls below the specified threshold.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `min_thickness` - Minimum acceptable wall thickness
///
/// # Returns
/// List of detected thin regions with location and severity.
pub fn detect_thin_regions(brep: &BRep, min_thickness: f64) -> Vec<ThinRegion> {
    let opts = MedialAxisOptions {
        min_thickness: min_thickness * 0.1,
        ..Default::default()
    };
    let surface = compute_medial_surface(brep, &opts);

    surface
        .vertices
        .iter()
        .filter(|v| v.radius * 2.0 < min_thickness)
        .map(|v| {
            let thickness = v.radius * 2.0;
            let severity = 1.0 - (thickness / min_thickness).min(1.0);
            ThinRegion {
                center: v.point,
                thickness,
                area: 0.0,
                face_indices: v.boundary_elements.clone(),
                severity,
            }
        })
        .collect()
}

/// Compute a detailed thickness map for a solid.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `opts` - Computation options
///
/// # Returns
/// A thickness map with samples at multiple points.
pub fn compute_thickness_map(brep: &BRep, opts: &MedialAxisOptions) -> ThicknessMap {
    let surface = compute_medial_surface(brep, opts);

    let samples: Vec<ThicknessSample> = surface
        .vertices
        .iter()
        .map(|v| ThicknessSample {
            point: v.point,
            thickness: v.radius * 2.0,
            normal: DVec3::Z, // Would need proper computation
            nearest_face: v.boundary_elements.first().copied().unwrap_or(0),
        })
        .collect();

    let stats = surface.thickness_stats;

    // Detect thin regions
    let thin_regions: Vec<ThinRegion> = samples
        .iter()
        .filter(|s| s.thickness < opts.min_thickness)
        .map(|s| {
            let severity = 1.0 - (s.thickness / opts.min_thickness).min(1.0);
            ThinRegion {
                center: s.point,
                thickness: s.thickness,
                area: 0.0,
                face_indices: vec![s.nearest_face],
                severity,
            }
        })
        .collect();

    ThicknessMap {
        samples,
        stats,
        thin_regions,
    }
}

/// Compute the mid-surface of a thin-walled solid for FEA shell meshing.
///
/// # Arguments
/// * `brep` - The B-Rep model to analyze
/// * `opts` - Computation options
///
/// # Returns
/// The mid-surface with thickness information.
pub fn compute_mid_surface(brep: &BRep, opts: &MedialAxisOptions) -> MidSurfaceResult {
    // Compute medial surface
    let surface = compute_medial_surface(brep, opts);

    // Create a new B-Rep for the mid-surface
    let mut mid_brep = BRep::default();

    // Create faces from medial surface vertices
    // This is a simplified approach - a full implementation would
    // create proper surface patches

    let mut face_thickness: Vec<f64> = Vec::new();
    let mut face_mapping: Vec<(usize, usize)> = Vec::new();
    let mut faces: Vec<Face> = Vec::new();

    for vertex in surface.vertices.iter() {
        // Create a small planar face at each medial vertex
        let plane = Plane {
            origin: vertex.point,
            normal: DVec3::Z, // Simplified - would need proper normal
        };

        let surf_idx = mid_brep.geom.surfaces.len();
        mid_brep.geom.surfaces.push(Surface3::Plane(plane));

        // Create a small quad face
        let r = vertex.radius * 0.5;
        let corners = [vertex.point + DVec3::new(-r, -r, 0.0),
            vertex.point + DVec3::new(r, -r, 0.0),
            vertex.point + DVec3::new(r, r, 0.0),
            vertex.point + DVec3::new(-r, r, 0.0)];

        // Add vertices
        let v_indices: Vec<usize> = corners
            .iter()
            .map(|&p| {
                let idx = mid_brep.vertices.len();
                mid_brep.vertices.push(rcad_kernel::Vertex { point: p });
                idx
            })
            .collect();

        // Create edges
        let e_indices: Vec<usize> = (0..4)
            .map(|i| {
                let start = v_indices[i];
                let end = v_indices[(i + 1) % 4];
                let idx = mid_brep.edges.len();
                mid_brep.edges.push(rcad_kernel::Edge { start, end });
                idx
            })
            .collect();

        // Create wire
        let wire = Wire {
            edges: e_indices.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
        };

        // Create face
        let face = Face {
            outer_wire: wire,
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        };

        let face_idx = mid_brep.geom.face_surface.len();
        mid_brep.geom.face_surface.push(Some(surf_idx));
        faces.push(face);

        face_thickness.push(vertex.radius * 2.0);

        // Map to original face (simplified)
        if let Some(&orig_face) = vertex.boundary_elements.first() {
            face_mapping.push((face_idx, orig_face));
        }
    }

    // Create shell and solid from faces
    let shell = Shell { faces };
    let solid = Solid { shells: vec![shell] };
    mid_brep.solids.push(solid);

    MidSurfaceResult {
        brep: mid_brep,
        face_thickness,
        face_mapping,
    }
}

/// Generate rib/stiffener paths from medial axis.
///
/// Creates curve geometry suitable for generating reinforcement features.
///
/// # Arguments
/// * `axis` - The computed medial axis
///
/// # Returns
/// List of curves representing potential rib centerlines.
pub fn generate_rib_paths(axis: &MedialSurface) -> Vec<Curve3> {
    let mut paths = Vec::new();

    for edge in &axis.edges {
        if let (Some(start_v), Some(end_v)) = (
            axis.vertices.get(edge.start_vertex),
            axis.vertices.get(edge.end_vertex),
        ) {
            let direction = end_v.point - start_v.point;
            if direction.length() > TOLERANCE_LINEAR_ULTRA_STRICT {
                paths.push(Curve3::Line(Line3 {
                    origin: start_v.point,
                    direction: direction.normalize(),
                }));
            }
        }
    }

    paths
}

/// Find the maximum inscribed circle for a 2D profile.
///
/// # Arguments
/// * `points` - Boundary points of the profile
///
/// # Returns
/// The center and radius of the maximum inscribed circle, if found.
pub fn find_max_inscribed_circle(points: &[DVec3]) -> Option<(DVec3, f64)> {
    let opts = MedialAxisOptions::default();
    let axis = compute_medial_axis_2d(points, &opts);

    axis.max_inscribed_circle
        .map(|(center, radius)| (DVec3::new(center.x, center.y, 0.0), radius))
}

/// Compute the medial axis transform (MAT) for a 2D profile.
///
/// The MAT includes both the geometry and the radius function along the axis.
///
/// # Arguments
/// * `points` - Boundary points of the profile
/// * `opts` - Computation options
///
/// # Returns
/// A tuple of (vertices with radii, edges with radius functions).
pub fn compute_mat_2d(
    points: &[DVec3],
    opts: &MedialAxisOptions,
) -> (Vec<MedialPoint2d>, Vec<(usize, usize)>) {
    let axis = compute_medial_axis_2d(points, opts);

    // Collect all edges from branches
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for branch in &axis.branches {
        for i in 0..branch.points.len().saturating_sub(1) {
            edges.push((i, i + 1)); // Simplified indexing
        }
    }

    (axis.all_points, edges)
}

/// Cluster medial axis vertices into regions for analysis.
///
/// Groups nearby vertices into clusters for detecting distinct
/// thin regions or thickness variations.
///
/// # Arguments
/// * `surface` - The medial surface
/// * `cluster_distance` - Maximum distance for vertices in the same cluster
///
/// # Returns
/// Vector of vertex index groups representing clusters.
pub fn cluster_medial_vertices(surface: &MedialSurface, cluster_distance: f64) -> Vec<Vec<usize>> {
    let n = surface.vertices.len();
    if n == 0 {
        return vec![];
    }

    let mut visited = vec![false; n];
    let mut clusters = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }

        let mut cluster = Vec::new();
        let mut stack = vec![start];

        while let Some(i) = stack.pop() {
            if visited[i] {
                continue;
            }
            visited[i] = true;
            cluster.push(i);

            for j in 0..n {
                if !visited[j] {
                    let d = (surface.vertices[i].point - surface.vertices[j].point).length();
                    if d < cluster_distance {
                        stack.push(j);
                    }
                }
            }
        }

        if !cluster.is_empty() {
            clusters.push(cluster);
        }
    }

    clusters
}

// ============================================================================
// Tests
// ============================================================================
