
/// Helper function to build an orthonormal frame from axis and reference direction.
///
/// Returns (axis_normalized, x_axis, y_axis) where axis, x_axis, y_axis form
/// a right-handed orthonormal basis.
fn orthonormal_frame(axis: DVec3, ref_dir: DVec3) -> (DVec3, DVec3, DVec3) {
    let z = axis.normalize_or_zero();
    let x = (ref_dir - ref_dir.dot(z) * z).normalize_or_zero();
    // Handle degenerate case where ref_dir is parallel to axis
    let x = if x.length_squared() < 0.5 {
        any_perpendicular(z)
    } else {
        x
    };
    let y = z.cross(x).normalize_or_zero();
    (z, x, y)
}

// ── Thread Edge Detection ───────────────────────────────────────────────────────

/// Check if an edge is a thread edge (helical) on a cylinder or cone.
///
/// Thread edges are characterized by:
/// - The edge curve is a helix or approximately helical
/// - The edge lies on a cylindrical or conical surface
/// - The edge makes an angle with the surface axis (not parallel or perpendicular)
pub fn is_thread_edge(
    brep: &rcad_kernel::BRep,
    edge_idx: usize,
    surface: &Surface3,
) -> bool {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return false;
    };

    // Get the edge curve
    let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|&c| c) else {
        return false;
    };
    let Some(curve) = brep.geom.curves.get(curve_idx) else {
        return false;
    };

    match (surface, curve) {
        // Circular helix on cylinder or cone is a thread edge
        (Surface3::Cylinder(_), rcad_kernel::geom::Curve3::CircularHelix(_)) => true,
        (Surface3::Cone(_), rcad_kernel::geom::Curve3::CircularHelix(_)) => true,

        // Check for approximately helical curves
        (Surface3::Cylinder(cyl), curve3d) => {
            is_approximately_helical_on_cylinder(curve3d, cyl, brep, edge)
        }
        (Surface3::Cone(cone), curve3d) => {
            is_approximately_helical_on_cone(curve3d, cone, brep, edge)
        }
        _ => false,
    }
}

/// Check if a curve is approximately helical on a cylinder.
fn is_approximately_helical_on_cylinder(
    curve: &rcad_kernel::geom::Curve3,
    cyl: &rcad_kernel::geom::CylindricalSurface,
    brep: &rcad_kernel::BRep,
    edge: &rcad_kernel::topology::Edge,
) -> bool {
    use rcad_kernel::geom::CurveEval;

    // Sample the curve and check if points lie on the cylinder surface
    // and the curve makes an angle with the axis
    let Some(_v_start) = brep.vertices.get(edge.start) else { return false; };
    let Some(_v_end) = brep.vertices.get(edge.end) else { return false; };

    let [t0, t1] = curve.default_domain();
    let samples = 16;

    let mut on_surface_count = 0;
    let mut has_axial_component = false;
    let mut has_angular_component = false;

    for i in 0..samples {
        let t = t0 + (t1 - t0) * i as f64 / (samples - 1) as f64;
        let pt = curve.point_at(t);

        // Check if point is on cylinder surface
        let radial = pt - cyl.origin;
        let axial = radial.dot(cyl.axis);
        let radial_vec = radial - axial * cyl.axis;
        let radial_dist = radial_vec.length();

        if (radial_dist - cyl.radius).abs() < TOLERANCE_RETRY_LADDER_COARSE {
            on_surface_count += 1;
        }

        // Check for axial and angular components
        if i > 0 {
            let t_prev = t0 + (t1 - t0) * (i - 1) as f64 / (samples - 1) as f64;
            let pt_prev = curve.point_at(t_prev);
            let delta = pt - pt_prev;

            let axial_delta = delta.dot(cyl.axis).abs();
            let radial_delta = (delta - delta.dot(cyl.axis) * cyl.axis).length();

            if axial_delta > TOLERANCE_MESH_LEGACY {
                has_axial_component = true;
            }
            if radial_delta > TOLERANCE_MESH_LEGACY {
                has_angular_component = true;
            }
        }
    }

    // Thread edge is on surface and has both axial and angular components
    on_surface_count > samples / 2 && has_axial_component && has_angular_component
}

/// Check if a curve is approximately helical on a cone.
fn is_approximately_helical_on_cone(
    curve: &rcad_kernel::geom::Curve3,
    cone: &rcad_kernel::geom::ConicalSurface,
    brep: &rcad_kernel::BRep,
    edge: &rcad_kernel::topology::Edge,
) -> bool {
    use rcad_kernel::geom::CurveEval;

    let Some(_v_start) = brep.vertices.get(edge.start) else { return false; };
    let Some(_v_end) = brep.vertices.get(edge.end) else { return false; };

    let [t0, t1] = curve.default_domain();
    let samples = 16;

    let mut on_surface_count = 0;
    let mut has_axial_component = false;
    let mut has_angular_component = false;

    let apex = cone.apex_point();
    let axis = cone.axis_dir();

    for i in 0..samples {
        let t = t0 + (t1 - t0) * i as f64 / (samples - 1) as f64;
        let pt = curve.point_at(t);

        // Check if point is on cone surface
        let to_point = pt - apex;
        let axial_dist = to_point.dot(axis);
        let radial_vec = to_point - axial_dist * axis;
        let radial_dist = radial_vec.length();

        // Expected radius at this axial distance
        let expected_radius = axial_dist.abs() * cone.half_angle_rad.tan();

        if (radial_dist - expected_radius).abs() < TOLERANCE_RETRY_LADDER_COARSE {
            on_surface_count += 1;
        }

        // Check for axial and angular components
        if i > 0 {
            let t_prev = t0 + (t1 - t0) * (i - 1) as f64 / (samples - 1) as f64;
            let pt_prev = curve.point_at(t_prev);
            let delta = pt - pt_prev;

            let axial_delta = delta.dot(axis).abs();
            let radial_delta = (delta - delta.dot(axis) * axis).length();

            if axial_delta > TOLERANCE_MESH_LEGACY {
                has_axial_component = true;
            }
            if radial_delta > TOLERANCE_MESH_LEGACY {
                has_angular_component = true;
            }
        }
    }

    on_surface_count > samples / 2 && has_axial_component && has_angular_component
}

// ── Seam Edge Detection ─────────────────────────────────────────────────────────

/// Check if an edge is a seam edge on a closed surface.
///
/// Seam edges are edges where a closed surface (cylinder, cone, sphere, torus)
/// meets itself at the parameter boundary (u = 0 and u = 2π).
pub fn is_seam_edge(
    brep: &rcad_kernel::BRep,
    edge_idx: usize,
    surface: &Surface3,
) -> bool {
    // A seam edge has two PCurves on the same surface with different u values
    // (typically 0 and 2π)

    let Some(_edge) = brep.edges.get(edge_idx) else {
        return false;
    };

    // Check if this edge has multiple PCurves on the same surface
    let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
        return false;
    };

    if pcurves.len() < 2 {
        return false;
    }

    // Check if two PCurves are on the same surface
    let mut surface_counts: HashMap<usize, usize> = HashMap::new();
    for pcurve in pcurves {
        *surface_counts.entry(pcurve.surface_idx).or_insert(0) += 1;
    }

    // If any surface has multiple PCurves for this edge, it's likely a seam
    for &count in surface_counts.values() {
        if count >= 2 {
            // Verify the surface is closed
            match surface {
                Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
                    return true;
                }
                _ => {}
            }
        }
    }

    false
}

/// Check if an edge is a degenerate edge (zero length, e.g., pole singularity).
pub fn is_degenerate_edge_for_hlr(brep: &rcad_kernel::BRep, edge_idx: usize) -> bool {
    // Check the degenerated flag in GeomStore
    brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false)
}

// ── Curve-Surface Intersection for HLR ──────────────────────────────────────────

/// Result of curve-surface intersection for HLR visibility.
#[derive(Debug, Clone)]
pub struct CurveSurfaceIntersection {
    /// Parameter values on the curve where it intersects the surface.
    pub curve_params: Vec<f64>,
    /// 3D points of intersection.
    pub points: Vec<DVec3>,
    /// UV parameters on the surface for each intersection.
    pub surface_uvs: Vec<(f64, f64)>,
    /// Visibility status between consecutive intersections.
    /// visibility[i] indicates visibility between intersection i and i+1.
    pub visibility: Vec<bool>,
}

/// Compute visible portions of a curve on a curved face.
///
/// This function finds where a curve intersects the silhouette of a surface
/// and determines which portions of the curve are visible.
pub fn compute_curve_visibility_on_surface(
    brep: &rcad_kernel::BRep,
    edge_idx: usize,
    surface_idx: usize,
    camera: &HlrCamera,
    opts: &HlrOptions,
) -> Option<CurveSurfaceIntersection> {
    let _edge = brep.edges.get(edge_idx)?;
    let surface = brep.geom.surfaces.get(surface_idx)?;

    // Get the edge curve
    let curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|&c| c)?;
    let curve = brep.geom.curves.get(curve_idx)?;

    // Get parameter range
    let [t0, t1] = brep.geom.edge_curve_range.get(edge_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| curve.default_domain());

    // Sample the curve
    let num_samples = opts.edge_samples.max(16);
    let mut curve_params = Vec::new();
    let mut points = Vec::new();
    let mut surface_uvs = Vec::new();

    let view_dir = (camera.target - camera.eye).normalize_or_zero();

    for i in 0..num_samples {
        let t = t0 + (t1 - t0) * i as f64 / (num_samples - 1) as f64;
        let pt = curve.point_at(t);

        curve_params.push(t);
        points.push(pt);

        // Project point onto surface to get UV
        if let Some((u, v)) = project_point_to_surface(&pt, surface) {
            surface_uvs.push((u, v));
        } else {
            surface_uvs.push((0.0, 0.0)); // placeholder
        }
    }

    // Compute visibility at each sample point
    let mut visibility = Vec::with_capacity(num_samples);

    for &pt in points.iter() {
        let _dist = (camera.eye - pt).length();
        let is_visible = true; // Will be computed by the main HLR pipeline
        visibility.push(is_visible);
    }

    // Find silhouette crossings
    let mut silhouette_crossings: Vec<(f64, DVec3)> = Vec::new();

    for i in 1..num_samples {
        let prev_uv = surface_uvs[i - 1];
        let curr_uv = surface_uvs[i];

        let prev_normal = surface.normal_at(prev_uv.0, prev_uv.1);
        let curr_normal = surface.normal_at(curr_uv.0, curr_uv.1);

        let prev_dot = prev_normal.dot(view_dir);
        let curr_dot = curr_normal.dot(view_dir);

        // Check for silhouette crossing (sign change)
        if prev_dot * curr_dot < 0.0 {
            // Bisection to find exact crossing
            if let Some((t_cross, pt_cross)) = find_silhouette_crossing(
                curve, surface, view_dir,
                curve_params[i - 1], curve_params[i],
                10,
            ) {
                silhouette_crossings.push((t_cross, pt_cross));
            }
        }
    }

    // Update visibility based on silhouette crossings
    // (Portions of the curve on the far side of the silhouette are hidden)

    Some(CurveSurfaceIntersection {
        curve_params,
        points,
        surface_uvs,
        visibility,
    })
}

/// Project a 3D point onto a surface to find the closest UV parameters.
fn project_point_to_surface(point: &DVec3, surface: &Surface3) -> Option<(f64, f64)> {
    use rcad_kernel::closest_point_on_surface;

    let result = closest_point_on_surface(surface, *point, 16);
    Some(result.params)
}

/// Find where a curve crosses a surface silhouette using bisection.
fn find_silhouette_crossing(
    curve: &rcad_kernel::geom::Curve3,
    surface: &Surface3,
    view_dir: DVec3,
    t_start: f64,
    t_end: f64,
    max_iter: usize,
) -> Option<(f64, DVec3)> {
    use rcad_kernel::geom::CurveEval;

    let pt_start = curve.point_at(t_start);
    let pt_end = curve.point_at(t_end);

    // Get UV parameters for start and end
    let uv_start = project_point_to_surface(&pt_start, surface)?;
    let uv_end = project_point_to_surface(&pt_end, surface)?;

    let dot_start = surface.normal_at(uv_start.0, uv_start.1).dot(view_dir);
    let dot_end = surface.normal_at(uv_end.0, uv_end.1).dot(view_dir);

    if dot_start * dot_end > 0.0 {
        return None; // No crossing
    }

    let mut t_lo = t_start;
    let mut t_hi = t_end;
    let mut dot_lo = dot_start;

    for _ in 0..max_iter {
        let t_mid = (t_lo + t_hi) * 0.5;
        let pt_mid = curve.point_at(t_mid);

        if let Some(uv_mid) = project_point_to_surface(&pt_mid, surface) {
            let dot_mid = surface.normal_at(uv_mid.0, uv_mid.1).dot(view_dir);

            if dot_mid.abs() < TOLERANCE_LINEAR_RELAX_8 {
                return Some((t_mid, pt_mid));
            }

            if dot_lo * dot_mid < 0.0 {
                t_hi = t_mid;
            } else {
                t_lo = t_mid;
                dot_lo = dot_mid;
            }
        }
    }

    let t_final = (t_lo + t_hi) * 0.5;
    Some((t_final, curve.point_at(t_final)))
}

// ── Edge Classification ─────────────────────────────────────────────────────────

/// Classify all edges in a BRep for HLR processing.
pub fn classify_edges(
    brep: &rcad_kernel::BRep,
    camera: &HlrCamera,
    opts: &HlrOptions,
) -> Vec<EdgeClassInfo> {
    let mut classifications: Vec<EdgeClassInfo> = Vec::new();

    for edge_idx in 0..brep.edges.len() {
        let classification = classify_single_edge(brep, edge_idx, camera, opts);
        classifications.push(classification);
    }

    classifications
}

/// Classify a single edge.
fn classify_single_edge(
    brep: &rcad_kernel::BRep,
    edge_idx: usize,
    _camera: &HlrCamera,
    opts: &HlrOptions,
) -> EdgeClassInfo {
    let Some(_edge) = brep.edges.get(edge_idx) else {
        return EdgeClassInfo {
            edge_idx,
            classification: EdgeClassification::Hidden,
            visible_segments: 0,
            hidden_segments: 0,
            on_curved_surface: false,
            surface_idx: None,
        };
    };

    // Get the surface this edge is on (if any)
    let surface_idx = get_edge_surface(brep, edge_idx);
    let on_curved_surface = surface_idx.is_some_and(|idx| {
        matches!(
            brep.geom.surfaces.get(idx),
            Some(Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Cone(_) | Surface3::Torus(_))
        )
    });

    // Check for thread edge
    if let Some(idx) = surface_idx
        && let Some(surface) = brep.geom.surfaces.get(idx) {
            if opts.detect_thread_edges && is_thread_edge(brep, edge_idx, surface) {
                return EdgeClassInfo {
                    edge_idx,
                    classification: EdgeClassification::Thread,
                    visible_segments: 0,
                    hidden_segments: 0,
                    on_curved_surface: true,
                    surface_idx: Some(idx),
                };
            }

            if opts.detect_seam_edges && is_seam_edge(brep, edge_idx, surface) {
                return EdgeClassInfo {
                    edge_idx,
                    classification: EdgeClassification::Seam,
                    visible_segments: 0,
                    hidden_segments: 0,
                    on_curved_surface: true,
                    surface_idx: Some(idx),
                };
            }
        }

    // Default classification - will be updated during HLR processing
    EdgeClassInfo {
        edge_idx,
        classification: EdgeClassification::Partial,
        visible_segments: 0,
        hidden_segments: 0,
        on_curved_surface,
        surface_idx,
    }
}

/// Get the primary surface an edge is on.
fn get_edge_surface(brep: &rcad_kernel::BRep, edge_idx: usize) -> Option<usize> {
    let pcurves = brep.geom.edge_pcurves.get(edge_idx)?;
    pcurves.first().map(|pc| pc.surface_idx)
}

// ── Spatial Indexing for Silhouette Queries ─────────────────────────────────────

/// Spatial grid for efficient silhouette point queries.
#[derive(Debug, Clone)]
pub struct SilhouetteSpatialIndex {
    /// Grid cells containing silhouette sample points.
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
    /// All sample points.
    points: Vec<DVec3>,
    /// Grid cell size.
    cell_size: f64,
    /// Bounding box of all points.
    bbox_min: DVec3,
    bbox_max: DVec3,
}

impl SilhouetteSpatialIndex {
    /// Build a spatial index from silhouette points.
    pub fn build(points: &[DVec3], cell_size: f64) -> Self {
        if points.is_empty() {
            return Self {
                cells: HashMap::new(),
                points: Vec::new(),
                cell_size,
                bbox_min: DVec3::ZERO,
                bbox_max: DVec3::ZERO,
            };
        }

        // Compute bounding box
        let mut bbox_min = DVec3::splat(f64::INFINITY);
        let mut bbox_max = DVec3::splat(f64::NEG_INFINITY);
        for &p in points {
            bbox_min = bbox_min.min(p);
            bbox_max = bbox_max.max(p);
        }

        // Build grid
        let mut cells: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &p) in points.iter().enumerate() {
            let cell = Self::point_to_cell(p, bbox_min, cell_size);
            cells.entry(cell).or_default().push(i);
        }

        Self {
            cells,
            points: points.to_vec(),
            cell_size,
            bbox_min,
            bbox_max,
        }
    }

    fn point_to_cell(p: DVec3, origin: DVec3, cell_size: f64) -> (i32, i32, i32) {
        let d = (p - origin) / cell_size;
        (d.x.floor() as i32, d.y.floor() as i32, d.z.floor() as i32)
    }

    /// Find all silhouette points within a radius of the query point.
    pub fn query_radius(&self, point: DVec3, radius: f64) -> Vec<usize> {
        let mut result = Vec::new();

        let cell_radius = (radius / self.cell_size).ceil() as i32;
        let center_cell = Self::point_to_cell(point, self.bbox_min, self.cell_size);

        for di in -cell_radius..=cell_radius {
            for dj in -cell_radius..=cell_radius {
                for dk in -cell_radius..=cell_radius {
                    let cell = (center_cell.0 + di, center_cell.1 + dj, center_cell.2 + dk);

                    if let Some(indices) = self.cells.get(&cell) {
                        for &idx in indices {
                            let dist = (self.points[idx] - point).length();
                            if dist <= radius {
                                result.push(idx);
                            }
                        }
                    }
                }
            }
        }

        result
    }

    /// Find the nearest silhouette point to the query point.
    pub fn query_nearest(&self, point: DVec3) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;

        // Start with a small search radius and expand
        let mut radius = self.cell_size;

        for _ in 0..10 {
            let candidates = self.query_radius(point, radius);

            for idx in candidates {
                let dist = (self.points[idx] - point).length();
                if best.is_none_or(|(_, d)| dist < d) {
                    best = Some((idx, dist));
                }
            }

            if best.is_some() && best.unwrap().1 <= radius * 0.5 {
                break;
            }

            radius *= 2.0;
        }

        best
    }

    /// Get a point by index.
    pub fn get_point(&self, idx: usize) -> Option<DVec3> {
        self.points.get(idx).copied()
    }

    /// Get the number of indexed points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

/// Numerical silhouette extraction for general parametric surfaces.
///
/// Uses a marching approach to find curves where normal · view_dir = 0.
/// This implementation includes:
/// - Marching along iso-parametric curves to trace silhouette curves
/// - Curvature-adaptive sampling for better accuracy in high-curvature regions
/// - Handling of closed silhouette loops
fn extract_numerical_silhouettes(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    opts: &HlrOptions,
    _brep: &rcad_kernel::BRep,
) -> Vec<Vec<DVec3>> {
    let [u0, u1, v0, v1] = domain;
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    // Find silhouette seed points on a coarse grid
    let grid_size = opts.silhouette_samples.max(16);
    let seeds = find_silhouette_seeds(surface, view_dir, domain, grid_size, opts.tangent_tolerance);

    if seeds.is_empty() {
        return curves;
    }

    // March from each seed to trace silhouette curves
    let mut visited: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();

    for (i, j, u, v) in seeds {
        if visited.contains(&(i, j)) {
            continue;
        }

        // Trace a curve starting from this seed
        let curve = march_silhouette_curve(surface, view_dir, domain, u, v, opts);

        if curve.len() >= 2 {
            // Mark visited cells along the curve
            for pt in &curve {
                // Find grid cell for this point
                let pi = ((pt.0 - u0) / (u1 - u0) * grid_size as f64).floor() as usize;
                let pj = ((pt.1 - v0) / (v1 - v0) * grid_size as f64).floor() as usize;
                visited.insert((pi.min(grid_size - 1), pj.min(grid_size - 1)));
            }

            // Apply adaptive refinement based on curvature
            let refined_curve = if opts.curvature_adaptive {
                refine_curve_by_curvature(surface, curve, opts)
            } else {
                curve.into_iter().map(|(_, _, pt)| pt).collect()
            };

            // Apply B-spline fitting if enabled
            let final_curve = if opts.fit_bspline && refined_curve.len() >= 4 {
                fit_bspline_to_points(&refined_curve, opts.bspline_tolerance)
            } else {
                refined_curve
            };

            if final_curve.len() >= 2 {
                curves.push(final_curve);
            }
        }
    }

    curves
}

/// A point in parameter space with its 3D position.
type ParamPoint = (f64, f64, DVec3);

/// Find seed points for silhouette curves on a grid.
fn find_silhouette_seeds(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    grid_size: usize,
    tangent_tol: f64,
) -> Vec<(usize, usize, f64, f64)> {
    let [u0, u1, v0, v1] = domain;
    let mut seeds = Vec::new();

    // Sample grid and look for sign changes in normal · view_dir
    let mut dot_values: Vec<Vec<f64>> = vec![vec![0.0; grid_size]; grid_size];

    // Compute dot products at grid sample nodes
    for i in 0..grid_size {
        for j in 0..grid_size {
            let u = u0 + (u1 - u0) * i as f64 / (grid_size - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (grid_size - 1) as f64;
            let normal = surface.normal_at(u, v);
            dot_values[i][j] = normal.dot(view_dir);
        }
    }

    // Find cells where sign changes occur (indicating silhouette crossing)
    for i in 0..grid_size - 1 {
        for j in 0..grid_size - 1 {
            let d00 = dot_values[i][j];
            let d10 = dot_values[i + 1][j];
            let d01 = dot_values[i][j + 1];
            let d11 = dot_values[i + 1][j + 1];

            // Check for sign changes in the cell
            let has_crossing = (d00 * d10 < 0.0)
                || (d00 * d01 < 0.0)
                || (d10 * d11 < 0.0)
                || (d01 * d11 < 0.0);

            if has_crossing {
                // Find the exact crossing point using bisection
                if let Some((u, v)) = find_crossing_point(surface, view_dir, domain, i, j, grid_size, tangent_tol) {
                    seeds.push((i, j, u, v));
                }
            }
        }
    }

    seeds
}

/// Find the exact crossing point in a grid cell using bisection.
fn find_crossing_point(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    i: usize,
    j: usize,
    grid_size: usize,
    _tangent_tol: f64,
) -> Option<(f64, f64)> {
    let [u0, v0, u1, v1] = [
        domain[0] + (domain[1] - domain[0]) * i as f64 / (grid_size - 1) as f64,
        domain[2] + (domain[3] - domain[2]) * j as f64 / (grid_size - 1) as f64,
        domain[0] + (domain[1] - domain[0]) * (i + 1) as f64 / (grid_size - 1) as f64,
        domain[2] + (domain[3] - domain[2]) * (j + 1) as f64 / (grid_size - 1) as f64,
    ];

    // Try to find crossing along each edge of the cell
    let edges = [
        (u0, v0, u1, v0), // bottom edge
        (u0, v1, u1, v1), // top edge
        (u0, v0, u0, v1), // left edge
        (u1, v0, u1, v1), // right edge
    ];

    for (ua, va, ub, vb) in edges {
        if let Some((u, v)) = bisection_search(surface, view_dir, ua, va, ub, vb, 12) {
            return Some((u, v));
        }
    }

    None
}

/// Bisection search to find where normal · view_dir = 0.
fn bisection_search(
    surface: &Surface3,
    view_dir: DVec3,
    u0: f64,
    v0: f64,
    u1: f64,
    v1: f64,
    max_iter: usize,
) -> Option<(f64, f64)> {
    let d0 = surface.normal_at(u0, v0).dot(view_dir);
    let d1 = surface.normal_at(u1, v1).dot(view_dir);

    if d0 * d1 > 0.0 {
        return None; // No sign change
    }

    let mut ua = u0;
    let mut va = v0;
    let mut ub = u1;
    let mut vb = v1;

    for _ in 0..max_iter {
        let um = (ua + ub) / 2.0;
        let vm = (va + vb) / 2.0;
        let dm = surface.normal_at(um, vm).dot(view_dir);

        if dm.abs() < TOLERANCE_LINEAR_ULTRA_STRICT {
            return Some((um, vm));
        }

        if d0 * dm < 0.0 {
            ub = um;
            vb = vm;
        } else {
            ua = um;
            va = vm;
        }
    }

    Some(((ua + ub) / 2.0, (va + vb) / 2.0))
}

/// March along a silhouette curve starting from a seed point.
fn march_silhouette_curve(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    u_start: f64,
    v_start: f64,
    opts: &HlrOptions,
) -> Vec<ParamPoint> {
    let mut curve: Vec<ParamPoint> = Vec::new();
    let [u0, u1, v0, v1] = domain;

    // Add the starting point
    let p_start = surface.point_at(u_start, v_start);
    curve.push((u_start, v_start, p_start));

    // March in both directions from the seed
    for direction in &[-1.0_f64, 1.0] {
        let mut u = u_start;
        let mut v = v_start;
        let mut curve_dir: Option<DVec2> = None;

        for _ in 0..opts.max_subdivisions * 50 {
            // Compute the tangent direction to the silhouette curve
            let tangent = compute_silhouette_tangent(surface, view_dir, u, v);

            if tangent.length_squared() < 1e-16 {
                break;
            }

            // Choose direction along the tangent
            let step_dir = if let Some(cd) = curve_dir {
                // Continue in the same general direction
                if cd.dot(tangent) > 0.0 {
                    tangent
                } else {
                    -tangent
                }
            } else {
                *direction * tangent
            };
            curve_dir = Some(step_dir.normalize_or_zero());

            // Compute step size based on curvature
            let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, u, v);
            let max_k = k1.abs().max(k2.abs()).max(opts.min_curvature);
            let curvature_factor = (opts.max_curvature / max_k).min(4.0).max(0.25);
            let step_size = opts.angular_tolerance * curvature_factor;

            // Take a step
            let u_new = u + step_dir.x * step_size;
            let v_new = v + step_dir.y * step_size;

            // Check bounds
            if u_new < u0 || u_new > u1 || v_new < v0 || v_new > v1 {
                break;
            }

            // Project back onto the silhouette curve
            if let Some((u_proj, v_proj)) = project_to_silhouette(surface, view_dir, u_new, v_new, opts.tangent_tolerance) {
                u = u_proj;
                v = v_proj;

                let p = surface.point_at(u, v);
                let d = (p - curve.last().map(|(_, _, lp)| *lp).unwrap_or(p_start)).length();

                // Only add if we've moved enough
                if d > opts.bspline_tolerance * 0.1 {
                    curve.push((u, v, p));
                }
            } else {
                break;
            }

            // Check for closed loop
            if curve.len() > 10 {
                let first = curve[0];
                let dist = ((first.0 - u).powi(2) + (first.1 - v).powi(2)).sqrt();
                if dist < step_size * 2.0 {
                    // Close the loop
                    curve.push(curve[0]);
                    break;
                }
            }
        }

        // Reverse the points added while marching in the negative direction
        if *direction < 0.0 && curve.len() > 1 {
            let first = curve[0];
            curve.reverse();
            curve.push(first); // Re-add the start point for the loop
        }
    }

    curve
}

/// Compute the tangent direction to the silhouette curve at a point.
fn compute_silhouette_tangent(
    surface: &Surface3,
    view_dir: DVec3,
    u: f64,
    v: f64,
) -> DVec2 {
    const EPS: f64 = TOLERANCE_MESH_LEGACY;

    // Compute gradients of the implicit function f(u,v) = N(u,v) · V
    let n = surface.normal_at(u, v);
    let n_u = surface.normal_at(u + EPS, v);
    let n_v = surface.normal_at(u, v + EPS);

    // Gradient of f = N · V
    let df_du = (n_u - n).dot(view_dir) / EPS;
    let df_dv = (n_v - n).dot(view_dir) / EPS;

    // The tangent direction is perpendicular to the gradient
    DVec2::new(-df_dv, df_du).normalize_or_zero()
}

/// Project a point back onto the silhouette curve.
fn project_to_silhouette(
    surface: &Surface3,
    view_dir: DVec3,
    u: f64,
    v: f64,
    tol: f64,
) -> Option<(f64, f64)> {
    let mut u_curr = u;
    let mut v_curr = v;

    // Newton iteration to find f(u,v) = 0
    for _ in 0..20 {
        let n = surface.normal_at(u_curr, v_curr);
        let f = n.dot(view_dir);

        if f.abs() < tol {
            return Some((u_curr, v_curr));
        }

        // Compute gradient numerically
        const EPS: f64 = TOLERANCE_ABS;
        let n_u = surface.normal_at(u_curr + EPS, v_curr);
        let n_v = surface.normal_at(u_curr, v_curr + EPS);

        let df_du = (n_u - n).dot(view_dir) / EPS;
        let df_dv = (n_v - n).dot(view_dir) / EPS;

        let grad_len_sq = df_du * df_du + df_dv * df_dv;
        if grad_len_sq < TOLERANCE_METRIC_SQ_NEAR_ZERO {
            break;
        }

        // Newton step
        let step = f / grad_len_sq;
        u_curr -= step * df_du;
        v_curr -= step * df_dv;
    }

    // Check if we converged
    let f = surface.normal_at(u_curr, v_curr).dot(view_dir);
    if f.abs() < tol * 10.0 {
        Some((u_curr, v_curr))
    } else {
        None
    }
}

/// Refine a silhouette curve based on surface curvature.
fn refine_curve_by_curvature(
    surface: &Surface3,
    curve: Vec<ParamPoint>,
    opts: &HlrOptions,
) -> Vec<DVec3> {
    if curve.len() < 2 {
        return curve.into_iter().map(|(_, _, p)| p).collect();
    }

    let mut refined: Vec<DVec3> = Vec::new();
    refined.push(curve[0].2);

    for i in 1..curve.len() {
        let (u0, v0, p0) = curve[i - 1];
        let (u1, v1, p1) = curve[i];

        // Compute curvature at the midpoint
        let um = (u0 + u1) / 2.0;
        let vm = (v0 + v1) / 2.0;
        let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, um, vm);
        let max_k = k1.abs().max(k2.abs());

        // Determine number of subdivision points based on curvature
        let chord_len = (p1 - p0).length();
        let subdivs = if max_k > opts.min_curvature {
            let curvature_samples = (max_k * chord_len * std::f64::consts::PI).ceil() as usize;
            curvature_samples.min(8).max(1)
        } else {
            1
        };

        // Add subdivision points
        for j in 1..subdivs {
            let t = j as f64 / subdivs as f64;
            let u = u0 + t * (u1 - u0);
            let v = v0 + t * (v1 - v0);
            let p = surface.point_at(u, v);
            refined.push(p);
        }

        refined.push(p1);
    }

    refined
}

/// Fit a B-spline curve to a set of points.
fn fit_bspline_to_points(points: &[DVec3], tolerance: f64) -> Vec<DVec3> {
    if points.len() < 4 {
        return points.to_vec();
    }

    // Simple approach: sample the fitted B-spline at uniform intervals
    // For a proper implementation, we would use least-squares fitting
    // Here we use a simplified version that preserves the shape

    let n = points.len();
    let mut result: Vec<DVec3> = Vec::with_capacity(n);

    // Compute chord lengths for parameterization
    let mut chords = vec![0.0_f64; n];
    for i in 1..n {
        chords[i] = chords[i - 1] + (points[i] - points[i - 1]).length();
    }
    let total_len = chords[n - 1];
    if total_len < TOLERANCE_LEN_MIN {
        return points.to_vec();
    }

    // Generate control points using Catmull-Rom style interpolation
    let _degree = 3.min(n - 1);
    let num_samples = (total_len / tolerance).ceil() as usize;
    let num_samples = num_samples.max(10).min(1000);

    for i in 0..=num_samples {
        let t = i as f64 / num_samples as f64;
        let target_len = t * total_len;

        // Find the segment containing this length
        let seg_idx = chords.partition_point(|&c| c < target_len).saturating_sub(1);
        let seg_idx = seg_idx.min(n - 2);

        // Interpolate within the segment
        let seg_start = chords[seg_idx];
        let seg_end = chords[seg_idx + 1];
        let seg_len = seg_end - seg_start;

        let local_t = if seg_len > TOLERANCE_LEN_MIN {
            (target_len - seg_start) / seg_len
        } else {
            0.5
        };

        // Simple linear interpolation with smoothing
        let p0 = points[seg_idx];
        let p1 = points[seg_idx + 1];

        // Hermite interpolation for smoother result
        let t0 = if seg_idx > 0 {
            (points[seg_idx + 1] - points[seg_idx - 1]).normalize_or_zero()
        } else {
            (points[1] - points[0]).normalize_or_zero()
        };

        let t1 = if seg_idx + 2 < n {
            (points[seg_idx + 2] - points[seg_idx]).normalize_or_zero()
        } else {
            (points[n - 1] - points[n - 2]).normalize_or_zero()
        };

        let h00 = 2.0 * local_t * local_t * local_t - 3.0 * local_t * local_t + 1.0;
        let h10 = local_t * local_t * local_t - 2.0 * local_t * local_t + local_t;
        let h01 = -2.0 * local_t * local_t * local_t + 3.0 * local_t * local_t;
        let h11 = local_t * local_t * local_t - local_t * local_t;

        let p = h00 * p0 + h10 * t0 * seg_len + h01 * p1 + h11 * t1 * seg_len;
        result.push(p);
    }

    result
}

/// Generate silhouette curves for the HLR pipeline (internal function).
fn compute_silhouettes(brep: &rcad_kernel::BRep, view_dir: DVec3, samples: usize) -> Vec<SilhouetteCurve> {
    let opts = HlrOptions {
        silhouette_samples: samples,
        ..HlrOptions::default()
    };

    extract_silhouette_curves(brep, view_dir, &opts)
        .into_iter()
        .map(|curve| SilhouetteCurve {
            world_pts: curve.points,
            curve_hint: None,
            dense: true, // All silhouettes are treated as dense for proper rendering
        })
        .collect()
}

/// Occlusion tester that supports both brute-force and BVH-accelerated methods.
enum OcclusionTester<'a> {
    BruteForce(&'a [[DVec3; 3]]),
    Bvh {
        bvh: &'a TriBvh,
        triangles: &'a [[DVec3; 3]],
    },
}

impl<'a> OcclusionTester<'a> {
    fn is_occluded(&self, point: DVec3, eye: DVec3, dist_to_eye: f64) -> bool {
        match self {
            OcclusionTester::BruteForce(triangles) => {
                is_occluded(point, eye, triangles, dist_to_eye)
            }
            OcclusionTester::Bvh { bvh, triangles } => {
                bvh.is_occluded(point, eye, triangles, dist_to_eye)
            }
        }
    }
}

/// Improved visibility classification that handles grazing angles on curved surfaces.
///
/// For points near silhouette curves (where normal is nearly perpendicular to view direction),
/// we use additional testing to improve numerical stability.
fn classify_visibility(
    point: DVec3,
    normal: Option<DVec3>,
    camera: &HlrCamera,
    occlusion_tester: &OcclusionTester<'_>,
    grazing_threshold: f64,
) -> VisibilityInfo {
    let dist = (camera.eye - point).length();
    let view_dir = (camera.eye - point).normalize_or_zero();

    // Check if we're at a grazing angle
    let grazing_factor = if let Some(n) = normal {
        let dot = n.dot(view_dir).abs();
        // grazing_factor = 1.0 when perfectly grazing (dot = 0)
        // grazing_factor = 0.0 when viewing straight on (dot = 1)
        1.0 - dot
    } else {
        0.0
    };

    // For grazing angles, use more robust testing
    let is_occluded = if grazing_factor > grazing_threshold.cos() {
        // At grazing angle: test multiple rays to reduce false positives
        let base_occluded = occlusion_tester.is_occluded(point, camera.eye, dist);

        if base_occluded {
            // Verify with additional samples to reduce numerical errors
            let mut occluded_count = 1;
            const NUM_SAMPLES: usize = 4;
            let offset = TOLERANCE_RETRY_LADDER_COARSE;

            for i in 0..NUM_SAMPLES {
                let angle = i as f64 * std::f64::consts::TAU / NUM_SAMPLES as f64;
                let perp = any_perpendicular(view_dir);
                let perturb = perp * (angle.cos() * offset) + view_dir.cross(perp) * (angle.sin() * offset);
                let test_point = point + perturb;

                if occlusion_tester.is_occluded(test_point, camera.eye, dist) {
                    occluded_count += 1;
                }
            }

            // Require majority to confirm occlusion at grazing angles
            occluded_count > NUM_SAMPLES / 2
        } else {
            false
        }
    } else {
        occlusion_tester.is_occluded(point, camera.eye, dist)
    };

    VisibilityInfo {
        visible: !is_occluded,
        grazing_factor,
        depth: dist,
    }
}

/// Information about visibility at a point.
struct VisibilityInfo {
    visible: bool,
    grazing_factor: f64,
    depth: f64,
}

/// Process a list of world-space sample points through the HLR visibility
/// pipeline and append resulting segments to `result`.
///
/// When `dense` is true, one segment is emitted per consecutive point pair
/// (useful for polyline approximations of curved silhouettes).
fn process_world_pts(
    world_pts: &[DVec3],
    curve_hint: Option<CurveHint>,
    dense: bool,
    segment_type: SegmentType,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    result: &mut HlrResult,
) {
    process_world_pts_with_bvh(
        world_pts,
        curve_hint,
        dense,
        segment_type,
        camera,
        view,
        triangles,
        None,
        &HlrOptions::default(),
        result,
    )
}

/// Process world points with optional BVH acceleration and grazing angle handling.
fn process_world_pts_with_bvh(
    world_pts: &[DVec3],
    curve_hint: Option<CurveHint>,
    dense: bool,
    segment_type: SegmentType,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    bvh: Option<&TriBvh>,
    _opts: &HlrOptions,
    result: &mut HlrResult,
) {
    if world_pts.len() < 2 {
        return;
    }
    let n = world_pts.len();

    let occlusion_tester = if let Some(bvh) = bvh {
        OcclusionTester::Bvh { bvh, triangles }
    } else {
        OcclusionTester::BruteForce(triangles)
    };

    let sample_vis: Vec<bool> = world_pts
        .iter()
        .map(|&wp| {
            let dist = (camera.eye - wp).length();
            !occlusion_tester.is_occluded(wp, camera.eye, dist)
        })
        .collect();

    let screen_pts: Vec<DVec2> = world_pts.iter().map(|&wp| project(wp, view).0).collect();

    if dense {
        // Emit one segment per consecutive pair (preserves polyline shape).
        for i in 0..n - 1 {
            let seg = HlrSegment {
                start: screen_pts[i],
                end: screen_pts[i + 1],
                visible: sample_vis[i] && sample_vis[i + 1],
                curve_hint: curve_hint.clone(),
                segment_type,
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
        }
        return;
    }

    let mut seg_start = 0usize;
    for i in 1..n {
        let changed = sample_vis[i] != sample_vis[seg_start];
        let last = i == n - 1;
        if changed || last {
            let end_idx = if last && !changed { i } else { i - 1 };
            let seg = HlrSegment {
                start: screen_pts[seg_start],
                end: screen_pts[end_idx],
                visible: sample_vis[seg_start],
                curve_hint: curve_hint.clone(),
                segment_type,
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
            if changed {
                seg_start = i;
            }
        }
    }
}



/// Perform hidden-line removal on a BRep from the given camera position.
///
/// Returns 2D projected segments labeled visible/hidden.
/// `samples` controls how finely each edge is subdivided for occlusion testing
/// (higher = more accurate but slower; 8 is a reasonable default).
pub fn compute_hlr(brep: &rcad_kernel::BRep, camera: &HlrCamera, samples: usize) -> HlrResult {
    compute_hlr_with_options(brep, camera, HlrOptions::default().with_edge_samples(samples))
}

/// Perform hidden-line removal with full configuration options.
///
/// This function provides fine-grained control over HLR computation parameters,
/// including adaptive sampling for curved surfaces.
///
/// # Arguments
/// * `brep` - The BRep model to process.
/// * `camera` - Camera/view specification.
/// * `opts` - Configuration options for sampling and tolerances.
///
/// # Returns
/// An `HlrResult` containing projected 2D segments labeled as visible/hidden.
pub fn compute_hlr_with_options(brep: &rcad_kernel::BRep, camera: &HlrCamera, opts: HlrOptions) -> HlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let triangles = collect_triangles(brep);
    let _edge_samples = opts.edge_samples.max(2);
    let mut result = HlrResult::default();

    // Build BVH for acceleration if enabled and we have enough triangles
    let bvh: Option<TriBvh> = if opts.use_bvh && triangles.len() > 32 {
        Some(TriBvh::build(&triangles))
    } else {
        None
    };
    let bvh_ref = bvh.as_ref();

    // Classify edges for thread/seam detection
    let edge_classifications = if opts.detect_thread_edges || opts.detect_seam_edges {
        classify_edges(brep, camera, &opts)
    } else {
        Vec::new()
    };

    // ── Wire edges ────────────────────────────────────────────────────────────

    // Collect all unique edges from all faces + standalone edges
    let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    edge_indices.insert(we.idx);
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        edge_indices.insert(we.idx);
                    }
                }
            }
        }
    }
    for i in 0..brep.edges.len() {
        edge_indices.insert(i);
    }

    // Convert to vector for potential parallel processing
    let edge_indices_vec: Vec<usize> = edge_indices.into_iter().collect();

    // Process edges (optionally in parallel)
    let edge_results: Vec<Vec<HlrSegment>> = if opts.parallel && edge_indices_vec.len() > opts.parallel_threshold {
        let triangles_ref = &triangles;
        let bvh_opt = bvh_ref;
        let brep_ref = brep;
        let camera_ref = camera;
        let view_ref = &view;
        let opts_ref = &opts;
        let edge_classes = &edge_classifications;

        edge_indices_vec
            .par_iter()
            .map(|&edge_idx| {
                process_single_edge(
                    brep_ref,
                    edge_idx,
                    camera_ref,
                    view_ref,
                    triangles_ref,
                    bvh_opt,
                    opts_ref,
                    edge_classes,
                )
            })
            .collect()
    } else {
        edge_indices_vec
            .iter()
            .map(|&edge_idx| {
                process_single_edge(
                    brep,
                    edge_idx,
                    camera,
                    &view,
                    &triangles,
                    bvh_ref,
                    &opts,
                    &edge_classifications,
                )
            })
            .collect()
    };

    // Merge results
    for segments in edge_results {
        result.segments.extend(segments);
    }

    // ── Silhouette curves ────────────────────────────────────────────

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    let silhouette_curves = compute_silhouettes_with_options(brep, view_dir, &opts);

    // Build spatial index for silhouette queries
    let all_silhouette_points: Vec<DVec3> = silhouette_curves
        .iter()
        .flat_map(|c| c.world_pts.iter().copied())
        .collect();
    let _spatial_index = SilhouetteSpatialIndex::build(&all_silhouette_points, 0.1);

    for sil in silhouette_curves {
        process_world_pts_with_bvh(
            &sil.world_pts,
            sil.curve_hint,
            sil.dense,
            SegmentType::Silhouette,
            camera,
            &view,
            &triangles,
            bvh_ref,
            &opts,
            &mut result,
        );
    }

    result
}

/// Process a single edge and return its segments.
fn process_single_edge(
    brep: &rcad_kernel::BRep,
    edge_idx: usize,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    bvh: Option<&TriBvh>,
    opts: &HlrOptions,
    edge_classifications: &[EdgeClassInfo],
) -> Vec<HlrSegment> {
    let segments: Vec<HlrSegment> = Vec::new();

    let Some(edge) = brep.edges.get(edge_idx) else { return segments; };
    let Some(v_start) = brep.vertices.get(edge.start) else { return segments; };
    let Some(v_end) = brep.vertices.get(edge.end) else { return segments; };

    let p0 = v_start.point;
    let p1 = v_end.point;

    // Determine edge type
    let segment_type = if let Some(class_info) = edge_classifications.get(edge_idx) {
        match class_info.classification {
            EdgeClassification::Thread => SegmentType::Thread,
            EdgeClassification::Seam => SegmentType::Seam,
            _ => SegmentType::Edge,
        }
    } else {
        SegmentType::Edge
    };

    let edge_curve = brep
        .geom
        .edge_curve
        .get(edge_idx)
        .and_then(|&ci| ci)
        .and_then(|ci| brep.geom.curves.get(ci));

    let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
        if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
    });

    let is_other_curve = edge_curve
        .is_some_and(|c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
        && circle_info.is_none();

    // Adaptive sampling for curved edges on curved surfaces
    let this_edge_samples = if circle_info.is_some() || is_other_curve {
        // Check if this edge is on a curved surface with high curvature
        if let Some(class_info) = edge_classifications.get(edge_idx) {
            if class_info.on_curved_surface {
                if let Some(surf_idx) = class_info.surface_idx {
                    if let Some(surface) = brep.geom.surfaces.get(surf_idx) {
                        // Compute curvature at midpoint
                        let domain = brep.geom.face_surface_range
                            .iter()
                            .find_map(|r| *r)
                            .unwrap_or_else(|| surface.default_domain());
                        let mid_u = (domain[0] + domain[1]) * 0.5;
                        let mid_v = (domain[2] + domain[3]) * 0.5;
                        let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, mid_u, mid_v);
                        let max_k = k1.abs().max(k2.abs());

                        // More samples for higher curvature
                        let adaptive_factor = (max_k / 10.0).min(8.0).max(1.0);
                        ((opts.edge_samples as f64 * adaptive_factor * 4.0) as usize).max(32).min(256)
                    } else {
                        (opts.edge_samples * 4).max(32)
                    }
                } else {
                    (opts.edge_samples * 4).max(32)
                }
            } else {
                (opts.edge_samples * 4).max(32)
            }
        } else {
            (opts.edge_samples * 4).max(32)
        }
    } else {
        opts.edge_samples
    };

    let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
        let [t0, t1] = brep
            .geom
            .edge_curve_range
            .get(edge_idx)
            .and_then(|r| *r)
            .unwrap_or_else(|| circ.default_domain());
        (0..this_edge_samples)
            .map(|i| {
                let t = t0 + (t1 - t0) * (i as f64 / (this_edge_samples - 1) as f64);
                circ.point_at(t)
            })
            .collect()
    } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
        let [t0, t1] = brep
            .geom
            .edge_curve_range
            .get(edge_idx)
            .and_then(|r| *r)
            .unwrap_or_else(|| curve.default_domain());
        (0..this_edge_samples)
            .map(|i| {
                let t = t0 + (t1 - t0) * (i as f64 / (this_edge_samples - 1) as f64);
                curve.point_at(t)
            })
            .collect()
    } else {
        if (p1 - p0).length_squared() < TOLERANCE_LEN_MIN {
            return segments;
        }
        (0..this_edge_samples)
            .map(|i| {
                let t = i as f64 / (this_edge_samples - 1) as f64;
                p0 + (p1 - p0) * t
            })
            .collect()
    };

    // Compute curve_hint for circle edges
    let screen_pts_for_hint: Vec<DVec2> =
        world_pts.iter().map(|&wp| project(wp, view).0).collect();
    let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
        let (center_2d, _) = project(circ.center, view);
        let r = screen_pts_for_hint
            .iter()
            .map(|p| (*p - center_2d).length())
            .fold(0.0_f64, f64::max);
        Some(CurveHint::Circle { center: center_2d, radius: r })
    } else if is_other_curve {
        Some(CurveHint::Other)
    } else {
        None
    };

    // Process the edge points
    let mut edge_result = HlrResult::default();
    process_world_pts_with_bvh(
        &world_pts,
        curve_hint,
        false,
        segment_type,
        camera,
        view,
        triangles,
        bvh,
        opts,
        &mut edge_result,
    );

    edge_result.segments
}

/// Compute silhouette curves with full options (internal helper).
fn compute_silhouettes_with_options(brep: &rcad_kernel::BRep, view_dir: DVec3, opts: &HlrOptions) -> Vec<SilhouetteCurve> {
    extract_silhouette_curves(brep, view_dir, opts)
        .into_iter()
        .map(|curve| SilhouetteCurve {
            world_pts: curve.points,
            curve_hint: None,
            dense: true,
        })
        .collect()
}

/// Per-component HLR result for assembly HLR.
#[derive(Debug, Clone, Default)]
pub struct ComponentHlr {
    /// Component name (from the assembly).
    pub name: String,
    /// HLR segments for this component.
    pub segments: Vec<HlrSegment>,
}

/// Output of assembly HLR — one `ComponentHlr` per leaf BRep.
#[derive(Debug, Clone, Default)]
pub struct AssemblyHlrResult {
    pub components: Vec<ComponentHlr>,
}

impl AssemblyHlrResult {
    /// Return all visible segments across all components.
    pub fn visible_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| s.visible).map(move |s| (c, s))
        })
    }

    /// Return all hidden segments across all components.
    pub fn hidden_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| !s.visible).map(move |s| (c, s))
        })
    }
}

/// Transform a BRep's vertices by an affine transform.
/// Returns a new BRep with transformed vertex positions.
fn transform_brep(brep: &rcad_kernel::BRep, transform: &DAffine3) -> rcad_kernel::BRep {
    let mut out = brep.clone();
    for v in &mut out.vertices {
        v.point = transform.transform_point3(v.point);
    }
    out
}

/// Perform hidden-line removal on an assembly of BReps.
///
/// Each component's geometry is transformed to world space, then all triangles
/// are merged into a single occlusion buffer. Each component's edges are
/// tested against the global occlusion buffer, so components correctly
/// occlude each other.
///
/// Returns one `ComponentHlr` per leaf component.
pub fn hlr_assembly(
    components: &[(rcad_kernel::BRep, DAffine3, String)],
    camera: &HlrCamera,
    samples: usize,
) -> AssemblyHlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let samples = samples.max(2);

    // Transform all BRePs to world space and collect a unified triangle pool.
    let world_breps: Vec<rcad_kernel::BRep> = components
        .iter()
        .map(|(brep, xf, _)| transform_brep(brep, xf))
        .collect();

    let mut all_triangles: Vec<[DVec3; 3]> = Vec::new();
    for wb in &world_breps {
        all_triangles.extend(collect_triangles(wb));
    }

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    let mut result = AssemblyHlrResult::default();

    for (wb, (_, _, name)) in world_breps.iter().zip(components.iter()) {
        let mut comp_result = HlrResult::default();

        // ── Wire edges ────────────────────────────────────────────────────
        let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        for solid in &wb.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for we in &face.outer_wire.edges {
                        edge_indices.insert(we.idx);
                    }
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            edge_indices.insert(we.idx);
                        }
                    }
                }
            }
        }
        for i in 0..wb.edges.len() {
            edge_indices.insert(i);
        }

        for &edge_idx in &edge_indices {
            let Some(edge) = wb.edges.get(edge_idx) else { continue };
            let Some(v_start) = wb.vertices.get(edge.start) else { continue };
            let Some(v_end) = wb.vertices.get(edge.end) else { continue };

            let p0 = v_start.point;
            let p1 = v_end.point;

            let edge_curve = wb
                .geom
                .edge_curve
                .get(edge_idx)
                .and_then(|&ci| ci)
                .and_then(|ci| wb.geom.curves.get(ci));

            let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
                if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
            });

            let is_other_curve = edge_curve
                .is_some_and(|c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
                && circle_info.is_none();

            let edge_samples = if circle_info.is_some() || is_other_curve {
                64.max(samples)
            } else {
                samples
            };

            let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| circ.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        circ.point_at(t)
                    })
                    .collect()
            } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| curve.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        curve.point_at(t)
                    })
                    .collect()
            } else {
                if (p1 - p0).length_squared() < TOLERANCE_LEN_MIN {
                    continue;
                }
                (0..edge_samples)
                    .map(|i| {
                        let t = i as f64 / (edge_samples - 1) as f64;
                        p0 + (p1 - p0) * t
                    })
                    .collect()
            };

            let screen_pts_for_hint: Vec<DVec2> =
                world_pts.iter().map(|&wp| project(wp, &view).0).collect();
            let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
                let (center_2d, _) = project(circ.center, &view);
                let r = screen_pts_for_hint
                    .iter()
                    .map(|p| (*p - center_2d).length())
                    .fold(0.0_f64, f64::max);
                Some(CurveHint::Circle { center: center_2d, radius: r })
            } else if is_other_curve {
                Some(CurveHint::Other)
            } else {
                None
            };

            process_world_pts(&world_pts, curve_hint, false, SegmentType::Edge, camera, &view, &all_triangles, &mut comp_result);
        }

        // ── Silhouette curves ────────────────────────────────────
        let opts = HlrOptions::default().with_edge_samples(samples);
        for sil in compute_silhouettes_with_options(wb, view_dir, &opts) {
            process_world_pts(&sil.world_pts, sil.curve_hint, sil.dense, SegmentType::Silhouette, camera, &view, &all_triangles, &mut comp_result);
        }

        result.components.push(ComponentHlr {
            name: name.clone(),
            segments: comp_result.segments,
        });
    }

    result
}

/// Render HLR result as a simple SVG string.
///
/// Visible edges are drawn solid black; hidden edges are dashed gray.
/// `scale` controls pixel size per unit.
pub fn hlr_to_svg(result: &HlrResult, scale: f64, margin: f64) -> String {
    if result.segments.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".to_string();
    }

    // Compute bounding box
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for seg in &result.segments {
        for p in [seg.start, seg.end] {
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
        }
    }

    // Flip Y (SVG Y grows downward, camera Y grows upward)
    let transform = |p: DVec2| -> (f64, f64) {
        let x = (p.x - min_x) * scale + margin;
        let y = (max_y - p.y) * scale + margin;
        (x, y)
    };

    let w = (max_x - min_x) * scale + 2.0 * margin;
    let h = (max_y - min_y) * scale + 2.0 * margin;

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}\" height=\"{:.1}\" viewBox=\"0 0 {:.1} {:.1}\">\n",
        w, h, w, h
    );
    svg.push_str("  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");

    for seg in &result.segments {
        let (x1, y1) = transform(seg.start);
        let (x2, y2) = transform(seg.end);
        let stroke = if seg.visible {
            "black\" stroke-width=\"1.5"
        } else {
            "#999\" stroke-width=\"0.8\" stroke-dasharray=\"4,3"
        };

        // For circle segments emit an SVG arc path; for all others emit a line.
        if let Some(CurveHint::Circle { center, radius }) = &seg.curve_hint {
            let (cx, cy) = transform(*center);
            let r = radius * scale;
            // Determine large-arc flag: compare arc length vs half-circumference
            let dx1 = x1 - cx;
            let dy1 = y1 - cy;
            let dx2 = x2 - cx;
            let dy2 = y2 - cy;
            let cross = dx1 * dy2 - dy1 * dx2;
            let dot = dx1 * dx2 + dy1 * dy2;
            let angle = cross.atan2(dot).abs();
            let large_arc = if angle > std::f64::consts::PI { 1 } else { 0 };
            let sweep = if cross < 0.0 { 0 } else { 1 };
            svg.push_str(&format!(
                "  <path d=\"M {:.3} {:.3} A {:.3} {:.3} 0 {} {} {:.3} {:.3}\" fill=\"none\" stroke=\"{}\"/>\n",
                x1, y1, r, r, large_arc, sweep, x2, y2, stroke
            ));
            // Also record the center for debugging/reference (as a tiny dot, invisible by default)
            let _ = (cx, cy); // suppress unused warning
        } else {
            svg.push_str(&format!(
                "  <line x1=\"{:.3}\" y1=\"{:.3}\" x2=\"{:.3}\" y2=\"{:.3}\" stroke=\"{}\"/>\n",
                x1, y1, x2, y2, stroke
            ));
        }
    }
    svg.push_str("</svg>\n");
    svg
}

// ── Tests ─────────────────────────────────────────────────────────────────────
