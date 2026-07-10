use std::sync::Arc;

// =============================================================================
// ENHANCED DEFEATURE: THROUGH-HOLE vs BLIND-HOLE DETECTION
// =============================================================================

/// Hole type classification based on geometry analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HoleType {
    /// Through-hole: the hole passes completely through the solid.
    ThroughHole,
    /// Blind hole: the hole has a closed bottom.
    BlindHole,
    /// Counterbore: a stepped hole with a larger diameter section.
    Counterbore,
    /// Countersink: a conical enlargement at the top of a hole.
    Countersink,
    /// Spotface: a shallow circular recess for washer/bolt head seating.
    Spotface,
    /// Unknown: unable to classify.
    Unknown,
}

/// Extended cylindrical feature with additional classification.
#[derive(Debug, Clone)]
pub struct CylindricalFeatureExtended {
    /// Base cylindrical feature.
    pub base: CylindricalFeature,
    /// Hole type classification.
    pub hole_type: HoleType,
    /// Whether the hole has a flat bottom (typical for blind holes).
    pub has_flat_bottom: bool,
    /// Whether the hole bottom is conical.
    pub has_conical_bottom: bool,
    /// Estimated depth for blind holes (0.0 for through-holes).
    pub blind_depth: f64,
    /// Face index of the bottom face (if blind hole).
    pub bottom_face_index: Option<usize>,
    /// Adjacent face indices at top and bottom openings.
    pub top_adjacent_faces: Vec<usize>,
    pub bottom_adjacent_faces: Vec<usize>,
}

/// Classify a cylindrical feature as through-hole or blind-hole.
///
/// This function analyzes the topology around a cylindrical feature to determine
/// whether it passes completely through the solid or has a closed bottom.
///
/// # Algorithm
///
/// 1. Find all faces adjacent to the cylindrical wall face(s) at each end
/// 2. Check if the adjacent faces at each end are planar (indicating a through-hole)
/// 3. Check for conical or spherical bottom faces (indicating blind hole)
/// 4. Analyze edge connectivity to determine hole termination
pub fn classify_hole_type(brep: &rcad_kernel::BRep, feature: &CylindricalFeature) -> CylindricalFeatureExtended {
    let si = 0;
    let shi = 0;

    // Walk TShape hierarchy to find the first shell
    let Some(shell_ts_idx) = brep.tshapes.iter().find_map(|ts| {
        if let TShape::Solid(sd) = ts.as_ref() {
            sd.shells.first().copied()
        } else {
            None
        }
    })
    .map(|sr| sr.index)
    else {
        return CylindricalFeatureExtended {
            base: feature.clone(),
            hole_type: HoleType::Unknown,
            has_flat_bottom: false,
            has_conical_bottom: false,
            blind_depth: 0.0,
            bottom_face_index: None,
            top_adjacent_faces: Vec::new(),
            bottom_adjacent_faces: Vec::new(),
        };
    };
    let TShape::Shell(shd) = &*brep.tshapes[shell_ts_idx] else {
        return CylindricalFeatureExtended {
            base: feature.clone(),
            hole_type: HoleType::Unknown,
            has_flat_bottom: false,
            has_conical_bottom: false,
            blind_depth: 0.0,
            bottom_face_index: None,
            top_adjacent_faces: Vec::new(),
            bottom_adjacent_faces: Vec::new(),
        };
    };

    // Build edge -> face adjacency map (edge tshape index -> [local face index])
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face_sr) in shd.faces.iter().enumerate() {
        if let TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
            if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                for edge_sr in &wd.edges {
                    edge_to_faces.entry(edge_sr.index).or_default().push(fi);
                }
            }
            for inner_sr in &fd.inner_wires {
                if let TShape::Wire(iwd) = &*brep.tshapes[inner_sr.index] {
                    for edge_sr in &iwd.edges {
                        edge_to_faces.entry(edge_sr.index).or_default().push(fi);
                    }
                }
            }
        }
    }

    // Find adjacent faces at each end of the cylinder
    let ax = feature.axis;
    let mut top_adjacent: Vec<usize> = Vec::new();
    let mut bottom_adjacent: Vec<usize> = Vec::new();

    // Collect all edges from the cylindrical wall faces
    let mut wall_edges: HashSet<usize> = HashSet::new();
    for &fi in &feature.face_indices {
        if let Some(face_sr) = shd.faces.get(fi) {
            if let TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                if let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                    for edge_sr in &wd.edges {
                        wall_edges.insert(edge_sr.index);
                    }
                }
            }
        }
    }

    // For each wall edge, find adjacent non-wall faces
    for &ei in &wall_edges {
        if let Some(adj_faces) = edge_to_faces.get(&ei) {
            for &afi in adj_faces {
                // Skip if this is a wall face
                if feature.face_indices.contains(&afi) {
                    continue;
                }

                // Determine if this face is at top or bottom of the cylinder
                // by analyzing the vertex positions of the shared edge
                if let TShape::Edge(ed) = &*brep.tshapes[ei] {
                    let mid_point = if let (Some(p1), Some(p2)) = (
                        brep.vertex_point(ed.first.index),
                        brep.vertex_point(ed.last.index),
                    ) {
                        (p1 + p2) * 0.5
                    } else {
                        continue;
                    };

                    // Project onto axis to determine position
                    let t = (mid_point - feature.origin).dot(ax);

                    if t > (feature.t_min + feature.t_max) * 0.5 {
                        top_adjacent.push(afi);
                    } else {
                        bottom_adjacent.push(afi);
                    }
                }
            }
        }
    }

    // Remove duplicates
    top_adjacent.sort();
    top_adjacent.dedup();
    bottom_adjacent.sort();
    bottom_adjacent.dedup();

    // Analyze adjacent faces to determine hole type
    let mut has_flat_bottom = false;
    let mut has_conical_bottom = false;
    let mut bottom_face_index: Option<usize> = None;

    // Check for planar bottom face (blind hole indicator)
    let check_faces = if top_adjacent.is_empty() && !bottom_adjacent.is_empty() {
        &bottom_adjacent
    } else if !top_adjacent.is_empty() && bottom_adjacent.is_empty() {
        &top_adjacent
    } else {
        &bottom_adjacent
    };

    for &afi in check_faces {
        if let Some(plane) = face_plane(brep, si, shi, afi) {
            // Check if the plane normal is opposite to the cylinder axis (bottom face)
            let dot = plane.normal.dot(ax);
            if dot.abs() > 0.9 {
                has_flat_bottom = true;
                bottom_face_index = Some(afi);
            }
        }
        if let Some(cone) = face_cone(brep, si, shi, afi) {
            // Conical bottom (drill point)
            let cone_axis = cone.axis.normalize_or_zero();
            if cone_axis.dot(ax).abs() > 0.9 {
                has_conical_bottom = true;
                bottom_face_index = Some(afi);
            }
        }
        if let Some(_sphere) = face_sphere(brep, si, shi, afi) {
            // Spherical bottom (ball-end drill)
            has_conical_bottom = true;
            bottom_face_index = Some(afi);
        }
    }

    // Determine hole type based on analysis
    let (hole_type, blind_depth) = if has_flat_bottom || has_conical_bottom {
        let depth = feature.height();
        (HoleType::BlindHole, depth)
    } else if top_adjacent.is_empty() && bottom_adjacent.is_empty() {
        (HoleType::ThroughHole, 0.0)
    } else if top_adjacent.len() > 1 && bottom_adjacent.len() > 1 {
        (HoleType::ThroughHole, 0.0)
    } else {
        (HoleType::ThroughHole, 0.0)
    };

    CylindricalFeatureExtended {
        base: feature.clone(),
        hole_type,
        has_flat_bottom,
        has_conical_bottom,
        blind_depth,
        bottom_face_index,
        top_adjacent_faces: top_adjacent,
        bottom_adjacent_faces: bottom_adjacent,
    }
}

/// Detect and classify all cylindrical features in a B-Rep.
///
/// Returns extended features with hole type classification.
pub fn detect_cylindrical_features_extended(
    brep: &rcad_kernel::BRep,
    max_hole_radius: f64,
    max_boss_radius: f64,
) -> Vec<CylindricalFeatureExtended> {
    let base_features = detect_cylindrical_features(brep, max_hole_radius, max_boss_radius);
    base_features
        .into_iter()
        .map(|f| classify_hole_type(brep, &f))
        .collect()
}

// =============================================================================
// POST-SUPPRESSION TOPOLOGY HEALING
// =============================================================================

/// Result of post-suppression healing.
#[derive(Debug, Clone, Default)]
pub struct PostSuppressionHealingReport {
    /// Number of gaps filled.
    pub gaps_filled: usize,
    /// Number of dangling edges removed.
    pub dangling_edges_removed: usize,
    /// Number of tolerance mismatches repaired.
    pub tolerance_repairs: usize,
    /// Number of vertices merged.
    pub vertices_merged: usize,
    /// Number of degenerate faces removed.
    pub degenerate_faces_removed: usize,
    /// Number of healing passes performed.
    pub passes_performed: usize,
    /// Whether healing succeeded.
    pub success: bool,
}

/// Options for post-suppression healing.
#[derive(Debug, Clone)]
pub struct PostSuppressionHealingOptions {
    /// Tolerance for gap detection.
    pub gap_tolerance: f64,
    /// Tolerance for vertex merging.
    pub merge_tolerance: f64,
    /// Minimum edge length (edges below this are candidates for removal).
    pub min_edge_length: f64,
    /// Maximum number of healing passes.
    pub max_passes: usize,
    /// Whether to attempt gap filling.
    pub fill_gaps: bool,
    /// Whether to remove dangling edges.
    pub remove_dangling_edges: bool,
    /// Whether to repair tolerance mismatches.
    pub repair_tolerances: bool,
    /// Tolerance growth factor for each pass.
    pub tolerance_growth: f64,
    /// Maximum tolerance cap.
    pub tolerance_cap: f64,
}

impl Default for PostSuppressionHealingOptions {
    fn default() -> Self {
        Self {
            gap_tolerance: TOLERANCE_ABS * 10.0,
            merge_tolerance: TOLERANCE_ABS * 5.0,
            min_edge_length: TOLERANCE_ABS * 2.0,
            max_passes: 5,
            fill_gaps: true,
            remove_dangling_edges: true,
            repair_tolerances: true,
            tolerance_growth: 1.5,
            tolerance_cap: TOLERANCE_ABS * 100.0,
        }
    }
}

impl PostSuppressionHealingOptions {
    /// Create aggressive healing options for difficult cases.
    pub fn aggressive() -> Self {
        Self {
            gap_tolerance: TOLERANCE_ABS * 50.0,
            merge_tolerance: TOLERANCE_ABS * 20.0,
            min_edge_length: TOLERANCE_ABS * 5.0,
            max_passes: 10,
            fill_gaps: true,
            remove_dangling_edges: true,
            repair_tolerances: true,
            tolerance_growth: 2.0,
            tolerance_cap: TOLERANCE_ABS * 500.0,
        }
    }

    /// Create conservative healing options for precise geometry.
    pub fn conservative() -> Self {
        Self {
            gap_tolerance: TOLERANCE_ABS * 5.0,
            merge_tolerance: TOLERANCE_ABS * 2.0,
            min_edge_length: TOLERANCE_ABS,
            max_passes: 3,
            fill_gaps: false,
            remove_dangling_edges: true,
            repair_tolerances: true,
            tolerance_growth: 1.2,
            tolerance_cap: TOLERANCE_ABS * 20.0,
        }
    }
}

/// Perform post-suppression topology healing.
///
/// This function repairs the topology after feature suppression operations,
/// addressing gaps, dangling edges, and tolerance mismatches.
pub fn heal_after_suppression(
    brep: &rcad_kernel::BRep,
    options: &PostSuppressionHealingOptions,
) -> (rcad_kernel::BRep, PostSuppressionHealingReport) {
    let mut current = brep.clone();
    let mut report = PostSuppressionHealingReport::default();

    for pass in 0..options.max_passes {
        let growth = options.tolerance_growth.powi(pass as i32);
        let current_merge_tol = (options.merge_tolerance * growth).min(options.tolerance_cap);
        let current_gap_tol = (options.gap_tolerance * growth).min(options.tolerance_cap);

        let mut changed = false;

        // Step 1: Merge close vertices
        if options.repair_tolerances {
            let (merged_brep, merged_count) =
                crate::brep_repair::merge_close_vertices(&current, current_merge_tol);
            if merged_count > 0 {
                current = merged_brep;
                report.vertices_merged += merged_count;
                report.tolerance_repairs += merged_count;
                changed = true;
            }
        }

        // Step 2: Remove small/dangling edges
        if options.remove_dangling_edges {
            let (cleaned_brep, removed_count) =
                crate::brep_repair::remove_small_edges(&current, options.min_edge_length);
            if removed_count > 0 {
                current = cleaned_brep;
                report.dangling_edges_removed += removed_count;
                changed = true;
            }
        }

        // Step 3: Attempt gap filling (if enabled)
        if options.fill_gaps {
            let (filled_brep, gaps_filled) = fill_topology_gaps(&current, current_gap_tol);
            if gaps_filled > 0 {
                current = filled_brep;
                report.gaps_filled += gaps_filled;
                changed = true;
            }
        }

        // Step 4: Remove degenerate faces
        let (cleaned_brep, degenerate_count) =
            crate::brep_repair::remove_degenerate_faces(&current);
        if degenerate_count > 0 {
            current = cleaned_brep;
            report.degenerate_faces_removed += degenerate_count;
            changed = true;
        }

        report.passes_performed = pass + 1;

        if !changed {
            break;
        }
    }

    report.success = true;
    (current, report)
}

/// Fill topology gaps by analyzing edge connectivity.
///
/// Gaps can occur after boolean operations when faces don't align perfectly.
fn fill_topology_gaps(brep: &rcad_kernel::BRep, tolerance: f64) -> (rcad_kernel::BRep, usize) {
    let mut gaps_filled = 0usize;
    let mut current = brep.clone();

    // Build edge->face adjacency and find boundary edges (read-only TShape walk)
    let (edge_face_count, boundary_edges) = {
        let Some(shell_ts_idx) = current.tshapes.iter().find_map(|ts| {
            if let TShape::Solid(sd) = ts.as_ref() {
                sd.shells.first().copied()
            } else {
                None
            }
        })
        .map(|sr| sr.index)
        else {
            return (current, 0);
        };
        let TShape::Shell(shd) = &*current.tshapes[shell_ts_idx] else {
            return (current, 0);
        };

        // Count face usage for each edge
        let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
        for face_sr in &shd.faces {
            if let TShape::Face(fd) = &*current.tshapes[face_sr.index] {
                if let TShape::Wire(wd) = &*current.tshapes[fd.outer_wire.index] {
                    for edge_sr in &wd.edges {
                        *edge_face_count.entry(edge_sr.index).or_default() += 1;
                    }
                }
                for inner_sr in &fd.inner_wires {
                    if let TShape::Wire(iwd) = &*current.tshapes[inner_sr.index] {
                        for edge_sr in &iwd.edges {
                            *edge_face_count.entry(edge_sr.index).or_default() += 1;
                        }
                    }
                }
            }
        }

        // Find boundary edges (used by only one face in a manifold solid)
        let boundary_edges: Vec<usize> = edge_face_count
            .iter()
            .filter(|(_, count)| **count == 1)
            .map(|(ei, _)| *ei)
            .collect();

        (edge_face_count, boundary_edges)
    };
    // shd dropped here -- OK to mut borrow current.tshapes below

    // Collect vertex merge operations
    let mut vertex_merges: Vec<(usize, DVec3)> = Vec::new();

    // For each boundary edge, try to find and close the gap
    for &ei in &boundary_edges {
        let TShape::Edge(ed) = &*current.tshapes[ei] else {
            continue;
        };
        let Some(p1) = current.vertex_point(ed.first.index) else {
            continue;
        };
        let Some(p2) = current.vertex_point(ed.last.index) else {
            continue;
        };

        // Look for other edges with close vertices
        for (&other_ei, &count) in &edge_face_count {
            if other_ei == ei || count != 1 {
                continue;
            }
            let TShape::Edge(other_ed) = &*current.tshapes[other_ei] else {
                continue;
            };
            let Some(op1) = current.vertex_point(other_ed.first.index) else {
                continue;
            };
            let Some(op2) = current.vertex_point(other_ed.last.index) else {
                continue;
            };

            // Check if vertices are close enough to merge
            let close_1_1 = (p1 - op1).length() < tolerance;
            let close_1_2 = (p1 - op2).length() < tolerance;
            let close_2_1 = (p2 - op1).length() < tolerance;
            let close_2_2 = (p2 - op2).length() < tolerance;

            if (close_1_1 || close_1_2) && (close_2_1 || close_2_2) {
                // Record vertex merges
                if close_1_1 || close_1_2 {
                    let target_v = if close_1_1 { other_ed.first.index } else { other_ed.last.index };
                    vertex_merges.push((target_v, p1));
                }
                if close_2_1 || close_2_2 {
                    let target_v = if close_2_1 { other_ed.first.index } else { other_ed.last.index };
                    vertex_merges.push((target_v, p2));
                }
                gaps_filled += 1;
            }
        }
    }

    // Apply vertex merges via Arc::get_mut on TShape::Vertex
    for (vi, new_point) in vertex_merges {
        if let Some(TShape::Vertex(vd)) = current.tshapes.get_mut(vi).and_then(|arc| Arc::get_mut(arc)) {
            vd.point = new_point;
        }
    }

    (current, gaps_filled)
}

// =============================================================================
// FEATURE INTERACTION ANALYSIS
// =============================================================================

/// Type of interaction between two features.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureInteraction {
    /// Features share an edge.
    ShareEdge,
    /// Features share a vertex.
    ShareVertex,
    /// Features overlap spatially.
    Overlap,
    /// One feature is contained within another.
    Contained,
    /// Features are adjacent (within tolerance).
    Adjacent,
    /// Features do not interact.
    None,
}

/// Analysis result for feature interactions.
#[derive(Debug, Clone)]
pub struct FeatureInteractionAnalysis {
    /// Index of first feature.
    pub feature_a: usize,
    /// Index of second feature.
    pub feature_b: usize,
    /// Type of interaction detected.
    pub interaction: FeatureInteraction,
    /// Distance between features (for adjacent features).
    pub distance: f64,
    /// Whether features should be processed together.
    pub should_process_together: bool,
}

/// Analyze interactions between cylindrical features.
///
/// This function identifies pairs of features that share edges, vertices,
/// or overlap spatially, which should be processed together for robust defeaturing.
pub fn analyze_feature_interactions(
    brep: &rcad_kernel::BRep,
    features: &[CylindricalFeature],
    tolerance: f64,
) -> Vec<FeatureInteractionAnalysis> {
    let mut analyses: Vec<FeatureInteractionAnalysis> = Vec::new();

    for i in 0..features.len() {
        for j in (i + 1)..features.len() {
            let fa = &features[i];
            let fb = &features[j];

            // Check for shared faces/edges
            let share_edge = fa.face_indices.iter().any(|&fi_a| {
                fb.face_indices.iter().any(|&fi_b| {
                    faces_share_edge(brep, fi_a, fi_b)
                })
            });

            if share_edge {
                analyses.push(FeatureInteractionAnalysis {
                    feature_a: i,
                    feature_b: j,
                    interaction: FeatureInteraction::ShareEdge,
                    distance: 0.0,
                    should_process_together: true,
                });
                continue;
            }

            // Check for shared vertices
            let share_vertex = fa.face_indices.iter().any(|&fi_a| {
                fb.face_indices.iter().any(|&fi_b| {
                    faces_share_vertex(brep, fi_a, fi_b)
                })
            });

            if share_vertex {
                analyses.push(FeatureInteractionAnalysis {
                    feature_a: i,
                    feature_b: j,
                    interaction: FeatureInteraction::ShareVertex,
                    distance: 0.0,
                    should_process_together: true,
                });
                continue;
            }

            // Check for spatial overlap/adjacency
            let dist = feature_distance(fa, fb);
            if dist < tolerance {
                let interaction = if dist < 0.0 {
                    FeatureInteraction::Overlap
                } else if dist < tolerance * 0.1 {
                    FeatureInteraction::Contained
                } else {
                    FeatureInteraction::Adjacent
                };

                analyses.push(FeatureInteractionAnalysis {
                    feature_a: i,
                    feature_b: j,
                    interaction,
                    distance: dist,
                    should_process_together: true,
                });
            }
        }
    }

    analyses
}

/// Check if two faces share an edge.
fn faces_share_edge(brep: &rcad_kernel::BRep, fi_a: usize, fi_b: usize) -> bool {
    let Some(fd_a) = get_face_data(brep, 0, 0, fi_a) else {
        return false;
    };
    let Some(fd_b) = get_face_data(brep, 0, 0, fi_b) else {
        return false;
    };

    let edges_a: HashSet<usize> = {
        let mut s = HashSet::new();
        if let TShape::Wire(wd) = &*brep.tshapes[fd_a.outer_wire.index] {
            for edge_sr in &wd.edges {
                s.insert(edge_sr.index);
            }
        }
        s
    };
    let edges_b: HashSet<usize> = {
        let mut s = HashSet::new();
        if let TShape::Wire(wd) = &*brep.tshapes[fd_b.outer_wire.index] {
            for edge_sr in &wd.edges {
                s.insert(edge_sr.index);
            }
        }
        s
    };

    !edges_a.is_disjoint(&edges_b)
}

/// Check if two faces share a vertex.
fn faces_share_vertex(brep: &rcad_kernel::BRep, fi_a: usize, fi_b: usize) -> bool {
    let Some(fd_a) = get_face_data(brep, 0, 0, fi_a) else {
        return false;
    };
    let Some(fd_b) = get_face_data(brep, 0, 0, fi_b) else {
        return false;
    };

    let mut vertices_a: HashSet<usize> = HashSet::new();
    if let TShape::Wire(wd) = &*brep.tshapes[fd_a.outer_wire.index] {
        for edge_sr in &wd.edges {
            if let TShape::Edge(ed) = &*brep.tshapes[edge_sr.index] {
                vertices_a.insert(ed.first.index);
                vertices_a.insert(ed.last.index);
            }
        }
    }

    if let TShape::Wire(wd) = &*brep.tshapes[fd_b.outer_wire.index] {
        for edge_sr in &wd.edges {
            if let TShape::Edge(ed) = &*brep.tshapes[edge_sr.index] {
                if vertices_a.contains(&ed.first.index) || vertices_a.contains(&ed.last.index) {
                    return true;
                }
            }
        }
    }

    false
}

/// Compute the distance between two cylindrical features.
///
/// Returns a negative distance if features overlap spatially.
fn feature_distance(fa: &CylindricalFeature, fb: &CylindricalFeature) -> f64 {
    // Compute distance between feature axes
    let origin_diff = fb.origin - fa.origin;

    // Project onto both axes
    let proj_a = origin_diff.dot(fa.axis);
    let proj_b = origin_diff.dot(fb.axis);

    // Closest points on each axis
    let closest_a = fa.origin + fa.axis * proj_a.clamp(fa.t_min, fa.t_max);
    let closest_b = fb.origin + fb.axis * proj_b.clamp(fb.t_min, fb.t_max);

    // Distance between axes
    let axis_dist = (closest_b - closest_a).length();

    // Adjust for radii
    let radius_sum = fa.radius + fb.radius;
    axis_dist - radius_sum
}

/// Build a feature processing order that respects interactions.
///
/// Features that interact should be processed together or in sequence.
pub fn build_processing_order(
    features: &[CylindricalFeature],
    interactions: &[FeatureInteractionAnalysis],
) -> Vec<Vec<usize>> {
    // Build adjacency from interactions
    let mut adjacency: HashMap<usize, HashSet<usize>> = HashMap::new();
    for interaction in interactions {
        if interaction.should_process_together {
            adjacency
                .entry(interaction.feature_a)
                .or_default()
                .insert(interaction.feature_b);
            adjacency
                .entry(interaction.feature_b)
                .or_default()
                .insert(interaction.feature_a);
        }
    }

    // Find connected components
    let mut visited = vec![false; features.len()];
    let mut groups: Vec<Vec<usize>> = Vec::new();

    for start in 0..features.len() {
        if visited[start] {
            continue;
        }

        let mut group = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(idx) = queue.pop_front() {
            group.push(idx);

            if let Some(neighbors) = adjacency.get(&idx) {
                for &n in neighbors {
                    if !visited[n] {
                        visited[n] = true;
                        queue.push_back(n);
                    }
                }
            }
        }

        groups.push(group);
    }

    groups
}

// =============================================================================
// ROBUSTNESS IMPROVEMENTS
// =============================================================================

/// Robustness options for defeaturing operations.
#[derive(Debug, Clone)]
pub struct RobustnessOptions {
    /// Maximum number of attempts for each operation.
    pub max_attempts: usize,
    /// Tolerance growth factor for each retry.
    pub tolerance_growth: f64,
    /// Maximum tolerance cap.
    pub max_tolerance: f64,
    /// Whether to use fuzzy boolean operations.
    pub use_fuzzy_boolean: bool,
    /// Whether to heal between operations.
    pub heal_between_operations: bool,
    /// Healing options for inter-operation healing.
    pub healing_options: PostSuppressionHealingOptions,
}

impl Default for RobustnessOptions {
    fn default() -> Self {
        Self {
            max_attempts: 3,
            tolerance_growth: 2.0,
            max_tolerance: TOLERANCE_ABS * 100.0,
            use_fuzzy_boolean: true,
            heal_between_operations: true,
            healing_options: PostSuppressionHealingOptions::default(),
        }
    }
}

/// Result of a robust feature suppression operation.
#[derive(Debug, Clone)]
pub struct RobustSuppressionResult {
    /// The resulting BRep.
    pub brep: rcad_kernel::BRep,
    /// Whether the operation succeeded.
    pub success: bool,
    /// Number of attempts made.
    pub attempts: usize,
    /// Final tolerance used.
    pub final_tolerance: f64,
    /// Whether healing was applied.
    pub healing_applied: bool,
    /// Healing report (if healing was applied).
    pub healing_report: Option<PostSuppressionHealingReport>,
}

/// Attempt to suppress a feature with robust error recovery.
///
/// This function wraps the boolean operation with multiple retry strategies
/// and inter-operation healing.
pub fn suppress_feature_robust(
    brep: &rcad_kernel::BRep,
    fill_solid: &rcad_kernel::BRep,
    is_hole: bool,
    options: &RobustnessOptions,
) -> RobustSuppressionResult {
    let mut current = brep.clone();
    let mut tolerance = TOLERANCE_ABS;
    let mut healing_applied = false;
    let mut healing_report: Option<PostSuppressionHealingReport> = None;

    let op = if is_hole {
        BooleanOpType::Union
    } else {
        BooleanOpType::Difference
    };

    for attempt in 0..options.max_attempts {
        // Try the boolean operation
        let result = if options.use_fuzzy_boolean && attempt > 0 {
            // Use fuzzy tolerance for retry
            let fuzzy_opts = BooleanOptions {
                fuzzy_tol: tolerance,
                use_glue: true,
                glue_tolerance: tolerance,
                ..Default::default()
            };
            boolean_op_with_options(op, &current, fill_solid, fuzzy_opts)
        } else {
            boolean_op(op, &current, fill_solid).map(|t| (t).clone())
        };

        match result {
            Ok(new_brep) => {
                return RobustSuppressionResult {
                    brep: new_brep,
                    success: true,
                    attempts: attempt + 1,
                    final_tolerance: tolerance,
                    healing_applied,
                    healing_report,
                };
            }
            Err(_) => {
                // Try healing before retry
                if options.heal_between_operations {
                    let (healed, heal_report) =
                        heal_after_suppression(&current, &options.healing_options);
                    current = healed;
                    healing_applied = true;
                    healing_report = Some(heal_report);
                }

                // Increase tolerance for next attempt
                tolerance = (tolerance * options.tolerance_growth).min(options.max_tolerance);
            }
        }
    }

    RobustSuppressionResult {
        brep: current,
        success: false,
        attempts: options.max_attempts,
        final_tolerance: tolerance,
        healing_applied,
        healing_report,
    }
}

/// Perform boolean operation with explicit options.
fn boolean_op_with_options(
    op: BooleanOpType,
    a: &rcad_kernel::BRep,
    b: &rcad_kernel::BRep,
    options: BooleanOptions,
) -> Result<rcad_kernel::BRep, crate::BooleanError> {
    // For now, delegate to the standard boolean with fuzzy tolerance
    // A full implementation would respect all options
    if options.fuzzy_tol > 0.0 {
        let robust_opts = BooleanRobustOptions {
            base: options,
            fuzzy_retry_ladder: vec![options.fuzzy_tol],
            retry_policy: BooleanRetryPolicy::AdaptiveByFailureClass,
            extreme_geometry: crate::ExtremeGeometryRetryConfig::default(),
        };
        boolean_op_robust(op, a, b, robust_opts).map(|(b, _)| b)
    } else {
        boolean_op(op, a, b).map(|t| (t).clone())
    }
}

// =============================================================================
// ENHANCED DEFEATURE WITH ALL IMPROVEMENTS
// =============================================================================

/// Enhanced defeaturing options with all improvements integrated.
#[derive(Debug, Clone)]
pub struct DefeaturingOptionsV2 {
    /// Base defeaturing options.
    pub base: DefeaturingOptions,
    /// Robustness options.
    pub robustness: RobustnessOptions,
    /// Post-suppression healing options.
    pub healing: PostSuppressionHealingOptions,
    /// Whether to classify hole types.
    pub classify_hole_types: bool,
    /// Whether to analyze feature interactions.
    pub analyze_interactions: bool,
    /// Whether to process interacting features together.
    pub process_interactions_together: bool,
    /// Interaction tolerance.
    pub interaction_tolerance: f64,
}

impl Default for DefeaturingOptionsV2 {
    fn default() -> Self {
        Self {
            base: DefeaturingOptions::default(),
            robustness: RobustnessOptions::default(),
            healing: PostSuppressionHealingOptions::default(),
            classify_hole_types: true,
            analyze_interactions: true,
            process_interactions_together: true,
            interaction_tolerance: TOLERANCE_ABS * 10.0,
        }
    }
}

impl DefeaturingOptionsV2 {
    /// Create options for simulation preprocessing.
    pub fn for_simulation() -> Self {
        Self {
            base: DefeaturingOptions {
                max_hole_radius: 5.0,
                max_boss_radius: 3.0,
                enable_conical_features: true,
                max_conical_hole_radius: 5.0,
                enable_retry: true,
                max_retries: 5,
                run_post_healing: false, // We use our own healing
                ..Default::default()
            },
            robustness: RobustnessOptions {
                max_attempts: 5,
                tolerance_growth: 1.5,
                heal_between_operations: true,
                ..Default::default()
            },
            healing: PostSuppressionHealingOptions::aggressive(),
            classify_hole_types: true,
            analyze_interactions: true,
            process_interactions_together: true,
            interaction_tolerance: 0.01,
        }
    }

    /// Create options for machining preparation.
    pub fn for_machining() -> Self {
        Self {
            base: DefeaturingOptions {
                max_hole_radius: 0.0, // Don't remove holes for machining
                max_boss_radius: 2.0,
                enable_blend_features: true,
                max_blend_radius: 1.0,
                max_chamfer_distance: 1.0,
                enable_retry: true,
                ..Default::default()
            },
            robustness: RobustnessOptions {
                max_attempts: 3,
                tolerance_growth: 1.2,
                heal_between_operations: false,
                ..Default::default()
            },
            healing: PostSuppressionHealingOptions::conservative(),
            classify_hole_types: true,
            analyze_interactions: false,
            process_interactions_together: false,
            interaction_tolerance: 0.001,
        }
    }
}

/// Enhanced report with all analysis details.
#[derive(Debug, Clone, Default)]
pub struct DefeaturingReportV2 {
    /// Base report.
    pub base: DefeaturingReport,
    /// Classified hole types.
    pub hole_types: Vec<(usize, HoleType)>,
    /// Feature interactions detected.
    pub interactions: Vec<FeatureInteractionAnalysis>,
    /// Processing groups.
    pub processing_groups: Vec<Vec<usize>>,
    /// Post-suppression healing report.
    pub healing_report: Option<PostSuppressionHealingReport>,
    /// Robustness statistics.
    pub total_attempts: usize,
    pub features_succeeded_on_retry: usize,
}

/// Perform enhanced defeaturing with all improvements.
///
/// This function integrates:
/// - Through-hole vs blind-hole classification
/// - Feature interaction analysis
/// - Robust error recovery
/// - Post-suppression topology healing
pub fn defeature_brep_v2(
    brep: &rcad_kernel::BRep,
    options: &DefeaturingOptionsV2,
) -> Result<(rcad_kernel::BRep, DefeaturingReportV2), DefeaturingError> {
    // Check that there's at least one solid with a non-empty shell
    let has_viable_solid = brep.tshapes.iter().any(|ts| {
        if let TShape::Solid(sd) = ts.as_ref() {
            !sd.shells.is_empty()
        } else {
            false
        }
    });
    if !has_viable_solid {
        return Err(DefeaturingError::EmptyInput);
    }

    let mut report = DefeaturingReportV2::default();
    let mut current = brep.clone();

    // Step 1: Detect cylindrical features
    let features = if options.base.max_hole_radius > 0.0 || options.base.max_boss_radius > 0.0 {
        detect_cylindrical_features(
            &current,
            options.base.max_hole_radius,
            options.base.max_boss_radius,
        )
    } else {
        Vec::new()
    };

    // Step 2: Classify hole types if requested
    let extended_features = if options.classify_hole_types {
        let extended: Vec<CylindricalFeatureExtended> = features
            .iter()
            .map(|f| classify_hole_type(&current, f))
            .collect();

        // Record classifications
        for (i, ext) in extended.iter().enumerate() {
            report.hole_types.push((i, ext.hole_type));
        }

        extended
    } else {
        features
            .iter()
            .map(|f| CylindricalFeatureExtended {
                base: f.clone(),
                hole_type: HoleType::Unknown,
                has_flat_bottom: false,
                has_conical_bottom: false,
                blind_depth: 0.0,
                bottom_face_index: None,
                top_adjacent_faces: Vec::new(),
                bottom_adjacent_faces: Vec::new(),
            })
            .collect()
    };

    // Step 3: Analyze feature interactions if requested
    let processing_groups = if options.analyze_interactions && !extended_features.is_empty() {
        let interactions = analyze_feature_interactions(
            &current,
            &extended_features.iter().map(|e| e.base.clone()).collect::<Vec<_>>(),
            options.interaction_tolerance,
        );
        report.interactions = interactions.clone();

        if options.process_interactions_together {
            build_processing_order(
                &extended_features.iter().map(|e| e.base.clone()).collect::<Vec<_>>(),
                &interactions,
            )
        } else {
            (0..extended_features.len()).map(|i| vec![i]).collect()
        }
    } else {
        (0..extended_features.len()).map(|i| vec![i]).collect()
    };
    report.processing_groups = processing_groups.clone();

    // Step 4: Process features with robust suppression
    let margin = if options.base.fill_margin > 0.0 {
        options.base.fill_margin
    } else {
        DEFAULT_FILL_MARGIN
    };

    for group in &processing_groups {
        for &idx in group {
            let ext_feature = &extended_features[idx];
            let feature = &ext_feature.base;

            // Determine if this feature should be processed
            let should_process = if feature.is_hole {
                options.base.max_hole_radius > 0.0 && feature.radius <= options.base.max_hole_radius
            } else {
                options.base.max_boss_radius > 0.0 && feature.radius <= options.base.max_boss_radius
            };

            if !should_process {
                continue;
            }

            // Build fill solid
            let fill_result = if feature.is_hole {
                make_fill_cylinder(feature, margin)
            } else {
                make_boss_cylinder(feature, margin)
            };

            let Ok(fill) = fill_result else {
                report.base.failed_features += 1;
                continue;
            };

            // Apply robust suppression
            let fill_old = (fill).clone();
            let result = suppress_feature_robust(&current, &fill_old, feature.is_hole, &options.robustness);

            report.total_attempts += result.attempts;
            if result.success {
                current = result.brep;
                if feature.is_hole {
                    report.base.holes_removed += 1;
                } else {
                    report.base.bosses_removed += 1;
                }
                if result.attempts > 1 {
                    report.features_succeeded_on_retry += 1;
                }
            } else {
                report.base.failed_features += 1;
            }
        }
    }

    // Step 5: Post-suppression healing
    let (healed, healing_report) = heal_after_suppression(&current, &options.healing);
    current = healed;
    report.healing_report = Some(healing_report);

    Ok((current, report))
}
